/*
 * All of the documentation and software included in the
 * Alchemy Software is copyrighted by Stanley Kok, Parag
 * Singla, Matthew Richardson, Pedro Domingos, Marc
 * Sumner, Hoifung Poon, Daniel Lowd, and Jue Wang.
 * 
 * Copyright [2004-09] Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang. All rights reserved.
 * 
 * Contact: Pedro Domingos, University of Washington
 * (pedrod@cs.washington.edu).
 * 
 * Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that
 * the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the
 * following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the
 * above copyright notice, this list of conditions and the
 * following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 * 
 * 3. All advertising materials mentioning features or use
 * of this software must display the following
 * acknowledgment: "This product includes software
 * developed by Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang in the Department of
 * Computer Science and Engineering at the University of
 * Washington".
 * 
 * 4. Your publications acknowledge the use or
 * contribution made by the Software to your research
 * using the following citation(s): 
 * Stanley Kok, Parag Singla, Matthew Richardson and
 * Pedro Domingos (2005). "The Alchemy System for
 * Statistical Relational AI", Technical Report,
 * Department of Computer Science and Engineering,
 * University of Washington, Seattle, WA.
 * http://alchemy.cs.washington.edu.
 * 
 * 5. Neither the name of the University of Washington nor
 * the names of its contributors may be used to endorse or
 * promote products derived from this software without
 * specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF WASHINGTON
 * AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OF WASHINGTON OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
#include "bpnode.h"

/*****************************************************************************/
// Functions for class BPNode
/*****************************************************************************/

  //add the factors with appropriate counts, also add the node to the
  //corresponding factor
void BPNode::addFactors(Array<BPFactor *> * const & allFactors,
                        LinkIdToTwoWayMessageMap* const & lidToTWMsg)
{
  BPFactor *factor;
  BPLink *link;
  double cnt;
  ClauseCounter * counter = (ClauseCounter *)superPred_->getClauseCounter();
  int numFactors = counter->getNumClauses();
  const Array<double> * cnts;

  for (int findex = 0; findex < numFactors; findex++)
  {
    int fid = counter->getClauseId(findex);  
    cnts = counter->getClauseCounts(findex);
    factor = (*allFactors)[fid];
      //predIndex is the predicate index in the clause/factor
    for (int predIndex = 0; predIndex < cnts->size(); predIndex++)
    {
      cnt = (*cnts)[predIndex]; 
      if ((*cnts)[predIndex] == 0)
        continue;

        //index where this node would be stored in the list of factors
      int reverseNodeIndex = factor->getNumLinks();
        //index where this factor would be stored in the list of nodes
      int reverseFactorIndex = getNumLinks();

      link = new BPLink(this, factor, reverseNodeIndex, reverseFactorIndex,
                        predIndex, cnt);

        //now find the messages from parent nodes
      LinkId *lid;
      TwoWayMessage *tmsg;
      double *nodeToFactorMsgs, *factorToNodeMsgs;
      LinkIdToTwoWayMessageMap::iterator lidToTMsgItr;
      int parentSuperPredId = getParentSuperPredId();
      int parentSuperClauseId = factor->getParentSuperClauseId();

      lid = new LinkId(predId_, parentSuperPredId, parentSuperClauseId,
                       predIndex);
      lidToTMsgItr = lidToTWMsg->find(lid);
      delete lid;

      if (lidToTMsgItr != lidToTWMsg->end())
      {
        tmsg = lidToTMsgItr->second;
        nodeToFactorMsgs = tmsg->getNodeToFactorMessage();
        factorToNodeMsgs = tmsg->getFactorToNodeMessage();
      }
      else
      {
        nodeToFactorMsgs = NULL;
        factorToNodeMsgs = NULL;
      }
      this->addLink(link, nodeToFactorMsgs);
      factor->addLink(link,factorToNodeMsgs);
    }
  }
}

  //send the messages to all the factor nodes connected to this node
void BPNode::sendMessage()
{
  BPLink *link;
  BPFactor *factor;
  double cnt;
  double *msgs;
  double *outMsgs = new double[2];

  for (int lindex = 0; lindex < links_->size(); lindex++)
  {
    link = (*links_)[lindex];
    factor = link->getFactor();
    cnt = link->getCount();
      //subtract the msg recieved from this factor node
    msgs = (*msgsArr_)[lindex];
    for (int i = 0; i < 2; i++)
    {
      outMsgs[i] = msgProds_[i] - msgs[i];
    }

      //Assumes pass by value copy of the messages
    factor->receiveMessage(outMsgs, link);
  }
  delete [] outMsgs;
}

  //send the messages to all the auxiliary factor nodes connected to this node
void BPNode::sendAuxMessage()
{
  BPLink *link;
  BPFactor *factor;
  double cnt;
  //double *msgs;
  double *outMsgs = new double[2];

  for (int lindex = 0; lindex < auxLinks_->size(); lindex++)
  {
    link = (*auxLinks_)[lindex];
    factor = link->getFactor();
    cnt = link->getCount();
    for (int i = 0; i < 2; i++)
    {
      outMsgs[i] = msgProds_[i];
    }

      //Assumes pass by value copy of the messages
    factor->receiveMessage(outMsgs, link);
  }
  delete [] outMsgs;
}


  //update the stored msgs and update the msgProduct
void BPNode::moveToNextStep()
{
  double * msgs;
  double * nextMsgs;
  double cnt;

    //initialize to zero
  for (int i = 0; i < 2; i++)
  {
    msgProds_[i] = 0;
  }

    //store the current messages in prevMessages array and
    //reinitialize the current message array
  for (int lindex = 0; lindex < links_->size(); lindex++)
  {
    msgs = (*msgsArr_)[lindex];
    nextMsgs = (*nextMsgsArr_)[lindex];

      // MS: Weight of clause, but should be count of nodes belonging to supernode
    if (superPred_)
      cnt = superPred_->getNumTuples();
    else
      cnt = 1;
    //cnt = (*links_)[lindex]->getCount();
    for (int i = 0; i < 2; i++)
    {
      msgs[i] = nextMsgs[i];
      nextMsgs[i] = 0;
        //since we store the logs, we need to take the sum
      msgProds_[i] = msgProds_[i] + cnt*msgs[i];
    }
  }
}

