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
#include "factor.h"
 
/*****************************************************************************/
// Functions for class Factor
/*****************************************************************************/
 
bool Factor::isApproxNetwork__ = false;
     
void Factor::setIsApproxNetwork(bool val)
{
  isApproxNetwork__ = val;
}


//the contribution of the factor itself
void Factor::initFactorMesssages()
{
  Predicate *pred;
  int numPreds = clause_->getNumPredicates();
  int stateCnt = (int)pow(2.0, numPreds);
  factorMsgs_ = new double[stateCnt];
  bool isSatisfied;
  for (int state = 0; state < stateCnt; state++)
  {
    isSatisfied = false;
    for (int predno = 0; predno < numPreds; predno++)
    {
      pred = clause_->getPredicate(predno);
      bool predBit = state & (1<<predno);
      if (pred->getSense() == predBit)
      {
        isSatisfied = true;
        break;
      }
    }

    if (isSatisfied)
    {
        // Always 1
        // MS: should be weight of clause in ground network (super or ground)
      factorMsgs_[state] = clause_->getWt();
    }
    else
    {
      factorMsgs_[state] = 0;
    }
  }
}

   
    //find the outgoing message for the given inpPredIndex
    double* Factor::multiplyMessagesAndSumOut(int inpPredIndex) {
             //cout<<"Doing for pred "<<inpPredIndex<<endl;
             int numPreds = clause_->getNumPredicates();
             int stateCnt = (int)pow(2.0,numPreds);
             double * prodMsgs = new double[stateCnt];
             
             Node *node;

			 double *gndNodeCnts = NULL;

			 /* This is done to handle the case when BP is run
             * on supernodes/superfeatures which have not yet
             * reached an equilibrium state */
             isApproxNetwork__ = false;
			 if(isApproxNetwork__) {
			  gndNodeCnts = new double[numPreds];
              for(int predno=0;predno<numPreds;predno++)
                  gndNodeCnts[predno] = 0;
             
              for(int lno=0;lno<links_->size();lno++) {
                  int predIndex = (*links_)[lno]->getPredIndex();
                  node = (*links_)[lno]->getNode();
                  gndNodeCnts[predIndex] += node->getGroundNodeCount();
              }
			 }

             //initialize the product
             for(int state=0;state<stateCnt;state++)
                  prodMsgs[state] = factorMsgs_[state];

             if(!action_)
             {
             for(int lno=0;lno<links_->size();lno++) {
                  int predIndex = (*links_)[lno]->getPredIndex();
                  if(predIndex == inpPredIndex)
                       continue;
                  
                  node = (*links_)[lno]->getNode();
                  double wt;

				  //wt must be equal to 1 in equilibrium. This step 
                  //averages (weighted) out the messages from various 
                  //supernodes (at this predIndex position).
				  
				  if(isApproxNetwork__) {
				   assert(gndNodeCnts[predIndex] != 0);
                   wt = node->getGroundNodeCount()/gndNodeCnts[predIndex];
				  } else { 
				   wt = 1; 
				  }
                  
				  for(int state=0;state<stateCnt;state++) {
                      bool predBit = state & (1<<predIndex);
                      if(predBit) {
                            prodMsgs[state] += (*msgsArr_)[lno][1]*wt;
                      } else
                            prodMsgs[state] += (*msgsArr_)[lno][0]*wt;
                 }
             }
             }
             
             /*
             for(int state=0;state<stateCnt;state++) {
                  cout<<"ProdMsgs["<<state<<"] = "<<prodMsgs[state]<<endl;
             }*/

             //caller is responsible for deleting it
             double *outMsgs = new double[2];    
             double maxMsgs[2];
             maxMsgs[0] = maxMsgs[1] = 0;
             bool firstTime[2];
             firstTime[0] = firstTime[1] = true;
             
             //now find the max messages
             /*
             for(int state=0;state<stateCnt;state++) {
                 if(maxMsg < prodMsgs[state]) {
                          maxMsg = prodMsgs[state];
                 }
             }*/

             //now find the max messages 
             for(int state=0;state<stateCnt;state++) {
                 bool predBit = state & (1<<inpPredIndex);
                 if(predBit && (maxMsgs[1] < prodMsgs[state] || firstTime[1])) {
                           firstTime[1] = false;
                           maxMsgs[1] = prodMsgs[state];
                 }

                 if(!predBit && (maxMsgs[0] < prodMsgs[state] || firstTime[0])) {
                           firstTime[0] = false;
                           maxMsgs[0] = prodMsgs[state];
                 }
             }

             outMsgs[0] = outMsgs[1] = 0;
             for(int state=0;state<stateCnt;state++) {
                 bool predBit = state & (1<<inpPredIndex);
                 /*
                 double msgDiff = prodMsgs[state] - maxMsg;
                 if(msgDiff < MINEXP)
                      msgDiff = MINEXP;
                  */
                 if(predBit) 
                      outMsgs[1] += expl(prodMsgs[state] - maxMsgs[1]);
                 else
                      outMsgs[0] += expl(prodMsgs[state] - maxMsgs[0]);
             }
             //outMsgs[1] = maxMsgs[1] + logl(outMsgs[1]);
             //outMsgs[0] = maxMsgs[0] + logl(outMsgs[0]);
             outMsgs[1] = maxMsgs[1] + logl(outMsgs[1]);
             outMsgs[0] = maxMsgs[0] + logl(outMsgs[0]);
             outMsgs[1] = outMsgs[1] - outMsgs[0];
             outMsgs[0] = 0;
             //cout<<"maxMsgs[0] = "<<maxMsgs[0]<<", maxMsgs[1] = "<<maxMsgs[1]<<endl;
             //cout<<"outMsgs[0] = "<<outMsgs[0]<<", outMsgs[1] = "<<outMsgs[1]<<endl;
             delete [] prodMsgs;
             if(isApproxNetwork__)
			  delete [] gndNodeCnts;
             
			 return outMsgs;
    }

   //send Message on all the links
void Factor::sendMessage()
{
  double *outMsgs = NULL;
  Link *link;
  Node *node;
  int predIndex;
  for (int lindex = 0; lindex < links_->size(); lindex++)
  {
    link = (*links_)[lindex];
    node = link->getNode();
    predIndex = link->getPredIndex();

    outMsgs = multiplyMessagesAndSumOut(predIndex);
      //Assumes pass by value copy of the messages
    node->receiveMessage(outMsgs, link);
    delete [] outMsgs;
  }
}

  //update the stored msgs and update the msgProduct
void Factor::moveToNextStep()
{
  double * msgs;
  double * nextMsgs;

    //store the current messages in prevMessages array and
    //reinitialize the current message array
  for (int lindex = 0;lindex < links_->size(); lindex++)
  {
    msgs = (*msgsArr_)[lindex];
    nextMsgs = (*nextMsgsArr_)[lindex];
    for (int i = 0; i < 2; i++)
    {
      msgs[i] = nextMsgs[i];
        // CHECK: commented out in dt
      nextMsgs[i] = 0;
    }
  }
}

void Factor::tagNeighborsSendReceive(set<Factor*> &sendFactors,
                                     set<Factor*> &receiveFactors,
                                     set<Node*> &sendNodes,
                                     set<Node*> &receiveNodes)
{
  //if(taggedSend_)
  //  return;
  //taggedSend_ = true;

  Node* node;
  for (int lindex = 0; lindex < links_->size(); lindex++)
  {
    node = (*links_)[lindex]->getNode();
    sendNodes.insert(node);
    receiveNodes.insert(node);
    node->tagNeighborsSend(sendFactors, receiveFactors, receiveNodes);
  }
}

void Factor::tagNeighborsSendReceive(list<Factor*> &sendFactors,
                                     list<Factor*> &receiveFactors,
                                     list<Node*> &sendNodes,
                                     list<Node*> &receiveNodes)
{
  //if(taggedSend_)
  //  return;
  //taggedSend_ = true;

  Node* node;
  for (int lindex = 0; lindex < links_->size(); lindex++)
  {
    node = (*links_)[lindex]->getNode();
    sendNodes.push_back(node);
    receiveNodes.push_back(node);
    node->tagNeighborsSend(sendFactors, receiveFactors, receiveNodes);
  }
}

void Factor::tagNeighborsReceive(set<Node*> &receiveNodes)
{
  //if(taggedReceive_)
  //  return;
  //taggedReceive_ = true;

  Node* node;
  for (int lindex = 0; lindex < links_->size(); lindex++)
  {
    node = (*links_)[lindex]->getNode();
    receiveNodes.insert(node);
  }
}

void Factor::tagNeighborsReceive(list<Node*> &receiveNodes)
{
  //if(taggedReceive_)
  //  return;
  //taggedReceive_ = true;

  Node* node;
  for (int lindex = 0; lindex < links_->size(); lindex++)
  {
    node = (*links_)[lindex]->getNode();
    receiveNodes.push_back(node);
  }
}



/**
 * Prints the factor as its variables separated by "/".
 */
ostream& Factor::print(ostream& out)
{
  /*
  if (superClause_ != NULL)
  {
    printArray(*(superClause_->getConstantTuple(0)), 1, out);
  }
  else
  {
    printArray(*constants_,1,out);
  }
  out << "// ";
  */
    
  for (int i = 0; i < links_->size(); i++)
  {
    (*links_)[i]->getNode()->print(out);
    if (i != links_->size() - 1) out << "/ ";
  }
  
  return out;
}

/**
 * Prints the weighted features for all states of the variables.
 */
ostream& Factor::printWts(ostream& out)
{
  for (int i = 0; i < (int)pow(2.0, clause_->getNumPredicates()); i++)
  {
    out << factorMsgs_[i] * outputWt_ << " ";
  }
  
  return out;
}

