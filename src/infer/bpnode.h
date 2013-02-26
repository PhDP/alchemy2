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
#ifndef _BPNODE_H_Jan_2008
#define _BPNODE_H_Jan_2008

#include <math.h>
#include "util.h"
#include "mrf.h"
#include "array.h"
#include "bplink.h"
#include "superpred.h"
#include "twowaymessage.h"
#include "bpfactor.h"

const double MIN = 1e-6;
const double MINEXP = logl(MIN);

const double SMOOTHTERM = 1e-6;

using namespace std;
using namespace __gnu_cxx;

class BPNode
{
 public:
  BPNode(int predId, SuperPred * const & superPred, 
         Array<int> * const & constants, Domain * const & domain)
  {
    predId_ = predId;
    superPred_ = superPred;
    constants_ = constants;
    domain_ = domain;
    links_ = new Array<BPLink *>();
    auxLinks_ = new Array<BPLink *>();
    msgsArr_ = new Array<double *>();
    nextMsgsArr_ = new Array<double *>();
    msgProds_[0] = msgProds_[1] = 0;
  }

  ~BPNode()
  {
    for (int i = 0; i < links_->size(); i++)
      delete (*links_)[i];
    for (int i = 0; i < auxLinks_->size(); i++)
      delete (*auxLinks_)[i];
    for (int i = 0; i < msgsArr_->size(); i++)
    {
      delete (*msgsArr_)[i];
      delete (*nextMsgsArr_)[i];
    }
    delete links_;
    delete auxLinks_;
    delete msgsArr_;
    delete nextMsgsArr_;
  }

  int getPredId() { return predId_;}

  int getSuperPredId()
  {
    if (superPred_) 
      return superPred_->getSuperPredId();
    else
      return -1;
  }

  int getParentSuperPredId()
  {
    if (superPred_) 
      return superPred_->getParentSuperPredId();
    else
      return -1;
  }

  SuperPred* getSuperPred() {return superPred_;}
  Array<int> * getConstants() {return constants_;}

  int getGroundNodeCount()
  {
    if (superPred_)
    {
      return superPred_->getNumTuples();
    }
    else
      return 1;
  }

  int getNumLinks() { return links_->size();}
  
  BPLink * getLink(int index) {return (*links_)[index];}
          
  void getMessage(int index, double msgs[])
  {
    msgs[0] = (*msgsArr_)[index][0];
    msgs[1] = (*msgsArr_)[index][1];
  }

 /**
  * Adds an auxiliary link (factor along with all the relevant information about
  * counts etc.).
  */
  void addAuxLink(BPLink *link)
  {
    auxLinks_->append(link);
  }

 /**
  * add the link (factor along with all the relevant information about
  * counts etc.).
  */
  void addLink(BPLink *link, double inpMsgs[2])
  {
    links_->append(link);
    double *msgs;
    msgs = new double[2];

    double cnt = link->getCount(); 
    for (int i = 0; i < 2; i++)
    {
      if (inpMsgs)
      {
        msgs[i] = inpMsgs[i];
      }
      else
      {
        msgs[i] = 0;
      }
      msgProds_[i] = msgProds_[i] + cnt*msgs[i];
    }

    msgsArr_->append(msgs);
    msgs = new double[2];
    nextMsgsArr_->append(msgs);
  }

    //add the factors with appropriate counts, also add the node to the
    //corresponding factor
  void addFactors(Array<BPFactor *> * const & allFactors,
                  LinkIdToTwoWayMessageMap* const & lidToTWMsg);

    //receive the message sent over a link
  void receiveMessage(double* inpMsgs, BPLink * const & link)
  {
    double *nextMsgs;
    int reverseFactorIndex = link->getReverseFactorIndex();
    nextMsgs = (*nextMsgsArr_)[reverseFactorIndex];
    nextMsgs[0] = inpMsgs[0];
    nextMsgs[1] = inpMsgs[1];
  }

  double getExp()
  {
    double exp = msgProds_[1] - msgProds_[0];
    return exp;
  }

  /**
   * Gets the probabilities for this node - it is simply the product of
   * messages coming at this node
   */
  double * getProbs(double * const & probs)
  {
    double exps[2];
    double norm, smooth;
    exps[0] = 0; //msgProds_[0];
    exps[1] = msgProds_[1] - msgProds_[0];

    for (int i = 0; i < 2; i++)
    {
      if(exps[i] < MINEXP)
        exps[i] = MINEXP;
      if(exps[i] > -MINEXP)
        exps[i] = -MINEXP;
      probs[i] = expl(exps[i]);
    }

    norm = probs[0] + probs[1];
    smooth = norm*SMOOTHTERM;
    norm += smooth;

    for (int i = 0; i < 2 ; i++)
    {
      probs[i] += smooth/2;
      probs[i] = probs[i]/norm;
    }
    return probs;	
  }

    //send the messages to all the factor nodes connected to this node
  void sendMessage();

    //send the messages to all the auxiliary factor nodes connected to this node
  void sendAuxMessage();

    //update the stored msgs and update the msgProduct
  void moveToNextStep();

  ostream& print(ostream& out)
  {
    out << predId_ << ": ";
    if (superPred_ != NULL)
      printArray(*(superPred_->getConstantTuple(0)),out);
    else
      printArray(*constants_,out);
    return out;
  }

 private:

  int predId_;

    //only one of the two below is used - superPred in case of lifted
    //inference and constants in case of ground inference
  SuperPred * superPred_;
  Array<int> *constants_;

  Domain * domain_;

  Array<BPLink *> * links_;
    // Links to auxiliary factors
  Array<BPLink *> * auxLinks_;
    //log of the actual message received from factor is stored
  Array<double *> * msgsArr_; 
  Array<double *> * nextMsgsArr_;

    //we store these so that we do not have to n^2 squared
    //computation while calculating the outgoing messages - we can 
    //just divide (subtract in the log domain) by the message from the 
    //factor node to which it is being sent
  double msgProds_[2];
};

#endif

