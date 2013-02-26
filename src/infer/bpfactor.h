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
#ifndef _BPFACTOR_H_Jan_2008
#define _BPFACTOR_H_Jan_2008

#include <math.h>
#include "util.h"
#include "mrf.h"
#include "array.h"
#include "bplink.h"
#include "bpnode.h"
#include "superclause.h"

using namespace std;
using namespace __gnu_cxx;

/**
 * Stores all relevant information about a factor.
 */
class BPFactor
{
 public:
  BPFactor(Clause * const & clause, SuperClause * const & superClause, 
           Array<int> * const & constants, Domain * const & domain,
           double outputWt)
  {
    clause_ = clause;
    superClause_ = superClause;
    constants_ = constants;
    domain_ = domain;

    msgsArr_ = new Array<double *>;
    nextMsgsArr_ = new Array<double *>;
    links_ = new Array<BPLink *>();
    outputWt_ = outputWt;
    initFactorMesssages();
  }
  
  BPFactor()
  {
    links_ = new Array<BPLink *>();
    msgsArr_ = new Array<double *>;
  }
		  
  ~BPFactor()
  {
    for (int i = 0; i < msgsArr_->size(); i++)
    {
      delete (*msgsArr_)[i];
      delete (*nextMsgsArr_)[i];
    }
    delete links_;
    delete msgsArr_;
    delete nextMsgsArr_;
  }

  void initFactorMesssages();
  
  int getSuperClauseId()
  {
    if (superClause_) 
      return superClause_->getSuperClauseId();
    else
      return -1;
  }
		  
  int getParentSuperClauseId()
  {
    if( superClause_) 
      return superClause_->getParentSuperClauseId();
    else
      return -1;
  }

  SuperClause* getSuperClause() {return superClause_;}
  Clause * getClause() {return clause_;}
  Domain *getDomain() {return domain_;}
  Array<int> * getConstants() {return constants_;}
  int getNumLinks() {return links_->size();}
		  
  void getMessage(int index, double msgs[])
  { 
    msgs[0] = (*msgsArr_)[index][0];
    msgs[1] = (*msgsArr_)[index][1];
  }

  void addLink(BPLink *link, double inpMsgs[2])
  {
    links_->append(link);
    double *msgs;
    msgs = new double[2];
			   
    if (inpMsgs)
    {
      msgs[0] = inpMsgs[0];
      msgs[1] = inpMsgs[1];
    }
    else
    {
      msgs[0] = msgs[1] = 0;
    }
    msgsArr_->append(msgs);
    msgs = new double[2];
    nextMsgsArr_->append(msgs);
  }

  /**
   * Stores the message sent over a particular link.
   */
  void receiveMessage(double* inpMsgs, BPLink *link)
  {
    double *nextMsgs;
    int reverseNodeIndex = link->getReverseNodeIndex();
    nextMsgs = (*nextMsgsArr_)[reverseNodeIndex];
    nextMsgs[0] = inpMsgs[0];
    nextMsgs[1] = inpMsgs[1];
  }
		  
    //find the outgoing message for the given predIndex - helper for sendMessage
  double* multiplyMessagesAndSumOut(int predIndex);

    //send Message on all the links
  void sendMessage();
         
    //update the stored msgs and update the msgProduct
  void moveToNextStep();

  ostream& print(ostream& out);

  ostream& printWts(ostream& out);

 protected:
  Clause * clause_;

    //only one of the two below is used - superClause in case of lifted
    //inference and constants in case of ground inference
  SuperClause * superClause_;
  Array<int> * constants_;
  Domain * domain_;
  Array<BPLink *> *links_;
    // Incoming messages from nodes
  Array<double *> *msgsArr_;
    // Incoming messages from nodes for next step
  Array<double *> *nextMsgsArr_;
    // Outgoing messages
  double *factorMsgs_;
  double outputWt_;
};

#endif

