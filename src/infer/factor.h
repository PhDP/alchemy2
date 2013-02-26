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
#ifndef _FACTOR_H_Jan_2008
#define _FACTOR_H_Jan_2008

#include <math.h>
#include "util.h"
#include "mrf.h"
#include "array.h"
#include "link.h"
#include "node.h"
#include "superclause.h"

const double ACTION_WT = 20.0;

using namespace std;
using namespace __gnu_cxx;

/**
 * Stores all relevant information about a factor in a factor graph.
 */
class Factor
{
 public:

  /**
   * Constructor.
   */
  Factor(Clause * const & clause, SuperClause * const & superClause, 
         Array<int> * const & constants, Domain * const & domain,
         double outputWt)
  {
    clause_ = clause;
    superClause_ = superClause;
    constants_ = constants;
    domain_ = domain;

    msgsArr_ = new Array<double *>;
    nextMsgsArr_ = new Array<double *>;
    links_ = new Array<Link *>();
    outputWt_ = outputWt;
    initFactorMesssages();
    
    action_ = false;
    actionValue_ = true;

    taggedSend_ = false;
    taggedReceive_ = false;
  }
  
  /**
   * Constructor.
   */
  Factor()
  {
    links_ = new Array<Link *>();
    msgsArr_ = new Array<double *>;

    action_ = false;
    actionValue_ = true;

    taggedSend_ = false;
    taggedReceive_ = false;
  }
		  
  /**
   * Destructor.
   */
  ~Factor()
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

  /**
   * Initializes outgoing factor messages. This is the contribution of the
   * factor itself.
   */
  void initFactorMesssages();
  
  /**
   * Returns the id of the super clause associated with this factor.
   */
  int getSuperClauseId()
  {
    if (superClause_) 
      return superClause_->getSuperClauseId();
    else
      return -1;
  }
		  
  /**
   * Returns the id of the parent super clause associated with this factor.
   */
  int getParentSuperClauseId()
  {
    if( superClause_) 
      return superClause_->getParentSuperClauseId();
    else
      return -1;
  }

  /**
   * Returns the super clause associated with this factor.
   */
  SuperClause* getSuperClause() {return superClause_;}

  /**
   * Returns the clause associated with this factor.
   */
  Clause * getClause() {return clause_;}

  /**
   * Returns the domain associated with this factor.
   */
  Domain *getDomain() {return domain_;}

  /**
   * Returns the constants associated with this factor.
   */
  Array<int> * getConstants() {return constants_;}

  /**
   * Returns the number of links of (or number of nodes attached to) this
   * factor.
   */
  int getNumLinks() {return links_->size();}
		  
  /**
   * Returns the number of links of (or number of nodes attached to) this
   * factor.
   */
  Link* getLink(const int& index) {return (*links_)[index];}

  /**
   * Gets the incoming message from a certain link.
   * 
   * @param index Index of the link
   * @param msgs msgs[0] contains the message from the node when it is false,
   * msgs[1] contains the message from the node when it is true.
   */
  void getMessage(int index, double msgs[])
  { 
    msgs[0] = (*msgsArr_)[index][0];
    msgs[1] = (*msgsArr_)[index][1];
  }

  /**
   * Adds a link to this factor.
   * 
   * @param link Link to be added.
   * @inpMsgs Input message from the node in the link. inpMsgs[0] contains the
   * message from the node when it is false, inpMsgs[1] contains the message
   * from the node when it is true. If not specified, the messages are set to 0.
   */
  void addLink(Link *link, double inpMsgs[2])
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
  void receiveMessage(double* inpMsgs, Link *link)
  {
    double *nextMsgs;
    int reverseNodeIndex = link->getReverseNodeIndex();
    nextMsgs = (*nextMsgsArr_)[reverseNodeIndex];
    nextMsgs[0] = inpMsgs[0];
    nextMsgs[1] = inpMsgs[1];
  }

  /**
   * Finds the outgoing message for the given predIndex - helper for sendMessage
   */
  double* multiplyMessagesAndSumOut(int predIndex);

  /**
   * Sends the messages to all the links.
   */
  void sendMessage();
   
  /**
   * Updates the stored msgs and updates the msgProduct
   */
  void moveToNextStep();

  ostream& print(ostream& out);

  bool isAction() { return action_; }

  bool getActionValue() { return actionValue_; }

  void flip()
  {
    actionValue_ = !actionValue_;
    factorMsgs_[actionValue_] = clause_->getWt();
    factorMsgs_[!actionValue_] = 0.0;
  }

    // Only makes sense for action factors!
  Node* getActionNode()
  {
    if (!action_)
    {
      cout << "NOT AN ACTION FACTOR." << endl;
      exit(1);
    }
    return (*links_)[0]->getNode();
  }

  void setAction()
  {
    //cout << "SETTING ACTION" << endl;
    action_ = true;
    actionValue_ = true;
    clause_->setWt(ACTION_WT);
    factorMsgs_[actionValue_] = clause_->getWt();
    factorMsgs_[!actionValue_] = 0.0;
  }

    // Add neighbors to send and receive sets.
  void tagNeighborsSendReceive(set<Factor*> &sendFactors,
                               set<Factor*> &receiveFactors,
                               set<Node*> &sendNodes,
                               set<Node*> &receiveNodes);
  void tagNeighborsSendReceive(list<Factor*> &sendFactors,
                               list<Factor*> &receiveFactors,
                               list<Node*> &sendNodes,
                               list<Node*> &receiveNodes);

    // Add neighbors to receive set.
  void tagNeighborsReceive(set<Node*> &receiveNodes);
  void tagNeighborsReceive(list<Node*> &receiveNodes);

  void resetTags()
  {
    taggedSend_ = false;
    taggedReceive_ = false;
  }
   
  static void setIsApproxNetwork(bool val);

  ostream& printWts(ostream& out);

 protected:
  bool action_;
  bool actionValue_;

  Clause * clause_;

    //only one of the two below is used - superClause in case of lifted
    //inference and constants in case of ground inference
  SuperClause * superClause_;
  Array<int> * constants_;
  Domain * domain_;
  Array<Link *> *links_;
    // Incoming messages from nodes
  Array<double *> *msgsArr_;
    // Incoming messages from nodes for next step
  Array<double *> *nextMsgsArr_;
    // Outgoing messages
  double *factorMsgs_;
  double outputWt_;
  static bool isApproxNetwork__;
  
  bool taggedSend_, taggedReceive_;
};

#endif

