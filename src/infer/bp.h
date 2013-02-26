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
#ifndef BP_H_
#define BP_H_

#include "inference.h"
#include "bpparams.h"
#include "twowaymessage.h"
#include "superclause.h"
#include "auxfactor.h"
#include "node.h"
#include "factorgraph.h"

#include <set>
#include <vector>
using namespace std;

const int bpdebug = false;

//const bool EFBP = true;
const bool VERBOSE = false;

const int MAX_ROUNDS = 100;
const int MIN_CONV_ITR = 10;
const double EVAL_THRESH = 1e-4;
const double SEARCH_THRESH = 1e-3;

//const double MAX_TIME = 259200.0; // 3 days
const double MAX_TIME = 288000.0; // 80 hours

/**
 * Class for belief propagation algorithm. This version of BP works on a factor
 * graph, which can be in a lifted representation or not.
 */
class BP : public Inference
{
 public:

  /**
   * Constructor. Requires a factor graph and a set of parameters for the
   * algorithm. Optionally, a set of query formulas is used.
   */
  BP(FactorGraph* factorGraph, BPParams* bpParams, const bool& efbp,
     Array<Array<Predicate* >* >* queryFormulas = NULL)
    : Inference(NULL, -1, false, queryFormulas)
  {
    factorGraph_ = factorGraph;
    maxSteps_ = bpParams->maxSteps;
    maxSeconds_ = bpParams->maxSeconds;
    convergenceThresh_ = bpParams->convergenceThresh;
    convergeRequiredItrCnt_ = bpParams->convergeRequiredItrCnt;
    outputNetwork_ = bpParams->outputNetwork;
    efbp_ = efbp;
  }

  /**
   * Destructor.
   */
  ~BP()
  {
  }
 
  /**
   * Initializes belief propagation. The factor graph is built.
   */
  void init()
  {
    Timer timer1;
    cout << "Initializing ";
    cout << "Belief Propagation..." << endl;

    factorGraph_->init();
    if (bpdebug)
    {
      cout << "[init] ";
      Timer::printTime(cout, timer1.time());
      cout << endl;
      timer1.reset();
    }
  }

  bool runBP(vector<Node*> &nodes, vector<Factor*> &factors, double conv_thresh)
  {
    if (VERBOSE)
      cout << "BP" << endl;

    int convItr = 0;
    for(int rd = 0; rd < MAX_ROUNDS; rd++)
    {
      if (VERBOSE)
        cout << "Round " << rd+1 << " of " << MAX_ROUNDS << endl;

      for (unsigned int i = 0; i < factors.size(); i++)
        factors[i]->sendMessage();
      for (unsigned int i = 0; i < nodes.size(); i++)
        nodes[i]->sendMessage();

      for (unsigned int i = 0; i < factors.size(); i++)
        factors[i]->moveToNextStep();
      bool conv = true;
      double maxDiff = 0.0;
      int numChanges = 0;
      for (unsigned int i = 0; i < nodes.size(); i++)
      {
        double oldProb = nodes[i]->getProb();
        nodes[i]->moveToNextStep();
        double newProb = nodes[i]->getProb();
        if(abs(oldProb - newProb) > conv_thresh)
        {
          if(VERBOSE)
          {
            cout << "\tNode " << i << " [";
            nodes[i]->print(cout);
            cout << "] changed from " << oldProb << " to " << newProb << endl;
          }
          conv = false;
          numChanges++;
        }
        if(abs(oldProb - newProb) > maxDiff)
          maxDiff = abs(oldProb - newProb);
        //if(VERBOSE)
        //  cout << nodes[i]->getId() << ": " << newProb << endl;
      }
      if(VERBOSE)
        cout << numChanges << " nodes changed." << endl;
      if(VERBOSE)
        cout << "Max diff: " << maxDiff << endl << endl;

      if(conv)
        convItr++;
      else
        convItr = 0;

      if(convItr >= MIN_CONV_ITR)
      {
        if(VERBOSE)
        {
          cout << "Converged in " << rd+1-MIN_CONV_ITR << " + "
            << MIN_CONV_ITR << " = " << rd+1 << " rounds." << endl;
        }
        break;
      }
    }
    if(convItr < MIN_CONV_ITR)
      return false;
    return true;
  }


  bool runEFBP(set<Node*> &sendNodes,
      set<Factor*> &sendFactors,
      set<Node*> &receiveNodes,
      set<Factor*> &receiveFactors,
      double conv_thresh)
  {
    if (VERBOSE)
      cout << "EFBP" << endl;

    set<Node*>::iterator nIter;
    set<Factor*>::iterator fIter;

    int convItr = 0;
    for(int rd = 0; rd < MAX_ROUNDS; rd++)
    {
      if(VERBOSE)
      {
        cout << "Round " << rd+1 << " of " << MAX_ROUNDS << endl;
        cout << "sendNodes: " << sendNodes.size() << endl;
        cout << "sendFactors: " << sendFactors.size() << endl;
        cout << "receiveNodes: " << receiveNodes.size() << endl;
        cout << "receiveFactors: " << receiveFactors.size() << endl;
      }

      //set<Node*> nextSendNodes = sendNodes;
      //set<Factor*> nextSendFactors = sendFactors;
      //set<Node*> nextReceiveNodes = receiveNodes;
      //set<Factor*> nextReceiveFactors = receiveFactors;
      list<Node*> nextSendNodes;
      list<Factor*> nextSendFactors;
      list<Node*> nextReceiveNodes;
      list<Factor*> nextReceiveFactors;

      for(fIter = sendFactors.begin(); fIter != sendFactors.end(); fIter++)
        (*fIter)->sendMessage();
      for(nIter = sendNodes.begin(); nIter != sendNodes.end(); nIter++)
        (*nIter)->sendMessage();

      for(fIter = receiveFactors.begin(); fIter != receiveFactors.end(); fIter++)
        (*fIter)->moveToNextStep();

      bool conv = true;
      double maxDiff = 0.0;
      int numChanges = 0;
      for(nIter = receiveNodes.begin(); nIter != receiveNodes.end(); nIter++)
      {
        double oldProb = (*nIter)->getProb();
        (*nIter)->moveToNextStep();
        double newProb = (*nIter)->getProb();
          
        if(abs(oldProb - newProb) > conv_thresh)
        {
          //cout << '\t' << receiveNodes[i]->getId() << " changed from " << oldProb << " to " << newProb << endl;
          //cout << "Tagging." << endl;
          numChanges++;
          conv = false;
          (*nIter)->tagNeighborsReceive(nextSendFactors, nextReceiveFactors, nextSendNodes, nextReceiveNodes);
          //cout << "Tagged." << endl;
        }
        if(abs(oldProb - newProb) > maxDiff)
          maxDiff = abs(oldProb - newProb);

        //if(VERBOSE)
        //  cout << (*nIter)->getId() << ": " << newProb << endl;
      }

      sendNodes.insert(nextSendNodes.begin(), nextSendNodes.end());
      sendFactors.insert(nextSendFactors.begin(), nextSendFactors.end());
      receiveNodes.insert(nextReceiveNodes.begin(), nextReceiveNodes.end());
      receiveFactors.insert(nextReceiveFactors.begin(), nextReceiveFactors.end());

      //if(VERBOSE)
      //  cout << "Max diff: " << maxDiff << endl << endl;

      if(conv)
        convItr++;
      else
        convItr = 0;

      if(convItr >= MIN_CONV_ITR)
      {
        if(VERBOSE)
        {
          cout << "Converged in " << rd+1-MIN_CONV_ITR << " + "
            << MIN_CONV_ITR << " = " << rd+1 << " rounds." << endl;
        }
        break;
      }
    }
    //for(nIter = receiveNodes.begin(); nIter != receiveNodes.end(); nIter++)
    //  (*nIter)->resetTags();
    //for(fIter = receiveFactors.begin(); fIter != receiveFactors.end(); fIter++)
    //  (*fIter)->resetTags();
    if(convItr < MIN_CONV_ITR)
      return false;
    return true;
  }


  void runDecisionBP()
  {
    Timer timer;
    vector<Node*> nodes;
    vector<Factor*> factors;

    for (int i = 0; i < factorGraph_->getNumNodes(); i++)
      nodes.push_back(factorGraph_->getNode(i));
    for (int i = 0; i < factorGraph_->getNumFactors(); i++)
      factors.push_back(factorGraph_->getFactor(i));

    for (unsigned int i = 0; i < factors.size(); i++)
    {
      if (factors[i]->isAction())
      {
        //cout << "ACTION" << '\t';
        //factors[i]->print(cout);
        //cout << endl;
        factors[i]->flip();
        actionFactors_.push_back(factors[i]);
      }
    }

    if (actionFactors_.size() == 0)
    {
      cout << "ERROR: No action predicates defined!" << endl;
      exit(1);
    }

    cout << "EFBP: " << efbp_ << endl;
    cout << "MIN_CONV_ITR: " << MIN_CONV_ITR << endl;
    cout << "EVAL_THRESH: " << EVAL_THRESH << endl;
    cout << "SEARCH_THRESH: " << SEARCH_THRESH << endl;
    cout << endl;

    cout << nodes.size() << " nodes." << endl;
    cout << factors.size() << " factors." << endl;
    cout << actionFactors_.size() << " actions." << endl;

    bestActionState_.resize(actionFactors_.size());

    //cout << "Pretagging neighbors." << endl;
    //for(unsigned int i = 0; i < nodes.size(); i++)
    //  nodes[i]->preTagNeighborsReceive();

    cout << "Making decisions." << endl;
    Timer::printTime(cout, timer.time());
    cout << endl;

    double util, prevUtil = -DBL_MAX;
    double maxUtil = -DBL_MAX;
    unsigned int maxUtilItr = (unsigned int) -1;
    bool reverse = false;
    double prevTime = timer.time();
    double initTime = timer.time();
    Factor* flipFactor = NULL;
    Factor* prevFactor = NULL;

    //int numFlips = 2 * actionFactors_.size() + 1;
    int numFlips = 100 * actionFactors_.size() + 1;

    for(int i = 0; i < numFlips; i++)
    {
      if (i%100 == 0 || VERBOSE)
        cout << i << "/" << numFlips << endl;

      set<Node*> sendNodes;
      set<Factor*> sendFactors;
      set<Node*> receiveNodes;
      set<Factor*> receiveFactors;
      set<Node*>::iterator nIter;
      set<Factor*>::iterator fIter;

      bool converge = true;
      if (efbp_ && i > 0)
      {
        sendNodes.insert(flipFactor->getActionNode());
        receiveNodes.insert(flipFactor->getActionNode());
        flipFactor->getActionNode()->tagNeighborsSend(sendFactors,
                                                      receiveFactors,
                                                      receiveNodes);
        if (reverse)
        {
          sendNodes.insert(prevFactor->getActionNode());
          receiveNodes.insert(prevFactor->getActionNode());
          prevFactor->getActionNode()->tagNeighborsSend(sendFactors,
                                                        receiveFactors,
                                                        receiveNodes);
        }
        converge = runEFBP(sendNodes, sendFactors, receiveNodes, receiveFactors,
                           SEARCH_THRESH);
        if (!converge)
          cout << "FAILED TO CONVERGE in iteration " << i << endl;

        util = prevUtil;
        for(nIter = receiveNodes.begin(); nIter != receiveNodes.end(); nIter++)
        {
          if(isnan((*nIter)->getProb()) || VERBOSE)
          {
            cout << "Node ";
            (*nIter)->print(cout);
            cout << ": prob = " << (*nIter)->getProb() << endl;
          }
          util += (*nIter)->getUtil() - (*nIter)->getPrevUtil();
        }
        // TODO: restore complex util factors.
        //for(fIter = receiveFactors.begin(); fIter != receiveFactors.end(); fIter++)
        //{
        //  util += (*fIter)->getUtil() - (*fIter)->getPrevUtil();
        //}
        //cout << receiveNodes.size() << " nodes updated." << endl;
      }
      else
      {
        if (i == 0)
          converge = runBP(nodes, factors, EVAL_THRESH);
        else
          converge = runBP(nodes, factors, SEARCH_THRESH);
        if (!converge)
          cout << "FAILED TO CONVERGE in iteration " << i << endl;

        if(VERBOSE)
          cout << "Done with BP" << endl;
        util = 0.0;
        for (unsigned int j = 0; j < nodes.size(); j++)
        {
          if (isnan(nodes[j]->getProb()) || VERBOSE)
          {
            cout << "Node " << j << " [";
            nodes[j]->print(cout);
            cout << "] : prob = " << nodes[j]->getProb() 
                 << ", util = " << nodes[j]->getUtil() << endl;
          }
          //nodes[j]->print();
          util += nodes[j]->getUtil();
        }
        //for(unsigned int j = 0; j < factors.size(); j++)
        //{
        //  util += factors[j]->getUtil();
        //}
      }
      if (!converge)
        util = -LDBL_MAX;

      if (VERBOSE)
        cout << "Util = " << util << endl;
      if (i%100 == 0)
      {
        double curTime = timer.time();
        Timer::printTime(cout, curTime);
        cout << endl;
        if (i > 0)
        {
          cout << "Total:  " << i/(curTime - initTime) << " flips/sec."
               << endl;
          cout << "Recent: " << 100/(curTime - prevTime) << " flips/sec."
               << endl;
        }
        else
          initTime = curTime;
        prevTime = curTime;
      }

      reverse = false;
      if (util < prevUtil)
      {
        prevFactor = actionFactors_[(i-1)%actionFactors_.size()];
        if(VERBOSE)
        {
          cout << "Reversing " << (i-1)%actionFactors_.size() << "; ";
          prevFactor->print(cout);
          cout << endl;
        }
        prevFactor->flip();
        reverse = true;
      }
      else
      {
        prevUtil = util;

        if (efbp_)
        {
          if (VERBOSE)
            cout << "Updating probabilities." << endl;
          if (i == 0)
          {
            for (unsigned int j = 0; j < nodes.size(); j++)
            {
              nodes[j]->updateProb();
            }
            //for(unsigned int j = 0; j < factors.size(); j++)
            //  factors[j]->updateProb();
          }
          else
          {
            for (nIter = receiveNodes.begin(); nIter != receiveNodes.end();
                 nIter++)
              (*nIter)->updateProb();
            //for(fIter = receiveFactors.begin(); fIter != receiveFactors.end(); fIter++)
            //  (*fIter)->updateProb();
          }
        }
      }

      if (util > maxUtil)
      {
        maxUtil = util;
        maxUtilItr = i;
        //cout << "New max util in round " << i << ": " << util << endl;

        for (unsigned int i = 0; i < actionFactors_.size(); i++)
          bestActionState_[i] = actionFactors_[i]->getActionValue();
      }

      flipFactor = actionFactors_[i%actionFactors_.size()];
      if (VERBOSE)
      {
        cout << "Flipping " << i%actionFactors_.size() << "; ";
        flipFactor->print(cout);
        cout << "[";
        flipFactor->getActionNode()->print(cout);
        cout << "]" << endl;
      }
      //actionFactors_[i%actionFactors_.size()]->flip();
      flipFactor->flip();

      if (i%actionFactors_.size() == 0)
      {
        cout << "Run " << i/actionFactors_.size() << ": best = " << maxUtil 
             << endl;
      }

      if (i - maxUtilItr > actionFactors_.size())
      {
        cout << "Action search has converged in iteration " << i << "." << endl;
        break;
      }
      if (timer.time() > MAX_TIME)
      {
        cout << "OK, I'm bored now. Stopping in iteration " << i << "." << endl;
        break;
      }
    }
    cout << "Max util: " << maxUtil << " in iteration " << maxUtilItr
         << endl << endl;

    cout << "Restoring optimal state:" << endl;
    for (unsigned int i = 0; i < actionFactors_.size(); i++)
    {
      if (bestActionState_[i] != actionFactors_[i]->getActionValue())
        actionFactors_[i]->flip();
      actionFactors_[i]->print(cout);
      cout << ": " << actionFactors_[i]->getActionValue() << endl;
    }
    runBP(nodes, factors, EVAL_THRESH);
    util = 0.0;
    for (unsigned int j = 0; j < nodes.size(); j++)
    {
      //if(nodes[j]->isUtil())
      //  nodes[j]->print(cout);
      util += nodes[j]->getUtil();
    }
    //for(unsigned int j = 0; j < factors.size(); j++)
    //{
    //  util += factors[j]->getUtil();
    //}
    cout << "Util = " << util << endl;
    maxUtil_ = util;

    //for(unsigned int i = 0; i < nodes.size(); i++)
    //  delete nodes[i];
    //for(unsigned int i = 0; i < factors.size(); i++)
    //  delete factors[i];
  }


  /**
   * Performs Belief Propagation inference.
   */
  void infer()
  {
//    vector<Node*> nodes;
//    vector<Factor*> factors;
//    vector<Factor*> actionFactors;
//
//    for (int i = 0; i < factorGraph_->getNumNodes(); i++)
//      nodes.push_back(factorGraph_->getNode(i));
//    for (int i = 0; i < factorGraph_->getNumFactors(); i++)
//      factors.push_back(factorGraph_->getFactor(i));
//
//    runBP(nodes, factors, EVAL_THRESH);
//    return;

    Timer timer1;

    double oldProbs[2];
    double newProbs[2];
    double diff;
    double maxDiff;
    int maxDiffNodeIndex;
    int convergeItrCnt = 0;
    bool converged = false;
    int numFactors = factorGraph_->getNumFactors();
    int numNodes = factorGraph_->getNumNodes();

    cout << "factorcnt = " << numFactors
         << ", nodecnt = " << numNodes << endl;

    if (bpdebug)
    {
      cout << "factors:" << endl;
      for (int i = 0; i < numFactors; i++)
      {
        factorGraph_->getFactor(i)->print(cout);
        cout << endl;
      }
      cout << "nodes:" << endl;
      for (int i = 0; i < numNodes; i++)
      {
        factorGraph_->getNode(i)->print(cout);
        cout << endl;
      }
    }
      // Pass around (send) the messages
    int itr;
    
      // Move to next step to transfer the message in the nextMsgsArr in the
      // beginning to the the msgsArr 
    for (itr = 1; itr <= maxSteps_; itr++)
    {
      if (bpdebug)
      {
        cout<<"*************************************"<<endl;
        cout<<"Performing Iteration "<<itr<<" of BP"<<endl;
        cout<<"*************************************"<<endl;
      }

      for (int i = 0; i < numFactors; i++)
      {
        if (bpdebug)
        {
          cout << "Sending messages for Factor: ";
          factorGraph_->getFactor(i)->print(cout); cout << endl;
        }
        factorGraph_->getFactor(i)->sendMessage();
      }
           
      for (int i = 0; i < numNodes; i++)
      {
        if (bpdebug)
        {
          cout << "Sending messages for Node: ";
          factorGraph_->getNode(i)->print(cout); cout << endl;
        }
        factorGraph_->getNode(i)->sendMessage();
      }

      for (int i = 0; i < numFactors; i++)
      {
        factorGraph_->getFactor(i)->moveToNextStep();
        if (bpdebug)
        {
          cout << "BP-Factor Iteration " << itr << " => ";
          factorGraph_->getFactor(i)->print(cout); cout << endl;
        }
      }
          
      maxDiff = -1;
      maxDiffNodeIndex = -1;
      for (int i = 0; i < numNodes; i++)
      {
        if (bpdebug)
        {
          cout<<"************************************"<<endl;
          cout<<"Node "<<i<<":"<<endl;
          cout<<"************************************"<<endl;
          cout<<"Getting Old Probabilities =>"<<endl; 
          cout<<endl;
          cout<<"Moving to next step "<<endl;
          cout<<endl;
          cout<<"Getting New Probabilities =>"<<endl; 
        }

        factorGraph_->getNode(i)->getProbs(oldProbs);
        factorGraph_->getNode(i)->moveToNextStep();
        factorGraph_->getNode(i)->getProbs(newProbs);

        diff = abs(newProbs[1] - oldProbs[1]);

        if (bpdebug)
        {
          cout<<endl<<endl<<"Final Probs : "<<endl;
          cout<<"Node "<<i<<": probs["<<0<<"] = "<<newProbs[0]
              <<", probs["<<1<<"] = "<<newProbs[1]<<endl;
          cout<<"BP-Node Iteration "<<itr<<": "<<newProbs[0]
              <<"  probs["<<1<<"] = "<<newProbs[1]<<endl;
          cout<<" : => ";
          factorGraph_->getNode(i)->print(cout);
          cout << endl;
        }
        
        if (maxDiff < diff)
        {
          maxDiff = diff;
          maxDiffNodeIndex = i;
        }
      }
           
      cout<<"At Iteration "<<itr<<": MaxDiff = "<<maxDiff<<endl;
      cout<<endl;
           
        //check if BP has converged
      if (maxDiff < convergenceThresh_)
        convergeItrCnt++;
      else
        convergeItrCnt = 0;

        // Check if for N continuous iterations, maxDiff has been below the
        // threshold
      if (convergeItrCnt >= convergeRequiredItrCnt_)
      {
        converged = true;
        //run all the way till end
		//break;    
      }
    }
      
    if (converged)
    {
      cout << "Converged in " << itr << " Iterations " << endl;
    }
    else
    {
      cout << "Did not converge in " << maxSteps_ << " (max allowed) Iterations"
           << endl;
    }
    
    if (queryFormulas_)
    {
      cout << "Computing probabilities of query formulas ..." << endl;
      for (int i = 0; i < numNodes; i++)
      {
        if (bpdebug)
        {
          cout << "Sending auxiliary messages for Node: ";
          factorGraph_->getNode(i)->print(cout); cout << endl;
        }
        factorGraph_->getNode(i)->sendAuxMessage();
          // Now, messages have been sent to the aux. factors
      }
      for (int j = 0; j < qfProbs_->size(); j++)
      {
        (*qfProbs_)[j] = factorGraph_->getAuxFactor(j)->getProb();
      }
    }
  }

  /**
   * Prints out the network.
   */
  void printNetwork(ostream& out)
  {
    factorGraph_->printNetwork(out);
  }
  
  /**
   * Prints the probabilities of each predicate to a stream.
   */
  void printProbabilities(ostream& out)
  {
    IntArrayHashArray *neqConstraints;
	Predicate *varPred;
	double probs[2];
    HyperCube *hyperCube;
	Array<Array<int>*> *constantsArr;
	Array<int>* constants;
    Predicate* pred;
    int predId;
    Node* node;
    double exp;
    Domain* domain = factorGraph_->getDomain();
    Database* db = domain->getDB();

	int predCnt = domain->getNumPredicates();

	//nodes stored by the predId
    Array<Array<Node *> *> *nodesByPredIdArr = 
      new Array<Array<Node *>*>(predCnt, NULL);
    Array<Node *>* nodesByPredId;
    PredicateHashArray *predGndingsHash;

	for (int i = 0; i < factorGraph_->getNumNodes(); i++)
    { 
      node = factorGraph_->getNode(i);
      predId = node->getPredId();
      if((*nodesByPredIdArr)[predId] == NULL) {
		   (*nodesByPredIdArr)[predId] = new Array<Node *>();
	  }
	  (*nodesByPredIdArr)[predId]->append(node);
	}
	   
    
	//while(true){}
    for(predId=0;predId<predCnt;predId++) {
	   nodesByPredId = (*nodesByPredIdArr)[predId];
	   if(nodesByPredId == NULL)
			  continue;
	   varPred = domain->createPredicate(predId,false);
	   cout<<"Processing the predicate ";
	   varPred->printWithStrVar(cout,domain);
	   varPred->print(cout,domain);
	   cout<<endl;
	   neqConstraints = factorGraph_->getNeqConstraints(varPred);
	   
	   int numTrue = db->getNumTrueGndPreds(predId);
	   int numFalse = db->getNumFalseGndPreds(predId);	 
	   
	   //If this is not a completely unknown predicate, then need
	   //to find out the unknown groundings explicitly and print
	   //their probabilities
	   //if(true) 
	   if(numTrue > 0 || numFalse > 0) {
        predGndingsHash = new PredicateHashArray();
	    factorGraph_->getUnknownGndingsHash(predGndingsHash,varPred,domain,db);
	   } else {
		predGndingsHash = NULL;
	   }
	   
	   //cout<<"Got "<<predGndingsHash->size()<<" gndings for the pred ";
	   
	   //while(true){}
	   for (int i = 0; i < nodesByPredId->size(); i++)
       { 
        node = (*nodesByPredId)[i];
	    node->getProbs(probs);         
        exp = node->getExp();
        SuperPred * superPred = node->getSuperPred();

        if (superPred)
        {
		 for (int index = 0; index < superPred->getNumHyperCubes(); index++)
         {
		  hyperCube = superPred->getHyperCube(index);
		  
		  if(factorGraph_->getUseCT()) {
		   constantsArr = hyperCube->getTuples(neqConstraints);
		  } else {
		   constantsArr = hyperCube->getTuples(NULL);
		  }

		  for(int tno=0;tno<constantsArr->size();tno++) {
			constants = (*constantsArr)[tno];
            //remove the constant for the dummy variable
			constants->removeItem(0);
			pred = domain->getPredicate(constants, predId);
			delete constants;
            
			if(predGndingsHash) {
			 //print the probabilities if this predicate actually belongs to the unknown list
			 int index = predGndingsHash->find(pred);
			 delete pred;
			 if(index < 0) 
				continue;
			 pred = predGndingsHash->removeItemFastDisorder(index);
			}
			pred->printWithStrVar(out, domain);
            out << " " << probs[1] << endl;
            //out<<" "<<exp<<endl;
			delete pred;
		   }
		  delete constantsArr;
		}
      }
	  else
      {
        constants = node->getConstants(); 
        assert(constants != NULL);
        pred = domain->getPredicate(constants, predId);
		if(predGndingsHash) {
	     int index = predGndingsHash->find(pred);
		 delete pred;
		 if(index < 0)
		  continue;
		 pred = predGndingsHash->removeItemFastDisorder(index);
		}
	    pred->printWithStrVar(out, domain);
        out << " " << probs[1] << endl;
        //out<<" "<<exp<<endl;
	    delete pred;
	   }
	}
	
	if(predGndingsHash) {
	 //if no noise was introduced, all unknown preds should be output by now
	 if (factorGraph_->getUseHC() && factorGraph_->getHcCreateNoise() == 0) {
		 assert(predGndingsHash->size() == 0);   
	 }
	 //print out any remaining predicates
     for(int i=0;i<predGndingsHash->size();i++) {
		 pred = (*predGndingsHash)[i];
		 pred->printWithStrVar(out, domain);
		 out<<" "<<0<<endl;
		 delete pred;
     }
	 delete predGndingsHash;
	}

	delete neqConstraints;
    delete varPred;
   }
 
   //clear up
   nodesByPredIdArr->deleteItemsAndClear();
   delete nodesByPredIdArr;
 }

  /**
   * Puts the predicates whose probability has changed with respect to the
   * reference vector oldProbs by more than probDelta in string form and the
   * corresponding probabilities of each predicate in two vectors. Currently
   * not implemented.
   * 
   * @param changedPreds Predicates whose probability have changed more than
   * probDelta are put here.
   * @param probs The probabilities corresponding to the predicates in
   * changedPreds are put here.
   * @param oldProbs Reference probabilities for checking for changes.
   * @param probDelta If probability of an atom has changed more than this
   * value, then it is considered to have changed.
   */
  void getChangedPreds(vector<string>& changedPreds, vector<float>& probs,
                       vector<float>& oldProbs, const float& probDelta)
  {
  }

  /**
   * Gets the probability of a ground predicate.
   * 
   * @param gndPred GroundPredicate whose probability is being retrieved.
   * @return Probability of gndPred if present in state, otherwise 0.
   */
  double getProbability(GroundPredicate* const& gndPred)
  {
/*
    double probs[2];
    Array<int>* constants;
    Predicate* pred;
    unsigned int predId;
    Node* node;
    Domain* domain = factorGraph_->getDomain();
    bool found = false;
    for (int i = 0; i < factorGraph_->getNumNodes(); i++)
    { 
      node = factorGraph_->getNode(i);
      predId = node->getPredId();
      if (predId != gndPred->getId()) continue;
      node->getProbs(probs);         
      SuperPred * superPred = node->getSuperPred();

      if (superPred)
      {        
        for (int index = 0; index < superPred->getNumTuples(); index++)
        {
          constants = superPred->getConstantTuple(index);
          pred = domain->getPredicate(constants, predId);
          if (!pred->same(gndPred))
          {
            delete pred;
            continue;
          }
          delete pred;
          found = true;
          return probs[1];
        }
      }
      else
      {
        constants = node->getConstants(); 
        assert(constants != NULL);
        pred = domain->getPredicate(constants, predId);
        if (!pred->same(gndPred))
        {
          delete pred;
          continue;
        }
        delete pred;
        found = true;
        return probs[1];
      }
    }
*/
    return 0.5;
  }

  /**
   * Gets the probability of a ground predicate. Currently not implemented.
   * 
   * @param gndPred GroundPredicate whose probability is being retrieved.
   * @return Probability of gndPred if present in state, otherwise 0.
   */
  double getProbabilityH(GroundPredicate* const& gndPred)
  {
    return 0.0;
  }

  /**
   * Prints each predicate with a probability of 0.5 or greater to a stream.
   * Currently not implemented.
   */
  void printTruePreds(ostream& out)
  {
  }
  
  /**
   * Prints each predicate with a probability of 0.5 or greater to a stream.
   * Currently not implemented.
   */
  void printTruePredsH(ostream& out)
  {
  }
  
  FactorGraph* getFactorGraph()
  {
    return factorGraph_;
  }
 
  void printDecisionResults(ostream& out)
  {
    Domain* domain = factorGraph_->getDomain();
    out << "Max. Utility: " << maxUtil_ << endl;
    for (unsigned int i = 0; i < actionFactors_.size(); i++)
    {
      Factor* factor = actionFactors_[i];
      Node* node = factor->getLink(0)->getNode();
      int predId = node->getPredId();
      Array<int>* constants = node->getConstants(); 
      assert(constants != NULL);
      Predicate* pred = domain->getPredicate(constants, predId);
      pred->printWithStrVar(out, domain);
      out << " " << factor->getActionValue() << endl;
      delete pred;
    }
  }
 
 private:
    // Network on which BP is run
  FactorGraph* factorGraph_;
    // Max. no. of BP iterations to perform
  int maxSteps_;
    // Max. no. of seconds BP should run
  int maxSeconds_;
  
    // Maximum difference between probabilities must be less than this
    // in order to converge
  double convergenceThresh_;
    // Convergence must last this number of iterations
  int convergeRequiredItrCnt_;
    // No inference is run, rather the factor graph is built
  bool outputNetwork_;
  
  bool efbp_;
  
  double maxUtil_;

  vector<Factor*> actionFactors_;
  vector<bool> bestActionState_;
};

#endif /*BP_H_*/
