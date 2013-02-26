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
#ifndef FACTORGRAPH_H_
#define FACTORGRAPH_H_

#include "clausefactory.h"
#include "twowaymessage.h"
#include "superclause.h"
#include "auxfactor.h"
#include "node.h"
#include "hypercube.h"
#include "hypercubeoperations.h"

const int fgdebug = true;

typedef hash_map<Predicate*, HyperCubeReverseIndex *, HashPredicate,
                 EqualPredicate> PredicateToHyperCubeReverseIndexMap;

typedef hash_map<Predicate*, Array<HyperCube*>*, HashPredicate, EqualPredicate> 
                PredicateToHyperCubeArrayMap;

typedef hash_map<int, IntArrayHashArray*, HashInt, EqualInt> 
                IntToIntArrayHashArrayMap;

extern IntHashArray* getPredicateIds(string & predsStr, Domain* const & domain);

/**
 * Class for a factor graph. It consists of an array of nodes and an array of
 * factors. The graph can be lifted or not, depending on if the lifted_ flag is
 * set or not. An optional array of auxiliary factors holds factors which are
 * attached to some nodes, but which do not send messages. These are used to
 * compute probabilities of query formulas (after inference is run, the
 * messages from the nodes attached to the auxiliary factors are to be sent.
 * Marginal probabilities are read directly from the nodes after inference is
 * run.
 */
class FactorGraph
{
 public:

  /**
   * Constructor. Data structures are initialized.
   */
  FactorGraph(bool lifted, bool useHC, bool useCT, bool implicitRep,
              HyperCubeCreateType hcCreateType, double hcCreateNoise,
              int lncIter, string noHCPredsStr, MLN* mln, Domain* domain,
              Array<Array<Predicate* >* >* queryFormulas = NULL)
  {
    lifted_ = lifted;
    useHC_ = useHC;
    useCT_ = useCT;
    implicitRep_ = implicitRep;
    hcCreateType_ = hcCreateType;
    hcCreateNoise_ = hcCreateNoise;
    lncIter_ = lncIter;
    noHCPredsStr_ = noHCPredsStr;
    //resultsFile_ = resultsFile;

    lidToTWMsg_ = new LinkIdToTwoWayMessageMap();
    superClausesArr_ = new Array<Array<SuperClause*>*>();
    factors_ = new Array<Factor*>();
    nodes_ = new Array<Node*>();
    mln_ = mln;
    domain_ = domain;

    auxFactors_ = NULL;
    if (queryFormulas)
    {
      auxFactors_ = new Array<AuxFactor*>();
      for (int i = 0; i < queryFormulas->size(); i++)
        auxFactors_->append(new AuxFactor((*queryFormulas)[i]));
    }
  }

  /**
   * Destructor. Data structures are destroyed.
   */
  ~FactorGraph()
  {
    delete lidToTWMsg_;
    delete superClausesArr_;
    delete factors_;
    delete nodes_;
    if (auxFactors_) delete auxFactors_;
  }
  
  /**
   * Builds the factor graph. If lifted_ is true, then the graph will be in a
   * lifted representation.
   */
  void init()
  {
    Timer timer1;

    removeTrivialClauses();

    cout << "Building ";
    if (lifted_) cout << "Lifted ";
    cout << "Factor Graph..." << endl;

    if (useHC_)
     cout << "Using Hyper Cubes " << endl;
    else
     cout << "NOT using Hyper Cubes " << endl;

    if (lifted_)
    {
      createSuper();
      createSuperNetwork();
    }
    else
    {
      createGround();
      createGroundNetwork();
    }
    
    if (fgdebug)
    {
      cout << "[init] ";
      Timer::printTime(cout, timer1.time());
      cout << endl;
      timer1.reset();
    }
  }

  /**
   * Prints out the network. The format is:
   * variables:
   * [var1]
   * [var2]
   * [var3]
   * ....
   * [varn]
   * factors:
   * [vari1] / [varj1] // [log table] /// [weight to vari1] [weight to varj1]
   * [vari2] / [varj2] / [vark2] // [log table] /// [weight to vari2]
   * [weight to varj2] [weight to vark2]
   * 
   * The first section of the file begins with a line containing only the text
   * "variables:".
   * Following which is list of of all the variable names. The only constraint
   * for variable names is that the variable name must not contain the '/'
   * symbol.
   *
   * The second section of file begins with a line containing only the text
   * "factors:".
   * Following which is list of of all the factors with one factor per line.
   * To demonstrate how the factor is stored, lets consider a factor F(X,Y,Z)
   * with binary variables.
   *
   * The part of the line to the left of the string '//' is a list of all the
   * variables in the factor, seperated by the '/' symbol.
   *
   * For instance, for the factor F(X,Y,Z), we will write
   * X / Y / Z // .....
   * The rest of the line is the table for log F(X,Y,Z) written in the iteration
   * order:
   * [log F(0,0,0)] [log F(1,0,0)] [log F(0,1,0)] [log F(1,1,0)] [log F(0,0,1)]
   * (basically increment the left most variable until it hits the limit, then
   * reset and increment the next variable and so on)
   * 
   * The part of each factor line to the right of '///' is a list of integer
   * weights which represent the weight of the edge connecting the factor to
   * each adjacent variable.
   * 
   */
  void printNetwork(ostream& out)
  {
    out << "variables:" << endl;
    for (int i = 0; i < nodes_->size(); i++)
    {
      (*nodes_)[i]->print(out);
      out << endl;
    }
    out << "factors:" << endl;
    for (int i = 0; i < factors_->size(); i++)
    {
      Factor* factor = (*factors_)[i];
      factor->print(out);
      out << "// ";
      factor->printWts(out);
      if (lifted_)
      {
        out << " ///";
        int numLinks = factor->getNumLinks();
        for (int lno = 0; lno < numLinks; lno++)
        {
          Link* link = factor->getLink(lno);
          Node* node = link->getNode();
          double cnt = node->getGroundNodeCount();
          out << " " << cnt;
        }
      }
      out << endl;
    }
  }

  /**
   * Returns the link to message map.
   */
  LinkIdToTwoWayMessageMap* getLinkIdToTwoWayMessageMap()
  {
    return lidToTWMsg_;
  }

  /**
   * Returns the number of nodes in the graph.
   */
  const int getNumNodes()
  {
    return nodes_->size();
  }
  
  /**
   * Returns the number of factors in the graph.
   */
  const int getNumFactors()
  {
    return factors_->size();
  }
  
  /**
   * Returns the number of auxiliary factors in the graph.
   */
  const int getNumAuxFactors()
  {
    if (auxFactors_ == NULL) return 0;
    return auxFactors_->size();
  }
  
  /**
   * Returns a node in the graph.
   * 
   * @param index Index of the node to be retrieved.
   */
  Node* getNode(const int& index)
  {
    return (*nodes_)[index];
  }

  /**
   * Returns a factor in the graph.
   * 
   * @param index Index of the factor to be retrieved.
   */
  Factor* getFactor(const int& index)
  {
    return (*factors_)[index];
  }

  /**
   * Returns an auxiliary factor in the graph.
   * 
   * @param index Index of the auxiliary factor to be retrieved.
   */
  AuxFactor* getAuxFactor(const int& index)
  {
    return (*auxFactors_)[index];
  }

  /**
   * Returns the domain on which the factor graph is built.
   */
  Domain* getDomain()
  {
    return domain_;
  }

 private:
 
 /**
  * Ground away the evidence from each clause. Then, create one 
  * superclause for each ground clause
  */
  void createGround()
  {
    MLN* mln = mln_;
    Domain* domain = domain_;
    Clause* mlnClause;

      //clause with all variables in it
    Clause *varClause;

    Array<Clause *> allClauses;
    Array<int> *mlnClauseTermIds;
  
    ClauseToSuperClauseMap *clauseToSuperClause;
    ClauseToSuperClauseMap::iterator clauseItr;
    int numClauses = mln->getNumClauses();
  
    PredicateTermToVariable *ptermToVar = NULL;
  
    double gndtime;
    Util::elapsed_seconds();
  
    clauseToSuperClause = new ClauseToSuperClauseMap();
    for (int i = 0; i < numClauses; i++)
    {
        //remove the unknown predicates
      mlnClause = (Clause*) mln->getClause(i);
      varClause = new Clause(*mlnClause);
      mlnClauseTermIds = varClause->updateToVarClause();
     
      clauseToSuperClause->clear();
      mlnClause->getConstantTuples(domain, domain->getDB(), mlnClauseTermIds, 
                                   varClause, ptermToVar, clauseToSuperClause,
                                   false);
      addSuperClauses(clauseToSuperClause);

        //can delete the var Clause now
      delete varClause;
    }

    gndtime = Util::elapsed_seconds();
    cout << "Grounding time is = " << gndtime << endl;

      //create the super preds  
    int totalGndClauseCnt = 0;
    for (int i = 0; i < superClausesArr_->size(); i++)
    {
      SuperClause *superClause = (*(*superClausesArr_)[i])[0];
      int gndCnt = superClause->getNumTuples();
      totalGndClauseCnt += gndCnt;
    }
    cout << "Total Number of Ground Clauses = " << totalGndClauseCnt << endl;
  }
 
 
  /**
   * Creates the ground network.
   */
  void createGroundNetwork()
  {
    Domain* domain = domain_;
    int domainPredCnt = domain->getNumPredicates();
    Array<IntArrayToIntMap*>* pconstantsToNodeIndexArr =
      new Array<IntArrayToIntMap*>();
    pconstantsToNodeIndexArr->growToSize(domainPredCnt);
     
      // Initialize the mappings from pred constants to the nodeIndices
    for(int i = 0; i < domainPredCnt; i++)
    {
      (*pconstantsToNodeIndexArr)[i] = new IntArrayToIntMap();
    }
     
    IntArrayToIntMap *pconstantsToNodeIndex;
    IntArrayToIntMap::iterator itr;
     
    Array<SuperClause*>* superClauses;
    SuperClause* superClause;

    Clause *clause;
    Predicate *pred;
    HyperCube *hyperCube;
    Array<Array<int> *> *constantsArr;
    Array<int>* constants;
    Array<int>* pconstants;
    
    Factor* factor;
    Node* node;

    for (int arrIndex = 0; arrIndex < superClausesArr_->size(); arrIndex++)
    {
      superClauses = (*superClausesArr_)[arrIndex];
      for (int scindex = 0; scindex < superClauses->size(); scindex++)
      {
        superClause = (*superClauses)[scindex];
        clause = superClause->getClause();
        int numHyperCubes = superClause->getNumHyperCubes();
        for (int hindex = 0; hindex < numHyperCubes; hindex++)
        {
          hyperCube = superClause->getHyperCube(hindex);
          double hcnt = superClause->getHyperCubeCount(hindex);
          constantsArr = hyperCube->getTuples();
          int numTuples = constantsArr->size();
          for (int tindex = 0; tindex < numTuples; tindex++)
          {
            // Number of times this tuple appears in the superclause
            // (because of counting elimination)
            double tcnt = hcnt;
            constants = (*constantsArr)[tindex];
            factor = new Factor(clause, NULL, constants, domain,
            					superClause->getOutputWt());  
            					
            if (clause->isActionFactor())
              factor->setAction();
            factors_->append(factor);
            for (int pindex = 0; pindex < clause->getNumPredicates(); pindex++)
            {
              pred = clause->getPredicate(pindex);
              pconstants = pred->getPredicateConstants(constants);
            
              int predId = domain->getPredicateId(pred->getName());
              pconstantsToNodeIndex = (*pconstantsToNodeIndexArr)[predId];
              itr = pconstantsToNodeIndex->find(pconstants);
              double util = clause->getUtil();
                 
              if (itr == pconstantsToNodeIndex->end())
              {
                int nodeIndex = nodes_->size();
                (*pconstantsToNodeIndex)[pconstants] = nodeIndex;
                    
                node = new Node(predId, NULL, pconstants, domain);
                node->setUtil(util);
                nodes_->append(node);
              }
              else
              {
                delete pconstants;
                int nodeIndex = itr->second;
                node = (*nodes_)[nodeIndex];
                util += node->getUtil();
                node->setUtil(util);
              }

                // Now, add the links to the node/factor
              int reverseNodeIndex = factor->getNumLinks();
                // Index where this factor would be stored in the list of nodes
              int reverseFactorIndex = node->getNumLinks();
              Link * link = new Link(node, factor, reverseNodeIndex,
                                     reverseFactorIndex, pindex, tcnt);
              node->addLink(link, NULL);
              factor->addLink(link, NULL);
            }
          }
        }
      }
    }
     
      // If query formulas present, make aux. links to the nodes
    if (auxFactors_)
    {
      for (int i = 0; i < auxFactors_->size(); i++)
      {
        AuxFactor* auxFactor = (*auxFactors_)[i];
        Array<Predicate* >* formula = auxFactor->getFormula();
        for (int pindex = 0; pindex < formula->size(); pindex++)
        {
          double z = 1.0;
          pred = (*formula)[pindex];
          pconstants = pred->getPredicateConstants();
          int predId = domain->getPredicateId(pred->getName());

          pconstantsToNodeIndex = (*pconstantsToNodeIndexArr)[predId];
          itr = pconstantsToNodeIndex->find(pconstants);

          if (itr == pconstantsToNodeIndex->end())
          {
            cout << "ERROR: couldn't find predicate ";
            pred->printWithStrVar(cout, domain);
            cout << " from query formula in factor graph" << endl;
            exit(-1);
          }
          else
          {
            delete pconstants;
            int nodeIndex = itr->second;
            node = (*nodes_)[nodeIndex];
          }

            // Now, add the links to the node/factor
          int reverseNodeIndex = auxFactor->getNumLinks();
            // Index where this factor would be stored in the list of nodes
          int reverseFactorIndex = node->getNumLinks();
          Link* link = new Link(node, auxFactor, reverseNodeIndex,
                                    reverseFactorIndex, pindex, z);
          node->addAuxLink(link);
          auxFactor->addLink(link, NULL);
        }
      }
    }

      // clean up
    for (int predId = 0; predId < domainPredCnt; predId++)
    {
      pconstantsToNodeIndex = (*pconstantsToNodeIndexArr)[predId];
      pconstantsToNodeIndex->clear();
      delete pconstantsToNodeIndex;
    }
    delete pconstantsToNodeIndexArr;
    cout << "Created Ground Network" << endl;
  }

 /**
  * Ground away the evidence from each clause. Then,
  * create the super clause/predicate network
  */
  void createSuper()
  {
    MLN* mln = mln_;
    Domain* domain = domain_;
    Clause* mlnClause;
    Array<Clause *> allClauses;

    ClauseToSuperClauseMap *clauseToSuperClause;
    ClauseToSuperClauseMap::iterator clauseItr;

    int numClauses = mln->getNumClauses();
       
    double hctime, jointime, refinetime, gndtime, bptime, setuptime;
    Util::elapsed_seconds();
    int predCnt = domain->getNumPredicates();
    Array<HyperCubeRefinement *> * hcRefinementArr = NULL;
    clauseToSuperClause = new ClauseToSuperClauseMap();
      
    /****************** HyperCubes ****************************/
    if (useHC_)
    {
      if (hcCreateType_ == Basic)
        cout<<"Creation Type is Basic..."<<endl;

      if (hcCreateType_ == DT)
        cout<<"Creation Type is DT..."<<endl;
     
      cout<<"Creation Noise = "<<hcCreateNoise_<<endl;
      cout<<"Number of Iterations to Run LNC for = "<<lncIter_<<endl;
     
      if (useCT_)
        cout<<"Using Constraints"<<endl;
      else
        cout<<"NOT Using Constraints"<<endl;

      Array<PredicateToHyperCubeReverseIndexMap*>* predToHCReverseIndexArr;
      if (useCT_)
        predToHCReverseIndexArr =
          createPredicateHyperCubesAndReverseIndexCT(mln, domain);
      else
        predToHCReverseIndexArr =
          createPredicateHyperCubesAndReverseIndex(mln, domain);

      hcRefinementArr = new Array<HyperCubeRefinement*>(predCnt, NULL);
      for (int predno = 0; predno < predCnt; predno++)
        (*hcRefinementArr)[predno] = new HyperCubeRefinement();
     
      hctime = Util::elapsed_seconds();

      cout<<"*****************************************************"<<endl;
      cout<<"                   Creating the Joins Now            "<<endl;
      cout<<"*****************************************************"<<endl;

      for (int i = 0; i < numClauses; i++)
      {
          //remove the unknown predicates
        mlnClause = (Clause*) mln->getClause(i);
        cout << "*************************************************"<<endl;
        cout << "Getting the join for the Clause "<<endl;
        mlnClause->print(cout,domain);
        cout<<endl;
        clauseToSuperClause->clear();
        getClauseToSuperClauseMap(mlnClause, predToHCReverseIndexArr,
                                  hcRefinementArr, domain, clauseToSuperClause);
        addSuperClauses(clauseToSuperClause);
      }
    
      jointime = Util::elapsed_seconds();
         
      cout<<"**********************************************************"<<endl;
      cout<<"               Done creating the Joins                    "<<endl;
      cout<<"**********************************************************"<<endl;

        /* delete the reverse index and the hypercubes in them */
      for (int tv = 0; tv < 3; tv++)
      {
        cout<<"delete indices for tv "<<tv<<endl;
        PredicateToHyperCubeReverseIndexMap *predToHCReverseIndex;
        PredicateToHyperCubeReverseIndexMap::iterator riItr;
        HyperCubeReverseIndex * hcReverseIndex;
       
        predToHCReverseIndex = (*predToHCReverseIndexArr)[tv];

        for(riItr = predToHCReverseIndex->begin();
            riItr != predToHCReverseIndex->end();
            riItr++)
        {
          Predicate *pred = riItr->first;
          pred->print(cout,domain);
          cout<<endl;
            
          hcReverseIndex = riItr->second;
          if (hcReverseIndex == NULL)
          {
            continue;
          }

          hcReverseIndex->deleteAllHyperCubes();
          delete hcReverseIndex;
        }
      }

      cout<<"Delete all the indices......."<<endl;
    
        /* HyperCube Refinement */
      cout<<endl<<"*****************************************************"<<endl;
      cout<<"************** Now Refining The HyperCubes ****************"<<endl;
      cout<<"***********************************************************"<<endl;
      for (int i = 0; i < hcRefinementArr->size(); i++)
      {
        cout<<"Processing the predicate = "<<domain->getPredicateName(i)<<endl;
        (*hcRefinementArr)[i]->createRefinedHyperCubeMapping();
      }
      refinetime = Util::elapsed_seconds();
     
      cout << "hctime = " << hctime << ", jointime = " << jointime
           << " refinetime = " << refinetime << endl;
      gndtime = hctime + jointime + refinetime;
    
     /******************** end code with hypercubes **************************/
    }
    else
    {
      /************* Original Code without hypercubes ***********************/

        //clause with all variables in it
      Clause *varClause;
      Array<int> *mlnClauseTermIds;
    
      PredicateTermToVariable *ptermToVar = NULL;
      if (implicitRep_)
      {
        ptermToVar = getPredicateTermToVariableMap(mln, domain);
        getIndexedConstants(ptermToVar, mln, domain);
      }
       
      for (int i = 0; i < numClauses; i++)
      {
          //remove the unknown predicates
        mlnClause = (Clause*) mln->getClause(i);
        varClause = new Clause(*mlnClause);
        mlnClauseTermIds = varClause->updateToVarClause();
    
        clauseToSuperClause->clear();
        mlnClause->getConstantTuples(domain, domain->getDB(), mlnClauseTermIds, 
                                     varClause, ptermToVar, clauseToSuperClause,
                                     implicitRep_);
        addSuperClauses(clauseToSuperClause);
          //can delete the var Clause now
        delete varClause;
      } 
      gndtime = Util::elapsed_seconds();
    }

    /****** End original code without hypercubes ********************/
      
      //set the hcrefinementarr - it is null for the code which does not use
      //hypercubes
    SuperPred::setHyperCubeRefinement(hcRefinementArr);

    cout<<"gnd time is = "<<gndtime<<endl;
    cout<<endl<<endl;
    cout<<"***************************************************"<<endl<<endl;
    cout<<"Now, starting the iterations of creating supernodes/superfeatures"
        <<endl;
      //now create the super preds corresponding to the current set of
      //superclauses

    SuperClause *superClause;
    int totalHyperCubeCnt = 0;
    int totalGndTupleCnt = 0;

    cout << "Counts in the beginning:" << endl;
    for (int i = 0; i < superClausesArr_->size(); i++)
    {
      superClause = (*(*superClausesArr_)[i])[0];
         
      int hcnt = superClause->getNumHyperCubes();
      int gndCnt = superClause->getNumTuples();
      double withDuplicateGndCnt = superClause->getNumTuplesWithDuplicates();

      Clause *clause = superClause->getClause();
      clause->print(cout,domain) << endl;
      cout << i << "=> hcnt :" << hcnt << ", ground cnt : " << gndCnt;
      cout << " gnd cnt (with Duplicates): " << withDuplicateGndCnt << endl;
      totalHyperCubeCnt += hcnt;
      totalGndTupleCnt += gndCnt;
    }
    cout<<"Total Number of Ground Tuples = "<<totalGndTupleCnt<<endl;
    cout<<"Total Number of HyperCubes Created = "<<totalHyperCubeCnt<<endl;

    /*************************************************************************
     * Start the Iterations now */
    /*************************************************************************/
       
    int newSuperClauseCnt = getNumArrayArrayElements(*superClausesArr_);
    int superClauseCnt = newSuperClauseCnt;
    int itr = 1;
    cout<<"*************************************************"<<endl<<endl;
    setuptime = 0;
    bptime = 0;

      //create the initial set of superpreds - one superpredicate for each
      //first order predicate
    initSuperPreds(superClausesArr_,domain);

    while (newSuperClauseCnt != superClauseCnt || itr <= 2)
    {      
      superClauseCnt = newSuperClauseCnt;
      cout<<"************************************************************"<<endl;
      cout<<"Iteration: "<<itr<<endl;
           
        //for iteration 1, superclauses have already been created
      if (itr > 1)
      {
        cout<<"Creating Super Clauses.. "<<endl;
        createSuperClauses(superClausesArr_, domain);
        newSuperClauseCnt = getNumArrayArrayElements(*superClausesArr_);
      }

      totalHyperCubeCnt = getHyperCubeCount(superClausesArr_);       
      cout<<"Total Number of HyperCubes Created = "<<totalHyperCubeCnt
          <<endl<<endl;
      
      cout<<"Number of superclauses after this iteration is = "
          <<newSuperClauseCnt<<endl;

      if(lncIter_ == itr)
      {
        cout<<endl<<"***************************************"<<endl;
        cout<<"Stopped Early At Iteration: "<<itr<<endl;
        itr++;
        bool isApproximateSuperPreds = true;
        cout<<"Creating Approximate Super Preds.. "<<endl;
        createSuperPreds(superClausesArr_,domain,isApproximateSuperPreds);
        break;
      }

      cout<<"Creating New Super Preds.. "<<endl;
      bool isApproximateSuperPreds = false;
      createSuperPreds(superClausesArr_, domain, isApproximateSuperPreds);
      setuptime += Util::elapsed_seconds(); 
      itr++;
    }
  
      //delete the hyperCube Refinement structures
    SuperPred::deleteHyperCubeRefinement();

    cout<<endl<<"***************************************"<<endl;
    cout<<"superclause setup time is = "<<setuptime<<endl;
    cout<<"total gnd time is = "<<(gndtime + setuptime)<<endl;

    cout<<"*********************************************"<<endl;
    cout<<"Number of Iterations For LNC = "<<--itr<<endl;
  
    Array<SuperClause *> *superClauses;
    int tupleCnt, hyperCubeCnt, constantCnt;

    tupleCnt = 0;
    superClauseCnt = 0;
    hyperCubeCnt = 0;
    constantCnt = 0;

    for(int i=0;i<superClausesArr_->size();i++)
    {
      superClauses = (*superClausesArr_)[i];
      for(int j=0;j<superClauses->size();j++)
      {
        superClause = (*superClauses)[j];
        tupleCnt += superClause->getNumTuples();
        superClauseCnt++;
        hyperCubeCnt += superClause->getNumHyperCubes();
        constantCnt += superClause->getNumConstants();
      }
    }
  
    cout<<"Total Number of Super Clauses = "<<superClauseCnt<<endl;
    cout<<"Total Number of HyperCubes = "<<hyperCubeCnt<<endl;
    cout<<"Total Number of Tuples = "<<tupleCnt<<endl;
    cout<<"Total Number of Constants = "<<constantCnt<<endl;
  
    cout<<endl<<"Memory-Results: ";
    cout<<"SuperClauseCnt "<<superClauseCnt<<" HyperCubeCnt "<<hyperCubeCnt;
    cout<<" TupleCnt "<<tupleCnt<<" ConstantCnt "<<constantCnt;
    cout<<endl<<endl;

    const PredicateTemplate *ptemplate;
    for(int predId=0;predId<predCnt;predId++)
    {
      ptemplate = domain->getPredicateTemplate(predId);
      if(ptemplate->isEqualPredicateTemplate())
        continue;
      int cnt = SuperPred::getSuperPredCount(predId);
      if(cnt > 0)
      {
        cout<<"SuperPred count for pred: ";
        ptemplate->print(cout);
        cout<<" = "<<cnt<<endl;
      }
    }
  }


  /**
   * Create the super network.
   */
  void createSuperNetwork()
  {
    Domain* domain = domain_;
    Array<SuperPred*> * superPreds;
    Array<SuperClause*> *superClauses;
  
    SuperClause *superClause;
    SuperPred *superPred;

    Factor *factor;
    Node *node;
    Clause *clause;
    Array<int>* constants = NULL;

    if (lncIter_ != 0)
    {
      cout<<"Working with an approximate network!"<<endl;
      Factor::setIsApproxNetwork(true);
    }

    //clean up - in case there were any previous runs
    for (int i = 0; i < nodes_->size(); i++)
      delete (*nodes_)[i];
    for (int i = 0; i < factors_->size(); i++)
      delete (*factors_)[i];
    nodes_->clear();
    factors_->clear();

      //create the factor (superclause) nodes
    for (int arrIndex = 0; arrIndex < superClausesArr_->size(); arrIndex++)
    {
      superClauses = (*superClausesArr_)[arrIndex];
      //cout<<arrIndex<<":"<<superClauses->size()<<endl;
      for (int scindex = 0; scindex < superClauses->size(); scindex++)
      {
        superClause = (*superClauses)[scindex];
        clause = superClause->getClause();
        factor = new Factor(clause, superClause, constants, domain,
                            superClause->getOutputWt());
        factors_->append(factor);
        if (clause->isActionFactor())
          factor->setAction();
      }
    }

      //create the variable (superpreds) nodes
    int predCnt = domain->getNumPredicates();
    for (int predId = 0;predId < predCnt; predId++)
    {
      double util = 0.0;
      Predicate* pred = domain->createPredicate(predId, false);
      if (pred)
      {
        Clause* unitClause = ClauseFactory::createUnitClause(pred, false);
        const Clause* newClause = mln_->findClause(unitClause);
        if (newClause)
        {
          util = newClause->getUtil();
          delete newClause;
        }
        delete pred;
      }
      
      superPreds = SuperPred::getSuperPreds(predId);
      //cout<<predId<<":"<<superPreds->size()<<endl;
      for (int spindex = 0; spindex < superPreds->size(); spindex++)
      {
        superPred = (*superPreds)[spindex];
        util = util * superPred->getNumHyperCubes();
        node = new Node(predId, superPred, constants, domain);
        node->addFactors(factors_, getLinkIdToTwoWayMessageMap());
        node->setUtil(util);
        nodes_->append(node);
      }
    }
  }
  
    //get the total number of hypercubes
  int getHyperCubeCount(Array<Array<SuperClause*>*>* const & superClausesArr_)
  {
    Array<SuperClause*>* superClauses;
    SuperClause* superClause;
    int totalHyperCubeCnt = 0;  
    for (int i = 0; i < superClausesArr_->size(); i++)
    {
      superClauses = (*superClausesArr_)[i];
      for(int j=0;j<superClauses->size();j++)
      {
        superClause = (*superClauses)[j];
        totalHyperCubeCnt += superClause->getNumHyperCubes();
      }
    }
    return totalHyperCubeCnt;
  }


  /**
   *  Manipulate the link id to message map
   */
  void updateLinkIdToTwoWayMessageMap()
  {
    Node *node;
    Link *link;
    Factor *factor;
    double nodeToFactorMsgs[2];
    double factorToNodeMsgs[2];

    LinkId *lid;
    TwoWayMessage *tmsg;
    LinkIdToTwoWayMessageMap::iterator lidToTMsgItr;
     
      // Delete the old values
    Array<LinkId*> keysArr;
    for (lidToTMsgItr = lidToTWMsg_->begin();
         lidToTMsgItr != lidToTWMsg_->end();
         lidToTMsgItr++)
    {
      keysArr.append(lidToTMsgItr->first);
      tmsg = lidToTMsgItr->second;
      delete tmsg;
    }
                       
    for (int i = 0; i < keysArr.size(); i++)
    {
      delete keysArr[i];
    }
    lidToTWMsg_->clear();

      // Now populate
    for (int i = 0; i < nodes_->size(); i++)
    {
      node = (*nodes_)[i];
      for (int j = 0; j < node->getNumLinks(); j++)
      {
        link = node->getLink(j);
        factor = link->getFactor();
              
        int predId = node->getPredId();
        int superPredId = node->getSuperPredId();
        int superClauseId = factor->getSuperClauseId();
        int predIndex = link->getPredIndex(); 
              
        lid = new LinkId(predId, superPredId, superClauseId, predIndex);

        int reverseFactorIndex = link->getReverseFactorIndex();
        node->getMessage(reverseFactorIndex, nodeToFactorMsgs);

        int reverseNodeIndex = link->getReverseNodeIndex();
        factor->getMessage(reverseNodeIndex, factorToNodeMsgs);
 
        tmsg = new TwoWayMessage(nodeToFactorMsgs,factorToNodeMsgs);
        (*lidToTWMsg_)[lid] = tmsg;
      }
    }
  }
  
  /**
   * Add the super clauses.
   */
  void addSuperClauses(ClauseToSuperClauseMap* const & clauseToSuperClause)
  {
    ClauseToSuperClauseMap::iterator clauseItr;
    Array<SuperClause *> * superClauses;
    SuperClause *superClause;
     
    for(clauseItr = clauseToSuperClause->begin();
        clauseItr != clauseToSuperClause->end(); 
        clauseItr++)
    {
      superClauses = new Array<SuperClause *>();
      superClausesArr_->append(superClauses);
      superClause = clauseItr->second;
      superClauses->append(superClause);
    }
  }
  
  
  /**
   * Merge the super clauses.
   */
  ClauseToSuperClauseMap*
  mergeSuperClauses(ClauseToSuperClauseMap* const & clauseToSuperClause)
  {
    Domain* domain = domain_;
    ClauseToSuperClauseMap *mergedClauseToSuperClause =
      new ClauseToSuperClauseMap();
    SuperClause *superClause, *mergedSuperClause;
    Clause *keyClause;
    ClauseToSuperClauseMap::iterator itr, mergedItr;
    for (itr = clauseToSuperClause->begin();
         itr != clauseToSuperClause->end();
         itr++)
    {
      superClause = itr->second;
      keyClause = superClause->getClause();
      mergedItr = mergedClauseToSuperClause->find(keyClause);
      if (mergedItr != mergedClauseToSuperClause->end())
      {
        mergedSuperClause = mergedItr->second;
        mergedSuperClause->merge(superClause);
        delete superClause;
      }
      else
      {
        (*mergedClauseToSuperClause)[keyClause] = superClause;
        keyClause->print(cout, domain);
        cout << endl;
      }
    }
    return mergedClauseToSuperClause;
  }
  
  /**
   * Creates the equivalence classes for the variables appearing in various
   * clauses.
   */
  PredicateTermToVariable* getPredicateTermToVariableMap(MLN * const & mln,
                                                         Domain* const & domain)
  {
    Clause *clause;
    Predicate *pred;
    const Term* term;
    Array<Variable *> *eqVars = new Array<Variable *>();
    int eqClassId = 0;
    Variable *var;
    const PredicateTemplate *ptemplate;
    const Array<int>* constants;

    for (int clauseno = 0; clauseno < mln->getNumClauses(); clauseno++)
    {
      clause = (Clause *)mln->getClause(clauseno);
      for (int predno = 0; predno < clause->getNumPredicates(); predno++)
      {
        pred = clause->getPredicate(predno);
        ptemplate = pred->getTemplate();
        for (int termno = 0; termno < pred->getNumTerms(); termno++)
        {
          term = pred->getTerm(termno);
          int varId = term->getId();
          int varTypeId = ptemplate->getTermTypeAsInt(termno);
          constants = domain->getConstantsByType(varTypeId); 
          var = new Variable(clause,varId,pred,termno,eqClassId,constants);
          eqVars->append(var);
          eqClassId++;
        }
      }
    }
    
    Variable  *var1, *var2;
    for (int i = 0; i < eqVars->size(); i++)
    {
      var1 = (*eqVars)[i];
      for (int j = i + 1; j < eqVars->size(); j++)
      {
        var2 = (*eqVars)[j];
        if (var1->same(var2))
        {
          var1->merge(var2); 
        }
      }
    }

      //now populate the map
    PredicateTermToVariable *ptermToVar = new PredicateTermToVariable();
    PredicateTermToVariable::iterator itr;
    PredicateTerm *pterm;
    Variable *tiedVar;
    int uniqueCnt = 0;
    for (int i = 0; i < eqVars->size(); i++)
    {
      var = (*eqVars)[i];
      if (var->isRepresentative())
      {
        uniqueCnt++;
        for (int j = 0; j < var->getNumTiedVariables(); j++)
        {
          tiedVar = var->getTiedVariable(j);
          int predId = tiedVar->getPredId();
          int termno = tiedVar->getTermno();
          pterm = new PredicateTerm(predId,termno);
          itr = ptermToVar->find(pterm);
          if (itr == ptermToVar->end())
          {
            (*ptermToVar)[pterm] = var;
          }
          else
          {
            delete pterm;
          }
        }
      }
    }

    cout << "size of PtermToVarMap is " << ptermToVar->size() << endl;
    cout << "count of Variable Eq Classes (Unique) is = " << uniqueCnt << endl;
    return ptermToVar;
  }


  void getIndexedConstants(PredicateTermToVariable * const & ptermToVar, 
                           MLN * const & mln, 
                           Domain * const & domain)
  {
    Predicate *pred, *gndPred;
    IntHashArray seenPredIds;
    const Clause *clause;
    const Term *term;
     
    Array<Predicate *> * indexedGndings;
    PredicateTerm *pterm;
    Database * db;
     
    int predId, termId, constantId;
    bool ignoreActivePreds = true;

    PredicateTermToVariable::iterator itr;
    Variable * var;

    indexedGndings = new Array<Predicate *>();
    db = domain->getDB();
    cout << "size of PtermToVarMap is " << ptermToVar->size() << endl;
     
    Clause *varClause;
    for (int clauseno = 0; clauseno < mln->getNumClauses(); clauseno++)
    {
      clause = mln->getClause(clauseno);    
      varClause = new Clause(*clause);
      varClause->updateToVarClause();

        //to make sure that we do not use clause
      clause = NULL;

      for (int predno = 0; predno < varClause->getNumPredicates(); predno++)
      {
        pred = varClause->getPredicate(predno);
        predId = pred->getId();

        if (seenPredIds.append(predId) < 0)
          continue;
        indexedGndings->clear();
          //Note: we assume that every predicate is indexable
        if (db->isClosedWorld(predId))
        {
            //precidate is closed world - rettrieve only true groundings
          db->getIndexedGndings(indexedGndings,pred,ignoreActivePreds,true);
        }
        else
        {
            //predicate is open world - retrieve both true and false groundings  
          db->getIndexedGndings(indexedGndings,pred,ignoreActivePreds,true);
          db->getIndexedGndings(indexedGndings,pred,ignoreActivePreds,false);
        }
        
        for (int gndno = 0; gndno < indexedGndings->size(); gndno++)
        {
          gndPred = (*indexedGndings)[gndno];

          for (int termno = 0; termno < gndPred->getNumTerms(); termno++)
          {
            pterm = new PredicateTerm(predId, termno);
            itr = ptermToVar->find(pterm);
            assert(itr != ptermToVar->end());
            var = itr->second;
            term = gndPred->getTerm(termno);
            constantId = term->getId();
            var->removeImplicit(constantId);
            delete pterm;
          }
          delete (*indexedGndings)[gndno];
        }
      }
      delete varClause;
    }
   
      //now explicitly handle the constants appearing in the clause
    for (int clauseno = 0; clauseno < mln->getNumClauses(); clauseno++)
    {
      clause = mln->getClause(clauseno);    
      for (int predno = 0; predno < clause->getNumPredicates(); predno++)
      {
        pred = clause->getPredicate(predno);
        predId = pred->getId();
        for (int termno = 0; termno < pred->getNumTerms(); termno++)
        {
          term = pred->getTerm(termno);
          termId = term->getId();
            // if it is a variable, nothing to do
          if (termId < 0) continue;
            // else, this constant also should be added the list of
            // indexed constants
          pterm = new PredicateTerm(predId, termno);
          itr = ptermToVar->find(pterm);
          assert(itr != ptermToVar->end());
          var = itr->second;
          var->removeImplicit(termId);
          delete pterm;
        }
      }
    }

    IntHashArray *seenEqClassIds = new IntHashArray();
    cout << "Implicit Set of constants are: " << endl;
    for (itr = ptermToVar->begin(); itr != ptermToVar->end(); itr++)
    {
      pterm = itr->first;
      var = itr->second;
      int eqClassId = var->getEqClassId();
      if (seenEqClassIds->find(eqClassId) >= 0)
        continue;
      seenEqClassIds->append(eqClassId);
      cout << "Implicit Constants for Eq class " << eqClassId << endl;
      cout << "Count =  " << var->getNumImplicitConstants() << " => " << endl;
      var->printImplicitConstants(cout, domain);
      cout << endl << endl << endl;
    }
    delete seenEqClassIds;
  }

  //remove those clauses from the MLN which only contain closed world predicates
  void removeTrivialClauses()
  {
    const Clause *clause;
    Predicate *pred;
    Array<const Clause*> *removeClauses = new Array<const Clause *>();
    Database *db;
    db = domain_->getDB();

    for(int clauseno=0;clauseno<mln_->getNumClauses();clauseno++)
    {
      clause = mln_->getClause(clauseno);    
      bool hasUnknown = false;
      for(int predno=0;predno<clause->getNumPredicates();predno++)
      {    
        pred = clause->getPredicate(predno);
        int predId = pred->getId();
        if(!db->isClosedWorld(predId))
        {
          hasUnknown = true;
        }
      }
      if(!hasUnknown)
        removeClauses->append(clause);
    }
       
    for(int i=0;i<removeClauses->size();i++)\
    {
      clause = (*removeClauses)[i];
      mln_->removeClause(clause);
    }
  }

/************* Methods for dealing with Hypercube Representation ********/

    //get the groundings of the equal predicate
  void getEqualPredicateIndexedGndings(Array<Predicate *> * const & predGndings,
                                       Predicate * const & pred,
                                       Domain * const & domain,
                                       bool areEqual)
  {
    cout<<endl<<endl<<"*******************************************"<<endl;
    cout<<"*******************************************************"<<endl;
    cout<<"For the Equal predicate (tv = "<<areEqual<<")";
    pred->print(cout,domain);
    cout<<endl;
    cout<<"Number of gndings = ";

    if(pred->isGrounded())
    {
      predGndings->append(new Predicate(*pred));
      cout<<predGndings->size()<<endl;
      return;
    }
       
    Predicate *predGnding;
    const PredicateTemplate * ptemplate = pred->getTemplate();
    assert(pred->getNumTerms() == 2);           
    assert(ptemplate->getTermTypeAsInt(0) == ptemplate->getTermTypeAsInt(1)); 
       
    int termId1 = pred->getTerm(0)->getId();
    int termId2 = pred->getTerm(1)->getId();
       
      //trivial case - variables are bound and they are "not equal" at the
      //same time
    if (termId1 == termId2 && !areEqual)
    {
      cout << predGndings->size() << endl;
      return;
    }
      
      //speical handling for constrained representation
    if (useCT_ && termId1 != termId2 && areEqual)
    {
      cout << predGndings->size() << endl;
      return;
    }

    IntArray *varConstants1 = new IntArray();
    IntArray *varConstants2 = new IntArray();
       
      //constants for first term
    int varTypeId = ptemplate->getTermTypeAsInt(0);
    if (termId1 < 0) 
      varConstants1->append(domain->getConstantsByType(varTypeId));
    else 
      varConstants1->append(termId1);
       
      //constants for second term
    if (termId2 < 0) 
      varConstants2->append(domain->getConstantsByType(varTypeId));
    else 
      varConstants2->append(termId2);

    Term *term;

    for (int index1 = 0; index1 < varConstants1->size(); index1++)
    {
      for (int index2 = 0; index2 < varConstants2->size(); index2++)
      {
        int constant1 = (*varConstants1)[index1];
        int constant2 = (*varConstants2)[index2];
             
        if (areEqual && constant1 != constant2)
          continue;
        if (!areEqual && (constant1 == constant2) && !useCT_) 
          continue;

        predGnding = new Predicate(*pred);
        term = (Term *)predGnding->getTerm(0);
        term->setId(constant1);
        term = (Term *)predGnding->getTerm(1);
        term->setId(constant2);
        predGndings->append(predGnding);
      }
    }

    cout << predGndings->size() << endl;
       
    delete varConstants1;
    delete varConstants2;
  }

public:

  //to get the unknown groundings of a predicate
  void getUnknownGndingsHash(PredicateHashArray * const & predGndingsHash,
                             Predicate * const pred,
                             Domain * const & domain,
                             Database * const & db)
  {
    Array<Predicate *> *predGndings = new Array<Predicate *>();
    getUnknownGndings(predGndings, pred, domain, db);
     
    for (int i = 0; i < predGndings->size(); i++)
    {
      predGndingsHash->append((*predGndings)[i]);
    }
    predGndings->clear();
    delete predGndings;
  }

private:

    //to get the unknown groundings of a predicate
  void getUnknownGndings(Array<Predicate*>* const & predGndings,
                         Predicate * const pred,
                         Domain * const & domain,
                         Database * const & db)
  {
    Array<Clause *> *unknownClauses = new Array<Clause *>();
    Clause *gndClause;
    Predicate *gndPred;
      //first create a clause from the predicate
    Clause *clause = new Clause();
    clause->appendPredicate(pred);

    clause->getUnknownClauses(domain, db, -1, NULL, NULL, NULL, unknownClauses);
    
    for (int i = 0; i < unknownClauses->size(); i++)
    {
      gndClause = (*unknownClauses)[i];
      gndPred = gndClause->getPredicate(0); 
      predGndings->append(gndPred);
      gndClause->removeAllPredicates();
      delete gndClause;
    }
    delete unknownClauses;
  }


  void populateTupleConstantsArray(Array<Predicate*>* const & predGndings,
                                   IntArrayHashArray* const & tupleConstantsArr)
  {
    Predicate *gndPred;
    int constantId;
    const Term *term;
    Array<int> * tupleConstants;
    tupleConstantsArr->clear();
    for (int gndno = 0; gndno < predGndings->size(); gndno++)
    {
      gndPred = (*predGndings)[gndno];
      tupleConstants = new Array<int>();
      tupleConstants->append(0); 
      for (int termno = 0; termno < gndPred->getNumTerms(); termno++)
      {
        term = gndPred->getTerm(termno);
        constantId = term->getId();
        tupleConstants->append(constantId);
      }
      tupleConstantsArr->append(tupleConstants);
    }
  }

  /**
   * Gets the hypercube corresponding to all the constants
   * in the domain of the predicate
   */
  HyperCube* getDomainHyperCube(Predicate* const & pred,
                                Domain* const & domain)
  {
    const PredicateTemplate * ptemplate = pred->getTemplate();
    int termCnt = ptemplate->getNumTerms();
    HyperCube * domainHyperCube = new HyperCube(termCnt);
    for (int termno = 0; termno < termCnt; termno++)
    {
      IntArray * varConstants = new IntArray();
      const Term *term = pred->getTerm(termno);
        //if this term is a constant, then append the constant to the list
        //of varConstants
      if (term->isConstant())
      {
        int constant = term->getId();
        varConstants->append(constant);
      }
      else
      { //else append the whole domain of the variable 
        int varTypeId = ptemplate->getTermTypeAsInt(termno);
        varConstants->append(domain->getConstantsByType(varTypeId)); 
      }
      domainHyperCube->setVarConstants(varConstants, termno + 1);
    }
          
      //if constrainted representation, duplicates variable occurrences should
      //be replaced by corresponding references
    if (useCT_)
    {
      IntArray *varConstants;
      IntToIntMap *refMap = getReferenceMap(pred);
      int varCnt = termCnt;
      for (int varId = 1; varId <= varCnt; varId++)
      {
        int refVarId = (*refMap)[varId];
        if (refVarId == varId) 
          continue;
        varConstants = new IntArray();
        varConstants->append(-refVarId);
        domainHyperCube->setVarConstantsAndDeleteOld(varConstants, varId);
      }
    }
    return domainHyperCube;
  }


  Array<PredicateToHyperCubeReverseIndexMap*>* 
  createPredicateHyperCubesAndReverseIndex(MLN* const & mln,
                                           Domain* const & domain)
  {
    IntHashArray *noHCPredIds = getPredicateIds(noHCPredsStr_,domain);
    cout << "NoHC Pred Ids Are =>" << endl;
    for (int i = 0; i < noHCPredIds->size(); i++)
    {
      cout << (*noHCPredIds)[i] << " ";
    }
    cout << "NOTE: NoHC Pred Ids are considered only for creating True "
         << "HyperCubes!"
         << endl;
     
    Array<PredicateToHyperCubeReverseIndexMap *> *predToHCReverseIndexArr =
      new Array<PredicateToHyperCubeReverseIndexMap *>(3, NULL);
    for (int i = 0; i < predToHCReverseIndexArr->size(); i++)
      (*predToHCReverseIndexArr)[i] = new PredicateToHyperCubeReverseIndexMap();
     
    int domainPredCnt = domain->getNumPredicates();
    Array<PredicateHashArray*>* canonicalPredsArr = 
      new Array<PredicateHashArray*>(domainPredCnt, NULL);
    for (int i = 0; i < domainPredCnt; i++)
    {
      (*canonicalPredsArr)[i] = new PredicateHashArray();
    }

    PredicateHashArray * trueOccurrencePreds = new PredicateHashArray();
    PredicateHashArray * canonicalPreds;

    Array<HyperCube *> *hyperCubes;
    PredicateToHyperCubeReverseIndexMap * predToHCReverseIndex;
    PredicateToHyperCubeReverseIndexMap::iterator riItr;
    HyperCubeReverseIndex *hcReverseIndex;

    Predicate *canonicalPred, *pred;
    Array<Predicate *> *predVariations;
    const Clause *clause;
    Database * db;
     
    bool ignoreActivePreds = true;
    db = domain->getDB();
     
      //first get all the unique occurrences of predicates as they appear in
      //the mln clauses (e.g. Friends(x,y) is treated different from
      //Friends(x,x))
    for (int clauseno = 0; clauseno < mln->getNumClauses(); clauseno++)
    {
      clause = mln->getClause(clauseno);    
       
      for (int predno = 0; predno < clause->getNumPredicates(); predno++)
      {
        pred = clause->getPredicate(predno);
        if (useCT_)
        {
          predVariations = getPredicateVariations(pred);
        }
        else
        {
          predVariations = new Array<Predicate *>();
          predVariations->append(new Predicate(*pred));
        }

        for (int i = 0; i < predVariations->size(); i++)
        {
          canonicalPred = (*predVariations)[i];
          canonicalPred->canonicalize();
          int predId = canonicalPred->getId();
          (*canonicalPredsArr)[predId]->append(canonicalPred);
         
          if (canonicalPred->getSense() == true)
            trueOccurrencePreds->append(canonicalPred);
        }
        delete predVariations;
      }
    }
     
    Array<Predicate *> * predGndings;
    IntArrayHashArray *tupleConstantsArr;
     
    predGndings = new Array<Predicate *>();
    tupleConstantsArr = new IntArrayHashArray();

    for (int predId = 0; predId < domainPredCnt; predId++)
    {
      canonicalPreds = (*canonicalPredsArr)[predId];
      for (int i = 0; i < canonicalPreds->size(); i++)
      {
        canonicalPred = (*canonicalPreds)[i];
          //the for loop is such that first true hypercubes are created
        for (int tv = 1; tv >= 0; tv--)
        {
            //if tv = 0, and this predicate does not appear anywhere in the MLN
            //as true, then do not need to create the false hypercubes
            //corresponding to this
          if (tv == 0)
          {
            if (trueOccurrencePreds->find(canonicalPred) < 0)
            {
              cout<<"Not creaing false hypercubes for pred "<<endl;
              canonicalPred->print(cout,domain);
              cout<<endl;
              continue;
            }
          }
               
          predGndings->clear();
          tupleConstantsArr->clear();

          if (canonicalPred->getTemplate()->isEqualPredicateTemplate())
          {
            getEqualPredicateIndexedGndings(predGndings, canonicalPred, domain,
                                            (bool)tv);
          }
          else
          {     
            db->getIndexedGndings(predGndings, canonicalPred, ignoreActivePreds,
                                  (bool)tv);
          }

            //nothing to do if there are no groundings for this combination
          if (predGndings->size() <= 0)
            continue;


            //special handling for closed world predicates, for the false
            //truth value the hypercubes are obtained by taking complement of
            //true hypercubes
            //- This is used only when predicate does not have any duplicate
            //variables. Otherwise, this method breaks down!!
          if (db->isClosedWorld(predId) && tv == 0 &&
              !containsDuplicateVariable(canonicalPred))
          {
            HyperCube * domainHyperCube;
            Array<HyperCube *> * trueHyperCubes;
               
            int tvTrue = 1;
            predToHCReverseIndex = (*predToHCReverseIndexArr)[tvTrue];
            riItr = predToHCReverseIndex->find(canonicalPred);
            if (riItr == predToHCReverseIndex->end())
            {
              trueHyperCubes = new Array<HyperCube *>();
            }
            else
            {
              hcReverseIndex = riItr->second;
              trueHyperCubes = hcReverseIndex->getAllHyperCubes();
            }
            domainHyperCube = getDomainHyperCube(canonicalPred, domain);
            hyperCubes = createComplementaryHyperCubes(trueHyperCubes,
                                                       domainHyperCube);
                
            cout<<"**************************************************"<<endl;
            cout<<"special handling for creating hypercubes .."<<endl;
            cout<<"size of hypercubes created = "<<hyperCubes->size()<<endl;
            delete domainHyperCube;
          }
          else
          {
            populateTupleConstantsArray(predGndings, tupleConstantsArr);
              
            cout<<"**************************************************"<<endl;
            cout << endl << endl << "Obtaining the Hypercubes for the "
                 << "predicate: (tv = " << tv << ")" << endl;
            canonicalPred->print(cout, domain);
            cout<<endl<<endl;
            cout<<"size of tuple constants = "<<tupleConstantsArr->size()<<endl;
            cout<<"**************************************************"<<endl;
            if (useCT_)
              replaceDuplicateConstantsByBindings(tupleConstantsArr,
                                                  canonicalPred);
              
              //currently done only for true values
            if (tv == TRUE && noHCPredIds->find(predId) >= 0)
            {
              cout<<"Getting the Default Ground HyperCubes!!"<<endl;
              hyperCubes = createGroundHyperCubes(tupleConstantsArr);  
            }
            else
            {
              hyperCubes = createHyperCubes(tupleConstantsArr, hcCreateType_,
                                            hcCreateNoise_);  
            }
            cout<<"size of hypercubes created = "<<hyperCubes->size()<<endl;
          }

          hcReverseIndex = new HyperCubeReverseIndex(hyperCubes);
          hcReverseIndex->createIndex();

          predToHCReverseIndex = (*predToHCReverseIndexArr)[tv];
          (*predToHCReverseIndex)[canonicalPred] = hcReverseIndex;
            
            //clear up
          predGndings->deleteItemsAndClear();
          tupleConstantsArr->deleteItemsAndClear();
        }
      } 
    }
     
     //now, populate the index for uknown groundings
     //- this is handled in slightly different way
     ////This is done to ensure disjoint property between
     //hypercubes coming from different partial groundings
     //of a predicate e.g. AdvisedBy(x,y), AdvisedBy(x,A) etc.
     //This property is not needed for evidence hypercubes
     int tvUnknown = 2;
     PredicateToHyperCubeArrayMap *predToHyperCubes;
     PredicateToHyperCubeArrayMap::iterator predItr;    
     
     cout<<"Getting unknown hypercubes.."<<endl;    
     for(int predId=0;predId<domainPredCnt;predId++)
     {
       if (db->isClosedWorld(predId))
         continue;
        
       canonicalPreds = (*canonicalPredsArr)[predId];
         //check the conditions for special handling:
         //All the predicate groundings are unknown
         //Few other conditions
       if (canonicalPreds->size() == 1)
       {
         bool unknownHCisDomainHC = true;
         canonicalPred = (*canonicalPreds)[0];
         if (containsDuplicateVariable(canonicalPred))
           unknownHCisDomainHC = false;
         for (int tv = 0; tv < 2; tv++)
         {
           predToHCReverseIndex = (*predToHCReverseIndexArr)[tv];
           riItr = predToHCReverseIndex->find(canonicalPred);
           if (riItr != predToHCReverseIndex->end())
             unknownHCisDomainHC = false;
         }
            
         if (unknownHCisDomainHC)
         {
           hyperCubes = new Array<HyperCube *>();
           HyperCube* domainHyperCube = getDomainHyperCube(canonicalPred,
                                                           domain);
           hyperCubes->append(domainHyperCube);
           hcReverseIndex = new HyperCubeReverseIndex(hyperCubes);
           hcReverseIndex->createIndex();
           predToHCReverseIndex = (*predToHCReverseIndexArr)[tvUnknown];
           (*predToHCReverseIndex)[canonicalPred] = hcReverseIndex;
           continue;
         }
      }

      predToHyperCubes = new PredicateToHyperCubeArrayMap();
      populateHyperCubes(canonicalPreds,predToHyperCubes,domain,db);

      for (predItr = predToHyperCubes->begin();
           predItr != predToHyperCubes->end();
           predItr++)
      {
        canonicalPred = predItr->first;
        hyperCubes = predItr->second;

        cout<<"For the predicate: "<<endl;
        canonicalPred->print(cout,domain);
        cout<<endl;
        cout<<"Number of hypercubes are = "<<hyperCubes->size();
        cout<<endl;

        hcReverseIndex = new HyperCubeReverseIndex(hyperCubes);
        hcReverseIndex->createIndex();
        predToHCReverseIndex = (*predToHCReverseIndexArr)[tvUnknown];
        (*predToHCReverseIndex)[canonicalPred] = hcReverseIndex;
      }
      delete predToHyperCubes;
    }

    delete predGndings;
    delete tupleConstantsArr;
    return predToHCReverseIndexArr;
  }

    //create the hypercubes for the given set of predicates (they come from
    //the same underlying predicate e.g. AdvisedBy - AdvisedBy(x,y),
    //AdvisedBy(x,A), AdvisedBy(x,x) etc.
  void populateHyperCubes(PredicateHashArray * const & canonicalPreds, 
                         PredicateToHyperCubeArrayMap* const & predToHyperCubes,
                          Domain* const & domain, Database * const & db)
  {
    Predicate *canonicalPred;
    Array<int> *tupleConstants;

    Array<Predicate *> * predGndings;
    IntArrayHashArray *tupleConstantsArr;
    IntArrayHashArray *allTupleConstantsArr;
    IntToIntArrayHashArrayMap *predIndexToTupleConstantsArr, 
                              *indexToTupleConstantsArr;
    IntToIntArrayHashArrayMap::iterator predIndexItr, indexItr;
    
    predIndexToTupleConstantsArr = new IntToIntArrayHashArrayMap();
    indexToTupleConstantsArr = new IntToIntArrayHashArrayMap();

    PredicateHashArray * newCanonicalPreds = new PredicateHashArray();  

    predGndings = new Array<Predicate *>();
    allTupleConstantsArr = new IntArrayHashArray();

    int predIndex = 0;

    cout<<"inside populate hypercube.."<<endl;
    for (int i = 0; i < canonicalPreds->size(); i++)
    {
      canonicalPred = (*canonicalPreds)[i];
      cout<<"Processing the Predicate ";
      canonicalPred->print(cout,domain);
      cout<<endl;
      predGndings->clear();
      getUnknownGndings(predGndings, canonicalPred, domain, db);
        
        //nothing to do if there are no groundings for this combination
      if (predGndings->size() <= 0)
        continue;
      cout<<"got non zero number of gndings.."<<endl;
      tupleConstantsArr = new IntArrayHashArray();
      newCanonicalPreds->append(canonicalPred);
      populateTupleConstantsArray(predGndings, tupleConstantsArr);
        //replace the duplicates
      if (useCT_)
      {
        replaceDuplicateConstantsByBindings(tupleConstantsArr, canonicalPred);
      }
      (*predIndexToTupleConstantsArr)[predIndex++] = tupleConstantsArr; 
      allTupleConstantsArr->append(tupleConstantsArr);
    }
    
    int index;
      //now place each tuple in the appropriate index
    for (int i = 0; i < allTupleConstantsArr->size(); i++)
    {
      tupleConstants = (*allTupleConstantsArr)[i];
      index = 0;
      for (predIndexItr = predIndexToTupleConstantsArr->begin();
           predIndexItr != predIndexToTupleConstantsArr->end();
           predIndexItr++)
      {
        predIndex = predIndexItr->first;
        tupleConstantsArr = predIndexItr->second;

        if (tupleConstantsArr->find(tupleConstants) >= 0)
        {
          index += (int)pow(2.0,predIndex);
        }
      }
      indexItr = indexToTupleConstantsArr->find(index);
      if (indexItr == indexToTupleConstantsArr->end())
      {
        tupleConstantsArr = new IntArrayHashArray();
        (*indexToTupleConstantsArr)[index] = tupleConstantsArr;
      }
      else
      {
        tupleConstantsArr = indexItr->second;
      }
      tupleConstantsArr->append(tupleConstants);
    }

    PredicateToHyperCubeArrayMap::iterator predItr;    
    Array<HyperCube *> * origHyperCubes, *hyperCubes;
     
      //now, create the hypercubes for tuple constants array in each index
    for (indexItr = indexToTupleConstantsArr->begin();
         indexItr != indexToTupleConstantsArr->end();
         indexItr++)
    {
      index = indexItr->first;
      tupleConstantsArr = indexItr->second;
      hyperCubes = createHyperCubes(tupleConstantsArr, hcCreateType_,
                                    hcCreateNoise_);  
      delete tupleConstantsArr;

      for (int predIndex = 0; predIndex < newCanonicalPreds->size(); 
           predIndex++)
      {
        if ((index>>predIndex & 1) > 0)
        {
          canonicalPred = (*newCanonicalPreds)[predIndex];
          predItr = predToHyperCubes->find(canonicalPred);
          if (predItr == predToHyperCubes->end())
          {
            origHyperCubes = new Array<HyperCube *>();
            (*predToHyperCubes)[canonicalPred] = origHyperCubes;
          }
          else
          {
            origHyperCubes = predItr->second;
          }

          for (int i = 0; i < hyperCubes->size(); i++)
          {
            origHyperCubes->append(new HyperCube((*hyperCubes)[i]));
          }
        }
      }

      for (int i = 0; i < hyperCubes->size(); i++)
      {
        (*hyperCubes)[i]->deleteVarConstants();
        delete (*hyperCubes)[i];
      }
    }
     
      //clean up
    for (predIndexItr = predIndexToTupleConstantsArr->begin();
         predIndexItr != predIndexToTupleConstantsArr->end();
         predIndexItr++)
    {
      tupleConstantsArr = predIndexItr->second;
      delete tupleConstantsArr;
    }

    delete predIndexToTupleConstantsArr;
    delete indexToTupleConstantsArr;

    delete newCanonicalPreds;
    delete predGndings;
    delete allTupleConstantsArr;
  }


    //for the constrained version - this one is much simpler
    //NOTE*****: May not work correctly if there are open world predicates
    //which appear with constants in the MLN clauses
    //********
  Array<PredicateToHyperCubeReverseIndexMap*>* 
  createPredicateHyperCubesAndReverseIndexCT(MLN* const & mln,
                                             Domain* const & domain)
  {
    assert(useCT_);
    
    IntHashArray *noHCPredIds = getPredicateIds(noHCPredsStr_, domain);
    cout<<"NoHC Pred Ids Are =>"<<endl;
    for (int i = 0; i < noHCPredIds->size(); i++)
    {
      cout<<(*noHCPredIds)[i]<<" ";
    }
    cout<<endl;

    Array<PredicateToHyperCubeReverseIndexMap*>* predToHCReverseIndexArr = 
      new Array<PredicateToHyperCubeReverseIndexMap *>(3, NULL);
    for (int i = 0; i < predToHCReverseIndexArr->size(); i++)
      (*predToHCReverseIndexArr)[i] = new PredicateToHyperCubeReverseIndexMap();
     
    int domainPredCnt = domain->getNumPredicates();
    Array<PredicateHashArray*>* canonicalPredsArr = 
      new Array<PredicateHashArray*>(domainPredCnt, NULL);
    for (int i = 0; i < domainPredCnt; i++)
      (*canonicalPredsArr)[i] = new PredicateHashArray();

    PredicateHashArray * trueOccurrencePreds = new PredicateHashArray();
    PredicateHashArray * canonicalPreds;

    Array<HyperCube *> *hyperCubes;
    PredicateToHyperCubeReverseIndexMap * predToHCReverseIndex;
    PredicateToHyperCubeReverseIndexMap::iterator riItr;
    HyperCubeReverseIndex *hcReverseIndex;

    Predicate *canonicalPred, *pred;
    Array<Predicate *> *predVariations;
    const Clause *clause;
    Database * db;
     
    bool ignoreActivePreds = true;

    db = domain->getDB();
     
      //first get all the unique occurrences of predicates as they appear in
      //the mln clauses (e.g. Friends(x,y) is treated different from
      //Friends(x,x))
    for (int clauseno = 0; clauseno < mln->getNumClauses(); clauseno++)
    {
      clause = mln->getClause(clauseno);    
       
      for (int predno = 0; predno < clause->getNumPredicates(); predno++)
      {    
        pred = clause->getPredicate(predno);
        predVariations = getPredicateVariations(pred);
         
        for (int i = 0; i < predVariations->size(); i++)
        {
          canonicalPred = (*predVariations)[i];
          canonicalPred->canonicalize();
          int predId = canonicalPred->getId();
          (*canonicalPredsArr)[predId]->append(canonicalPred);
         
          if (canonicalPred->getSense() == true)
            trueOccurrencePreds->append(canonicalPred);
        }
        delete predVariations;
      }
    }

    Array<Predicate *> * predGndings;
    IntArrayHashArray *tupleConstantsArr;
     
    predGndings = new Array<Predicate *>();
    tupleConstantsArr = new IntArrayHashArray();

    int tvSize = 3;
    IntArray *orderedTVs = new IntArray();
    IntArray *counts = new IntArray();
     
    orderedTVs->growToSize(tvSize);
    counts->growToSize(tvSize);
    Array<HyperCube *> *otherTVHyperCubes = new Array<HyperCube *>();
    HyperCube *domainHyperCube;

    for (int predId = 0; predId < domainPredCnt; predId++)
    {
      (*counts)[TRUE]  = db->getNumTrueGndPreds(predId);
      (*counts)[FALSE] = db->getNumFalseGndPreds(predId);     
      (*counts)[UNKNOWN] = db->getNumUnknownGndPreds(predId);
      for (int tv = 0; tv < 2;tv++)
      {
        if ((*counts)[tv] < 0)
        {
          cout<<"For Predicate "<<domain->getPredicateName(predId)<<", ";
          cout<<"Counts for tv = "<<tv<<" Overflow. Substituting Max Int "
              <<endl<<endl;
          (*counts)[tv] = MAX_INT;
        }
      }
         
      sortAscending(orderedTVs,counts);
      cout<<endl<<endl<<endl;
      cout<<"=========================================================="<<endl;
      cout<<"For Predicate "<<domain->getPredicateName(predId)<<endl;
      for (int index = 0; index < orderedTVs->size(); index++)
      {
        int tv = (*orderedTVs)[index];
        cout<<"orderedTVs["<<index<<"] = "<<tv<<" ** count = "
            <<(*counts)[tv]<<endl;
      }

      canonicalPreds = (*canonicalPredsArr)[predId];
      for (int i = 0; i < canonicalPreds->size(); i++)
      {
        canonicalPred = (*canonicalPreds)[i];
        cout<<endl<<endl;   
        otherTVHyperCubes->clear();
        for (int index = 0; index <= 2;index++)
        {
          int tv = (*orderedTVs)[index];
          cout<<"Processing ";
          canonicalPred->print(cout,domain);
          cout<<" (tv = "<<tv<<")"<<endl;
               
            //no need to process false groundings if there is no MLN clause
            //which would use them
          if (index == 2 && tv == FALSE && 
              trueOccurrencePreds->find(canonicalPred) < 0)
          {
            cout<<"**************************************************"<<endl;
            continue;
          }

            //if it is the last tv, then create hypercubes as complement of 
            //the hypercubes for other two truth values - if ground hypercubes
            //need to be created for this predicate, then do not create using 
            //complements
          if (index == 2 && noHCPredIds->find(predId) < 0)
          {
            cout<<"Obtaining the HyperCubes through complement "<<endl;
            domainHyperCube = getDomainHyperCube(canonicalPred, domain);
            hyperCubes = createComplementaryHyperCubes(otherTVHyperCubes,
                                                       domainHyperCube);
            cout<<"size of hypercubes created = "<<hyperCubes->size()<<endl;
            cout<<"**************************************************"<<endl;
            delete domainHyperCube;
          }
          else
          {
            predGndings->clear();
            tupleConstantsArr->clear();
               
            if (tv == UNKNOWN)
            {
              if (db->isClosedWorld(predId))
              {
                cout<<"**********************************************"<<endl;
                continue;
              }
              getUnknownGndings(predGndings,canonicalPred,domain,db);
            }
            else
            {
              if (canonicalPred->getTemplate()->isEqualPredicateTemplate())
                getEqualPredicateIndexedGndings(predGndings, canonicalPred,
                                                domain, (bool)tv);
              else
                db->getIndexedGndings(predGndings, canonicalPred,
                                      ignoreActivePreds, (bool)tv);
            } 

              //nothing to do if there are no groundings for this combination
            if (predGndings->size() <= 0)
            {
              continue;
            }

            populateTupleConstantsArray(predGndings, tupleConstantsArr);
             
            cout<<endl<<"Obtaining the Hypercubes "<<endl;
            canonicalPred->print(cout,domain);
            cout<<endl<<endl;
            cout<<"size of tuple constants = "<<tupleConstantsArr->size()
                <<endl;
               
            replaceDuplicateConstantsByBindings(tupleConstantsArr,
                                                canonicalPred);
               
            if (noHCPredIds->find(predId) >= 0)
            {
              cout<<"Getting the Default Ground HyperCubes!!"<<endl;
              hyperCubes = createGroundHyperCubes(tupleConstantsArr);  
            }
            else
            {
              hyperCubes = createHyperCubes(tupleConstantsArr, hcCreateType_,
                                            hcCreateNoise_);  
            }
               
            cout<<"size of hypercubes created = "<<hyperCubes->size()<<endl;
            cout<<"**************************************************"<<endl;
            otherTVHyperCubes->append(hyperCubes);
          }
           
          hcReverseIndex = new HyperCubeReverseIndex(hyperCubes);
          hcReverseIndex->createIndex();
               
          predToHCReverseIndex = (*predToHCReverseIndexArr)[tv];
          (*predToHCReverseIndex)[canonicalPred] = hcReverseIndex;
             
            //clear up
          predGndings->deleteItemsAndClear();
          tupleConstantsArr->deleteItemsAndClear();
        } 
      } 
    } 
    
    delete otherTVHyperCubes; 
    delete predGndings;
    delete tupleConstantsArr;
    return predToHCReverseIndexArr;
  }

  void getClauseToSuperClauseMap(Clause* const & mlnClause,
           Array<PredicateToHyperCubeReverseIndexMap*>* predToHCReverseIndexArr,
                           Array<HyperCubeRefinement*>* const & hcRefinementArr,
                                 Domain* const & domain,
                            ClauseToSuperClauseMap* const & clauseToSuperClause)
  {
    PredicateToHyperCubeReverseIndexMap::iterator riItr;
    bool predBit;
    int predCnt = mlnClause->getNumPredicates();
    Clause *varClause, *clause;
    Predicate *pred, *mlnPred, *canonicalPred;
    int tv;
    SuperClause *superClause;
    Array<int> *varIdToCanonicalVarId;
    Array<HyperCubeReverseIndex *> *hcReverseIndexArr =
      new Array<HyperCubeReverseIndex *>(predCnt, NULL);
    PredicateToHyperCubeReverseIndexMap *predToHCReverseIndex;

    int configCnt = (int)pow(2.0, predCnt);
    varClause = new Clause(*mlnClause);
    varClause->updateToVarClause();
    IntHashArray * unknownPredIndices = new IntHashArray();

    for (int n = 1; n < configCnt; n++)
    {
      clause = new Clause();
      unknownPredIndices->clear();
      for (int predno = 0; predno < predCnt; predno++)
      {
        pred = varClause->getPredicate(predno);
        bool sense = pred->getSense();
             
        predBit = n>>predno & 1;
        if (predBit)
        {
            //this predicate should be treated as unknown
          tv = 2;
          clause->appendPredicate(new Predicate(*pred));
          unknownPredIndices->append(predno);
        }
        else
        {
          tv = sense ? 0: 1;
        }
              
        predToHCReverseIndex = (*predToHCReverseIndexArr)[tv];

        mlnPred = mlnClause->getPredicate(predno);
        canonicalPred = new Predicate(*mlnPred);
        canonicalPred->canonicalize();

        riItr = predToHCReverseIndex->find(canonicalPred);
        delete canonicalPred;
        if (riItr == predToHCReverseIndex->end())
        {
          delete clause;
          clause = NULL;
          break;
        }
        else
        {
          (*hcReverseIndexArr)[predno] =  riItr->second;
        }
      }
        
      if (!clause || clause->getNumPredicates() <= 0)
      {
        continue;
      }

         /*
         cout<<"Getting the map for the simplified clause : "<<endl;
         clause->print(cout,domain);
         cout<<endl;
         */

         //make it a variable clause
         double wt = varClause->getWt();
         clause->setWt(wt); 

         int maxVarId = getMaxVarId(varClause); 
         varIdToCanonicalVarId = new Array<int>(maxVarId+1,-1);
         clause->canonicalize(varIdToCanonicalVarId);
         
         /*
         cout<<"Var clause is:"<<endl;
         clause->print(cout,domain);
         cout<<endl;
         */

         //try all the variations of this clause
         Clause *mlnVariationClause;
         ClauseHashArray *mlnVariationClauses;
         
         IntArrayHashArray *neqConstraints, *canonicalNeqConstraints;
         
         neqConstraints = NULL;
         canonicalNeqConstraints = NULL;
         if(useCT_) {
            neqConstraints = new IntArrayHashArray();
            mlnVariationClauses = getClauseVariations(mlnClause, unknownPredIndices,neqConstraints);
            
            canonicalNeqConstraints = getCanonicalConstraints(neqConstraints, varIdToCanonicalVarId);
            cout<<"constraints are :"<<endl;
            for(int i=0;i<canonicalNeqConstraints->size();i++) {
              cout<<i<<":";
              printArray(*(*canonicalNeqConstraints)[i], cout);
              cout<<endl;
          }
         } else {
            mlnVariationClauses = new ClauseHashArray();
            mlnVariationClauses->append(new Clause(*mlnClause));
         }
    
         superClause = new SuperClause(clause, canonicalNeqConstraints,
                                       varIdToCanonicalVarId, useCT_, wt);

         cout<<"MLN clause is "<<endl;
         mlnClause->print(cout,domain);
         cout<<endl<<" Variations Are :"<<endl;
         for(int i=0;i<mlnVariationClauses->size();i++) {
          mlnVariationClause = (*mlnVariationClauses)[i];
          //mlnVariationClause->canonicalize();
          cout<<i<<":"<<endl;
          mlnVariationClause->print(cout,domain);
          cout<<endl;
          
          bool hasZeroGndings = false;

          for(int predno=0;predno<predCnt;predno++) {
             mlnPred = mlnVariationClause->getPredicate(predno);
             bool sense = mlnPred->getSense();
             
             predBit = n>>predno & 1;
             if(predBit) {
                  tv = 2;
             } else {
                  tv = sense?0:1;
             }
             
             mlnPred->canonicalize();
             /*
             mlnPred->print(cout,domain);
             cout<<" ** tv = "<<tv<<endl;
             */

             predToHCReverseIndex = (*predToHCReverseIndexArr)[tv];

             canonicalPred = new Predicate(*mlnPred);
             canonicalPred->canonicalize();
             
             /*
             canonicalPred->print(cout,domain);
             cout<<" ** tv = "<<tv<<endl;
             */

             //cout<<"came here 1 ..."<<endl;
             riItr = predToHCReverseIndex->find(canonicalPred);
             //delete canonicalPred;
             //assert(riItr != predToHCReverseIndex->end());
             if(riItr == predToHCReverseIndex->end()) {
                  canonicalPred->print(cout,domain);
                  cout<<" ** tv = "<<tv<<endl;
                  cout<<"Could not find this one in index!!"<<endl;    
                  hasZeroGndings = true;
             } else {
             (*hcReverseIndexArr)[predno] =  riItr->second;
           }
            delete canonicalPred;
          }
          if(hasZeroGndings)
               continue;
          
          getClauseHyperCubes(varClause, superClause, hcReverseIndexArr, hcRefinementArr,unknownPredIndices,
                              neqConstraints, useCT_, domain);
         
         }

         cout<<"Done with the join :"<<endl;
         cout<<"Number of hypercubes is = "<<superClause->getNumHyperCubes()<<endl;
        
        
         /*
         for(int i=0;i<superClause->getNumHyperCubes();i++) {
              superClause->getHyperCube(i)->print(cout);
         }
         */
         cout<<endl;
         cout<<"Tuple Count is = "<<superClause->getNumTuples()<<endl;
         cout<<"Clause is = "<<endl;
         clause->print(cout,domain);
         cout<<endl;


         if(superClause->getNumTuples() == 0) {
              delete superClause;
              continue;
         } 

         //test for memory leak
        
         /*
         Array<HyperCube *> *hcArr;
         while(true) {
             cout<<"doing one more time for the clause "<<endl;
             clause->print(cout,domain);
             cout<<endl;
             hcArr = getClauseHyperCubes(varClause, NULL, hcReverseIndexArr, hcRefinementArr,unknownPredIndices,domain);
             cout<<"Number of hyperClauses = "<<hcArr->size()<<endl;
             hcArr->deleteItemsAndClear();
             delete hcArr;
         }*/

         SuperClause *existingSuperClause;
         ClauseToSuperClauseMap::iterator scItr = clauseToSuperClause->find(clause);
         //if this clause already exists, merge with corresponding superClause
         //otherwise create a new mapping
         if(scItr == clauseToSuperClause->end()) {
          (*clauseToSuperClause)[clause] = superClause;
         } else {
              cout<<"came to merge the clauses!!!"<<endl;
              existingSuperClause = scItr->second;
              existingSuperClause->merge(superClause);
              delete superClause;
         }
       }

       delete varClause;
       delete unknownPredIndices;
       delete hcReverseIndexArr;
  }

public:

//get all the possible not equal constraints for the terms in the
//given predicate - assumes that predicate has been variablized
IntArrayHashArray * getNeqConstraints(Predicate * const & pred) {
      const Term * term;
      IntToIntMap *varIdToTypeId = new IntToIntMap();
      IntToIntMap::iterator itr1, itr2;
      IntArrayHashArray *neqConstraints = new IntArrayHashArray();
      IntArray *constraint;

      int termCnt = pred->getNumTerms();
      for(int termno=0;termno<termCnt;termno++) {
           term = pred->getTerm(termno);
           int varId = -term->getId();
           assert(varId > 0);
           int typeId = pred->getTermTypeAsInt(termno);
           (*varIdToTypeId)[varId] = typeId;
      }

      for(itr1 = varIdToTypeId->begin();itr1 != varIdToTypeId->end(); itr1++) {
           int varId1 = itr1->first;
           int typeId1 = itr1->second;
       for(itr2 = varIdToTypeId->begin();itr2 != varIdToTypeId->end(); itr2++) {
           int varId2 = itr2->first;
           int typeId2 = itr2->second;
           if(typeId1 != typeId2 || varId1 >= varId2)
                continue;
           constraint = new IntArray();
           constraint->append(-varId1);
           constraint->append(-varId2);
           neqConstraints->append(constraint);
       }
      }
  delete varIdToTypeId;
  return neqConstraints;
}

private:
  
//get the mapping from original terms to the new terms, for the variation
//specified the SatCode and the constraints. 
//Note that we also need the constIds, because the constraints are defined 
//over the variablized structure
IntToIntMap *getVariationIdMap(IntArrayHashArray * const & neqConstraints, 
                     int constraintSatCode,
                     //variables which actually map to constants
                     IntHashArray * const & constIds,
                     //total number of ids
                     int idCnt) {
         
         IntToIntMap * origToVariationIds = new IntToIntMap();
         //initially each id maps to itself
         for(int i=1;i<=idCnt;i++) {
              (*origToVariationIds)[-i] = -i;
         }

         Array<IntHashArray *> *eqIdsArr = new Array<IntHashArray *>();
         IntHashArray *eqIds;
         
         IntArray *variationIdArr = new IntArray();
         IntArray * constraint;
         int constraintCnt = neqConstraints->size();
         for(int cno=0;cno<constraintCnt;cno++) {
             bool cbit = constraintSatCode>>cno & 1;
             if(cbit)
                  continue;
             constraint = (*neqConstraints)[cno];
             int id1 = (*constraint)[0];
             int id2 = (*constraint)[1];
             int foundPos = -1;
             for(int i=0;i<eqIdsArr->size();i++) {
                  eqIds = (*eqIdsArr)[i];
                  if(eqIds->find(id1) >= 0 || eqIds->find(id2) >= 0) {
                       eqIds->append(id1);
                       eqIds->append(id2);
                       foundPos = i;
                       break;
                  }
             }
             if(foundPos < 0) {
                  eqIds = new IntHashArray();
                  eqIds->append(id1);
                  eqIds->append(id2);
                  eqIdsArr->append(eqIds);
                  
                  //one of them is representative id
                  variationIdArr->append(id1);
                  foundPos = eqIdsArr->size()-1; 
             }
             
             //if it is a constant, then it should be mapped id
             if(constIds->find(id1) >= 0 ) {
                  (*variationIdArr)[foundPos] = id1;
             }
             if(constIds->find(id2) >= 0 ) {
                  (*variationIdArr)[foundPos] = id2;
             }
         }   
             
         for(int i=0;i<eqIdsArr->size();i++) {
                  eqIds = (*eqIdsArr)[i];
                  for(int j=0;j<eqIds->size();j++) {
                   int origId = (*eqIds)[j];
                   (*origToVariationIds)[origId] = (*variationIdArr)[i];
                  }
         }
        
         eqIdsArr->deleteItemsAndClear();
         
         delete eqIdsArr;
         delete variationIdArr;
         return origToVariationIds;  
}

Array<Predicate *> * getPredicateVariations(Predicate * const & pred) {

       //first create the variable Predicate
       IntArray *mlnTermIds;
       Clause *clause = new Clause();
       clause->appendPredicate(new Predicate(*pred));
       mlnTermIds = clause->updateToVarClause();
       Predicate *varPred = clause->getPredicate(0);
       
       clause->removeAllPredicates();
       delete clause;
       
       IntHashArray *constIds = new IntHashArray();
       //cout<<"orig term ids are = "<<endl;
       for(int i=1;i<mlnTermIds->size();i++) {
          //cout<<(*mlnTermIds)[i]<<" : "<<endl;
          if((*mlnTermIds)[i] >= 0)
               constIds->append(i);
       }

       //cout<<"trying to get the eq constraints.."<<endl;
       IntArrayHashArray * neqConstraints = getNeqConstraints(varPred);
       //cout<<"got neq constraints, Size = "<<neqConstraints->size()<<endl;
       
       IntToIntMap *origToVariationIds;
       Predicate *variationPred;
       Array<Predicate *> *variationPreds = new Array<Predicate *>();
       
       //for checking
       /*
       variationPreds->append(pred);
       return variationPreds;
       */

       int constraintCnt = neqConstraints->size();
       int configCnt= (int)pow(2.0,constraintCnt);

       //cout<<"config cnt is = "<<configCnt<<endl;
       int idCnt = mlnTermIds->size()-1;
       IntToIntMap::iterator itr;
       for(int n=0;n<configCnt;n++) {
         int constraintSatCode = n;
         origToVariationIds = getVariationIdMap(neqConstraints, constraintSatCode,constIds, idCnt);
        
         /*
         cout<<"mapping obtained is =>"<<endl;
         for(itr=origToVariationIds->begin();itr != origToVariationIds->end();itr++) {
             cout<<itr->first<<" : "<<itr->second<<endl;
         }*/

         variationPred = new Predicate(*varPred);
         int termCnt = varPred->getNumTerms();
         Term *term;
         for(int termno=0;termno<termCnt;termno++) {
          term = (Term *)variationPred->getTerm(termno);
          int origId = term->getId();
          int variationId = (*origToVariationIds)[origId];
          
          //if it is a constant
          if((*mlnTermIds)[-variationId] >= 0) {
               variationId = (*mlnTermIds)[-variationId];
          }
          term->setId(variationId);
         }
         variationPreds->append(variationPred);
         delete origToVariationIds;
       }

       delete mlnTermIds;
       delete constIds;
       delete varPred;
       
       neqConstraints->deleteItemsAndClear();
       delete neqConstraints;
       
       return variationPreds;
}


ClauseHashArray * getClauseVariations(Clause * const & mlnClause,
                                      IntHashArray * const & unknownPredIndices,
                                      IntArrayHashArray * const & neqConstraints) {
      //first create the variable Predicate
      IntArray *mlnTermIds;
      Clause *varClause = new Clause(*mlnClause); 
      Predicate *pred, *origPred;
      mlnTermIds = varClause->updateToVarClause();
       
      IntHashArray *constIds = new IntHashArray();
      for(int i=0;i<mlnTermIds->size();i++) {
          if((*mlnTermIds)[i] >= 0)
               constIds->append(i);
      }
      
      //Array<IntArray *> * neqConstraints = new Array<IntArray *>();
      IntArrayHashArray *predNeqConstraints;
      int predCnt = varClause->getNumPredicates();
      for(int predno=0;predno<predCnt;predno++) {
           pred = varClause->getPredicate(predno);
           
           //constraints imposed only for unknown preds and equality preds
           if((unknownPredIndices->find(predno) >= 0) ||
               pred->getTemplate()->isEqualPredicateTemplate()) {
            pred = varClause->getPredicate(predno);
            predNeqConstraints = getNeqConstraints(pred);
            neqConstraints->append(predNeqConstraints);
            delete predNeqConstraints;
           }
      }

      IntToIntMap *origToVariationIds;
      Predicate *variationPred;
      Clause *variationClause;
      ClauseHashArray *variationClauses = new ClauseHashArray();
      
      int constraintCnt = neqConstraints->size();
      /*
      cout<<"got neq constraints, Size = "<<neqConstraints->size()<<endl;
      cout<<"constraints are :"<<endl;
      for(int i=0;i<neqConstraints->size();i++) {
           printArray((*(*neqConstraints)[i]),cout);
           cout<<endl;
      }*/

      int configCnt=(int)pow(2.0,constraintCnt);

      int idCnt = mlnTermIds->size()-1;
      IntToIntMap::iterator itr;
      
      for(int n=0;n<configCnt;n++) {
         int constraintSatCode = n;
         origToVariationIds = getVariationIdMap(neqConstraints, constraintSatCode,constIds,idCnt);
         
         /*
         cout<<"mapping obtained (SatCode = "<<constraintSatCode<<") is =>"<<endl;
         for(itr=origToVariationIds->begin();itr != origToVariationIds->end();itr++) {
             cout<<itr->first<<" : "<<itr->second<<endl;
         }*/

         variationClause = new Clause();
         for(int predno=0;predno<predCnt;predno++) {
          origPred = varClause->getPredicate(predno);    
          
          variationPred = new Predicate(*origPred);
          int termCnt = origPred->getNumTerms();
          Term *term;
          for(int termno=0;termno<termCnt;termno++) {
           term = (Term *)variationPred->getTerm(termno);
           int origId = term->getId();
           int variationId = (*origToVariationIds)[origId];
           //if it is a constant
           if((*mlnTermIds)[-variationId] >= 0) {
               variationId = (*mlnTermIds)[-variationId];
           }
           term->setId(variationId);
         }
         variationClause->appendPredicate(variationPred);
         }
         variationClauses->append(variationClause);
         delete origToVariationIds;
       }
       
       delete mlnTermIds;
       delete constIds;
       delete varClause;
       
       //neqConstraints->deleteItemsAndClear();
       //delete neqConstraints;
       
       return variationClauses;
}


IntToIntMap *getReferenceMap(Predicate * const & canonicalPred) {
     int termCnt = canonicalPred->getNumTerms();
     IntToIntMap *refMap = new IntToIntMap();
     IntToIntMap * varIdToTermNo = new IntToIntMap();
     IntToIntMap::iterator itr;
     const Term * term;
     for(int termno=0;termno<termCnt;termno++) {
        term = canonicalPred->getTerm(termno);
        int varId = -term->getId();
        itr = varIdToTermNo->find(varId);
        if(itr == varIdToTermNo->end()) {
            (*varIdToTermNo)[varId] = termno;
        }
        (*refMap)[termno+1] = (*varIdToTermNo)[varId]+1;
     }
     return refMap;
}


void replaceDuplicateConstantsByBindings(IntArrayHashArray * tupleConstantsArr,
                                         Predicate * const & canonicalPred) {
     
     IntToIntMap * varIdToTermNo = new IntToIntMap();
     IntToIntMap::iterator itr;
     IntArray *tupleConstants;
     const Term * term;
     int termCnt = canonicalPred->getNumTerms();
     for(int termno=0;termno<termCnt;termno++) {
        term = canonicalPred->getTerm(termno);
        int varId = -term->getId();
        itr = varIdToTermNo->find(varId);
        if(itr == varIdToTermNo->end()) {
            (*varIdToTermNo)[varId] = termno;
        }
     }
    
     IntArrayHashArray *newTupleConstantsArr = new IntArrayHashArray();
     for(int i=0;i<tupleConstantsArr->size();i++) {
          tupleConstants = (*tupleConstantsArr)[i];
          assert(termCnt+1 == tupleConstants->size());
          for(int termno=0;termno<termCnt;termno++) {
               term = canonicalPred->getTerm(termno);
               int varId = -term->getId();
               int refTermNo = (*varIdToTermNo)[varId];
               //store the reference if it is not the original term
               if(refTermNo != termno) {
                 (*tupleConstants)[termno+1] = -(refTermNo+1);
                 //cout<<"need to replace by ref "<<refTermNo<<endl;
                 //printArray(*tupleConstants,cout);
                 //cout<<endl;
               } 
          }
               //printArray(*tupleConstants,cout);
               //cout<<endl;
               newTupleConstantsArr->append(tupleConstants);
     }
    
     tupleConstantsArr->clear();
     for (int i=0;i<newTupleConstantsArr->size();i++) {
          tupleConstantsArr->append((*newTupleConstantsArr)[i]);
     }
     newTupleConstantsArr->clear();
     delete newTupleConstantsArr;
     delete varIdToTermNo;
}



IntArrayHashArray* getCanonicalConstraints(IntArrayHashArray * const & constraints, 
                             IntArray * const & varIdToCanonicalVarId) {
     
     IntArrayHashArray *canonicalConstraints = new IntArrayHashArray();
     IntArray *constraint;
     IntArray *canonicalConstraint;
     int varId1, varId2, canonicalVarId1, canonicalVarId2;

     for(int i=0;i<constraints->size();i++) {
          constraint = (*constraints)[i];
          varId1 = -(*constraint)[0];
          varId2 = -(*constraint)[1];
          canonicalVarId1 = (*varIdToCanonicalVarId)[varId1];
          canonicalVarId2 = (*varIdToCanonicalVarId)[varId2];
          if(canonicalVarId1 < 0 || canonicalVarId2 < 0)
               continue;
          canonicalConstraint = new IntArray();
          if(canonicalVarId1 > canonicalVarId2) {
               int tmp = canonicalVarId1;
               canonicalVarId1 = canonicalVarId2;
               canonicalVarId2 = tmp;
          }
          canonicalConstraint->append(-canonicalVarId1);
          canonicalConstraint->append(-canonicalVarId2);
          canonicalConstraints->append(canonicalConstraint);
     }
     return canonicalConstraints;
}
  
 public:

  bool getUseHC()
  {
    return useHC_;
  }

  bool getUseCT()
  {
    return useCT_;
  }

  double getHcCreateNoise()
  {
    return hcCreateNoise_;
  }  
  
  MLN* getMLN()
  {
    return mln_;
  }
  
 private:
    // Indicates if lifted inference will be run
  bool lifted_;

    //whether to use the hypercube representation
  bool useHC_;

    //use the constrainted hypercube representation
  bool useCT_;
  
    //type of method used for hypercube creation
  HyperCubeCreateType hcCreateType_;

    //amt of noise tolerated during hypercube creation
  double hcCreateNoise_;

    //number of iterations of LNC to be performed
  int lncIter_;

    //list of predicates for which hypercubes should not be created
  string noHCPredsStr_;

    // Indicates if implicit representation is to be used
  bool implicitRep_;

  LinkIdToTwoWayMessageMap* lidToTWMsg_;

  Array<Array<SuperClause*>*>* superClausesArr_;

    // Factors in the graph  
  Array<Factor*>* factors_;
    // Nodes in the graph
  Array<Node*>* nodes_;

    // MLN from which the factor graph is built
  MLN* mln_;
    // Domain containing the constants from which the factor graph is built
  Domain* domain_;
  
    // Stores auxiliary factors used for query formulas
  Array<AuxFactor*>* auxFactors_;
};

#endif /*FACTORGRAPH_H_*/
