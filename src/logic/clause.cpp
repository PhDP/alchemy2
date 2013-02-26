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
#include "clause.h"
#include "mrf.h"
#include "variable.h"
#include "superclause.h"

ClauseSampler* Clause::clauseSampler_ = NULL;
double Clause::fixedSizeB_ = -1;
double AuxClauseData::fixedSizeB_ = -1;


  //check if addition of a constant with the given implicit index still keeps
  //the tuple a representative one. A tuple is a representative tuple if the
  //first occurrence of any given implicit constant is preceded by all the
  //constants which have index less than this.
inline bool isRepresentativePartialTuple(Array<int>* const & constants,
                                         int & implicitIndex,
                                         Array<Variable*>* const & eqVars,
                                         int varId)
{
  IntHashArray * seenConstants = new IntHashArray();
  Variable *var = (*eqVars)[-varId];
  for (int i = 0; i < constants->size(); i++)
  {
    int constantId = (*constants)[i];
    if(((*eqVars)[i] == var) && (constantId >= 0) && 
       var->isImplicit(constantId))
      seenConstants->append(constantId);
  }

    //it is sufficient to check the unique count as we know the input tuple is
    //representative
  bool representative = (implicitIndex <= seenConstants->size());
  delete seenConstants;
  return representative;
}

template <typename Type>
inline Array<Type>* getCanonicalArray(Array<Type>* const & arr,
                                     Array<int> * const & varIdToCanonicalVarId)
{
  Array<Type> * canonicalArr = new Array<Type>();
  int canonicalVarId;
  Type val;
  for (int varId = 0; varId < arr->size(); varId++)
  {
    val = (*arr)[varId];
    canonicalVarId = (*varIdToCanonicalVarId)[varId];
    if (canonicalVarId < 0)
      continue;
    if (canonicalArr->size() < canonicalVarId + 1)
      canonicalArr->growToSize(canonicalVarId + 1);
    (*canonicalArr)[canonicalVarId] = val;
  }
  return canonicalArr;
}

  //returns true if the (ground) clause has two literals with opposite sense
  //i.e. the clause is satisfied; otherwise returns false
bool Clause::createAndAddUnknownClause(
                                Array<GroundClause*>* const& unknownGndClauses,
                                Array<Clause*>* const& unknownClauses,
                                double* const & numUnknownClauses,
                                const AddGroundClauseStruct* const & agcs,
                                const Database* const & db)
{ 
  PredicateSet predSet; // used to detect duplicates
  PredicateSet::iterator iter;
  
  Clause* clause = NULL;
  for (int i = 0; i < predicates_->size(); i++)
  {
	Predicate* predicate = (*predicates_)[i];
    assert(predicate->isGrounded());
    if (db->getValue(predicate) == UNKNOWN)
    {
      if ( (iter=predSet.find(predicate)) != predSet.end() )
      {
          // the two gnd preds are of opp sense, so clause must be satisfied
        if ((*iter)->getSense() !=  predicate->getSense())
        {
          if (clause) delete clause;
          return true;
        }
        // since the two gnd preds are identical, no point adding a dup
        continue;
      }
      else
        predSet.insert(predicate);
      
      if (clause == NULL) clause = new Clause();
      Predicate* pred = new Predicate(*predicate, clause);
      clause->appendPredicate(pred);
    }
  }
  
  if (clause) 
  {
    if (numUnknownClauses) (*numUnknownClauses)++;

    clause->setWt(wt_);
    clause->canonicalizeWithoutVariables();
    
    if (agcs)
    {
      if (clausedebug >= 2)
      {
        cout << "Appending unknown clause to MRF ";
        clause->print(cout, db->getDomain());
        cout << endl;
      }
      MRF::addUnknownGndClause(agcs, this, clause, isHardClause_);
    }
    
    // MARC: The case with unknownGndClauses appears to be obsolete!
	if (unknownGndClauses)
    {
      if (clausedebug >= 2)
      {
        cout << "Appending unknown ground clause ";
        clause->print(cout, db->getDomain());
        cout << endl;
      }
      unknownGndClauses->append(new GroundClause(clause, agcs->gndPreds));
	  if (isHardClause_) unknownGndClauses->lastItem()->setWtToHardWt();
    }
    else if (unknownClauses)
    {
      if (clausedebug >= 2)
      {
        cout << "Appending unknown clause ";
        clause->print(cout, db->getDomain());
        cout << endl;
      }
      unknownClauses->append(clause);
      if (isHardClause_) clause->setIsHardClause(true);
    }
    if (unknownClauses == NULL) delete clause;
  }
  return false;
}

void addPredicateToHash(const Clause* const & c,
                        PredicateHashArray* const & predHashArray)
{
  int numPreds = c->getNumPredicates();
  // for each predicate in clause c
  for (int i = 0; i < numPreds; i++)
  {
    Predicate* pred = new Predicate(*(c->getPredicate(i)));
	int index = predHashArray->find(pred);
	if(index < 0 )
	{
	  index = predHashArray->append(pred) + 1;
	}
	else
	{
	  delete pred;
	  index++;
	}
  }	 
}


/**
 * Creates and adds a ground active clause.
 * 
 * @param activeGroundClauses If not NULL, then active GroundClauses are
 * accumulated here.
 * @param seenGndPreds GroundPredicates which have been seen before. Used when
 * accumulating GroundClauses.
 * @param db Database used to check truth values and evidence status of
 * predicates.
 * @param getSatisfied If true, satisfied clauses are also retrieved.
 */
bool Clause::createAndAddActiveClause(
                           Array<GroundClause *> * const & activeGroundClauses,
                           GroundPredicateHashArray* const& seenGndPreds,
		                   const Database* const & db,
                           bool const & getSatisfied)
{
  bool accumulateClauses = activeGroundClauses;
  Predicate *cpred;
  PredicateSet predSet; // used to detect duplicates
  PredicateSet::iterator iter;
 
  GroundClause *groundClause;
  
  Clause* clause = NULL;
  bool isEmpty = true;
  for (int i = 0; i < predicates_->size(); i++)
  {
    Predicate* predicate = (*predicates_)[i];
    assert(predicate); 
	assert(predicate->isGrounded());
    if ( (iter = predSet.find(predicate)) != predSet.end() )
    {
        // The two gnd preds are of opp sense, so clause must be satisfied
        // and no point in adding it 
	  if (wt_ >= 0 && !getSatisfied &&
          (*iter)->getSense() !=  predicate->getSense())
      {
        if (clause) delete clause;
		return false;
      }

        // Since the two gnd preds are identical, no point adding a dup
      continue;
	}
    else
      predSet.insert(predicate);
      
	bool isEvidence = db->getEvidenceStatus(predicate);
    
    if (clausedebug >= 2)
    {
      cout << "isEvidence " << isEvidence << endl;
      predicate->printWithStrVar(cout, db->getDomain());
      cout << endl;
    }
	if (!isEvidence)
      isEmpty = false;
	  
      // True evidence in a neg. clause: Discard clause
    if (wt_ < 0 && isEvidence && !getSatisfied &&
        db->sameTruthValueAndSense(db->getValue(predicate),
                                   predicate->getSense()))
    {
      if (clause) delete clause;
      return false;
    }
    
	  // Add only non-evidence prdicates
	if (accumulateClauses && !isEvidence)
	{
	  if (!clause) clause = new Clause();
        
	  cpred = new Predicate(*predicate, clause);
      assert(cpred);
      if (clausedebug >= 2)
      {
        cout << "Appending pred ";
        predicate->printWithStrVar(cout, db->getDomain());
        cout << " to clause ";
        clause->print(cout, db->getDomain());
        cout << endl;
      }
	  clause->appendPredicate(cpred);
      if (clausedebug >= 2) cout << "Appended pred to clause" << endl;
	}
  }

    // If the clause is empty after taking evidence into account, it should 
    // be discarded
  if (isEmpty)
  {
    if (clausedebug >= 2) cout << "Clause is empty" << endl;
	assert(!clause);
    return false;
  }
    // Came to this point means that the clause is active (and hence nonEmpty)
  else
  {
   	  // Add the corresponding ground clause if accumulateClauses is true
   	if (accumulateClauses)
   	{
      assert(clause);	
      if (clausedebug >= 2) cout << "Canonicalizing clause" << endl;
      clause->canonicalizeWithoutVariables();

      groundClause = new GroundClause(clause, seenGndPreds);
      if (isHardClause_)
        groundClause->setWtToHardWt();
      if (clausedebug >= 2) cout << "Appending ground clause to active set" << endl;
      activeGroundClauses->appendUnique(groundClause);
      delete clause;
    }
    return true;
  } 
}


double Clause::getConstantTuples(const Domain* const & domain,
                                 const Database* const & db,
                                 Array<int>* const & mlnClauseTermIds,
                                 const Clause* const & varClause,
                                 PredicateTermToVariable * const & ptermToVar,
                            ClauseToSuperClauseMap* const & clauseToSuperClause,
                                 bool useImplicit)
{
  bool ignoreActivePreds = true;
  double numTrueGndings = 0;

    //initialize the constants tuple with the var Ids - we will fill them in
    //with constants as we go along
  Array<int> *constants = new Array<int>(*mlnClauseTermIds);

    //Copy the literals so that their original order in the clause is
    //not affected by the subsequent sorting
    
  cout<<"***************************************************************"<<endl;
  cout<<"Came to process the clause : "<<endl;
  print(cout, domain);
  cout<<endl;

  createVarIdToVarsGroundedType(domain);

  Array<Predicate*>* origClauseLits = new Array<Predicate*>(*predicates_);

    // Array of partially grounded clauses achieved by using the inverted
    // index
  Array<Array<Predicate*>* > partGroundedClauses;

    //Assign to each variable in clause, the corresponding variable class
    //(stores the list of indexed constants). It is accessed through the
    //predicates appearing in the clause
  PredicateTermToVariable::iterator itr;
  PredicateTerm *pterm;
  Predicate *pred;
  const Term *term;
  Array<Variable *> *eqVars = new Array<Variable *>();
  eqVars->growToSize(mlnClauseTermIds->size());
    
    //initialize with NULL
  for (int i = 0; i < eqVars->size(); i++)
    (*eqVars)[i] = NULL;

  if (useImplicit)
  {
    for (int predno = 0; predno < varClause->getNumPredicates(); predno++)
    {
      pred = varClause->getPredicate(predno);
      int predId = pred->getId();
      for (int termno = 0; termno < pred->getNumTerms(); termno++)
      {
        term = pred->getTerm(termno);
        assert(!term->isConstant());
        int termId = term->getId();
        pterm = new PredicateTerm(predId, termno);
        itr = ptermToVar->find(pterm);
        assert(itr != ptermToVar->end());
        assert(-termId < eqVars->size());
        (*eqVars)[-termId] = itr->second;
        delete pterm;
      }
    }
  }

  if (useInverseIndex)
  {
      // Put the indexable literals first and ground them
    sortLiteralsByNegationAndArity(*origClauseLits, ignoreActivePreds, db);
    groundIndexableLiterals(domain, db, *origClauseLits, partGroundedClauses,
                            ignoreActivePreds);
  }
  else
  {
      //Sort preds in decreasing order of #TrueGndOfLiteral/#numOfGroundings.
      //The larger the number of true groundings of a literal, the more likely
      //it is to be true, so put it in front so that we can decide whether the
      //clause is true early.The larger the number of groundings of the
      //literal, the larger the savings when we decide that preceding literals
      //are true.
    sortLiteralsByTrueDivTotalGroundings(*origClauseLits, domain, db);
      // Put the original clause as the only clause into partGroundedClauses
    Array<Predicate*>* clauseLitsCopy = new Array<Predicate*>;
    clauseLitsCopy->growToSize(origClauseLits->size());
    for (int i = 0; i < origClauseLits->size(); i++)
      (*clauseLitsCopy)[i] = new Predicate(*(*origClauseLits)[i]);
    partGroundedClauses.append(clauseLitsCopy);
  }

    // At this point partGroundedClauses holds the nodes of the branch and
    // bound algorithm. This means nothing more is indexed and we must ground
    // out the rest of the predicates
  if (clausedebug)
  {
    cout << "Partially grounded clauses to be completed: " << endl;
    for (int pgcIdx = 0; pgcIdx < partGroundedClauses.size(); pgcIdx++)
    {
      cout << "\t";
      for (int i = 0; i < partGroundedClauses[pgcIdx]->size(); i++)
      {
        (*partGroundedClauses[pgcIdx])[i]->printWithStrVar(cout, domain);
        cout << " ";
      }
      cout << endl;
    }
  }

  bool skip;
    // Go through each clause in partGroundedClauses (nodes of the branch and
    // bound algorithm if using inverted index; otherwise, the original
    // clause), ground them out and check truth values
  for (int pgcIdx = 0; pgcIdx < partGroundedClauses.size(); pgcIdx++)
  {
      //intially, the list of constants is simply the mln term ids 
    constants->copyFrom(*mlnClauseTermIds);   

    skip = false;
      // clauseLits is a sorted copy of predicates_
    Array<Predicate*> clauseLits = *(partGroundedClauses[pgcIdx]);
    assert(clauseLits.size() == origClauseLits->size());
      // Set the var to groundings in this clause to be those in clauseLits
    Array<int>* origVarIds = new Array<int>;
      
    for (int i = 0; i < clauseLits.size(); i++)
    {
      assert(clauseLits[i]->getNumTerms() == 
             (*origClauseLits)[i]->getNumTerms());
        // Ground variables throughout clause
      for (int j = 0; j < (*origClauseLits)[i]->getNumTerms(); j++)
      {
        const Term* oldTerm = (*origClauseLits)[i]->getTerm(j);
        const Term* newTerm = clauseLits[i]->getTerm(j);
          
        if (oldTerm->getType() == Term::VARIABLE)
        {
          int varId = oldTerm->getId();
          origVarIds->append(varId);
          if (newTerm->getType() == Term::CONSTANT)
          {
            int constId = newTerm->getId();
            assert(constId >= 0);
            Array<Term*>& vars = (*varIdToVarsGroundedType_)[-varId]->vars;
            assert(constants->size() >= (-varId+1));

            if (useImplicit)
            {
              int implicitIndex = 
                (*eqVars)[-varId]->getImplicitIndex(constId);
              if (implicitIndex < 0)
              {
                (*constants)[-varId] = constId;
              }
              else
              {
                if (isRepresentativePartialTuple(constants, implicitIndex,
                                                 eqVars, varId))
                {
                  (*constants)[-varId] = constId;
                }
                else
                {
                  skip = true;
                }
              }
            }
            else
            {
              (*constants)[-varId] = constId;
            }

            for (int k = 0; k < vars.size(); k++) vars[k]->setId(constId);
          }
        }
      }
        // Set the preds in clauseLits to point to the original predicates_
      delete clauseLits[i];
      clauseLits[i] = (*origClauseLits)[i];
    }
      
    if (!skip)
    {
        //simulate a stack, back/front corresponds to top/bottom
        //ivg stands for index, varIds, groundings
      Array<LitIdxVarIdsGndings*> ivgArr;
      createAllLitIdxVarsGndings(clauseLits, ivgArr, domain, true);
      int ivgArrIdx = 0; //store current position in ivgArr
      bool lookAtNextLit = false;
    
        // while stack is not empty
      while (ivgArrIdx >= 0)
      {
          // get variable groundings at top of stack
        LitIdxVarIdsGndings* ivg = ivgArr[ivgArrIdx];
        Predicate* lit = (*origClauseLits)[ivg->litIdx];
        
        Array<int>& varIds = ivg->varIds;
        ArraysAccessor<int>& varGndings = ivg->varGndings;
        bool& litUnseen = ivg->litUnseen;
        bool hasComb;

        if (clausedebug)
        {
          cout << "Looking at lit: ";
          lit->printWithStrVar(cout, domain);
          cout << endl;
        }

        bool gotoNextComb = false;
          // while there are groundings of literal's variables
        while ((hasComb = varGndings.hasNextCombination()) || litUnseen)
        {
            // there may be no combinations if the literal is fully grounded
          if (litUnseen) litUnseen = false;

          if (hasComb)
          {
              //replace back the variables into the constants array
            for (int v = 0; v < varIds.size(); v++)
            {
              (*constants)[-varIds[v]] = varIds[v];
            }
              //ground the literal's variables throughout the clause
            int constId;
            int v = 0; // index of varIds
              //for each variable in literal
            gotoNextComb = false;
            while (varGndings.nextItemInCombination(constId))
            {
              int varId = varIds[v];
              Array<Term*>& vars = (*varIdToVarsGroundedType_)[-varId]->vars;
              
                //store the assignment to this variable
              assert(constants->size() >= (-varId+1));
   
              if (useImplicit)
              {
                int implicitIndex =
                  (*eqVars)[-varId]->getImplicitIndex(constId);
                if (implicitIndex < 0)
                {
                  (*constants)[-varId] = constId;
                }
                else
                {
                    //in case of implicit constant, proceed only if
                    //substitution of this constant
                    //will form a representative tuple - see the body of the
                    //function for the details of 
                    //what a representative tuple is
                  if (isRepresentativePartialTuple(constants, implicitIndex,
                                                   eqVars, varId))
                  {
                    (*constants)[-varId] = constId;
                  }
                  else
                  {
                    gotoNextComb = true;
                  }
                }
              }
              else
              {
                (*constants)[-varId] = constId;
              }
              v++;
              for (int i = 0; i < vars.size(); i++) vars[i]->setId(constId);
            }

              //a sanity check - for debugging purposes   
            assert(varIds.size() == v);
          }

          //removeRedundantPredicates();

          if (clausedebug)
          {
            cout << "Clause is now: ";
            printWithWtAndStrVar(cout, domain);
            cout << endl;
          }
          
          if (gotoNextComb)
            continue;
            
          if (literalOrSubsequentLiteralsAreTrue(lit, ivg->subseqGndLits, db))
          {
            if (clausedebug)
              cout << "Clause satisfied" << endl;
              
              //count the number of combinations of remaining variables
            double numComb = 1;
            for (int i = ivgArrIdx + 1; i < ivgArr.size(); i++)
            {
              int numVar = ivgArr[i]->varGndings.getNumArrays();
              for (int j = 0; j < numVar; j++)
                numComb *= ivgArr[i]->varGndings.getArray(j)->size();
            }
            numTrueGndings += numComb;
          }
          else
          {
              // if there are more literals
            if (ivgArrIdx + 1 < ivgArr.size())
            {
              if (clausedebug) cout << "Moving to next literal" << endl;
              lookAtNextLit = true;
              ivgArrIdx++; // move up stack
              break;
            }
              //At this point all the literals are grounded, and they are
              //either unknown or false (have truth values opposite of their
              //senses).

              //- this so that it matches the hypercube representation.
              //There is no real need to simplify the clause!
            if (false) 
            {
              ++numTrueGndings;
            }
            else
            {
                //Create a new constant tuple
              addConstantTuple(domain, db, varClause, constants, eqVars,
                               clauseToSuperClause, useImplicit);
            }
          }
        } //while there are groundings of literal's variables
        
          //if we exit the while loop in order to look at next literal 
          //(i.e. without considering all groundings of current literal)
        if (lookAtNextLit) { lookAtNextLit = false; }
          //mv down stack
        else
        { 
          varGndings.reset();
          litUnseen = true; 
          ivgArrIdx--; 
            //replace back the variables into the constants array
          for (int v = 0; v < varIds.size(); v++)
          {
            (*constants)[-varIds[v]] = varIds[v];
          }
        }
      } // while stack is not empty
      deleteAllLitIdxVarsGndings(ivgArr); 
    } //end skip
      
      // Restore variables
    for (int i = 0; i < origVarIds->size(); i++)
    {
      int varId = (*origVarIds)[i];
      assert(varId < 0);
      Array<Term*>& vars = (*varIdToVarsGroundedType_)[-varId]->vars;
      for (int j = 0; j < vars.size(); j++) vars[j]->setId(varId);
      (*varIdToVarsGroundedType_)[-varId]->isGrounded = false;
    }
      
    delete origVarIds;
    delete partGroundedClauses[pgcIdx];
  }
  delete origClauseLits;
  delete constants;
  return numTrueGndings;
}


//add a constant tuple to the list of constant tuples
inline void Clause::addConstantTuple(const Domain* const & domain,
                                     const Database* const & db,
                                     const Clause * const & varClause,
                                     Array<int> * const & constants,
                                     Array<Variable *> * const & eqVars,
                            ClauseToSuperClauseMap* const & clauseToSuperClause,
                                     bool useImplicit)
{
    // used to detect duplicates
  PredicateSet predSet; 
  
  PredicateSet::iterator iter;
  Predicate *predicate;
  
  Clause *clause = new Clause();

    //first make sure the clause is not satisfied (because of presence of two
    //identical gnd preds with opposite senses)
  for (int i = 0; i < predicates_->size(); i++)
  {
    predicate = (*predicates_)[i];
    assert(predicate->isGrounded());
	
	/*
    if (db->getValue(predicate) == UNKNOWN)
    {
      clause->appendPredicate(varClause->getPredicate(i));
      if ( (iter=predSet.find(predicate)) != predSet.end() )
      {
          // the two gnd preds are of opp sense, so clause must be satisfied
        if ((*iter)->getSense() !=  predicate->getSense())
        {
          clause->removeAllPredicates();
          delete clause;
          return;
        }
        // since the two gnd preds are identical, no point adding a dup
        continue;
      }
      else
        predSet.insert(predicate);
    } */
    
    //- this so that it matches the hypercube representation.
    //There is no real need to simplify the clause!
	if (db->getValue(predicate) == UNKNOWN)
    {
      clause->appendPredicate(varClause->getPredicate(i));
    }
  }
  
  SuperClause *superClause;
  ClauseToSuperClauseMap::iterator itr;
  Clause *keyClause;
  Array<int> * varIdToCanonicalVarId;
  Array<int> * canonicalConstants;
  Array<Variable *> * canonicalEqVars;

  //Slightly different Imeplementation for SuperClause
  bool useCT = false;
  IntArrayHashArray *neqConstraints = NULL;

    //nothing to do if the clause has no predicate
  if (clause->getNumPredicates() == 0)
  {
    delete clause;
    return;
  }

  itr = clauseToSuperClause->find(clause);
  if (itr == clauseToSuperClause->end())
  {
      //it is best to create a completely new copy of the clause
      //at this point
    keyClause = new Clause(*clause);
    keyClause->setWt(this->getWt());
    keyClause->setUtil(this->getUtil());
    keyClause->setAction(this->isActionFactor());
    //keyClause->setWt(1);

    int varCnt = constants->size();
    varIdToCanonicalVarId = new Array<int>(varCnt, -1);
    keyClause->canonicalize(varIdToCanonicalVarId);
    canonicalEqVars = getCanonicalArray(eqVars, varIdToCanonicalVarId);
      //Slightly different Imeplementation for SuperClause
    superClause = new SuperClause(keyClause, neqConstraints,
                                  varIdToCanonicalVarId, useCT, wt_);
   
    (*clauseToSuperClause)[clause] = superClause;
    delete canonicalEqVars;
    //delete varIdToCanonicalVarId;
  }
  else
  {
    superClause = itr->second;

      //sort of a hack, but it works: 
      //may need to increase the size of varIdToCanonicalVarId map
    int varCnt = constants->size();
    varIdToCanonicalVarId = superClause->getVarIdToCanonicalVarId();
    varIdToCanonicalVarId->growToSize(varCnt,-1);
      //delete the clause. Predicates should not be deleted, as
      //they belong to the varClause
    clause->removeAllPredicates();
    delete clause;
  }
  
  varIdToCanonicalVarId = superClause->getVarIdToCanonicalVarId();
  canonicalConstants = getCanonicalArray(constants, varIdToCanonicalVarId);
   //superClause->addNewConstantsAndIncrementCount(canonicalConstants,this->getWt());
  superClause->addNewConstantsAndIncrementCount(canonicalConstants,1);
    
    //clean up
  delete canonicalConstants;
}


//replaces all the constants in the clause by a new variable. Also, returns an 
//array of terms persent in the original clause (first the variable ids and 
//then the constant ids)
Array<int> * Clause::updateToVarClause()
{
  Array<int> * termIds = new Array<int>();
    //dummy values at position 0 - varIds start from 1
  termIds->append(0);

  const Predicate *pred;
     
    //first populate the variable ids to get a count
  for (int i = 0; i < predicates_->size(); i++)
  {
    pred = (*predicates_)[i];
    for (int j = 0; j < pred->getNumTerms(); j++)
    {
      const Term* t = pred->getTerm(j);
      if (t->getType() == Term::VARIABLE)
      {
        int id = t->getId();
        assert(id < 0);
        termIds->growToSize(-id+1);
        (*termIds)[-id] = id;
      }
    }
  }
    
    //now, convert all constants to variables while storing
    //their ids (constant values)
  for (int i = 0; i < predicates_->size(); i++)
  {
    pred = (*predicates_)[i];
    for (int j = 0; j < pred->getNumTerms(); j++)
    {
      Term* t = (Term *) pred->getTerm(j);
        //if its a variable - do not need to do anything
      if (t->getType() == Term::VARIABLE)
        continue;
      int constantId = t->getId();
      assert(constantId >= 0);
      termIds->append(constantId);
        
        //now convert this term to a new variable
      int varId = termIds->size()-1;
      t->setId(-varId);
    }
  }
  return termIds;
}

