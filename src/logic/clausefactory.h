/*
 * All of the documentation and software included in the
 * Alchemy Software is copyrighted by Stanley Kok, Parag
 * Singla, Matthew Richardson, Pedro Domingos, Marc
 * Sumner, Hoifung Poon, and Daniel Lowd.
 * 
 * Copyright [2004-07] Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, and Daniel Lowd. All rights reserved.
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
 * Poon, and Daniel Lowd in the Department of Computer Science and
 * Engineering at the University of Washington".
 * 
 * 4. Your publications acknowledge the use or
 * contribution made by the Software to your research
 * using the following citation(s): 
 * Stanley Kok, Parag Singla, Matthew Richardson and
 * Pedro Domingos (2005). "The Alchemy System for
 * Statistical Relational AI", Technical Report,
 * Department of Computer Science and Engineering,
 * University of Washington, Seattle, WA.
 * http://www.cs.washington.edu/ai/alchemy.
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
#ifndef CLAUSEFACTORY_H_NOV_10_2005
#define CLAUSEFACTORY_H_NOV_10_2005

#include <climits>
#include "clause.h"
#include "domain.h"


class ClauseFactory
{
 public:
  ClauseFactory() : maxNumVars_(INT_MAX), maxNumPredicates_(INT_MAX), 
                    domain_(NULL) {}
  ClauseFactory(const int& maxNumVars, const int& maxNumPredicates,
                const Domain* const & domain) 
      : maxNumVars_(maxNumVars), maxNumPredicates_(maxNumPredicates),
        domain_(domain) {}
  ~ClauseFactory() {}
  

    //Caller is responsible for deleting preds (and its contents), clause, and
    //the contents of newClauses.
  void addPredicateToClause(const Array<Predicate*>& preds, 
                            Clause* const & clause,
                            const int& op, const int& removeClauseIdx,
                            const bool& clauseHasBeenCanonicalized,
                            ClauseOpHashArray& newClauses,
                            const bool& addTruePredOnly)
  {
    if (!clauseHasBeenCanonicalized) clause->canonicalize();

    int smallestVarId;
    Array<Array<int>*>* typeIdToVarIds 
      = clause->getTypeIdToVarIdsMapAndSmallestVarId(smallestVarId);

    for (int i = 0; i < preds.size(); i++)
    {
      bool truePredOnly =
      	(preds[i]->isEqualPred() || preds[i]->isInternalPred()) ? true : addTruePredOnly;

      addPredicateToClause(preds[i], clause, op, removeClauseIdx, newClauses, 
                           typeIdToVarIds, smallestVarId, truePredOnly, NULL);
    }

    typeIdToVarIds->deleteItemsAndClear();
    delete typeIdToVarIds;
  }


    //Caller is responsible for deleting pred, clause and newClauses' contents.
  void addPredicateToClause(const Predicate* const & pred, 
                            Clause* const & clause, 
                            const int& op, const int& removeClauseIdx,
                            const bool& clauseHasBeenCanonicalized,
                            ClauseOpHashArray& newClauses,
                            const bool& addTruePredOnly)
  {
    Array<Predicate*> arr; 
    arr.append((Predicate*)pred);
      bool truePredOnly =
      	(pred->isEqualPred() || pred->isInternalPred()) ? true : addTruePredOnly;

    addPredicateToClause(arr, clause, op, removeClauseIdx, 
                         clauseHasBeenCanonicalized, newClauses, truePredOnly);
  }


    //Caller is responsible for deleting clause, and the contents of newClauses
  void removePredicateFromClause(const int& predIdx, 
                                 const Clause* const & clause,
                                 const int& op, const int& removeClauseIdx,
                                 ClauseOpHashArray& newClauses)
  {
    if (clause->getNumPredicates() == 1) return;
    Clause* newClause = new Clause(*clause);
    Predicate* p = newClause->removePredicate(predIdx);
    if (p == NULL) return;

    string prevClauseStr, addPredStr;
    if (domain_)
    {
      ostringstream oss1, oss2; 
      clause->printWithoutWtWithStrVar(oss1, domain_); prevClauseStr=oss1.str();
      p->printWithStrVar(oss2, domain_); addPredStr = oss2.str();
    }

    delete p;
    newClause->canonicalize();

    if (newClause->getNumVariablesAssumeCanonicalized() > maxNumVars_ ||
        newClause->getNumPredicates() > maxNumPredicates_) 
    { delete newClause; return; }
    if (!newClause->checkPredsAreConnected()) { delete newClause; return; }
    if (!validClause(newClause)) { delete newClause; return; }

    newClause->setAuxClauseData(new AuxClauseData(0, op, removeClauseIdx, 
                                                  prevClauseStr, addPredStr,
                                                  predIdx));
    if (((Clause*)clause)->containsConstants()) newClause->trackConstants();
    if (newClauses.append(newClause) < 0) { delete newClause; return; }
    newClause->compress();
  }


    //Caller is responsible for deleting clause, and the contents of newClauses
  void removePredicateFromClause(const Clause* const & clause,
                                 const int& op, const int& removeClauseIdx,
                                 ClauseOpHashArray& newClauses)
  {
    for (int i = 0; i < clause->getNumPredicates(); i++)
      removePredicateFromClause(i, clause, op, removeClauseIdx, newClauses);
  }


  void flipSensesInClause(const Clause* const & clause, const int& op, 
                          const int& removeClauseIdx,
                          ClauseOpHashArray& newClauses,
                          const bool& canonicalizeNewClauses)
  {
    Array<Array<bool>*> newedArrays; // for clean up later

    ArraysAccessor<bool> acc;
    for (int i = 0; i < clause->getNumPredicates(); i++)
    {
      Array<bool>* bArr = new Array<bool>(2); 
      bArr->append(false); bArr->append(true);
      newedArrays.append(bArr);
      acc.appendArray(bArr);      
    }

    Array<bool> bArr;
    while (acc.getNextCombination(bArr))
    {
      assert(bArr.size() == clause->getNumPredicates());
      Clause* newClause = new Clause(*clause);
      for (int i = 0; i < bArr.size(); i++)
      {
        Predicate* pred = (Predicate*) newClause->getPredicate(i);
        pred->setSense(bArr[i]);
      }
      
      int numVarsInNewClause;
      if (canonicalizeNewClauses)
      {
        newClause->canonicalize();
        numVarsInNewClause = newClause->getNumVariablesAssumeCanonicalized();
      }
      else        
        numVarsInNewClause = newClause->getNumVariables();

      if (numVarsInNewClause > maxNumVars_ ||
          newClause->getNumPredicates() > maxNumPredicates_)
      { delete newClause; continue; }
      if (!validClause(newClause)) { delete newClause; continue; }

      string prevClauseStr;
      if (domain_)
      {
        ostringstream oss; 
        clause->printWithoutWtWithStrVar(oss, domain_); 
        prevClauseStr = oss.str();
      }


      newClause->setAuxClauseData(new AuxClauseData(0, op, removeClauseIdx,
                                                    prevClauseStr, "", -1));
      if (((Clause*)clause)->containsConstants()) newClause->trackConstants();
      if (newClauses.append(newClause) < 0) { delete newClause; continue; }
      newClause->compress();
    }

    for (int i = 0; i < newedArrays.size(); i++) delete newedArrays[i];
  }


    // Caller is responsible for deleting the returned clause and ppred 
  static Clause* createUnitClause(const Predicate* const & predicate,
                                  const bool& allowEqualPred)
  {
    if (!allowEqualPred && 
    	(predicate->isEqualPred() || predicate->isInternalPred()))
      return NULL;
    Clause* clause = new Clause();
    Predicate* pred = new Predicate(*predicate);
    pred->setSense(true);
    pred->setParent(clause);
    clause->appendPredicate(pred);
    clause->canonicalize();
    return clause;
  }


    //Caller is responsible for deleting the contents of unitClauses and preds
  static void createUnitClauses(Array<Clause*>& unitClauses,
                                Array<Predicate*>& preds,
                                const bool& allowEqualPred)
  {
    for (int i = 0; i < preds.size(); i++)
    {
      Clause* uclause = createUnitClause(preds[i], allowEqualPred);
      if (uclause) unitClauses.append(uclause);
    }
  }


    //Caller is responsible for deleting pred and the contents of newClauses
  void createUnitClausesWithDiffCombOfVar(const Predicate* const & pred,
                                          const int& op, 
                                          const int& removeClauseIdx,
                                          ClauseOpHashArray& newClauses)
  {
    Clause* unitClause = createUnitClause(pred, true);

     Array<Predicate*> newPreds;
    int smallestVarId;
    Array<Array<int>*>* typeIdToVarIds 
      = unitClause->getTypeIdToVarIdsMapAndSmallestVarId(smallestVarId);
    ClauseOpHashArray unused;
    addPredicateToClause(pred, unitClause, OP_ADD, -1, unused, typeIdToVarIds, 
                         smallestVarId, true, &newPreds);
    typeIdToVarIds->deleteItemsAndClear();
    delete typeIdToVarIds;
    delete unitClause;
    
    for (int i = 0; i < newPreds.size(); i++)
    {
      Clause* uclause = createUnitClause(newPreds[i], true);
      delete newPreds[i];
      if (uclause == NULL) continue;
      if (uclause->getNumVariablesAssumeCanonicalized() > maxNumVars_ ||
          uclause->getNumPredicates() > maxNumPredicates_) 
      { delete uclause; continue; }

      uclause->setAuxClauseData(new AuxClauseData(0, op, removeClauseIdx));
      uclause->trackConstants();
      if (newClauses.append(uclause) < 0) { delete uclause; continue;}
      uclause->compress();
      
    }
  }


  void createUnitClausesWithDiffCombOfVar(const Array<Predicate*>& preds,
                                          const int& op, 
                                          const int& removeClauseIdx,
                                          ClauseOpHashArray& newClauses)
  {
    for (int i = 0; i < preds.size(); i++)
      createUnitClausesWithDiffCombOfVar(preds[i], op, removeClauseIdx, 
                                         newClauses);
  }


    //NOTE: add code to filter out clauses here
  bool validClause(const Clause* const & c)
  {
    if (hasEqualPredWithFreeVar(c)) return false;
    if (hasRedundantEqualPreds(c)) return false;
    if (c->getNumPredicates() > maxNumPredicates_) return false;

      //for testing UW-CSE domain
    //if (hasSamePredWithFreeVarSameVarWrongOrderVar(c)) return false; 
    //if (hasRedundantSamePreds(c)) return false;
    
    return true;
  }


 private:
    // If the returned predicate is added to a clause, its parent should be set
    // as the clause
  Predicate* createPredicate(const Predicate* const pred,
                             const Array<int>& termIds)
  {
    Predicate* newPred = new Predicate(*pred);
    assert(newPred->getNumTerms() == termIds.size());
    for (int i = 0; i < termIds.size(); i++)
      ((Term*) newPred->getTerm(i))->setId(termIds[i]);
    newPred->setSense(true);
    return newPred;
  }
  

  bool allNewVars(const Array<int>& termIds, const int& smallestVarId)
  {
    for (int i = 0; i < termIds.size(); i++)
      if (termIds[i] >=  smallestVarId) return false;
    return true;
  }


  void addPredicateToClause(const Predicate* const & pred, 
                            const Clause* const & clause, 
                            const int& op, const int& removeClauseIdx,
                            ClauseOpHashArray& newClauses,
                            const Array<Array<int>*>* const & typeIdToVarIds,
                            const int& smallestVarId,
                            const bool& addTruePredOnly, 
                            Array<Predicate*>* const returnNewPredOnly)
  {
    Array<Array<int>*> newedArrays; //for clean up later

    int newVarId = smallestVarId;
    Array<int>* intArr;
    ArraysAccessor<int> acc;
    for (int i = 0; i < pred->getNumTerms(); i++)
    {
      if (pred->getTerm(i)->getType() == Term::VARIABLE)
      {        
        int typeId = pred->getTermTypeAsInt(i);
        if (typeId < typeIdToVarIds->size() && (*typeIdToVarIds)[typeId]) 
          intArr = new Array<int>( *((*typeIdToVarIds)[typeId]) );
        else  
          intArr = new Array<int>; 
        newedArrays.append(intArr);
        intArr->append(--newVarId);
        acc.appendArray(intArr);
      }
      else
      {
        assert(pred->getTerm(i)->getType() == Term::CONSTANT);
        intArr = new Array<int>; 
        newedArrays.append(intArr);
        intArr->append(pred->getTerm(i)->getId());
        acc.appendArray(intArr);
      }
    }

    int numComb = acc.numCombinations();
    int comb = 0;
    Array<int> termIds;
    while(acc.getNextCombination(termIds))
    {
        //prevent creating both P(a,b) v a = b and P(a,b) v b = a
      if (pred->isEqualPred() && termIds[0] <= termIds[1]) continue;

        //in last combination, pred contains all new variables, so ignore it
      if (++comb==numComb) {assert(allNewVars(termIds,smallestVarId));continue;}

      Predicate* newPred = createPredicate(pred, termIds);
      if (returnNewPredOnly) { returnNewPredOnly->append(newPred); continue; }

      Predicate* fnewPred = new Predicate(*newPred);

        //check newPred is not in clause
      if (clause->containsPredicate(newPred)) 
        { delete newPred; delete fnewPred; continue; }
      Clause* newClause = new Clause(*clause); 
      newClause->appendPredicate(newPred);
      newPred->setParent(newClause);
      newClause->canonicalize();

      if (newClause->getNumVariablesAssumeCanonicalized() > maxNumVars_ ||
          newClause->getNumPredicates() > maxNumPredicates_) 
      { delete newClause; delete fnewPred; continue; }

      if (!validClause(newClause)) { delete newClause;delete fnewPred;continue;}

      string prevClauseStr, addPredStr;
      if (domain_)
      {
        ostringstream oss1, oss2; 
        clause->printWithoutWtWithStrVar(oss1, domain_); 
        prevClauseStr = oss1.str(); 
        pred->printWithStrVar(oss2, domain_); 
        addPredStr = oss2.str();
      }


      newClause->setAuxClauseData(new AuxClauseData(0, op, removeClauseIdx,
                                                    prevClauseStr, addPredStr,
                                                    -1));
      if (((Clause*)clause)->containsConstants()) newClause->trackConstants();
      if (newClauses.append(newClause) < 0) { delete newClause; }
      else                                  { newClause->compress(); }

      if (addTruePredOnly) { delete fnewPred; continue; }

      //Predicate* fnewPred = new Predicate(*newPred);
      fnewPred->setSense(false);
      newClause = new Clause(*clause);
      newClause->appendPredicate(fnewPred);
      fnewPred->setParent(newClause);
      newClause->canonicalize();

      newClause->setAuxClauseData(new AuxClauseData(0, op, removeClauseIdx,
                                                    prevClauseStr, addPredStr,
                                                    -1));
      if (((Clause*)clause)->containsConstants()) newClause->trackConstants();
      if (newClauses.append(newClause) < 0) { delete newClause; }
      else                                  { newClause->compress(); }
    }

    newedArrays.deleteItemsAndClear();
  }


 private:
  //commented out: the problem of same var and wrong order is taken care of 
  //               in addPredicateToClause().
  /*
    // returns true if clause c contains an '=' pred with the same 
    // variable/constant (e.g., a = a), with a free var 
    // (e.g., person(a) v same(a,b)), or with the vars in the wrong order 
    // (e.g. advisedBy(a1,a2) v a2=a1). The last check prevents us from 
    // creating duplicate clauses like advisedBy(a1,a2) v a1=a2  and
    // advisedBy(a1,a2) v a2=a1.
  bool hasEqualPredWithFreeVarSameVarWrongOrderVar(const Clause* const & c)
  {
    for (int i = 0; i < c->getNumPredicates(); i++)
    {
      const Predicate* p = c->getPredicate(i);
      if (p->isEqualPred())
      {
        if (p->getTerm(0)->getId() >= p->getTerm(1)->getId()) return true;

        for (int j = 0; j < p->getNumTerms(); j++)
        {
          int tid = p->getTerm(j)->getId();
          bool bound = false;
          for (int k = 0; k < c->getNumPredicates(); k++)
          {
            if (i == k) continue;
            const Predicate* pp = c->getPredicate(k);
            for (int l = 0; l < pp->getNumTerms(); l++)
              if (pp->getTerm(l)->getId() == tid) { bound = true; break; }
            if (bound) break;
          }
          if (!bound) return true;
        }
      }
    }
    return false;
  }
  */

    // returns true if clause c contains an '=' pred with a free var 
    // (e.g., person(a) v a = b).
  bool hasEqualPredWithFreeVar(const Clause* const & c)
  {
    for (int i = 0; i < c->getNumPredicates(); i++)
    {
      const Predicate* p = c->getPredicate(i);
      if (p->isEqualPred())
      {
        for (int j = 0; j < p->getNumTerms(); j++)
        {
          int tid = p->getTerm(j)->getId();
          bool bound = false;
          for (int k = 0; k < c->getNumPredicates(); k++)
          {
            if (i == k) continue;
            const Predicate* pp = c->getPredicate(k);
            for (int l = 0; l < pp->getNumTerms(); l++)
              if (pp->getTerm(l)->getId() == tid) { bound = true; break; }
            if (bound) break;
          }
          if (!bound) return true;
        }
      }
    }
    return false;
  }
  
  
    //returns true if clause c contains two predicates of the form
    //(a = b) v (b = a)
  bool hasRedundantEqualPreds(const Clause* const & c)
  {
    for (int i = 0; i < c->getNumPredicates(); i++)
    {
      const Predicate* p = c->getPredicate(i);
      if (p->isEqualPred())
      {
        for (int j = i+1; j < c->getNumPredicates(); j++)
        {
          const Predicate* pp = c->getPredicate(j);
          if (p->getId() == pp->getId() && 
              p->getTerm(0)->getId() == pp->getTerm(1)->getId() &&
              p->getTerm(1)->getId() == pp->getTerm(0)->getId()) return true;
        }
      }
    }
    return false;
  }


    //for testing UW-CSE domain. 
    // returns true if clause c contains an 'same' pred with the same 
    // variable/constant (e.g., samePerson(a,a) or a pred with a free var 
    // (e.g., person(a) v samePerson(a,b)), or with the vars in the wrong order 
    // (e.g. advisedBy(a1,a2) v same(a2,a1)). The last check prevents us from 
    // creating duplicate clauses like advisedBy(a1,a2) v same(a1,a2)  and
    // advisedBy(a1,a2) v same(a2,a1).
  bool hasSamePredWithFreeVarSameVarWrongOrderVar(const Clause* const & c)
  {
    for (int i = 0; i < c->getNumPredicates(); i++)
    {
      const Predicate* p = c->getPredicate(i);
      const char* predName = domain_->getPredicateName(p->getId());
      if (strncmp(predName, "same", 4) == 0)
      {
        if (p->getTerm(0)->getId() >= p->getTerm(1)->getId()) return true;

        for (int j = 0; j < p->getNumTerms(); j++)
        {
          int tid = p->getTerm(j)->getId();
          bool bound = false;
          for (int k = 0; k < c->getNumPredicates(); k++)
          {
            if (i == k) continue;
            const Predicate* pp = c->getPredicate(k);
            for (int l = 0; l < pp->getNumTerms(); l++)
              if (pp->getTerm(l)->getId() == tid) { bound = true; break; }
            if (bound) break;
          }
          if (!bound) return true;
        }
      }
    }
    return false;
  }


    //for testing UW-CSE domain.  
    //returns true if clause c contains two predicates of the form
    //samePerson(a,b) v (a = b),  samePerson(a,b) v (b = a), 
    //samePerson(a,b) v samePerson(b,a), (a = b) v (b = a),
  bool hasRedundantSamePreds(const Clause* const & c)
  {
    for (int i = 0; i < c->getNumPredicates(); i++)
    {
      const Predicate* p = c->getPredicate(i);
      if (p->isEqualPred())
      {
        for (int j = i+1; j < c->getNumPredicates(); j++)
        {
          const Predicate* pp = c->getPredicate(j);
          if (p->getId() == pp->getId())
          {
            if (p->getTerm(0)->getId() == pp->getTerm(1)->getId() &&
                p->getTerm(1)->getId() == pp->getTerm(0)->getId()) return true;
          }
          else
          {
            const char* predName = domain_->getPredicateName(pp->getId());
            if (strncmp(predName, "same", 4) == 0)
            {
              if ( (p->getTerm(0)->getId() == pp->getTerm(0)->getId() &&
                    p->getTerm(1)->getId() == pp->getTerm(1)->getId()) ||
                   (p->getTerm(0)->getId() == pp->getTerm(1)->getId() &&
                    p->getTerm(1)->getId() == pp->getTerm(0)->getId()) )
                return true;
            }
          }
        }        
      }
      else
      {
        const char* predName = domain_->getPredicateName(p->getId());
        if (strncmp(predName, "same", 4) == 0)
        {
          for (int j = i+1; j < c->getNumPredicates(); j++)
          {
            const Predicate* pp = c->getPredicate(j);
            if (p->getId() == pp->getId())
            {
              if (p->getTerm(0)->getId() == pp->getTerm(1)->getId() &&
                  p->getTerm(1)->getId() == pp->getTerm(0)->getId())return true;
            }
            else
            if (pp->isEqualPred())
            {
              if ( (p->getTerm(0)->getId() == pp->getTerm(0)->getId() &&
                    p->getTerm(1)->getId() == pp->getTerm(1)->getId()) ||
                   (p->getTerm(0)->getId() == pp->getTerm(1)->getId() &&
                    p->getTerm(1)->getId() == pp->getTerm(0)->getId()) )
                return true;
            }
          }
        }
      }
    }
    return false;
  }


 private:
  int maxNumVars_;
  int maxNumPredicates_;
  const Domain* domain_; //not owned by this object, do not delete

};


#endif
