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
#ifndef MLN_H_JUN_26_2005
#define MLN_H_JUN_26_2005

#include <cfloat>
#include <ext/hash_map>
#include <ext/hash_set>
using namespace __gnu_cxx;
#include "clause.h"
#include "equalstr.h"
#include "mlnhelper.h"


class MLN
{
 public: 
  MLN() : clauses_(new ClauseHashArray),clauseInfos_(new Array<MLNClauseInfo*>),
          formAndClausesArray_(new FormulaAndClausesArray), 
          predIdToClausesMap_(new Array<Array<IndexClause*>*>),
          externalClause_(new Array<bool>), hybridClauses_(new Array<Clause*>),
          hybridFormulaStrings_(new Array<string>),
          numericPreds_(new Array<Predicate*>),
          numericMeans_(new Array<double>){}

  ~MLN()
  {
    if (clauses_)
    {
      clauses_->deleteItemsAndClear();
      delete clauses_;
    }
    
    if (clauseInfos_)
    {
      clauseInfos_->deleteItemsAndClear();
      delete clauseInfos_;
    }
    
    if (formAndClausesArray_)
    {
      formAndClausesArray_->deleteItemsAndClear();
      delete formAndClausesArray_;
    }
    
    if (predIdToClausesMap_)
    {
      for (int i = 0; i < predIdToClausesMap_->size(); i++)
        if ((*predIdToClausesMap_)[i]) 
        {
          (*predIdToClausesMap_)[i]->deleteItemsAndClear();
          delete (*predIdToClausesMap_)[i];
        }
      delete predIdToClausesMap_;
    }

    if (externalClause_) delete externalClause_;

    if (hybridClauses_)
    {
      hybridClauses_->deleteItemsAndClear();
      delete hybridClauses_;
    }

    if (hybridFormulaStrings_) delete hybridFormulaStrings_;

    if (numericPreds_)
    {
      numericPreds_->deleteItemsAndClear();
      delete numericPreds_;
    }

    if (numericMeans_) delete numericMeans_;
  }

  int getNumClauses() const { return clauses_->size(); }

  int getNumHybridClauses() const { return hybridClauses_->size(); }

  int getNumHardClauses() const
  {
  	int numHardClauses = 0;
  	for (int i = 0; i < clauses_->size(); i++)
  	{
  	  if ((*clauses_)[i]->isHardClause()) numHardClauses++;
  	}
  	return numHardClauses;
  }


  bool containsClause(const Clause* const & c) const
  { return clauses_->contains((Clause*)c); }

  /**
   * Append a clause to the mln and mark it as external. It is assumed that
   * the clause is in CNF.
   */
  bool appendExternalClause(const string& formulaString, const bool& hasExist,
                            Clause* const& c, const Domain* const & domain,
                            const bool& tiedClauses,
                            const bool& hasWeightFullStop,
                            const bool& isConjunction)
  {
    int idx;
    bool app = appendClause(formulaString, hasExist, c, c->getWt(),
                            c->isHardClause(), idx, tiedClauses,
                            hasWeightFullStop, isConjunction, 0.0);
    if (app)
    {
      setFormulaNumPreds(formulaString, c->getNumPredicates());
      setFormulaIsHard(formulaString, c->isHardClause());
      setFormulaPriorMean(formulaString, c->getWt());
      (*externalClause_)[idx] = true;
    }
    return app;   
  }

    //MLN owns clause c and is responsible for its deletion.
    //hasExist is true if formulaString contains an existential quantifier
    //during the CNF conversion process.
    //Call setHardClausesWts() after this function to set the weights of hard
    //clauses to twice the maximum of soft clause wts.
    //Returns true of the clause is appended; otherwise returns false.
    //retIdx is the index of the clause in MLN's array of clauses
  bool appendClause(const string& formulaString, const bool& hasExist, 
                    Clause* const& c, const double& wt,
                    const bool& isHardClause, int& retClauseIdx,
                    const bool& tiedClauses, const bool& hasWeightFullStop,
                    const bool& conjunction, const double& util)
  {
    assert(c);
    if (!tiedClauses) c->canonicalize();
    bool isAppended;
    Clause* clause;
    if ((retClauseIdx = clauses_->find(c)) >= 0) // if clauses_ contains c
    {
      clause = (*clauses_)[retClauseIdx];
      delete c;
      isAppended = false;
    }
    else
    {
        // clause does not exist in MLN, so add it
      retClauseIdx = clauses_->append(c);
      clause = c;
      isAppended = true;
      clause->setWt(0); //start accumulating weight from 0
      externalClause_->append(false);
    }
    clause->addWt(wt);
    clause->addUtil(util);
    if (util != 0.0) clause->setAction(true);
    if (isHardClause) clause->setIsHardClause(isHardClause);
    if (hasWeightFullStop) clause->setStaticWt(hasWeightFullStop);

 
    FormulaAndClauses* fac = new FormulaAndClauses(formulaString, 
                                                   formAndClausesArray_->size(),
                                                   hasExist, tiedClauses,
                                                   conjunction);
    int fidx;
      // if this is the first time we see the formulaString
    if ( (fidx = formAndClausesArray_->find(fac)) < 0 )
    {
      fidx = formAndClausesArray_->append(fac);
      assert(fac->index == fidx);
    }
    else
    {
      delete fac;
      fac = (*formAndClausesArray_)[fidx];
      assert(fac->index == fidx);
    }
    IndexClause* ic = new IndexClause(fac->indexClauses->size(), clause);
    int cidx = fac->indexClauses->append(ic);
    if (cidx < 0) { delete ic; ic = NULL; }
    else          assert(ic->index == cidx);

    if (isAppended) 
    {
        // add to predIdToClausesMap_, and append to clauseInfos_
      appendClauseInfo(clause, retClauseIdx, &(fac->index), &(ic->index));
    }
    else 
    if (cidx >= 0) 
    {   
        //a clause was appended to an entry in formAndClausesArray_, so update 
        //the formulaClauseIndexes of the MLNClauseInfo corresponding to clause
      updateClauseInfo(retClauseIdx, &(fac->index), &(ic->index));
    }

    assert(clauses_->size() == clauseInfos_->size());
    return isAppended;
  } //appendClause()


    //MLN owns clause c and is responsible for its deletion.
    //Returns true of the clause is appended; otherwise returns false.
    //retIdx is the index of the clause in MLN's array of hybrid clauses
  bool appendHybridClause(const string& formulaString, Clause* const& c,
                          const double& wt, int& retClauseIdx,
                          Predicate* const& pred, const double& mean)
  {
    assert(c);
    c->canonicalize();
    bool isAppended;
    Clause* clause;
      // if hybridClauses_ contains c
    if ((retClauseIdx = hybridClauses_->find(c)) >= 0 &&
        strcmp(formulaString.c_str(),
               (*hybridFormulaStrings_)[retClauseIdx].c_str()) == 0)
    {
      clause = (*hybridClauses_)[retClauseIdx];
      delete c;
      isAppended = false;
    }
    else
    {
        // hybrid clause does not exist in MLN, so add it
      retClauseIdx = hybridClauses_->append(c);
      hybridFormulaStrings_->append(formulaString);

      //numericTerms_->append(continuousString);
      numericPreds_->append(pred);
      numericMeans_->append(mean);
      clause = c;
      isAppended = true;
      clause->setWt(0); //start accumulating weight from 0
    }
    clause->addWt(wt);
 
    return isAppended;
  } //appendHybridClause()


    // idx is the index into the clauses_ array
  Clause* removeClause(const int& remIdx)
  {
      //remove the clause and its corresponding MLNClauseInfo
    Clause* r = clauses_->removeItemFastDisorder(remIdx);
    MLNClauseInfo* ci = clauseInfos_->removeItemFastDisorder(remIdx);
    assert(ci->index == remIdx);

    externalClause_->removeItemFastDisorder(remIdx);

    if (clauseInfos_->size() != remIdx) //if we didn't just remove the last item
    {
        //(*clauseInfos_)[remIdx] is the last item moved to the position of the
        //removed IndexClause
      (*clauseInfos_)[remIdx]->index = remIdx;
    }

      //update MLNClauseInfo's predIdsClauseIndexes
    Array<PredIdClauseIndex*>& piciArr = ci->predIdsClauseIndexes;
    for (int i = 0; i < piciArr.size(); i++)
    {
      int predId = piciArr[i]->predId;
      int remIdx = *(piciArr[i]->clauseIndex);

      Array<IndexClause*>* indexesClauses = (*predIdToClausesMap_)[predId];
      IndexClause* ic = indexesClauses->removeItemFastDisorder(remIdx);
      assert(ic->index == remIdx && ic->clause == r);
      delete ic;

        //if pred's clause array isn't empty & we didn't just removed last item
      if (indexesClauses->size() != remIdx)
      {
          //(*predClauses)[remIdx] is the last item moved to the position of 
          //removed clause, so update its position to remIdx
        (*indexesClauses)[remIdx]->index = remIdx;        
      }

      if (indexesClauses->empty())
      {
        delete indexesClauses;
        (*predIdToClausesMap_)[predId] = NULL;
      }
    }

      //update MLNClauseInfo's formulaClauseIndexes
    Array<FormulaClauseIndexes*>& fciArr = ci->formulaClauseIndexes;
    for (int i = 0; i < fciArr.size(); i++)
    {
      int fidx = *(fciArr[i]->formulaIndex);
      int remIdx = *(fciArr[i]->clauseIndex);
      FormulaAndClauses* fac = (*formAndClausesArray_)[fidx];
      IndexClauseHashArray* indexClauses = fac->indexClauses;
      IndexClause* ic = indexClauses->removeItemFastDisorder(remIdx);
      assert(ic->index == remIdx);
      delete ic;
      if (indexClauses->size() != remIdx)
        (*indexClauses)[remIdx]->index = remIdx;
      
      if (indexClauses->empty())
      {
        fac = formAndClausesArray_->removeItemFastDisorder(fidx);
        assert(fac->index == fidx);
        delete fac;
        if (formAndClausesArray_->size() != fidx)
          (*formAndClausesArray_)[fidx]->index = fidx;
      }
    }

    delete ci;
    assert(clauses_->size() == clauseInfos_->size());
    return r;
  }//removeClause()


    // returns the removed clause or NULL if c is not in the MLN
  Clause* removeClause(const Clause* const & c)
  {
    int idx = findClauseIdx(c);
    if (idx < 0) return NULL;
    return removeClause(idx);
  }


  void removeAllClauses(Array<Clause*>* const & clauses)
  {
    while (!clauses_->empty()) 
    {
      Clause* c = removeClause(0);
      if (clauses) clauses->append(c);
      else delete c;
    }
  }


    // Clauses  are reinserted into clauses_, and the indexClauses arrays of 
    // each FormulaAndClauses in formAndClausesArrays_. Call this when the 
    // clauses' hash codes change, e.g., when constants are reordered. You have
    // to do this because hash_map does not register the hash values of its
    // contents if the values change after the contents are inserted.
  void rehashClauses()
  {
    ClauseHashArray* newha =  new ClauseHashArray;
    for (int i = 0; i < clauses_->size(); i++) 
    {
      int a = newha->append((*clauses_)[i]);
      assert(a >= 0); a = 0; //avoid compilation warning
    }
    delete clauses_;
    clauses_ = newha;
 

    for (int i = 0; i < formAndClausesArray_->size(); i++)
    {
      FormulaAndClauses* fac = (*formAndClausesArray_)[i];
      IndexClauseHashArray* newArr = new IndexClauseHashArray;
      IndexClauseHashArray* curArr = fac->indexClauses;
      for (int j = 0; j < curArr->size(); j++)
      {
        int a = newArr->append((*curArr)[j]);
        assert((*curArr)[j]->index == a); a = 0; //avoid compilation warning
      }
      curArr->clear();
      delete curArr;
      fac->indexClauses = newArr;
    }
  }


    //returns the maximum absolute value of soft weights (weights of clauses
    //that are not hard and not belonging to existentially and uniquely quant.
    //formulas).
  double getMaxAbsSoftWt()
  {
      //find the max of the (absolute) soft clause wts
    double maxSoftWt = 0;
    for (int i = 0; i < clauses_->size(); i++)
    {
      Clause* c = (*clauses_)[i];
      if (!c->isHardClause() && !isExistUniqueClause(c))
      {
        double abWt = fabs(c->getWt());
        if (abWt > maxSoftWt) maxSoftWt = abWt;
      }
    }
    return maxSoftWt;
  }


    // Returns the index into clauses_ where c is found
  int findClauseIdx(const Clause* const & c) const 
  { return clauses_->find((Clause*)c); }


  const Clause* findClause(const Clause* const & c) const 
  { 
cout << "m1" << endl;
    int i = findClauseIdx(c);
cout << "m2" << endl;
    if (i < 0) return NULL;
cout << "m3" << endl;
    return (*clauses_)[i];
  }


    //returns clause at position of clauses_ 
    //or NULL if i is not a valid index of clauses_     
  const Clause* getClause(const int& i) const
  {
    if (i < 0 || i >= clauses_->size()) return NULL;
    return (*clauses_)[i];
  }

    //returns clause at position of hybridClauses_ 
    //or NULL if i is not a valid index of hybridClauses_     
  const Clause* getHybridClause(const int& i) const
  {
    if (i < 0 || i >= hybridClauses_->size()) return NULL;
    return (*hybridClauses_)[i];
  }

  const Predicate* getNumericPred(const int& i) const
  {
    if (i < 0 || i >= numericPreds_->size()) return NULL;
    return (*numericPreds_)[i];
  }

  const double getNumericMean(const int& i) const
  {
    return (*numericMeans_)[i];
  }

  double getHybridClauseVariance(const int& i,
                                 const Domain* const & domain) const
  {
    double sum = 0.0;
    double mean = (*numericMeans_)[i];
    
    Predicate* pred = (*numericPreds_)[i];
    Array<Predicate*> predArr;
    pred->createAllGroundings(domain, predArr);

    int numPreds = predArr.size();
    for (int i = 0; i < numPreds; i++)
    {
      Predicate* newPred = predArr[i];
      double rv = domain->getDB()->getRealValue(newPred);
      sum += pow((rv - mean),2.0);
      delete newPred;
    }
    double variance = (sum / ((double)numPreds));
    return variance;
  }

  bool isExternalClause(const int& i) const
  {
    assert(0 <= i && i < clauses_->size());
    return (*externalClause_)[i];
  }

  bool isExternalClause(const Clause* const & c) const
  {
    int i = findClauseIdx(c);
    assert(i >= 0);
    return isExternalClause(i);
  }

    // returns true if the ith clause is in the CNF of an existentially
    // quantified formula
  bool isExistClause(const int& i) const
  {
    assert(0 <= i && i < clauses_->size());
    Array<FormulaClauseIndexes*>& fciArr 
      = (*clauseInfos_)[i]->formulaClauseIndexes;
    for (int i = 0; i < fciArr.size(); i++)
      if ((*formAndClausesArray_)[*(fciArr[i]->formulaIndex)]->hasExist) 
        return true;
    return false;
  }


    // returns true if clause c is in the CNF of an existentially
    // quantified formula
  bool isExistClause(const Clause* const & c) const
  {
    int i = findClauseIdx(c);
    assert(i >= 0);
    return isExistClause(i);
  }


    // returns true if the ith clause is in the CNF of an existentially and
    // uniquely quantified formula
  bool isExistUniqueClause(const int& i) const
  {
    assert(0 <= i && i < clauses_->size());
    Array<FormulaClauseIndexes*>& fciArr 
      = (*clauseInfos_)[i]->formulaClauseIndexes;
    for (int i = 0; i < fciArr.size(); i++)
      if ((*formAndClausesArray_)[*(fciArr[i]->formulaIndex)]->isExistUnique) 
        return true;
    return false;
  }


    // returns true if clause c is in the CNF of an existentially and
    // uniquely quantified formula
  bool isExistUniqueClause(const Clause* const & c) const
  {
    int i = findClauseIdx(c);
    assert(i >= 0);
    return isExistUniqueClause(i);
  }


    // returns true if the ith clause is in the CNF of a non-existentially
    // quantified formula
  bool clauseInNonExistAndNonExistUniqueFormulaCNF(const int& i) const
  {
    assert(0 <= i && i < clauses_->size());
    Array<FormulaClauseIndexes*>& fciArr 
      = (*clauseInfos_)[i]->formulaClauseIndexes;
    for (int i = 0; i < fciArr.size(); i++)
    {
      FormulaAndClauses*fnc=(*formAndClausesArray_)[*(fciArr[i]->formulaIndex)];
      if (!fnc->hasExist && !fnc->isExistUnique) return true;
    }
    return false;
  }


    // returns true if the ith clause is in the CNF of a non-existentially
    // quantified formula
  bool clauseInNonExistAndNonExistUniqueFormulaCNF(const Clause* const & c)
  const
  {
    int i = findClauseIdx(c);
    assert(i >= 0);
    return clauseInNonExistAndNonExistUniqueFormulaCNF(i);
  }


    // The list and its contents should not be modified.
  const ClauseHashArray* getClauses() const { return clauses_; }

    // The list and its contents should not be modified.
  const Array<Clause*>* getHybridClauses() const { return hybridClauses_; }

  int getSoftClauseSize()
  {
	int clauseNum = clauses_->size();
	for (int i = 0; i < clauseNum; i++)
    {
      Clause* clause = (*clauses_)[i];
      if (clause->isHardClause())
      {
        clauseNum --;   
      }
    }
    return clauseNum;
  }

  void getClauses(Array<Clause*>* const & clauses) const
  { for (int i = 0; i < clauses_->size();i++) clauses->append((*clauses_)[i]); }

  void setClauses(ClauseHashArray* const & clauses)
  { clauses_ = clauses; }

  void replaceClauses(ClauseHashArray* const & clauses)
  {
    if (clauses_)
    {
      clauses_->deleteItemsAndClear();
      delete clauses_;
    }
    clauses_ = clauses;
  }

  const MLNClauseInfo* getMLNClauseInfo(const int& i) const
  {
    if (i < 0 || i >= clauseInfos_->size()) return NULL;
    return (*clauseInfos_)[i];
  }


  int* getMLNClauseInfoIndexPtr(const int& i) const
  {
    assert(0 <= i && i < clauseInfos_->size());
    return &((*clauseInfos_)[i]->index);
  }

  const Array<MLNClauseInfo*>* getMLNClauseInfos() const
  { return clauseInfos_; }

  void setMLNClauseInfos(Array<MLNClauseInfo*>* const & clauseInfos)
  { clauseInfos_ = clauseInfos; }

  void replaceMLNClauseInfos(Array<MLNClauseInfo*>* const & clauseInfos)
  {
    if (clauseInfos_)
    {
      clauseInfos_->deleteItemsAndClear();
      delete clauseInfos_;
    }
    clauseInfos_ = clauseInfos;
  }

  
    // Caller should not modify the returned array or its contents
  const Array<Array<IndexClause*>*>* getPredIdToClausesMap() const
  { return predIdToClausesMap_; }
  
  void setPredIdToClausesMap(Array<Array<IndexClause*>*>* const & predIdToClausesMap)
  { predIdToClausesMap_ = predIdToClausesMap; }

  void replacePredIdToClausesMap(Array<Array<IndexClause*>*>* const & predIdToClausesMap)
  {
    if (predIdToClausesMap_)
    {
      for (int i = 0; i < predIdToClausesMap_->size(); i++)
        if ((*predIdToClausesMap_)[i]) 
        {
          (*predIdToClausesMap_)[i]->deleteItemsAndClear();
          delete (*predIdToClausesMap_)[i];
        }
      
      delete predIdToClausesMap_;
    }
    
    predIdToClausesMap_ = predIdToClausesMap;
  }
  
  const FormulaAndClausesArray* getFormulaAndClausesArray() const
  { return formAndClausesArray_; }
  
  void setFormulaAndClausesArray(FormulaAndClausesArray* const & formAndClausesArray)
  { formAndClausesArray_ = formAndClausesArray; }

  void replaceFormulaAndClausesArray(FormulaAndClausesArray* const & formAndClausesArray)
  {
    if (formAndClausesArray_)
    {
      formAndClausesArray_->deleteItemsAndClear();
      delete formAndClausesArray_;
    }
    formAndClausesArray_ = formAndClausesArray;
  }
  
  const Array<bool>* getExternalClause() const
  { return externalClause_; }

  void setExternalClause(Array<bool>* const & externalClause)
  {
    externalClause_ = externalClause;
  }

  void replaceExternalClause(Array<bool>* const & externalClause)
  {
    if (externalClause_)
    {
      delete externalClause_;
    }
    externalClause_ = externalClause;
  }


  void setClauseInfoPriorMeansToClauseWts()
  {
    assert(clauses_->size() == clauseInfos_->size());
    for (int i = 0; i < clauseInfos_->size(); i++)
      (*clauseInfos_)[i]->priorMean = (*clauses_)[i]->getWt();
  }


  void getClauseWts(Array<double>& wts) const
  {
    wts.clear();
    for (int i = 0; i < clauses_->size(); i++) 
      wts.append((*clauses_)[i]->getWt());
  }


  void setClauseWts(Array<double>& wts) 
  {
    assert (wts.size() == clauses_->size());
    for (int i = 0; i < clauses_->size(); i++) 
      (*clauses_)[i]->setWt(wts[i]);
  }

  
  const IndexClauseHashArray* 
  getClausesOfFormula(const string& formulaStr) const
  {
    FormulaAndClauses tmp(formulaStr, 0, false, false, false);    
    int i = formAndClausesArray_->find(&tmp);
    if (i < 0) return NULL;
    FormulaAndClauses* fnc = (*formAndClausesArray_)[i];
    return fnc->indexClauses;
  }


    //Returns true if formulaStr is in mln and its numPred is set; otherwise,
    //returns false.
  bool setFormulaNumPreds(const string& formulaStr, const int& numPreds)
  {
    FormulaAndClauses tmp(formulaStr, 0, false, false, false);    
    int i = formAndClausesArray_->find(&tmp);
    if (i < 0) return false;
    FormulaAndClauses* fnc = (*formAndClausesArray_)[i];
    fnc->numPreds = numPreds;
    return true;
  }


    //Returns true if formulaStr is in mln and its isHard is set; otherwise,
    //returns false.
  bool setFormulaIsHard(const string& formulaStr, const bool& isHard)
  {
    FormulaAndClauses tmp(formulaStr, 0, false, false, false);    
    int i = formAndClausesArray_->find(&tmp);
    if (i < 0) return false;
    FormulaAndClauses* fnc = (*formAndClausesArray_)[i];
    fnc->isHard = isHard;
    return true;
  }


    //Returns true if formulaStr is in mln and its priorMean is set; otherwise,
    //returns false.
  bool setFormulaPriorMean(const string& formulaStr, const double& priorMean)
  {
    FormulaAndClauses tmp(formulaStr, 0, false, false, false);
    int i = formAndClausesArray_->find(&tmp);
    if (i < 0) return false;
    FormulaAndClauses* fnc = (*formAndClausesArray_)[i];
    fnc->priorMean = priorMean;
    return true;
  }

    //Returns true if formulaStr is in mln and its wt is set; otherwise,
    //returns false.
  bool setFormulaWt(const string& formulaStr, const double& wt)
  {
    FormulaAndClauses tmp(formulaStr, 0, false, false, false);    
    int i = formAndClausesArray_->find(&tmp);
    if (i < 0) return false;
    FormulaAndClauses* fnc = (*formAndClausesArray_)[i];
    fnc->wt = wt;
    return true;
  }


    //Returns true if formulaStr is in mln and its isExistUnique is set; 
    //otherwise returns false.
  bool setFormulaIsExistUnique(const string& formulaStr, 
                               const bool& isExistUnique)
  {
    FormulaAndClauses tmp(formulaStr, 0, false, false, false);
    int i = formAndClausesArray_->find(&tmp);
    if (i < 0) return false;
    FormulaAndClauses* fnc = (*formAndClausesArray_)[i];
    fnc->isExistUnique = isExistUnique;
    return true;
  }


    //formulaStr must be in MLN
  double getFormulaWt(const string& formulaStr)
  {
    FormulaAndClauses tmp(formulaStr, 0, false, false, false);    
    int i = formAndClausesArray_->find(&tmp);
    FormulaAndClauses* fnc = (*formAndClausesArray_)[i];
    return fnc->wt;
  }


  const Array<IndexClause*>* getClausesContainingPred(const int& predId) const
  { 
    if (predId < predIdToClausesMap_->size())
      return (*predIdToClausesMap_)[predId];
    return NULL;
  }

  /**
   * Returns a parent formula of a clause.
   * 
   * @param clauseIdx Index of clause for which parent formula is returned
   * @param formulaIdx Index of parent formula which is returned
   * 
   * @return Parent formula as string 
   */
  const string getParentFormula(const int& clauseIdx, const int& formulaIdx) const 
  {
    FormulaClauseIndexes* fcIdxs =
      (*clauseInfos_)[clauseIdx]->formulaClauseIndexes[formulaIdx];
    return (*formAndClausesArray_)[*fcIdxs->formulaIndex]->formula;
  }

  /**
   * Returns a parent formula of a clause.
   * 
   * @param c Clause for which parent formula is returned
   * @param formulaIdx Index of parent formula which is returned
   * 
   * @return Parent formula as string 
   */
  const string getParentFormula(const Clause* const & c, const int& formulaIdx) const 
  {
    int i = findClauseIdx(c);
    if (i < 0) return NULL;
    return getParentFormula(i, formulaIdx);
  }


  void compress()
  {
    clauses_->compress();
    for (int i = 0; i < clauseInfos_->size(); i++) 
      (*clauseInfos_)[i]->compress();
    clauseInfos_->compress();
    for (int i = 0; i < formAndClausesArray_->size(); i++)
      (*formAndClausesArray_)[i]->indexClauses->compress();
    formAndClausesArray_->compress();
    for (int i = 0; i < predIdToClausesMap_->size(); i++)
      if ((*predIdToClausesMap_)[i]) (*predIdToClausesMap_)[i]->compress();
    predIdToClausesMap_->compress();
  }


    //print formulas (commented out) followed by the clauses in its CNF
  void printMLN(ostream& out, const Domain* const & domain)
  {
    int outprec = 6;
    out.precision(outprec);
    out.setf(ios_base::left, ios_base::adjustfield);

    const FormulaAndClausesArray* fncArr = formAndClausesArray_;
    
    Array<int> fncArrIdxs;
    for (int i = 0; i < fncArr->size(); i++)
      if (!(*fncArr)[i]->hasExist && !(*fncArr)[i]->isExistUnique) 
        fncArrIdxs.append(i);

    for (int i = 0; i < fncArr->size(); i++)
      if ((*fncArr)[i]->isExistUnique) 
        fncArrIdxs.append(i);
    
    for (int i = 0; i < fncArr->size(); i++)
      if ((*fncArr)[i]->hasExist) 
        fncArrIdxs.append(i);

      // for each user-specified formula
    for (int ii = 0; ii < fncArrIdxs.size(); ii++)  
    {
      int i = fncArrIdxs[ii];
      double totalWt = 0;
      IndexClauseHashArray* indexClauses = (*fncArr)[i]->indexClauses;
      for (int j = 0; j < indexClauses->size(); j++)
      {
        Clause* c = (*indexClauses)[j]->clause;
        totalWt += c->getWt()/getNumParentFormulas(c);
      }
        // output the original formula and its weight
      out.width(0); out << "// "; out.width(outprec); 
      if ((*fncArr)[i]->isHard) 
        out << (*fncArr)[i]->formula << endl;
      else
      {
        if ((*fncArr)[i]->isConjunction) totalWt = -totalWt;
        out << totalWt << "  " << (*fncArr)[i]->formula << endl;
      }
      
      if ((*fncArr)[i]->hasExist || (*fncArr)[i]->isExistUnique)
      {
        if ((*fncArr)[i]->isHard) 
          out << (*fncArr)[i]->formula << "." << endl;
        else
          out << totalWt << "  " << (*fncArr)[i]->formula <<endl;
      }
      else
      {
        if ((*fncArr)[i]->tiedClauses)
        {
          out.width(outprec); 
          if ((*fncArr)[i]->isHard) 
            out << (*fncArr)[i]->formula << "." << endl;
          else
            out << totalWt << "  " << (*fncArr)[i]->formula << endl;
        }
        else
        {
            // output clauses derived from the original formula and their
            // weights
          for (int j = 0; j < indexClauses->size(); j++)
          {
            Clause* c = (*indexClauses)[j]->clause;    
            if ((*fncArr)[i]->isHard)
            {
              c->printWithoutWtWithStrVar(out, domain);
              out << ".";
            }
            else
            {
              out.width(outprec); 
              out << c->getWt()/getNumParentFormulas(c) << "  "; 
              c->printWithoutWtWithStrVar(out, domain);
            }
            out << endl;
          }
        }
      }
      out << endl;
    }
    
    for (int i = 0; i < hybridClauses_->size(); i++)
    {
      Clause* c = (*hybridClauses_)[i];

        // output the original formula and its weight
      out.width(0); out << "// "; out.width(outprec); 
      if (!c->isHardClause()) out << c->getWt() << "  ";
      out << (*hybridFormulaStrings_)[i];
      out << endl;

      out.width(outprec);
      if (!c->isHardClause()) out << c->getWt() << "  ";
      c->printWithoutWtWithStrVar(out, domain);
      out << " * (";
      (*numericPreds_)[i]->printWithStrVar(out, domain);
      out << " = " << (*numericMeans_)[i] << ")";
      if (c->isHardClause()) out << ".";
      out << endl;
    }
  }

  
    //Print non-existential formulas (commented out) followed by the clauses in 
    //its CNF. Clause weights are divided among non-existential formulas.
  void printMLNNonExistFormulas(ostream& out, const Domain* const & domain)
  {
    int outprec = 6;
    out.precision(outprec);
    out.setf(ios_base::left, ios_base::adjustfield);

      // for each user-specified formula
    const FormulaAndClausesArray* fncArr = formAndClausesArray_;
    for (int i = 0; i < fncArr->size(); i++)
    {
      if ((*fncArr)[i]->hasExist || (*fncArr)[i]->isExistUnique) continue;
    
      double totalWt = 0;
      IndexClauseHashArray* indexClauses = (*fncArr)[i]->indexClauses;
      for (int j = 0; j < indexClauses->size(); j++)
      {
        Clause* c = (*indexClauses)[j]->clause;
        //assert(clauseOrdering->find(c) >= 0);
        totalWt 
          += c->getWt()/getNumNonExistNonExistUniqueParentFormulas(c);
      }
        // output the original formula and its weight
      out.width(0); out << "// "; out.width(outprec); 
      if ((*fncArr)[i]->isHard) 
        out << (*fncArr)[i]->formula << endl;
      else
      {
        if ((*fncArr)[i]->isConjunction) totalWt = -totalWt;
        out << totalWt << "  " << (*fncArr)[i]->formula << endl;
      }

      if ((*fncArr)[i]->tiedClauses)
      {
        out.width(outprec); 
        if ((*fncArr)[i]->isHard)
          out << (*fncArr)[i]->formula << "." << endl;
        else
          out << totalWt << "  " << (*fncArr)[i]->formula << endl;
      }
      else
      {
          // output the clauses derived from the original formula and their
          // weights
        for (int j = 0; j < indexClauses->size(); j++)
        {
          Clause* c = (*indexClauses)[j]->clause;    
          out.width(outprec); 

          if (!(*fncArr)[i]->isHard)
          {
            out << c->getWt()/getNumNonExistNonExistUniqueParentFormulas(c)
                << "  ";
          }
          c->printWithoutWtWithStrVar(out, domain);
          if ((*fncArr)[i]->isHard) out << ".";
          out << endl;
        }
      }
      out << endl;
    } 
  }


    //Print each non-existential clause as a single formula, followed by
    //existential formulas
  void printMLNClausesFormulas(ostream& out, const Domain* const & domain,
                               const bool& includeIdx)
  {
    int idx = 0;
    int* startIdx = (includeIdx) ? &idx : NULL;
    bool includeExistClauses = false;
    bool divideWtAmongExistFormulas = true;
    bool sortByLen = true;
    printClausesWithWeights(out, domain, startIdx, includeExistClauses,
                            sortByLen, divideWtAmongExistFormulas);
    printExistOrExistUniqueFormulasWithWeights(out, startIdx);
  }


    //If includeExistClauses is false, we exclude clauses that ONLY appear
    //in the CNFs of existential formulas.
  void printClausesWithWeights(ostream& out, const Domain* const & domain,
                               int* const & startIdx = NULL,
                               const bool& includeExistClauses=true,
                               const bool& sortByLen=false,
                               const bool& divWtAmongExistFormulas=false) const
  {
    Array<Clause*> ca;
    for (int i = 0; i < clauses_->size(); i++) 
    {
      if ( !includeExistClauses && 
           !clauseInNonExistAndNonExistUniqueFormulaCNF(i) ) continue;
      ca.append((*clauses_)[i]);
    }

    if (sortByLen) Clause::sortByLen(ca);
    ClauseHashArray cha;
    for (int i = 0; i < ca.size(); i++) cha.append(ca[i]);

    printClausesWithWeights(out, domain, &cha,startIdx,divWtAmongExistFormulas);
  }

  
  void printClausePriorMeans(ostream& out, const Domain* const & domain)
  {
    out.setf(ios_base::left, ios_base::adjustfield);    
    out.precision(6);
    assert(clauseInfos_->size() == clauses_->size());
    for (int i = 0; i < clauses_->size(); i++)
    {
      out << i << ":  "; out.width(14); 
      out << (*clauseInfos_)[i]->priorMean; out.width(0); out << " "; 
      (*clauses_)[i]->printWithoutWtWithStrVar(out, domain); 
      out << endl; 
    }
    out.width(0);
  }

  
  void printFormulaPriorMeans(ostream& out)
  {
    out.setf(ios_base::left, ios_base::adjustfield);
    out.precision(6);
    for (int i = 0; i < formAndClausesArray_->size(); i++)
    {
      FormulaAndClauses* fnc = (*formAndClausesArray_)[i];
      out << i << ":  "; out.width(14); 
      out << fnc->priorMean; out.width(0); 
      out << " " << fnc->formula << endl;       
    }
    out.width(0);
  }


 private:
  void appendClauseInfo(const Clause* const & clause, const int& clauseIdx, 
                        int* const& formulaIdx, int* const& clauseIdxForFormula)
  {
    assert(clauseInfos_->size() == clauseIdx);

    MLNClauseInfo* ci = new MLNClauseInfo(clauseInfos_->size());
    clauseInfos_->append(ci);

      // keep track of appended clause's positions in predIdToClausesMap_
    Array<PredIdClauseIndex*>& piciArr = ci->predIdsClauseIndexes;
    hash_set<int> seenPredIds;
    const Array<Predicate*>* preds = clause->getPredicates();
    for (int i = 0; i < preds->size(); i++)
    {
      int predId = (*preds)[i]->getId();
        //if we have already noted that clause contains pred
      if (seenPredIds.find(predId) != seenPredIds.end()) continue;
      seenPredIds.insert(predId);
      
      if (predId >= predIdToClausesMap_->size())
        predIdToClausesMap_->growToSize(predId+1, NULL);
      Array<IndexClause*>*& icArr = (*predIdToClausesMap_)[predId];
      if (icArr == NULL) icArr = new Array<IndexClause*>;
      IndexClause* ic = new IndexClause(icArr->size(), (Clause*)clause);
      icArr->append(ic);
      
      piciArr.append(new PredIdClauseIndex(predId, &(ic->index)));
    }

    updateClauseInfo(clauseIdx, formulaIdx, clauseIdxForFormula);
  }

  
  void updateClauseInfo(const int& clauseIdx, int* const & formulaIdx, 
                        int* const & clauseIdxForFormula)
  {
    MLNClauseInfo* clauseInfo = (*clauseInfos_)[clauseIdx];
      // keep track of appendedClause's position in formAndClauseArray_    
    assert(*clauseIdxForFormula >= 0);
    Array<FormulaClauseIndexes*>& fciArr = clauseInfo->formulaClauseIndexes;

    fciArr.append(new FormulaClauseIndexes(formulaIdx, clauseIdxForFormula));
  }


    //Returns -1 if c is not found in MLN
  int getNumParentFormulas(const Clause* const & c) const 
  {
    int i = findClauseIdx(c);
    if (i < 0) return -1;
    return (*clauseInfos_)[i]->formulaClauseIndexes.size();
  }


  int getNumNonExistNonExistUniqueParentFormulas(const Clause* const & c) const 
  {
    int i = findClauseIdx(c);
    if (i < 0) return -1;
    Array<FormulaClauseIndexes*>& fciArr 
      = (*clauseInfos_)[i]->formulaClauseIndexes;
    int n = 0;
    for (int i = 0; i < fciArr.size(); i++)
    {
      int fidx = *(fciArr[i]->formulaIndex);
      FormulaAndClauses* fnc = (*formAndClausesArray_)[fidx];
      if (!fnc->hasExist && !fnc->isExistUnique) n++;
    }
    return n;
  }


  int getNumExistOrExistUniqueParentFormulas(const Clause* const & c) const 
  {
    int i = findClauseIdx(c);
    if (i < 0) return -1;
    Array<FormulaClauseIndexes*>& fciArr 
      = (*clauseInfos_)[i]->formulaClauseIndexes;
    int n = 0;
    for (int i = 0; i < fciArr.size(); i++)
    {
      int fidx = *(fciArr[i]->formulaIndex);
      FormulaAndClauses* fnc = (*formAndClausesArray_)[fidx];
      if (fnc->hasExist || fnc->isExistUnique) n++;
    }
    return n;
  }


  void printClausesWithWeights(ostream& out, const Domain* const & domain,
                               const ClauseHashArray* const & clauses,
                               int* const & startIdx,
                               const bool& divideWithAmongExistFormulas) const
  {
    out.setf(ios_base::left, ios_base::adjustfield);
    int i;
    for (i = 0; i < clauses->size(); i++)
    {      
      if (startIdx) { out << (*(startIdx))++ << ":  "; out.width(14); }
      else          { out.width(10); }
      
      double wt = (*clauses)[i]->getWt();
      if (divideWithAmongExistFormulas) 
        wt /= 1+getNumExistOrExistUniqueParentFormulas((*clauses)[i]);

      out << wt; out.width(0); out << " ";
      (*clauses)[i]->printWithoutWtWithStrVar(out, domain); 
      out << endl; 
    }
    out.width(0);
  }


  void printExistOrExistUniqueFormulasWithWeights(ostream& out, 
                                                  int* const & startIdx=NULL)
  {
    int i;
    for (i = 0; i < formAndClausesArray_->size(); i++)
    {
      if (!(*formAndClausesArray_)[i]->hasExist && 
          !(*formAndClausesArray_)[i]->isExistUnique) continue;
      double totalWt = 0;
      IndexClauseHashArray* indexClauses 
        = (*formAndClausesArray_)[i]->indexClauses;
      for (int j = 0; j < indexClauses->size(); j++)
      {
        Clause* c = (*indexClauses)[j]->clause;
        totalWt += c->getWt()/getNumParentFormulas(c);
      }
      
      if (startIdx) { out << (*(startIdx))++ << ":  "; out.width(14); }
      else            { out.width(10); }
      out << totalWt; out.width(0); out << " " 
          << (*formAndClausesArray_)[i]->formula <<endl;
    }
  }


 private:
  ClauseHashArray* clauses_;
  Array<MLNClauseInfo*>* clauseInfos_;
  FormulaAndClausesArray* formAndClausesArray_;

    //predIdToClausesMap_[p] maps pred id p to an array of IndexClause.
    //The clause in each IndexClause contains a pred with id p, and the 
    //index is the index of the IndexClause in predIdToClausesMap_[p].
    //predIdToClausesMap_[p] may be NULL.
  Array<Array<IndexClause*>*>* predIdToClausesMap_;

    // externalClause_[c] indicates if clause c is external, i.e. it contains
    // constants not present in this mln (can occur when per-constant rules are
    // used with multiple dbs).
  Array<bool>* externalClause_;
  
  Array<Clause*>* hybridClauses_;
  Array<string>* hybridFormulaStrings_;
  Array<Predicate*>* numericPreds_;
  Array<double>* numericMeans_;
};

#endif
