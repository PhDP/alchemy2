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
#ifndef CLAUSE_H_JUN_26_2005
#define CLAUSE_H_JUN_26_2005

#include <ext/hash_set>
using namespace __gnu_cxx;
#include <ostream>
#include <sstream>
using namespace std;
#include <climits>
#include "domain.h"
#include "hashstring.h"
#include "hashlist.h"
#include "arraysaccessor.h"
#include "powerset.h"
#include "multdarray.h"
#include "hashint.h"
#include "database.h"
#include "predicate.h"
#include "groundpredicate.h"
#include "clausesampler.h"
#include "clausehelper.h"

const bool useInverseIndex = true;
  // 0 = no ouput, 1 = some output, 2 = full output
const int clausedebug = 0;


/******************************************************************************/
/** Need to define this type in the beginning as it is used in this file */
class HashClause;
class EqualClause;
typedef HashArray<Clause*, HashClause, EqualClause> ClauseHashArray;

class SuperClause;
class PredicateTerm;
class Variable;
class HashPredicateTerm;
class EqualPredicateTerm;

typedef hash_map<Clause*, SuperClause*, HashClause, EqualClause> 
  ClauseToSuperClauseMap;

typedef hash_map<PredicateTerm*, Variable*, HashPredicateTerm,
                 EqualPredicateTerm> PredicateTermToVariable;
/******************************************************************************/

class AddGroundClauseStruct;


  // Ensure that the dirty_ bit is consistently updated.
class Clause
{
  friend class ClauseSampler;

 public: 
  Clause() 
    : wt_(0), predicates_(new Array<Predicate*>), intArrRep_(NULL),
      hashCode_(0), dirty_(true), isHardClause_(false), locked_(false),
      varIdToVarsGroundedType_(NULL), auxClauseData_(NULL), staticWt_(false),
      action_(false), util_(0.0) {}

  Clause(const double& wt) 
    : wt_(wt), predicates_(new Array<Predicate*>), intArrRep_(NULL),
      hashCode_(0), dirty_(true), isHardClause_(false), locked_(false),
      varIdToVarsGroundedType_(NULL), auxClauseData_(NULL), staticWt_(false),
      action_(false), util_(0.0) {}


  Clause(const Clause& c)
  {
    wt_ = c.wt_;
    locked_ = c.locked_;
    predicates_ = new Array<Predicate*>;
    Array<Predicate*>* cpredicates = c.predicates_;
    for (int i = 0; i < cpredicates->size(); i++)
    {
      Predicate* p = (*cpredicates)[i];
      predicates_->append(new Predicate(*p, this));      
    }
    dirty_ = c.dirty_;

    if (!dirty_) { assert(noDirtyPredicates()); }

    if (c.intArrRep_)  intArrRep_ = new Array<int>(*(c.intArrRep_));
    else               intArrRep_ = NULL;

    hashCode_ = c.hashCode_;
    
    isHardClause_ = c.isHardClause_;

    staticWt_ = c.staticWt_;
    
    action_ = c.action_;
    
    util_ = c.util_;

    if (c.varIdToVarsGroundedType_)
    {
      varIdToVarsGroundedType_ = new Array<VarsGroundedType*>;
      for (int i = 0; i < c.varIdToVarsGroundedType_->size(); i++)
      {
        VarsGroundedType* vgt = (*(c.varIdToVarsGroundedType_))[i];
        varIdToVarsGroundedType_->append(new VarsGroundedType(*vgt));
      }
    }
    else
      varIdToVarsGroundedType_ = NULL;

    if (c.auxClauseData_)
    {
      auxClauseData_ = new AuxClauseData(c.auxClauseData_->gain,
                                         c.auxClauseData_->op,
                                         c.auxClauseData_->removedClauseIdx,
                                         c.auxClauseData_->hasBeenExpanded,
                                         c.auxClauseData_->lastStepExpanded,
                                       c.auxClauseData_->lastStepOverMinWeight);
      if (c.auxClauseData_->constTermPtrs) trackConstants();
      if (c.auxClauseData_->cache)
      {
        Array<Array<Array<CacheCount*>*>*>* cache, * otherCache;
        cache  = new Array<Array<Array<CacheCount*>*>*>;
        otherCache = c.auxClauseData_->cache;

        cache->growToSize(otherCache->size(),NULL);

        for (int i = 0; i < otherCache->size(); i++)
        {
          (*cache)[i] = new Array<Array<CacheCount*>*>;
          (*cache)[i]->growToSize((*otherCache)[i]->size(), NULL);
          for (int j = 0; j < (*otherCache)[i]->size(); j++)
          {
            Array<CacheCount*>* ccArr = (*(*otherCache)[i])[j];
            if (ccArr == NULL) continue;
            (*(*cache)[i])[j] = new Array<CacheCount*>;
            (*(*cache)[i])[j]->growToSize(ccArr->size());
            for (int k = 0; k < ccArr->size(); k++)
            {
              (*(*(*cache)[i])[j])[k] = new CacheCount((*ccArr)[k]->g,
                                                       (*ccArr)[k]->c,
                                                       (*ccArr)[k]->cnt); 
            }
          }
        }
        auxClauseData_->cache = cache;
      }
      
      auxClauseData_->prevClauseStr = c.auxClauseData_->prevClauseStr;
      auxClauseData_->addedPredStr = c.auxClauseData_->addedPredStr;
      auxClauseData_->removedPredIdx = c.auxClauseData_->removedPredIdx;
    }
    else
      auxClauseData_ = NULL;
  }


  ~Clause() 
  {
    for (int i = 0; i < predicates_->size(); i++)
      delete (*predicates_)[i];
    delete predicates_;
    if (intArrRep_) delete intArrRep_;
    if (varIdToVarsGroundedType_) deleteVarIdToVarsGroundedType();
    if (auxClauseData_) delete auxClauseData_;
  }

  
    //returns approximate size in MB, mainly due to cache in auxClauseData_
  double sizeMB() const
  {
    double sizeMB = (fixedSizeB_ + intArrRep_->size()*sizeof(int) + 
                     predicates_->size()*sizeof(Predicate*)) /1000000.0;
    for (int i = 0; i < predicates_->size(); i++)
      sizeMB += (*predicates_)[i]->sizeMB();
    if (auxClauseData_ != NULL) sizeMB += auxClauseData_->sizeMB(); 
    return sizeMB;
  }


  static void computeFixedSizeB()
  {
    fixedSizeB_ = sizeof(Clause) + sizeof(Array<Predicate*>) +
                  sizeof(Array<int>);
    //+ sizeof(Array<VarsGroundedType*>); // this is transient
  }


  void compress()
  {
    predicates_->compress();
    for (int i = 0; i < predicates_->size(); i++) (*predicates_)[i]->compress();
    if (intArrRep_) intArrRep_->compress();
    if (varIdToVarsGroundedType_) varIdToVarsGroundedType_->compress();
    if (auxClauseData_) auxClauseData_->compress();
  }


  bool same(Clause* const & c)
  {
    if (this == c)  return true;
    const Array<int>* cArr  = c->getIntArrRep();
    const Array<int>* myArr = getIntArrRep();
    if (myArr->size() != cArr->size()) return false;
    const int* cItems  = c->getIntArrRep()->getItems();
    const int* myItems = getIntArrRep()->getItems();
    return (memcmp(myItems, cItems, myArr->size()*sizeof(int))==0);
  }

  size_t hashCode() 
  {
    if (dirty_) computeAndStoreIntArrRep();
    return hashCode_;
  }


  int getNumPredicates() const { return predicates_->size(); }

  double getWt() const { return wt_; }

  const double* getWtPtr() const { return &wt_; }

    // not setting dirty bit because it does not affect the clause
  void setWt(const double& wt) { if (!locked_) wt_ = wt; }  
  void addWt(const double& wt) { if (!locked_) wt_ += wt; }
 
  double getUtil() const { return util_; }

    // not setting dirty bit because it does not affect the clause
  void setUtil(const double& util) { util_ = util; }  
  void addUtil(const double& util) { util_ += util; }

  void setDirty() { dirty_ = true; }
  bool isDirty() const { return dirty_; }

    // Caller should not delete returned Predicate*.
  Predicate* getPredicate(const int& idx) const { return (*predicates_)[idx]; }
  
    // Caller should not delete returned array nor modify its contents.
  const Array<Predicate*>* getPredicates() const { return predicates_; }


  bool containsPredicate(const Predicate* const & pred) const
  {
    for (int i = 0; i < predicates_->size(); i++)
      if ((*predicates_)[i]->same((Predicate*)pred)) return true;
    return false;
  }

  
  int getNumVariables() const
  {
    hash_set<int> intset;
    for (int i = 0; i < predicates_->size(); i++)
      for (int j = 0; j < (*predicates_)[i]->getNumTerms(); j++)
      {
        if ((*predicates_)[i]->getTerm(j)->getType() == Term::VARIABLE)
          intset.insert((*predicates_)[i]->getTerm(j)->getId());
      }
    return intset.size();
  }


  int getNumVariablesAssumeCanonicalized() const
  {
    int minVarId = 0;
    for (int i = 0; i < predicates_->size(); i++)
      for (int j = 0; j < (*predicates_)[i]->getNumTerms(); j++)
      {
        if ((*predicates_)[i]->getTerm(j)->getType() == Term::VARIABLE &&
            (*predicates_)[i]->getTerm(j)->getId() < minVarId)
          minVarId = (*predicates_)[i]->getTerm(j)->getId();

      }
    return -minVarId;
  }


  bool isHardClause() const { return isHardClause_; }
  void setIsHardClause(const bool& b) { isHardClause_ = b; }

  bool isStaticWt() const { return staticWt_; }
  void setStaticWt(const bool& b) { staticWt_ = b; }

  bool isAction() const { return action_; }
  void setAction(const bool& b) { action_ = b; }

  bool isLocked() const { return locked_; }
  void lock()   { locked_ = true; }
  void unlock() { locked_ = false; }


    // p is stored in MLN and the caller of this function should not delete it.
  void appendPredicate(Predicate* const& p) {predicates_->append(p);setDirty();}

 /**
  * Removes a predicate from this clause. After removing a predicate, the clause
  * is no longer canonicalized, so canonicalize() must be called.
  * 
  * @param i Index of predicate to be removed.
  */
  Predicate* removePredicate(const int& i) 
  { 
    if (0 <= i && i < predicates_->size()) 
    { setDirty(); return predicates_->removeItemFastDisorder(i); }
    return NULL;
  }


 /**
  * returns true if this clause contains redundant predicates.
  */
  bool hasRedundantPredicates()
  {
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* ip = (*predicates_)[i];
      for (int j = i+1; j < predicates_->size(); j++)
      {
        Predicate* jp = (*predicates_)[j];
        if (jp->same(ip) && jp->getSense() == ip->getSense()) return true;
      }
    }
    return false;
  }


    //returns true if redundant predicates were removed
  bool removeRedundantPredicates()
  {
    bool changed = false;
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* ip = (*predicates_)[i];
      for (int j = i+1; j < predicates_->size(); j++)
      {
        Predicate* jp = (*predicates_)[j];
        if (jp->same(ip) && jp->getSense() == ip->getSense())
        {
          predicates_->removeItemFastDisorder(j);
          changed = true;
          j--;
        }
      }
    }
    return changed;
  }

  
  void removeRedundantPredicatesAndCanonicalize()
  { if (removeRedundantPredicates()) canonicalize(); }

  
  //There is an additional (optional) parameter in which the mapping
  //from old variable ids to the new canonicalized variable ids can be stored
  void canonicalize(Array<int> * const & oldVarIdToNewVarId=NULL)
  {
      // rename the variables in decreasing order (negative numbers) 
    Array<VarsGroundedType*>* vgtArr = new Array<VarsGroundedType*>;
    createVarIdToVarsGroundedType(vgtArr);

    sortPredicatesByIdAndSenseAndTerms(0, predicates_->size()-1, vgtArr);

    IntHashArray varAppearOrder;
    getVarOrder(varAppearOrder);

    int newVarId = 0;
    for (int i = 0; i < varAppearOrder.size(); i++)
    {
      int varId = varAppearOrder[i];
      ++newVarId;
      Array<Term*>& vars = (*vgtArr)[varId]->vars;
      for (int j = 0; j < vars.size(); j++)  vars[j]->setId(-newVarId);
      
	  if(oldVarIdToNewVarId)
      {
        (*oldVarIdToNewVarId)[varId] = newVarId;
      }
    }
	
    if(oldVarIdToNewVarId)
      (*oldVarIdToNewVarId)[0] = 0;
    assert(newVarId == varAppearOrder.size());

    int i = 0, j = 0;
    while (i < predicates_->size())
    {
      while (++j < predicates_->size() 
             && (*predicates_)[i]->getId() == (*predicates_)[j]->getId() 
             && (*predicates_)[i]->getSense() == (*predicates_)[j]->getSense());
      if (i < j-1) sortPredicatesByTerms(i, j-1);
      i = j;
    }
    
    deleteVarIdToVarsGroundedType(vgtArr);
    setDirty();
  }


  static bool moveTermsFromUnseenToSeen(Array<Term*>* const & terms,
                                 PredicateHashArray& unseenPreds,
                                 Array<Predicate*>& seenPreds)
  {
    for (int j = 0; j < terms->size(); j++)
    {
      bool parIsPred;
      Predicate* parent = (Predicate*) (*terms)[j]->getParent(parIsPred);      
      assert(parIsPred);
      int a;
      if ((a=unseenPreds.find(parent)) >= 0)
      {
        Predicate* par = unseenPreds.removeItemFastDisorder(a);
        if (unseenPreds.empty()) return true;
        assert(par == parent);
        seenPreds.append(par);
      }
    }
    return false;
  }


    //Returns true if there is a path of shared variables between any two preds
  bool checkPredsAreConnected()
  {
    if (predicates_->size() <= 1) return true;
    Array<Array<Term*>*>* varIdToTerms = new Array<Array<Term*>*>;
    hash_map<int,Array<Term*>*>* constIdToTerms=new hash_map<int,Array<Term*>*>;
    createVarConstIdToTerms(varIdToTerms, constIdToTerms);
    PredicateHashArray unseenPreds;
    Array<Predicate*> seenPreds;
    hash_set<int> seenIds;

    seenPreds.append((*predicates_)[0]);
    for (int i = 1; i < predicates_->size(); i++) 
      unseenPreds.append((*predicates_)[i]);

    while (!seenPreds.empty())
    {
      Predicate* curPred = seenPreds.removeLastItem();
      for (int i = 0; i < curPred->getNumTerms(); i++) 
      {
        const Term* t = curPred->getTerm(i);
        int id = t->getId();
        if (seenIds.find(id) != seenIds.end()) continue;
        seenIds.insert(id);

        if (t->getType() == Term::VARIABLE)
        {
          Array<Term*>* terms = (*varIdToTerms)[-id];
          for (int j = 0; j < terms->size(); j++)
          {
            bool parIsPred;
            Predicate* parent = (Predicate*) (*terms)[j]->getParent(parIsPred);
            assert(parIsPred);
            int a;
            if ((a=unseenPreds.find(parent)) >= 0)
            {
              Predicate* par = unseenPreds.removeItemFastDisorder(a);
              if (unseenPreds.empty()) 
              {
                deleteVarConstIdToTerms(varIdToTerms, constIdToTerms);
                return true;
              }
              seenPreds.append(par);

              //commented out: not true if there a duplicate preds in clause
              //               e.g. !isMale(p) v isMale(p)
              //assert(par == parent);
            }
          }
        }
        else
        if (t->getType() == Term::CONSTANT)
        {
          Array<Term*>* terms = (*constIdToTerms)[id];
          for (int j = 0; j < terms->size(); j++)
          {
            bool parIsPred;
            Predicate* parent = (Predicate*) (*terms)[j]->getParent(parIsPred);
            assert(parIsPred);
            int a;
            if ((a=unseenPreds.find(parent)) >= 0)
            {
              Predicate* par = unseenPreds.removeItemFastDisorder(a);
              if (unseenPreds.empty()) 
              {
                deleteVarConstIdToTerms(varIdToTerms, constIdToTerms);
                return true;
              }
              seenPreds.append(par);
                //commented out: not true if there a duplicate preds in clause
                //               e.g. !isMale(p) v isMale(p)
              //assert(par == parent);
            }
          }
        }
      } // for each term of curPred
    } // while (!seenPreds.empty())

    assert(!unseenPreds.empty());
    deleteVarConstIdToTerms(varIdToTerms, constIdToTerms);
    return false;
  }


  void canonicalizeWithoutVariables()
  {
    sortPredicatesByIdAndSenseAndTerms(0, predicates_->size()-1, NULL);
    int i = 0, j = 0;
    while (i < predicates_->size())
    {
      while (++j < predicates_->size() 
             && (*predicates_)[i]->getId() == (*predicates_)[j]->getId() 
             && (*predicates_)[i]->getSense() == (*predicates_)[j]->getSense());
      if (i < j-1) sortPredicatesByTerms(i, j-1);
      i = j;
    }
    setDirty();
  }

  
  AuxClauseData* getAuxClauseData() const { return auxClauseData_; }


  void setAuxClauseData(AuxClauseData* const& acd ) 
  { 
    if (auxClauseData_ == acd) return;
    if (auxClauseData_) delete auxClauseData_;
    auxClauseData_ = acd;
  }


  void newAuxClauseData() 
  {if (auxClauseData_) delete auxClauseData_; auxClauseData_=new AuxClauseData;}


  static void setClauseSampler(ClauseSampler* const & cs) 
  { if (clauseSampler_) delete clauseSampler_; clauseSampler_ = cs; }


  static const ClauseSampler* getClauseSampler() { return clauseSampler_; }
  

  bool containsConstants() const
  {
    return (auxClauseData_ && auxClauseData_->constTermPtrs && 
            auxClauseData_->constTermPtrs->size() > 0);
  }


    //auxClauseData_ must have been set to a valid AuxClauseData object
  void trackConstants()
  {
    if (auxClauseData_ == NULL) auxClauseData_ = new AuxClauseData; 
    Array<Term*>*& constTermPtrs = auxClauseData_->constTermPtrs;
    if (constTermPtrs) constTermPtrs->clear();
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* p = (*predicates_)[i];
      for (int j = 0; j < p->getNumTerms(); j++)
      {
        const Term* t = p->getTerm(j);
        if (t->getType() == Term::CONSTANT)
        {
          if (constTermPtrs == NULL) constTermPtrs = new Array<Term*>;
          constTermPtrs->append((Term*)t);
        }
      }
    }
  }


    //auxClauseData_ must have been set to a valid AuxClauseData object  
  void newCache(const int& numDomains, const int& numPreds)
  {
    assert(auxClauseData_);
    auxClauseData_->deleteCache();
    auxClauseData_->cache = new Array<Array<Array<CacheCount*>*>*>;
    auxClauseData_->cache->growToSize(numDomains, NULL);
    for (int i = 0; i < numDomains; i++)
    {
      (*auxClauseData_->cache)[i] = new Array<Array<CacheCount*>*>;
      (*auxClauseData_->cache)[i]->growToSize(numPreds, NULL);
    }
  }

  
  void translateConstants(const Domain* const & orig, const Domain* const& nnew)
  {
    if (auxClauseData_ == NULL || auxClauseData_->constTermPtrs == NULL) return;

    Array<Term*>* constTermPtrs = auxClauseData_->constTermPtrs;
    for (int i = 0; i < constTermPtrs->size(); i++)
    {
      Term* t = (*constTermPtrs)[i];
      int newId = nnew->getConstantId(orig->getConstantName(t->getId()));
      t->setId(newId);
      if (newId < 0)
      {
        cout << "ERROR: in Clause::translateConstants(). Failed to translate "
             <<orig->getConstantName(t->getId())<<" from one domain to another."
             << "Check that the constant exists in all domains." << endl;
        exit(-1);                                      
      }
    }
  }

    //Returns a map from typeId to variable ids. ReturnedArray[typeId] is NULL
    //if there are no variables for the corresponding typeId.
    //Caller is responsible for deleting the return Array and its contents.
  Array<Array<int>*>* getTypeIdToVarIdsMapAndSmallestVarId(int& smallestVarId)
  {
    smallestVarId = 0;
    hash_set<int> intSet;
    Array<Array<int>*>* arr = new Array<Array<int>*>;
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* pred = (*predicates_)[i];
      for (int j = 0; j < pred->getNumTerms(); j++)
      {
        const Term* t = pred->getTerm(j);
        if (t->getType() == Term::VARIABLE)
        {
          int varId = t->getId();
          if (intSet.find(varId) != intSet.end()) continue;
          intSet.insert(varId);

          int typeId = pred->getTermTypeAsInt(j); 
          if (typeId >= arr->size()) arr->growToSize(typeId+1, NULL);
          if ((*arr)[typeId] == NULL) (*arr)[typeId] = new Array<int>;
          (*arr)[typeId]->append(varId); 
          
          if (varId < smallestVarId) smallestVarId = varId;
        }
      }
    }
    return arr;
  }


  static void sortByLen(Array<Clause*>& ca) { sortByLen(ca, 0, ca.size()-1); } 


  static string getOpAsString(const int& op)
  {
    if (op == OP_ADD)             return "OP_ADD";
    if (op == OP_REPLACE_ADDPRED) return "OP_REPLACE_ADDPRED";
    if (op == OP_REPLACE_REMPRED) return "OP_REPLACE_REMPRED";
    if (op == OP_REMOVE)          return "OP_REMOVE";
    if (op == OP_REPLACE)         return "OP_REPLACE";
    if (op == OP_NONE)            return "OP_NONE";
    return "";
  }


  ostream& print(ostream& out, const Domain* const& domain, 
                 const bool& withWt, const bool& asInt, 
                 const bool& withStrVar) const
  {
    if (withWt) out << wt_ << " ";

    Array<Predicate*> eqPreds;
    Array<Predicate*> internalPreds;
    for (int i = 0; i < predicates_->size(); i++)
    {
      if ((*predicates_)[i]->isEqualPred()) 
        eqPreds.append((*predicates_)[i]);
      else
      if ((*predicates_)[i]->isInternalPred()) 
        internalPreds.append((*predicates_)[i]);
      else
      {
        if (asInt)           (*predicates_)[i]->printAsInt(out);
        else if (withStrVar) (*predicates_)[i]->printWithStrVar(out, domain);
        else                 (*predicates_)[i]->print(out,domain);
        if (i < predicates_->size()-1 || !eqPreds.empty() ||
        	!internalPreds.empty()) out << " v ";
      }
    }

    for (int i = 0; i < eqPreds.size(); i++)
    {
      if (asInt)           eqPreds[i]->printAsInt(out);
      else if (withStrVar) eqPreds[i]->printWithStrVar(out,domain);
      else                 eqPreds[i]->print(out,domain);
      out << ((i != eqPreds.size()-1 || !internalPreds.empty())?" v ":"");      
    }

    for (int i = 0; i < internalPreds.size(); i++)
    {
      if (asInt)           internalPreds[i]->printAsInt(out);
      else if (withStrVar) internalPreds[i]->printWithStrVar(out,domain);
      else                 internalPreds[i]->print(out,domain);
      out << ((i!=internalPreds.size()-1)?" v ":"");      
    }

    return out;
  }


  ostream& printAsInt(ostream& out) const
  { return print(out, NULL, true, true, false); }

  ostream& printWithoutWt(ostream& out, const Domain* const & domain) const
  { return print(out, domain, false, false, false); }

  ostream& 
  printWithoutWtWithStrVar(ostream& out, const Domain* const & domain) const
  { return print(out, domain, false, false, true); }

  ostream& printWithWtAndStrVar(ostream& out, const Domain* const& domain) const
  { return print(out, domain, true, false, true); }

  ostream& print(ostream& out, const Domain* const& domain) const
  { return print(out, domain, true, false, false); }
    
  ostream& printWithoutWtWithStrVarAndPeriod(ostream& out, 
                                              const Domain* const& domain) const
  {
    printWithoutWtWithStrVar(out, domain);
    if (isHardClause_) out << ".";
    return out;
  }


 private:
  static int comparePredicatesByIdAndSenseAndTerms(
                             const Predicate* const & l, 
                                             const Predicate* const & r,
                                 const Array<VarsGroundedType*>* const & vgtArr)
  {
    if (l->getId() < r->getId()) return -1;
    if (l->getId() > r->getId()) return 1;
    if (l->getSense() && !(r->getSense())) return -1;
    if (!(l->getSense()) && r->getSense()) return 1;
    assert(l->getNumTerms() == r->getNumTerms());
    for (int i = 0; i < l->getNumTerms(); i++)
    {
      int lid = l->getTerm(i)->getId();
      int rid = r->getTerm(i)->getId();
      if (vgtArr && lid < 0 && rid < 0)
      {
        bool lbound = (*vgtArr)[-lid]->vars.size() > 1;
        bool rbound = (*vgtArr)[-rid]->vars.size() > 1;
        if (lbound && !rbound) return -1;
        if (!lbound && rbound) return 1;
      }
      if (lid > rid) return -1;
      if (lid < rid) return 1;
    }
    return 0;
  }


  void sortPredicatesByIdAndSenseAndTerms(const int& l, const int& r,
                                 const Array<VarsGroundedType*>* const & vgtArr)
  {
    if (l >= r) return;
    Predicate* tmp = (*predicates_)[l];
    (*predicates_)[l] = (*predicates_)[(l+r)/2];
    (*predicates_)[(l+r)/2] = tmp;

    int last = l;
    for (int i = l+1; i <= r; i++)
      if (comparePredicatesByIdAndSenseAndTerms((*predicates_)[i],
                                                (*predicates_)[l], vgtArr) < 0)
      {
        ++last;
        tmp = (*predicates_)[last];
        (*predicates_)[last] = (*predicates_)[i];
        (*predicates_)[i] = tmp;
      }
    
    tmp = (*predicates_)[l];
    (*predicates_)[l] = (*predicates_)[last];
    (*predicates_)[last] = tmp;
    sortPredicatesByIdAndSenseAndTerms(l, last-1, vgtArr);
    sortPredicatesByIdAndSenseAndTerms(last+1, r, vgtArr);  
  }


    //l and r must have the same id and sense
  static int comparePredicatesByTerms(const Predicate* const & l,
                                      const Predicate* const & r)
  {
    assert(l->getId() == r->getId());
    assert(l->getSense() == r->getSense());
    for (int i = 0; i < l->getNumTerms(); i++)
    {
      const Term* lTerm =l->getTerm(i);
      const Term* rTerm =r->getTerm(i);
      int lTermType = lTerm->getType();
      int rTermType = rTerm->getType();
      if (lTermType == Term::VARIABLE && rTermType == Term::VARIABLE)
      {
        if (lTerm->getId() > rTerm->getId()) return -1;
        if (lTerm->getId() < rTerm->getId()) return 1;
      }
      else
      if (lTermType == Term::CONSTANT && rTermType == Term::CONSTANT)
      {
        if (lTerm->getId() < rTerm->getId()) return -1;
        if (lTerm->getId() > rTerm->getId()) return 1;
      }
      else
      if (lTermType == Term::VARIABLE && rTermType == Term::CONSTANT)
        return -1;
      else
      if (lTermType == Term::CONSTANT && rTermType == Term::VARIABLE)
        return 1;
      else
      {
        assert(false);
      }
    }
    return 0;
  }


    //The lth to rth predicates must have the same id and sense 
  void sortPredicatesByTerms(const int& l, const int& r)
  {
    if (l >= r) return;
    Predicate* tmp = (*predicates_)[l];
    (*predicates_)[l] = (*predicates_)[(l+r)/2];
    (*predicates_)[(l+r)/2] = tmp;
    
    int last = l;
    for (int i = l+1; i <= r; i++)
      if (comparePredicatesByTerms((*predicates_)[i],(*predicates_)[l]) < 0)
      {
        ++last;
        tmp = (*predicates_)[last];
        (*predicates_)[last] = (*predicates_)[i];
        (*predicates_)[i] = tmp;
      }
    
    tmp = (*predicates_)[l];
    (*predicates_)[l] = (*predicates_)[last];
    (*predicates_)[last] = tmp;
    sortPredicatesByTerms(l, last-1);
    sortPredicatesByTerms(last+1, r);  
  }


  bool noDirtyPredicates() const
  {
    for (int i = 0; i < predicates_->size(); i++)
      if ((*predicates_)[i]->isDirty()) return false;
    return true;
  }

  const Array<int>* getIntArrRep() 
  { if (dirty_) computeAndStoreIntArrRep(); return intArrRep_; }


  void computeAndStoreIntArrRep()
  {
    dirty_ = false;
    if (intArrRep_ == NULL) intArrRep_ = new Array<int>;
    else                    intArrRep_->clear();

    int numPred = predicates_->size();
    for (int i = 0; i < numPred; i++)
    {
      if ((*predicates_)[i]->getSense()) intArrRep_->append(1);
      else                               intArrRep_->append(0);
      (*predicates_)[i]->appendIntArrRep(*intArrRep_);      
    }
    hashCode_ = Hash::hash(*intArrRep_);
  }

  
  ///////////// functions for counting number of true groundings ////////////
 public:
  void addUnknownClauses(const Domain* const & domain, 
                         const Database* const& db, const int& gndPredIdx, 
                         const GroundPredicate* const & groundPred,
                         const AddGroundClauseStruct* const & agcs)
  {
    if (clausedebug >= 2)
    {
      cout << "In Clause::addUnknownClauses()" << endl;
    }
    createVarIdToVarsGroundedType(domain);
    if (gndPredIdx >= 0) groundPredVars(gndPredIdx, groundPred);
    countNumTrueGroundings(domain, db, true, false, gndPredIdx, groundPred, 
                           NULL, NULL, NULL, NULL, agcs, NULL, NULL, NULL,
                           true, false);
    restoreVars();
    deleteVarIdToVarsGroundedType();
  }


  void getUnknownClauses(const Domain* const & domain, 
                         const Database* const& db, const int& gndPredIdx, 
                         const GroundPredicate* const & groundPred,
                         const Predicate* const & gndPred,
                         Array<GroundClause*>* const& unknownGndClauses,
                         Array<Clause*>* const& unknownClauses)
  {
    createVarIdToVarsGroundedType(domain);
    if (gndPredIdx >= 0) 
    {
      if (groundPred)          groundPredVars(gndPredIdx, groundPred);
      else  { assert(gndPred); groundPredVars(gndPredIdx, (Predicate*)gndPred);}
    }
    countNumTrueGroundings(domain, db, true, false, gndPredIdx, groundPred, 
                           gndPred, unknownGndClauses, unknownClauses, NULL,
                           NULL, NULL, NULL, NULL, true, false);
    restoreVars();
    deleteVarIdToVarsGroundedType();
  }


  bool isSatisfiable(const Domain* const & domain, const Database* const& db,
                     const bool& hasUnknownPreds)
  {
    createVarIdToVarsGroundedType(domain);

    double i = countNumTrueGroundings(domain, db, hasUnknownPreds ,true,
                                      -1, NULL, NULL, NULL, NULL, NULL, NULL,
                                      NULL, NULL, NULL, true, false);
    restoreVars();
    deleteVarIdToVarsGroundedType();
    return (i > 0);
  }


    //count the number of groundings of clause variables
  double getNumGroundings(const Domain* const & domain)
  {
    createVarIdToVarsGroundedType(domain);
    double n = countNumGroundings();
    deleteVarIdToVarsGroundedType();
    return n;
  }


  double getNumTrueGroundings(const Domain* const & domain,
                              const Database* const & db,
                              const bool& hasUnknownPreds)
  {
    createVarIdToVarsGroundedType(domain);
    double n = countNumTrueGroundings(domain, db, hasUnknownPreds, false, 
                                      -1, NULL, NULL, NULL, NULL, NULL, NULL,
                                      NULL, NULL, NULL, true, false);
    restoreVars();
    deleteVarIdToVarsGroundedType();
    return n;
  }


  void getNumTrueUnknownGroundings(const Domain* const & domain,
                                   const Database* const & db,
                                   const bool& hasUnknownPreds,
                                   double& numTrue, double& numUnknown)
  {
    createVarIdToVarsGroundedType(domain);
    numTrue = countNumTrueGroundings(domain, db, hasUnknownPreds, false, 
                                     -1, NULL, NULL, NULL, NULL, &numUnknown,
                                     NULL, NULL, NULL, NULL, true, false);
    restoreVars();
    deleteVarIdToVarsGroundedType();
  }


  double getNumUnknownGroundings(const Domain* const & domain,
                                 const Database* const & db,
                                 const bool& hasUnknownPreds)
  {
    double numTrue, numUnknown;
    getNumTrueUnknownGroundings(domain, db, hasUnknownPreds,numTrue,numUnknown);
    return numUnknown;
  }
  
  void getNumTrueFalseUnknownGroundings(const Domain* const & domain,
                                        const Database* const & db,
                                        const bool& hasUnknownPreds,
                                        double& numTrue, double& numFalse,
                                        double& numUnknown)
  {
    getNumTrueUnknownGroundings(domain, db, hasUnknownPreds,numTrue,numUnknown);
    numFalse = getNumGroundings(domain) - numTrue - numUnknown;
    assert(numTrue >= 0);
    assert(numUnknown >= 0);
    assert(numFalse >= 0);
  }


    //Count the difference between the number of true groundings of the clause
    //when gndPred is held to the opposite of its actual value and when held to
    //its actual value. 
  double countDiffNumTrueGroundings(Predicate* const & gndPred, 
                                    const Domain* const & domain,
                                    Database* const & db,
                                    const bool& hasUnknownPreds,
                                    const bool& sampleClauses,
                                    const int& combo,
                                    Array<Clause*>* const & tiedClauses)
  {
    assert(gndPred->isGrounded());

      //store the indexes of the predicates that can be grounded as gndPred
    Array<int> gndPredIndexes;

    for (int i = 0; i < predicates_->size(); i++)
    {
      if ((*predicates_)[i]->canBeGroundedAs(gndPred)) gndPredIndexes.append(i);
    }
      //create mapping of variable ids (e.g. -1) to variable addresses,
      //note whether they have been grounded, and store their types
    createVarIdToVarsGroundedType(domain); 

    TruthValue actual = db->getValue(gndPred);
    assert(actual == TRUE || actual == FALSE);
    TruthValue opp = (actual == TRUE) ? FALSE : TRUE;
    bool flipped = false;

      //count # true groundings when gndPred is held to its actual value
    double numTrueGndActual = 
      countNumTrueGroundingsForAllComb(gndPredIndexes, gndPred, actual, flipped,
                                       domain, hasUnknownPreds, sampleClauses,
                                       tiedClauses);
    //cout << "numTrueGndActual = " << numTrueGndActual << endl;
    
      //count # true groundings when gndPred is held to opposite value
    double numTrueGndOpp = 0.0;
    
    int blockIdx = domain->getBlock(gndPred);
    if (blockIdx >= 0)
    {
        // Pred in block: We have to look at combination c
        // of other preds in block
      assert(combo < domain->getBlockSize(blockIdx));
      
      const Predicate* oldTruePred = domain->getTruePredInBlock(blockIdx);

      int oldTrueOne = 0;
      if (oldTruePred)
        oldTrueOne = domain->getIndexOfPredInBlock(oldTruePred, blockIdx);
//cout << "oldTrueOne " << oldTrueOne << endl;
      assert(oldTrueOne > -1);

      int newTrueOne = (oldTrueOne <= combo) ? combo + 1 : combo;
//cout << "newTrueOne " << newTrueOne << endl;
      const Predicate* newTruePred =
        domain->getPredInBlock(newTrueOne, blockIdx);

      bool isOldTrue =
        oldTruePred ? gndPred->same((Predicate*)oldTruePred) : false;
      bool isNewTrue = gndPred->same((Predicate*)newTruePred);      
      TruthValue newTV =
        ((!isOldTrue && isNewTrue) || (isOldTrue && !isNewTrue)) ? opp : actual;
      
      if (oldTruePred) assert(db->getValue(oldTruePred) == TRUE);
      assert(db->getValue(newTruePred) == FALSE);
      
      if (oldTruePred) db->setValue(oldTruePred, FALSE);
      db->setValue(newTruePred, TRUE);
      
      numTrueGndOpp +=
        countNumTrueGroundingsForAllComb(gndPredIndexes, gndPred, newTV, 
                                         newTV != actual, domain,
                                         hasUnknownPreds, sampleClauses,
                                         tiedClauses);

      //numTrueGndOpp +=
      //  countNumTrueGroundingsForAllComb(gndPredIndexes, (*block)[oldTrueOne],
      //                                   opp, flipped, domain, hasUnknownPreds,
      //                                   sampleClauses, tiedClauses);
      //numTrueGndOpp +=
      //  countNumTrueGroundingsForAllComb(gndPredIndexes, (*block)[newTrueOne],
      //                                   opp, flipped, domain, hasUnknownPreds,
      //                                   sampleClauses, tiedClauses);
      
      if (oldTruePred) db->setValue(oldTruePred, TRUE);
      db->setValue(newTruePred, FALSE);
      
      delete newTruePred;
    }
    else
    {
        // Pred not in block: Count gndings with pred flipped
      flipped = true;

        //set gndPred to have the opposite of its actual value
      db->setValue(gndPred, opp);

      numTrueGndOpp +=
        countNumTrueGroundingsForAllComb(gndPredIndexes, gndPred, opp, flipped,
                                        domain, hasUnknownPreds, sampleClauses,
                                        tiedClauses);

      db->setValue(gndPred, actual);
    }
    //cout << "numTrueGndOpp    = " << numTrueGndOpp << endl;

    deleteVarIdToVarsGroundedType();
    return numTrueGndOpp - numTrueGndActual;
  }


 private:
  static void addVarId(Term* const & t, const int& typeId, 
                const Domain* const & domain,                
                Array<VarsGroundedType*>*& vgtArr)
  {
    int id = -(t->getId());
    assert(id > 0);
    if (id >= vgtArr->size()) vgtArr->growToSize(id+1,NULL);
    VarsGroundedType*& vgt = (*vgtArr)[id];
    if (vgt == NULL) 
    {
      vgt = new VarsGroundedType; 
      // vgt->isGrounded init to false
      vgt->typeId = typeId;
      assert(vgt->typeId >= 0);
      vgt->numGndings = domain->getNumConstantsByType(vgt->typeId);
      assert(vgt->numGndings > 0);
    }
    assert(typeId == vgt->typeId);
    vgt->vars.append(t);
  }
  
  void createVarIdToVarsGroundedType(const Domain* const & domain,
                                     Array<VarsGroundedType*>*& vgtArr) const
  {    
      //for each predicate 
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* p = (*predicates_)[i];
        //for each variable of the predicate
      for (int j = 0; j < p->getNumTerms(); j++)
      {
        Term* t = (Term*) p->getTerm(j);
        if (t->getType() == Term::VARIABLE)
          addVarId(t, p->getTermTypeAsInt(j), domain, vgtArr);
        else
        if (t->getType() == Term::FUNCTION)
        {
          cout << "Clause::createVarIdToVarsGroundedType() not expecting a "
               << "FUNCTION term" << endl;
          exit(-1);
        }
      }// for each variable of the predicate
    } // for each predicate
  }


  void createVarIdToVarsGroundedType(const Domain* const & domain)
  {    
    deleteVarIdToVarsGroundedType();
    varIdToVarsGroundedType_ = new Array<VarsGroundedType*>;
    createVarIdToVarsGroundedType(domain, varIdToVarsGroundedType_);
  }


  static void deleteVarIdToVarsGroundedType(Array<VarsGroundedType*>*& vgtArr)
  {
    for (int i = 0; i < vgtArr->size(); i++)
      if ((*vgtArr)[i]) delete (*vgtArr)[i];
    delete vgtArr;
    vgtArr = NULL;
  }


  void deleteVarIdToVarsGroundedType()
  {
    if (varIdToVarsGroundedType_)
      deleteVarIdToVarsGroundedType(varIdToVarsGroundedType_);
  }


  void getVarOrder(IntHashArray& varAppearOrder) const
  {
    varAppearOrder.clear();
    for (int i = 0; i < predicates_->size(); i++)
    {
      const Predicate* p = (*predicates_)[i];
      for (int j = 0; j < p->getNumTerms(); j++)
      {
        const Term* t = p->getTerm(j);
        if (t->getType() == Term::VARIABLE)
        {
          int id = -(t->getId());
          assert(id > 0);
          varAppearOrder.append(id);
        }
      }
    }
  }


  void createVarIdToVarsGroundedType(Array<VarsGroundedType*>*& vgtArr)
  {    
      //for each predicate 
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* p = (*predicates_)[i];
        //for each variable of the predicate
      for (int j = 0; j < p->getNumTerms(); j++)
      {
        Term* t = (Term*) p->getTerm(j);
        if (t->getType() == Term::VARIABLE)
        {
          assert(t->getId()<0);
          int id = -(t->getId());
          if (id >= vgtArr->size()) vgtArr->growToSize(id+1,NULL);
          VarsGroundedType*& vgt = (*vgtArr)[id];
          if (vgt == NULL)  vgt = new VarsGroundedType; 
          vgt->vars.append(t);
        }
        assert(t->getType() != Term::FUNCTION);
      }// for each variable of the predicate
    } // for each predicate    
  }
  

  void createVarConstIdToTerms(Array<Array<Term*>*>*& varIdToTerms,
                               hash_map<int,Array<Term*>*>*& constIdToTerms)
  {    
    for (int i = 0; i < predicates_->size(); i++) //for each predicate 
    {
      Predicate* p = (*predicates_)[i];
      for (int j = 0; j < p->getNumTerms(); j++) //for each term of predicate
      {
        Term* t = (Term*) p->getTerm(j);
        if (t->getType() == Term::VARIABLE)
        {
          if (varIdToTerms == NULL) continue;
          int id = -(t->getId());
          if (id >= varIdToTerms->size()) varIdToTerms->growToSize(id+1,NULL);
          Array<Term*>*& termArr = (*varIdToTerms)[id];
          if (termArr == NULL) termArr = new Array<Term*>;             
          termArr->append(t);
        }
        else
        if (t->getType() == Term::CONSTANT)
        {
          if (constIdToTerms == NULL) continue;
          int id = t->getId();
          Array<Term*>* termArr;
          hash_map<int,Array<Term*>*>::iterator it = constIdToTerms->find(id);
          if (it == constIdToTerms->end()) 
          {
            termArr = new Array<Term*>;
            (*constIdToTerms)[id] = termArr;
          }
          else
            termArr = (*it).second;
          termArr->append(t);
        }
        else
        if (t->getType() == Term::FUNCTION)
        {
          cout << "Clause::createVarIdToTerms() not expecting a "
               << "FUNCTION term" << endl;
          exit(-1);
        }
      }// for each variable of the predicate
    } // for each predicate
  }


  void deleteVarConstIdToTerms(Array<Array<Term*>*>*& varIdToTerms,
                               hash_map<int,Array<Term*>*>*&  constIdToTerms)
  {
    if (varIdToTerms)
    {
      for (int i = 0; i < varIdToTerms->size(); i++)
        if ((*varIdToTerms)[i]) delete (*varIdToTerms)[i];
      delete varIdToTerms;
      varIdToTerms = NULL;
    }
    if (constIdToTerms)
    {
      hash_map<int,Array<Term*>*>::iterator it = constIdToTerms->begin();
      for (; it != constIdToTerms->end(); it++) delete (*it).second;
      delete constIdToTerms;
      constIdToTerms = NULL;
    }
  }


    //Predicate at position predIdx must be able to be grounded as gndPred.
    //Call Predicate::canBeGroundedAs() to check this.
    //After one or more invocation of this function restoreVars() should
    //be called to restore the variables in the clause to their original 
    //values. Since the variables will be restored later, the dirty_ bit is 
    //not set.
  void groundPredVars(const int& predIdx, Predicate* const& gndPred,
                      Array<VarsGroundedType*>*& vgtArr) const
  {
    assert((*predicates_)[predIdx]->canBeGroundedAs(gndPred));
    assert(predIdx < predicates_->size());
    assert(gndPred->isGrounded());

    const Predicate* pred = (*predicates_)[predIdx];
    for (int i = 0; i < pred->getNumTerms(); i++)
    {
      const Term* t = pred->getTerm(i);
      if (t->getType() == Term::VARIABLE)
      {
        int constId = gndPred->getTerm(i)->getId();
        VarsGroundedType* vgt = (*vgtArr)[-t->getId()];
        Array<Term*>& vars = vgt->vars;
        for (int j = 0; j < vars.size(); j++)  vars[j]->setId(constId);
        vgt->isGrounded = true;
      }
      assert(t->getType() != Term::FUNCTION);
    }
  }


    //Predicate at position predIdx must be able to be grounded as gndPred.
    //Call Predicate::canBeGroundedAs() to check this.
    //After one or more invocation of this function restoreVars() should
    //be called to restore the variables in the clause to their original 
    //values. Since the variables will be restored later, the dirty_ bit is 
    //not set.
  void groundPredVars(const int& predIdx, 
                      const GroundPredicate* const& gndPred) const
  {
    assert(varIdToVarsGroundedType_);
    assert((*predicates_)[predIdx]->canBeGroundedAs(gndPred));
    assert(predIdx < predicates_->size());

    const Predicate* pred = (*predicates_)[predIdx];
    for (int i = 0; i < pred->getNumTerms(); i++)
    {
      const Term* t = pred->getTerm(i);
      if (t->getType() == Term::VARIABLE)
      {
        int constId = gndPred->getTermId(i);
        VarsGroundedType* vgt = (*varIdToVarsGroundedType_)[-t->getId()];
        Array<Term*>& vars = vgt->vars;
        for (int j = 0; j < vars.size(); j++)  vars[j]->setId(constId);
        vgt->isGrounded = true;
      }
      assert(t->getType() != Term::FUNCTION);
    }
  }


    //Predicate at position predIdx must be able to be grounded as gndPred.
    //Call Predicate::canBeGroundedAs() to check this.
    //After one or more invocation of this function restoreVars() should
    //be called to restore the variables in the clause to their original 
    //values. Since the variables will be restored later, the dirty_ bit is 
    //not set.
  void groundPredVars(const int& predIdx, Predicate* const& gndPred)
  {
    assert(varIdToVarsGroundedType_);
    groundPredVars(predIdx, gndPred, varIdToVarsGroundedType_); 
  }


    // restore variables to original values
  static void restoreVars(Array<VarsGroundedType*>* const & vgtArr)
  {    
    for (int i = 1; i < vgtArr->size(); i++)
    {
      if ((*vgtArr)[i] == NULL) continue;
      Array<Term*>& vars = (*vgtArr)[i]->vars;
      for (int j = 0; j < vars.size(); j++)  vars[j]->setId(-i);
      (*vgtArr)[i]->isGrounded = false;
    }
  }
  

    // restore variables to original values
  void restoreVars()
  {    
    assert(varIdToVarsGroundedType_);
    restoreVars(varIdToVarsGroundedType_);
  }


  /**
   * Adds variable ids and groundings to a LitIdxVarIdsGndings structure.
   * 
   * @param varId Id of variable being added.
   * @param varType Type of variable being added.
   * @param domain Domain in which this all takes place.
   * @param ivg LitIdxVarIdsGndings to which the ids and groundings are
   * being added.
   */
  static void addVarIdAndGndings(const int& varId, const int& varType,
                                 const Domain* const & domain,
                                 LitIdxVarIdsGndings* const & ivg)
  {
    ivg->varIds.append(varId);
    const Array<int>* constants = domain->getConstantsByType(varType);
    ivg->varGndings.appendArray(constants);
  }

  /**
   * Constructs a LitIdxVarIdsGndings structure for one literal.
   * Caller should delete the returned LitIdxVarsGndings*.
   * Returns NULL if the literal at litIdx is grounded.
   * 
   * @param lit Literal whose groundings are being constructed.
   * @param litIdx Index of literal in the clause
   * @param domain Domain in which the clause occurs
   * 
   * @return Structure containing the index, variable ids and groundings
   * for the given literal.
   */
  static LitIdxVarIdsGndings* createLitIdxVarIdsGndings(Predicate* const & lit,
                                                     const unsigned int& litIdx,
                                                   const Domain* const & domain)
  {
    LitIdxVarIdsGndings* ivg = new LitIdxVarIdsGndings;
    ivg->litIdx = litIdx;
    for (int i = 0; i < lit->getNumTerms(); i++)
    {
      const Term* t = lit->getTerm(i);
      int termType = t->getType();
      if (termType == Term::VARIABLE)
      {
        int varId = t->getId();
        if (!ivg->varIds.contains(varId))
          addVarIdAndGndings(varId, lit->getTermTypeAsInt(i), domain, ivg);
      }
      assert(t->getType() != Term::FUNCTION);
      assert(ivg->varIds.size() == ivg->varGndings.getNumArrays());
    }
    ivg->litUnseen = true;
    return ivg;
  }

  /**
   * Constructs a LitIdxVarIdsGndings structure for each literal in the
   * clause.
   * 
   * @param clauseLits Array of literals for which the groundings are being
   * constructed.
   * @param ivgArr Array where the structures are put.
   * @param domain Domain in which the clause occurs.
   * @param setGroundedClausesToNull If true, grounded clauses are set to
   * NULL; otherwise, they are not.
   */
  void createAllLitIdxVarsGndings(Array<Predicate*>& clauseLits, 
                                  Array<LitIdxVarIdsGndings*>& ivgArr,
                                  const Domain* const & domain,
                                  const bool& setGroundedClausesToNull) const
  {
    assert(varIdToVarsGroundedType_); // this must already be created
    
      // for each literal
    for (unsigned int i = 0; i < (unsigned int) clauseLits.size(); i++)
    {
        //the literal was grounded when a previous literal was grounded
      if (clauseLits[i] == NULL) continue;

      if (clausedebug >= 2)
      {
        cout << "createAllLitIdxVarsGndings lit before: ";
        clauseLits[i]->printWithStrVar(cout, domain);
        cout << endl;
      }
            
      ivgArr.append(createLitIdxVarIdsGndings(clauseLits[i], i, domain));
      
        //ground variables of the last literal we looked at throughout clause
      ArraysAccessor<int>& varGndings = ivgArr.lastItem()->varGndings;
      Array<int>& varIds = ivgArr.lastItem()->varIds;
      for (int j = 0; j < varIds.size(); j++)
      {
          // get the first constant that can be used to ground the var
        int constId = varGndings.getArray(j)->item(0);
        
          // ground all occurrences of var
        Array<Term*>& vars = (*varIdToVarsGroundedType_)[-varIds[j]]->vars;
        for (int k = 0; k < vars.size(); k++) vars[k]->setId(constId);
      }
    
        //store subsequent literals that are grounded when literal i is grounded
      Array<Predicate*>& subseqGndLits = ivgArr.lastItem()->subseqGndLits;

      for (int j = i + 1; j < clauseLits.size(); j++)
      {
        Predicate* subseqLit = clauseLits[j];
        if (subseqLit == NULL) continue;
        if (subseqLit->isGrounded()) 
        {
          subseqGndLits.append(subseqLit);
          if (setGroundedClausesToNull) clauseLits[j] = NULL;
        }
      } //for each subsequent literal

      if (clausedebug >= 2)
      {
        cout << "createAllLitIdxVarsGndings lit after: ";
        clauseLits[i]->printWithStrVar(cout, domain);
        cout << endl;
      }
    } //for each literal
  }


    // Also sets to -1 the ids of the parent terms of functions in ivgArr[i]. 
  static void deleteAllLitIdxVarsGndings(Array<LitIdxVarIdsGndings*>& ivgArr)
  { 
    for (int i = 0; i < ivgArr.size(); i++)
    {
      //ivgArr[i]->varGndings.deleteArraysAndClear();
      delete ivgArr[i];
    }
  }


  static void quicksortLiterals(pair<double,Predicate*> items[], 
                         const int& l, const int& r)
  {
    if (l >= r) return;
    pair<double,Predicate*> tmp = items[l];
    items[l] = items[(l+r)/2];
    items[(l+r)/2] = tmp;

    int last = l;
    for (int i = l+1; i <= r; i++)
      if (items[i].first > items[l].first)
      {
        ++last;
        tmp = items[last];
        items[last] = items[i];
        items[i] = tmp;
      }
    
    tmp = items[l];
    items[l] = items[last];
    items[last] = tmp;
    quicksortLiterals(items, l, last-1);
    quicksortLiterals(items, last+1, r);  
  }


    // sort literals in decreasing order of (# true groundings/# groundings)
  void sortLiteralsByTrueDivTotalGroundings(Array<Predicate*>& clauseLits,
                                            const Domain* const & domain,
                                            const Database* const & db) const
  {
    assert(predicates_->size() == clauseLits.size());

    Array<pair<double, Predicate*> > arr;
    for (int i = 0; i < clauseLits.size(); i++)
    {
      Predicate* lit = clauseLits[i];
    
        // put all the grounded literal in the front
      if (lit->isGrounded()) 
      {
        arr.append(pair<double,Predicate*>(DBL_MAX, lit));
        continue;
      }

        //estimate how likely the literal is true
      double numTrue = (lit->getSense())? db->getNumTrueGndPreds(lit->getId())
                                         :db->getNumFalseGndPreds(lit->getId());
      double numTotal = lit->getNumGroundingsIfAllVarDiff(domain);

        //get number of groundings of the literal
      double numGnd = 1;

      Array<int> varIds; //used to check unique var ids. A hash_set is slower.
      varIds.growToSize(lit->getNumTerms(),1);
      for (int i = 0; i < lit->getNumTerms(); i++)
      {
        const Term* t = lit->getTerm(i);
        if (t->getType() == Term::VARIABLE)
        {
          int tid = t->getId();
          if (!varIds.contains(tid))
          {
            varIds.append(tid);
            numGnd *= domain->getNumConstantsByType(lit->getTermTypeAsInt(i));
          }
        }
        assert(t->getType() != Term::FUNCTION);
      }

      arr.append(pair<double,Predicate*>(numTrue/numTotal/numGnd, lit));
    }
  
    quicksortLiterals((pair<double,Predicate*>*) arr.getItems(),0,arr.size()-1);
    assert(arr.size() == clauseLits.size());
    for (int i = 0; i < arr.size(); i++) clauseLits[i] = arr[i].second;
  }


  static bool literalOrSubsequentLiteralsAreTrue(Predicate* const & lit,
                                          const Array<Predicate*>& subseqLits,
                                          const Database* const & db)
  {
    TruthValue tv = db->getValue(lit);
    lit->setTruthValue(tv);
    if (db->sameTruthValueAndSense(tv, lit->getSense())) return true;
    for (int i = 0; i < subseqLits.size(); i++)
    {
      tv = db->getValue(subseqLits[i]);
      subseqLits[i]->setTruthValue(tv);
      if (db->sameTruthValueAndSense(tv,subseqLits[i]->getSense())) return true;
    }
    return false;
  }


  bool hasTwoLiteralsWithOppSense(const Database* const & db) const
  {    
    PredicateSet predSet; // used to detect duplicates
    PredicateSet::iterator iter;
    
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* predicate = (*predicates_)[i];
      assert(predicate->isGrounded());
      if (db->getValue(predicate) == UNKNOWN)
      {
        if ( (iter = predSet.find(predicate)) != predSet.end() )
        {
            // the two gnd preds are of opp sense, so clause must be satisfied
          if ((*iter)->getSense() !=  predicate->getSense()) return true;
        }
        else
          predSet.insert(predicate);        
      }
    }
      
    return false;
  }


    //returns true if the (ground) clause has two literals with opposite sense
    //i.e. the clause is satisfied; otherwise returns false
  bool createAndAddUnknownClause(Array<GroundClause*>* const& unknownGndClauses,
                                 Array<Clause*>* const& unknownClauses,
                                 double* const & numUnknownClauses,
                                 const AddGroundClauseStruct* const & agcs,
                                 const Database* const & db);

  bool groundPredicates(const Array<int>* const & set,
                        const Array<int>& gndPredIndexes,
                        Predicate* const & gndPred,
                        const TruthValue& gndPredTV,
                        const Database* const & db,
                        bool& sameTruthValueAndSense,
                        bool& gndPredPosHaveSameSense)
  {
    sameTruthValueAndSense = false;
    gndPredPosHaveSameSense = true;
    bool prevSense = false;
    for (int i = 0; i < set->size(); i++)
    {
      int gndPredIdx = gndPredIndexes[(*set)[i]];
      Predicate* pred = (*predicates_)[gndPredIdx];
        // if inconsistent grounding of variables, proceed to next combination
      if (!pred->canBeGroundedAs(gndPred)) return false;
      groundPredVars(gndPredIdx, gndPred);
      if (db->sameTruthValueAndSense(gndPredTV, pred->getSense()))
        sameTruthValueAndSense = true;
      if (i == 0) 
        prevSense = pred->getSense();
      else
      if (prevSense != pred->getSense())
        gndPredPosHaveSameSense = false;
    }
    return true;
  }


    //count the number of groundings of clause variables
  double countNumGroundings()
  {
    double n = 1;
    for (int i = 1; i < varIdToVarsGroundedType_->size(); i++)
    {
      if (!((*varIdToVarsGroundedType_)[i]->isGrounded))
        n *= (*varIdToVarsGroundedType_)[i]->numGndings;
    }
    return n;
  }
  

  static double addCountToCombination(bool inComb[], const int& inCombSize,
                               const Array<int>* const & set,
                               MultDArray<double>& gndedPredPosArr,
                               const double& count)
  {
    memset(inComb, false, inCombSize*sizeof(bool));
    for (int i = 0; i < set->size(); i++)  inComb[(*set)[i]] = true;
    Array<int> multDArrIndexes(inCombSize); //e.g. [0][1][0][1]
    for (int i = 0; i < inCombSize; i++)
    {
      if (inComb[i]) multDArrIndexes.append(1);
      else           multDArrIndexes.append(0);    
    }
    gndedPredPosArr.addItem(&multDArrIndexes, count);   
    return gndedPredPosArr.getItem(&multDArrIndexes);
  }

  
  static void minusRepeatedCounts(bool inComb[], const int& inCombSize,
                           const Array<int>& inCombIndexes,
                           const Array<int>* const & set,
                           const Array<int>* const & falseSet,
                           MultDArray<double>& gndedPredPosArr,
                           const double& count)
  {
    memset(inComb, false, inCombSize*sizeof(bool));
    for (int i = 0; i < set->size(); i++)  inComb[(*set)[i]] = true;
    
    for (int i = 0; i < falseSet->size(); i++)
    {
      int idx = inCombIndexes[(*falseSet)[i]];
      assert(inComb[idx] == true);
      inComb[idx] = false;
    }
    
    Array<int> multDArrIndexes(inCombSize); //e.g. [0][1][0][1]
    for (int i = 0; i < inCombSize; i++)
    {
      if (inComb[i]) multDArrIndexes.append(1);
      else           multDArrIndexes.append(0);    
    }      
    
      //subtract count from that of a smaller combination that includes it
    gndedPredPosArr.addItem(&multDArrIndexes, -count);
  }  


  ////////////////////// for getting unknown clauses ////////////////////////
  void getBannedPreds(Array<Predicate*>& bannedPreds,
                      const int& gndPredIdx,
                      const GroundPredicate* const & groundPred,
                      const Predicate* const & gndPred) const
  {
    assert(gndPredIdx < predicates_->size());
    for (int i = 0; i < gndPredIdx; i++)
    {
      if ( (groundPred && (*predicates_)[i]->canBeGroundedAs(groundPred)) ||
           (gndPred && (*predicates_)[i]->canBeGroundedAs((Predicate*)gndPred)))
        bannedPreds.append((*predicates_)[i]);
    }
  }



  static void createBannedPreds(Array<Predicate*>& clauseLits,
                                Array<LitIdxVarIdsGndings*>& ivgArr,
                                Array<Predicate*>& bannedPreds)
  {
    if (bannedPreds.size() == 0) return;

    int a;
    for (int i = 0; i < ivgArr.size(); i++)
    {
      LitIdxVarIdsGndings* ivg = ivgArr[i];
      a = bannedPreds.find(clauseLits[ivg->litIdx]);
      if (a >= 0) ivg->bannedPreds.append(bannedPreds[a]);
      for (int j = 0; j < ivg->subseqGndLits.size(); j++)
      {
        a = bannedPreds.find(ivg->subseqGndLits[j]);
        if (a >= 0) ivg->bannedPreds.append(bannedPreds[a]);        
      }
    }   
  }


    //the array parameter should be LitIdxVarIdsGndings.bannedGndPreds
  static bool bannedPredsAreGndedAsGndPred(
                                const Array<Predicate*>& bannedPreds,
                                    const GroundPredicate* const & groundPred,
                                    const Predicate* const & gndPred)
  {
    for (int i = 0; i < bannedPreds.size(); i++)
    {
      if ( (groundPred && bannedPreds[i]->same(groundPred)) ||
           (gndPred    && bannedPreds[i]->same((Predicate*)gndPred)) ) 
        return true;
    }
    return false;
  }


  bool containsGndPredBeforeIdx(const int& gndPredIdx, 
                                const GroundPredicate* const & groundPred,
                                const Predicate* const & gndPred)
  {
    assert(gndPredIdx < predicates_->size());
 
    if (gndPredIdx < 0) return false;

    if (groundPred)
    {
      for (int i = 0; i < gndPredIdx; i++)
        if ((*predicates_)[i]->same(groundPred)) return true;
    }
    else
    {
      assert(gndPred);
      for (int i = 0; i < gndPredIdx; i++)
        if ((*predicates_)[i]->same((Predicate*)gndPred)) return true; 
    }
    return false;
  }


  ////////////////////////////////////////////////////////////////////////////

    //Even though it is more intuitive to use a recursive function to count
    //the number of true groundings, we are not doing so in order to allow the 
    //compiler to inline it.
    //Returns the number of true groundings, unknown clauses, number of unknown
    //clauses, and satisfiability
    //If gndPredId >= 0, the returned clauses do not contain groundPred/gndPred 
    //before position gndPredIdx
    //No more than one of the array parameters can be non-NULL.
    //No more than one of the groundPred/gndPred parameters can be non-NULL  
  double countNumTrueGroundings(const Domain* const & domain,
                                const Database* const & db,
                                const bool& hasUnknownPreds,
                                const bool& checkSatOnly,
                                  // params below: find unknown clauses
                                const int& gndPredIdx,
                                const GroundPredicate* const & groundPred,
                                const Predicate* const & gndPred,
                                Array<GroundClause*>* const & unknownGndClauses,
                                Array<Clause*>* const & unknownClauses,
                                double* const & numUnknownClauses,
                                  // params below: add unknown clauses to MRF
                                const AddGroundClauseStruct* const & agcs,
                                  // params below: get active clauses and count
                            Array<GroundClause *> * const & activeGroundClauses,
                                int* const & activeClauseCnt,
                                GroundPredicateHashArray* const& seenGndPreds,
                                bool const & ignoreActivePreds,
                                bool const & getSatisfied)                                
  {
    assert(unknownGndClauses == NULL || unknownClauses == NULL);
    assert(groundPred == NULL || gndPred == NULL);
      // Assert if activeGroundClauses isn't NULL, then others are and
      // vice versa
    if (activeGroundClauses)
    {
      assert(unknownGndClauses == NULL && unknownClauses == NULL &&
             agcs == NULL);
    }
    if (unknownGndClauses || unknownClauses || agcs)
    {
      assert(activeGroundClauses == NULL);
    }
    bool getActiveClauses = (activeGroundClauses || activeClauseCnt);

    if (activeClauseCnt) *activeClauseCnt = 0;
    if (numUnknownClauses) *numUnknownClauses = 0;

    bool findUnknownClauses = (unknownGndClauses || unknownClauses || 
                               numUnknownClauses || agcs);
      // these predicates must not be grounded as groundPred/gndPred
    Array<Predicate*> bannedPreds;
    if (findUnknownClauses)
      getBannedPreds(bannedPreds, gndPredIdx, groundPred, gndPred);

    double numTrueGndings = 0;
    
      //Copy the literals so that their original order in the clause is
      //not affected by the subsequent sorting
    Array<Predicate*>* origClauseLits = new Array<Predicate*>(*predicates_);

      // Array of partially grounded clauses achieved by using the inverted
      // index
    Array<Array<Predicate*>* > partGroundedClauses;

    if ((findUnknownClauses || getActiveClauses) && useInverseIndex)
    {
      bool newIgnoreActivePreds = ignoreActivePreds;
        // If in a neg. clause, then we can only index on evidence
      if (wt_ < 0) newIgnoreActivePreds = true;
      
        // Put the indexable literals first and ground them
      sortLiteralsByNegationAndArity(*origClauseLits, newIgnoreActivePreds, db);
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
    if (clausedebug >= 2)
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

    bool stillActivating = true;
      // Go through each clause in partGroundedClauses (nodes of the branch and
      // bound algorithm if using inverted index; otherwise, the original
      // clause), ground them out and check truth values
    for (int pgcIdx = 0; pgcIdx < partGroundedClauses.size(); pgcIdx++)
    {
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
              for (int k = 0; k < vars.size(); k++) vars[k]->setId(constId);
            }
          }
        }
          // Set the preds in clauseLits to point to the original predicates_
        delete clauseLits[i];
        clauseLits[i] = (*origClauseLits)[i];
      }
      
        //simulate a stack, back/front corresponds to top/bottom
        //ivg stands for index, varIds, groundings
      Array<LitIdxVarIdsGndings*> ivgArr;
      createAllLitIdxVarsGndings(clauseLits, ivgArr, domain, true);
      if (findUnknownClauses)
        createBannedPreds(clauseLits, ivgArr, bannedPreds);
      int ivgArrIdx = 0; //store current position in ivgArr
      bool lookAtNextLit = false;
    
        // while stack is not empty
      while (ivgArrIdx >= 0 && stillActivating)
      {
          // get variable groundings at top of stack
        LitIdxVarIdsGndings* ivg = ivgArr[ivgArrIdx];
        //Predicate* lit = clauseLits[ivg->litIdx];
        Predicate* lit = (*origClauseLits)[ivg->litIdx];
        Array<int>& varIds = ivg->varIds;
        ArraysAccessor<int>& varGndings = ivg->varGndings;
        bool& litUnseen = ivg->litUnseen;
        bool hasComb;

        if (clausedebug >= 2)
        {
          cout << "Looking at lit: ";
          lit->printWithStrVar(cout, domain);
          cout << endl;
        }

          // while there are groundings of literal's variables
        while ((hasComb = varGndings.hasNextCombination()) || litUnseen)
        {
            // there may be no combinations if the literal is fully grounded
          if (litUnseen) litUnseen = false;

          if (hasComb)
          {
              //ground the literal's variables throughout the clause
            int constId;
            int v = 0; // index of varIds
              //for each variable in literal
            while (varGndings.nextItemInCombination(constId))
            {
              Array<Term*>& vars =
                (*varIdToVarsGroundedType_)[-varIds[v++]]->vars;
              for (int i = 0; i < vars.size(); i++) vars[i]->setId(constId);
            }
          }

          if (clausedebug >= 2)
          {
            cout << "Clause is now: ";
            printWithWtAndStrVar(cout, domain);
            cout << endl;
          }

            // Getting active clauses / counts
          if (getActiveClauses)
          {
              // Check if memory is available: Linux specific
              // TODO: Check available memory in Windows and Mac
#ifdef _SC_AVPHYS_PAGES
            if (sysconf(_SC_AVPHYS_PAGES) < 100)
            {
              stillActivating = false;
              cout << "Stopped activating" << endl;
            }
#endif
            if (!stillActivating) break;

              // proceed further only if:
              // 1. positive weight and partially gnded clause is unsatisfied or
              // 2. negative weight and partially gnded clause is satisfied
            bool proceed = true;
            if (wt_ >= 0 && !getSatisfied)
              proceed = isUnsatisfiedGivenActivePreds(lit, ivg->subseqGndLits,
                                                      db, ignoreActivePreds);

            if (clausedebug >= 2)
            {
              cout << " proceed " << proceed << endl;
            }

            if (proceed)
            {
                // if there are more literals
              if (ivgArrIdx + 1 < ivgArr.size())
              {
                lookAtNextLit = true;
                ivgArrIdx++; // move up stack
                break;
              }

                // Now we can check neg. clauses: if not satisfied (no true
                // literals) or satisfied with evidence atom, then do not
                // activate
              if (wt_ < 0 && !getSatisfied &&
                  !isSatisfiedGivenActivePreds(db, ignoreActivePreds))
              {
                if (clausedebug >= 2) cout << "continuing..." << endl;
                continue;
              }
          
                // At this point all the literals are grounded
                // and does not have any true literal. To make sure that
                // it is active, need to check the following two conditions:
                // 1. It may have the same literal appearing in opposite senses
                // => satisfied (and hence not active)
                // 2. It may be empty when evidence is pruned away => not active
              bool active;
              bool accumulateClauses = activeGroundClauses;
              if (!accumulateClauses)
                active = isActive(db);
              else
                active = createAndAddActiveClause(activeGroundClauses,
                                                  seenGndPreds, db,
                                                  getSatisfied);
            
              if (clausedebug >= 2) cout << "Active " << active << endl;
              if (active) (*activeClauseCnt)++;
            }
          }
            // Getting unknown clauses / counts
          else
          {
              // if literal or subsequent grounded literals are true,
            if (literalOrSubsequentLiteralsAreTrue(lit, ivg->subseqGndLits, db))
            {
              if (clausedebug >= 2)
                cout << "Clause satisfied" << endl;
              
              if (checkSatOnly) return 1;
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
            if (findUnknownClauses && 
                bannedPredsAreGndedAsGndPred(ivg->bannedPreds, groundPred,
                                             gndPred))
            {
              if (clausedebug >= 2)
                cout << "Banned preds are grounded as ground pred" << endl;
              //do nothing, will move down stack later
            }
            else
            {
                // if there are more literals
              if (ivgArrIdx + 1 < ivgArr.size())
              {
                if (clausedebug >= 2) cout << "Moving to next literal" << endl;
                lookAtNextLit = true;
                ivgArrIdx++; // move up stack
                break;
              }
                //At this point all the literals are grounded, and they are
                //either unknown or false (have truth values opposite of their
                //senses).
              bool twoLitWithOppSense = false;
              if (hasUnknownPreds)
              {
                if (hasTwoLiteralsWithOppSense(db)) 
                {
                  twoLitWithOppSense = true;
                  ++numTrueGndings;
                  if (checkSatOnly) return 1;
                }
              }

              if (!twoLitWithOppSense && findUnknownClauses)
              {
                assert(!containsGndPredBeforeIdx(gndPredIdx, groundPred,
                                                 gndPred));

                  //Create a new clause by appending unknown predicates to it.
                createAndAddUnknownClause(unknownGndClauses, unknownClauses, 
                                          numUnknownClauses, agcs, db);
              }
            }
          }
        } //while there are groundings of literal's variables

          // If we have stopped due to no memory, then exit the while loop
        if (!stillActivating) break;
        
          //if we exit the while loop in order to look at next literal 
          //(i.e. without considering all groundings of current literal)
        if (lookAtNextLit) { lookAtNextLit = false; }
          //mv down stack
        else { varGndings.reset(); litUnseen = true; ivgArrIdx--; }

      } // while stack is not empty

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

      deleteAllLitIdxVarsGndings(ivgArr);
      delete partGroundedClauses[pgcIdx];
    }
    delete origClauseLits;

    if (getActiveClauses)
    {
      if (stillActivating) return 1;
      else return 0;
    }
    return numTrueGndings;
  }


  double countNumTrueGroundingsForAllComb(const Array<int>& gndPredIndexes,
                                          Predicate* const & gndPred,
                                          const TruthValue& gndPredTV,
                                          const bool& gndPredFlipped,
                                          const Domain* const & domain,
                                          const bool& hasUnknownPreds,
                                          const bool& sampleClauses,
                                          Array<Clause*>* const & tiedClauses)
  {
    assert(varIdToVarsGroundedType_);
    const Database* db = domain->getDB();
    double count = 0;
    int inCombSize = gndPredIndexes.size();
    bool* inComb = new bool[inCombSize];

      //initialize the number of true groundings for each way grounding the
      //predicates represented by gndPredIndexes.
      //e.g. gndedPredPosArr[0][0][1][1] is the number of true groundings when
      //the 1st and 2nd predicates in gndPredIndexes are not grounded and
      //the 3rd and 4th are grounded
    Array<int> dim;
    for (int i = 0; i < gndPredIndexes.size(); i++) dim.append(2);
      //WARNING: this may take up a lot of memory when gndPred can be grounded
      //at many positions in the clause
    MultDArray<double> gndedPredPosArr(&dim);
    const Array<double>* oneDArray = gndedPredPosArr.get1DArray();
    double* oneDArrayItems = (double*) oneDArray->getItems();
    memset(oneDArrayItems, 0, oneDArray->size()*sizeof(double));

      //for each possible combination of grounding the predicates
      //(represented by gndPredIndexes) as gndPred
    PowerSet* ps = PowerSet::getPowerSet();
      //the combinations are accessed in order of decreasing size
    ps->prepareAccess(gndPredIndexes.size(), false);

    const Array<int>* set;
    while (ps->getNextSet(set))
    {
        //ground the predicates in current combination
      bool sameTruthValueAndSense; //a ground pred has the same tv and sense
      bool gndPredPosSameSense;
      bool valid = groundPredicates(set, gndPredIndexes, gndPred, gndPredTV, db,
                                    sameTruthValueAndSense,gndPredPosSameSense);

        //If it is not possible to ground the predicates according to current
        //combination or the grounded predicates are not all of the same sense,
        //then skip the combination.
        //We can ignore the combination when the grounded predicates are not all
        //of the same sense because the counts when gndPred is held to true and
        //false cancel one another out exactly
      if (!valid || !gndPredPosSameSense) { restoreVars(); continue; }
      
        //count number of true groundings
      double cnt, numGndings = countNumGroundings();
      if (tiedClauses)
      {
        for (int i = 0; i < tiedClauses->size(); i++)
        {
          double gndings = (*tiedClauses)[i]->getNumGroundings(domain);
          cnt += gndings;
          numGndings += gndings;
        }
      }

        //if any of the grounded predicates has the same truth value and sense
      if (sameTruthValueAndSense)
        cnt = numGndings;
      else
      {
        double samp; int np;
        bool toSampleClauses = (sampleClauses && (np=predicates_->size()) > 1 &&
                                (samp = clauseSampler_->computeNumSamples(np))
                                < numGndings);
        //commented out: for testing sampling only
        //toSampleClauses = (sampleClauses && np > 1);

        if (toSampleClauses)
        {
          Predicate* flippedGndPred = (gndPredFlipped) ? gndPred : NULL;
          cnt = clauseSampler_->estimateNumTrueGroundings(this, flippedGndPred,
                                                          domain, samp);
          if (cnt > numGndings) cnt = numGndings;
          else if (cnt < 0)     cnt = 0;
        }
        else
          cnt = countNumTrueGroundings(domain, db, hasUnknownPreds, false,
                                       -1, NULL, NULL, NULL, NULL, NULL, NULL,
                                       NULL, NULL, NULL, true, false);
      }

        // add cnt to that of current combination
      double cntDueToThisComb
        = addCountToCombination(inComb,inCombSize,set,gndedPredPosArr,cnt);
      count += cntDueToThisComb;

      //for (int i = 0; i < inCombSize; i++)  cout << ((int) inComb[i]) << " ";
      //cout << " = " << cntDueToThisComb << "/" << cnt << endl;
    
        //find the indexes that are in this combination
      Array<int> inCombIndexes;
      for (int i = 0; i < set->size(); i++)  inCombIndexes.append((*set)[i]);

        // subtract all the repeated counts of cntDueToThisComb
      PowerSetInstanceVars psInstVars;
      ps->prepareAccess(inCombIndexes.size(), psInstVars);
      const Array<int>* falseSet;
      while (ps->getNextSet(falseSet, psInstVars))
      {
          // at least one of the predicates must be gnded as gndPred
        if (falseSet->size() == set->size()) continue;
        minusRepeatedCounts(inComb, inCombSize, inCombIndexes, set, falseSet,
                            gndedPredPosArr, cntDueToThisComb);
      }
      restoreVars();
    } //for each possible combination of grounding the predicates as gndPred
    delete [] inComb;
    return count;
  }


    // Returns true if the ground clause was active 
  bool createAndAddActiveClause(
                            Array<GroundClause *> * const & activeGroundClauses,
                                GroundPredicateHashArray* const& seenGndPreds,
                                const Database* const & db,
                                bool const & getSatisfied);
  
  bool isActive(const Database* const & db)
  {
    PredicateSet predSet; // used to detect duplicates
    PredicateSet::iterator iter;
    bool isEmpty = true;
 
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* predicate = (*predicates_)[i];
      assert(predicate); 
	  assert(predicate->isGrounded());
      if ( (iter=predSet.find(predicate)) != predSet.end() )
      {
          // the two gnd preds are of opp sense, so clause must be satisfied
	    if (wt_ >= 0 && (*iter)->getSense() !=  predicate->getSense())
		  return false;
        continue;
	  }
      else
        predSet.insert(predicate);
	  
      bool isEvidence = db->getEvidenceStatus(predicate);
	  if (!isEvidence)
        isEmpty = false;
    }
    return !isEmpty;
  }

    // Assumption: clause has pos. weight
  bool isUnsatisfiedGivenActivePreds(Predicate* const & lit,
                                     const Array<Predicate*>& subseqLits,
                                     const Database* const & db,
                                     bool const & ignoreActivePreds)
  {
//cout << "ignoreActivePreds " << ignoreActivePreds << endl;
//cout << "lit "; lit->printWithStrVar(cout, db->getDomain()); cout << endl;
      // If atom has been deactivated, then we don't want any clauses that
      // it's in
    if (db->getDeactivatedStatus(lit)) return false;
    bool active = false;
    if (!ignoreActivePreds)
      active = db->getActiveStatus(lit);
//cout << "active " << active << endl;
    TruthValue tv = db->getValue(lit);
    lit->setTruthValue(tv);
//cout << "tv  " << tv << endl;
    if (!active && db->sameTruthValueAndSense(tv, lit->getSense()))
      return false;

    for (int i = 0; i < subseqLits.size(); i++)
    {
//cout << "subseqLit " << i << " ";
//subseqLits[i]->printWithStrVar(cout, db->getDomain());
//cout << endl;
    if (!ignoreActivePreds)
      active = db->getActiveStatus(subseqLits[i]);
//cout << "active " << active << endl;
    tv = db->getValue(subseqLits[i]);
    subseqLits[i]->setTruthValue(tv);
//cout << "tv  " << tv << endl;
    if (!active && db->sameTruthValueAndSense(tv,subseqLits[i]->getSense()))
      return false;
    }
    return true;
  }


    // Assumption: clause has neg. weight
  bool isSatisfiedGivenActivePreds(const Database* const & db,
                                   bool const & ignoreActivePreds)
  {
    bool isSatisfied = false;
    for (int i = 0; i < predicates_->size(); i++)
    {
      Predicate* lit = getPredicate(i);
      assert(lit->isGrounded());
        // If atom has been deactivated, then we don't want any clauses
        // that it's in
      if (db->getDeactivatedStatus(lit)) return false;
      bool active = false;
      bool evidence = false;
      if (!ignoreActivePreds)
        active = db->getActiveStatus(lit);
    
      TruthValue tv = db->getValue(lit);
      lit->setTruthValue(tv);
        // If true evidence atom is in clause, then clause is always satisfied
        // and we want to ignore it.
      evidence = db->getEvidenceStatus(lit);

//cout << "Lit: ";
//lit->printWithStrVar(cout, db->getDomain());
//cout << " ev. " << evidence << " act. " << active << " stvas "
//<< db->sameTruthValueAndSense(tv, lit->getSense()) << endl;
      if (evidence && db->sameTruthValueAndSense(tv, lit->getSense()))
        return false;
        // Any active atom in clause or any true non-active atom
        // means it is candidate for active clause
      if (active || db->sameTruthValueAndSense(tv, lit->getSense()))
        isSatisfied = true;
    }
    return isSatisfied;
  }


  /**
   * Sort literals for use by the inverted index.
   * For pos. clauses: negative literals first, then of increasing arity.
   * For neg. clauses: positive literals first, then of increasing arity.
   * 
   * @param clauseLits Array of Predicate* to be sorted.
   * @param ignoreActivePreds If true, active preds are not indexed, which
   * means only evidence preds are indexed (sorted to the front). This is used
   * in eager inference.
   */
  void sortLiteralsByNegationAndArity(Array<Predicate*>& clauseLits,
                                      const bool& ignoreActivePreds,
                                      const Database* const & db) const
  {
    assert(predicates_->size() == clauseLits.size());

    Array<pair<double, Predicate*> > arr;
    for (int i = 0; i < clauseLits.size(); i++)
    {
      Predicate* lit = clauseLits[i];
      bool isIndexable = lit->isIndexable(wt_ >= 0);
      
        // If lit is indexable, put it at the beginning
      if (isIndexable)
      {
          // Evidence first
        if (db->isClosedWorld(lit->getId()))
        {
            // Put neg. literals first because we assume sparseness
          if (lit->getSense())
            arr.append(pair<double, Predicate*>(5 + 
                                      (1.0 / (double)lit->getNumTerms()), lit));
          else
            arr.append(pair<double, Predicate*>(7 +
                                      (1.0 / (double)lit->getNumTerms()), lit));
        }
          // Then non-evidence if performing lazy
        else if (!ignoreActivePreds)
        {
            // Put neg. literals first because we assume sparseness
          if (lit->getSense())
            arr.append(pair<double, Predicate*>(1 + 
                                      (1.0 / (double)lit->getNumTerms()), lit));
          else
            arr.append(pair<double, Predicate*>(3 +
                                      (1.0 / (double)lit->getNumTerms()), lit));
        }
        else
          arr.append(pair<double, Predicate*>(-(double)lit->getNumTerms(), lit));        
      }
      else
        arr.append(pair<double, Predicate*>(-(double)lit->getNumTerms(), lit));
    }
  
    quicksortLiterals((pair<double,Predicate*>*) arr.getItems(),0,arr.size()-1);
    assert(arr.size() == clauseLits.size());
    for (int i = 0; i < arr.size(); i++)
    {
  	  clauseLits[i] = arr[i].second;
  	  //cout << clauseLits[i]->getSense() << " " << clauseLits[i]->getName()
      //     << endl;
    }
  }


  /**
   * Uses the inverted index to ground literals starting from the first literal
   * up to the last indexable literal.
   * Assumes clauseLits is ordered so that indexable lits are in the front.
   * 
   * @param domain Domain in which this takes place.
   * @param db Database containing the truth values / activity of the preds
   * @param clauseLits Array of Predicate*s representing the first-order
   * clause.
   * @param resultingClauses Partially grounded clauses resulting from grounding
   * indexable literals.
   */
  void groundIndexableLiterals(const Domain* const & domain,
                               const Database* const & db,
                               Array<Predicate*>& clauseLits,
                               Array<Array<Predicate*>* >& resultingClauses,
                               bool const & ignoreActivePreds)
  {
    if (clausedebug >= 1)
    {
      cout << "Grounding indexable literals for ";
      for (int i = 0; i < clauseLits.size(); i++)
      {
        clauseLits[i]->printWithStrVar(cout, domain);
        cout << " ";
      }
      cout << endl;
    }
    assert(clauseLits.size() > 0);
    bool posClause = (wt_ >= 0) ? true : false;

      // Initially, resultingClauses holds only copies of the original clause
    Array<Predicate*>* clauseLitsCopy = new Array<Predicate*>;
    clauseLitsCopy->growToSize(clauseLits.size());
    for (int i = 0; i < clauseLits.size(); i++)
      (*clauseLitsCopy)[i] = new Predicate(*clauseLits[i]);
    resultingClauses.append(clauseLitsCopy);
      // Go through each literal and, if indexable, branch and bound on it
    for (int litIdx = 0; litIdx < clauseLits.size(); litIdx++)
    {
        // Indexable lit
      if (clauseLits[litIdx]->isIndexable(posClause) &&
          (!ignoreActivePreds||db->isClosedWorld(clauseLits[litIdx]->getId())))
      {
        if (clausedebug >= 1)
          cout << "Looking at literal " << litIdx << endl;
          
          // tmpClauses holds the partially grounded clauses in one level of the
          // tree
        Array<Array<Predicate*>* > tmpClauses;
        if (clausedebug >= 1)
        {
          cout << "Resulting clauses size before: " << resultingClauses.size()
               << endl;
        }
          // pgc = partially grounded clause
        for (int pgcIdx = 0; pgcIdx < resultingClauses.size(); pgcIdx++)
        {
          if (clausedebug >= 2)
          {
            cout << "Resulting clause before " << pgcIdx << ": ";
            for (int i = 0; i < resultingClauses[pgcIdx]->size(); i++)
            {
              (*resultingClauses[pgcIdx])[i]->printWithStrVar(cout, domain);
              cout << " ";
            }
            cout << endl;
          }
	        // Branch and bound on indexable literals
	      Array<Predicate *>* indexedGndings = new Array<Predicate *>;
	        // Bound: Get true indexed gndings if negated literal,
            // otherwise false groundings
          bool trueGndings = !(clauseLits[litIdx]->getSense());
          
          // Fix: otherwise fail to retrieve neg-wt clauses if activating
          // positive query atom
          if (wt_ < 0 && !((Database *)db)->
                isClosedWorld((*resultingClauses[pgcIdx])[litIdx]->getId())) 
            trueGndings = clauseLits[litIdx]->getSense();
              
	      ((Database *)db)->getIndexedGndings(indexedGndings,
                                            (*resultingClauses[pgcIdx])[litIdx],
                                              ignoreActivePreds, trueGndings);
          if (clausedebug >= 2)
          {
            cout << "indexedGndings: " << endl;
            for (int i = 0; i < indexedGndings->size(); i ++)
            {
              cout << "\t";
              (*indexedGndings)[i]->printWithStrVar(cout, domain);
              cout << endl;
            }
          }
            // Branch on the indexed groundings: No. of gndings is the no.
            // of new branches
          int oldTmpSize = tmpClauses.size();
          tmpClauses.growToSize(oldTmpSize + indexedGndings->size());
          if (clausedebug >= 2)
          {
            cout << "tmpClause size old: " << oldTmpSize << " new: "
                 << tmpClauses.size() << endl;
          }
          
          for (int i = 0; i < indexedGndings->size(); i++)
          {
              // Reference clause is what was branched on. First, copy this
              // into tmpClauses. Later, variables are replaced with the
              // constants from the indexed grounding
            tmpClauses[oldTmpSize + i] =
              new Array<Predicate*>(resultingClauses[pgcIdx]->size());
            tmpClauses[oldTmpSize + i]->
              growToSize(resultingClauses[pgcIdx]->size());
            if (clausedebug >= 2)
              cout << "Tmp clause before " << oldTmpSize + i << ": ";
            for (int predNo = 0; predNo < tmpClauses[oldTmpSize + i]->size();
                 predNo++)
            {
              (*tmpClauses[oldTmpSize + i])[predNo] =
                new Predicate(*((*resultingClauses[pgcIdx])[predNo]));
              if (clausedebug >= 2)
              {
                (*tmpClauses[oldTmpSize + i])[predNo]->
                  printWithStrVar(cout, domain);
                cout << " ";
              }
            }
            if (clausedebug >= 2) cout << endl;

              // Ground variables throughout clause
            for (int j = 0; 
                 j < (*resultingClauses[pgcIdx])[litIdx]->getNumTerms(); j++)
            {
              const Term* term =
                (*resultingClauses[pgcIdx])[litIdx]->getTerm(j);
                // Was a variable that has now been grounded
              if (term->getType() == Term::VARIABLE)
              {
                int varId = term->getId();
                int constId = (*indexedGndings)[i]->getTerm(j)->getId();
                assert(constId >= 0);
                  // Ground this var in lit and subsequent lits
                for (int k = litIdx; k < tmpClauses[oldTmpSize + i]->size();
                     k++)
                {
                  Predicate* pred = (*tmpClauses[oldTmpSize + i])[k];
                  for (int l = 0; l < pred->getNumTerms(); l++)
                  {
                      // Matching variable
                    if (pred->getTerm(l)->getId() == varId)
                      pred->setTermToConstant(l, constId);
                  }
                }
              }
            }
              // Indexed grounding can be deleted
            delete (*indexedGndings)[i];
            if (clausedebug >= 2)
            {
              cout << "Tmp clause after " << oldTmpSize + i << ": ";
              for (int d = 0; d < tmpClauses[oldTmpSize + i]->size(); d++)
              {
                (*tmpClauses[oldTmpSize + i])[d]->printWithStrVar(cout, domain);
                cout << " ";
              }
              cout << endl;
            }
          }
          delete indexedGndings;
        }
          // Replace the clauses in previous level with the new level
          // First, get rid of the old ones
        for (int pgcIdx = 0; pgcIdx < resultingClauses.size(); pgcIdx++)
        {
          for (int i = 0; i < resultingClauses[pgcIdx]->size(); i++)
            delete (*resultingClauses[pgcIdx])[i];
          delete resultingClauses[pgcIdx];
        }
        
        if (clausedebug >= 1)
        {
          cout << "tmpClauses size: " << tmpClauses.size() << endl;
        }
        
          // Expand or shrink resultingClauses depending on the no. of new ones
        if (tmpClauses.size() > resultingClauses.size())
        {
          resultingClauses.growToSize(tmpClauses.size());
        }
        else if (resultingClauses.size() > tmpClauses.size())
        {
          resultingClauses.shrinkToSize(tmpClauses.size());
        }
          // Copy the new ones into resultingClauses
        for (int i = 0; i < tmpClauses.size(); i++)
        {
          resultingClauses[i] = tmpClauses[i];
        }

        if (clausedebug >= 1)
        {
          cout << "Resulting clauses size after: " << resultingClauses.size()
               << endl;
        }
      }
        // Reached first non-indexable literal
      else
      {
        break;
      }
    }  // For all literals

    if (clausedebug >= 2)
    {
      cout << "Resulting clauses returned from grounding indexables: " << endl;
      for (int pgcIdx = 0; pgcIdx < resultingClauses.size(); pgcIdx++)
      {
        cout << "\t";
        for (int i = 0; i < resultingClauses[pgcIdx]->size(); i++)
        {
          (*resultingClauses[pgcIdx])[i]->printWithStrVar(cout, domain);
          cout << " ";
        }
        cout << endl;
      }
    }

    return;
  }

  /**
   * get Active Clauses unifying with the given predicate - if ignoreActivePreds 
   * is true, this is equivalent to getting all the unsatisfied clauses
   * 
   * @return true if memory is still available to keep activating. If no more
   * memory is available to activate further, false is returned.
   */
  bool getActiveClausesAndCnt(Predicate*  const & gndPred,
                              const Domain* const & domain,
                              Array<GroundClause *>* const & activeGroundClauses,
                              int & activeClauseCnt,
                              GroundPredicateHashArray* const& seenGndPreds,
                              bool const & ignoreActivePreds,
                              bool const & getSatisfied)
  {  
      //create mapping of variable ids (e.g. -1) to variable addresses,
      //note whether they have been grounded, and store their types   
    createVarIdToVarsGroundedType(domain); 
    double stillActivating = 1;
    if (gndPred == NULL)
    {
      stillActivating = countNumTrueGroundings(domain, domain->getDB(), false,
                                               false, -1, NULL, NULL, NULL, NULL,
                                               NULL, NULL, activeGroundClauses,
                                               &activeClauseCnt, seenGndPreds,
                                               ignoreActivePreds, getSatisfied);
      restoreVars();
    }
    else
    {
	  assert(gndPred->isGrounded());

      	//store the indexes of the predicates that can be grounded as gndPred
      Array<int> gndPredIndexes;
      for (int i = 0; i < predicates_->size(); i++)
	    if ((*predicates_)[i]->canBeGroundedAs(gndPred))
          gndPredIndexes.append(i);
    
      const Database* db = domain->getDB();
      Array<int> unarySet;
	  unarySet.append(-1);

      activeClauseCnt = 0;
      for (int i = 0; i < gndPredIndexes.size(); i++)
      {
          //ground the predicates in current combination
        bool sameTruthValueAndSense; //a ground pred has the same tv and sense
        bool gndPredPosSameSense;
        unarySet[0] = i; //gndPredIndexes[i];
          //cout<<"size of unary predicate set "<<unarySet.size()<<endl;
          //cout<<"Element[0] = "<<unarySet[0]<<endl;
        groundPredicates(&unarySet, gndPredIndexes, gndPred,
                         db->getValue(gndPred), db, sameTruthValueAndSense,
                         gndPredPosSameSense);
        int cnt;
        stillActivating = countNumTrueGroundings(domain, domain->getDB(), false,
                                                 false, -1, NULL, NULL, NULL,
                                                 NULL, NULL, NULL,
                                                 activeGroundClauses, & cnt,
                                                 seenGndPreds, ignoreActivePreds,
                                                 getSatisfied);
        activeClauseCnt += cnt;
        restoreVars();
      }
    }

    assert(!activeGroundClauses ||
           activeGroundClauses->size() == activeClauseCnt);
    return (stillActivating == 1);
  }


public:

  /**
   * Retrieves active clauses unifying with the given predicate - if
   * ignoreActivePreds is true, this is equivalent to getting all the
   * unsatisfied clauses. Returns the groundClauses in activeGroundClauses
   * 
   * @param gndPred Predicate with which clauses must unify.
   * @param domain Domain in which the clauses occur
   * @param activeGroundClauses Array to hold the retrieved clauses.
   * @param seenGndPreds HashArray of seen ground predicates.
   * @param ignoreActivePreds If true active predicates are ignored and only
   * unsatisfied clauses are retrieved.
   * 
   * @return true if memory is still available to keep activating. If no more
   * memory is available to activate further, false is returned.
   */
  bool getActiveClauses(Predicate* const & gndPred, 
                        const Domain* const & domain,
                        Array<GroundClause *>* const & activeGroundClauses,
                        GroundPredicateHashArray * const & seenGndPreds,
                        bool const & ignoreActivePreds)
  {
    int cnt = 0;
    bool getSatisfied = false;
    return getActiveClausesAndCnt(gndPred, domain, activeGroundClauses, cnt,
                                  seenGndPreds, ignoreActivePreds, getSatisfied);
  }


  //get the count of Active Clauses unifying with the given predicate
  //- if ignoreActivePreds is true, this is equivalent to getting the
  //count of all the unsatisfied clauses
  int getActiveClauseCnt(Predicate* const & gndPred, 
                         const Domain* const & domain,
                         bool const & ignoreActivePreds)
  {
    int cnt = 0;
    Array<GroundClause *> *activeGroundClauses = NULL;
    GroundPredicateHashArray* const & seenGndPreds = NULL;
    bool getSatisfied = false;
    getActiveClausesAndCnt(gndPred, domain, activeGroundClauses, cnt,
                           seenGndPreds, ignoreActivePreds, getSatisfied);

    return cnt;
  }


private:
  
  static void sortByLen(Array<Clause*>& clauses, const int& l, const int& r)
  {
    Clause** items = (Clause**) clauses.getItems();
    if (l >= r) return;
    Clause* tmp = items[l];
    items[l] = items[(l+r)/2];
    items[(l+r)/2] = tmp;

    int last = l;
    for (int i = l+1; i <= r; i++)
      if (items[i]->getNumPredicates() < items[l]->getNumPredicates())
      {
        ++last;
        Clause* tmp = items[last];
        items[last] = items[i];
        items[i] = tmp;
      }
    
    tmp = items[l];
    items[l] = items[last];
    items[last] = tmp;
    sortByLen(clauses, l, last-1);
    sortByLen(clauses, last+1, r); 
  }


public:

  double getConstantTuples(const Domain* const & domain,
                           const Database* const & db,
                           Array<int>* const & mlnClauseTermIds,
                           const Clause* const & varClause,
                           PredicateTermToVariable * const & ptermToVar,
                           ClauseToSuperClauseMap* const & clauseToSuperClause,
                           bool useImplicit);

  void addConstantTuple(const Domain* const & domain,
                        const Database* const & db,
                        const Clause * const & origClauseLits,
                        Array<int> * const & constants,
                        Array<Variable *> * const & eqVars,
                        ClauseToSuperClauseMap* const & clauseToSuperClause,
                        bool useImplicit);

  Array<int>* updateToVarClause();

  void removeAllPredicates() { setDirty(); predicates_->clear();}

  bool isActionFactor()
  {
    if (predicates_->size() != 1)
      return false;
    if (action_)
      return true;
    return false;
  }

  void flip()
  {
    if (!isActionFactor())
    {
      cout << "ERROR: This isn't an action factor." << endl;
      cout << predicates_->size() << "\t" << isHardClause_ << endl;
      exit(1);
    }
    (*predicates_)[0]->invertSense();
  }


 private:
  double wt_;
  Array<Predicate*>* predicates_;
  Array<int>* intArrRep_;
  size_t hashCode_;
  bool dirty_;
  bool isHardClause_;
  bool locked_;

    // (*varIdToVarsGroundedType_)[v] is the VarsGroundType of variable -v
    // start accessing this array from index 1
  Array<VarsGroundedType*>* varIdToVarsGroundedType_ ;

  AuxClauseData* auxClauseData_;
  bool staticWt_;
  static ClauseSampler* clauseSampler_;
  static double fixedSizeB_;
  bool action_;
  double util_;
};


////////////////////////////////// hash /////////////////////////////////

class HashClause
{
 public:
  size_t operator()(Clause* const & c) const  { return c->hashCode(); }
};


class EqualClause
{
 public:
  bool operator()(Clause* const & c1, Clause* const & c2) const
  { return c1->same(c2); }
};


class EqualClauseOp
{
 public:
    //the auxClauseData_ of c1 and c2 must be NON-NULL
  bool operator()(Clause* const & c1, Clause* const & c2) const
  { 
    AuxClauseData* acd1 = c1->getAuxClauseData();
    AuxClauseData* acd2 = c2->getAuxClauseData();
    if (acd1->op != acd2->op) return false;
    if (acd1->op == OP_ADD) return c1->same(c2);
    if (acd1->op == OP_REMOVE) 
      return acd1->removedClauseIdx == acd2->removedClauseIdx;
      //acd1->op is OP_REPLACE || OP_REPLACE_ADDPRED || OP_REPLACE_REMPRED
    return acd1->removedClauseIdx == acd2->removedClauseIdx && c1->same(c2);
  }
};


/////////////////////////// containers /////////////////////////////////

typedef hash_set<Clause*, HashClause, EqualClause> ClauseSet;
typedef hash_set<Clause*, HashClause, EqualClauseOp> ClauseOpSet;
typedef hash_map<Clause*, Array<Clause*>*, HashClause, EqualClause> 
  ClauseToClausesMap;

typedef HashList<Clause*, HashClause, EqualClause> ClauseHashList;
//typedef HashArray<Clause*, HashClause, EqualClause> ClauseHashArray;
typedef HashArray<Clause*, HashClause, EqualClauseOp> ClauseOpHashArray;

#endif

