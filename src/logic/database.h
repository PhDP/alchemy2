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
#ifndef DATABASE_H_JUN_28_2005
#define DATABASE_H_JUN_28_2005

#include <set>
#include "hash.h"
#include "groundpreds.h"
#include "domain.h"
#include "arraysaccessor.h"
#include "groundpredicate.h"

  // Used by Database::paramMultByPred_. 
  // First value is the multiplier and the second is the class id.
typedef pair<unsigned long long, unsigned int> MultAndType;

const int dbdebug = 0;
//const int dbdebug = 1;

// T_INDEX = True ev., F_INDEX = False ev., A_INDEX = Active ev.
enum IndexType { T_INDEX = 0, F_INDEX = 1, A_INDEX = 2, INDEX_TYPE_COUNT = 3 };

struct NumTrueFalse 
{ 
  NumTrueFalse(): numTrue(0), numFalse(0) {}  
  int numTrue; 
  int numFalse; 
};


////////////////////////////////// hash /////////////////////////////////
class HashLongLong
{
 public:
  size_t operator()(const unsigned long long ld) const
  {
    return (size_t)ld;
  }
};


class EqualLongLong
{
 public:
  bool operator()(const unsigned long long ld1, const unsigned long long ld2)
  const
  { return ld1 == ld2; }
};

////////////////////////////////// containers /////////////////////////////////

typedef hash_set<unsigned long long, HashLongLong, EqualLongLong>
    LongLongHashSet;

typedef hash_map<pair<unsigned int, unsigned int>, set<unsigned long long>,
                 IntPairHash>
    LongLongHashMap;

typedef hash_map<unsigned long long, double, HashLongLong, EqualLongLong>
    LongLongToDoubleMap;

class Database
{
 public:
    //ASSUMPTION: constants of the same type are ordered consecutively
    //            should have called d's reorderConstants()
  Database(Domain* const & d, const Array<bool>& closedWorld, 
           const bool &storeGndPreds) : domain_(d)
  {
	lazyFlag_ = false;
	performingInference_ = false;
	assert(closedWorld.size() == domain_->getNumPredicates());
    int numFOPreds = domain_->getNumPredicates();

    closedWorld_.growToSize(numFOPreds);
    memcpy((void*)closedWorld_.getItems(), closedWorld.getItems(),
           closedWorld.size()*sizeof(bool));
	
    termMultByPred_.growToSize(numFOPreds, NULL);

    truePredIdxSet_ = new Array<LongLongHashSet>(numFOPreds);
    falsePredIdxSet_ = new Array<LongLongHashSet>(numFOPreds);
    numberOfGroundings_.growToSize(numFOPreds);
    truePredIdxSet_->growToSize(numFOPreds);
    falsePredIdxSet_->growToSize(numFOPreds);

      //Initialize inverted index
    trueEvIndex_ = new Array<LongLongHashMap>(numFOPreds);
    trueEvIndex_->growToSize(numFOPreds);
    falseEvIndex_ = new Array<LongLongHashMap>(numFOPreds);
    falseEvIndex_->growToSize(numFOPreds);
    activeIndex_ = new Array<LongLongHashMap>(numFOPreds);
    activeIndex_->growToSize(numFOPreds);

    for (int i = 0; i < numFOPreds; i++)
    {
        // if this is a '=' pred, leave termMultByPred_[i]
        // as NULL
      const PredicateTemplate* t = domain_->getPredicateTemplate(i);
	  if (t->isEqualPredWithType()) continue;
		// Internal predicates are included!
		
	  const Array<int>* termTypes = domain_->getPredicateTermTypesAsInt(i);
      Array<MultAndType>* matArr = new Array<MultAndType>;
      int numTermTypes = termTypes->size(); 
      matArr->growToSize(numTermTypes);
      unsigned long long curMult = 1;
      for (int j = numTermTypes-1; j >= 0; j--)
      {
        int typeId = (*termTypes)[j];
        (*matArr)[j] = MultAndType(curMult,typeId);
        curMult *= domain_->getNumConstantsByTypeWithExt(typeId);
      }
      termMultByPred_[i] = matArr;

      int numGnd = 1;
      for (int j = 0; j < t->getNumTerms(); j++)
      {
        int numConstByType =
          domain_->getNumConstantsByType(t->getTermTypeAsInt(j));
        numGnd *= numConstByType;
      }
      numberOfGroundings_[i] = numGnd;
    }
  
    predIdToNumTF_ = new Array<NumTrueFalse>(numFOPreds);
    predIdToNumTF_->growToSize(numFOPreds);

    activePredIdxSet_ = NULL;
    evidencePredIdxSet_ = NULL;
    deactivatedPredIdxSet_ = NULL;
  }

    //Copy constructor
  Database(const Database& db)
  {
    domain_ = db.domain_;
  	closedWorld_.growToSize(db.closedWorld_.size());
    memcpy((void*)closedWorld_.getItems(), db.closedWorld_.getItems(),
           db.closedWorld_.size()*sizeof(bool));
    //closedWorld_ = db.closedWorld_;
    lazyFlag_ = db.lazyFlag_;
    performingInference_ = db.performingInference_;
    numberOfGroundings_ = db.numberOfGroundings_;
    
    oppEqGndPreds_ = db.oppEqGndPreds_;

	termMultByPred_.growToSize(db.termMultByPred_.size(), NULL);
	for (int i = 0; i < db.termMultByPred_.size(); i++)
	{
	  if (db.termMultByPred_[i])
	  {
	    Array<MultAndType>* matArr = new Array<MultAndType>;
      	matArr->growToSize(db.termMultByPred_[i]->size());
      	for (int j = 0; j < db.termMultByPred_[i]->size(); j++)
      	{
		  MultAndType mat;
		  mat.first = (*db.termMultByPred_[i])[j].first;
		  mat.second = (*db.termMultByPred_[i])[j].second;
          (*matArr)[j] = mat;
      	}
      	termMultByPred_[i] = matArr;		
	  }
	}

    if (db.trueEvIndex_)
    {
      trueEvIndex_ = new Array<LongLongHashMap>(db.trueEvIndex_->size());
      trueEvIndex_->growToSize(db.trueEvIndex_->size());
      for (int i = 0; i < db.trueEvIndex_->size(); i++)
      {
        (*trueEvIndex_)[i] = (*db.trueEvIndex_)[i];
      }
    }

    if (db.falseEvIndex_)
    {
      falseEvIndex_ = new Array<LongLongHashMap>(db.falseEvIndex_->size());
      falseEvIndex_->growToSize(db.falseEvIndex_->size());
      for (int i = 0; i < db.falseEvIndex_->size(); i++)
      {
        (*falseEvIndex_)[i] = (*db.falseEvIndex_)[i];
      }
    }

    if (db.activeIndex_)
    {
      activeIndex_ = new Array<LongLongHashMap>(db.activeIndex_->size());
      activeIndex_->growToSize(db.activeIndex_->size());
      for (int i = 0; i < db.activeIndex_->size(); i++)
      {
        (*activeIndex_)[i] = (*db.activeIndex_)[i];
      }
    }

	if (db.predIdToNumTF_)
	{
	  predIdToNumTF_ = new Array<NumTrueFalse>(db.predIdToNumTF_->size());
	  predIdToNumTF_->growToSize(db.predIdToNumTF_->size());
	  for (int i = 0; i < db.predIdToNumTF_->size(); i++)
	  {
	  	(*predIdToNumTF_)[i] = (*db.predIdToNumTF_)[i];
	  }
	}
	
	if (db.truePredIdxSet_)
	{
	  truePredIdxSet_ = new Array<LongLongHashSet>(db.truePredIdxSet_->size());
	  truePredIdxSet_->growToSize(db.truePredIdxSet_->size());
	  for (int i = 0; i < db.truePredIdxSet_->size(); i++)
	  {
	  	(*truePredIdxSet_)[i] = (*db.truePredIdxSet_)[i];
	  }
	}
	
	if (db.falsePredIdxSet_)
	{	
	  falsePredIdxSet_ = new Array<LongLongHashSet>();
	  falsePredIdxSet_->growToSize(db.falsePredIdxSet_->size());
	  for (int i = 0; i < db.falsePredIdxSet_->size(); i++)
	  	(*falsePredIdxSet_)[i] = (*db.falsePredIdxSet_)[i];
	}
	
	if (db.activePredIdxSet_)
	{
	  activePredIdxSet_ = new Array<LongLongHashSet>();
	  activePredIdxSet_->growToSize(db.activePredIdxSet_->size());
	  for (int i = 0; i < db.activePredIdxSet_->size(); i++)
	  	(*activePredIdxSet_)[i] = (*db.activePredIdxSet_)[i];
	}
	
	if (db.evidencePredIdxSet_)
	{
	  evidencePredIdxSet_ = new Array<LongLongHashSet>();
	  evidencePredIdxSet_->growToSize(db.evidencePredIdxSet_->size());
	  for (int i = 0; i < db.evidencePredIdxSet_->size(); i++)
	  	(*evidencePredIdxSet_)[i] = (*db.evidencePredIdxSet_)[i];
	}

	if (db.deactivatedPredIdxSet_)
	{
	  deactivatedPredIdxSet_ = new Array<LongLongHashSet>();
	  deactivatedPredIdxSet_->growToSize(db.deactivatedPredIdxSet_->size());
	  for (int i = 0; i < db.deactivatedPredIdxSet_->size(); i++)
	  	(*deactivatedPredIdxSet_)[i] = (*db.deactivatedPredIdxSet_)[i];
	}
  }
  
  ~Database()
  {
    for (int i = 0; i < termMultByPred_.size(); i++) 
      if (termMultByPred_[i]) delete termMultByPred_[i];
    termMultByPred_.clearAndCompress();

    for (int i = 0; i < truePredIdxSet_->size(); i++)
      (*truePredIdxSet_)[i].clear();
    delete truePredIdxSet_;
	  	
    for (int i = 0; i < falsePredIdxSet_->size(); i++)
      (*falsePredIdxSet_)[i].clear();
    delete falsePredIdxSet_;

    if (lazyFlag_)
    {
      for (int i = 0; i < activePredIdxSet_->size(); i++)
        (*activePredIdxSet_)[i].clear();
      delete activePredIdxSet_;

      for (int i = 0; i < evidencePredIdxSet_->size(); i++)
        (*evidencePredIdxSet_)[i].clear();
      delete evidencePredIdxSet_;
	  	
      for (int i = 0; i < deactivatedPredIdxSet_->size(); i++)
        (*deactivatedPredIdxSet_)[i].clear();
      delete deactivatedPredIdxSet_;
    }

    if (predIdToNumTF_)  { delete predIdToNumTF_; predIdToNumTF_ = NULL; }

    for (int i = 0; i < trueEvIndex_->size(); i++)
      (*trueEvIndex_)[i].clear();
    delete trueEvIndex_;

    for (int i = 0; i < falseEvIndex_->size(); i++)
      (*falseEvIndex_)[i].clear();
    delete falseEvIndex_;

    for (int i = 0; i < activeIndex_->size(); i++)
      (*activeIndex_)[i].clear();
    delete activeIndex_;

  }
  
  void compress()
  {
    for (int i = 0; i < termMultByPred_.size(); i++) 
      if (termMultByPred_[i]) termMultByPred_[i]->compress();
    termMultByPred_.compress();
  }

  void printInfo()
  {
    cout << "GNDINGS " << endl;
    for (int i = 0; i < numberOfGroundings_.size(); i++)
      cout << i << ": " << numberOfGroundings_[i] << endl;
      
    cout << "TRUE " << truePredIdxSet_->size() << endl;
    for (int i = 0; i < truePredIdxSet_->size(); i++)
    {
      LongLongHashSet hs = (*truePredIdxSet_)[i];
      cout << i << ": " << hs.size() << endl;
    }
    cout << "FALSE " << falsePredIdxSet_->size() << endl;
    for (int i = 0; i < falsePredIdxSet_->size(); i++)
    {
      LongLongHashSet hs = (*falsePredIdxSet_)[i];
      cout << i << ": " << hs.size() << endl;
    }
    cout << "ACTIVE " << activePredIdxSet_->size() << endl;
    cout << "EVIDENCE " << evidencePredIdxSet_->size() << endl;
    cout << "DEACTIVE " << deactivatedPredIdxSet_->size() << endl;
  }
  
  // Set the lazy flag, update the structures accordingly
  void setLazyFlag()
  {
    setLazyFlag(true);
  }

  void setLazyFlag(const bool& lf)
  {
      // If lazyFlag_ is already set this way, do nothing
    if ((lf && lazyFlag_) || (!lf && !lazyFlag_)) return;

      // If changing lazy flag, our lazy index sets become invalid
    if (activePredIdxSet_)
    {
      for (int i = 0; i < activePredIdxSet_->size(); i++)
        (*activePredIdxSet_)[i].clear();
      delete activePredIdxSet_;
      activePredIdxSet_ = NULL;
    }
    if (evidencePredIdxSet_)
    {
      for (int i = 0; i < evidencePredIdxSet_->size(); i++)
        (*evidencePredIdxSet_)[i].clear();
      delete evidencePredIdxSet_;
      evidencePredIdxSet_ = NULL;
    }
    if (deactivatedPredIdxSet_)
    {
      for (int i = 0; i < deactivatedPredIdxSet_->size(); i++)
        (*deactivatedPredIdxSet_)[i].clear();
      delete deactivatedPredIdxSet_;
      deactivatedPredIdxSet_ = NULL;
    }
    if (lf)
    {
      lazyFlag_ = true;
      int numFOPreds = domain_->getNumPredicates();

      activePredIdxSet_ = new Array<LongLongHashSet>(numFOPreds);
      evidencePredIdxSet_ = new Array<LongLongHashSet>(numFOPreds);
      deactivatedPredIdxSet_ = new Array<LongLongHashSet>(numFOPreds);
      activePredIdxSet_->growToSize(numFOPreds);
      evidencePredIdxSet_->growToSize(numFOPreds);
      deactivatedPredIdxSet_->growToSize(numFOPreds);
    }
    if (!lf)
    {
      lazyFlag_ = false;    
    }
  }
       
  // Set the performing inference flag. Returns the previous setting
  bool setPerformingInference(bool pi)
  {
    bool previous = performingInference_;
    performingInference_ = pi;
    return previous;
  }

  const Domain* getDomain() const { return domain_; }

  static bool sameTruthValueAndSense(const TruthValue& tv, const bool& sense)
  { return (tv==TRUE && sense) || (tv==FALSE && !sense); }

  /**
   * Get the truth value of a predicate in this database.
   */
  TruthValue getValue(const Predicate* const& pred) const
  {
  	if (dbdebug >= 1) cout << "Calling database::getValue" << endl;	  	
    assert(((Predicate*)pred)->isGrounded());
    Predicate* ppred = (Predicate*) pred;
    int predId = pred->getId();

      //if this is a '=' predicate
    if (pred->isEqualPredWithType())
    {
      // If the two const ids are the same, actual truth value is TRUE,
      // else FALSE
      TruthValue actual 
        = (pred->getTerm(0)->getId()==pred->getTerm(1)->getId()) ? TRUE : FALSE;

        //if pred is in oppEqGndPreds_, return opposite of actual value
      for (int i = 0; i < oppEqGndPreds_.size(); i++)
        if (ppred->same(oppEqGndPreds_[i]))
        {
          if (actual == TRUE)
          {
          	if (dbdebug >= 1) cout << "Returning FALSE" << endl;
          	return FALSE;
          }
          if (dbdebug >= 1) cout << "Returning TRUE" << endl;
          return TRUE;
        }

      if (dbdebug >= 1) cout << "Returning " << actual << endl;
      return actual;
    }
    
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return getValue(idx, predId);
  }
  
  /**
   * Get the truth value of a predicate (in GroundPredicate representation)
   * in this database.
   */
  TruthValue getValue(const GroundPredicate* const& pred) const
  {
    if (dbdebug >= 1) cout << "Calling database::getValue" << endl;     
    int predId = pred->getId();
    
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return getValue(idx, predId);
  }
 
  /**
   * Get the active status of a predicate in this database.
   */
  bool getActiveStatus(const Predicate* const& pred) const
  {
  	if (dbdebug >= 1) cout << "Calling database::getActiveStatus" << endl;
	if (pred->isEqualPredWithType()) return false;
	assert(lazyFlag_);
	assert(((Predicate*)pred)->isGrounded());
    
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return getActiveStatus(idx, pred->getId());
  }
  
  /**
   * Get the active status of a predicate (in GroundPredicate representation)
   * in this database.
   */
  bool getActiveStatus(const GroundPredicate* const& pred) const
  {
    if (dbdebug >= 1) cout << "Calling database::getActiveStatus" << endl;
    assert(lazyFlag_);

    unsigned long long idx = getIdxOfGndPredValues(pred);
    return getActiveStatus(idx, pred->getId());
  }
  

  /**
   * Get the deactivated status of a predicate in this database.
   */  
  bool getDeactivatedStatus(const Predicate* const& pred) const
  {
  	if (dbdebug >= 1) cout << "Calling database::getDeactivatedStatus" << endl;
  	if (pred->isEqualPredWithType()) return false;
	assert(lazyFlag_);
	assert(((Predicate*)pred)->isGrounded());
    
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return getDeactivatedStatus(idx, pred->getId());
  }
  
  /**
   * Get the deactivated status of a predicate (in GroundPredicate
   * representation) in this database.
   */  
  bool getDeactivatedStatus(const GroundPredicate* const& pred) const
  {
    if (dbdebug >= 1) cout << "Calling database::getDeactivatedStatus" << endl;
    assert(lazyFlag_);
    
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return getDeactivatedStatus(idx, pred->getId());
  }
    
  /**
   * Get the evidence status of a predicate in this database.
   */
  bool getEvidenceStatus(const Predicate* const& pred) const
  {	
  	if (dbdebug >= 1) cout << "Calling database::getEvidenceStatus" << endl;
	if (pred->isEqualPredWithType()) return true;
    assert(lazyFlag_);
	assert(((Predicate*)pred)->isGrounded());

    unsigned long long idx = getIdxOfGndPredValues(pred);
    return getEvidenceStatus(idx, pred->getId());
  }
  
  /**
   * Get the evidence status of a predicate (in GroundPredicate representation)
   * in this database.
   */
  bool getEvidenceStatus(const GroundPredicate* const& pred) const
  { 
    if (dbdebug >= 1) cout << "Calling database::getEvidenceStatus" << endl;
    assert(lazyFlag_);

    unsigned long long idx = getIdxOfGndPredValues(pred);
    return getEvidenceStatus(idx, pred->getId());
  }
  
  string getValueAsString(const Predicate* const& pred) const
  {
    TruthValue tv = getValue(pred);
    if (tv == TRUE)    return "TRUE";
    if (tv == FALSE)   return "FALSE";
    if (tv == UNKNOWN) return "UNKNOWN";
    assert(false); 
    return "UNKNOWN";
  }


  /**
   * Sets the value of a predicate in this database. Caller is responsible for
   * deleting pred.
   * 
   * @param flip If true, value of predicate is flipped.
   * @param ttv TruthValue to which the predicate is set (This is ignored
   * if flip is true)
   * 
   * @return Previous value of the predicate in this database.
   */
  TruthValue setValue(const Predicate* const & pred, const TruthValue& tv)
  { return setValueHelper(pred, false, tv); }

  /**
   * Sets the value of a predicate (in GroundPredicate representation)
   * in this database. Caller is responsible for deleting pred.
   * 
   * @param flip If true, value of predicate is flipped.
   * @param ttv TruthValue to which the predicate is set (This is ignored
   * if flip is true)
   * 
   * @return Previous value of the predicate in this database.
   */
  TruthValue setValue(const GroundPredicate* const & pred, const TruthValue& tv)
  { 
    return setValueHelper(pred, false, tv);
  }

  /**
   * Flips the value of a predicate in this database. Caller is responsible for
   * deleting pred.
   * 
   * @return Previous value of the predicate in this database.
   */
  TruthValue flipValue(const Predicate* const & pred)
  { return setValueHelper(pred, true, UNKNOWN); }
    
  /**
   * Flips the value of a predicate (in GroundPredicate representation)
   * in this database. Caller is responsible for deleting pred.
   * 
   * @return Previous value of the predicate in this database.
   */
  TruthValue flipValue(const GroundPredicate* const & pred)
  { return setValueHelper(pred, true, UNKNOWN); }
    
  /**
   * Sets the active status of a predicate in this database. Caller is
   * responsible for deleting pred.
   * 
   * @return Previous active status of predicate.
   */
  bool setActiveStatus(const Predicate* const & pred, const bool& as)
  { 
    if (dbdebug >= 1) cout << "Calling database::setActiveStatus" << endl;
    if (pred->isEqualPredWithType())
    {
      if (dbdebug >= 1) cout << "Returning false" << endl;
      return false;
    }
      // Active preds can not be evidence
    if (as) assert(!getEvidenceStatus(pred));
    
    assert(lazyFlag_);
    assert(((Predicate*)pred)->isGrounded());
    
    int predId = pred->getId();
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return setActiveStatus(idx, predId, as, pred, NULL);
  }
  
  /**
   * Sets the active status of a predicate (in GroundPredicate representation)
   * in this database. Caller is responsible for deleting pred.
   * 
   * @return Previous active status of predicate.
   */
  bool setActiveStatus(const GroundPredicate* const & pred, const bool& as)
  { 
    if (dbdebug >= 1) cout << "Calling database::setActiveStatus" << endl;
      // Active preds can not be evidence
    if (as) assert(!getEvidenceStatus(pred));
    assert(lazyFlag_);
    
    int predId = pred->getId();
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return setActiveStatus(idx, predId, as, NULL, pred);
  }
  
  void resetActiveStatus()
  {
  	assert(lazyFlag_);
    for (int i = 0; i < activePredIdxSet_->size(); i++)
      (*activePredIdxSet_)[i].clear();
  }

  /**
   * Sets the deactivated status of a predicate in this database. Caller is
   * responsible for deleting pred.
   * 
   * @return Previous deactivated status of predicate.
   */
  bool setDeactivatedStatus(const Predicate* const & pred, const bool& das)
  { 
  	if (dbdebug >= 1) cout << "Calling database::setDeactivatedStatus" << endl;
	if (pred->isEqualPredWithType())
	{
	  if (dbdebug >= 1) cout << "Returning false" << endl;
	  return false;
	}
	//deactivated preds can not be evidence
	if (das) assert(!getEvidenceStatus(pred));
	
	assert(lazyFlag_);
    assert(((Predicate*)pred)->isGrounded());
	
    int predId = pred->getId();
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return setDeactivatedStatus(idx, predId, das);
  }
  
  /**
   * Sets the deactivated status of a predicate (in GroundPredicate
   * representation) in this database. Caller is responsible for deleting pred.
   * 
   * @return Previous deactivated status of predicate.
   */
  bool setDeactivatedStatus(const GroundPredicate* const & pred,
                            const bool& das)
  { 
    if (dbdebug >= 1) cout << "Calling database::setDeactivatedStatus" << endl;
      // Deactivated preds can not be evidence
    if (das) assert(!getEvidenceStatus(pred));
    assert(lazyFlag_);
    
    int predId = pred->getId();
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return setDeactivatedStatus(idx, predId, das);
  }
  
  void resetDeactivatedStatus()
  {
  	assert(lazyFlag_);
    for (int i = 0; i < deactivatedPredIdxSet_->size(); i++)
      (*deactivatedPredIdxSet_)[i].clear();
  }

  /**
   * Sets the evidence status of a predicate in this database. Caller is
   * responsible for deleting pred.
   * 
   * @return Previous evidence status of predicate.
   */
  bool setEvidenceStatus(const Predicate* const & pred, const bool& es)
  {
    if (dbdebug >= 1) cout << "Calling database::setEvidenceStatus" << endl;
    if (pred->isEqualPredWithType())
    {
      if (dbdebug >= 1) cout << "Returning true" << endl;
      return true;
    }
      //active preds can not be evidence
    if (es) assert(!getActiveStatus(pred));
			  
    assert(lazyFlag_);
    assert(((Predicate*)pred)->isGrounded());

    int predId = pred->getId();
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return setEvidenceStatus(idx, predId, es);
  }
  
  /**
   * Sets the evidence status of a predicate (in GroundPredicate representation)
   * in this database. Caller is responsible for deleting pred.
   * 
   * @return Previous evidence status of predicate.
   */
  bool setEvidenceStatus(const GroundPredicate* const & pred, const bool& es)
  {
    if (dbdebug >= 1) cout << "Calling database::setEvidenceStatus" << endl;
      // Active preds can not be evidence
    if (es) assert(!getActiveStatus(pred));
    assert(lazyFlag_);

    int predId = pred->getId();
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return setEvidenceStatus(idx, predId, es);
  }
  
  void setValuesToUnknown(const Array<Predicate*>* const & gndPreds,
                          Array<TruthValue>* const & prevValues)
  {
    if (prevValues) prevValues->clear();
    for (int i = 0; i < gndPreds->size(); i++)
    {
      TruthValue prev = setValue((*gndPreds)[i], UNKNOWN);
      if (prevValues) prevValues->append(prev);
    }
  }

  void setValuesToUnknown(const GroundPredicateHashArray* const & gndPreds,
                          Array<TruthValue>* const & prevValues)
  {
    if (prevValues) prevValues->clear();
    for (int i = 0; i < gndPreds->size(); i++)
    {
      TruthValue prev = setValue((*gndPreds)[i], UNKNOWN);
      if (prevValues) prevValues->append(prev);
    }
  }


  void setValuesToGivenValues(const Array<Predicate*>* const & gndPreds,
                              const Array<TruthValue>* const & values)
  {
    assert(values);
    assert(values->size() == gndPreds->size());
    for (int i = 0; i < gndPreds->size(); i++)
      setValue((*gndPreds)[i], (*values)[i]);
  }

  void setValuesToGivenValues(const GroundPredicateHashArray* const & gndPreds,
                              const Array<TruthValue>* const & values)
  {
    assert(values);
    assert(values->size() == gndPreds->size());
    for (int i = 0; i < gndPreds->size(); i++)
      setValue((*gndPreds)[i], (*values)[i]);
  }

    //Change all ground predicates with toBeAlteredVal to newVal. 
    //gndPredValues contains the previous values of gndPreds.
  void alterTruthValue(const Array<Predicate*>* const & gndPreds, 
                       const TruthValue& tobeAlteredVal,
                       const TruthValue& newVal, 
                       Array<TruthValue>* const & gndPredValues) 
  {
    for (int i = 0; i < gndPreds->size(); i++) 
    {
	  TruthValue val = getValue((*gndPreds)[i]);
	  gndPredValues->append(val);
	  if (val == tobeAlteredVal) 
      {
        setValue((*gndPreds)[i], newVal);
	  }
    }
  }


	// Called only for evidence atoms
    // Caller should delete pred if required
  void addEvidenceGroundPredicate(Predicate* const & pred)
  {
    if (dbdebug >= 1)
      cout << "Calling database::addEvidenceGroundPredicate" << endl;
    setValue(pred, pred->getTruthValue());
    if (pred->getRealValue() > numeric_limits<double>::min())
    {
      setRealValue(pred, pred->getRealValue());
    }
      // Evidence status only needs to be set for open-world preds
    if (!closedWorld_[pred->getId()])
    {
      setEvidenceStatus(pred, true);
    }
    unsigned long long idx = getIdxOfGndPredValues(pred);
      // Add to evidence index
    if (pred->getTruthValue() == TRUE)
      addToInvertedIndex(pred, idx, T_INDEX);
    //else if (pred->getTruthValue() == FALSE)
    //  addToInvertedIndex(pred, idx, F_INDEX);
    if (dbdebug >= 1) cout << "Returning" << endl;
  }

  /**
   * Adds a grounding to the evidence inverted index. Assumption is that
   * evidence groundings have already been added to the db so their truth values
   * can bee looked up.
   * 
   * @param pred Predicate* being added to the inverted index.
   * @param sense If true, grounding occurs as pos. literal in clause;
   * otherwise, it occurs as neg. literal.
   */
  void addPredToEvidenceIndex(Predicate* const& pred, const bool& sense)
  {
    TruthValue tv = getValue(pred);
    unsigned long long idx = getIdxOfGndPredValues(pred);
      // Occurs as pos. lit. and false evidence
    if (sense && tv == FALSE)
      addToInvertedIndex(pred, idx, F_INDEX);
      // Occurs as neg. lit. and true evidence
    else if (!sense && tv == TRUE)
      addToInvertedIndex(pred, idx, T_INDEX);    
  }

  /**
   * Sets the real value of a predicate in this database.
   * 
   * @param rv Real value to which the predicate is set
   */
  void setRealValue(const Predicate* const & pred, const double& rv)
  {
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return setRealValue(idx, rv);
  }

  /**
   * Sets the real value of a predicate in this database.
   * 
   * @param rv Real value to which the predicate is set
   */
  void setRealValue(const unsigned long long& idx, const double& rv)
  {
    realValues_[idx] = rv;
  }

  double getRealValue(const Predicate* const & pred)
  {
    unsigned long long idx = getIdxOfGndPredValues(pred);
    return realValues_[idx];    
  }

  int getNumGroundings(const int& predId) const
  {
      //if it is a '=' predicate
    const PredicateTemplate* t = domain_->getPredicateTemplate(predId);    
    if (t->isEqualPredWithType())
    {
      int nc = domain_->getNumConstantsByType(t->getTermTypeAsInt(0));
      return nc*nc;
    }
    
    return numberOfGroundings_[predId];
  }

  int getNumEvidenceGndPreds(const int& predId) const
  {
      // All '=' predicates are evidence
    const PredicateTemplate* t = domain_->getPredicateTemplate(predId);    
    if (t->isEqualPredWithType())
    {
      int nc = domain_->getNumConstantsByType(t->getTermTypeAsInt(0));
      return nc*nc;
    }

      // All closed-world preds are evidence
    if (closedWorld_[predId])
    {
      return numberOfGroundings_[predId];
    }
    return (*evidencePredIdxSet_)[predId].size();
  }

  int getNumTrueGndPreds(const int& predId) const
  {
      //if it is a '=' predicate
    const PredicateTemplate* t = domain_->getPredicateTemplate(predId);    
    if (t->isEqualPredWithType())
    {
      int nc = domain_->getNumConstantsByType(t->getTermTypeAsInt(0));
      int minus = 0;
      for (int i = 0; i < oppEqGndPreds_.size(); i++)
      {
        Predicate* pred = oppEqGndPreds_[i];
        if (pred->getId() == predId &&
            pred->getTerm(0)->getId() == pred->getTerm(1)->getId())
          minus++;
      }

      return nc-minus;
    }

    return (*predIdToNumTF_)[predId].numTrue;
  }


  int getNumFalseGndPreds(const int& predId) const
  {
      //if it is a '=' predicate
    const PredicateTemplate* t = domain_->getPredicateTemplate(predId);    
    if (t->isEqualPredWithType())
    {
      int nc = domain_->getNumConstantsByType(t->getTermTypeAsInt(0));
      int minus = 0;
      for (int i = 0; i < oppEqGndPreds_.size(); i++)
      {
        Predicate* pred = oppEqGndPreds_[i];
        if (pred->getId() == predId &&
            pred->getTerm(0)->getId() != pred->getTerm(1)->getId())
          minus++;
      }

      return nc*(nc-1)-minus;
    }

    if (closedWorld_[predId])
      return numberOfGroundings_[predId]-(*predIdToNumTF_)[predId].numTrue;
    else
      return (*predIdToNumTF_)[predId].numFalse;
  }


  int getNumUnknownGndPreds(const int& predId) const
  {
      //if it is a '=' predicate
    const PredicateTemplate* t = domain_->getPredicateTemplate(predId);    
    if (t->isEqualPredWithType()) return 0;

    //assert(!closedWorld_[predId]);
    if (closedWorld_[predId]) return 0;
    
    return numberOfGroundings_[predId] -(*predIdToNumTF_)[predId].numTrue
                                       -(*predIdToNumTF_)[predId].numFalse;
  }


  int getNumPredIds() const 
  {
    return predIdToNumTF_->size();
  }

  bool isClosedWorld(const int& predId) const { return closedWorld_[predId]; }
  void setClosedWorld(const int& predId, const bool& b) 
  { closedWorld_[predId] = b; }
  const Array<bool>& getClosedWorld() const { return closedWorld_; }

  /**
   * Gets true groundings of a predicate.
   * 
   * @param predId Id of predicate whose true groundings are being retrieved.
   * @param indexedGndings True groundings are put here
   */ 
  void getTrueGndings(int predId, Array<Predicate*>* const & indexedGndings)
  {
      // True gndings from truePredIdxSet_
    LongLongHashSet::const_iterator it = (*truePredIdxSet_)[predId].begin();
    for (; it != (*truePredIdxSet_)[predId].end(); it++)
    {
      Predicate* p =
        getPredFromIdx((*it), domain_->getPredicateTemplate(predId));
      indexedGndings->append(p);
    }
  }

  /**
   * Retrieves indexed groundings of a predicate. The indexed groundings are
   * the true (and false and active) groundings of the (possibly partially
   * grounded) predicate. Caller should delete indexedGndings and elements of
   * indexedGndings, if necessary.
   * 
   * @param indexedGndings Groundings found are appended here.
   * @param pred Predicate whose groundings are found.
   * @param ignoreActivePreds If true, only true groundings are retrieved;
   * otherwise, groundings which are false and active are retrieved as well.
   * @param trueGndings If true, true indexed groundings are retrieved;
   * otherwise, false ones are retrieved.
   */
  void getIndexedGndings(Array<Predicate*>* const & indexedGndings,
  						 Predicate* const & pred,
                         bool const & ignoreActivePreds,
                         bool const & trueGndings)
  {
    int predId = pred->getId();
      // All terms are grounded: just return itself
    if (pred->isGrounded())
    {
        // Check if pred is true/false or, if getting active preds, active
      TruthValue tv = getValue(pred);
      if ((trueGndings && tv == TRUE) ||
          (!trueGndings && tv == FALSE) || 
          (!ignoreActivePreds && getActiveStatus(pred)))
      {
        Predicate* p = new Predicate(*pred);
        indexedGndings->append(p);
      }
      return;
    }
    
      // No terms are grounded: Return activePredIdxSet_ and
      // truePredIdxSet_ / falsePredIdxSet_
    if (!pred->containsConstants())
    {
      if (trueGndings)
      {
          // True gndings from truePredIdxSet_
        LongLongHashSet::const_iterator it = (*truePredIdxSet_)[predId].begin();
        for (; it != (*truePredIdxSet_)[predId].end(); it++)
        {
          Predicate* p = getPredFromIdx((*it), pred->getTemplate());
          if (pred->canBeGroundedAs(p))
            indexedGndings->append(p);
        }
      }
      else
      {
        if (closedWorld_[predId])
        {
            // False groundings is the complement of truePredIdxSet_
            // Generate all groundings and omit true ones
          Array<Predicate*> predArr;
          Predicate::createAllGroundings(predId, domain_, predArr);
          int numPreds = predArr.size();
          for (int i = 0; i < numPreds; i++)
          {
            Predicate* newPred = predArr[i];
              // Check if in true groundings
            LongLongHashSet::const_iterator it;
            if ((it = (*truePredIdxSet_)[predId].find(i)) !=
                (*truePredIdxSet_)[predId].end())
              delete newPred;
            else if (pred->canBeGroundedAs(newPred))
              indexedGndings->append(newPred);
          }
        }
        else
        {
          LongLongHashSet::const_iterator it =
            (*falsePredIdxSet_)[predId].begin();
          for (; it != (*falsePredIdxSet_)[predId].end(); it++)
          {
            Predicate* p = getPredFromIdx((*it), pred->getTemplate());
            if (pred->canBeGroundedAs(p))
              indexedGndings->append(p);
          }
        }        
      }
        // active gndings from activePredIdxSet_
      if (!ignoreActivePreds)
      {
          // Active gndings from activePredIdxSet_
        LongLongHashSet::const_iterator ait =
          (*activePredIdxSet_)[predId].begin();
        for (; ait != (*activePredIdxSet_)[predId].end(); ait++)
        {
          Predicate* p = getPredFromIdx((*ait), pred->getTemplate());
          indexedGndings->append(p);
        }
      }
      return;
    }

      // Not all terms are grounded: Compute intersection of index sets
    assert(predId < trueEvIndex_->size());
    assert(predId < falseEvIndex_->size());
      // Only true OR false gndings are retrieved
    set<unsigned long long> trueOrFalseIntersection;
    set<unsigned long long> activeIntersection;
    bool initial = true;
      // For each constant, retrieve the groundings and merge them
    for (int term = 0; term < pred->getNumTerms(); term++)
    {
      int constantId = pred->getTerm(term)->getId();
        // If at a variable, do nothing
      if (constantId < 0) continue;
      pair<unsigned int, unsigned int> placeAndConstant(term, constantId);
        // True index
      if (initial)
      {
        if (trueGndings)
          trueOrFalseIntersection = ((*trueEvIndex_)[predId])[placeAndConstant];
        else
          trueOrFalseIntersection =
            ((*falseEvIndex_)[predId])[placeAndConstant];
      }
      else
      {
        set<unsigned long long> tmpTrueOrFalseIntersection;
        if (trueGndings)
          set_intersection(trueOrFalseIntersection.begin(),
                           trueOrFalseIntersection.end(),
                           ((*trueEvIndex_)[predId])[placeAndConstant].begin(),
                           ((*trueEvIndex_)[predId])[placeAndConstant].end(),
                           inserter(tmpTrueOrFalseIntersection,
                                    tmpTrueOrFalseIntersection.begin()));
        else
          set_intersection(trueOrFalseIntersection.begin(),
                           trueOrFalseIntersection.end(),
                           ((*falseEvIndex_)[predId])[placeAndConstant].begin(),
                           ((*falseEvIndex_)[predId])[placeAndConstant].end(),
                           inserter(tmpTrueOrFalseIntersection,
                                    tmpTrueOrFalseIntersection.begin()));
        trueOrFalseIntersection = tmpTrueOrFalseIntersection;
	//cout << "db size " << trueOrFalseIntersection.size() << endl;
      }
        // Active index
      if (!ignoreActivePreds)
      {
        if (initial)
          activeIntersection = ((*activeIndex_)[predId])[placeAndConstant];
        else
        {
          set<unsigned long long> tmpActiveIntersection;
          set_intersection(activeIntersection.begin(),
                           activeIntersection.end(),
                           ((*activeIndex_)[predId])[placeAndConstant].begin(),
                           ((*activeIndex_)[predId])[placeAndConstant].end(),
                           inserter(tmpActiveIntersection,
                                    tmpActiveIntersection.begin()));
          activeIntersection = tmpActiveIntersection;
        }
      }
      initial = false;
    }    

      // Now, intersection sets contain the intersections of groundings
      // Make preds and append them to return array
      // True intersection
    set<unsigned long long>::const_iterator it =
      trueOrFalseIntersection.begin();
    for (; it != trueOrFalseIntersection.end(); it++)
    {
      Predicate* p = getPredFromIdx((*it), pred->getTemplate());
      indexedGndings->append(p);
    }
      // active intersection
    if (!ignoreActivePreds)
    {
      it = activeIntersection.begin();
      for (; it != activeIntersection.end(); it++)
      {
        Predicate* p = getPredFromIdx((*it), pred->getTemplate());
        indexedGndings->append(p);
      }
    }
    return;
  }


  /**
   * Reindexes all predicates with new constant ids. This is necessary after
   * constant ids have been changed.
   * 
   * @param oldToNewConstIds Mapping of old to new constant ids.
   */
  void changeConstantsToNewIds(hash_map<int,int>& oldToNewConstIds)
  {
      // Change in all Arrays of hash_sets
      // - compute new first-const/multi, but do not overwrite old struc yet
      // - update predIdSet 
      //   . calls getIdOfGndValues/getPredFromIdx w. old first-const/multi
      //   . compute new hashcode w. new first-const/multi
      //   . update hashcode
      // - update first-const/multi
      
    // compute new termMultByPred_
    int numFOPreds = domain_->getNumPredicates();
    Array<Array<MultAndType>*> newtermMultByPred;
    newtermMultByPred.growToSize(numFOPreds, NULL);
    for (int i = 0; i < numFOPreds; i++)
    {
        // if this is a '=' pred, leave termMultByPred_[i] as NULL
      const PredicateTemplate* t = domain_->getPredicateTemplate(i);
      if (t->isEqualPredWithType()) continue;
        
      const Array<int>* termTypes = domain_->getPredicateTermTypesAsInt(i);
      Array<MultAndType>* matArr = new Array<MultAndType>;
      int numTermTypes = termTypes->size(); 
      matArr->growToSize(numTermTypes);
      unsigned long long curMult = 1;
      for (int j = numTermTypes - 1; j >= 0; j--)
      {
        int typeId = (*termTypes)[j];
        (*matArr)[j] = MultAndType(curMult, typeId);
        curMult *= domain_->getNumConstantsByTypeWithExt(typeId);
        //curMult *= domain_->getNumConstantsByType(typeId);
      }
      //delete termMultByPred_[i];
      newtermMultByPred[i] = matArr;
    }
/*
      // update predIdSet: retrieve w. old DS, compute hashcode w. new DS
    changeConstantsToNewIds(truePredIdxSet_, oldToNewConstIds,
                            newtermMultByPred);
    changeConstantsToNewIds(falsePredIdxSet_, oldToNewConstIds,
                            newtermMultByPred);
    changeConstantsToNewIds(activePredIdxSet_, oldToNewConstIds,
                            newtermMultByPred);
    changeConstantsToNewIds(evidencePredIdxSet_, oldToNewConstIds,
                            newtermMultByPred);
    changeConstantsToNewIds(deactivatedPredIdxSet_, oldToNewConstIds,
                            newtermMultByPred);
    changeConstantsToNewIds(trueEvIndex_, oldToNewConstIds,
                            newtermMultByPred);
    changeConstantsToNewIds(falseEvIndex_, oldToNewConstIds,
                            newtermMultByPred);
    changeConstantsToNewIds(activeIndex_, oldToNewConstIds,
                            newtermMultByPred);
*/
      // update first-const/multi
    for (int i = 0; i < numFOPreds; i++)
    {
      delete termMultByPred_[i];
    }
    termMultByPred_ = newtermMultByPred;
  }



 private:

  /**
   * Retrieves a set of true or false ground predicates.
   * 
   * @param predId Id of predicate whose groundings should be retrieved.
   * @param truthValue If true, true groundings are retrieved, otherwise
   * false groundings.
   * @param withEvidence If false, only non-evidence groundings are retrieved,
   * otherwise evidence and non-evidence are retrieved.
   */ 
  const PredicateHashArray* getGndPreds(const int& predId,
                                        const bool& truthValue,
                                        const bool& withEvidence) const 
  {
    PredicateHashArray* preds = new PredicateHashArray();
    LongLongHashSet predIdxSet;
      // True or false preds?
    if (truthValue) predIdxSet = (*truePredIdxSet_)[predId];
    else predIdxSet = (*falsePredIdxSet_)[predId];

    LongLongHashSet::const_iterator it = predIdxSet.begin();
    for (; it != predIdxSet.end(); it++)
    {
        // Check if in evidence if only non-evidence are to be included
      if (withEvidence ||
          (*evidencePredIdxSet_)[predId].find(*it) ==
          (*evidencePredIdxSet_)[predId].end())
      {
        Predicate* p = getPredFromIdx((*it),
                                      domain_->getPredicateTemplate(predId));
        preds->append(p);
      }
    }
    return preds;
  }

  /**
   * Get the index of the ground predicate in Predicate representation. This
   * does exactly the same thing as getIdxOfGndPredValues(GroundPredicate*).
   * The index is computed from the number of terms, the id of the predicate
   * and the ids of the terms which are identical in both representations.
   * This should be fixed so that the two functions share code, i.e. Predicate
   * and GroundPredicate should be a subclass of a common class. This, however,
   * slows down the code because it prevents inlining of the functions.
   * 
   * @param pred Predicate for which the index is computed.
   * 
   * @return Index of the Predicate.
   */
  unsigned long long getIdxOfGndPredValues(const Predicate* const & pred,
                    const Array<Array<MultAndType>*> & currtermMultByPred) const
  {
    int numTerms = pred->getNumTerms();
    Array<MultAndType>* multAndTypes = currtermMultByPred[pred->getId()];
    unsigned long long idx = 0;
    for (int i = 0; i < numTerms; i++)
    {
      int constId = pred->getTerm(i)->getId();
      assert(constId >= 0);
      int constIdx =
        domain_->getConstantIndexInType(constId, (*multAndTypes)[i].second);
      assert(constIdx >= 0);
        // idx += mutliplier * num of constants belonging to type d]);
      idx += (*multAndTypes)[i].first * constIdx;
    }
    return idx;
  }

    // default: use member vars
  unsigned long long getIdxOfGndPredValues(const Predicate* const & pred) const
  {
    return getIdxOfGndPredValues(pred, termMultByPred_);
  }

  /**
   * Get the index of the ground predicate in GroundPredicate representation.
   * This does exactly the same thing as getIdxOfGndPredValues(Predicate*).
   * The index is computed from the number of terms, the id of the predicate
   * and the ids of the terms which are identical in both representations.
   * This should be fixed so that the two functions share code, i.e. Predicate
   * and GroundPredicate should be a subclass of a common class. This, however,
   * slows down the code because it prevents inlining of the functions.
   * 
   * @param pred GroundPredicate for which the index is computed.
   * 
   * @return Index of the GroundPredicate.
   */  
  unsigned long long getIdxOfGndPredValues(const GroundPredicate* const & pred,
                    const Array<Array<MultAndType>*> & currtermMultByPred) const
  {
    int numTerms = pred->getNumTerms();
    Array<MultAndType>* multAndTypes = currtermMultByPred[pred->getId()];
    unsigned long long idx = 0;
    for (int i = 0; i < numTerms; i++)
    {
      int constId = pred->getTermId(i);
      assert(constId >= 0);
        // idx += mutliplier * num of constants belonging to type d]);
      idx += (*multAndTypes)[i].first 
        * domain_->getConstantIndexInType(constId, (*multAndTypes)[i].second);
    }
    return idx;
  }

    // default: use member vars
  unsigned long long getIdxOfGndPredValues(const GroundPredicate* const & pred)
  const
  {
    return getIdxOfGndPredValues(pred, termMultByPred_);
  }

  /**
   * Create a Predicate based on the index in this database and its template.
   * 
   * @param idx Index of the predicate in this database.
   * @param predTemplate PredicateTemplate of the predicate being created.
   * 
   * @return Predicate created.
   */
  Predicate* getPredFromIdx(unsigned long long idx,
                            const PredicateTemplate* const & predTemplate,
                    const Array<Array<MultAndType>*> & currtermMultByPred) const
  {
    Predicate* p = new Predicate(predTemplate);
    Array<unsigned long long>* auxIdx =
      new Array<unsigned long long>(predTemplate->getNumTerms());
    auxIdx->growToSize(predTemplate->getNumTerms());
    Array<MultAndType>* multAndTypes = currtermMultByPred[p->getId()];
//    for (int i = 0; i < predTemplate->getNumTerms(); i++)
//    {
//      unsigned long long aux = 0;
//      for (int j = 0; j < i; j++)
//      {
//        aux += (*auxIdx)[j] * ((*multAndTypes)[j].first);
//      }
//      (*auxIdx)[i] = (idx - aux) / (*multAndTypes)[i].first;
//      int constId =
//        (int)(*auxIdx)[i] + currfirstConstIdByType[(*multAndTypes)[i].second];
//      p->setTermToConstant(i, constId);
//    }
//    delete auxIdx;
    for (int i = 0; i < predTemplate->getNumTerms(); i++)
    {
      unsigned long long aux = 0;
      for (int j = 0; j < i; j++)
      {
        aux += (*auxIdx)[j] * ((*multAndTypes)[j].first);
      }
      (*auxIdx)[i] = (idx - aux) / (*multAndTypes)[i].first;
      const Array<int>* cbt =
        domain_->getConstantsByType((*multAndTypes)[i].second);
      int constId = (*cbt)[(int)(*auxIdx)[i]];
      p->setTermToConstant(i, constId);
    }
    delete auxIdx;
    return p;
  }

  Predicate* getPredFromIdx(unsigned long long idx,
                            const PredicateTemplate* const & predTemplate) const
  {
    return getPredFromIdx(idx, predTemplate, termMultByPred_);
  }

  /**
   * Reindexes all predicates in an Array of hash_sets with new constant ids.
   * This is necessary after constant ids have been changed.
   * 
   * @param predIdxSets Array of hash_sets containing the pred indices
   * @param oldToNewConstIds Mapping of old to new constant ids.
   */
  void changeConstantsToNewIds(Array<LongLongHashSet >* predIdxSets,
                               hash_map<int,int>& oldToNewConstIds,
                          const Array<Array<MultAndType>*> & currtermMultByPred)
  {
    for (int i = 0; i < predIdxSets->size(); i++)
    {
      LongLongHashSet newHashSet;
      LongLongHashSet& predIdxSet = (*predIdxSets)[i];
      LongLongHashSet::const_iterator it = predIdxSet.begin();
      for (; it != predIdxSet.end(); it++)
      {
        Predicate* p = getPredFromIdx((*it), domain_->getPredicateTemplate(i));
        
        for (int j = 0; j < p->getNumTerms(); j++)
        {
          Term* t = (Term*) p->getTerm(j);
          if (t->getType() == Term::CONSTANT)
          {
            int oldId = t->getId();
            assert(oldToNewConstIds.find(oldId) != oldToNewConstIds.end());
            t->setId(oldToNewConstIds[oldId]);
          }
        }

        unsigned long long d = getIdxOfGndPredValues(p, currtermMultByPred);
        newHashSet.insert(d);
        delete p;
      }
      assert(predIdxSet.size() == newHashSet.size());
      (*predIdxSets)[i] = newHashSet;
    }
  }

  /**
   * Reindexes all predicates in an Array of hash_maps with new constant ids.
   * This is necessary after constant ids have been changed.
   * 
   * @param predIdxMaps Array of hash_maps containing the pred indices
   * @param oldToNewConstIds Mapping of old to new constant ids.
   */
  void changeConstantsToNewIds(Array<LongLongHashMap>* predIdxMaps,
                               hash_map<int,int>& oldToNewConstIds,
                          const Array<Array<MultAndType>*> & currtermMultByPred)
  {
    for (int i = 0; i < predIdxMaps->size(); i++)
    {
      LongLongHashMap newHashMap;
      LongLongHashMap& predIdxMap = (*predIdxMaps)[i];
      LongLongHashMap::iterator it = predIdxMap.begin();
      for (; it != predIdxMap.end(); it++)
      {
        set<unsigned long long> newSet;
          // The second term of the pair in *it contains the set of groundings
        set<unsigned long long> gndings = (*it).second;
        set<unsigned long long>::const_iterator setit = gndings.begin();
          // Change the constant ids in each grounding
        for (; setit != gndings.end(); setit++)
        {
          Predicate* p = getPredFromIdx((*setit),
                                        domain_->getPredicateTemplate(i));
        
          for (int j = 0; j < p->getNumTerms(); j++)
          {
            Term* t = (Term*) p->getTerm(j);
            if (t->getType() == Term::CONSTANT)
            {
              int oldId = t->getId();
              assert(oldToNewConstIds.find(oldId) != oldToNewConstIds.end());
              t->setId(oldToNewConstIds[oldId]);
            }
          }

          unsigned long long d = getIdxOfGndPredValues(p, currtermMultByPred);
          delete p;
          newSet.insert(d);
        }
        //(*it).second = newSet;
          // The constant id must be changed in the firste element of *it
        pair<unsigned int, unsigned int> placeAndConstant = (*it).first;
        pair<unsigned int, unsigned int>
          newPlaceAndConstant(placeAndConstant.first,
                              oldToNewConstIds[placeAndConstant.second]);
        newHashMap[newPlaceAndConstant] = newSet;
      }
      assert(predIdxMap.size() == newHashMap.size());
      (*predIdxMaps)[i] = newHashMap;
    }
  }


  /**
   * Get the truth value of a predicate in this database.
   * 
   * @param idx Index of the predicate in this database.
   * @param predId Id of the predicate.
   */
  TruthValue getValue(const unsigned long long& idx, const int& predId) const
  {
    LongLongHashSet::const_iterator it = (*truePredIdxSet_)[predId].find(idx);
    if (it != (*truePredIdxSet_)[predId].end())
    {
      if (dbdebug >= 1) cout << "returning TRUE" << endl;
      return TRUE;
    }

    if (!closedWorld_[predId] && !performingInference_)
    {
      if((*falsePredIdxSet_)[predId].find(idx) !=
         (*falsePredIdxSet_)[predId].end())
      {
        if (dbdebug >= 1) cout << "returning FALSE" << endl;
        return FALSE;
      }
      else
      {
        if (dbdebug >= 1) cout << "returning UNKNOWN" << endl;
        return UNKNOWN;
      }
    }
    if (dbdebug >= 1) cout << "returning FALSE" << endl;
    return FALSE;
  }
  
  /**
   * Get the active status of a predicate in this database.
   * 
   * @param idx Index of the predicate in this database.
   * @param predId Id of the predicate.
   */
  bool getActiveStatus(const unsigned long long& idx, const int& predId) const
  { 
      // Evidence atoms can not be active
    if (getEvidenceStatus(idx, predId))
    {
      if (dbdebug >= 1) cout << "Returning false(1)" << endl;       
      return false;
    }

    if ((*activePredIdxSet_)[predId].find(idx) !=
        (*activePredIdxSet_)[predId].end())
    {
      if (dbdebug >= 1) cout << "Returning true" << endl;
      return true;
    }
    if (dbdebug >= 1) cout << "Returning false(2)" << endl;     
    return false;
  }

  /**
   * Get the deactivated status of a predicate in this database.
   * 
   * @param idx Index of the predicate in this database.
   * @param predId Id of the predicate.
   */  
  bool getDeactivatedStatus(const unsigned long long& idx, const int& predId) const
  {
      // Evidence atoms can not be active or deactivated
    if (getEvidenceStatus(idx, predId))
    {
      if (dbdebug >= 1) cout << "Returning false(1)" << endl;       
      return false;
    }

    if ((*deactivatedPredIdxSet_)[predId].find(idx) !=
        (*deactivatedPredIdxSet_)[predId].end())
    {
      if (dbdebug >= 1) cout << "Returning true" << endl;       
      return true;
    }
    if (dbdebug >= 1) cout << "Returning false(2)" << endl;     
    return false;
  }
  
  /**
   * Get the evidence status of a predicate in this database.
   * 
   * @param idx Index of the predicate in this database.
   * @param predId Id of the predicate.
   */
  bool getEvidenceStatus(const unsigned long long& idx, const int& predId) const
  {
      // All closed-world preds are evidence
    if (closedWorld_[predId])
    {
      if (dbdebug >= 1) cout << "Returning true(1)" << endl;        
      return true;
    }
    if ((*evidencePredIdxSet_)[predId].find(idx) !=
        (*evidencePredIdxSet_)[predId].end())
    {
      if (dbdebug >= 1) cout << "Returning true(2)" << endl;        
      return true;
    }
    if (dbdebug >= 1) cout << "Returning false" << endl;        
    return false;
  }

  /**
   * Sets the value of a predicate in this database.
   * 
   * @param flip If true, value of predicate is flipped.
   * @param ttv TruthValue to which the predicate is set (This is ignored
   * if flip is true)
   * 
   * @return Previous value of the predicate in this database.
   */
  TruthValue setValueHelper(const Predicate* const & pred, const bool& flip,
                            const TruthValue& ttv)
  {
    if (dbdebug >= 1) cout << "Calling database::setValueHelper" << endl;       
    Predicate* ppred = (Predicate *)pred;
    assert(ppred->isGrounded());
      
      // If this is a '=' predicate
    if (pred->isEqualPredWithType())
    {
      // If the two const ids are the same, actual truth value is TRUE,
      // otherwise FALSE
      TruthValue actual
        = (pred->getTerm(0)->getId() == pred->getTerm(1)->getId()) ?
          TRUE : FALSE;
      TruthValue opposite = (actual==TRUE) ? FALSE : TRUE;

        //if pred is in oppEqGndPreds_, return opposite of actual value
      for (int i = 0; i < oppEqGndPreds_.size(); i++)
        if (ppred->same(oppEqGndPreds_[i]))
        {
          if (flip || ttv == actual)
          {
            Predicate* pr = oppEqGndPreds_[i];
            oppEqGndPreds_.removeItemFastDisorder(i);
            delete pr;
          }
          if (dbdebug >= 1) cout << "returning " << opposite << endl;
          return opposite;
        }

      if (flip || ttv == opposite) oppEqGndPreds_.append(new Predicate(*pred));
      if (dbdebug >= 1) cout << "returning " << actual << endl;
      return actual;
    }

    unsigned long long idx = getIdxOfGndPredValues(pred);
    int predId = pred->getId();
    return setValueHelper(idx, predId, flip, ttv, pred, NULL);
  }

  /**
   * Sets the value of a predicate (in GroundPredicate representation)
   * in this database.
   * 
   * @param flip If true, value of predicate is flipped.
   * @param ttv TruthValue to which the predicate is set (This is ignored
   * if flip is true)
   * 
   * @return Previous value of the predicate in this database.
   */
  TruthValue setValueHelper(const GroundPredicate* const & pred,
                            const bool& flip, const TruthValue& ttv)
  {
    if (dbdebug >= 1) cout << "Calling database::setValueHelper" << endl;      
    unsigned long long idx = getIdxOfGndPredValues(pred);
    int predId = pred->getId();
    return setValueHelper(idx, predId, flip, ttv, NULL, pred);
  }

  /**
   * Sets the value of a predicate in this database.
   * 
   * Right now, pred and gndPred are being sent as parameters in order to
   * utilize the inverted index. This needs to be fixed. Exactly one of pred
   * and gndPred must be NULL.
   * 
   * @param idx Index of the predicate in this database.
   * @param predId Id of the predicate.
   * @param flip If true, value of predicate is flipped.
   * @param ttv TruthValue to which the predicate is set (This is ignored
   * if flip is true)
   * 
   * @return Previous value of the predicate in this database.
   */
  TruthValue setValueHelper(const unsigned long long& idx, const int& predId,
                            const bool& flip, const TruthValue& ttv,
                            const Predicate* const & pred,
                            const GroundPredicate* const & gndPred)
  {
    assert((pred && !gndPred) || (!pred && gndPred));
    TruthValue oldtv = getValue(idx, predId);
    //TruthValue oldtv = FALSE;
    TruthValue tv = ttv;

    if (flip) 
    {
      assert(oldtv == TRUE || oldtv == FALSE);
      if (oldtv == FALSE) tv = TRUE;
      else                tv = FALSE;
    }

    if (oldtv != tv)
    {
      if (closedWorld_[predId])
      {
        if (oldtv == TRUE)  //TRUE->FALSE
        {
          assert(tv == FALSE);
          (*predIdToNumTF_)[predId].numTrue--;
          (*truePredIdxSet_)[predId].erase(idx);
        }
        else 
        { //FALSE->TRUE
          assert(oldtv == FALSE && tv == TRUE);
          (*predIdToNumTF_)[predId].numTrue++;
          (*truePredIdxSet_)[predId].insert(idx);
        }
      }
      else
      { //open world
        if (oldtv == UNKNOWN)
        {
          if (tv == TRUE)
          { //UNKNOWN->TRUE
            (*predIdToNumTF_)[predId].numTrue++;
            (*truePredIdxSet_)[predId].insert(idx);
          }
          else
          { //UKNOWN->FALSE
            (*predIdToNumTF_)[predId].numFalse++;
            if (!performingInference_) (*falsePredIdxSet_)[predId].insert(idx);
          }
        }
        else
        if (oldtv == TRUE)
        {
          (*predIdToNumTF_)[predId].numTrue--;
          (*truePredIdxSet_)[predId].erase(idx);
          if (tv == FALSE) 
          { //TRUE->FALSE
            (*predIdToNumTF_)[predId].numFalse++;
            if (!performingInference_) (*falsePredIdxSet_)[predId].insert(idx);
          }
        }
        else
        {
          assert(oldtv == FALSE);
          (*predIdToNumTF_)[predId].numFalse--;
          if (!performingInference_) (*falsePredIdxSet_)[predId].erase(idx);
          if (tv == TRUE)
          { //FALSE->TRUE
            (*predIdToNumTF_)[predId].numTrue++;
            (*truePredIdxSet_)[predId].insert(idx);
          }
        }          
      }
    } // if (oldtv != tv)
    if (dbdebug >= 1) cout << "returning " << oldtv << endl;
    return oldtv;
  }

  /**
   * Sets the active status of a predicate in this database. Caller is
   * responsible for deleting pred.
   * 
   * Right now, pred and gndPred are being sent as parameters in order to
   * utilize the inverted index. This needs to be fixed. Exactly one of pred
   * and gndPred must be NULL.
   * 
   * @return Previous active status of predicate.
   */
  bool setActiveStatus(const unsigned long long& idx, const int& predId,
                       const bool& as, const Predicate* const & pred,
                       const GroundPredicate* const & gndPred)
  {
    assert((pred && !gndPred) || (!pred && gndPred));
    bool oldas = getActiveStatus(idx, predId);
    if (oldas != as)
    {
      if (as)
      {
        (*activePredIdxSet_)[predId].insert(idx);
      }
      else
      {
        (*activePredIdxSet_)[predId].erase(idx);
      }
      
        // Add to/remove from inverted index
      if (lazyFlag_)
      {
        if (as)
        {
          if (pred) addToInvertedIndex(pred, idx, A_INDEX);
          else if (gndPred) addToInvertedIndex(gndPred, idx, A_INDEX);
        }
        else
        {
          if (pred) removeFromInvertedIndex(pred, idx, A_INDEX);
          if (gndPred) removeFromInvertedIndex(gndPred, idx, A_INDEX);
        }
      }
    }
    if (dbdebug >= 1) cout << "Returning " << oldas << endl;
    return oldas;
  }

  /**
   * Sets the deactivated status of a predicate in this database. Caller is
   * responsible for deleting pred.
   * 
   * @return Previous deactivated status of predicate.
   */
  bool setDeactivatedStatus(const unsigned long long& idx, const int& predId,
                            const bool& das)
  {
    bool olddas = getDeactivatedStatus(idx, predId);
    if (olddas != das)
    {
      if (das)
      {
        (*deactivatedPredIdxSet_)[predId].insert(idx);
      }
      else
      {
        (*deactivatedPredIdxSet_)[predId].erase(idx);
      }
    }
    if (dbdebug >= 1) cout << "Returning " << olddas << endl;
    return olddas;
  }

  /**
   * Sets the evidence status of a predicate in this database. Caller is
   * responsible for deleting pred.
   * 
   * @return Previous evidence status of predicate.
   */
  bool setEvidenceStatus(const unsigned long long& idx, const int& predId, const bool& es)
  {  
    bool oldes = getEvidenceStatus(idx, predId);
    if (oldes != es)
    {
      if (es)
      {
        (*evidencePredIdxSet_)[predId].insert(idx);
      }
      else
      {
        (*evidencePredIdxSet_)[predId].erase(idx);
      }
    }
    if (dbdebug >= 1) cout << "Returning " << oldes << endl;
    return oldes;
  }

  /**
   * Add a predicate to the given index type.
   */
  void addToInvertedIndex(const Predicate* const & pred,
                          const unsigned long long& predIdx, IndexType idxType)
  {
      // Assumption is: Predicate is grounded
    assert(((Predicate*)pred)->isGrounded());
    int predId = pred->getId();
    assert(predId < trueEvIndex_->size());
    for (int term = 0; term < pred->getNumTerms(); term++)
    {
      int constantId = pred->getTerm(term)->getId();
      assert(constantId >= 0);
        // Insert this grounding for each constant term
      pair<unsigned int, unsigned int> placeAndConstant(term, constantId);
      pair<set<unsigned long long>::iterator, bool> prevElement;
        // True index
      if (idxType == T_INDEX)
      {
        prevElement =
          ((*trueEvIndex_)[predId])[placeAndConstant].insert(predIdx);
          // prevElement.second is false if key already in set (assert that this
          // is the same element as was tried to be inserted
        if (!prevElement.second)
          assert(predIdx == *(prevElement.first));
      }
        // False index
      else if (idxType == F_INDEX)
      {
        prevElement =
          ((*falseEvIndex_)[predId])[placeAndConstant].insert(predIdx);
          // prevElement.second is false if key already in set (assert that this
          // is the same element as was tried to be inserted
        if (!prevElement.second)
          assert(predIdx == *(prevElement.first));
      }
        // Active index
      else if (idxType == A_INDEX)
      {
        prevElement =
          ((*activeIndex_)[predId])[placeAndConstant].insert(predIdx);
          // prevElement.second is false if key already in set (assert that this
          // is the same element as was tried to be inserted
        if (!prevElement.second)
          assert(predIdx == *(prevElement.first));
      }
    }    
  }

  /**
   * Add a predicate to the given index type.
   */
  void addToInvertedIndex(const GroundPredicate* const & pred,
                          const unsigned long long& predIdx, IndexType idxType)
  {
    int predId = pred->getId();
    assert(predId < trueEvIndex_->size());
    for (unsigned int term = 0; term < pred->getNumTerms(); term++)
    {
      int constantId = pred->getTermId(term);
      assert(constantId >= 0);
        // Insert this grounding for each constant term
      pair<unsigned int, unsigned int> placeAndConstant(term, constantId);
      pair<set<unsigned long long>::iterator, bool> prevElement;
        // True index
      if (idxType == T_INDEX)
      {
        prevElement =
          ((*trueEvIndex_)[predId])[placeAndConstant].insert(predIdx);
          // prevElement.second is false if key already in set (assert that this
          // is the same element as was tried to be inserted
        if (!prevElement.second)
          assert(predIdx == *(prevElement.first));
      }
        // False index
      else if (idxType == F_INDEX)
      {
        prevElement =
          ((*falseEvIndex_)[predId])[placeAndConstant].insert(predIdx);
          // prevElement.second is false if key already in set (assert that this
          // is the same element as was tried to be inserted
        if (!prevElement.second)
          assert(predIdx == *(prevElement.first));
      }
        // Active index
      else if (idxType == A_INDEX)
      {
        prevElement =
          ((*activeIndex_)[predId])[placeAndConstant].insert(predIdx);
          // prevElement.second is false if key already in set (assert that this
          // is the same element as was tried to be inserted
        if (!prevElement.second)
          assert(predIdx == *(prevElement.first));
      }
    }    
  }

  /**
   * Remove a predicate from the given index type.
   */
  void removeFromInvertedIndex(const Predicate* const & pred,
                               const unsigned long long& predIdx,
                               IndexType idxType)
  {
      // Assumption is: Predicate is grounded
    assert(((Predicate*)pred)->isGrounded());
    int predId = pred->getId();
    assert(predId < trueEvIndex_->size());
    for (int term = 0; term < pred->getNumTerms(); term++)
    {
      int constantId = pred->getTerm(term)->getId();
      assert(constantId >= 0);
        // Remove this grounding for each possible constant term
      pair<unsigned int, unsigned int> placeAndConstant(term, constantId);
      pair<set<unsigned long long>::iterator, bool> prevElement;
      int numErased = 0;
        // True index
      if (idxType == T_INDEX)
      {
        numErased = ((*trueEvIndex_)[predId])[placeAndConstant].erase(predIdx);
        assert(numErased <= 1);
      }
        // False index
      else if (idxType == F_INDEX)
      {
        numErased = ((*falseEvIndex_)[predId])[placeAndConstant].erase(predIdx);
        assert(numErased <= 1);
      }
        // Active index
      else if (idxType == A_INDEX)
      {
        numErased =
          ((*activeIndex_)[predId])[placeAndConstant].erase(predIdx);
        assert(numErased <= 1);
      }
    }
  }

  /**
   * Remove a predicate from the given index type.
   */
  void removeFromInvertedIndex(const GroundPredicate* const & pred,
                               const unsigned long long& predIdx,
                               IndexType idxType)
  {
    int predId = pred->getId();
    assert(predId < trueEvIndex_->size());
    for (unsigned int term = 0; term < pred->getNumTerms(); term++)
    {
      int constantId = pred->getTermId(term);
      assert(constantId >= 0);
        // Remove this grounding for each possible constant term
      pair<unsigned int, unsigned int> placeAndConstant(term, constantId);
      pair<set<unsigned long long>::iterator, bool> prevElement;
      int numErased = 0;
        // True index
      if (idxType == T_INDEX)
      {
        numErased = ((*trueEvIndex_)[predId])[placeAndConstant].erase(predIdx);
        assert(numErased <= 1);
      }
        // False index
      else if (idxType == F_INDEX)
      {
        numErased = ((*falseEvIndex_)[predId])[placeAndConstant].erase(predIdx);
        assert(numErased <= 1);
      }
        // Active index
      else if (idxType == A_INDEX)
      {
        numErased =
          ((*activeIndex_)[predId])[placeAndConstant].erase(predIdx);
        assert(numErased <= 1);
      }
    }
  }

 private:
  const Domain* domain_; // not owned by Database

    // closedWorld_[p]==true: the closed world assumption is made for pred p
  Array<bool> closedWorld_; 

  	// Flag is set when performing lazy inference
  bool lazyFlag_;
  
    // Flag set when performing inference. If set, hash table of false ground
    // atoms is not used (all atoms are true or false).
  bool performingInference_;
  
    //Suppose a predicate with id i has three terms, and the number of 
    //groundings of the first, second and third terms are 
    //4, 3, 2 respectively. Then termMultByPred_[i][0].first = 3x2,
    //termMultByPred_[i][1].first = 2, and termMultByPred_[i][2].first = 1.
    //termMultByPred_[i][j].first is the product of the number of groundings
    //of terms j > i;
  Array<Array<MultAndType>*> termMultByPred_;

    // predIdToNumTF_[p] contains the number of TRUE/FALSE groundings of pred
    // with id p  
  Array<NumTrueFalse>* predIdToNumTF_;


    //'=' ground preds that are set to opposite of their actual values
    //We are using a simple array because we are setting one ground predicate 
    //at a time to the opposite of its actual value and then reverting it to 
    //its actual value,so there's no efficiency penalty when searching the array
  Array<Predicate*> oppEqGndPreds_;
  
    // Inverted index of true evidence ground atoms: For each predicate p,
    // trueEvIndex_[p] contains a hash map mapping a pair of ints (argument
    // number and constant id) to a set of true groundings.
  Array<LongLongHashMap>* trueEvIndex_;

    // Inverted index of false evidence ground atoms: For each predicate p,
    // trueEvIndex_[p] contains a hash map mapping a pair of ints (argument
    // number and constant id) to a set of false groundings.
  Array<LongLongHashMap>* falseEvIndex_;

    // Inverted index of active ground atoms: For each predicate p,
    // activeIndex_[p] contains a hash map mapping a pair of ints (argument
    // number and constant id) to a set of false, active groundings.
  Array<LongLongHashMap>* activeIndex_;

    // Hash set of true ground atoms
  Array<LongLongHashSet>* truePredIdxSet_;
  
    // Hash set of false ground atoms - only used for open world preds
  Array<LongLongHashSet>* falsePredIdxSet_;

    // Hash set of active ground atoms
  Array<LongLongHashSet>* activePredIdxSet_;

    // Hash set of evidence ground atoms
  Array<LongLongHashSet>* evidencePredIdxSet_;
  
    // Hash set of deactivated ground atoms
  Array<LongLongHashSet>* deactivatedPredIdxSet_;

    // Total number of groundings for each predicate
  Array<unsigned long long> numberOfGroundings_;
  
    // Map from pred idx to real value
  LongLongToDoubleMap realValues_;
};



#endif
