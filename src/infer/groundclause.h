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
#ifndef GROUNDCLAUSE_H_JUN_26_2005
#define GROUNDCLAUSE_H_JUN_26_2005

#include <cfloat>
#include <ext/hash_set>
using namespace __gnu_cxx;
#include "hasharray.h"
#include "array.h"
#include "hash.h"
#include <map>
#include <cstring>

using namespace std;

// containers    /////////////////////////////////////////////
typedef map<int, pair<int,bool> > IntBoolPair;
typedef IntBoolPair::iterator IntBoolPairItr;
//////////////////////////////////////////////////////////////

// Constants
//const double HARD_GROUNDCLAUSE_WT = DBL_MAX;
const double HARD_GROUNDCLAUSE_WT = 10;
const bool gcdebug = false;

// Forward declarations
class MLN;
class Domain;
class Clause;
class GroundPredicate;
class HashGroundPredicate;
class EqualGroundPredicate;

// Typedefs
typedef HashArray<GroundPredicate*, HashGroundPredicate, EqualGroundPredicate>
 GroundPredicateHashArray;

/**
 * Represents a grounded clause.
 */
class GroundClause
{
 public:
  GroundClause(const Clause* const & c,
               GroundPredicateHashArray* const & gndPredHashArray);

  ~GroundClause()
  {
    if (gndPredIndexes_) delete gndPredIndexes_;
    if (foClauseFrequencies_) delete foClauseFrequencies_;
  }

  void deleteFoClauseFrequencies()
  {
    if (foClauseFrequencies_) delete foClauseFrequencies_;
    foClauseFrequencies_ = NULL;
  }

  void addWt(const double& wt)
  { if (wt_ == HARD_GROUNDCLAUSE_WT) return; wt_ += wt; }

  void setWt(const double& wt)
  { if (wt_ == HARD_GROUNDCLAUSE_WT) return; wt_ = wt; }

  double getWt() const { return wt_; }

  void setWtToHardWt() { wt_ = HARD_GROUNDCLAUSE_WT; }
  bool isHardClause() const { return (wt_ == HARD_GROUNDCLAUSE_WT); }

  int getNumGroundPredicates() const { return gndPredIndexes_->size(); }

  const GroundPredicate* getGroundPredicate(const int& i,
                      GroundPredicateHashArray* const & gndPredHashArray) const
  {
    return (*gndPredHashArray)[abs((*gndPredIndexes_)[i]) - 1];
  }

  /**
   * Appends this GroundClause to all GroundPredicates in it.
   *
   * @param gndPredHashArray Reference HashArray containing the
   * GroundPredicates indexed in the GroundClause.
   */
  void appendToGndPreds(GroundPredicateHashArray* const & gndPredHashArray);

  bool getGroundPredicateSense(const int& i) const
  { return ((*gndPredIndexes_)[i] > 0); }

  void setGroundPredicateSense(const int& i, const bool& sense)
  {
      // Already the sense being set to, then return
    if ((sense && (*gndPredIndexes_)[i] > 0) ||
        (!sense && (*gndPredIndexes_)[i] < 0))
      return;

    (*gndPredIndexes_)[i] = -(*gndPredIndexes_)[i];
    rehash();
  }

  void setGroundPredicateIndex(const int& i, const int& gndPredIdx)
  { (*gndPredIndexes_)[i] = gndPredIdx; }

  int getGroundPredicateIndex(const int& i) const
  { return (*gndPredIndexes_)[i]; }

  const Array<int>* getGndPredIndexes() const
  {
    return gndPredIndexes_;
  }

  /**
   * The weight of this ground clause is set to the sum of its parent weights.
   * If the weight has been inverted from the parent, this is taken into account.
   *
   * @param mln Reference MLN to which the clause indices in foClauseFrequencies_
   * correspond.
   */
  void setWtToSumOfParentWts(const MLN* const & mln);

  IntBoolPair *getClauseFrequencies()
  {
    return foClauseFrequencies_;
  }

  int getClauseFrequency(int clauseno)
  {
	if (!foClauseFrequencies_) return 0;
	IntBoolPairItr itr = foClauseFrequencies_->find(clauseno);
	if (itr == foClauseFrequencies_->end())
	  return 0;
	else
	  return itr->second.first;
  }

  void incrementClauseFrequency(int clauseno, int increment, bool invertWt)
  {
	if (!foClauseFrequencies_)
      foClauseFrequencies_ = new IntBoolPair;
	IntBoolPairItr itr = foClauseFrequencies_->find(clauseno);
	if (itr == foClauseFrequencies_->end())
	  foClauseFrequencies_->
        insert(make_pair(clauseno, make_pair(increment, invertWt)));
	else
	  itr->second.first += increment;
  }

  /**
   * Removes a ground predicate from this ground clause. The ground predicate
   * at the given index is removed and the new hash code is stored.
   *
   * @param gndPred The ground predicate to be removed.
   */
  void removeGndPred(const int& gndPred)
  {
    for (int i = 0; i < gndPredIndexes_->size(); i++)
    {
      if (gndPred == (*gndPredIndexes_)[i])
      {
        gndPredIndexes_->removeItem(i);
        gndPredIndexes_->compress();
        rehash();
        break;
      }
    }
  }

  /**
   * Changes the index of a ground predicate in this ground clause. The new hash
   * code is stored.
   *
   * @param oldIdx Index of the ground predicate to be changed.
   * @param newIdx New index of the ground predicate.
   */
  void changeGndPredIndex(const int& oldIdx, const int& newIdx)
  {
    for (int i = 0; i < gndPredIndexes_->size(); i++)
    {
      if (oldIdx == (*gndPredIndexes_)[i])
      {
        (*gndPredIndexes_)[i] = newIdx;
        rehash();
        break;
      }
    }
  }

  size_t hashCode() { return hashCode_; }

  bool same(const GroundClause* const & gc)
  {
    if (this == gc) return true;
    if (gndPredIndexes_->size() != gc->getGndPredIndexes()->size())
    {
      return false;
    }
    return (memcmp(gndPredIndexes_->getItems(),
                   gc->getGndPredIndexes()->getItems(),
                   (gndPredIndexes_->size())*sizeof(int)) == 0);
  }

  void printWithoutWt(ostream& out) const;
  void print(ostream& out) const;

  ostream& print(ostream& out, const Domain* const& domain,
                 const bool& withWt, const bool& asInt,
                 const bool& withStrVar,
                 const GroundPredicateHashArray* const & predHashArray) const;

  ostream& printWithoutWt(ostream& out, const Domain* const & domain,
                  const GroundPredicateHashArray* const & predHashArray) const;

  ostream&
  printWithoutWtWithStrVar(ostream& out, const Domain* const & domain,
                  const GroundPredicateHashArray* const & predHashArray) const;

  ostream& printWithWtAndStrVar(ostream& out, const Domain* const& domain,
                  const GroundPredicateHashArray* const & predHashArray) const;

  ostream& print(ostream& out, const Domain* const& domain,
                 const GroundPredicateHashArray* const & predHashArray) const;

  ostream& printWithoutWtWithStrVarAndPeriod(ostream& out,
                                             const Domain* const& domain,
                  const GroundPredicateHashArray* const & predHashArray) const;

  double sizeKB();

 private:

  /**
   * Computes the hash code and stores it.
   */
  void rehash()
  {
    Array<unsigned int>* intArrRep = new Array<unsigned int>;

      // For each predicate
    for (int i = 0; i < gndPredIndexes_->size(); i++)
    {
        // For each pred 1 (if pos.) or 0 (if neg.) is appended to intArrRep
      if ((*gndPredIndexes_)[i] > 0)
        intArrRep->append(1);
      else
        intArrRep->append((unsigned int)0);
      intArrRep->append(abs((*gndPredIndexes_)[i]));
    }

    hashCode_ = Hash::hash(*intArrRep);
    delete intArrRep;
  }

 private:
    // Hash code of this ground clause
  size_t hashCode_;
  Array<int>* gndPredIndexes_; // 4 + 4*n bytes (n is no. of preds)

    // overloaded to indicate whether this is a hard clause
    // if this is a hard clause, wt_ is set to HARD_GROUNDCLAUSE_WT
  double wt_; // 8 bytes

    // Number of first-order clauses this clause corresponds to. Also stores
    // if the weight has been flipped from each parent clause
  IntBoolPair* foClauseFrequencies_;

};


////////////////////////////////// hash /////////////////////////////////


class HashGroundClause
{
 public:
  size_t operator()(GroundClause* const & gc) const  { return gc->hashCode(); }
};


class EqualGroundClause
{
 public:
  bool operator()(GroundClause* const & c1, GroundClause* const & c2) const
    { return c1->same(c2); }
};


/////////////////////////////// containers  /////////////////////////////////

typedef HashArray<GroundClause*, HashGroundClause, EqualGroundClause>
  GroundClauseHashArray;

typedef hash_set<GroundClause*, HashGroundClause, EqualGroundClause>
  GroundClauseSet;


#endif
