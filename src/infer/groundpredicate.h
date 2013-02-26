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
#ifndef GROUNDPREDICATE_H_JUN_26_2005
#define GROUNDPREDICATE_H_JUN_26_2005

#include <cmath>
#include "domain.h"
#include "predicate.h"
#include "groundclause.h"

class GroundPredicate
{
 public:
  GroundPredicate(Predicate* const & pred)
    : negGndClauses_(new Array<GroundClause*>),
      posGndClauses_(new Array<GroundClause*>),
      gndClauseSet_(new GroundClauseSet),
      truthValue_(false), wtWhenFalse_(0), wtWhenTrue_(0)
  {
    assert(pred->isGrounded());
    int numTerms = pred->getNumTerms();    
    intArrRep_ = new Array<unsigned int>(numTerms+1);
    intArrRep_->growToSize(numTerms+1);
    (*intArrRep_)[0] = pred->getId();
    assert((*intArrRep_)[0] >= 0);
    for (int i = 0; i < numTerms; i++)
    {
      (*intArrRep_)[i+1] = pred->getTerm(i)->getId();
      assert((*intArrRep_)[i+1] >= 0);
    }
    intArrRep_->compress();
    hashCode_ = Hash::hash(*intArrRep_);
  }

  
  ~GroundPredicate()
  { 
    if (intArrRep_)    delete intArrRep_; 
    if (gndClauseSet_) delete gndClauseSet_;
    delete negGndClauses_;
    delete posGndClauses_;
  }


  void deleteGndClauseSet() 
  { if (gndClauseSet_) delete gndClauseSet_; gndClauseSet_ = NULL; }

  /**
   * Removes the ground clauses from this ground predicate. The gndClauseSet_
   * is deleted and the gndClauses_ Array is reset.
   */
  void removeGndClauses() 
  {
    deleteGndClauseSet();
    delete negGndClauses_;
    delete posGndClauses_;
    negGndClauses_ = new Array<GroundClause*>;
    posGndClauses_ = new Array<GroundClause*>;
  }

  /**
   * Creates a Predicate* equivalent to this GroundPredicate.
   * The caller is responsible for deleting the returned Predicate*.
   */
  Predicate* createEquivalentPredicate(const Domain* const & domain) const
  {
    const PredicateTemplate* pt = domain->getPredicateTemplate(getId());
    assert(pt);
    Predicate* pred = new Predicate(pt);
    for (int j = 0; j < pt->getNumTerms(); j++)
      pred->appendTerm(new Term(getTermId(j), (void*)pred, true));
    return pred;
  }

  unsigned int getId() const { return (*intArrRep_)[0]; }
  unsigned int getTermId(const int& idx) const { return (*intArrRep_)[idx+1]; }
  unsigned int getNumTerms() const { return intArrRep_->size() - 1; }

  bool getTruthValue() const { return truthValue_; }
  void setTruthValue(const bool& tv) { truthValue_ = tv; }

  double getWtWhenFalse() const { return wtWhenFalse_; }
  void setWtWhenFalse(const double& wt) {wtWhenFalse_ = wt;}
  void addWtWhenFalse(const double& wt) {wtWhenFalse_ += wt;}

  double getWtWhenTrue() const { return wtWhenTrue_; }
  void setWtWhenTrue(const double& wt) {wtWhenTrue_ = wt;}
  void addWtWhenTrue(const double& wt) {wtWhenTrue_ += wt;}

    // Get the probability that the ground predicate is true
  double getProbability()
  { return 1.0 / ( 1.0 + exp(wtWhenFalse_ - wtWhenTrue_)); }


  void compress()
  {
    negGndClauses_->compress();
    posGndClauses_->compress();
  }


  size_t hashCode() { return hashCode_; }
  
  const Array<unsigned int>* getIntArrRep() const { return intArrRep_; }

    // add a copy of this GroundPredicate's intArrRep_ to rep
  void appendIntArrRep(Array<unsigned int>& rep) { rep.append(*intArrRep_); }


  bool appendGndClause(GroundClause* const & gc, const bool& senseInGndClause) 
  { 
    if (gndClauseSet_== NULL || gndClauseSet_->find(gc) == gndClauseSet_->end())
    {
      if (gndClauseSet_) gndClauseSet_->insert(gc);
      //gndClauses_->append(gc);
      if (senseInGndClause) posGndClauses_->append(gc);
      else                  negGndClauses_->append(gc);
      //senseInGndClauses_->append(senseInGndClause);
      return true;
    }
    return false;
  } 

  const Array<GroundClause*>* getNegGndClauses() const { return negGndClauses_;}
  const Array<GroundClause*>* getPosGndClauses() const { return posGndClauses_;}

  int getNumGndClauses() const 
  { return negGndClauses_->size() + posGndClauses_->size(); }


  bool same(const GroundPredicate* const & gp)
  {
    if (intArrRep_->size() != gp->getIntArrRep()->size()) return false;
    return (memcmp(intArrRep_->getItems(), gp->getIntArrRep()->getItems(), 
                   intArrRep_->size()*sizeof(unsigned int)) == 0);
  }

  string getPredicateStr(const Domain* const& domain)
  {
	  const char* predName = domain->getPredicateName((*intArrRep_)[0]);
	  string strPredName = predName;
	  strPredName = strPredName + "(";
	  int size = intArrRep_->size();
	  for (int i = 1; i < size; i++)
	  {
		  string name = domain->getConstantName((*intArrRep_)[i]);
		  string::size_type at = name.rfind("@");
		  if (at != string::npos) name = name.substr(at+1, name.length()-at-1);
		  strPredName = strPredName + name;
		  if (i < size-1) strPredName = strPredName + ",";
		  else strPredName = strPredName + ")";
	  }	
	  return strPredName;
  }
  string getPredName(const Domain* const& domain) 
  {
	  string s = string(domain->getPredicateName((*intArrRep_)[0]));
	  return s;
  }

  string getPredString(const Domain* const & domain)
  {
	  string strPred = string(domain->getPredicateName((*intArrRep_)[0]));
	  strPred = strPred + "(";
	  int size = intArrRep_->size();
	  for (int i = 1; i < size; i++)
	  {
		  string name = domain->getConstantName((*intArrRep_)[i]);
		  string::size_type at = name.rfind("@");
		  if (at != string::npos) name = name.substr(at+1, name.length()-at-1);
		  strPred = strPred + name;
		  if (i < size-1) strPred = strPred + ",";
		  else            strPred = strPred + ")";
	  }
	  
	 return strPred; 
  }

  void print(ostream& out, const Domain* const & domain) const
  {
    const char* predName = domain->getPredicateName((*intArrRep_)[0]);
    out << predName;
      // If dealing with a propositional variable
    if (strcmp(domain->getConstantName((*intArrRep_)[1]),
               Domain::PROPOSITIONAL_CONSTANT) == 0)
    {
      return;
    }    
    
    out << "(";
    int size = intArrRep_->size();
    for (int i = 1; i < size; i++)
    {
      string name = domain->getConstantName((*intArrRep_)[i]);
      string::size_type at = name.rfind("@");
      if (at != string::npos) name = name.substr(at+1, name.length()-at-1);
      out << name;
      if (i < size-1) out << ",";
      else            out << ")";
    }
  }


  void print(ostream& out) const
  {
    out << (*intArrRep_)[0] << "(";
    int size = intArrRep_->size();
    for (int i = 1; i < size; i++)
    {
      out << (*intArrRep_)[i];
      if (i < size-1) out << ",";
      else            out << ")";
    }
  }
  
  /**
   * Computes and returns the size of this ground predicate.
   */
  double sizeKB()
  {
    double size = 0;
      //intArrRep_
    if (intArrRep_)
      size += (intArrRep_->size()*sizeof(unsigned int) / 1024.0);
      // hashCode_
    size += (sizeof(size_t) / 1024.0);

    if (negGndClauses_)
    {
      for (int i = 0; i < negGndClauses_->size(); i++)
        size += ((*negGndClauses_)[i]->sizeKB());
    }
    if (posGndClauses_)
    {
      for (int i = 0; i < posGndClauses_->size(); i++)
        size += ((*posGndClauses_)[i]->sizeKB());
    }

      // truthValue_
    size += (sizeof(bool) / 1024.0);
      // wtWhenFalse_
    size += (sizeof(double) / 1024.0);
      // wtWhenTrue_
    size += (sizeof(double) / 1024.0);
    
    return size;
  }

 private:
    // intArrRep_[0] is the pred id; intArrRep_[i] is the constant id
    // where i > 0
  Array<unsigned int>* intArrRep_; // 4 bytes + 4*n bytes (n is the arity)
  size_t hashCode_; // 4 bytes
    // gnd pred is a neg lit in these clauses
  Array<GroundClause*>* negGndClauses_; // 4 bytes + 4*n bytes
    // gnd pred is a pos lit in these clauses 
  Array<GroundClause*>* posGndClauses_; // 4 bytes + 4*n bytes
    // Pointers to ground clauses in which this pred occurs
  GroundClauseSet* gndClauseSet_; // 4 bytes + 4*n bytes

    // Truth value
  bool truthValue_;
    // Wt when false
  double wtWhenFalse_;
    // Wt when true
  double wtWhenTrue_;
 
};


////////////////////////////////// hash /////////////////////////////////


class HashGroundPredicate
{
 public:
  size_t operator()(GroundPredicate* const & gp) const  {return gp->hashCode();}
};


class EqualGroundPredicate
{
 public:
  bool operator()(GroundPredicate* const & p1, GroundPredicate* const& p2) const
  { return p1->same(p2); }
};

////////////////////////////////// containers /////////////////////////////////

typedef hash_set<GroundPredicate*, HashGroundPredicate, EqualGroundPredicate> 
  GroundPredicateSet;

typedef hash_map<GroundPredicate*, int,HashGroundPredicate,EqualGroundPredicate>
  GroundPredicateToIntMap;



#endif
