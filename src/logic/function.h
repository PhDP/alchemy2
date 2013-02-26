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
#ifndef FUNCTION_H_JUN_26_2005
#define FUNCTION_H_JUN_26_2005

#include <ext/hash_set>
using namespace __gnu_cxx;
#include "functiontemplate.h"
#include "hash.h"
#include "term.h"


  //Ensure that the dirty_ bit is consistently updated.
class Function
{
 public:

  Function(const FunctionTemplate* const & ft) 
    : template_(ft), terms_(new Array<Term*>), retConstId_(-1), 
      intArrRep_(NULL), hashCode_(0), dirty_(true), parent_(NULL) {}


  Function(const FunctionTemplate* const & ft, Term* const & parent) 
    : template_(ft), terms_(new Array<Term*>), retConstId_(-1), 
    intArrRep_(NULL), hashCode_(0), dirty_(true), parent_(NULL) {}
    

  Function(const Function& f) { parent_ = NULL;  copy(f); }
  Function(const Function& f, Term* const & parent) { parent_=parent; copy(f); }
  

  ~Function() 
  {
    for (int i = 0; i < terms_->size(); i++)
      delete (*terms_)[i];
    delete terms_;
    if (intArrRep_) delete intArrRep_;
  }


    //term is owned by Function. Caller should not delete it.
  void appendTerm(Term* const & term)
  {
    if (template_->getNumTerms() == terms_->size()) 
    {
      cout << "Error: In Function::appendTerm. Tried to add more terms than "
           << "the declared num of " << template_->getNumTerms() << endl;
      exit(-1);
    }
    terms_->append(term);
    setDirty();
  }

  /**
   * Removes the last term in the function. The template is not changed and terms_
   * is compressed. Returned Term* no longer belongs to the function and caller is
   * responsible for deleting it.
   */
  Term* removeLastTerm()
  {
    Term* term = terms_->removeLastItem();
    terms_->compress();
    setDirty();
    return term;
  }
  
  int getNumTerms() const { return terms_->size(); }

  const Term* getTerm(const int& idx) const { return (*terms_)[idx]; }

  void setTemplate(FunctionTemplate* const & t) { template_ = t; }

  const FunctionTemplate* getTemplate() const { return template_; }

    // Caller should not delete the returned const char* 
  const char* getName() const { return template_->getName(); }

  int getId() const { return template_->getId(); }


    // retConstId refers to the id of the constant returned by the function
  void setRetConstId(const int& constId) { retConstId_ = constId; setDirty();}
  int getRetConstId() const { return retConstId_; }

  void setDirty();
  bool isDirty() const { return dirty_; }

  void setParent(Term* const parent) { parent_ = parent; setDirty(); }
  Term* getParent() const { return parent_; }


    // Caller should not delete the returned const char* nor modify its
    // contents. Returns NULL if idx is larger than the possible number of 
    // terms
  const char* getTermTypeAsStr(const int& idx) const 
  {
    if (idx >= template_->getNumTerms()) return NULL;
    return template_->getTermTypeAsStr(idx); 
  }

    // Returns -1 if idx is larger than the possible number of terms
  int getTermTypeAsInt(const int& idx) const 
  {
    if (idx >= template_->getNumTerms()) return -1;
    return template_->getTermTypeAsInt(idx); 
  }


  int getRetTypeId() const { return template_->getRetTypeId(); }
  
    // Caller should not delete the returned const char*
  const char* getRetTypeName() const { return template_->getRetTypeName(); }

  bool same(const Function* const & f) const
  {
    if (this == f) return true;
    //commented out: no need to check return values are the same
    //if (retConstId_ != f->getRetConstId()) return false;
    if (template_->getId() != f->getId()) return false;
    int numTerms = terms_->size();
    if (numTerms != f->getNumTerms()) return false;
    for (int i = 0; i < numTerms; i++)
      if ((*terms_)[i]->getId() != f->getTerm(i)->getId()) return false;
    return true;
  }

  
    // represent function as an Array<int> and append it to rep
  void appendIntArrRep(Array<int>& rep) 
  {
    if (dirty_) computeAndStoreIntArrRep();
    rep.append(*intArrRep_);
  }


  size_t hashCode()
  {
    if (dirty_) computeAndStoreIntArrRep();
    return hashCode_;
  }


  bool same(Function* const & f)
  {
    if (this==f)  return true;
    const Array<int>* fArr  = f->getIntArrRep();
    const Array<int>* myArr = getIntArrRep();
    if (myArr->size() != fArr->size()) return false;
    const int* fItems  = f->getIntArrRep()->getItems();
    const int* myItems = getIntArrRep()->getItems();
    return (memcmp(myItems, fItems, myArr->size()*sizeof(int))==0); 
  }


  ostream& printAsInt(ostream& out) const
  {
    out << template_->getId() << "(";
    for (int i = 0; i < terms_->size(); i++)
    {
      (*terms_)[i]->printAsInt(out); 
      out << ((i!=terms_->size()-1)?",":")");
    }
    return out;
  }


  ostream& printAsIntWithRetConstId(ostream& out) const
  {
    printAsInt(out);
    out << "=" << retConstId_ << endl;
    return out;
  }


  ostream& print(ostream& out, const Domain* const & domain) const
  {
    out << template_->getName() << "(";
    for (int i = 0; i < terms_->size(); i++)
    {
      (*terms_)[i]->print(out, domain); 
      out << ((i!=terms_->size()-1)?",":")");
    }
    return out;
  }


  ostream& 
  printWithRetConstName(ostream& out, const Domain* const & domain) const;


 private:
  void copy(const Function& f)
  {
    template_ = f.template_;

    terms_ = new Array<Term*>;
    Array<Term*>* fterms = f.terms_;
    for (int i = 0; i < fterms->size(); i++)
    {
      Term* t = (*fterms)[i];
      terms_->append(new Term(*t, (void*)this, false));
    }
    dirty_ = f.dirty_;

    if (!dirty_) { assert(noDirtyTerms()); }

    retConstId_ = f.retConstId_;

    if (f.intArrRep_)  intArrRep_ = new Array<int>(*(f.intArrRep_));
    else               intArrRep_ = NULL;

    hashCode_ = f.hashCode_; 

    if (parent_) parent_->setDirty();
  }

  
  bool noDirtyTerms()
  {
    for (int i = 0; i < terms_->size(); i++)
      if ((*terms_)[i]->isDirty()) return false;
    return true;
  }


  const Array<int>* getIntArrRep() 
  { if (dirty_) computeAndStoreIntArrRep(); return intArrRep_; }


  void computeAndStoreIntArrRep()
  {
      //note that retConstId_ is not used to identify a function
    dirty_ = false;
    if (intArrRep_ == NULL) intArrRep_ = new Array<int>;
    else                    intArrRep_->clear();
    intArrRep_->append(template_->getId());
    
    int numTerms = terms_->size();
    for (int i = 0; i < numTerms; i++)
      (*terms_)[i]->appendIntArrRep(*intArrRep_);

    intArrRep_->compress();
    hashCode_ = Hash::hash(*intArrRep_);
  }


 private:
  const FunctionTemplate* template_; // not owned by Function
  Array<Term*>* terms_;
  int retConstId_;
  Array<int>* intArrRep_;
  size_t hashCode_;
  bool dirty_;
  Term* parent_;  // not owned by Function, so do not delete in destructor
};

inline
ostream& operator<<(ostream& out, const Function& f) {return f.printAsInt(out);}


////////////////////////////////// hash  /////////////////////////////////

class HashFunction
{
 public:
  size_t operator()(Function* const & f) const  { return f->hashCode(); }
};


class EqualFunction
{
 public:
  bool operator()(const Function* const & f1, const Function* const & f2) const
  { return f1->same(f2); }
};


///////////////////////////// container ////////////////////////////////////

typedef hash_set<Function*, HashFunction, EqualFunction> FunctionSet;


#endif
