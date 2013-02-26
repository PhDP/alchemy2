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
#ifndef PREDICATE_H_JUN_26_2005
#define PREDICATE_H_JUN_26_2005

#include <limits>
#include <ext/hash_set>
using namespace __gnu_cxx;
#include "predicatetemplate.h"
#include "term.h"
#include "hash.h"
#include "hasharray.h"
#include "function.h"

const int MAX_VAR = 10;
enum TruthValue { FALSE = 0, TRUE = 1, UNKNOWN = 2 };

struct VarsTypeId { Array<Term*> vars; int typeId; }; 

class Clause;
class GroundPredicate;
class Domain;

  //Ensure that the dirty_ bit is consistently updated.
class Predicate
{
 public:
  Predicate(const PredicateTemplate* const & pt) 
    : template_(pt), terms_(new Array<Term*>), sense_(true), 
      truthValue_(UNKNOWN), allTermsAreDiffVars_(false), isGrounded_(false), 
      intArrRep_(NULL), hashCode_(0), dirty_(true), parent_(NULL),
      realValue_(std::numeric_limits<double>::min()) {}

  Predicate(const PredicateTemplate* const & pt, Clause* const & parent)
    : template_(pt), terms_(new Array<Term*>), sense_(true), 
      truthValue_(UNKNOWN), allTermsAreDiffVars_(false), isGrounded_(false), 
      intArrRep_(NULL), hashCode_(0), dirty_(true), parent_(parent),
      realValue_(std::numeric_limits<double>::min()) {}

  Predicate(const Predicate& p)  { parent_ = NULL; copy(p); }
  Predicate(const Predicate& p,  Clause* const & par) { parent_= par; copy(p); }

  
  ~Predicate()
  {
    for (int i = 0; i < terms_->size(); i++)  delete (*terms_)[i];
    delete terms_;
    if (intArrRep_) delete intArrRep_;
  }


  double sizeMB() const
  {
    double sizeMB = (fixedSizeB_ + intArrRep_->size()*sizeof(int) +
                     terms_->size()*sizeof(Term*))/1000000.0;
    for (int i = 0; i < terms_->size(); i++) sizeMB += (*terms_)[i]->sizeMB();
    return sizeMB;
  }


  static void computeFixedSizeB()
  { fixedSizeB_ = (sizeof(Predicate)+sizeof(Array<Term*>)+sizeof(Array<int>)); }


  void compress() 
  { 
    for (int i = 0; i < terms_->size(); i++)  (*terms_)[i]->compress();
    terms_->compress(); 
    intArrRep_->compress(); 
  }


    //returns true if this predicate represents "="
  bool isEqualPred() const  { return template_->isEqualPred(); }
  
  bool isEqualPredWithType() const  { return template_->isEqualPredWithType(); }

    //returns true if this predicate is empty predicate
  bool isEmptyPred() const  { return template_->isEmptyPredicateTemplate(); }
 
     //returns true if this predicate is an internal predicate
  bool isInternalPred() const  { return template_->isInternalPredicateTemplate(); }
  
  bool isInternalPredWithoutType() const
  { return template_->isInternalPredicateTemplateWithoutType(); }

  void canonicalize()
  {
    Array<int> varIdToNewVarId;
    int newVarId = 0;
    for (int i = 0; i < terms_->size(); i++)
    {
      if ((*terms_)[i]->getType() == Term::VARIABLE)
      {
        int id = -((*terms_)[i]->getId());
        assert(id > 0);
        if (id >= varIdToNewVarId.size()) varIdToNewVarId.growToSize(id+1,0);
          //if a new var id has not been assigned to old var id
        if (varIdToNewVarId[id] >= 0) varIdToNewVarId[id] = --newVarId;
        (*terms_)[i]->setId(varIdToNewVarId[id]);
      }
      assert((*terms_)[i]->getType() != Term::FUNCTION);
    }
  }


    //term is owned by Predicate. Caller should not delete it
  void appendTerm(Term* const & term)
  {
    if (template_->getNumTerms() == terms_->size()) 
    {
      cout << "Error: In Predicate::appendTerm. Tried to add more terms than "
           << "the declared num of " << template_->getNumTerms() << endl;
      exit(-1);
    }
    terms_->append(term);
    setDirty();
  }

  /**
   * Removes the last term in the predicate. The template is not changed and terms_
   * is compressed. Returned Term* no longer belongs to the predicate and caller is
   * responsible for deleting it.
   */
  Term* removeLastTerm()
  {
    Term* term = terms_->removeLastItem();
    terms_->compress();
    setDirty();
    return term;
  }

  
  bool getSense() const { return sense_; }
  void setSense(const bool& s) { sense_ = s; setDirty();}
  void invertSense() { sense_ = (sense_) ? false : true; setDirty(); }

    //get and set the definitive truth value of the predicate
  TruthValue getTruthValue() const { return truthValue_; }
  void setTruthValue(const TruthValue& tv) { truthValue_ = tv; }
  string getTruthValueAsStr() const 
  {
    if (truthValue_ == TRUE)    return "TRUE";
    if (truthValue_ == FALSE)   return "FALSE";
    if (truthValue_ == UNKNOWN) return "UNKNOWN";
    assert(false); 
    return "UNKNOWN"; //avoid compilation warning
  }

  double getRealValue() const { return realValue_; }
  void setRealValue(const double rv) { realValue_ = rv; }

  double getNumGroundingsIfAllVarDiff(const Domain* const & domain) const;
 
  
    //Caller is responsible for deleting the Predicate* in returnArray
  static void createAllGroundings(const int& predId, 
                                  const Domain* const & domain,
                                  Array<Predicate*>& returnArray);

  
  //Caller is responsible for deleting the Predicate* in returnArray
  static void createAllGroundingsUnifyingWithTerm(const int& predId, 
                                  const Domain* const & domain,
                                  Array<Predicate*>& returnArray,
                                  int termTypeId, int termVal);


  void createAllGroundingsIfAllVarDiff(const Domain* const & domain,
                                       Array<Predicate*>& returnArray)
  { createAllGroundings(getId(), domain, returnArray); }

     //create groundings of predicate, taking into account vars may be shared
  void createAllGroundings(const Domain* const & domain,
                           Array<Predicate*>& returnArray)
  { createAllGroundings(domain, &returnArray, NULL, -1); }


    //create groundings of predicate, taking into account vars may be shared
  void createAllGroundings(const Domain* const & domain,
                           Array<int*>& returnArray)
  { createAllGroundings(domain, NULL, &returnArray, -1); }

  /**
   * Generates exactly one grounding of this predicate.
   * 
   * @param domain Domain in which this predicate exists.
   * @param returnArray The one grounding produced is placed as the only element
   * in this Array.
   * @param grounding The index of the grounding to be produced. This is the
   * grounding-th combination of constants.
   */
  void getGroundingNumber(const Domain* const & domain,
                          Array<Predicate*>& returnArray, const int& grounding)
  { createAllGroundings(domain, &returnArray, NULL, grounding); }

  int getNumTerms() const { return terms_->size(); }

  const Term* getTerm(const int& idx) const { return (*terms_)[idx]; }


  void setTermToConstant(const int& termNum, const int& constId)
  {
    assert(termNum < template_->getNumTerms());
    assert(constId >= 0);
    if (termNum >= terms_->size()) terms_->growToSize(termNum+1,NULL);
    if ((*terms_)[termNum]) delete (*terms_)[termNum];
    (*terms_)[termNum] = new Term(constId, (void*)this, true);
    setDirty();
  }

  
  bool containsConstant(const int& constId) const
  {
    assert(constId >= 0);
    for (int i = 0; i < terms_->size(); i++)
    {
      Term* t = (*terms_)[i];
      if (t->getType() == Term::CONSTANT && t->getId() == constId) return true;
    }
    return false;
  }

  /**
   * Checks if this predicate contains at least one constant.
   * 
   * @return True if and only if at least one term is a constant
   */
  bool containsConstants() const
  {
    for (int i = 0; i < terms_->size(); i++)
    {
      Term* t = (*terms_)[i];
      if (t->getType() == Term::CONSTANT) return true;
    }
    return false;    
  }

  void setTemplate(PredicateTemplate* const & t) { template_ = t; }

  const PredicateTemplate* getTemplate() const { return template_; }


    // Caller should not delete the returned const char* 
  const char* getName() const { return template_->getName(); }

  int getId() const { return template_->getId(); }
  

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


  bool allTermsAreDiffVars()
  {
    if (dirty_) computeAndStoreIntArrRep();
    return allTermsAreDiffVars_;
  }
  

  bool checkAllTermsAreDiffVars()
  {
    Array<int> vars;
    allTermsAreDiffVars_ = true;
    for (int i = 0; i < terms_->size(); i++)
    {
      if ((*terms_)[i]->getType() == Term::VARIABLE)
      {
        int tid = (*terms_)[i]->getId();
        if (vars.contains(tid)) { allTermsAreDiffVars_ = false; return false;}
        else                    vars.append(tid);
      }
      else
      { 
        assert((*terms_)[i]->getType() != Term::FUNCTION);
          //is a constant
        allTermsAreDiffVars_ = false;
        return false;
      }
    }
    return true;
  }


  bool isGrounded()
  {
    if (dirty_) computeAndStoreIntArrRep();
    return isGrounded_;    
  }


  /**
   * Checks if this predicate can be grounded as another (potentially)
   * partially grounded predicate.
   * 
   * @param partGndPred Partially grounded predicate to which this predicate
   * is being compared.
   * @return true, if this predicate can be grounded as partGndPred.
   */  
  bool canBeGroundedAs(Predicate* const & partGndPred)
  {
    //assert(partGndPred->isGrounded());
    if (template_->getId() != partGndPred->getId()) return false;
    if (allTermsAreDiffVars()) return true;
    if (isGrounded() && partGndPred->isGrounded()) return same(partGndPred);
    
    int varGndings[MAX_VAR];
    memset(varGndings, -1, MAX_VAR*sizeof(int));
    for (int i = 0; i < terms_->size(); i++) 
    {
      int termType = (*terms_)[i]->getType();
      
      if (termType == Term::CONSTANT) 
      {
        if ((*terms_)[i]->getId() != partGndPred->getTerm(i)->getId()) 
          return false;
      }
      else 
      if (termType == Term::VARIABLE) 
      {
        int varId = -(*terms_)[i]->getId();
        assert(varId > 0);
        assert(varId < MAX_VAR);
        if (varGndings[varId] < 0) // if variable has not been grounded
          varGndings[varId] = partGndPred->getTerm(i)->getId();
        else 
          if (varGndings[varId] != partGndPred->getTerm(i)->getId())
            return false;
      }
      else
      {
        assert(false);
      }
    }
    return true;
  }


  bool canBeGroundedAs(const GroundPredicate* const & gndPred);

  Array<int> * getPredicateConstants(Array<int> * const & constants)
  {
    Array<int> * pconstants = new Array<int>();
    for (int termno = 0; termno < getNumTerms(); termno++)
    {
      const Term *term = getTerm(termno);
      int id = term->getId();
      assert(id < 0);
      pconstants->append((*constants)[-id]);
    }
    return pconstants;
  }

  Array<int> * getPredicateConstants()
  {
    Array<int> * pconstants = new Array<int>();
    for (int termno = 0; termno < getNumTerms(); termno++)
    {
      const Term *term = getTerm(termno);
      int id = term->getId();
      assert(id >= 0);
      pconstants->append(id);
    }
    return pconstants;
  }

  void setDirty();
  bool isDirty() const { return dirty_; }

  void setParent(Clause* const & parent) { parent_ = parent; setDirty(); }
  Clause* getParent() const { return parent_; }

	// True, if the predicate can be used in the inverted index
    // If in pos. clause, then it must be negated
    // If in neg. clause, then the inverted index can not be used
  bool isIndexable(bool posClause)
  {
      // At most one of the two terms is grounded
      // and equality pred is never indexed
    //if (isGrounded() || isEqualPred()) return false;
    if (isEqualPred()) return false;
    
    //return (!getSense());
    return true;
  }

  void createVarsTypeIdArr(Array<VarsTypeId*>*& varsTypeIdArr)
  {
    if (varsTypeIdArr == NULL) varsTypeIdArr = new Array<VarsTypeId*>;

      //for each variable of the predicate
    for (int j = 0; j < terms_->size(); j++)
    {
      Term* t = (*terms_)[j];
      if (t->getType() == Term::VARIABLE)
      {
        int id = -(t->getId());
        assert(id > 0);
        if (id >= varsTypeIdArr->size()) varsTypeIdArr->growToSize(id+1, NULL);
        VarsTypeId*& vti = (*varsTypeIdArr)[id];
        if (vti == NULL) 
        {
          vti = new VarsTypeId;
          vti->typeId = getTermTypeAsInt(j);
          assert(vti->typeId >= 0);
        }
        assert(getTermTypeAsInt(j) == vti->typeId);
        vti->vars.append(t);
      }
      assert(t->getType() != Term::FUNCTION);
    }// for each variable of the predicate
  }


  void deleteVarsTypeIdArr(Array<VarsTypeId*>*& varsTypeIdArr)
  {
    for (int i = 0; i < varsTypeIdArr->size(); i++)
      if ((*varsTypeIdArr)[i]) delete (*varsTypeIdArr)[i];
    delete varsTypeIdArr;
    varsTypeIdArr = NULL; 
  }


    //Does not consider the sense of the predicate.
  bool same(Predicate* const & p)
  {
    if (this == p)  return true;
    const Array<int>* pArr  = p->getIntArrRep();
    const Array<int>* myArr = getIntArrRep();
    if (myArr->size() != pArr->size()) return false;
    const int* pItems  = p->getIntArrRep()->getItems();
    const int* myItems = getIntArrRep()->getItems();
    return (memcmp(myItems, pItems, myArr->size()*sizeof(int))==0);
  }


  bool same(const GroundPredicate* const & gp);

  
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


  ostream& printAsInt(ostream& out) const
  {
    if(isEqualPred()) return printAsIntEqualPred(out);
    if(isInternalPred()) return printAsIntInternalPred(out);

    if (!sense_) out << "!";
    out << template_->getId() << "(";
    for (int i = 0; i < terms_->size(); i++)
    {
      (*terms_)[i]->printAsInt(out); 
      out << ((i!=terms_->size()-1)?",":")");
    }
    return out;
  }

  ostream& printWithStrVar(ostream& out, const Domain* const & domain) const;

  ostream& print(ostream& out, const Domain* const & domain) const;


 private:
 
 /**
   * Create groundings of predicate, taking into account vars may be shared.
   * Exactly one of predReturnArray or constReturnArray must be NULL. If
   * constReturnArray is NULL, then the generated predicates are stored in
   * predReturnArray; otherwise arrays of the constant ids are stored in
   * constReturnArray.
   * Callers are responsible for deleting the parameters and their contents.
   * 
   * @param domain Domain in which this predicate exists.
   * @param predReturnArray The groundings produced are placed in this Array.
   * @param constReturnArray The constants produced are placed in this Array.
   * @param grounding If non-negative, then only the grounding with this index
   * is produced. This is the grounding-th combination of constants starting at
   * 0.
   */
  void createAllGroundings(const Domain* const & domain,
                           Array<Predicate*>* const & predReturnArray,
                           Array<int*>* const & constReturnArray,
                           const int& grounding);
 
  void copy(const Predicate& p)
  {
    template_ = p.template_;

    terms_ = new Array<Term*>;
    Array<Term*>* tterms = p.terms_;
    for (int i = 0; i < tterms->size(); i++)
    {
      Term* t = (*tterms)[i];
      terms_->append(new Term(*t, (void*)this, true));      
    }
    dirty_ = p.dirty_;

    if (!dirty_) { assert(noDirtyTerms()); }

    sense_ = p.sense_;
    truthValue_ = p.truthValue_;
    allTermsAreDiffVars_ = p.allTermsAreDiffVars_;
    isGrounded_ = p.isGrounded_;
    
    if (p.intArrRep_)  intArrRep_ = new Array<int>(*(p.intArrRep_));
    else               intArrRep_ = NULL;    

    hashCode_ = p.hashCode_;

    setParentDirty();
  }


  void setParentDirty();

  bool noDirtyTerms()
  {
    for (int i = 0; i < terms_->size(); i++)
      if ((*terms_)[i]->isDirty()) return false;
    return true;
  }


  const Array<int>* getIntArrRep() 
  { if (dirty_) computeAndStoreIntArrRep(); return intArrRep_; }


    //Also finds out whether all of its terms are variables/constants.
  void computeAndStoreIntArrRep()
  {
    dirty_ = false;
    if (intArrRep_ == NULL) intArrRep_ = new Array<int>;
    else                    intArrRep_->clear();
      //commented out: predicate comparison ignores sense
    //intArrRep_->append(sense_?1:0);
    intArrRep_->append(template_->getId());
    int numTerms = terms_->size();
    Array<int> vars(numTerms);
    allTermsAreDiffVars_ = true;
    isGrounded_ = true;
    for (int i = 0; i < numTerms; i++)
    {
      (*terms_)[i]->appendIntArrRep(*intArrRep_);
      int termType = (*terms_)[i]->getType();
      
      if (termType == Term::VARIABLE)
      {
        if (allTermsAreDiffVars_)
        {
          int tid = (*terms_)[i]->getId();
          if (vars.contains(tid)) allTermsAreDiffVars_ = false;
          else                    vars.append(tid);
        }
        isGrounded_ = false;
      }
      else
      if (termType == Term::CONSTANT)
      {
        allTermsAreDiffVars_ = false;
        if ((*terms_)[i]->getType() != Term::CONSTANT) isGrounded_ = false;
      }
      else
      {
        assert(false);
      }
    }
    hashCode_ = Hash::hash(*intArrRep_);
  }


  ostream& printAsIntEqualPred(ostream& out) const
  {
    assert(isEqualPred());
    assert(terms_->size()==2);
    if (!sense_) out << "!(";
    (*terms_)[0]->printAsInt(out); 
    out << " = ";
    (*terms_)[1]->printAsInt(out); 
    if (!sense_) out << ")";
    return out;
  }


  ostream& printEqualPred(ostream& out, const Domain* const & domain) const
  {
    assert(isEqualPred());
    assert(terms_->size()==2);
    if (!sense_) out << "!(";
    (*terms_)[0]->print(out, domain); 
    out << " = ";
    (*terms_)[1]->print(out, domain); 
    if (!sense_) out << ")";
    return out;
  }


  ostream& 
  printEqualPredWithStrVar(ostream& out, const Domain* const & domain) const
  {
    assert(isEqualPred());
    assert(terms_->size()==2);
    if (!sense_) out << "!(";
    (*terms_)[0]->printWithStrVar(out, domain); 
    out << " = ";
    (*terms_)[1]->printWithStrVar(out, domain); 
    if (!sense_) out << ")";
    return out;
  }
    
    
  ostream& printAsIntInternalPred(ostream& out) const
  {
    assert(isInternalPred());
    assert(terms_->size()==2);
    
      //No infix
	if (strncmp(getName(), PredicateTemplate::SUBSTR_NAME,
				strlen(PredicateTemplate::SUBSTR_NAME)) == 0)
	{
	  if (!sense_) out << "!";
      out << template_->getId() << "(";
      for (int i = 0; i < terms_->size(); i++)
      {
      	(*terms_)[i]->printAsInt(out); 
        out << ((i!=terms_->size()-1)?",":")");
      }
      return out; 
	}
	else // Infix
	{
      if (!sense_) out << "!(";
      (*terms_)[0]->printAsInt(out); 
	  
	  if (strncmp(getName(), PredicateTemplate::GT_NAME,
				  strlen(PredicateTemplate::GT_NAME)) == 0)
      	out << " > ";
      else
	  if (strncmp(getName(), PredicateTemplate::LT_NAME,
				  strlen(PredicateTemplate::LT_NAME)) == 0)
      	out << " < ";
      else
      if (strncmp(getName(), PredicateTemplate::GTEQ_NAME,
				  strlen(PredicateTemplate::GTEQ_NAME)) == 0)
      	out << " >= ";
      else
	  if (strncmp(getName(), PredicateTemplate::LTEQ_NAME,
				  strlen(PredicateTemplate::LTEQ_NAME)) == 0)
      	out << " <= ";
      
      (*terms_)[1]->printAsInt(out); 
      if (!sense_) out << ")";
      return out;
	}
  }


  ostream& printInternalPred(ostream& out, const Domain* const & domain) const
  {
    assert(isInternalPred());
    assert(terms_->size()==2);

	  //No infix
	if (strncmp(getName(), PredicateTemplate::SUBSTR_NAME,
				strlen(PredicateTemplate::SUBSTR_NAME)) == 0)
	{
      if (!sense_) out << "!";
      out << PredicateTemplate::SUBSTR_NAME << "(";
      for (int i = 0; i < terms_->size(); i++)
      {
        (*terms_)[i]->print(out, domain); 
        out << ((i!=terms_->size()-1)?",":")");
      }
      return out;
	}
	else // Infix
	{
	  if (!sense_) out << "!(";
      (*terms_)[0]->print(out, domain); 
	  
	  if (strncmp(getName(), PredicateTemplate::GT_NAME,
				  strlen(PredicateTemplate::GT_NAME)) == 0)
      	out << " > ";
      else
	  if (strncmp(getName(), PredicateTemplate::LT_NAME,
				  strlen(PredicateTemplate::LT_NAME)) == 0)
      	out << " < ";
      else
      if (strncmp(getName(), PredicateTemplate::GTEQ_NAME,
				  strlen(PredicateTemplate::GTEQ_NAME)) == 0)
      	out << " >= ";
      else
	  if (strncmp(getName(), PredicateTemplate::LTEQ_NAME,
				  strlen(PredicateTemplate::LTEQ_NAME)) == 0)
      	out << " <= ";
      
      (*terms_)[1]->print(out, domain); 
      if (!sense_) out << ")";
      return out;
	}
  }


  ostream& 
  printInternalPredWithStrVar(ostream& out, const Domain* const & domain) const
  {
    assert(isInternalPred());
    assert(terms_->size()==2);

	  //No infix
	if (strncmp(getName(), PredicateTemplate::SUBSTR_NAME,
				strlen(PredicateTemplate::SUBSTR_NAME)) == 0)
	{
      if (!sense_) out << "!";
      out << PredicateTemplate::SUBSTR_NAME << "(";
      for (int i = 0; i < terms_->size(); i++)
      {
        (*terms_)[i]->printWithStrVar(out, domain); 
        out << ((i!=terms_->size()-1)?",":")");
      }
      return out;
	}
	else // Infix
	{
      if (!sense_) out << "!(";
      (*terms_)[0]->printWithStrVar(out, domain);
	  
	  if (strncmp(getName(), PredicateTemplate::GT_NAME,
				  strlen(PredicateTemplate::GT_NAME)) == 0)
      	out << " > ";
      else
	  if (strncmp(getName(), PredicateTemplate::LT_NAME,
				  strlen(PredicateTemplate::LT_NAME)) == 0)
      	out << " < ";
      else
      if (strncmp(getName(), PredicateTemplate::GTEQ_NAME,
				  strlen(PredicateTemplate::GTEQ_NAME)) == 0)
      	out << " >= ";
      else
	  if (strncmp(getName(), PredicateTemplate::LTEQ_NAME,
				  strlen(PredicateTemplate::LTEQ_NAME)) == 0)
      	out << " <= ";
      
      (*terms_)[1]->printWithStrVar(out, domain); 
      if (!sense_) out << ")";
      return out;
	}
  }
  
   
 private:
  const PredicateTemplate* template_; // not owned by Predicate
  Array<Term*>* terms_;
  bool sense_; // indicate whether the predicate is non-negated/negated
  TruthValue truthValue_;

  bool allTermsAreDiffVars_; // all of its terms are different variables
  bool isGrounded_;

  Array<int>* intArrRep_;
  size_t hashCode_;
  bool dirty_;
  Clause* parent_; //not owned by Predicate, so no deleted in destructor

  static double fixedSizeB_;
  double realValue_;
};

inline
ostream& operator<<(ostream& out, const Predicate& p){return p.printAsInt(out);}


////////////////////////////////// hash /////////////////////////////////

class HashPredicate
{
 public:
  size_t operator()(Predicate* const & p) const  { return p->hashCode(); }
};


class EqualPredicate
{
 public:
  bool operator()(Predicate* const & p1, Predicate* const & p2) 
    const { return p1->same(p2); }
};

/////////////////////////////// containers  /////////////////////////////////

typedef hash_set<Predicate*, HashPredicate, EqualPredicate> PredicateSet;

typedef HashArray<Predicate*, HashPredicate, EqualPredicate> PredicateHashArray;

#endif
