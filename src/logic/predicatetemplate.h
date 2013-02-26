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
#ifndef PREDICATETEMPLATE_JUN_22_2005
#define PREDICATETEMPLATE_JUN_22_2005

#include <ctype.h>
#include "array.h"
#include <cstring>

class Domain;


class PredicateTemplate
{
 public:
  static const char* EMPTY_NAME; // name of empty template
  //Internal predicate names
  static const char* EQUAL_NAME; // name of template representing '='
  static const char* GT_NAME; // name of template representing '>'
  static const char* LT_NAME; // name of template representing '<'
  static const char* GTEQ_NAME; // name of template representing '>='
  static const char* LTEQ_NAME; // name of template representing '<='
  static const char* SUBSTR_NAME; // name of template representing 'substring'
  static const char* ANY_TYPE_NAME;   // name of any type
  static const char* INT_TYPE_NAME;   // name of type representing integers
  static const char* STRING_TYPE_NAME;   // name of type representing strings
  // Start of the predicate name used for replacing functions
  static const char* ZZ_RETURN_PREFIX;
  
  // Internal function names
  static const char* SUCC_NAME; // name of template representing 'succ'
  static const char* PLUS_NAME; // name of template representing '+'
  static const char* MINUS_NAME; // name of template representing '-'
  static const char* TIMES_NAME; // name of template representing '*'
  static const char* DIVIDEDBY_NAME; // name of template representing '/'
  static const char* MOD_NAME; // name of template representing '%'
  static const char* CONCAT_NAME; // name of template representing 'concat'
 private:
  enum { NONE = 0, EQUAL_PRED = 1, EQUAL_PRED_WITH_TYPE = 2 };
  enum { INT_PRED = 1, INT_PRED_WITH_TYPE = 2 };

 public:
  PredicateTemplate() : id_(-1), name_(NULL), termTypesAsInt_(new Array<int>), 
                        termTypesAsStr_(new Array<const char*>),
                        termsUnique_(new Array<bool>), numGnd_(-1.0), 
                        uniqueVarIndexes_(new Array<int>),
                        equalPred_(NONE), intPred_(NONE)
  {}
  
  virtual ~PredicateTemplate() 
  { 
    delete [] name_;
    for (int i = 0; i < termTypesAsStr_->size(); i++)
      delete [] (*termTypesAsStr_)[i];
    delete termTypesAsInt_; 
    delete termTypesAsStr_; 
    delete termsUnique_;
    delete uniqueVarIndexes_;
  }


    // returns EQUAL_NAME append with _typeName
  static string createEqualPredTypeName(const char* const & typeName)
  {
    string ptName;
    ptName.append(PredicateTemplate::EQUAL_NAME).append("_").append(typeName);
    return ptName;
  }

    // returns predName append with _typeName
  static string createInternalPredTypeName(const char* const & predName,
  										   const char* const & typeName)
  {
    string ptName;
    ptName.append(predName).append("_").append(typeName);
    return ptName;
  }
  
  	// returns array with internal predicate names appended with _typeName
  static Array<string> createInternalPredTypeNames(const char* const & typeName)
  {
  	Array<string> predTypeNames;
  	
  	string gtName;
  	gtName.append(PredicateTemplate::GT_NAME).append("_").append(typeName);
  	predTypeNames.append(gtName);
  	
  	string ltName;
  	ltName.append(PredicateTemplate::LT_NAME).append("_").append(typeName);
  	predTypeNames.append(ltName);
  	
  	string gteqName;
  	gteqName.append(PredicateTemplate::GTEQ_NAME).append("_").append(typeName);
  	predTypeNames.append(gteqName);
  	
  	string lteqName;
  	lteqName.append(PredicateTemplate::LTEQ_NAME).append("_").append(typeName);
  	predTypeNames.append(lteqName);
  	
  	string substrName;
  	substrName.append(PredicateTemplate::SUBSTR_NAME).append("_").append(typeName);
  	predTypeNames.append(substrName);
  	
  	return predTypeNames;
  }

  static bool isEqualPredName(const char* predName)
  { return (strncmp(predName, EQUAL_NAME, strlen(EQUAL_NAME))==0); }

  static bool isInternalPredicateTemplateName(const char* predName)
  { 
  	return (strncmp(predName, GT_NAME, strlen(GT_NAME)) == 0 || 
  			strncmp(predName, LT_NAME, strlen(LT_NAME)) == 0 ||
  			strncmp(predName, GTEQ_NAME, strlen(GTEQ_NAME)) == 0 ||
  			strncmp(predName, LTEQ_NAME, strlen(LTEQ_NAME)) == 0 ||
  			strncmp(predName, SUBSTR_NAME, strlen(SUBSTR_NAME)) == 0);
  }


  bool isEqualPredicateTemplate() const 
  { return (equalPred_ == EQUAL_PRED || equalPred_ == EQUAL_PRED_WITH_TYPE); }

  bool isEmptyPredicateTemplate() const 
  { return (strcmp(name_, EMPTY_NAME) == 0); }

  bool isInternalPredicateTemplate() const 
  { 
  	return (strncmp(name_, GT_NAME, strlen(GT_NAME)) == 0 || 
  			strncmp(name_, LT_NAME, strlen(LT_NAME)) == 0 ||
  			strncmp(name_, GTEQ_NAME, strlen(GTEQ_NAME)) == 0 ||
  			strncmp(name_, LTEQ_NAME, strlen(LTEQ_NAME)) == 0 ||
  			strncmp(name_, SUBSTR_NAME, strlen(SUBSTR_NAME)) == 0);
  }

  bool isInternalPredicateTemplateWithoutType() const 
  { 
  	return (strcmp(name_, GT_NAME) == 0 || 
  			strcmp(name_, LT_NAME) == 0 ||
  			strcmp(name_, GTEQ_NAME) == 0 ||
  			strcmp(name_, LTEQ_NAME) == 0 ||
  			strcmp(name_, SUBSTR_NAME) == 0);
  }

  bool isPredicateTemplateFromInternalFunction() const
  {
  	return (strcmp(name_, string(ZZ_RETURN_PREFIX).append(string(SUCC_NAME)).c_str()) == 0 ||
  			strcmp(name_, string(ZZ_RETURN_PREFIX).append(string(PLUS_NAME)).c_str()) == 0 ||
  			strcmp(name_, string(ZZ_RETURN_PREFIX).append(string(MINUS_NAME)).c_str()) == 0 ||
  			strcmp(name_, string(ZZ_RETURN_PREFIX).append(string(TIMES_NAME)).c_str()) == 0 ||
  			strcmp(name_, string(ZZ_RETURN_PREFIX).append(string(DIVIDEDBY_NAME)).c_str()) == 0 ||
  			strcmp(name_, string(ZZ_RETURN_PREFIX).append(string(MOD_NAME)).c_str()) == 0 ||
  			strcmp(name_, string(ZZ_RETURN_PREFIX).append(string(CONCAT_NAME)).c_str()) == 0);
  }
  
  bool isPredicateTemplateFromFunction() const
  {
  	if (strlen(name_) > strlen(ZZ_RETURN_PREFIX))
  	  return (strncmp(name_, ZZ_RETURN_PREFIX, strlen(ZZ_RETURN_PREFIX)) == 0);
  	return false;
  }
  
  static string translateEqualPredicateName(const string& eqPredName)
  {
    if (strncmp(eqPredName.c_str(), "same", 4) != 0) return "";
    string typeStr = eqPredName.substr(4, eqPredName.length()-4);
    string newName = PredicateTemplate::EQUAL_NAME;
    newName.append("_").append(typeStr);
    return newName;
  }


  void setId(const int& id)  { id_ = id; }

  int getId() const { return id_; }


    // Caller is responsible for deleting name if required.
  void setName(const char* const & name) 
  { 
    if (name_) delete [] name_;
    name_ = new char[strlen(name)+1]; 
    strcpy(name_, name);

    unsigned int len = strlen(EQUAL_NAME);
    if (strncmp(name_, EQUAL_NAME, len)==0 && strlen(name_) > len)
      equalPred_ = EQUAL_PRED_WITH_TYPE;
    else
    if (strncmp(name_,EQUAL_NAME, len)==0) 
      equalPred_ = EQUAL_PRED;
    else 
      equalPred_ = NONE;
  }

  
  bool isEqualPred() const 
  { return (equalPred_ == EQUAL_PRED || equalPred_ == EQUAL_PRED_WITH_TYPE); }


  bool isEqualPredWithType() const
  { return (equalPred_ == EQUAL_PRED_WITH_TYPE); }

    // Caller should not delete returned const char*.
  const char* const getName() const { return name_; }


  int getNumTerms() const 
  {
    assert(termTypesAsStr_->size() == termTypesAsInt_->size());
    return termTypesAsInt_->size(); 
  }


    // Caller is responsible for deleting typeName if required.
    // returns true if termType is successfully added; false otherwise
  bool appendTermType(const char* const & typeName, const bool& isUnique,
                      const Domain* const & domain);

    // returns true if termType is successfully added; false otherwise
  bool appendTermType(const int& typeId, const bool& isUnique, 
                      const Domain* const & domain);
                      

    // Caller should not delete the returned Array<const char&>* nor modify its
    // contents.
  const Array<const char*>* getTermTypesAsStr() const 
  { return termTypesAsStr_; }


    // Caller should not delete the returned const char* nor modify its
    // contents.
  const char* getTermTypeAsStr(const int& idx) const 
  { return (*termTypesAsStr_)[idx]; }


    // Caller should not delete the returned Array<int>* nor modify its
    // contents.
  const Array<int>* getTermTypesAsInt() const 
  { return termTypesAsInt_; }


  int getTermTypeAsInt(const int& idx) const 
  { return (*termTypesAsInt_)[idx]; }
  

    // Caller should not delete the returned Array<bool>* nor modify its
    // contents.
  const Array<bool>* getTermsUnique() const { return termsUnique_; }


  const bool termIsUnique(const int& idx) const 
  { return (*termsUnique_)[idx]; }

    // Caller should not delete the returned Array<bool>* nor modify its
    // contents.
  const Array<int>* getUniqueVarIndexes() const { return uniqueVarIndexes_; }

  /**
   * Marks a term as a unique variable index (a blocked variable).
   * 
   * @param idx Index of the term to be marked as a unique variable index.
   */
  void addUniqueVarIndex(const int& idx) const 
  {
    assert(idx < getNumTerms());
    uniqueVarIndexes_->append(idx);
  }

  virtual ostream& print(ostream& out) const
  {
    out << name_ << "(";
    for (int i = 0; i < termTypesAsStr_->size(); i++)
    {
      out << (*termTypesAsStr_)[i];
      if (uniqueVarIndexes_->contains(i))
        out << "!";
      out << ((i!=termTypesAsStr_->size()-1)?",":")");
    }
    return out;
  }


  virtual ostream& printWithStrVar(ostream& out) const
  {
    out << name_;
      // If dealing with a propositional variable
    if (strcmp(getTermTypeAsStr(0), "AlchemyPropositionalType") == 0)
    {
      return out;
    }
    
    out << "(";
    for (int i = 0; i < termTypesAsStr_->size(); i++)
    {
      out << "a" << i+1; 
      if (uniqueVarIndexes_->contains(i))
        out << "!";
      out << ((i != termTypesAsStr_->size() - 1) ? "," : ")");
    }
    return out;
  }


 protected:
  int id_;     // id of predicate
  char* name_; // name of predicate
  Array<int>* termTypesAsInt_;
  Array<const char*>* termTypesAsStr_;
  Array<bool>* termsUnique_;
  double numGnd_;
    // Indices of blocked variables
  Array<int>* uniqueVarIndexes_;
  
 private:
  int equalPred_;
  int intPred_;
};


inline
ostream& operator<<(ostream& o, const PredicateTemplate& p) {return p.print(o);}


#endif
