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
#include "predicatetemplate.h"
#include "domain.h"

const char* PredicateTemplate::EMPTY_NAME = "SK_EMPTY";
const char* PredicateTemplate::EQUAL_NAME = "SK_EQUAL";
const char* PredicateTemplate::GT_NAME = "greaterThan";
const char* PredicateTemplate::LT_NAME = "lessThan";
const char* PredicateTemplate::GTEQ_NAME = "greaterThanEq";
const char* PredicateTemplate::LTEQ_NAME = "lessThanEq";
const char* PredicateTemplate::SUBSTR_NAME = "substr";
const char* PredicateTemplate::ANY_TYPE_NAME = "SK_ANY";
const char* PredicateTemplate::INT_TYPE_NAME = "int";
const char* PredicateTemplate::STRING_TYPE_NAME = "string";
const char* PredicateTemplate::ZZ_RETURN_PREFIX = "isReturnValueOf";

const char* PredicateTemplate::SUCC_NAME = "succ";
const char* PredicateTemplate::PLUS_NAME = "plus";
const char* PredicateTemplate::MINUS_NAME = "minus";
const char* PredicateTemplate::TIMES_NAME = "times";
const char* PredicateTemplate::DIVIDEDBY_NAME = "dividedBy";
const char* PredicateTemplate::MOD_NAME = "mod";
const char* PredicateTemplate::CONCAT_NAME = "concat";

  // Caller is responsible for deleting typeName if required.
  // returns true if termType is successfully added; false otherwise
  // isUnique is false by default 
bool PredicateTemplate::appendTermType(const char* const & typeName,
                                       const bool& isUnique,
                                       const Domain* const & domain) 
{
  int id = domain->getTypeId(typeName);
  if (id < 0) 
  {
    cout << "Warning: In PredicateTemplate::addTermType(). Type " << typeName 
         << " has not been declared." << endl;
    return false;
  }
  
  char* tname = new char[strlen(typeName)+1]; 
  strcpy(tname, typeName);
  termTypesAsStr_->append(tname);
  termTypesAsInt_->append(id);
  termsUnique_->append(isUnique);
  assert(termTypesAsStr_->size() == termTypesAsInt_->size());
  return true;
}


  // returns true if termType is successfully added; false otherwise
  // isUnique is false by default 
bool PredicateTemplate::appendTermType(const int& typeId, const bool& isUnique,
                                       const Domain* const & domain)
{
  const char* typeName = domain->getTypeName(typeId);
  if (typeName==NULL) 
  {
    cout << "Error: In PredicateTemplate::addTermType(). TypeId " << typeId 
         << " does not exist." << endl;
    return false;
  }
  
  char* tname = new char[strlen(typeName)+1]; 
  strcpy(tname, typeName);
  termTypesAsStr_->append(tname);
  termTypesAsInt_->append(typeId);
  termsUnique_->append(isUnique);
  assert(termTypesAsStr_->size() == termTypesAsInt_->size());
  return true;
}


