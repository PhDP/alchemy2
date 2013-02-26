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
#ifndef FUNCTIONTEMPLATE_JUN_22_2005
#define FUNCTIONTEMPLATE_JUN_22_2005

#include "predicatetemplate.h"

class FunctionTemplate : public PredicateTemplate
{
 public:
  static const char* EMPTY_FTEMPLATE_NAME; // name of empty template
 public:
  FunctionTemplate() : PredicateTemplate(), retTypeId_(-1),retTypeName_(NULL) {}
  virtual ~FunctionTemplate() { delete [] retTypeName_; }


    // Sets name and id of return type.
  void setRetTypeId(const int& id, const Domain* const & domain);

    // Sets name and id of return type.
    // Caller is responsible for deleteing name if required
  void setRetTypeName(const char* const & name, const Domain* const & domain);

  
  int getRetTypeId() const { return retTypeId_; }
  
    // Caller should not delete the returned const char*
  const char* getRetTypeName() const { return retTypeName_; }

    // returns funcName append with _typeNames
  static string createInternalFuncTypeName(const char* const & funcName,
  										   const Array<string>& typeNames)
  {
    string ftName;
    ftName.append(funcName);
    for (int i = 0; i < typeNames.size(); i++)
    {
	  ftName.append("_").append(typeNames[i]);
    }
    return ftName;
  }
  
  static bool isInternalFunctionTemplateName(const char* funcName)
  { 
  	return (strncmp(funcName, SUCC_NAME, strlen(SUCC_NAME)) == 0 || 
  			strncmp(funcName, PLUS_NAME, strlen(PLUS_NAME)) == 0 ||
  			strncmp(funcName, MINUS_NAME, strlen(MINUS_NAME)) == 0 ||
  			strncmp(funcName, TIMES_NAME, strlen(TIMES_NAME)) == 0 ||
  			strncmp(funcName, DIVIDEDBY_NAME, strlen(DIVIDEDBY_NAME)) == 0 ||
  			strncmp(funcName, MOD_NAME, strlen(MOD_NAME)) == 0 ||
  			strncmp(funcName, CONCAT_NAME, strlen(CONCAT_NAME)) == 0);
  }

  static bool isInternalFunctionUnaryTemplateName(const char* funcName)
  { 
  	return (strncmp(funcName, SUCC_NAME, strlen(SUCC_NAME)) == 0);
  }


  bool isInternalFunctionTemplate() const 
  { 
  	return (strncmp(name_, SUCC_NAME, strlen(SUCC_NAME)) == 0 || 
  			strncmp(name_, PLUS_NAME, strlen(PLUS_NAME)) == 0 ||
  			strncmp(name_, MINUS_NAME, strlen(MINUS_NAME)) == 0 ||
  			strncmp(name_, TIMES_NAME, strlen(TIMES_NAME)) == 0 ||
  			strncmp(name_, DIVIDEDBY_NAME, strlen(DIVIDEDBY_NAME)) == 0 ||
  			strncmp(name_, MOD_NAME, strlen(MOD_NAME)) == 0 ||
  			strncmp(name_, CONCAT_NAME, strlen(CONCAT_NAME)) == 0);
  }

  ostream& print(ostream& out) const
  {
  	out << retTypeName_ << " ";
    out << name_ << "(";
    for (int i = 0; i < termTypesAsStr_->size(); i++)
    {
      out << (*termTypesAsStr_)[i]; 
      out << ((i!=termTypesAsStr_->size()-1)?",":")");
    }
    return out;
  }


  ostream& printWithStrVar(ostream& out) const
  {
  	out << retTypeName_ << " ";
    out << name_ << "(";
    for (int i = 0; i < termTypesAsStr_->size(); i++)
    {
      out << "a" << i+1; 
      out << ((i!=termTypesAsStr_->size()-1)?",":")");
    }
    return out;
  }
  
  
 private:
  int retTypeId_;
  char* retTypeName_;
};

#endif
