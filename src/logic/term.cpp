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
#include "term.h"
#include "function.h"
#include "domain.h"
#include <cstring>

double Term::fixedSizeB_ = -1;


Term::Term(const Term& t)
{
  parent_ = NULL;
  parentIsPred_ = true;
  copy(t);
}


Term::Term(const Term& t, void* const & parent, const bool& parentIsPred)
{
  parent_ = parent;
  parentIsPred_ = parentIsPred;
  copy(t);
}


Term::~Term() 
{ 
  if (function_) delete function_; 
  if (intArrRep_) delete intArrRep_; 
}


void Term::copy(const Term& t)
{
  id_ = t.id_;

  if (t.function_)  function_ = new Function(*(t.function_), this);
  else              function_ = NULL;

  if (t.intArrRep_)  intArrRep_ = new Array<int>(*(t.intArrRep_));
  else               intArrRep_ = NULL;  
  
  dirty_ = t.dirty_;

  if (!dirty_) { assert(noDirtyFunction()); }

  if (parent_ != NULL)
  {
    if (parentIsPred_) ((Predicate*)parent_)->setDirty();
    else               ((Function*)parent_)->setDirty();
  }
}


bool Term::noDirtyFunction()
{
  if (function_ == NULL) return true;
  return !(function_->isDirty());
}


void Term::setDirty() 
{ 
  dirty_ = true;
  if (parent_ != NULL)
  {
    if (parentIsPred_) ((Predicate*)parent_)->setDirty();
    else               ((Function*)parent_)->setDirty();
  }
}


void Term::computeAndStoreIntArrRep()
{
  dirty_ = false;

  if (intArrRep_ == NULL) intArrRep_ = new Array<int>;
  else                    intArrRep_->clear();
  
  if (function_ == NULL)
    intArrRep_->append(id_);
  else
  {
    if (id_ >= 0) intArrRep_->append(id_);
    else          function_->appendIntArrRep(*intArrRep_);
  }
}


void Term::printAsInt(ostream& out) const
{
  if (function_ == NULL) 
    out << id_;
  else 
  { //function_ != NULL;
    if (id_ >= 0) out << id_;
    else          function_->printAsInt(out);
  }
}


void Term::printWithStrVar(ostream& out, const Domain* const & domain) const
{
  if (function_ == NULL)
  {
    if (id_ < 0) 
      out << "a" << -id_;  // variable
    else
    {
      string cn = domain->getConstantName(id_);
      string::size_type at = cn.rfind("@");
      if (at != string::npos) cn = cn.substr(at+1, cn.length()-at-1);
      out << cn; // constant
    }
  }
  else
  { //function_ != NULL;
    if (id_ >= 0) 
    {
      string cn = domain->getConstantName(id_);
      string::size_type at = cn.rfind("@");
      if (at != string::npos) cn = cn.substr(at+1, cn.length()-at-1);
      out << cn; // return value is known
    }
    else          
      function_->print(out,domain); 
  }
}


void Term::print(ostream& out, const Domain* const & domain) const
{
  if (function_ == NULL)
  {
    if (id_ < 0) out << id_;  // variable
    else         
    {
      string cn = domain->getConstantName(id_);
      string::size_type at = cn.rfind("@");
      if (at != string::npos) cn = cn.substr(at+1, cn.length()-at-1);
      out << cn; // constant
    }
  }
  else
  { //function_ != NULL;
    if (id_ >= 0) 
    {
      string cn = domain->getConstantName(id_);
      string::size_type at = cn.rfind("@");
      if (at != string::npos) cn = cn.substr(at+1, cn.length()-at-1);
      out << cn; // return value is known
    }
    else          
      function_->print(out, domain); 
  }
}
