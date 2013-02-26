/*
 * All of the documentation and software included in the
 * Alchemy Software is copyrighted by Stanley Kok, Parag
 * Singla, Matthew Richardson, Pedro Domingos, Marc
 * Sumner, Hoifung Poon, and Daniel Lowd.
 * 
 * Copyright [2004-08] Stanley Kok, Parag Singla, Matthew
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
//NOTE: A better name for this class is needed. Variable may be a little
//misleading
#ifndef VARIABLE_H_JAN_2008
#define VARIABLE_H_JAN_2008

#include "util.h"
#include "predicate.h"


/**
 * This class represents a predicate term (stores the predicate id and term
 * no inside the predicate)
 */
class PredicateTerm
{
 public:
  PredicateTerm(int & predId, int & termno)
  {
    predId_ = predId;
    termno_ = termno;
  }

  int getPredId() { return predId_;}

  int getTermno() { return termno_;}

 private:
  int predId_;
  int termno_;
};


/**
 * This class stores all the information about a variable appearing a predicate
 * in a class. Also, stores the variable equivalence class id (defined by the
 * tied(v1,v2) relation to which this variable belongs.
 */
class Variable
{
 public:
  Variable(Clause* const & clause, int varId, Predicate* const & pred, 
           int& termno, int& eqClassId, const Array<int>* const & constants)
  {
    clause_ = clause;
    pred_ = pred;
    varId_ = varId;
    termno_ = termno;
    eqClassId_ = eqClassId; 
    tiedVars_ = new Array<Variable*>();
    tiedVars_->append(this);
    representative_ = true;
    implicitConstants_ = new IntHashArray();

      //all constants are initially implicit
    for (int i = 0; i < constants->size(); i++)
    {
      implicitConstants_->append((*constants)[i]);
    }
  }

  Array<Variable*>* getTiedVariables() {return tiedVars_;} 
          
  Clause * getClause(){return clause_;}
  int getPredId() {return pred_->getId();}
  int getVarId() {return varId_;}
  int getTermno() {return termno_;}
  int getEqClassId() {return eqClassId_;}
  bool isRepresentative() {return representative_;}
  int getNumTiedVariables() {return tiedVars_->size();}
  Variable * getTiedVariable(int pos) {return (*tiedVars_)[pos];}

    //update the clausevarids tied to this object
  void setTiedVars(Array<Variable *> * const & newTiedVars)
  {
    tiedVars_ = newTiedVars;
  }

  void setEqClassId(int id) { eqClassId_ = id;}

  void setRepresentative(bool b) { representative_ = b;}

  bool hasSameEqClass(Variable * const & var)
  {
    return (eqClassId_ == var->getEqClassId());
  }

  bool same(Variable * const & var)
  {
    return
      //either they are the same variable (not the constant) in a clause
	(clause_->same(var->getClause()) &&
     varId_ == var->getVarId() && varId_ < 0) || 
      //or some variable in a predicate
    ((pred_->getId() == var->getPredId()) && (termno_ == var->getTermno()));
  }

    //merge two Variable objects (and the Vars tied to them)
  void merge(Variable * const & var)
  {
    if (hasSameEqClass(var))
      return;
    tiedVars_->append(var->tiedVars_);
    delete var->tiedVars_;

    Variable *linkedVar;
      //copy over stuff to all the variables linked to his variable
    for (int i = 0; i < tiedVars_->size(); i++)
    {
      linkedVar = (*tiedVars_)[i];
      linkedVar->setEqClassId(eqClassId_);
      linkedVar->setTiedVars(tiedVars_);
      linkedVar->setRepresentative(false);
    }
    representative_ = true;
  }

  void print(ostream& out, const Domain* const & domain)
  {
    out<<"Class id: "<<eqClassId_<<endl;
    out<<pred_->getName()<<" :"<<termno_<<endl;
     cout<<endl;
  }

  void printAll(ostream& out, const Domain* const & domain)
  {
    for (int i = 0; i < tiedVars_->size(); i++)
    {
      (*tiedVars_)[i]->print(out,domain);
    }
  }

/****************************************************************************/
// Functions relating to manipulation of constants
/****************************************************************************/

  bool isImplicit(int constantId)
  {
    return implicitConstants_->find(constantId) >= 0;
  }

  int getNumImplicitConstants() { return implicitConstants_->size();}

  int getImplicitIndex(int constantId)
  {
    return implicitConstants_->find(constantId);
  }

    //remove this constant from the list of implicit constants
  void removeImplicit(int constantId)
  {
    int index = implicitConstants_->find(constantId);
    if (index >= 0)
      implicitConstants_->removeItem(index);
  }

  int getImplicitConstant(int index)
  {
    assert(implicitConstants_->size() > index);
    return (*implicitConstants_)[index]; 
  }

    //first constant of the implicit array is taken as the place holder for
    //the whole set
  bool isImplicitPlaceHolder(int constantId)
  {
    if (implicitConstants_->size() <= 0)
      return false;
    return constantId == (*implicitConstants_)[0];
  }

  void printImplicitConstants(ostream & out, Domain * const & domain)
  {
    for (int i = 0; i < implicitConstants_->size(); i++)
    {
      out<<(*implicitConstants_)[i]<<" ";
    }
  }

 private:
  Clause *clause_;
  Predicate *pred_;
  int varId_;
  int termno_;
  Array<Variable *> *tiedVars_;
  int eqClassId_;
  bool representative_;

    //constants which are represented implicitly
  IntHashArray * implicitConstants_;
};


class HashPredicateTerm
{
 public:
  size_t operator()(PredicateTerm *pt) const
  {
    return 31*pt->getPredId()+ pt->getTermno();
  }
};


class EqualPredicateTerm
{
 public:
  bool operator()(PredicateTerm * const & pt1,
                  PredicateTerm * const & pt2) const
  {
    return (pt1->getPredId() == pt2->getPredId()) &&
           (pt1->getTermno() == pt2->getTermno());
  }
};

#endif

