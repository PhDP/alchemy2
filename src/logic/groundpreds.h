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
#ifndef GROUNDPREDS_H_JUN_28_2005
#define GROUNDPREDS_H_JUN_28_2005

#include "predicate.h"


class GroundPreds
{
 public:
  GroundPreds() : trueGndPreds_(new PredicateHashArray), 
                  falseGndPreds_(new PredicateHashArray),
                  unknownGndPreds_(new PredicateHashArray) {}

  ~GroundPreds()
  {
    for (int i=0; i<trueGndPreds_->size();i++)    delete (*trueGndPreds_)[i];
    for (int i=0; i<falseGndPreds_->size();i++)   delete (*falseGndPreds_)[i];
    for (int i=0; i<unknownGndPreds_->size();i++) delete (*unknownGndPreds_)[i];
    delete trueGndPreds_;
    delete falseGndPreds_;
    delete unknownGndPreds_;
  }


  bool contains(Predicate* const & p)
  {
    return (trueGndPreds_->contains(p) ||  falseGndPreds_->contains(p) || 
            unknownGndPreds_->contains(p));
  }


    // GroundPreds owns pred and is responsible for deleting it.
  void appendTrueGndPred(Predicate* const & p)  
  { p->setTruthValue(TRUE); trueGndPreds_->append(p); }

    // GroundPreds owns pred and is responsible for deleting it.
  void appendFalseGndPred(Predicate* const & p) 
  { p->setTruthValue(FALSE); falseGndPreds_->append(p);}

    // GroundPreds owns pred and is responsible for deleting it.
  void appendUnknownGndPred(Predicate* const & p)
  { p->setTruthValue(UNKNOWN); unknownGndPreds_->append(p); }


    // Caller should not delete returned array nor modify its contents.
  const PredicateHashArray* getTrueGndPreds() const  { return trueGndPreds_; }
  const PredicateHashArray* getFalseGndPreds() const { return falseGndPreds_; }
  const PredicateHashArray* getUnknownGndPreds() const 
  {return unknownGndPreds_; }


  int getNumTrueGndPreds() const    { return trueGndPreds_->size(); }
  int getNumFalseGndPreds() const   { return falseGndPreds_->size(); }
  int getNumUnknownGndPreds() const { return unknownGndPreds_->size(); }


    // caller may delete p if it is not originally obtained from GroundPreds
  void changeGndPredTruthValue(Predicate* const & pred, 
                               const TruthValue& oldValue,
                               const TruthValue& newValue)
  {
    Predicate* ppred = (Predicate*) pred;
    Predicate* p = NULL;
    
      //modified to remove the ambiguity from removeItem(const int & index)
      //when type is also an int
    if (oldValue==TRUE)
      p = trueGndPreds_->removeInputItemFastDisorder(ppred);
    else
    if (oldValue==FALSE)
      p = falseGndPreds_->removeInputItemFastDisorder(ppred);
    else 
    if (oldValue==UNKNOWN)
      p = unknownGndPreds_->removeInputItemFastDisorder(ppred);

    assert(p);
    if (newValue == TRUE)         trueGndPreds_->append(p);
    else if (newValue == FALSE)   falseGndPreds_->append(p);
    else if (newValue == UNKNOWN) unknownGndPreds_->append(p);
    else assert(false);
  }
  

  void compress()
  {
    trueGndPreds_->compress();
    falseGndPreds_->compress();
    unknownGndPreds_->compress();
  }


 private:
  PredicateHashArray* trueGndPreds_;
  PredicateHashArray* falseGndPreds_;
  PredicateHashArray* unknownGndPreds_;
  
};


#endif
