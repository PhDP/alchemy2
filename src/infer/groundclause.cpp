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
#include "groundclause.h"
#include "groundpredicate.h"
#include "clause.h"
#include "mln.h"

GroundClause::GroundClause(const Clause* const & c, 
                           GroundPredicateHashArray* const & gndPredHashArray) 
  : wt_(c->getWt()), foClauseFrequencies_(NULL)
{
  if (gcdebug) cout << "Constructing GroundClause" << endl;
  if (gcdebug)
  {
    cout << "Seen ground preds: " << endl;
    for (int i = 0; i < gndPredHashArray->size(); i++)
    {
      (*gndPredHashArray)[i]->print(cout);
      cout << endl;
    }
  }

  int numPreds = c->getNumPredicates();
  if (gcdebug) cout << "Number of preds: " << numPreds << endl;
  gndPredIndexes_ = new Array<int>;
  Array<unsigned int>* intArrRep = new Array<unsigned int>;

    // For each predicate in clause c
  for (int i = 0; i < numPreds; i++)
  {
    if (gcdebug) cout << "Looking at pred " << i << endl;
    bool deletePred = false;
    Predicate* pred = c->getPredicate(i);

    if (gcdebug) cout << "Constructing GroundPredicate" << endl;
    GroundPredicate* gndPred = new GroundPredicate(pred);
    assert(gndPred);
    if (gcdebug)
    {
      cout << "Constructed GroundPredicate ";
      gndPred->print(cout);
      cout << endl;
    }
    
    if (gcdebug) cout << "Finding GroundPredicate" << endl;
    int index = gndPredHashArray->find(gndPred);

    if (gcdebug) cout << "Index of GroundPredicate: " << index << endl;
    if (index < 0 )
    { // Predicate not seen before: add it to hash array
      if (gcdebug) cout << "Pred " << i << " not seen before" << endl;
      index = gndPredHashArray->append(gndPred) + 1;
      if (gcdebug) cout << "Appended it" << endl;
    }
    else
    { // Predicate seen before: must be deleted later
      if (gcdebug) cout << "Pred " << i << " seen before" << endl;
      deletePred = true;
      index++;
    }
      // index is the index in hash array + 1    
    int wlit;
    if (pred->getSense())
    {
      wlit = index;
      intArrRep->append(1);
    }
    else
    {
      wlit = -index;
      intArrRep->append((unsigned int)0);
    }
    intArrRep->append(index);
    gndPredIndexes_->append(wlit);

    if (deletePred) delete gndPred;
  }
  
  hashCode_ = Hash::hash(*intArrRep);
  delete intArrRep;
  gndPredIndexes_->compress();
}

/**
 * Appends this GroundClause to all GroundPredicates in it.
 * 
 * @param gndPredHashArray Reference HashArray containing the
 * GroundPredicates indexed in the GroundClause.
 */
void GroundClause::appendToGndPreds(
                      GroundPredicateHashArray* const & gndPredHashArray)
{ 
    // For each ground pred in this clause
  for (int i = 0; i < gndPredIndexes_->size(); i++)
  {
    bool sense = ((*gndPredIndexes_)[i] > 0);
    int index = abs((*gndPredIndexes_)[i]) - 1;
      // Tell the ground pred that it occurs in this ground clause
    (*gndPredHashArray)[index]->appendGndClause(this, sense);
    //assert(ok); ok = true;  // avoid compilation warning
  }
}

/**
 * The weight of this ground clause is set to the sum of its parent weights.
 * If the weight has been inverted from the parent, this is taken into account.
 * 
 * @param mln Reference MLN to which the clause indices in foClauseFrequencies_
 * correspond.
 */
void GroundClause::setWtToSumOfParentWts(const MLN* const & mln)
{
  wt_ = 0;

  IntBoolPairItr itr;
  for (itr = foClauseFrequencies_->begin();
       itr != foClauseFrequencies_->end(); itr++)
  {
    int clauseno = itr->first;
    int frequency = itr->second.first;
    bool invertWt = itr->second.second;
    double parentWeight = mln->getClause(clauseno)->getWt();
    if (invertWt) wt_ -= parentWeight*frequency;
    else wt_ += parentWeight*frequency;
  }
}

void GroundClause::printWithoutWt(ostream& out) const
{
  for (int i = 0; i < gndPredIndexes_->size(); i++)
  {
    out << (*gndPredIndexes_)[i];
    if (i < gndPredIndexes_->size() - 1) out << " v ";
  }
}

void GroundClause::print(ostream& out) const
{ out << wt_ << " "; printWithoutWt(out); }

ostream&
GroundClause::print(ostream& out, const Domain* const& domain, 
                    const bool& withWt, const bool& asInt, 
                    const bool& withStrVar,
                   const GroundPredicateHashArray* const & predHashArray) const
{
  if (withWt) out << wt_ << " ";

  Array<Predicate*> eqPreds;
  Array<Predicate*> internalPreds;
  for (int i = 0; i < gndPredIndexes_->size(); i++)
  {
    int index = (*gndPredIndexes_)[i];
    Predicate* pred =
      (*predHashArray)[abs(index)-1]->createEquivalentPredicate(domain);
    //Predicate* pred = 
      //new Predicate((*(*predHashArray)[abs(index)-1]));
    assert(pred);
    if (index < 0) pred->setSense(false);
    else pred->setSense(true);

    if (pred->isEqualPred()) 
      eqPreds.append(pred);
    else
    if (pred->isInternalPred()) 
      internalPreds.append(pred);
    else
    {
      if (asInt)           pred->printAsInt(out);
      else if (withStrVar) pred->printWithStrVar(out, domain);
      else                 pred->print(out,domain);
      if (i < gndPredIndexes_->size() - 1 || !eqPreds.empty() ||
          !internalPreds.empty()) out << " v ";
      delete pred;
    }
  }

  for (int i = 0; i < eqPreds.size(); i++)
  {
    if (asInt)           eqPreds[i]->printAsInt(out);
    else if (withStrVar) eqPreds[i]->printWithStrVar(out,domain);
    else                 eqPreds[i]->print(out,domain);
    out << ((i != eqPreds.size()-1 || !internalPreds.empty())?" v ":"");      
    delete eqPreds[i];
  }

  for (int i = 0; i < internalPreds.size(); i++)
  {
    if (asInt)           internalPreds[i]->printAsInt(out);
    else if (withStrVar) internalPreds[i]->printWithStrVar(out,domain);
    else                 internalPreds[i]->print(out,domain);
    out << ((i!=internalPreds.size()-1)?" v ":"");      
    delete internalPreds[i];
  }

  return out;
}


ostream&
GroundClause::printWithoutWt(ostream& out, const Domain* const & domain,
                   const GroundPredicateHashArray* const & predHashArray) const
{ return print(out, domain, false, false, false, predHashArray); }

ostream& 
GroundClause::printWithoutWtWithStrVar(ostream& out,
                                       const Domain* const & domain,
                   const GroundPredicateHashArray* const & predHashArray) const
{ return print(out, domain, false, false, true, predHashArray); }

ostream&
GroundClause::printWithWtAndStrVar(ostream& out, const Domain* const& domain,
                   const GroundPredicateHashArray* const & predHashArray) const
{ return print(out, domain, true, false, true, predHashArray); }

ostream&
GroundClause::print(ostream& out, const Domain* const& domain,
                   const GroundPredicateHashArray* const & predHashArray) const
{ return print(out, domain, true, false, false, predHashArray); }
    
ostream&
GroundClause::printWithoutWtWithStrVarAndPeriod(ostream& out, 
                                                const Domain* const& domain,
                   const GroundPredicateHashArray* const & predHashArray) const
{
  printWithoutWtWithStrVar(out, domain, predHashArray);
  if (isHardClause()) out << ".";
  return out;
}


/**
 * Computes and returns the size of this ground clause.
 */
double GroundClause::sizeKB()
{
  double size = 0;
    // gndPredIndexes_
  if (gndPredIndexes_)
    size += (gndPredIndexes_->size()*sizeof(int) / 1024.0);
    // wt_
  size += (sizeof(double) / 1024.0);
    // foClauseFrequencies_
  if (foClauseFrequencies_)
    size += (foClauseFrequencies_->size()*
             (2*sizeof(int) + sizeof(bool)) / 1024.0);
  
  return size;
}

