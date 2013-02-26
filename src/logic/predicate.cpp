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
#include "predicate.h"
#include "clause.h"
#include "groundpredicate.h"

double Predicate::fixedSizeB_ = -1;


void Predicate::setDirty() { dirty_ = true; if (parent_) parent_->setDirty(); }


void Predicate::setParentDirty() { if (parent_) parent_->setDirty(); }


bool Predicate::canBeGroundedAs(const GroundPredicate* const & gndPred)
{
  if ( ((unsigned int) template_->getId()) != gndPred->getId()) return false;
  if (allTermsAreDiffVars()) return true;
    
  int varGndings[MAX_VAR];
  memset(varGndings, -1, MAX_VAR*sizeof(int));
  for (int i = 0; i < terms_->size(); i++) 
  {
    int termType = (*terms_)[i]->getType();
      
    if (termType == Term::CONSTANT) 
    {
      assert((*terms_)[i]->getId() >= 0);
      if ( ((unsigned int) (*terms_)[i]->getId()) != gndPred->getTermId(i)) 
        return false;
    }
    else 
    if (termType == Term::VARIABLE) 
    {
      int varId = -(*terms_)[i]->getId();
      assert(varId > 0);
      assert(varId < MAX_VAR);
      if (varGndings[varId] < 0) // if variable has not been grounded
        varGndings[varId] = gndPred->getTermId(i);
      else 
      if (varGndings[varId] != (int) gndPred->getTermId(i))
        return false;
    }
    else
    {
      assert(false);
    }
  }
  return true;
}

  //Does not consider the sense of the predicate.
bool Predicate::same(const GroundPredicate* const & gp)
{
  if ( getId() != (int) gp->getId() ) return false;
  if ( terms_->size() != (int) gp->getNumTerms() ) return false;
  for (int i = 0; i < terms_->size(); i++)
    if ( (*terms_)[i]->getId() != (int) gp->getTermId(i) ) return false;
  return true;
}

double 
Predicate::getNumGroundingsIfAllVarDiff(const Domain* const & domain) const
{ 
  int numGnd = 1;
  for (int i = 0; i < template_->getNumTerms(); i++)
    numGnd *= domain->getNumConstantsByType(template_->getTermTypeAsInt(i));
  return numGnd;
}

    //Caller is responsible for deleting the Predicate* in returnArray
void Predicate::createAllGroundings(const int& predId, 
                                    const Domain* const & domain,
                                    Array<Predicate*>& returnArray)
{
  const PredicateTemplate* pt = domain->getPredicateTemplate(predId);
  if (pt == NULL)
  {
    cout << "ERROR in Predicate::createAllGroundings: no predicate with id " 
         << predId << " in domain" << endl;
    return;
  }
  if (pt->isEqualPredicateTemplate())
  {
    cout << "ERROR in Predicate::createAllGroundings: cannot create " 
         << "groundings for '=' predicate" << endl;
    return;      
  }
  
  ArraysAccessor<int> acc;
  for (int i = 0; i < pt->getNumTerms(); i++)
  {
    int typeId = pt->getTermTypeAsInt(i);
    assert(typeId >= 0);
    const Array<int>* constArr = domain->getConstantsByType(typeId);
    assert(constArr->size() > 0);
    acc.appendArray(constArr);
  }

  while (acc.hasNextCombination())
  {
    Predicate* pred = new Predicate(pt);
    int constId;
    while (acc.nextItemInCombination(constId))
      pred->appendTerm(new Term(constId, (void*)pred, true));
    returnArray.append(pred);    
  }
}

//get all the groundings unifying with the term given 
//Caller is responsible for deleting the Predicate* in returnArray
void Predicate::createAllGroundingsUnifyingWithTerm(const int& predId, 
                                    const Domain* const & domain,
                                    Array<Predicate*>& returnArray,
									int termTypeId,int termConstId)
{
  Array<int>* tmpArr;
  const PredicateTemplate* pt = domain->getPredicateTemplate(predId);
  if (pt == NULL)
  {
    cout << "ERROR in Predicate::createAllGroundings: no predicate with id " 
         << predId << " in domain" << endl;
    return;
  }
  if (pt->isEqualPredicateTemplate())
  {
    cout << "ERROR in Predicate::createAllGroundings: cannot create " 
         << "groundings for '=' predicate" << endl;
    return;      
  }
  
  ArraysAccessor<int> acc;
  for (int i = 0; i < pt->getNumTerms(); i++)
  {
    int typeId = pt->getTermTypeAsInt(i);
    assert(typeId >= 0);
    const Array<int>* constArr;
    if (typeId != termTypeId)
    {
      constArr = domain->getConstantsByType(typeId);
	}
    else
    {
      tmpArr = new Array<int>;
      tmpArr->append(termConstId);
      constArr = tmpArr;
	}
    assert(constArr->size() > 0);
    acc.appendArray(constArr);
  }
  
  while (acc.hasNextCombination())
  {
    Predicate* pred = new Predicate(pt);
    int constId;
    while (acc.nextItemInCombination(constId))
      pred->appendTerm(new Term(constId, (void*)pred, true));
    returnArray.append(pred);
  }
}


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
void Predicate::createAllGroundings(const Domain* const & domain,
                                    Array<Predicate*>* const & predReturnArray,
                                    Array<int*>* const & constReturnArray,
                                    const int& grounding)
{
  assert(predReturnArray == NULL || constReturnArray == NULL);

  if (isGrounded())
  {
    if (grounding >= 0) assert(grounding == 0);
    if (predReturnArray)
      predReturnArray->append(new Predicate(*this));
    else
    {
      assert(constReturnArray);
      int* constArr = new int[terms_->size()];
      for (int j = 0; j < terms_->size(); j++) 
      {
        constArr[j] = (*terms_)[j]->getId();
        assert(constArr[j] >= 0);
      }
      constReturnArray->append(constArr);
    }
    return;
  }

  Array<VarsGroundedType*> vgtArr;
  for (int i = 0; i < terms_->size(); i++)
  {
    const Term* t = (*terms_)[i];
    if (t->getType() == Term::VARIABLE)
    {
      int id = -(t->getId());
      assert(id > 0);
      if (id >= vgtArr.size()) vgtArr.growToSize(id + 1,NULL);
      VarsGroundedType*& vgt = vgtArr[id];
      if (vgt == NULL) 
      {
        vgt = new VarsGroundedType; 
        vgt->typeId = getTermTypeAsInt(i);
        assert(vgt->typeId >= 0);        
      }
      assert(getTermTypeAsInt(i) == vgt->typeId);
      vgt->vars.append((Term*)t);
    }
  }

  Array<int> negVarId; 
  ArraysAccessor<int> acc;
  for (int i = 0; i < vgtArr.size(); i++)
  {
    if (vgtArr[i] == NULL) continue;
    negVarId.append(i);
    acc.appendArray(domain->getConstantsByType(vgtArr[i]->typeId));
  }

  int combination = 0;
  while (acc.hasNextCombination())
  {
    int i = 0;
    int constId;
    while (acc.nextItemInCombination(constId))
    {
      Array<Term*>& vars = vgtArr[negVarId[i]]->vars;
      for (int j = 0; j < vars.size(); j++) vars[j]->setId(constId);
      i++;
    }
    
      // If collecting all preds or we are at the right combination
    if (grounding < 0 || grounding == combination)
    {
      if (predReturnArray)
        predReturnArray->append(new Predicate(*this));
      else
      {
        assert(constReturnArray);
        int* constArr = new int[terms_->size()];
        for (int j = 0; j < terms_->size(); j++) 
        {
          constArr[j] = (*terms_)[j]->getId();
          assert(constArr[j] >= 0);
        }
        constReturnArray->append(constArr);
      }
      if (grounding == combination) break;
    }
    combination++;
  }
  if (grounding >= 0) assert(combination == grounding);
  
    // Restore variables
  for (int i = 0; i < vgtArr.size(); i++) 
  {
    if (vgtArr[i] == NULL) continue;
    Array<Term*>& vars  = vgtArr[i]->vars;
    for (int j = 0; j < vars.size(); j++) vars[j]->setId(-i);
    delete vgtArr[i];
  }  
}


ostream& Predicate::printWithStrVar(ostream& out,
                                    const Domain* const & domain) const
{
  if (isEqualPred()) return printEqualPredWithStrVar(out,domain);
  if (isInternalPred()) return printInternalPredWithStrVar(out,domain);

  if (!sense_) out << "!";
  out << template_->getName();
    // If dealing with a propositional variable
  if (strcmp(getTermTypeAsStr(0), Domain::PROPOSITIONAL_TYPE) == 0)
    return out;

  out << "(";
  for (int i = 0; i < terms_->size(); i++)
  {
    (*terms_)[i]->printWithStrVar(out, domain); 
    out << ((i != terms_->size() - 1) ? "," : ")");
  }
  return out;
}

ostream& Predicate::print(ostream& out, const Domain* const & domain) const
{
  if (isEqualPred()) return printEqualPred(out, domain);
  if (isInternalPred()) return printInternalPred(out, domain);

  if (!sense_) out << "!";
  out << template_->getName();
    // If dealing with a propositional variable
  if (strcmp(getTermTypeAsStr(0), Domain::PROPOSITIONAL_TYPE) == 0)
    return out;

  out << "(";
  for (int i = 0; i < terms_->size(); i++)
  {
    (*terms_)[i]->print(out, domain); 
    out << ((i != terms_->size() - 1) ? "," : ")");
  }
  return out;
}

