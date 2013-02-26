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
#include "clause.h"
#include "domain.h"
#include "database.h"
#include "mln.h"
#include "truefalsegroundingsstore.h"

const char* Domain::PROPOSITIONAL_TYPE = "AlchemyPropositionalType";
const char* Domain::PROPOSITIONAL_CONSTANT = "AlchemyPropositionalConstant";

Domain::~Domain()
{
  if (typeDualMap_) delete typeDualMap_;
  if (constDualMap_) delete constDualMap_;
  
  if (strToPredTemplateMap_)
  {
      // strToPredTemplateMap_'s keys (predicate names) are shared with 
      // predDualMap_. The latter will delete the keys.
    StrToPredTemplateMap::iterator it = strToPredTemplateMap_->begin(); 
    for (; it != strToPredTemplateMap_->end(); it++)
      delete (*it).second; //delete PredTemplate*;
    delete strToPredTemplateMap_;
  }
  if (predDualMap_) delete predDualMap_;

  if (strToFuncTemplateMap_)
  {
      // strToFuncTemplateMap_'s keys (function names) are shared with 
      // funcDualMap_. The latter will delete the keys.
    StrToFuncTemplateMap::iterator it2 = strToFuncTemplateMap_->begin(); 
    for (; it2 != strToFuncTemplateMap_->end(); it2++)
      delete (*it2).second; //delete FuncTemplate*;
    delete strToFuncTemplateMap_;
  }
  if (funcDualMap_) delete funcDualMap_;
  
  if (equalPredTemplate_) delete equalPredTemplate_;
  
  if (emptyPredTemplate_) delete emptyPredTemplate_;
  
  if (emptyFuncUnaryTemplate_) delete emptyFuncUnaryTemplate_;
  
  if (emptyFuncBinaryTemplate_) delete emptyFuncBinaryTemplate_;
  
  if (constantsByType_)
  {
    for (int i = 0; i < constantsByType_->size(); i++)
      delete (*constantsByType_)[i];
    delete constantsByType_;
  }

  if (externalConstantsByType_)
  {
    for (int i = 0; i < externalConstantsByType_->size(); i++)
      delete (*externalConstantsByType_)[i];
    delete externalConstantsByType_;
  }

  if (predBlocks_)
  {
    for (int i = 0; i < predBlocks_->size(); i++)
      delete (*predBlocks_)[i];
    delete predBlocks_;
  }

  if (truePredsInBlock_)
  {
    for (int i = 0; i < truePredsInBlock_->size(); i++)
      delete (*truePredsInBlock_)[i];
    delete truePredsInBlock_;
  }

  if (blockSizes_) delete blockSizes_;

  if (blockEvidence_) delete blockEvidence_;

  if (db_) delete db_;

  if (trueFalseGroundingsStore_) delete trueFalseGroundingsStore_;
  
  if (funcSet_)
  {
    FunctionSet::iterator fit;
    while (!funcSet_->empty())
    { 
      fit = funcSet_->begin();
      funcSet_->erase(fit);
      delete *fit;
    }
    delete funcSet_;
  }

  //if (externalConstant_) delete externalConstant_;
  
  if (numNonEvidAtomsPerPred_) delete numNonEvidAtomsPerPred_;
  if (numTrueNonEvidGndingsPerClause_) delete numTrueNonEvidGndingsPerClause_;
  if (numFalseNonEvidGndingsPerClause_) delete numFalseNonEvidGndingsPerClause_;
}


void Domain::deleteDB() { if (db_) delete db_; db_ = NULL; }


void Domain::newTrueFalseGroundingsStore() 
{ trueFalseGroundingsStore_ = new TrueFalseGroundingsStore(this); }

  
void Domain::compress() 
{
  typeDualMap_->compress();
  constDualMap_->compress();
  predDualMap_->compress();
  funcDualMap_->compress();
  constantsByType_->compress();
  for (int i = 0; i < constantsByType_->size(); i++)
    (*constantsByType_)[i]->compress();
  externalConstantsByType_->compress();
  for (int i = 0; i < externalConstantsByType_->size(); i++)
    (*externalConstantsByType_)[i]->compress();
  db_->compress();
}

/**
 * The constants of this domain are replaced by new constants. The new set of
 * constants must be a superset of the original one. New constants are marked
 * as external.
 */
void Domain::updatePerOldToNewIds(MLN* const & mln,
                                  hash_map<int,int> & oldToNewConstIds)
{
    //ensure that the constants in MLN clauses have the new ids 
  const ClauseHashArray* clauses = mln->getClauses();
  for (int i = 0; i < clauses->size(); i++)
  {
    Clause* c = (*clauses)[i];
    for (int i = 0; i < c->getNumPredicates(); i++)
      changePredTermsToNewIds(c->getPredicate(i), oldToNewConstIds);
  }

  const Array<Clause*>* hclauses = mln->getHybridClauses();
  for (int i = 0; i < hclauses->size(); i++)
  {
    Clause* c = (*hclauses)[i];
    for (int i = 0; i < c->getNumPredicates(); i++)
      changePredTermsToNewIds(c->getPredicate(i), oldToNewConstIds);
  }

    // Change the const ids in the pred blocks
  for (int i = 0; i < predBlocks_->size(); i++)
  {
    changePredTermsToNewIds((*predBlocks_)[i], oldToNewConstIds);
  }
  for (int i = 0; i < truePredsInBlock_->size(); i++)
  {
    changePredTermsToNewIds((*truePredsInBlock_)[i], oldToNewConstIds);
  }

    // Change the const ids in the database
  if (db_)
  {
    db_->changeConstantsToNewIds(oldToNewConstIds);
  }
}

// for parser
void Domain::reorderConstants(MLN* const & mln,
                           hash_map<int, PredicateHashArray*>& predIdToPredsMap)
{
  hash_map<int,int> oldToNewConstIds;

  ConstDualMap* newConstDualMap = new ConstDualMap;
  Array<Array<int>*>* newConstantsByType = new Array<Array<int>*>;
  Array<Array<int>*>* newExternalConstantsByType = new Array<Array<int>*>;

  int prevNewConstId = -1;
  bool constChanged = false;
  for (int i = 0; i < constantsByType_->size(); i++)
  {
    newConstantsByType->append(new Array<int>);
    Array<int>* constIds = (*constantsByType_)[i];
    for (int j = 0; j < constIds->size(); j++)
    {
      int constId = (*constIds)[j];
      const char* constName = getConstantName(constId);
      int newConstId = newConstDualMap->insert(constName, i);
      prevNewConstId = newConstId;
      (*newConstantsByType)[i]->append(newConstId);
      oldToNewConstIds[constId] = newConstId;
      if (constId != newConstId) constChanged = true;
    }
  }
  for (int i = 0; i < externalConstantsByType_->size(); i++)
  {
    newExternalConstantsByType->append(new Array<int>);
    Array<int>* constIds = (*externalConstantsByType_)[i];
    for (int j = 0; j < constIds->size(); j++)
    {
      int constId = (*constIds)[j];
      const char* constName = getConstantName(constId);
      int newConstId = newConstDualMap->insert(constName, i);
      prevNewConstId = newConstId;
      (*newExternalConstantsByType)[i]->append(newConstId);
      oldToNewConstIds[constId] = newConstId;
      if (constId != newConstId) constChanged = true;
    }
  }

  if (!constChanged)
  {
    delete newConstDualMap;
    for (int i = 0; i < newConstantsByType->size(); i++)
      delete (*newConstantsByType)[i];
    delete newConstantsByType;
    for (int i = 0; i < newExternalConstantsByType->size(); i++)
      delete (*newExternalConstantsByType)[i];
    delete newExternalConstantsByType;
    return;
  }

  delete constDualMap_;
  for (int i = 0; i < constantsByType_->size(); i++)
    delete (*constantsByType_)[i];
  delete constantsByType_;
  for (int i = 0; i < externalConstantsByType_->size(); i++)
    delete (*externalConstantsByType_)[i];
  delete externalConstantsByType_;

  constDualMap_ = newConstDualMap;
  constantsByType_ = newConstantsByType;
  externalConstantsByType_ = newExternalConstantsByType;

  constantsByType_->compress();
  for (int i = 0; i < constantsByType_->size(); i++)
    (*constantsByType_)[i]->compress();
  externalConstantsByType_->compress();
  for (int i = 0; i < externalConstantsByType_->size(); i++)
    (*externalConstantsByType_)[i]->compress();
  constDualMap_->compress();

    // update
  updatePerOldToNewIds(mln, oldToNewConstIds);

	// predIdToPredMap
  hash_map<int, PredicateHashArray*>::iterator mit = predIdToPredsMap.begin();
  for (; mit != predIdToPredsMap.end(); mit++)
  {
    PredicateHashArray* pha = (*mit).second;
    for (int i = 0; i < pha->size(); i++ )
      changePredTermsToNewIds((*pha)[i], oldToNewConstIds);
  }

    //clauses and hence their hash values have changed, and they need to be 
    //inserted into the MLN
  mln->rehashClauses();
}

  // for domain0 - postprocess
void Domain::reorderConstants(MLN* const & mln)
{
  hash_map<int,int> oldToNewConstIds;

  ConstDualMap* newConstDualMap = new ConstDualMap;
  Array<Array<int>*>* newConstantsByType = new Array<Array<int>*>;
  Array<Array<int>*>* newExternalConstantsByType = new Array<Array<int>*>;

  int prevNewConstId = -1;
  bool constChanged = false;
  assert(constantsByType_->size() == externalConstantsByType_->size());
  for (int i = 0; i < constantsByType_->size(); i++)
  {
    newConstantsByType->append(new Array<int>);
    Array<int>* constIds = (*constantsByType_)[i];
    for (int j = 0; j < constIds->size(); j++)
    {
      int constId = (*constIds)[j];
      const char* constName = getConstantName(constId);
      int newConstId = newConstDualMap->insert(constName, i);
      prevNewConstId = newConstId;
      (*newConstantsByType)[i]->append(newConstId);
      oldToNewConstIds[constId] = newConstId;
      if (constId != newConstId) constChanged = true;
    }
    newExternalConstantsByType->append(new Array<int>);
    Array<int>* extConstIds = (*externalConstantsByType_)[i];
    for (int j = 0; j < extConstIds->size(); j++)
    {
      int constId = (*extConstIds)[j];
      const char* constName = getConstantName(constId);
      int newConstId = newConstDualMap->insert(constName, i);
      prevNewConstId = newConstId;
      (*newExternalConstantsByType)[i]->append(newConstId);
      oldToNewConstIds[constId] = newConstId;
      if (constId != newConstId) constChanged = true;
    }
  }

  if (!constChanged)
  {
    delete newConstDualMap;
    for (int i = 0; i < newConstantsByType->size(); i++)
      delete (*newConstantsByType)[i];
    delete newConstantsByType;
    for (int i = 0; i < newExternalConstantsByType->size(); i++)
      delete (*newExternalConstantsByType)[i];
    delete newExternalConstantsByType;
    return;
  }

  delete constDualMap_;
  for (int i = 0; i < constantsByType_->size(); i++)
    delete (*constantsByType_)[i];
  delete constantsByType_;
  for (int i = 0; i < externalConstantsByType_->size(); i++)
    delete (*externalConstantsByType_)[i];
  delete externalConstantsByType_;

  constDualMap_ = newConstDualMap;
  constantsByType_ = newConstantsByType;
  externalConstantsByType_ = newExternalConstantsByType;

  constantsByType_->compress();
  for (int i = 0; i < constantsByType_->size(); i++)
    (*constantsByType_)[i]->compress();
  externalConstantsByType_->compress();
  for (int i = 0; i < externalConstantsByType_->size(); i++)
    (*externalConstantsByType_)[i]->compress();
  constDualMap_->compress();

    // update
  updatePerOldToNewIds(mln, oldToNewConstIds);

    //clauses and hence their hash values have changed, and they need to be 
    //inserted into the MLN
  mln->rehashClauses();
}

  // for domain1-N - postprocess
  // Assumption is: map is a superset of this domain's ConstDualMap
void Domain::reorderConstants(ConstDualMap* const & map,
                              Array<Array<int>*>* const & cbt,
                              Array<Array<int>*>* const & ecbt,
                              MLN* const & mln)
{
    // Generate oldtoNew
  hash_map<int,int> oldToNewConstIds;

    // Mark constants not currently in domain as external
  for (int i = 0; i < externalConstantsByType_->size(); i++)
    delete (*externalConstantsByType_)[i];
  delete externalConstantsByType_;
  externalConstantsByType_ = new Array<Array<int>*>;
  for (int i = 0; i < cbt->size(); i++)
    externalConstantsByType_->append(new Array<int>);

  int numConstants = map->getNumInt();
  for (int i = 0; i < numConstants; i++)
  {
    const char* name = map->getStr(i);
    if (isConstant(name)) 
    {
        // Constant is present: note the new id to replace later
      oldToNewConstIds[constDualMap_->getInt(name)] = i;
    }
    else
    {
/*
      bool continueOuter = false;
        // Add as external constant
        // Constant is in reference domain
      for (int j = 0; j < cbt->size(); j++)
      {
        bool breakOuter = false;
        for (int k = 0; k < (*cbt)[j]->size(); k++)
        {
          if (map->getInt(name) == (*(*cbt)[j])[k])
          {
              // Add as external constant
            (*externalConstantsByType_)[j]->append((*(*cbt)[j])[k]);
            breakOuter = true;
            continueOuter = true;
            break;
          }
        }
        if (breakOuter) break;
      }
      if (continueOuter) continue;
        // Constant is external in reference domain (comes from another domain)
      for (int j = 0; j < ecbt->size(); j++)
      {
        bool breakOuter = false;
        for (int k = 0; k < (*ecbt)[j]->size(); k++)
        {
          if (map->getInt(name) == (*(*ecbt)[j])[k])
              // Add as external constant
          {
            (*externalConstantsByType_)[j]->append((*(*ecbt)[j])[k]);
            breakOuter = true;
            break;
          }
        }
        if (breakOuter) break;
      }
*/
      int mi = map->getInt(name);
      Array<int>* ints = map->getInt2(mi);
      for (int i = 0; i < ints->size(); i++)
        (*externalConstantsByType_)[(*ints)[i]]->append(mi);
    }
  }
  
    // Replace old ids with new ids in constantsByType_
  for (int i = 0; i < constantsByType_->size(); i++)
  {
    for (int j = 0; j < (*constantsByType_)[i]->size(); j++)
    {
      int oldId = (*(*constantsByType_)[i])[j];
      assert(oldToNewConstIds.find(oldId) != oldToNewConstIds.end());
      (*(*constantsByType_)[i])[j] = oldToNewConstIds[oldId];
    }
    (*constantsByType_)[i]->quicksort();
  }

    // swap map/type
  delete constDualMap_;  
  constDualMap_ = map;

    // update
  updatePerOldToNewIds(mln, oldToNewConstIds);

    //clauses and hence their hash values have changed, and they need to be 
    //inserted into the MLN
  mln->rehashClauses();
}


void Domain::changePredTermsToNewIds(Predicate* const & p,
                                     hash_map<int,int>& oldToNewConstIds)
{
    // p could be NULL
  if (p)
  {
    for (int j = 0; j < p->getNumTerms(); j++)
    {
      Term* t = (Term*) p->getTerm(j);
      if (t->getType() == Term::CONSTANT)
      {
        int oldId = t->getId();
        assert(oldToNewConstIds.find(oldId) != oldToNewConstIds.end());
        t->setId(oldToNewConstIds[oldId]);
      } 
    }
  }
} 


  //Caller is responsible for deleting returned pointer
Predicate* Domain::createPredicate(const int& predId, 
                                   const bool& includeEqualPreds) const
{
  const PredicateTemplate* pt = getPredicateTemplate(predId);
  if (!includeEqualPreds && pt->isEqualPredicateTemplate()) return NULL;
  Predicate* pred = new Predicate(pt);
  pred->setSense(true);
  for (int j = 0; j < pt->getNumTerms(); j++)
    pred->appendTerm(new Term(-(j+1), (void*)pred, true));
  return pred;
}


  //Caller is responsible for deleting Array and its contents
void Domain::createPredicates(Array<Predicate*>* const & preds,
                              const bool& includeEqualPreds) const
{
  for (int i = 0; i < getNumPredicates(); i++)
  {
    Predicate* p = createPredicate(i, includeEqualPreds);
    if (p) preds->append(p);
  }
}


  //Caller is responsible for deleting Array and its contents
void Domain::createPredicates(Array<Predicate*>* const & preds,
                              const Array<string>* const & predNames)
{
  for (int i = 0; i < predNames->size(); i++)
  {
    int predId = getPredicateId((*predNames)[i].c_str());    
    Predicate* p = createPredicate(predId, true);
    if (p) preds->append(p);
  }
}

int Domain::getNumNonEvidenceAtoms() const
{
  int numAtoms = 0;
  for (int i = 0; i < getNumPredicates(); i++)
  {
    numAtoms += (*numNonEvidAtomsPerPred_)[i];
  }
  return numAtoms;
}

/*
 * Caller is responsible for deleting returned Predicate* if necessary
 */
Predicate* Domain::getNonEvidenceAtom(const int& index) const
{
  int predId = -1;
  int numAtomsPerPred;
  int numAtoms = 0;
  for (int i = 0; i < getNumPredicates(); i++)
  {
    numAtomsPerPred = (*numNonEvidAtomsPerPred_)[i];
    if (numAtoms + numAtomsPerPred >= index + 1)
    {
      predId = i;
      break;
    }
    numAtoms += numAtomsPerPred;
  }
  assert(predId >= 0);

    // Get the newIndex-th grounding of f.o. pred with id predId   
  Predicate* pred = createPredicate(predId, false);
    // Not all groundings of pred are non-evidence, so we need the while loop
  bool foundNE = false;
  while(!foundNE)
  {
    for (int i = 0; i < pred->getNumTerms(); i++)
    {
      int termType = pred->getTermTypeAsInt(i);
      const Array<int>* constantsByType = getConstantsByType(termType);
      int constIdx = random() % constantsByType->size();
      pred->setTermToConstant(i, (*constantsByType)[constIdx]);
    }
    assert(pred->isGrounded());
    if (!db_->getEvidenceStatus(pred)) foundNE = true;
  }
  return pred;
}

int Domain::addPredBlock(Predicate* const & predBlock) const
{
  int idx = predBlocks_->append(predBlock);
    // Record number of groundings
  Array<Predicate*> predArr;
  predBlock->createAllGroundings(this, predArr);
  blockSizes_->append(predArr.size());
  predArr.deleteItemsAndClear();

    // If nothing known, then no pred is set to true and no evidence in block
  Predicate* gndPred = NULL;
  truePredsInBlock_->append(gndPred);
  blockEvidence_->append(false);
  return idx;
}

/**
 * Returns the index of the block which a ground predicate
 * is in. If not in any, returns -1.
 * 
 * @param pred Predicate whose block is being searched for.
 * @return Index of the block found, or -1 if none.
 */
int Domain::getBlock(const Predicate* const & pred) const
{
  assert(((Predicate*)pred)->isGrounded());
  for (int i = 0; i < predBlocks_->size(); i++)
  {
    if ((*predBlocks_)[i]->canBeGroundedAs((Predicate*)pred))
      return i;
  }
  return -1;
}

/**
 * Returns the index of the block which a ground predicate
 * is in. If not in any, returns -1.
 * 
 * @param pred Predicate whose block is being searched for.
 * @return Index of the block found, or -1 if none.
 */
int Domain::getBlock(const GroundPredicate* const & pred) const
{
  for (int i = 0; i < predBlocks_->size(); i++)
  {
    if ((*predBlocks_)[i]->canBeGroundedAs(pred))
      return i;
  }
  return -1;
}

int Domain::getNumPredBlocks() const
{
  return predBlocks_->size();
}

const Array<Predicate*>* Domain::getPredBlocks() const
{
  return predBlocks_;
}

const Predicate* Domain::getPredBlock(const int& index) const
{
  return (*predBlocks_)[index];
}

const Array<bool>* Domain::getBlockEvidenceArray() const
{
  return blockEvidence_;
}

const bool Domain::getBlockEvidence(const int& index) const
{
  return (*blockEvidence_)[index];
}

void Domain::setBlockEvidence(const int& index, const bool& value) const
{
  (*blockEvidence_)[index] = value;
}

/**
 * Computes the number of non-evidence atoms for each first-order predicate.
 */
void Domain::computeNumNonEvidAtoms()
{
  if (numNonEvidAtomsPerPred_) delete numNonEvidAtomsPerPred_;
  numNonEvidAtomsPerPred_ = new Array<int>(getNumPredicates());
  numNonEvidAtomsPerPred_->growToSize(getNumPredicates());
  for (int i = 0; i < getNumPredicates(); i++)
  {
    (*numNonEvidAtomsPerPred_)[i] =
      (db_->getNumGroundings(i) - db_->getNumEvidenceGndPreds(i));
  }
}

