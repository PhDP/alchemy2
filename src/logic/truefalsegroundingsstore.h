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
#ifndef TRUEFALSEGROUNDINGSSTORE_H_NOV_16_2005
#define TRUEFALSEGROUNDINGSSTORE_H_NOV_16_2005

#include <ext/hash_map>
using namespace __gnu_cxx;
#include "random.h"
#include "domain.h"


struct TrueFalseGndings 
{
  TrueFalseGndings() {}
  ~TrueFalseGndings() 
  {
    for (int i = 0; i < trueGndings.size(); i++) delete [] trueGndings[i];
    for (int i = 0; i < falseGndings.size(); i++) delete [] falseGndings[i];
  }

  Array<int*> trueGndings; 
  Array<int*> falseGndings; 
};


typedef hash_map<Predicate*, TrueFalseGndings*, HashPredicate,EqualPredicate>
  PredToGndingsMap;



class TrueFalseGroundingsStore
{
 public:
  TrueFalseGroundingsStore(const Domain* const & domain) : domain_(domain) 
  {
    predToGndingsMaps_.growToSize(domain_->getNumPredicates(), NULL);
    random_.init(-1);
  }


  ~TrueFalseGroundingsStore() 
  {
    for (int i = 0; i < predToGndingsMaps_.size(); i++)
    {
      PredToGndingsMap* ptgMap = predToGndingsMaps_[i];
      if (ptgMap == NULL) continue;
      Array<Predicate*> toBeDeleted;      
      for (PredToGndingsMap::iterator it=ptgMap->begin();it!=ptgMap->end();it++)
      {
        toBeDeleted.append((*it).first);
        delete (*it).second;
      }
      delete ptgMap;
      toBeDeleted.deleteItemsAndClear();
    }
  }


    //Caller is responsible for deleting returned predicate
  int getNumTrueLiteralGroundings(const Predicate* const & pred,
                                  const Predicate* const & flippedGndPred)
  {
    if (((Predicate*)pred)->isGrounded())
    {
      TruthValue tv = domain_->getDB()->getValue(pred);
      if (Database::sameTruthValueAndSense(tv, pred->getSense())) return 1;
      return 0;
    }

    TrueFalseGndings* tfGndings = getTrueFalseGndings(pred, flippedGndPred);
    bool sense = pred->getSense();
    int n = (sense) ?  tfGndings->trueGndings.size() 
                     : tfGndings->falseGndings.size();

      // check whether flippedGndPred should be included in numGndings
    if (flippedGndPred && 
        ((Predicate*)pred)->canBeGroundedAs((Predicate*)flippedGndPred))
    {
      TruthValue flippedtv = domain_->getDB()->getValue(flippedGndPred);
      TruthValue tv = FALSE;
      if (flippedtv == TRUE)       tv = FALSE;
      else if (flippedtv == FALSE) tv = TRUE;
      else assert(false);
      if (Database::sameTruthValueAndSense(tv, pred->getSense())) n--;
    }
    return n;
  }


    //Caller is responsible for deleting returned predicate and the input params
  Predicate* getRandomTrueLiteralGrounding(const Predicate* const & pred,
                                           const Predicate* const & flippedPred)
  {
    if (((Predicate*)pred)->isGrounded())
    {
      TruthValue tv = domain_->getDB()->getValue(pred);
      if (Database::sameTruthValueAndSense(tv, pred->getSense())) 
        return new Predicate(*pred);
      return NULL;
    }


    TrueFalseGndings* tfGndings = getTrueFalseGndings(pred, flippedPred);
    Predicate* randPred = new Predicate(*pred);

    bool sense = pred->getSense();
    bool firstTry = true;
    if (sense)
    {
      int numGndings = tfGndings->trueGndings.size();
      while (firstTry || 
             (flippedPred && randPred->same((Predicate*)flippedPred)))
      {
        firstTry = false;
        int* constIds=tfGndings->trueGndings[ random_.randomOneOf(numGndings) ];
        for (int i = 0; i < randPred->getNumTerms(); i++)
          ((Term*)randPred->getTerm(i))->setId( constIds[i] );
      }
      return randPred;
    }
    else       
    {
      int numGndings = tfGndings->falseGndings.size();
      while (firstTry || 
             (flippedPred && randPred->same((Predicate*)flippedPred)))
      {
        firstTry = false;
        int* constIds=tfGndings->falseGndings[ random_.randomOneOf(numGndings)];
        for (int i = 0; i < randPred->getNumTerms(); i++)
          ((Term*)randPred->getTerm(i))->setId( constIds[i] );
      }
      return randPred;
    }
  }


 private:
  void fillInTrueFalseGndings(const Predicate* const & pred, 
                              TrueFalseGndings* const & tfGndings,
                              const Predicate* const & flippedPred)
  {
    bool flip = false;
    if (flippedPred && 
        ((Predicate*)pred)->canBeGroundedAs((Predicate*)flippedPred))
    { flip = true; domain_->getDB()->flipValue((Predicate*)flippedPred); }

    Predicate* tmpPred = new Predicate(*pred);
    Array<int*> constArrays;
    tmpPred->createAllGroundings(domain_, constArrays);
    for (int i = 0; i < constArrays.size(); i++)
    {
      int* constIds = constArrays[i];
      for (int j = 0; j < tmpPred->getNumTerms(); j++)
        ((Term*)tmpPred->getTerm(j))->setId(constIds[j]);

      TruthValue tv = domain_->getDB()->getValue(tmpPred); 
      if (tv == TRUE) tfGndings->trueGndings.append(constIds);
      else if (tv == FALSE) tfGndings->falseGndings.append(constIds);
    }
    tfGndings->trueGndings.compress();
    tfGndings->falseGndings.compress();
    delete tmpPred;

    if (flip) domain_->getDB()->flipValue((Predicate*)flippedPred);
  }


    //Caller is responsible for deleting pred
  TrueFalseGndings* newEntryInPredToGndingsMaps(const Predicate* const & pred,
                                         const Predicate* const & flippedPred)
  {
    PredToGndingsMap* ptgMap = predToGndingsMaps_[pred->getId()];
    TrueFalseGndings* tfGndings = new TrueFalseGndings;
    Predicate* newPred = new Predicate(*pred);
    (*ptgMap)[newPred] = tfGndings;
    fillInTrueFalseGndings(newPred, tfGndings, flippedPred);
    return tfGndings;
  }


  TrueFalseGndings* getTrueFalseGndings(const Predicate* const & pred,
                                        const Predicate* const & flippedPred)
  {
    int predId = pred->getId();
    if (predToGndingsMaps_[predId] == NULL) 
      predToGndingsMaps_[predId] = new PredToGndingsMap;
    PredToGndingsMap* ptgMap = predToGndingsMaps_[predId];

    PredToGndingsMap::iterator it;
    if ((it=ptgMap->find((Predicate*)pred)) == ptgMap->end())
      return newEntryInPredToGndingsMaps(pred, flippedPred);
    return (*it).second;
  }

  
 private:
  const Domain* domain_; // not owned by TrueFalseGroundingsStore; don't delete

    //PredToGndingsMaps[p] is the PredToGndingsMap of predicate with id p
  Array<PredToGndingsMap*> predToGndingsMaps_;
  Random random_;

};

#endif
