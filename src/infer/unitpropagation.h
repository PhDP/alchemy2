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
#ifndef UNITPROPAGATION_H_
#define UNITPROPAGATION_H_

#include "inference.h"

const int updebug = false;

class UnitPropagation : public Inference
{
 public:
  
  /**
   * There is nothing to do in the constructor for unit propagation.
   */
  UnitPropagation(VariableState* state, long int seed,
                  const bool& trackClauseTrueCnts)
    : Inference(state, seed, trackClauseTrueCnts) {}

  ~UnitPropagation() {}
  
  /**
   * There is nothing to initialize in unit propagation.
   */
  void init() { return; }
  
  /**
   * Perform unit propagation on the clauses in the state.
   */
  void infer()
  {
    if (updebug) cout << "Entering UnitPropagation::infer" << endl;
    bool done = false;
      // We are done when no more propagation occurs
    while (!done)
    {
      if (updebug) cout << endl << endl;
      done = true;
      for (int i = 0; i < state_->getNumClauses(); i++)
      {
          // Don't look at dead clauses
        if (state_->isDeadClause(i)) continue;
          // Simplify clause i using the atoms already fixed
        Array<int>* atoms = state_->simplifyClauseFromFixedAtoms(i);
        if (atoms->size() == 0)
        {
          delete atoms;
          continue;
        }
        if (updebug) cout << "Clause " << i << " ";
          // If unit clause or negative weight, then fix the atom(s)
        long double cost = state_->getClauseCost(i);
        if (atoms->size() == 1 || cost < 0)
        {
          if (updebug) cout << "unit or negative" << endl;
            // We are fixing atoms, so we have possible propagation
          done = false;
            // Fix each atom in the clause
          for (int j = 0; j < atoms->size(); j++)
          {
            int lit = (*atoms)[j];
              // Set value to true if pos. lit. and pos. clause
              // or neg. lit. and neg. clause
            bool value = ((lit > 0 && cost > 0) || (lit < 0 && cost < 0))
              ? true : false;
            state_->fixAtom(abs(lit), value);
          }
        }
        else
        {
          if (updebug) cout << endl;
        }
        delete atoms;
      }
    }
      // Done with UP, save this state
    state_->saveLowState();
    if (updebug) cout << "Leaving UnitPropagation::infer" << endl;
  }
  
  /**
   * Prints the best state found.
   */
  void printProbabilities(ostream& out)
  {
    for (int i = 0; i < state_->getNumAtoms(); i++)
    {
      state_->printGndPred(i, out);
      out << " " << state_->getValueOfLowAtom(i + 1) << endl;
    }
  }

  /**
   * Puts the predicates whose truth value has changed with respect to the
   * reference vector oldProbs in string form and the corresponding
   * probabilities of each predicate (1 or 0) in two vectors.
   * 
   * @param changedPreds Predicates whose truth values have changed are put
   * here.
   * @param probs The probabilities corresponding to the predicates in
   * nonZeroPreds are put here (the number 1 or 0).
   * @param oldProbs Reference truth values for checking for changes.
   * @param probDelta This parameter is ignored for MAP inference (either the
   * truth value has changed or it hasn't).
   */
  void getChangedPreds(vector<string>& changedPreds, vector<float>& probs,
                       vector<float>& oldProbs, const float& probDelta)
  {
    changedPreds.clear();
    probs.clear();
    int numAtoms = state_->getNumAtoms();
      // Atoms may have been added to the state, previous tv was 0
    oldProbs.resize(numAtoms, 0);
    for (int i = 0; i < numAtoms; i++)
    {
      int tv = state_->getValueOfLowAtom(i + 1);
      if (tv != oldProbs[i])
      {
          // Truth value has changed: Store new value in oldProbs and add to
          // two return vectors
        oldProbs[i] = tv;
        ostringstream oss(ostringstream::out);
        state_->printGndPred(i, oss);
        changedPreds.push_back(oss.str());
        probs.push_back(tv);
      }
    }
  }

  /**
   * Gets the truth value of a ground predicate in the best state found.
   */
  double getProbability(GroundPredicate* const& gndPred)
  {
    int idx = state_->getGndPredIndex(gndPred);
    int truthValue = 0;
    if (idx >= 0) truthValue = state_->getValueOfLowAtom(idx + 1);
    return truthValue;
  }

  /**
   * Prints the predicates set to true in the best state to a stream.
   */
  void printTruePreds(ostream& out)
  {
    for (int i = 0; i < state_->getNumAtoms(); i++)
    {
      if (state_->getValueOfLowAtom(i + 1))
      {
        state_->printGndPred(i, out);
        out << endl;
      }
    }
  }
  
 private:
 
   
};


#endif /*UNITPROPAGATION_H_*/
