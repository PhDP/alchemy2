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
#ifndef SIMULATEDTEMPERING_H_
#define SIMULATEDTEMPERING_H_

#include "mcmc.h"
#include "simulatedtemperingparams.h"
#include "maxwalksat.h"
#include "convergencetest.h"
#include "gelmanconvergencetest.h"

/**
 * Simulated Tempering algorithm.
 */
class SimulatedTempering : public MCMC
{
 public:

  /**
   * Constructor.
   */  
  SimulatedTempering(VariableState* state, long int seed,
                     const bool& trackClauseTrueCnts, 
                     SimulatedTemperingParams* stParams)
    : MCMC(state, seed, trackClauseTrueCnts, stParams)
  {
      // User-set parameters
    subInterval_ = stParams->subInterval;
    numST_ = stParams->numST;
    numSwap_ = stParams->numSwap;
      // Number of chains is determined here
    numChains_ = numSwap_*numST_;        
    // ------------------------------------------ //
    // Chained method
    //  10 chains: i and i+1 swap attempt at
    //      selInterval*k + selInterval/10*i
    // ------------------------------------------ //
      // 9 possible swaps out of 10 chains
    selInterval_ = subInterval_*(numSwap_ - 1);

      // invTemp for chain chainIds_[i]
    invTemps_ = new double*[numST_];
      // curr chainId for ith temperature
    chainIds_ = new int*[numST_];
      // curr tempId for ith chain
    tempIds_ = new int*[numST_];
      // We don't need to track clause true counts in mws
    mws_ = new MaxWalkSat(state_, seed, false, stParams->mwsParams);
  }

  /**
   * Destroys allocated memory.
   */
  ~SimulatedTempering()
  {
    for (int i = 0; i < numST_; i++)
    {
      delete [] invTemps_[i];
      delete [] chainIds_[i];
      delete [] tempIds_[i];
    }
    delete [] invTemps_;
    delete [] chainIds_;
    delete [] tempIds_;
    delete mws_;
  }
  
  /**
   * Initializes Simulated Tempering.
   */
  void init()
  {
      // Initialize gndPreds' truthValues & wts
    //state_->initTruthValuesAndWts(numChains_, start);
    initTruthValuesAndWts(numChains_);

      // Initialize with MWS
    cout << "Initializing Simulated Tempering with MaxWalksat" << endl;
    state_->eliminateSoftClauses();
      // Set num. of solutions temporarily to 1
    int numSolutions = mws_->getNumSolutions();
    mws_->setNumSolutions(1);
    for (int c = 0; c < numChains_; c++)
    {
      cout << "for chain " << c << "..." << endl;
        // Initialize with MWS
      mws_->init();
      mws_->infer();
      saveLowStateToChain(c);
    }
    mws_->setNumSolutions(numSolutions);
    state_->resetDeadClauses();

    // *** Initialize temperature schedule ***
    double maxWt = state_->getMaxClauseWeight();
    double maxWtForEvenSchedule = 100.0;
    double base = log(maxWt) / log((double)numSwap_);
    double* divs = new double[numSwap_];
    divs[0] = 1.0;

    for (int i = 1; i < numSwap_; i++)
    {
      divs[i] = divs[i - 1] / base;
    }

    for (int i = 0; i < numST_; i++)
    {
      invTemps_[i] = new double[numSwap_];
      chainIds_[i] = new int[numSwap_];
      tempIds_[i]  = new int[numSwap_];
      for (int j = 0; j < numSwap_; j++)
      {         
        chainIds_[i][j] = j;
        tempIds_[i][j] = j;
          // log vs even
        if (maxWt > maxWtForEvenSchedule)
        {
          invTemps_[i][j] = divs[j];
        }
        else
        {
          invTemps_[i][j] = 1.0-((double)j)/((double) numSwap_);
        }
      }
    }
    delete [] divs;
      
      // Initialize gndClauses' number of satisfied literals
    //int start = 0;
    initNumTrueLits(numChains_);
  }

  /**
   * Runs Simulated Tempering.
   */
  void infer()
  {
    initNumTrue();
    Timer timer;
      // Burn-in only if burnMaxSteps positive
    bool burningIn = (burnMaxSteps_ > 0) ? true : false;
    double secondsElapsed = 0;
    double startTimeSec = timer.time();
    double currentTimeSec;
    int samplesPerOutput = 100;

      // If keeping track of true clause groundings, then init to zero
    if (trackClauseTrueCnts_)
      for (int clauseno = 0; clauseno < clauseTrueCnts_->size(); clauseno++)
        (*clauseTrueCnts_)[clauseno] = 0;

      // Holds the ground preds which have currently been affected
    GroundPredicateHashArray affectedGndPreds;
    Array<int> affectedGndPredIndices;

    int numAtoms = state_->getNumAtoms();
    for (int i = 0; i < numAtoms; i++)
    {
      affectedGndPreds.append(state_->getGndPred(i), numAtoms);
      affectedGndPredIndices.append(i);
    }
    for (int c = 0; c < numChains_; c++)
      updateWtsForGndPreds(affectedGndPreds, affectedGndPredIndices, c);
    affectedGndPreds.clear();
    affectedGndPredIndices.clear();

    cout << "Running Simulated Tempering sampling..." << endl;
      // Sampling loop
    int sample = 0;
    int numSamplesPerPred = 0;
    bool done = false;
    while (!done)
    {
      ++sample;

      if (sample % samplesPerOutput == 0)
      { 
        currentTimeSec = timer.time();
        secondsElapsed = currentTimeSec-startTimeSec;
        cout << "Sample (per pred per chain) " << sample << ", time elapsed = ";
        Timer::printTime(cout, secondsElapsed); cout << endl;
      }

        // Attempt to swap temperature
      if ((sample % selInterval_) % subInterval_ == 0)
      {
        int attemptTempId = (sample % selInterval_) / subInterval_;
        if (attemptTempId < numSwap_ - 1)
        {
          double wl, wh, itl, ith;
          for (int i = 0; i < numST_; i++)
          {
            int lChainId = chainIds_[i][attemptTempId];
            int hChainId = chainIds_[i][attemptTempId + 1];
              // compute w_low, w_high: e = -w
              // swap acceptance ratio=e^(0.1*(w_h-w_l))
            wl = getWeightSum(i*numSwap_ + lChainId);
            wh = getWeightSum(i*numSwap_ + hChainId);
            itl = invTemps_[i][attemptTempId];
            ith = invTemps_[i][attemptTempId + 1];

            if (wl <= wh || random() <= RAND_MAX*exp((itl - ith)*(wh - wl)))
            {
              chainIds_[i][attemptTempId] = hChainId;
              chainIds_[i][attemptTempId+1] = lChainId;
              tempIds_[i][hChainId] = attemptTempId;
              tempIds_[i][lChainId] = attemptTempId + 1;
            }
          }
        }
      }

        // Generate new truth value based on temperature
      for (int c = 0; c < numChains_; c++) 
      {
          // For each block: select one to set to true
        for (int i = 0; i < state_->getDomain()->getNumPredBlocks(); i++)
        {
            // If evidence atom exists, then all others stay false
          if (state_->getDomain()->getBlockEvidence(i)) continue;
 
          double invTemp =
            invTemps_[c/numSwap_][tempIds_[c/numSwap_][c%numSwap_]];
            // chosen is index in the block, block[chosen] is index in gndPreds_
          int chosen = gibbsSampleFromBlock(c, i, invTemp);

          const Predicate* pred =
            state_->getDomain()->getPredInBlock(chosen, i);
          GroundPredicate* gndPred = new GroundPredicate((Predicate*)pred);
          int idx = state_->getIndexOfGroundPredicate(gndPred);

          delete gndPred;
          delete pred;
      
            // If gnd pred in state:
          if (idx >= 0)
          {
            bool truthValue = truthValues_[idx][c];
              // If chosen pred was false, then need to set previous true
              // one to false and update wts
            if (!truthValue)
            {
              int blockSize = state_->getDomain()->getBlockSize(i);
              for (int j = 0; j < blockSize; j++)
              {
                const Predicate* otherPred = 
                  state_->getDomain()->getPredInBlock(j, i);
                GroundPredicate* otherGndPred =
                  new GroundPredicate((Predicate*)otherPred);
                int otherIdx = state_->getIndexOfGroundPredicate(gndPred);

                delete otherGndPred;
                delete otherPred;
      
                  // If gnd pred in state:
                if (otherIdx >= 0)
                {
                  bool otherTruthValue = truthValues_[otherIdx][c];
                  if (otherTruthValue)
                  {
                    truthValues_[otherIdx][c] = false;
              
                    affectedGndPreds.clear();
                    affectedGndPredIndices.clear();
                    gndPredFlippedUpdates(otherIdx, c, affectedGndPreds,
                                          affectedGndPredIndices);
                    updateWtsForGndPreds(affectedGndPreds,
                                         affectedGndPredIndices, c);
                  }
                }
              }
                // Set truth value and update wts for chosen atom
              truthValues_[idx][c] = true;
              affectedGndPreds.clear();
              affectedGndPredIndices.clear();
              gndPredFlippedUpdates(idx, c, affectedGndPreds,
                                    affectedGndPredIndices);
              updateWtsForGndPreds(affectedGndPreds, affectedGndPredIndices, c);
            }

              // If in actual sampling phase, track the num of times
              // the ground predicate is set to true
            if (!burningIn && tempIds_[c/numSwap_][c%numSwap_] == 0)
              numTrue_[idx]++;
          }
        }

          // Now go through all preds not in blocks
        for (int i = 0; i < state_->getNumAtoms(); i++) 
        {
            // Predicates in blocks have been handled above
          if (state_->getBlockIndex(i) >= 0) continue;
            // Calculate prob
          double invTemp =
            invTemps_[c/numSwap_][tempIds_[c/numSwap_][c%numSwap_]];
          double p = getProbabilityOfPred(i, c, invTemp);

            // Flip updates
          bool newAssignment = genTruthValueForProb(p);
          //if (newAssignment != pred->getTruthValue(c))
          if (newAssignment != truthValues_[i][c])
          {
            //pred->setTruthValue(c, newAssignment);
            truthValues_[i][c] = newAssignment;
            affectedGndPreds.clear();
            affectedGndPredIndices.clear();
            gndPredFlippedUpdates(i, c, affectedGndPreds,
                                  affectedGndPredIndices);
            updateWtsForGndPreds(affectedGndPreds, affectedGndPredIndices, c);
          }

            // if in actual sim. tempering phase, track the num of times
            // the ground predicate is set to true
          if (!burningIn && newAssignment &&
              tempIds_[c/numSwap_][c%numSwap_] == 0)
            //pred->incrementNumTrue();
            numTrue_[i]++;
        }
      }
      if (!burningIn) numSamplesPerPred += numST_;

        // If keeping track of true clause groundings
      if (!burningIn && trackClauseTrueCnts_)
        state_->getNumClauseGndings(clauseTrueCnts_, true);

      if (burningIn) 
      {
        if (   (burnMaxSteps_ >= 0 && sample >= burnMaxSteps_)
            || (maxSeconds_ > 0 && secondsElapsed >= maxSeconds_))
        {
          cout << "Done burning. " << sample << " samples per chain " << endl;
          burningIn = false;
          sample = 0;
        }
      }
      else 
      {
        if (   (maxSteps_ >= 0 && sample >= maxSteps_)
            || (maxSeconds_ > 0 && secondsElapsed >= maxSeconds_)) 
        {
          cout << "Done simulated tempering sampling. " << sample
               << " samples per chain" << endl;
          done = true;
        }
      }
      cout.flush();
    } // while (!done)
    
    cout<< "Time taken for Simulated Tempering sampling = "; 
    Timer::printTime(cout, timer.time() - startTimeSec); cout << endl;

      // update gndPreds probability that it is true
    for (int i = 0; i < state_->getNumAtoms(); i++)
    {
      //GroundPredicate* gndPred = state_->getGndPred(i);
      //gndPred->setProbTrue(gndPred->getNumTrue() / numSamplesPerPred);
      setProbTrue(i, numTrue_[i] / numSamplesPerPred);
    }
    
      // If keeping track of true clause groundings
    if (trackClauseTrueCnts_)
    {
        // Set the true counts to the average over all samples
      for (int i = 0; i < clauseTrueCnts_->size(); i++)
        (*clauseTrueCnts_)[i] = (*clauseTrueCnts_)[i] / numSamplesPerPred;
    }
  }
  
 private:
 
  /**
   * Calculates energy as the sum of weights of satisfied pos. clauses and
   * unsatisfied neg. clauses.
   * 
   * @param chainIdx Index of chain for which weight sum is calculated.
   * @return Sum of weights of sat. pos. clauses and unsat. neg. clauses
   */
  long double getWeightSum(const int& chainIdx)
  {
    long double w = 0;
    for (int i = 0; i < state_->getNumClauses(); i++)
    {
      long double wt = state_->getClauseCost(i);
      if ((wt > 0 && numTrueLits_[i][chainIdx] > 0) ||
          (wt < 0 && numTrueLits_[i][chainIdx] == 0))
        w += abs(wt);
    }
    return w;
  }
 
 private:
 
    // User-set parameters:
    // Selection interval between swap attempts
  int subInterval_;
    // Number of simulated tempering runs
  int numST_;
    // Number of swapping chains
  int numSwap_;

    // MaxWalksat is used for initialization
  MaxWalkSat* mws_;  

    // 9 possible swaps out of 10 chains
  int selInterval_;
    // invTemp for chain chainIds_[i]
  double** invTemps_;
    // curr chainId for ith temperature
  int** chainIds_;
    // curr tempId for ith chain
  int** tempIds_; 
};

#endif /*SIMULATEDTEMPERING_H_*/
