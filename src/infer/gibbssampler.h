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
#ifndef GIBBSSAMPLER_H_
#define GIBBSSAMPLER_H_

#include "mcmc.h"
#include "gibbsparams.h"
#include "maxwalksat.h"
#include "convergencetest.h"
#include "gelmanconvergencetest.h"

  // Initialize values randomly or with MaxWalkSat?
enum WalksatType { NONE = 0, MAXWALKSAT = 1 };
  // Set to true for more output
const bool gibbsdebug = false;

/**
 * Gibbs sampling algorithm.
 */
class GibbsSampler : public MCMC
{
 public:

  /**
   * Constructor. User-set parameters are set and an instance of MaxWalksat is
   * created for initial values (not used if initialized randomly).
   * 
   * @see MCMC#Constructor(VariableState*, long int, const bool&, MCMCParams*)
   */
  GibbsSampler(VariableState* state, long int seed,
               const bool& trackClauseTrueCnts, GibbsParams* gibbsParams)
    : MCMC(state, seed, trackClauseTrueCnts, gibbsParams)
  {
      // User-set parameters
    gamma_ = gibbsParams->gamma;
    epsilonError_ = gibbsParams->epsilonError;
    fracConverged_ = gibbsParams->fracConverged;
    walksatType_ = gibbsParams->walksatType;
    testConvergence_ = gibbsParams->testConvergence;
    samplesPerTest_ = gibbsParams->samplesPerTest;
    
      // We don't need to track clause true counts in MWS
    mws_ = new MaxWalkSat(state_, seed, false, gibbsParams->mwsParams);
  }

  /**
   * Destructor. Convergence tests and the instance of MaxWalksat are deleted.
   */
  ~GibbsSampler()
  {
    deleteConvergenceTests(burnConvergenceTests_, gibbsConvergenceTests_,
                           state_->getNumAtoms());
    delete mws_;
  }
  
  /**
   * Initializes the Gibbs sampler.
   */
  void init()
  {
      // Initialize gndPreds' truthValues & wts
    initTruthValuesAndWts(numChains_);
     printNetwork(cout);
    cout << "Initializing Gibbs sampling " ;
      // Initialize with MWS
    if (walksatType_ == 1)
    {
      cout << "with MaxWalksat" << endl;
      for (int c = 0; c < numChains_; c++)
      {
        cout << "for chain " << c << "..." << endl;
        mws_->init();
        mws_->infer();
        saveLowStateToChain(c);
      }
      if (numChains_ == 1) state_->saveLowStateToGndPreds();
    }
      // Initialize randomly
    else
    {
      cout << "randomly" << endl;
      randomInitGndPredsTruthValues(numChains_);
    }
      
      // Initialize gndClauses' number of satisfied literals
    //int start = 0;
    initNumTrueLits(numChains_);

    int numGndPreds = state_->getNumAtoms();
      // Initialize convergence test
    initConvergenceTests(burnConvergenceTests_, gibbsConvergenceTests_,
                         gamma_, epsilonError_, numGndPreds, numChains_);
  }

  /**
   * Runs Gibbs sampling.
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

    cout << "Running Gibbs sampling..." << endl;
      // Sampling loop
    int sample = 0;
    int numSamplesPerPred = 0;
    bool done = false;
    while (!done)
    {
      ++sample;

      if (sample % samplesPerTest_ == 0)
      { 
        currentTimeSec = timer.time();
        secondsElapsed = currentTimeSec-startTimeSec;
        cout << "Sample (per pred per chain) " << sample << ", time elapsed = ";
        Timer::printTime(cout, secondsElapsed); cout << endl;
      }

        // For each chain, for each node, generate the node's new truth value
      for (int c = 0; c < numChains_; c++) 
      {
        performGibbsStep(c, burningIn, affectedGndPreds,
                         affectedGndPredIndices);
        if (!burningIn) numSamplesPerPred++;
      }
  
        // Add current truth values to the convergence testers
      for (int i = 0; i < state_->getNumAtoms(); i++) 
      {
        const bool* vals = truthValues_[i].getItems();
          //WARNING: implicit cast from bool* to double*
        if (burningIn) burnConvergenceTests_[i]->appendNewValues(vals);
        else           gibbsConvergenceTests_[i]->appendNewValues(vals);
      }

      if (sample % samplesPerTest_ != 0) continue;      
      if (burningIn) 
      {
          // Use convergence criteria stated in "Probability and Statistics",
          // DeGroot and Schervish
        bool burnConverged = false;
        
        if (testConvergence_)
          burnConverged = 
            GelmanConvergenceTest::checkConvergenceOfAll(burnConvergenceTests_,
                                                         state_->getNumAtoms(),
                                                         true);
        if (   (sample >= burnMinSteps_ && burnConverged)
            || (burnMaxSteps_ >= 0 && sample >= burnMaxSteps_)
            || (maxSeconds_ > 0 && secondsElapsed >= maxSeconds_))
        {
          cout << "Done burning. " << sample << " samples per pred per chain";
          if (testConvergence_)
          {
            cout << " (" << (burnConverged? "converged":"didn't converge") 
                 <<" at total of " << numChains_*sample << " samples per pred)";
          }
          cout << endl;
          burningIn = false;
          sample = 0;          
        }
      }
      else
      {  // Doing actual gibbs sampling
        bool gibbsConverged = false;
        
        if (testConvergence_)
          gibbsConverged =
            ConvergenceTest::checkConvergenceOfAtLeast(gibbsConvergenceTests_, 
                                                       state_->getNumAtoms(),
                                                       sample, fracConverged_,
                                                       true);

        if (   (sample >= minSteps_ && gibbsConverged) 
            || (maxSteps_ >= 0 && sample >= maxSteps_)
            || (maxSeconds_ > 0 && secondsElapsed >= maxSeconds_)) 
        {
          cout << "Done Gibbs sampling. " << sample 
               << " samples per pred per chain";
          if (testConvergence_)
          {
            cout << " (" << (gibbsConverged? "converged":"didn't converge") 
                 <<" at total of " << numSamplesPerPred << " samples per pred)";
          }
          cout << endl;
          done = true;
        }
      }
      cout.flush();
    } // while (!done)
    
    cout<< "Time taken for Gibbs sampling = "; 
    Timer::printTime(cout, timer.time() - startTimeSec); cout << endl;

      // update gndPreds probability that it is true
    for (int i = 0; i < state_->getNumAtoms(); i++)
    {
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
   * Initializes convergence tests for burning in and sampling.
   */
   void initConvergenceTests(GelmanConvergenceTest**& burnConvergenceTests,
                            ConvergenceTest**& gibbsConvergenceTests, 
                            const double& gamma, const double& epsilonFrac, 
                            const int& numGndPreds, const int& numChains)
  {
    burnConvergenceTests = new GelmanConvergenceTest*[numGndPreds];
    gibbsConvergenceTests = new ConvergenceTest*[numGndPreds];
    for (int i = 0; i < numGndPreds; i++) 
    {
      burnConvergenceTests[i]  = new GelmanConvergenceTest(numChains);
      gibbsConvergenceTests[i] = new ConvergenceTest(numChains, gamma,
                                                     epsilonFrac);
    }
  }

  /**
   * Deletes convergence tests.
   */
  void deleteConvergenceTests(GelmanConvergenceTest**& burnConvergenceTests,
                              ConvergenceTest**& gibbsConvergenceTests, 
                              const int& numGndPreds)
  {
    for (int i = 0; i < numGndPreds; i++) 
    {
      delete burnConvergenceTests[i];
      delete gibbsConvergenceTests[i];
    }
    delete [] burnConvergenceTests;
    delete [] gibbsConvergenceTests;
  }  
  
 private:
    // Gamma used by convergence test
  double gamma_;
    // Epsilon used by convergence test
  double epsilonError_;
    // Fraction of samples needed to converge
  double fracConverged_;
    // 0 = Initialize randomly, 1 = initialize with MaxWalksat
  int walksatType_;
    // If true, test for convergence, otherwise do not test
  int testConvergence_;
    // Number of samples between checking for convergence
  int samplesPerTest_;
    // Convergence test for burning in
  GelmanConvergenceTest** burnConvergenceTests_;
    // Convergence test for sampling
  ConvergenceTest** gibbsConvergenceTests_;
    
    // MaxWalksat is used for initialization if walksatType_ = 1
  MaxWalkSat* mws_;  
};

#endif /*GIBBSSAMPLER_H_*/
