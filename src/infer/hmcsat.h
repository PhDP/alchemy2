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
#ifndef HMCSAT_H_					   
#define HMCSAT_H_

#include "mcmc.h"
#include "mcsatparams.h"
#include "hmaxwalksat.h"
#include "hvariablestate.h"

const int hmsdebug = false;

/**
* HMC-SAT is an MCMC inference algorithm designed to deal efficiently with
* probabilistic and deterministic dependencies (See Poon and Domingos, 2006).
* It wraps a procedure around SampleSat, thus enabling it to sample nearly
* uniform.
*/
class HMCSAT : public MCMC
{
 public:
 /**
  * Constructor: Constructs unit propagation and SampleSat.
  */  
  HMCSAT(HVariableState* state, long int seed, const bool& trackClauseTrueCnts,
         MCSatParams* mcsatParams) 
    : MCMC(state, seed, trackClauseTrueCnts, mcsatParams)
  {
	Timer timer1;

      // We don't need to track clause true counts in up and ss
    mws_ = new HMaxWalkSat(hstate_, seed, false, mcsatParams->mwsParams);
    bMaxWalkSat_ = false;
    contSamples_.growToSize(hstate_->contAtomNum_);
    print_vars_per_sample_ = false;
    
    if (hmsdebug)
    {
      cout << "[MCSAT] ";
      Timer::printTime(cout, timer1.time());
      cout << endl;
      timer1.reset();
    }
  }

 /**
  * Destructor: Instances of unit propagation and MaxWalksat are deleted.
  */
  ~HMCSAT()
  {
    delete mws_;
  }

  void initNumTrueTotal()
  {
    int numDisPreds = hstate_->getNumAtoms();
    numTrue_.growToSize(numDisPreds, 0);
  }

 /**
  * Initializes MC-SAT with MaxWalksat on hard clauses.
  */
  void init()
  {
    Timer timer1;
    assert(numChains_ == 1);

    cout << "Initializing MC-SAT with MaxWalksat on hard clauses..." << endl;
    
    initNumTrueTotal();
    hstate_->eliminateSoftClauses();
    hstate_->setInferenceMode(hstate_->MODE_HARD);
      // Set num. of solutions temporarily to 1
    int numSolutions = mws_->getNumSolutions();
    mws_->setNumSolutions(1);
      // Initialize with MWS
/*    
    if (!bMaxWalkSat_)
    {
      hstate_->makeUnitCosts();
    }

      // get the initialization by satisfying only hard clauses
      // Since the hybrid constraints' threshold are intialized to very large negative values, they are all satisfied and wound't be involved in inference here.
    mws_->setHeuristic(TABU);
*/
    mws_->init();
    mws_->infer();

    if (hmsdebug) 
    {
      cout << "Low state:" << endl;
      hstate_->printLowState(cout);
    }

    hstate_->saveLowStateToGndPreds();
    hstate_->saveCurrentAsLastAssignment();
    hstate_->UpdateHybridConstraintTh();
      // Set heuristic to SampleSat (Initialization could have been different)
    mws_->setHeuristic(SS);
    mws_->setNumSolutions(numSolutions);
    mws_->setTargetCost(0.0);
    hstate_->resetDeadClauses();
    sample_ = 0;

	// state_->makeUnitCosts();
    hstate_->setInferenceMode(hstate_->MODE_SAMPLESAT);	

    if (hmsdebug)
    {
      cout << "[MCSAT.init] ";
      Timer::printTime(cout, timer1.time());
      cout << endl;
      timer1.reset();
    }
  }

	void printContSamples()
	{
		for(int i = 0; i < contSamples_.size(); i ++)
		{
			for (int j = 0; j < contSamples_[i].size(); j++)
			{
				contSampleLog_ << contSamples_[i][j] << " ";
			}
			contSampleLog_ << endl;
		}
	}

	/**
	* Performs MC-SAT inference.
	*/
	void infer()
	{
  	Timer timer1;
		// Burn-in only if burnMaxSteps positive
		bool burningIn = (burnMaxSteps_ > 0) ? true : false;
		double secondsElapsed = 0;
    //upSecondsElapsed_ = 0;
		ssSecondsElapsed_ = 0;
		double startTimeSec = timer1.time();
		double currentTimeSec;
		int samplesPerOutput = 100;

		// Holds the ground preds which have currently been affected
		GroundPredicateHashArray affectedGndPreds;
		Array<int> affectedGndPredIndices;
		
		// Update the weights for Gibbs step
		int numAtoms = hstate_->getNumAtoms();
		for (int i = 0; i < numAtoms; i++)
		{
			affectedGndPreds.append(hstate_->getGndPred(i), numAtoms);
			affectedGndPredIndices.append(i);
		}
		updateWtsForGndPredsH(affectedGndPreds, affectedGndPredIndices, 0);
		affectedGndPreds.clear();
		affectedGndPredIndices.clear();

    if (hmsdebug)
    {	
      cout << "[MCSAT.infer.prep] "; 
      Timer::printTime(cout, timer1.time());
      cout << endl;
      timer1.reset();
    }

		cout << "Running MC-SAT sampling..." << endl;
		// Sampling loop
		int numSamplesPerPred = 0;
		bool done = false;
		while (!done)
		{
			++sample_;
			if (sample_ % samplesPerOutput == 0)
			{ 
				currentTimeSec = timer1.time();
				secondsElapsed = currentTimeSec - startTimeSec;
        cout << "Sample (per pred) " << sample_ << ", time elapsed = ";
        Timer::printTime(cout, secondsElapsed);
        cout << ", num. preds = " << hstate_->getNumAtoms();
		cout << ", num. clauses = " << hstate_->getNumClauses();
		cout << endl;		
			}

          performMCSatStep(burningIn);

			if (!burningIn) numSamplesPerPred++;

			if (burningIn) 
			{
				if ((burnMaxSteps_ >= 0 && sample_ >= burnMaxSteps_) || (maxSeconds_ > 0 && secondsElapsed >= maxSeconds_))
				{
					cout << "Done burning. " << sample_ << " samples." << endl;
					burningIn = false;
					sample_ = 0;
				}
			}
			else 
			{
				if (   (maxSteps_ >= 0 && sample_ >= maxSteps_)
					|| (maxSeconds_ > 0 && secondsElapsed >= maxSeconds_)) 
				{
					cout << "Done MC-SAT sampling. " << sample_ << " samples."             
						<< endl;
					done = true;
				}
			}
			cout.flush();
		} // while (!done)

		cout<< "Time taken for MC-SAT sampling = "; 
		Timer::printTime(cout, timer1.time() - startTimeSec); cout << endl;

    //cout<< "Time taken for unit propagation = "; 
    //Timer::printTime(cout, upSecondsElapsed_); cout << endl;

		cout<< "Time taken for SampleSat = "; 
		Timer::printTime(cout, ssSecondsElapsed_); cout << endl;

		// Update gndPreds probability that it is true
		for (int i = 0; i < hstate_->getNumAtoms(); i++)
		{
			setProbTrue(i, numTrue_[i] / numSamplesPerPred);
		}
/*
		// If keeping track of true clause groundings
		if (trackClauseTrueCnts_)
		{
			// Set the true counts to the average over all samples
			for (int i = 0; i < clauseTrueCnts_->size(); i++)
				(*clauseTrueCnts_)[i] = (*clauseTrueCnts_)[i] / numSamplesPerPred;
			for (int i = 0; i < clauseTrueCntsCont_->size(); i++)
				(*clauseTrueCntsCont_)[i] = (*clauseTrueCntsCont_)[i] / numSamplesPerPred;
		}
*/		
	}


	void SetContSampleFile(const char* contsamplelog)
	{
		contSampleLog_.open(contsamplelog);
		if (!contSampleLog_.is_open())
		{
			//cout << "cont sample log file is not specified." << endl;
		}
	}



	void SetPrintVarsPerSample(bool pps) {
            print_vars_per_sample_ = pps; 
        }


	bool bMaxWalkSat_;  // Variable indicating whether this inference is for hybrid MaxWalkSAT or hybrid MCSAT.
	
private:
	/**
	* Performs one step of MC-SAT.
	* 
	* @param burningIn Indicates if we are in a burning-in phase or not.
	*/
	void performMCSatStep(const bool& burningIn)
	{
		if (hmsdebug)
		{
			cout << "Num of clauses " << hstate_->getNumClauses() << endl;
            cout << "Num of false clauses " << hstate_->getNumFalseClauses() << endl;
			cout << "Num of dead clauses " << hstate_->getNumDeadClauses() << endl;
		}
		Timer timer;
		double startTime;
		if (hmsdebug) cout << "Entering MC-SAT step" << endl;

      // Clause selection
    hstate_->setUseThreshold(true);
    hstate_->updatePrevSatisfied();

    int start = 0;
    
    hstate_->resetMakeBreakCostWatch();

    hstate_->killClauses(start);

		// SampleSat on the clauses
		startTime = timer.time();
		mws_->init();
		mws_->infer();
		ssSecondsElapsed_ += (timer.time() - startTime);
		//here we save the lowest state as the current state
		hstate_->saveLowStateToGndPreds();
		if (print_vars_per_sample_) {
		    cout << "HMCS Iter#" << sample_ << ":\t";
		    for (int i = 0; i < hstate_->getNumAtoms(); ++i) {
		        cout << hstate_->atom_[i + 1] << "\t";
		    }

		    for (int i = 0; i < hstate_->getNumContAtoms(); ++i) {
		        cout << hstate_->contAtoms_[i+1] << "\t";
		    }
		    cout << endl;
		}
		if (hmsdebug && hstate_->costOfTotalFalseConstraints_ < mws_->getTargetCost() + SMALLVALUE)
		{
			hstate_->printLowStateAll(cout);
		}
		if (hmsdebug && hstate_->costOfTotalFalseConstraints_ > mws_->getTargetCost() + SMALLVALUE)
		{
			cout << "not all satisfied, at sample " << sample_ << ". " << endl;
			cout << "Error at sample: " << sample_ << ", " << hstate_->costOfTotalFalseConstraints_ << endl;
			cout << "False cont constraints: " << hstate_->costHybridFalseConstraint_ << ", false discrete: " << hstate_->getCostOfFalseClauses() << endl;
			//print out false clauses & constraints
			hstate_->printFalseClauses(cout);
		}
		// Reset parameters needed for HMCSAT step		
		int numAtoms = hstate_->getNumAtoms();
		for (int i = 0; i < numAtoms; i++)
		{
			GroundPredicate* gndPred = hstate_->getGndPred(i);
			bool newAssignment = hstate_->getValueOfLowAtom(i + 1);
			// No need to update weight but still need to update truth/NumSat
			if (newAssignment != gndPred->getTruthValue()) {
				gndPred->setTruthValue(newAssignment);
				updateClauses(i);
			}
		}
		// Write the current sample of continuous variables to file.
		for(int i = 0; i < hstate_->contAtomNum_; i++)
		{
			contSamples_[i].append(hstate_->contAtoms_[i+1]);
			contSampleLog_ << hstate_->contAtoms_[i+1] << " ";
		}
		contSampleLog_ << endl;

		for( int i = 0; i < hstate_->getNumAtoms(); i++)
		{
			bool newAssignment = hstate_->getValueOfAtom(i+1);
			if (!burningIn && newAssignment) numTrue_[i]++;
		}	
		hstate_->resetFixedAtoms();
		hstate_->resetDeadClauses(); // update true lit num for each dis clause
		hstate_->saveCurrentAsLastAssignment(); // backup
		// Continuous thresholds are updated after each sample round.
		hstate_->UpdateHybridConstraintTh();
		hstate_->setUseThreshold(false);
		// If keeping track of true clause groundings
		if (!burningIn && trackClauseTrueCnts_)
		{
			hstate_->getNumClauseGndings(clauseTrueCnts_, true);
			hstate_->getContClauseGndings(clauseTrueCntsCont_);			
		}		
		if (hmsdebug) cout << "Leaving MC-SAT step" << endl;
	}

	/**
	* Updates the number of satisfied literals in the clauses in which a ground
	* predicate occurs.
	* 
	* @param gndPredIdx Index of ground predicate which occurs in the ground
	* clauses being updated.
	*/
	void updateClauses(const int& gndPredIdx)
	{
		if (hmsdebug) cout << "Entering updateClauses" << endl;
		GroundPredicate* gndPred = hstate_->getGndPred(gndPredIdx);
		Array<int>& negGndClauses =
			hstate_->getNegOccurenceArray(gndPredIdx + 1);
		Array<int>& posGndClauses =
			hstate_->getPosOccurenceArray(gndPredIdx + 1);
		int gndClauseIdx;
		bool sense;

		for (int i = 0; i < negGndClauses.size() + posGndClauses.size(); i++)
		{
			if (i < negGndClauses.size())
			{
				gndClauseIdx = negGndClauses[i];
				sense = false;
			}
			else
			{
				gndClauseIdx = posGndClauses[i - negGndClauses.size()];
				sense = true;
			}

			if (gndPred->getTruthValue() == sense)
				hstate_->incrementNumTrueLits(gndClauseIdx);
			else
				hstate_->decrementNumTrueLits(gndClauseIdx);
		}
		if (hmsdebug) cout << "Leaving updateClauses" << endl;
	}

private:
	ofstream contSampleLog_;
	int sample_;
	bool print_vars_per_sample_;

	Array<Array<double> > contSamples_;

	// Unit propagation is performed in MC-SAT  
	//HUnitPropagation* up_;
	// The base algorithm is SampleSat (MaxWalkSat with SS parameters)
	HMaxWalkSat* mws_;

	// Time spent on UnitPropagation
  //double upSecondsElapsed_;
	// Time spent on SampleSat
	double ssSecondsElapsed_;
};

#endif /*HMCSAT_H_*/
