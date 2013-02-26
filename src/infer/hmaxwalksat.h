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
#ifndef HMAXWALKSAT_H
#define HMAXWALKSAT_H

#include <cfloat>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <cstring>
#include <sstream>
#include <fstream>
#include "array.h"
#include "timer.h"
#include "util.h"
#include "sat.h"
#include "logichelper.h"
#include "maxwalksatparams.h"


const bool hmwsdebug = false;
#define BIGSTEP 100
/**
* The HMaxWalkSat algorithm. This code is based on the HMaxWalkSat
* package of Kautz et al. and SampleSat by Wei et al.
* 
* Walksat is achieved by using unit weights (1 and -1) in the clauses
* in the state and setting the target cost to 0.
* 
* SampleSat is achieved by using unit weights (1 and -1) in the clauses
* in the state and setting the target cost to 0 and the heuristic to SS.
*/
class HMaxWalkSat : public SAT
{
public:
	/**
	* Constructor: user-set parameters are set.
	*/
	HMaxWalkSat(HVariableState* state, int seed, const bool& trackClauseTrueCnts,
		MaxWalksatParams* params)
		: SAT(state, seed, trackClauseTrueCnts)
	{
		int numAtoms = hstate_->getNumAtoms();
		saSteps_ = 0;
		saInterval_ = 100;
		//falseDisNum_ = 0;
		// User-set parameters
		maxSteps_ = params->maxSteps;
		maxTries_ = params->maxTries;
		targetCost_ = params->targetCost;
		hard_ = params->hard;
      // Set this in the state
    hstate_->setBreakHardClauses(hard_);
		numSolutions_ = params->numSolutions;
		heuristic_ = params->heuristic;									   

		tabuLength_ = params->tabuLength;
		lazyLowState_ = params->lazyLowState;
		// For SampleSat
		saRatio_ = params->ssParams->saRatio;
		saTemp_ = params->ssParams->saTemp;
		saTempInit_ = saTemp_;
		lateSa_ = params->ssParams->lateSa;

		// Info from SAT class
		if (heuristic_ == TABU || heuristic_ == SS || heuristic_ == HMWS)
			changed_.growToSize(numAtoms + 1);
		numFlips_ = 0;

		// Info from HMaxWalkSat
		numerator_ = 50; // User-set?
		denominator_ = 100; // User-set?
		// We assume we are not in a block
		inBlock_ = false;
		flipInBlock_ = NOVALUE;
		// Initialize array of function pointers
		int idx = 0;
		pickcode[idx++] = &HMaxWalkSat::pickRandom;
		pickcode[idx++] = &HMaxWalkSat::pickBest;
		pickcode[idx++] = &HMaxWalkSat::pickHMCS;
		pickcode[idx++] = &HMaxWalkSat::pickSS;
		pickcode[idx++] = &HMaxWalkSat::pickHMWS;
		pickcode[idx++] = &HMaxWalkSat::pickSA;
		maxSeconds_ = -1;
	}

	/**
	* Destructor (destructors from superclasses are called)
	*/
	~HMaxWalkSat() { }


	void SetMaxSeconds(double maxSeconds)
	{
		maxSeconds_ = maxSeconds;
	}

	void SetNoisePra(int numerator, int denominator)
	{
		numerator_ = numerator;
		denominator_ = denominator;
	}

	void SetSAInterval(int saInterval)
	{
		saInterval_ = saInterval;
	}
	/**
	* Initializes the HMaxWalkSat algorithm. A random assignment is given to the
	* atoms and the state is initialized.
	*/
	void init()
	{
		if (heuristic_ == TABU || heuristic_ == SS || heuristic_ == HMWS)
		{
			int numAtoms = hstate_->getNumAtoms();
			if (changed_.size() != numAtoms + 1)
			{
				if (int(changed_.size()) < numAtoms + 1) {
					changed_.growToSize(numAtoms + 1);
				} else {
					changed_.shrinkToSize(numAtoms + 1);
				}
			}
			for(int i = 1; i <= numAtoms; ++i)
				changed_[i] = -(tabuLength_ + 1);
		}

		// Random initialize the variables 
		hstate_->initRandom();		
		//reinitialize the SATemp here
		saTemp_ = saTempInit_;
	}

	double getTargetCost()
	{
		return targetCost_;
	}

	void setSATempDownRatio(double saTempDownRatio)
	{
		saTempDownRatio_ = saTempDownRatio;
	}


	/**
	* Performs the given number of HMaxWalkSat inference tries.
	*/
	void infer()
	{
		int numtry = 0;
		int numsuccesstry = 0;
		Timer timer;
		double startTimeSec = timer.time();
		
		double currentTimeSec, secondsElapsed;
		// Whether the expected time is consumed by the inference.
		bool timeUp = false;
		double lowCost;
		if (!hstate_->bMaxOnly_)
		{
			lowCost = BigValue;
		} else {
			lowCost = -BigValue;
		}

		int totaltrynum = 0;
		// If keeping track of true clause groundings, then init to zero
		if (trackClauseTrueCnts_)
		{
			for (int clauseno = 0; clauseno < clauseTrueCnts_->size(); clauseno++)
				(*clauseTrueCnts_)[clauseno] = 0;
			for (int clauseno = 0; clauseno < clauseTrueCntsCont_->size(); clauseno++)
				(*clauseTrueCntsCont_)[clauseno] = 0;
		}

		int numhardfalse = 0;
		if (hstate_->bMaxOnly_)
		{
			for(int i = 0; i < hstate_->getNumFalseClauses(); i++)
			{
				GroundClause* clause = (*(hstate_->gndClauses_))[hstate_->falseClause_[i]];
				if (clause->isHardClause())
				{
					numhardfalse++;
				}
			}
			cout << "before HMWS, false dis clause#: " << hstate_->getNumFalseClauses() << " false hard clause#: " << numhardfalse<< endl;
		}

		// Perform up to maxTries tries or numSolutions successful tries
		while (numsuccesstry < numSolutions_ &&	numtry < maxTries_ && !timeUp) //default 1*10
		{
			numtry++;
			numFlips_ = 0;
			bool lowStateInThisTry = false;

			// Make random assignment in subsequent tries
			// (for first try, init() should have been called)
			if (numtry > 1 && numsuccesstry == 0)
			{
				init();
				startTimeSec = timer.time();
			}

			while (((numFlips_ < maxSteps_ && maxSeconds_ < 0) || (maxSeconds_ > 0)) 
				   && numsuccesstry < numSolutions_)
			{	 
				currentTimeSec = timer.time();
				secondsElapsed = currentTimeSec - startTimeSec;
				if (maxSeconds_ > 0 && secondsElapsed > maxSeconds_) {
					cout << "Time is up for current HMaxWalkSAT try. "	<< endl;
					timeUp = true;
					break;
				}
				numFlips_++;
				int atomIdx = (this->*(pickcode[heuristic_]))();

				if (hmwsdebug)
				{
					cout << "HMWS picked atom: " << atomIdx << "." << endl;
				}
				
				//we need to decide whether to record the continuous sample here.
				if (atomIdx < 0) // Continuous variables
				{
					atomIdx = -atomIdx;
					if (!hstate_->bMaxOnly_)  // HMCSAT
					{
						if (pickedContValue_ == UNSOLVABLE) {
							continue;
						}
						hstate_->UpdateHybridClauseInfoByContAtom(atomIdx, pickedContValue_);
					}  // For HMCS ends
					else  // HMWS
					{
						if (pickedContValue_ == UNSOLVABLE) {
							continue;
						}
						hstate_->UpdateHybridClauseWeightSumByContAtom(atomIdx, pickedContValue_);
					}  // For HMWS ends
				}  // Continuous variable ends
				else if (atomIdx > 0 && atomIdx != NOVALUE)  // Discrete variables
				{
					if (!hstate_->bMaxOnly_)  // HMCS
					{
						if (hmwsdebug)
						{
							cout << "flipping dis atom: " << atomIdx << " before flip:"; hstate_->printDisAtom(atomIdx, cout); cout << endl;
						}

						flipAtom(atomIdx);
						if (hmwsdebug)
						{
							cout << "after flip in dis: ";hstate_->printDisAtom(atomIdx, cout); cout<< endl;
						}

						bool atomValue = hstate_->getValueOfAtom(atomIdx);
						hstate_->UpdateHybridClauseInfoByDisAtom(atomIdx, atomValue);
						if (hmwsdebug)
						{
							cout << "after flip in cont: ";hstate_->printDisAtom(atomIdx, cout); cout<< endl;
						}
					}  // For HMCS ends
					else  // HMWS
					{	
						flipAtom(atomIdx);
						bool atomValue = hstate_->getValueOfAtom(atomIdx);
						hstate_->UpdateHybridClauseWeightSumByDisAtom(atomIdx, atomValue);
					}  // For HMWS ends
				}  // Discrete variable ends
				/*else {
					continue;
				}*/

				double discost = 0, hybridcost = 0;
				long double costOfFalseConstraints = 0;
				if (!hstate_->bMaxOnly_)  // For HMCS
				{
					costOfFalseConstraints = hstate_->getCostOfTotalFalseConstraints();
				}
				else  // For HMWS
				{
					costOfFalseConstraints =  hstate_->getCostOfAllClauses(discost, hybridcost);
				}

				//long double lowCost = hstate_->getLowCost();
				// If the cost of false clauses is lower than the low cost of all
				// tries, then save this as the low state.
				if(!hstate_->bMaxOnly_)  // HMCS
				{
					if (costOfFalseConstraints <= lowCost) // <= saves the last satisfying assignment
					{
						if (hmwsdebug)
						{
							cout << "Cost of false clauses: " << costOfFalseConstraints
								<< " less than lowest cost: " << lowCost << endl;
						}
						lowCost = costOfFalseConstraints;
						hstate_->saveLowStateAll();
					}

					// If successful try
					// Add SMALLVALUE to targetCost_ because of double precision error
					// This needs to be fixed
					if (costOfFalseConstraints <= targetCost_ + SMALLVALUE)
						numsuccesstry++; 
				}
				else
				{
					if (costOfFalseConstraints >= lowCost) // <= saves the last satisfying assignment
					{
						lowCost = costOfFalseConstraints;
						hstate_->saveLowStateAll();
					}
				}					
			}  // End inner while.

			if (!hstate_->bMaxOnly_)  // HMCS
			{
				totaltrynum += numFlips_;
				if (lowStateInThisTry)
					reconstructLowState();
				hstate_->saveLowStateToDB(); // only for discrete
				if (hmwsdebug)
				{
					cout<< "In the end of try " << numtry << ": " << endl;
					cout<< "Lowest num. of false constraints: " << hstate_->getLowBadAll() << endl;
					cout<< "Lowest cost of false constraints: " << hstate_->getLowCostAll() << endl;
					cout<< "Number of flips: " << numFlips_ << endl;
				}
			}				
		}  // End outter while
		
		if (!hstate_->bMaxOnly_)  // HMCS
		{
			if (hmwsdebug)
			{
				cout << "Total flips num: " << totaltrynum << ". Succeed:" << numsuccesstry << " times." << endl;
			}
			
			if (numsuccesstry == 0) //no satisfying assignment found, should reset last time's assignmetn here to low state 
			{
				//print out false clauses & constraints
				hstate_->saveLastAsCurrentAssignment();
				hstate_->updateCost();
				hstate_->saveLowStateAll();
			}

			hstate_->saveLowToCurrentAll();
			// If keeping track of true clause groundings
			if (trackClauseTrueCnts_)
			{
				hstate_->updateNumTrueLits();
				hstate_->getNumClauseGndings(clauseTrueCnts_, true);   
				hstate_->getContClauseGndings(clauseTrueCntsCont_);
			}
		}  // HMCS ends.
		else  // HMWS
		{
			cout << "finish hybrid MaxWalkSAT." << endl;
			cout << "Lowest state cost:" << lowCost << endl; //			
			numhardfalse = 0;
			for(int i = 0; i < hstate_->getNumFalseClauses(); i++)
			{
				GroundClause* clause = (*(hstate_->gndClauses_))[hstate_->falseClause_[i]];
				if (clause->isHardClause())
				{
					numhardfalse++;
				}
			}
			cout << "after HMWS, false dis clause#: " 
                 << hstate_->getNumFalseClauses() << " false hard clause#: "
                 << numhardfalse << endl;
			hstate_->saveLowToCurrentAll();
			if (trackClauseTrueCnts_)
			{
				hstate_->updateNumTrueLits();
				hstate_->getNumClauseGndings(clauseTrueCnts_, true);
				hstate_->getContClauseGndings(clauseTrueCntsCont_);
			}
		}
	}

	const int getHeuristic()
	{
		return heuristic_;
	}

	void setHeuristic(const int& heuristic)
	{
		heuristic_ = heuristic;
	}
	
	virtual const Array<double>* getClauseTrueCnts() 
	{ 
		//assert(clauseTrueCnts_->size());
		if (!trackClauseTrueCnts_)
		{
			hstate_->updateNumTrueLits();
			hstate_->getNumClauseGndings(clauseTrueCnts_, true);
		}		
		return clauseTrueCnts_; 
	}

	virtual const Array<double>* getClauseTrueCntsCont() 
	{ 
		if (!trackClauseTrueCnts_)
		{
			hstate_->getContClauseGndings(clauseTrueCntsCont_);
		}
		return clauseTrueCntsCont_; 
	}

protected:

	/**
	* Flips the truth value of an atom and marks this iteration as the
	* last time it was flipped.
	* 
	* @param toFlip Index of atom to flip.
	*/
	void flipAtom (int toFlip)
	{
		//if (hmwsdebug) cout << "Entering HMaxWalkSat::flipAtom" << endl;
		if (toFlip == NOVALUE)
			return;

		// Flip the atom in the state
		hstate_->flipAtom(toFlip, -1);
		// Mark this flip as last time atom was changed if tabu is used
		if (heuristic_ == TABU || heuristic_ == SS || heuristic_ == HMWS)
		{
			changed_[toFlip] = numFlips_;
			changed_.growToSize(hstate_->getNumAtoms() + 1, -(tabuLength_ + 1));
		}

		if (lazyLowState_)
		{
			// Mark this variable as flipped since last low state. If already
			// present, then it has been flipped back to its value in the low
			// state and we can remove it.
			if (varsFlippedSinceLowState_.find(toFlip) !=
				varsFlippedSinceLowState_.end())
			{
				varsFlippedSinceLowState_.erase(toFlip);
			}
			else
			{
				varsFlippedSinceLowState_.insert(toFlip);
			}
		}
      if (hmwsdebug) cout << "Leaving HMaxWalkSat::flipAtom" << endl;
	}

	/**
	* Pick a random atom in a random unsatisfied pos. clause or a random
	* true literal in a random satsified neg. clause to flip.
	* 
	* @return Index of atom picked.
	*/
	int pickRandom()
	{
		return 0;
	}

	/**
	* Pick the best atom (clauses made satisfied - clauses made unsatisfied)
	* in a random unsatisfied clause.
	* 
	* @return Index of atom picked.
	*/
	int pickBest()
	{
		return 0;
	}

	/**
	* Pick the atom idx according to hybrid WalkSAT strategy. 
	* 
	* @return Index of atom picked.
	*/
	int pickHMCS()
	{
		if (hmwsdebug)
		{
			cout << "HMCS" << endl;
		}

		assert(!inBlock_);
		// Clause to fix picked at random
		// need to be changed to deal with hybrid case
		bool noisyPick = (numerator_ > 0 && random() % denominator_ < numerator_); 
		// Pick a clause.
		// For H-MCSAT, the picked clause is either a false discrete clause, or a hybrid one not satisfying the inequality constraint.
		// For HMWS, the picked clause is either a false discrete clause or a hybrid one.
		int toFix = hstate_->getRandomFalseClauseIndexHMCS();  // HMCS version.
		if (toFix ==	 NOVALUE) // No false clause.
		{
			return NOVALUE;
		}
		else if (toFix > 0) //discrete clause, clause idx = toFix - 1;
		{
			toFix --;
			if (hmwsdebug)
			{
				cout << "Pick false dis clause: " << toFix << ", "; hstate_->printDisClause(toFix, cout); cout << "Cost: " << hstate_->clauseCost_[toFix] << endl;
			}
			long double improvement;
			int clauseSize = hstate_->getClauseSize(toFix);
			long double cost = hstate_->getClauseCost(toFix);
			// Holds the best atoms so far
			Array<int> best;
			// How many are tied for best
			register int numbest = 0;
			// Best value so far
			long double bestvalue = -LDBL_MAX;
			// With prob. do a noisy pick, prob = denominator / numerator.
			// if random move, exclude illegitimate ones and place others in a lottery
			if (noisyPick)
			{
				for (int i = 0; i < clauseSize; ++i)
				{       
					register int lit = hstate_->getAtomInClause(i, toFix);
					// Neg. clause: Only look at true literals
					if (cost < 0 && !hstate_->isTrueLiteral(lit)) continue;
					best.append(abs(lit));
					numbest++;
				}
			}  // Noisy pick ends.
			else  // Greedy: pick the best value
			{
				for (int i = 0; i < clauseSize; ++i)
				{
					register int lit = hstate_->getAtomInClause(i, toFix);
					// Neg. clause: Only look at true literals
					if (cost < 0 && !hstate_->isTrueLiteral(lit)) continue;
					// var is the index of the atom
					register int var = abs(lit);
					improvement = calculateImprovementDisHMCS(var);
					if (mwsdebug) {
						cout << "Improvement of var " << var << " is: " << improvement << endl;
					}

					// TABU: If pos. improvement, then always add it to best
					if (improvement > 0 && improvement >= bestvalue)
					{ // Tied or better than previous best
						if (improvement > bestvalue)
						{
							numbest = 0;
							best.clear();
						}
						bestvalue = improvement;
						best.append(var);
						numbest++;
					} else if (improvement > -LDBL_MAX && tabuLength_ < numFlips_ - changed_[var]) {
						if (improvement >= bestvalue)
						{ // Tied or better than previous best
							if (improvement > bestvalue)
							{
								numbest = 0;
								best.clear();
							}
							bestvalue = improvement;
							best.append(var);
							numbest++;
						}
					}
				}
			}

			if (numbest == 0) 
			{
				if (hmwsdebug) cout << "Leaving MaxWalkSat::pickHMCS (NOVALUE)" << endl;
				return NOVALUE;
			}
			int toFlip = NOVALUE;
			// Exactly one best atom
			if (numbest == 1) {
				toFlip = best[0];
			} else {
				// Choose one of the best at random
				toFlip = best[random()% best.size()];
			}

			if (hmwsdebug) 
			{
				cout << "Improvement:" << bestvalue << endl;
			}
			return toFlip;
		}  // Discrete clause
		else  // Hybrid clause, clause idx = - (toFix + 1); toFix < 0
		{
			toFix++;
			int contClauseIdx = -toFix;
			// If the weight of the hybrid constraint is 0, then we skip the current pick.
			if (hstate_->hybridWts_[contClauseIdx] == 0) return NOVALUE;

			if (hmwsdebug) {
				cout << "Pick false hybrid clause: " << contClauseIdx << ", constraint:"; hstate_->printHybridConstraint(contClauseIdx, cout);
			}
			int contVarNum = hstate_->hybridContClause_[contClauseIdx].size();
			int disVarNum = hstate_->hybridDisClause_[contClauseIdx].size();
			int totalVarNum = contVarNum + disVarNum;
			if (noisyPick)  // Random pick and then flip, p = numerator_/denominator_;
			{
				set<int> pickedContVars;
				int inIdx = random() % totalVarNum;
				if (inIdx < contVarNum)  // Continuous variable
				{
					int contAtomIdx = hstate_->hybridContClause_[contClauseIdx][inIdx];
					// Solve the constraint in regarding with the chosen continuous variable.
					// Then compute the improvement in term of #false clause reduction.
					pickedContValue_ = hstate_->SolveHybridConstraintToContVarSample(contClauseIdx, inIdx);
					if (pickedContValue_ == FLIPDIS) {
						// Generate a random solution to the discrete part of the hybrid clause, so that making it true.
						Array<bool> solution;
						GetClauseRandomSolution(hstate_->hybridDisClause_[contClauseIdx].size(), true, hstate_->hybridConjunctionDisjunction_[contClauseIdx], &solution);
						// Update each discrete variable.
						for (int i = 0; i < solution.size(); ++i) {
							int atomIdx = hstate_->hybridDisClause_[contClauseIdx][i];
							if (solution[i] != hstate_->getValueOfAtom(atomIdx)) {
								flipAtom(atomIdx);
								bool atomValue = hstate_->getValueOfAtom(atomIdx);
								hstate_->UpdateHybridClauseInfoByDisAtom(atomIdx, atomValue);
							}
						}
						// Then check if the continuous part's value satisfies the constraint or not.
						double cont_part_value = hstate_->HybridClauseContPartValue(contClauseIdx);
						if (cont_part_value < hstate_->hybridConstraints_[contClauseIdx].vThreshold_) 						{
							// If not, solve the continuous part for the given variable, set both solutions and update status.
							double cont_var_value = hstate_->SolveConstraintAndRandomSample(contClauseIdx, inIdx);
							if (cont_var_value == UNSOLVABLE)
								return NOVALUE;
							else {
								// Update by the continuous variable.
								hstate_->UpdateHybridClauseInfoByContAtom(contAtomIdx, cont_var_value);
							}
						}
	
						// Return the status that no action needs to be performed.
						return NOVALUE;
					}  // FLIP_DIS ends

					return -contAtomIdx;
				}
				else  // discrete variable
				{
					// Check if flipping this discrete variable could satisfy this constraint.
					// If so, return it, if not, do nothing.
					int atomIdx = abs(hstate_->hybridDisClause_[contClauseIdx][inIdx - contVarNum]);
					hstate_->atom_[atomIdx] = !hstate_->atom_[atomIdx];
					double hybrid_clause_val = hstate_->HybridClauseValue(contClauseIdx);
					hstate_->atom_[atomIdx] = !hstate_->atom_[atomIdx];
					if (hstate_->isSatisfied(hstate_->hybridConstraints_[contClauseIdx], hybrid_clause_val)) {
						return atomIdx;
					} else {
						return NOVALUE;
					}
				}
			}  // Noisy pick
			else //greedy try to flip every variable and find the one get biggest improvement
			{
				long double bestvalue = -LDBL_MAX;
				Array<int> best;
				int numbest = 0;
				map<int, double> pickedContValues;
				for (int i = 0; i < contVarNum; ++i) {
					int contAtomIdx = hstate_->hybridContClause_[contClauseIdx][i];
					// Solve the constraint in regarding with the chosen continuous variable.
					// Then compute the improvement in term of #false clause reduction.
					double cont_var_val = 0;
					double improvement =
							hstate_->GetImprovementBySolvingHybridConstraintToContVar(contClauseIdx, i, cont_var_val);
					// Tied or better than previous best
					if (improvement != CANNOTSATCURRENT && improvement >= bestvalue) {
						if (improvement > bestvalue) {
							numbest = 0;
							best.clear();
						}
						bestvalue = improvement;
						best.append(-contAtomIdx);
						numbest++;
						pickedContValues[-contAtomIdx] = cont_var_val;
					}
				}

				for (int i = 0; i < disVarNum; ++i)
                {
					int disAtomIdx = abs(hstate_->hybridDisClause_[contClauseIdx][i]);
					double improvement = calculateImprovementDisHMCS(disAtomIdx);

					if (improvement >= bestvalue)
					{ // Tied or better than previous best
						if (improvement > bestvalue)
						{
							numbest = 0;
							best.clear();
						}
						bestvalue = improvement;
						best.append(disAtomIdx);
						numbest++;
					}
				}

				int atomIdx = NOVALUE;
				// Exactly one best atom
				if (numbest == 1)
					atomIdx = best[0];
				else
					// Choose one of the best at random
					atomIdx = best[random()%numbest];

				if (atomIdx < 0)  // Picked variable is a continuous variable.
				{
					pickedContValue_ = pickedContValues[atomIdx];
				}
				return atomIdx;
			}  // Greedy
		}  // Hybrid clause.
		return NOVALUE;
	}

	void checkAr(const Array<int>& ar)
	{
		int v = ar[0];
		for(int i = 0;i < ar.size(); i++)
		{
			if(ar[i] != v)
			{
				cout << v << "\t" << ar[i] << ". Error!!!!!" <<endl;
				break;
			}
		}
	}
	
	struct FlipDisVarIdxCont
	{
		Array<int> disIdx_;
	};

	int pickHMWS()
	{
		if (hmwsdebug)
		{
			cout << "HMWSSTEP" << endl;
		}

		assert(!inBlock_);
		// Clause to fix picked at random    
		//need to be changed to deal with hybrid case
		bool noisyPick = (numerator_ > 0 && random()%denominator_ < numerator_); 

		// Pick a clause.
		// For H-MCSAT, the picked clause is either a false discrete clause, or a hybrid one not satisfying the inequality constraint.
		// For HMWS, the picked clause is either a false discrete clause or a hybrid one.
		int toFix = hstate_->getRandomFalseClauseIndexHMWS();  // HMCS version.

		if (toFix == NOVALUE) // no false clause
		{
			return NOVALUE;
		} else if (toFix > 0) {  // Discrete clause, clause idx = toFix - 1;
			toFix--;
			if (hmwsdebug)
			{
				cout << "Pick false dis clause: " << toFix << ", "; hstate_->printDisClause(toFix, cout);
			}
			long double improvement;
			int clauseSize = hstate_->getClauseSize(toFix);
			long double cost = hstate_->getClauseCost(toFix);
			// Holds the best atoms so far
			Array<int> best;
			// How many are tied for best
			register int numbest = 0;
			// Best value so far
			long double bestvalue = -LDBL_MAX;
			// With prob. do a noisy pick, prob = denominator / numerator.
			// if random move, exclude illegitimate ones and place others in a lottery
			if (noisyPick)
			{
				for (int i = 0; i < clauseSize; i++)
				{       
					register int lit = hstate_->getAtomInClause(i, toFix);
					// Neg. clause: Only look at true literals
					if (cost < 0 && !hstate_->isTrueLiteral(lit)) continue;
					best.append(abs(lit));
					numbest++;
				}
			}
			// greedy: pick the best value
			else {
				for (int i = 0; i < clauseSize; i++)
				{
					register int lit = hstate_->getAtomInClause(i, toFix);
					// Neg. clause: Only look at true literals
					if (cost < 0 && !hstate_->isTrueLiteral(lit)) continue;
					// var is the index of the atom
					register int var = abs(lit);
					improvement = calculateImprovementDisHMWS(var);
					if (mwsdebug) {
						cout << "Improvement of var " << var << " is: " << improvement << endl;
					}

					// TABU: If pos. improvement, then always add it to best
					if (improvement > 0 && improvement >= bestvalue)
					{ // Tied or better than previous best
						if (improvement > bestvalue)
						{
							numbest = 0;
							best.clear();
						}
						bestvalue = improvement;
						best.append(var);
						numbest++;
					} else if (improvement > -LDBL_MAX && tabuLength_ < numFlips_ - changed_[var]) {
						if (improvement >= bestvalue)
						{ // Tied or better than previous best
							if (improvement > bestvalue)
							{
								numbest = 0;
								best.clear();
							}
							bestvalue = improvement;
							best.append(var);
							numbest++;
						}
					}
				}
			}

			if (numbest == 0)
			{
				if (hmwsdebug) cout << "Leaving HMaxWalkSat::pickHMWS (NOVALUE)" << endl;
				return NOVALUE;
			}

			int toFlip = NOVALUE;
			// Exactly one best atom
			if (numbest == 1)
				toFlip = best[0];
			else
				// Choose one of the best at random
				toFlip = best[random()%numbest];

			if (hmwsdebug) 
			{
				cout << "Picked dis atom to flip: " << toFlip << ",";hstate_->printDisAtom(toFlip, cout);
				cout << ". Improvement:" <<bestvalue <<endl;
			}

			return toFlip;
		}  // Discrete clause
		else  // Hybrid clause, clause idx = - (toFix + 1); toFix < 0
		{
			toFix++;
			int contClauseIdx = -toFix;
			if (hstate_->hybridWts_[contClauseIdx] == 0) return NOVALUE;
			if (hmwsdebug) {
				cout << "Pick hybrid clause: " << contClauseIdx << ", constraint:";
				hstate_->printHybridConstraint(contClauseIdx, cout);
			}
			
			int contVarNum =  hstate_->hybridContClause_[contClauseIdx].size();
			int disVarNum = hstate_->hybridDisClause_[contClauseIdx].size();
			int totalVarNum = contVarNum + disVarNum;
			if (noisyPick) //random pick and then flip, p = numerator_/denominator_;
			{
				set<int> pickedContVars;
				int inIdx = random() % totalVarNum;
				if (inIdx < contVarNum)  // continuous variable
				{
					int contAtomIdx = hstate_->hybridContClause_[contClauseIdx][inIdx];
					// 1. If the discrete part of this hybrid constraint is 0, then do nth.
					// 2. If the discrete part of this hybrid constraint is 1, then return the continuous variable (Optimize this hybrid constraint along this variable direction w./ L-BFGS).
					if (hstate_->HybridClauseDisPartValue(contClauseIdx)) {
						// TODO: optimize the constraint along this variable to its optimal state.
						// We need to compute the assignment to the continuous variable corresponding to the optimal state.
						pickedContValue_ = hstate_->OptimizeHybridClauseToContVar(contClauseIdx, inIdx);
						return -contAtomIdx;
					}
					else return NOVALUE;
				}
				else  // discrete variable
				{
					// Check if flipping this discrete variable could improve the current hybrid clause's weighted contribution;
					// If so, return it, if not, do nothing.
					int disAtomIdx = abs(hstate_->hybridDisClause_[contClauseIdx][inIdx - contVarNum]);					
					double hybridClauseValueBeforeFlipping = hstate_->HybridClauseValue(contClauseIdx);
					hstate_->atom_[disAtomIdx] = !hstate_->atom_[disAtomIdx];
					double hybridClauseValueAfterFlipping = hstate_->HybridClauseValue(contClauseIdx);
					hstate_->atom_[disAtomIdx] = !hstate_->atom_[disAtomIdx];
					if (hybridClauseValueAfterFlipping - hybridClauseValueBeforeFlipping >= 0) return disAtomIdx;
					else return NOVALUE;
				}
			}  // Noisy pick
			else //greedy try to flip every variable and find the one get biggest improvement
			{
				long double bestvalue = -LDBL_MAX;
				Array<int> best;
				int numbest = 0;
				map<int, double> pickedContValues;

				for (int i = 0; i < contVarNum; ++i)
				{
					int contAtomIdx = hstate_->hybridContClause_[contClauseIdx][i];
					double assign = 0;
					// compute the improvement by optimizing the contribution of contAtomIdx.
					double improvement = hstate_->ReduceClauseAndOptimize(contAtomIdx, &assign);
					if (improvement == NOVALUE) continue;
					if (improvement >= bestvalue)
					{ // Tied or better than previous best
						if (improvement > bestvalue)
						{
							numbest = 0;
							best.clear();
						}
						bestvalue = improvement;
						best.append(-contAtomIdx);
						numbest++;
						pickedContValues.insert(map<int, double>::value_type(-contAtomIdx, assign));
					}
				}
				for (int i = 0; i < disVarNum; ++i) {
					int disAtomIdx = hstate_->hybridDisClause_[contClauseIdx][i];
					double improvement = calculateImprovementDisHMWS(disAtomIdx);
					if (improvement >= bestvalue)
					{ // Tied or better than previous best
						if (improvement > bestvalue)
						{
							numbest = 0;
							best.clear();
						}
						bestvalue = improvement;
						best.append(disAtomIdx);
						numbest++;
					}
				}

				int atomIdx = NOVALUE;
				// Exactly one best atom
				if (numbest == 1)
					atomIdx = best[0];
				else
					// Choose one of the best at random
					atomIdx = best[random()%numbest];

				if (atomIdx < 0)  // Picked variable is a continuous variable.
				{
					pickedContValue_ = pickedContValues[atomIdx];
				}
				return atomIdx;
			}  // Greedy
		}  // Hybrid clause.
		return NOVALUE;
	}

	int pickSA()
	{
		//	if (hmwsdebug) cout << "Entering HMaxWalkSat::pickSA" << endl;
		if (hmwsdebug)
		{
			cout << "SA" << endl;

		}

		// Choose a random atom to flip
		int toFlip = hstate_->getIndexOfRandomAtom(); //change to hybrid, pos, discrete, neg, continuous
		if (toFlip > 0)//discrete atom
		{
			if (hmwsdebug)
			{
				cout << "Choose to flip dis atom " << toFlip << ":";hstate_->printDisAtom(toFlip, cout);
			}
			long double improvement = hstate_->getImprovementInWeightSumByFlipping(toFlip);
			
			if (hmwsdebug)
			{
				cout << "By flipping, improvement is " << improvement << endl;
			}

			// If pos. or no improvement, then the atom will be flipped
			// Or neg. improvement: According to temp. flip atom anyway
			if ((improvement >= 0) ||
				(random() <= exp(improvement/(saTemp_)) * RAND_MAX))
			{
				// If atom can not be flipped, then null flip
				saSteps_ ++;
				if (saSteps_ % saInterval_ == 0)
				{
					saTemp_ = saTemp_ * saTempDownRatio_;
				}
				return toFlip;
			}
			else
			{
				return NOVALUE;
			}
		}
		else //continuous atom
		{
			//generate the new value for the continuous variable according to the proposal distribution
			int contAtomIdx = -toFlip;
			double vContAtomFlipped = hstate_->getContVarValueByRandomMove(contAtomIdx);	
			//compute the proposed cont value here

			long double improvement = hstate_->getImprovementRandomMoveCont(contAtomIdx, vContAtomFlipped); // change to hybrid

			if (hmwsdebug)
			{
				cout << "Choose to flip cont atom " << contAtomIdx << ":";hstate_->printContAtom(contAtomIdx, cout);
			}
			if (hmwsdebug)
			{
				cout << "By flipping to " << vContAtomFlipped << ", improvement is " << improvement << endl;
			}
			if ((improvement >= 0) || (random() <= exp(improvement/(saTemp_)) * RAND_MAX))
			{	
				pickedContValue_ = vContAtomFlipped; //the value picked.
				saSteps_ ++;
				if (saSteps_ % saInterval_ == 0)
				{
					saTemp_ = saTemp_ * saTempDownRatio_;
				}
				return toFlip; //neg, indicate cont variable.
			} else {
				return NOVALUE;
			}
		}	
		//need to update the temperature here
	}

 /**
  * Pick an atom according to the SampleSat heuristic. This means sometimes
  * sim. annealing, sometimes HMaxWalkSat.
  * 
  * @return Index of atom picked.
  */
  int pickSS()
  {
    long double costOfFalseClauses = hstate_->getCostOfTotalFalseConstraints();
      // If we have already reached a solution or if in an SA step,
      // then perform SA
    if (costOfFalseClauses <= targetCost_ + SMALLVALUE ||
        (random() % 100 < saRatio_ && !lateSa_))
    {
      if (hmwsdebug)
        cout << "SASTEP" << endl;

        // Choose a random atom to flip
        // change to hybrid, pos, discrete, neg, continuous
      int toFlip = hstate_->getIndexOfRandomAtom();
      if (toFlip > 0)  // Discrete atom.
      {
        if (hmwsdebug)
        {
          cout << "Choose to flip dis atom " << toFlip << ":";
          hstate_->printDisAtom(toFlip, cout);
          cout << endl;
        }
        long double improvement = calculateImprovementDisHMCS(toFlip); 

        if (hmwsdebug)
          cout << "By flipping, improvement is " << improvement << endl;

          // If pos. or no improvement, then the atom will be flipped
          // Or neg. improvement: According to temp. flip atom anyway
        if ((improvement >= 0) || 
            (random() <= exp(improvement/(saTemp_/100.0)) * RAND_MAX))
        {
          return toFlip;
        }
        else
        {
          return NOVALUE;
        }
      }
      else  // Continuous variable.
      {
          //generate the new value for the continuous variable according to the
          //proposal distribution
        int contAtomIdx = -toFlip;
        double vContAtomFlipped;				
        long double improvement =
          hstate_->GetImprovementByMovingContVar(contAtomIdx, vContAtomFlipped); 
				
        if (hmwsdebug)
        {
          cout << "Choose to move cont atom " << contAtomIdx << ":";
          hstate_->printContAtom(contAtomIdx, cout);
          cout << " To " << vContAtomFlipped << ", improvement is "
               << improvement << endl;
        }

        if ((improvement >= 0) ||
            (random() <= exp(improvement/(saTemp_/100.0)) * RAND_MAX))
        {	
          pickedContValue_ = vContAtomFlipped; //the value picked.
          return toFlip; //neg, indicate cont variable.
        }	  
        else
        {
          return NOVALUE;
        }
      }	
    }
      // Not in a solution or SA step: perform Hybrid WalkSAT step
    else
    {
      return pickHMCS();  // TODO: juewang
    }
  }

	/**
	* Calculates the improvement (makecost - breakcost) by flipping an atom.
	* If the atom is in a block, then its index is saved in candidateBlockedVars
	* and the index of another atom chosen to be flipped in the block is
	* appended to othersToFlip. If the atom is in a block with evidence, then
	* its index is appended to canNotFlip.
	* 
	* @param atomIdx Index of atom for which the improvement is calculated.
	* @param canNotFlip Holds indices of atoms which can not be flipped due
	* to evidence in a block.
	* @param candidateBlockedVars If dealing with an atom in a block, then its
	* index is appended here.
	* @param othersToFlip If dealing with an atom in a block, then the index of
	* the other atom chosen to be flipped is appended here.
	* 
	*/
	long double calculateImprovementDisHMCS(const int& atomIdx)	//should return the number of constraints satisfied by flipping minus the number of constraints violated by flipping
	{
		return hstate_->getImprovementByFlippingDisHMCS(atomIdx);
	}


	// No consideration of handling lazy and blocking is taken into account.
	// The improvement is computed based on the actual weighted feature values.
	long double calculateImprovementDisHMWS(int atomIdx)
	{
		if (atomIdx > 0) //discrete
		{
			return hstate_->getImprovementByFlippingDisHMWS(atomIdx);
		}
		return 0;	
	}

	/**
	* Reconstructs the state with the lowest cost by flipping atoms back to
	* their value in this state.
	*/
	void reconstructLowState()
	{
		assert(lazyLowState_);
		if (hmwsdebug) cout << "Reconstructing low state ..." << endl;
		hash_set<int>::const_iterator it = varsFlippedSinceLowState_.begin();
		for (; it != varsFlippedSinceLowState_.end(); it++)
		{
			if (hmwsdebug)
			{
				cout << "Flipping atom " << (*it) << endl;
			}
			hstate_->flipAtomValue(*it, -1);
		}
		hstate_->saveLowState();
		if (hmwsdebug) cout << "Done reconstructing low state ..." << endl;
	}

private:
	int saInterval_;
	unsigned int saSteps_;
	double maxSeconds_;
	double saTempDownRatio_;
	double saTempInit_;
	double pickedContValue_;
	int pickedContIdx_;
	//Array<int> pick2VarIdx_;
	Array<int> pickDisIdx_;
	//Array<double> pick2VarValue_;

	////////// BEGIN: User parameters ///////////
	// Heuristic to be used to pick an atom
	int heuristic_;
	// At least this many flips must occur before flipping the same atom
	int tabuLength_;

	// BEGIN: SampleSat parameters
	// Percent of sim. annealing steps
	int saRatio_;
	// Sim. annealing temperature
	double saTemp_;
	// Run sim. annealing only at a plateur
	bool lateSa_;
	// END: SampleSat parameters
	////////// END: User parameters ///////////

	// Function pointer holds which function is to be used to pick an atom
	// = {pickRandom, pickBest, pickTabu, pickSS};
	int (HMaxWalkSat::*pickcode[15])(void);
	// If true, never break a highest cost clause
	bool hard_;

	// Make random flip with numerator/denominator probability
	int numerator_;
	int denominator_;

	// Which other atom to flip in block
	bool inBlock_;
	int flipInBlock_;

	// If false, the naive way of saving low states (each time a low state is
	// found, the whole state is saved) is used; otherwise, a list of variables
	// flipped since the last low state is kept and the low state is
	// reconstructed. This can be much faster for very large data sets.
	bool lazyLowState_;
	// List of variables flipped since last low state was found. This is used
	// to reconstruct the low state when lazyLowState_ is set.
	hash_set<int> varsFlippedSinceLowState_;
	//unsigned falseDisNum_;
};

#endif
