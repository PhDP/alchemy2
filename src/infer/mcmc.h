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
#ifndef MCMC_H_
#define MCMC_H_

#include "inference.h"
#include "mcmcparams.h"

  // Set to true for more output
const bool mcmcdebug = false;

/**
 * Superclass of all MCMC inference algorithms. This class does not implement
 * all pure virtual functions of Inference and is thus an abstract class.
 */
class MCMC : public Inference
{
 public:

  /**
   * Constructor. User-set parameters are set.
   * 
   * @see Inference#Constructor(VariableState*, long int, const bool&,
   *                            GibbsParams*)
   */  
  MCMC(VariableState* state, long int seed, const bool& trackClauseTrueCnts,
       MCMCParams* params, Array<Array<Predicate* >* >* queryFormulas = NULL)
    : Inference(state, seed, trackClauseTrueCnts, queryFormulas)
  {
      // User-set parameters
    numChains_ = params->numChains;
    burnMinSteps_ = params->burnMinSteps;
    burnMaxSteps_ = params->burnMaxSteps;
    minSteps_ = params->minSteps;
    maxSteps_ = params->maxSteps;
    maxSeconds_ = params->maxSeconds;
  }


  MCMC(HVariableState* state, long int seed, const bool& trackClauseTrueCnts,
	  MCMCParams* params)
	  : Inference(state, seed, trackClauseTrueCnts)
  {
	  // User-set parameters
	  numChains_ = params->numChains;
	  burnMinSteps_ = params->burnMinSteps;
	  burnMaxSteps_ = params->burnMaxSteps;
	  minSteps_ = params->minSteps;
	  maxSteps_ = params->maxSteps;
	  maxSeconds_ = params->maxSeconds;
  }

  /**
   * Destructor.
   */
  ~MCMC() {}

  /*double computeHybridClauseValue(int contClauseIdx)
  {
	  return hstate_->hybridClauseCost_[contClauseIdx] * HybridClauseContPartValue(contClauseIdx) * HybridClauseDisPartValue(contClauseIdx);
  }*/

  double computeHybridClauseValue(int clauseIdx, int c)
  {
	  double contClauseContPartValue = HybridClauseContPartValue(clauseIdx, c);
	  double contClauseDisPartValue = HybridClauseDisPartValue(clauseIdx, c);
	  return hstate_->hybridWts_[clauseIdx] * contClauseContPartValue * contClauseDisPartValue;
  }

  double HybridClauseContPartValue(int contClauseIdx, int c)
  {
	  PolyNomial& pl = hstate_->GetHybridClausePolynomial(contClauseIdx);

	  assert(hstate_->hybridContClause_[contClauseIdx].size() == pl.GetVarNum());

	  Array<double> arVar;
	  for(int i = 0; i < hstate_->hybridContClause_[contClauseIdx].size(); ++i)
	  {
		  arVar.append(truthValuesCont_[hstate_->hybridContClause_[contClauseIdx][i] - 1][c]); // the order of the cont atom here should be the same as they are in the Pl
	  }
	  double v = pl.ComputePlValue(arVar);
	  
	  return v;
  }

  double HybridClauseDisPartValue(int contClauseIdx, int c)
  {
	  bool bAndOr = hstate_->hybridConjunctionDisjunction_[contClauseIdx];	  
	  int numTrueLits = 0;
	  int numFalseLits = 0;
	  for(int j = 0; j < hstate_->hybridDisClause_[contClauseIdx].size(); ++j)
	  {
		  int atomIdx = hstate_->hybridDisClause_[contClauseIdx][j];
		  if((atomIdx > 0) == truthValues_[abs(atomIdx)-1][c]) // true literal
		  {
			  numTrueLits ++;
			  if(!bAndOr) // disjunctive clause
			  {
				  break;
			  }
		  }
		  else // false literal
		  {
			  numFalseLits ++;
			  if(bAndOr) // conjunctive clause
			  {
				  break;
			  }
		  }
	  }
	  if(!bAndOr) // disjunctive clause
	  {
		  return (numTrueLits > 0)?1.0:0.0;
	  }
	  else // conjunctive clause
	  {
		  return (numFalseLits > 0)?0.0:1.0;
	  }
	  return 0.0;
  }

  void updateDisPredValue(int predIdx, int chainIdx, bool updateValue)
  {
	  bool bBak = truthValues_[predIdx][chainIdx];
	  Array<int>& occDisPos = hstate_->getPosOccurenceArray(predIdx + 1);
	  Array<int>& occDisNeg = hstate_->getNegOccurenceArray(predIdx + 1);
	  Array<int>& occCont = hstate_->getContDisOccurenceArray(predIdx + 1) ;

	  if(updateValue)
	  {
		  if(!bBak)
		  {
			  // update numTrueLits of discrete clauses where it appears 
			  for(int j = 0; j < occDisPos.size(); ++j)
			  {
				  int disClauseIdx = occDisPos[j];
				  numTrueLits_[disClauseIdx][chainIdx] += 1;
			  }

			  for( int j = 0; j < occDisNeg.size(); ++j)
			  {
				  int disClauseIdx = occDisNeg[j];
				  numTrueLits_[disClauseIdx][chainIdx] -= 1;
			  }
		  }
	  }
	  else
	  {
		  if(bBak)
		  {
			  for( int j = 0; j < occDisPos.size(); ++j)
			  {
				  int disClauseIdx = occDisPos[j];
				  numTrueLits_[disClauseIdx][chainIdx] -= 1;
			  }

			  for( int j = 0; j < occDisNeg.size(); ++j)
			  {
				  int disClauseIdx = occDisNeg[j];
				  numTrueLits_[disClauseIdx][chainIdx] += 1;
			  }
		  }	
	  }

	  truthValues_[predIdx][chainIdx] = updateValue;

	  if(bBak != updateValue)
	  {
		  for( int j = 0; j < occCont.size(); ++j)
		  {
			  int hybridClauseIdx = occCont[j];
			  hybridClauseDisPartValueMCMC_[hybridClauseIdx] = HybridClauseDisPartValue(hybridClauseIdx, chainIdx)==1.0?true:false;
		  }
	  }  
  }

  void updateProposalContValue(int contPredIdx, int chainIdx)
  {
	  Array<int>& occContCont = hstate_->hybridContOccurrence_[contPredIdx + 1];
	  PolyNomial pl;
	  for(int j = 0; j < occContCont.size(); j++)
	  {
		  int contClauseIdx = occContCont[j];
		  PolyNomial pltmp = hstate_->hybridPls_[contClauseIdx];
		  //find the in idx here
		  int inIdx = -1;
		  Array<double> arVars;
		  for(int k = 0; k < hstate_->hybridContClause_[contClauseIdx].size(); k++)
		  {
			  arVars.append(truthValuesCont_[hstate_->hybridContClause_[contClauseIdx][k]-1][chainIdx]);
			  if (hstate_->hybridContClause_[contClauseIdx][k] == contPredIdx+1)
			  {
				  //pltmp.ReduceToOneVar(,k)
				  inIdx = k;
			  }			  
		  }
		  if (inIdx == -1)
		  {
			  cout <<  "faint" << endl;
		  }			  

		  pltmp.ReduceToOneVar(arVars, inIdx);
		  pl.AddPl(pltmp);
	  }
	  
	  double miu = 0, stdev = 0;
	  pl.GetGaussianPara(&miu, &stdev);
	  truthValuesCont_[contPredIdx][chainIdx] = ExtRandom::gaussRandom(miu, stdev);

	  for(int j = 0; j < occContCont.size(); j++)
	  {
		  int hybridClauseIdx = occContCont[j];
		  hybridClauseContPartValueMCMC_[hybridClauseIdx][chainIdx] = HybridClauseContPartValue(hybridClauseIdx,chainIdx);
	  }  
  }

  // Computes the probability of a ground discrete predicate in a chain, adapted to accommodate the hybrid case.
  double getProbabilityOfPredH(const int& predIdx, const int& chainIdx, const double& invTemp)
  {
	  // Different for multi-chain
	  if (numChains_ > 1)
	  {
		  double wtDisAsTrue = 0, wtDisAsFalse = 0;
		  bool bBak = truthValues_[predIdx][chainIdx];

		  Array<int>& occDisPos = hstate_->getPosOccurenceArray(predIdx + 1);
		  Array<int>& occDisNeg = hstate_->getNegOccurenceArray(predIdx + 1);
		  Array<int>& occCont = hstate_->getContDisOccurenceArray(predIdx + 1) ;

		  // for clauses where the variable appears as a true literal
		  for(int j = 0; j < occDisPos.size(); j++)
		  {
			  int disClauseIdx = occDisPos[j];
			  wtDisAsTrue += hstate_->clauseCost_[disClauseIdx];

			  if(numTrueLits_[disClauseIdx][chainIdx] > 1)
			  {
				  wtDisAsFalse += hstate_->clauseCost_[disClauseIdx];
			  }
			  else if(numTrueLits_[disClauseIdx][chainIdx] == 1 && !bBak) // there was only 1 true literal and the only true literal is some one else
			  {
				  wtDisAsFalse += hstate_->clauseCost_[disClauseIdx];
			  }
		  }

		  // for clauses where the variable appears as a false literal
		  for(int j = 0; j < occDisNeg.size(); j++)
		  {
			  int disClauseIdx = occDisNeg[j];
			  if(numTrueLits_[disClauseIdx][chainIdx] > 1) // there were more than one true literals 
			  {
				  wtDisAsTrue += hstate_->clauseCost_[disClauseIdx];
			  }
			  else if(numTrueLits_[disClauseIdx][chainIdx] == 1 && bBak) // there was only 1 true literal and the only true literal is some one else
			  {
				  wtDisAsTrue += hstate_->clauseCost_[disClauseIdx];
			  }
			  wtDisAsFalse += hstate_->clauseCost_[disClauseIdx];
		  }

		  for (int j = 0; j < occCont.size(); j++) {
			  int hybridClauseIdx = occCont[j];
			  double contPart = hybridClauseContPartValueMCMC_[hybridClauseIdx][chainIdx];
			  double wt = hstate_->hybridWts_[hybridClauseIdx];

			  truthValues_[predIdx][chainIdx] = true;
			  double disPartAsTrue = HybridClauseDisPartValue(hybridClauseIdx,chainIdx);
			  truthValues_[predIdx][chainIdx] = false;
			  double disPartAsFalse = HybridClauseDisPartValue(hybridClauseIdx,chainIdx);
			  truthValues_[predIdx][chainIdx] = bBak;

			  wtDisAsTrue += wt*disPartAsTrue*contPart;
			  wtDisAsFalse += wt*disPartAsFalse*contPart;
		  }

		  // get probabilities
		  double wtDiff = (wtDisAsFalse - wtDisAsTrue) * invTemp;
		  double prob;
		  if (wtDiff > 403.429)
		  {
			  prob = 0;
		  }
		  else if (wtDiff < -403.429)
		  {
			  prob = 1;
		  }
		  else
		  {
			  prob = 1 / ( 1 + exp(wtDiff));
		  }
		  return prob;
	  }  // For chain id ends.
	  else
	  {
		  GroundPredicate* gndPred = hstate_->getGndPred(predIdx);
		  double wtDisAsTrue = 0, wtDisAsFalse = 0;
		  bool bBak = gndPred->getTruthValue();

		  Array<int>& occDisPos = hstate_->getPosOccurenceArray(predIdx + 1);
		  Array<int>& occDisNeg = hstate_->getNegOccurenceArray(predIdx + 1);
		  Array<int>& occCont = hstate_->getContDisOccurenceArray(predIdx + 1) ;

		  gndPred->setTruthValue(true);
		  //hstate_->setValueOfAtom(i+1, true);		  
		  for(int j = 0; j < occDisPos.size(); j++)
		  {
			  int disClauseIdx = occDisPos[j];
			  wtDisAsTrue += hstate_->clauseCost_[disClauseIdx];			  
		  }

		  for(int j = 0; j < occDisNeg.size(); j++)
		  {
			  int disClauseIdx = occDisNeg[j];
			  for(int k = 0; k < hstate_->clause_[disClauseIdx].size(); k++)
			  {
				  int lit = hstate_->clause_[disClauseIdx][k];
				  GroundPredicate* gndPredtmp =	 hstate_->getGndPred(abs(lit)-1);
				  if ((lit > 0) == gndPredtmp->getTruthValue()) //true literal
				  {
					  wtDisAsTrue += hstate_->clauseCost_[disClauseIdx];
					  break;
				  }
			  }
		  }

		  for (int j = 0; j < occCont.size(); j++)
		  {
			  int contClauseIdx = occCont[j];
			  wtDisAsTrue += hstate_->HybridClauseValue(contClauseIdx);
				  // computeHybridClauseValue(contClauseIdx);
		  }


		  gndPred->setTruthValue(false);
		  for(int j = 0; j < occDisNeg.size(); j++)
		  {
			  int disClauseIdx = occDisNeg[j];
			  wtDisAsFalse += hstate_->clauseCost_[disClauseIdx];
		  }

		  for(int j = 0; j < occDisPos.size(); j++)
		  {
			  int disClauseIdx = occDisPos[j];
			  for(int k = 0; k < hstate_->clause_[disClauseIdx].size(); k++)
			  {
				  int lit = hstate_->clause_[disClauseIdx][k];
				  GroundPredicate* gndPredtmp =	 hstate_->getGndPred(abs(lit)-1);
				  if ((lit > 0) == gndPredtmp->getTruthValue()) //true literal
				  {
					  wtDisAsFalse += hstate_->clauseCost_[disClauseIdx];
					  break;
				  }
			  }
		  }

		  for (int j = 0; j < occCont.size(); j++)
		  {
			  int contClauseIdx = occCont[j];
			  wtDisAsFalse += hstate_->HybridClauseValue(contClauseIdx);
		  }

		  gndPred->setTruthValue(bBak);

		  return 1.0 / ( 1.0 + exp((wtDisAsFalse - wtDisAsTrue) * invTemp));
	  }
  }

  /**
   * Prints out the network.
   */
  virtual void printNetwork(ostream& out)
  {
    state_->printNetwork(out);
  } 

  /**
   * Prints the probabilities of each predicate to a stream.
   */
  void printProbabilities(ostream& out)
  {
    for (int i = 0; i < state_->getNumAtoms(); i++)
    {
      double prob = getProbTrue(i);

        // Uniform smoothing
      prob = (prob*10000 + 1/2.0)/(10000 + 1.0);
      state_->printGndPred(i, out);
      out << " " << prob << endl;
    }    
  }

  /**
   * Puts the predicates whose probability has changed with respect to the
   * reference vector oldProbs by more than probDelta in string form and the
   * corresponding probabilities of each predicate in two vectors.
   * 
   * @param changedPreds Predicates whose probability have changed more than
   * probDelta are put here.
   * @param probs The probabilities corresponding to the predicates in
   * changedPreds are put here.
   * @param oldProbs Reference probabilities for checking for changes.
   * @param probDelta If probability of an atom has changed more than this
   * value, then it is considered to have changed.
   */
  void getChangedPreds(vector<string>& changedPreds, vector<float>& probs,
                       vector<float>& oldProbs, const float& probDelta)
  {
    changedPreds.clear();
    probs.clear();
    int numAtoms = state_->getNumAtoms();
      // Atoms may have been added to the state, previous prob. was 0
    oldProbs.resize(numAtoms, 0);
    for (int i = 0; i < numAtoms; i++)
    {
      double prob = getProbTrue(i);
      if (abs(prob - oldProbs[i]) > probDelta)
      {
          // Truth value has changed: Store new value (not smoothed) in oldProbs
          // and add to two return vectors
        oldProbs[i] = prob;
          // Uniform smoothing
        prob = (prob*10000 + 1/2.0)/(10000 + 1.0);
        ostringstream oss(ostringstream::out);
        state_->printGndPred(i, oss);
        changedPreds.push_back(oss.str());
        probs.push_back(prob);
      }
    }    
  }


  double getProbabilityH(GroundPredicate* const& gndPred)
  {
	  int idx = hstate_->getGndPredIndex(gndPred);
	  double prob = 0.0;
	  if (idx >= 0) prob = getProbTrue(idx);
	  // Uniform smoothing
	  return (prob*10000 + 1/2.0)/(10000 + 1.0);
  }

  /**
   * Gets the probability of a ground predicate.
   * 
   * @param gndPred GroundPredicate whose probability is being retrieved.
   * @return Probability of gndPred if present in state, otherwise 0.
   */
  double getProbability(GroundPredicate* const& gndPred)
  {
    int idx = state_->getGndPredIndex(gndPred);
    double prob = 0.0;
    if (idx >= 0) prob = getProbTrue(idx);
      // Uniform smoothing
    return (prob*10000 + 1/2.0)/(10000 + 1.0);
  }

  /**
   * Prints each predicate with a probability of 0.5 or greater to a stream.
   */
  void printTruePreds(ostream& out)
  {
    for (int i = 0; i < state_->getNumAtoms(); i++)
    {
      double prob = getProbTrue(i);

        // Uniform smoothing
      prob = (prob*10000 + 1/2.0)/(10000 + 1.0);
      if (prob >= 0.5) state_->printGndPred(i, out);
    }    
  }

  void printTruePredsH(ostream& out)
  {
	  for (int i = 0; i < hstate_->getNumAtoms(); i++)
	  {
		  double prob = getProbTrue(i);

		  // Uniform smoothing
		  prob = (prob * 10000 + 1/2.0) / (10000 + 1.0);
		  if (prob >= 0.5) hstate_->printGndPred(i, out);
	  }    
  }


 protected:

  /**
   * Initializes truth values and weights in the ground preds.
   * 
   * @param numChains Number of chains for which the values should be
   * initialized.
   * @param startIdx All predicates with index greater than or equal to this
   * will be initialized.
   */
  void initTruthValuesAndWts(const int& numChains)
  {
    int numPreds = state_->getNumAtoms();
    truthValues_.growToSize(numPreds);
    wtsWhenFalse_.growToSize(numPreds);
    wtsWhenTrue_.growToSize(numPreds);
    for (int i = 0; i < numPreds; i++)
    {
      truthValues_[i].growToSize(numChains, false);
      wtsWhenFalse_[i].growToSize(numChains, 0);
      wtsWhenTrue_[i].growToSize(numChains, 0);
    }
    
    int numClauses = state_->getNumClauses();
    numTrueLits_.growToSize(numClauses);
    for (int i = 0; i < numClauses; i++)
    {
      numTrueLits_[i].growToSize(numChains, 0);
    }
  }

  /**
   * Initializes structure for holding number of times a predicate was set
   * to true.
   */
  void initNumTrue()
  {
    int numPreds = state_->getNumAtoms();
    numTrue_.growToSize(numPreds);
    for (int i = 0; i < numTrue_.size(); i++)
      numTrue_[i] = 0;
  }

  /**
   * Initializes the number of true lits in each clause in each chain.
   * 
   * @param numChains Number of chains for which the initialization should
   * take place.
   */
  void initNumTrueLits(const int& numChains)
  {
      // Single chain
    if (numChains == 1) state_->resetMakeBreakCostWatch();
    for (int i = 0; i < state_->getNumClauses(); i++)
    {
      GroundClause* gndClause = state_->getGndClause(i);
      for (int j = 0; j < gndClause->getNumGroundPredicates(); j++)
      {
        const int atomIdx = abs(state_->getAtomInClause(j, i)) - 1;
        const bool sense = gndClause->getGroundPredicateSense(j);
        if (numChains > 1)
        {
          for (int c = 0; c < numChains; c++)
          {
            if (truthValues_[atomIdx][c] == sense)
            {
              numTrueLits_[i][c]++;
              assert(numTrueLits_[i][c] <= state_->getNumAtoms());
            }
          }
        }
        else
        { // Single chain
          GroundPredicate* gndPred = state_->getGndPred(atomIdx);
          if (gndPred->getTruthValue() == sense)
            state_->incrementNumTrueLits(i);
          assert(state_->getNumTrueLits(i) <= state_->getNumAtoms());
        }        
      }
    }
  }
 
  /**
   * Randomly initializes the ground predicate truth values, taking blocks
   * into account.
   * 
   * @param numChains Number of chains for which the initialization should
   * take place.
   */
  void randomInitGndPredsTruthValues(const int& numChains)
  {
    for (int c = 0; c < numChains; c++)
    {
      if (mcmcdebug) cout << "Chain " << c << ":" << endl;
        // For each block: select one to set to true
      for (int i = 0; i < state_->getDomain()->getNumPredBlocks(); i++)
      {
          // If evidence atom exists, then all others are false
        if (state_->getDomain()->getBlockEvidence(i))
        {
            // If 2nd argument is -1, then all are set to false
          setOthersInBlockToFalse(c, -1, i);
          continue;
        }

        bool ok = false;
        while (!ok)
        {
          const Predicate* pred = state_->getDomain()->getRandomPredInBlock(i);
          GroundPredicate* gndPred = new GroundPredicate((Predicate*)pred);
          int idx = state_->getIndexOfGroundPredicate(gndPred);
            
          delete gndPred;
          delete pred;

          if (idx >= 0)
          {
              // Truth values are stored differently for multi-chain
            if (numChains_ > 1)
              truthValues_[idx][c] = true;
            else
            {
              GroundPredicate* gndPred = state_->getGndPred(i);
              gndPred->setTruthValue(true);
            }
            setOthersInBlockToFalse(c, idx, i);
            ok = true;
          }
        }
      }
      
        // Random tv for all not in blocks
      for (int i = 0; i < truthValues_.size(); i++)
      {
          // Predicates in blocks have been handled above
        if (state_->getBlockIndex(i) == -1)
        {
          bool tv = genTruthValueForProb(0.5);
            // Truth values are stored differently for multi-chain
          if (numChains_ > 1)
            truthValues_[i][c] = tv;
          else
          {
            GroundPredicate* gndPred = state_->getGndPred(i);
            gndPred->setTruthValue(tv);
          }
          if (mcmcdebug) cout << "Pred " << i << " set to " << tv << endl;
        }
      }
    }
  }

  /**
   * Generates a truth value based on a probability.
   * 
   * @param p Number between 0 and 1. With probability p, truth value will be
   * true and with probability 1 - p, it will be false.
   */
  bool genTruthValueForProb(const double& p)
  {
    if (p == 1.0) return true;
    if (p == 0.0) return false;
    bool r = random() <= p*RAND_MAX;
    return r;
  }

  /**
   * Computes the probability of a ground predicate in a chain.
   * 
   * @param predIdx Index of predicate.
   * @param chainIdx Index of chain.
   * @param invTemp InvTemp used in simulated tempering.
   * 
   * @return Probability of predicate.
   */
  double getProbabilityOfPred(const int& predIdx, const int& chainIdx,
                              const double& invTemp)
  {
      // Different for multi-chain
    if (numChains_ > 1)
    {
      return 1.0 /
             ( 1.0 + exp((wtsWhenFalse_[predIdx][chainIdx] - 
                          wtsWhenTrue_[predIdx][chainIdx]) *
                          invTemp));      
    }
    else
    {
      GroundPredicate* gndPred = state_->getGndPred(predIdx);
      return 1.0 /
             ( 1.0 + exp((gndPred->getWtWhenFalse() - 
                          gndPred->getWtWhenTrue()) *
                          invTemp));
    }
  }
 
  /**
   * Sets the truth values of all atoms for a given chain in a block except
   * for the one given.
   * 
   * @param chainIdx Index of chain for which atoms are being set.
   * @param atomIdx Index of atom in block exempt from being set to false.
   * @param blockIdx Index of block whose atoms are set to false.
   */
  void setOthersInBlockToFalse(const int& chainIdx, const int& atomIdx,
                               const int& blockIdx)
  {
    int blockSize = state_->getDomain()->getBlockSize(blockIdx);
    for (int i = 0; i < blockSize; i++)
    {
      const Predicate* pred = state_->getDomain()->getPredInBlock(i, blockIdx);
      GroundPredicate* gndPred = new GroundPredicate((Predicate*)pred);
      int idx = state_->getIndexOfGroundPredicate(gndPred);

      delete gndPred;
      delete pred;

        // Pred is in the state
      if (idx >= 0 && idx != atomIdx)
        truthValues_[idx][chainIdx] = false;
    }
  }

  /**
   * Performs one step of Gibbs sampling in one chain.
   * 
   * @param chainIdx Index of chain in which the Gibbs step is performed.
   * @param burningIn If true, burning-in is occuring. Otherwise, false.
   * @param affectedGndPreds Used to store GroundPredicates which are affected
   * by changing truth values.
   * @param affectedGndPredIndices Used to store indices of GroundPredicates
   * which are affected by changing truth values.
   */
  void performGibbsStep(const int& chainIdx, const bool& burningIn,
                        GroundPredicateHashArray& affectedGndPreds,
                        Array<int>& affectedGndPredIndices)
  {
    if (mcmcdebug) cout << "Gibbs step" << endl;

      // For each block: select one to set to true
    for (int i = 0; i < state_->getDomain()->getNumPredBlocks(); i++)
    {
        // If evidence atom exists, then all others stay false
      if (state_->getDomain()->getBlockEvidence(i)) continue;

      int chosen = gibbsSampleFromBlock(chainIdx, i, 1);
        // Truth values are stored differently for multi-chain
      bool truthValue;
      const Predicate* pred = state_->getDomain()->getPredInBlock(chosen, i);
      GroundPredicate* gndPred = new GroundPredicate((Predicate*)pred);
      int idx = state_->getIndexOfGroundPredicate(gndPred);

      delete gndPred;
      delete pred;
      
        // If gnd pred in state:
      if (idx >= 0)
      {
        gndPred = state_->getGndPred(idx);
        if (numChains_ > 1) truthValue = truthValues_[idx][chainIdx];
        else truthValue = gndPred->getTruthValue();
          // If chosen pred was false, then need to set previous true
          // one to false and update wts
        if (!truthValue)
        {
          int blockSize = state_->getDomain()->getBlockSize(i);
          for (int j = 0; j < blockSize; j++)
          {
              // Truth values are stored differently for multi-chain
            bool otherTruthValue;
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
              otherGndPred = state_->getGndPred(otherIdx);
              if (numChains_ > 1)
                otherTruthValue = truthValues_[otherIdx][chainIdx];
              else
                otherTruthValue = otherGndPred->getTruthValue();
              if (otherTruthValue)
              {
                  // Truth values are stored differently for multi-chain
                if (numChains_ > 1)
                  truthValues_[otherIdx][chainIdx] = false;
                else
                  otherGndPred->setTruthValue(false);
              
                affectedGndPreds.clear();
                affectedGndPredIndices.clear();
                gndPredFlippedUpdates(otherIdx, chainIdx, affectedGndPreds,
                                      affectedGndPredIndices);
                updateWtsForGndPreds(affectedGndPreds, affectedGndPredIndices,
                                     chainIdx);
              }
            }
          }
            // Set truth value and update wts for chosen atom
            // Truth values are stored differently for multi-chain
          if (numChains_ > 1) truthValues_[idx][chainIdx] = true;
          else gndPred->setTruthValue(true);
          affectedGndPreds.clear();
          affectedGndPredIndices.clear();
          gndPredFlippedUpdates(idx, chainIdx, affectedGndPreds,
                                affectedGndPredIndices);
          updateWtsForGndPreds(affectedGndPreds, affectedGndPredIndices,
                               chainIdx);
        }
          // If in actual gibbs sampling phase, track the num of times
          // the ground predicate is set to true
        if (!burningIn) numTrue_[idx]++;
      }
    }

      // Now go through all preds not in blocks
    for (int i = 0; i < state_->getNumAtoms(); i++)
    {
        // Predicates in blocks have been handled above
      if (state_->getBlockIndex(i) >= 0) continue;

      if (mcmcdebug)
      {
        cout << "Chain " << chainIdx << ": Probability of pred "
             << i << " is " << getProbabilityOfPred(i, chainIdx, 1) << endl;
      }
      
      bool newAssignment
        = genTruthValueForProb(getProbabilityOfPred(i, chainIdx, 1));

        // Truth values are stored differently for multi-chain
      bool truthValue;
      GroundPredicate* gndPred = state_->getGndPred(i);
      if (numChains_ > 1) truthValue = truthValues_[i][chainIdx];
      else truthValue = gndPred->getTruthValue();
        // If gndPred is flipped, do updates & find all affected gndPreds
      if (newAssignment != truthValue)
      {
        if (mcmcdebug)
        {
          cout << "Chain " << chainIdx << ": Changing truth value of pred "
               << i << " to " << newAssignment << endl;
        }
        
        if (numChains_ > 1) truthValues_[i][chainIdx] = newAssignment;
        else gndPred->setTruthValue(newAssignment);
        affectedGndPreds.clear();
        affectedGndPredIndices.clear();
        gndPredFlippedUpdates(i, chainIdx, affectedGndPreds,
                              affectedGndPredIndices);
        updateWtsForGndPreds(affectedGndPreds, affectedGndPredIndices,
                             chainIdx);
      }

        // If in actual gibbs sampling phase, track the num of times
        // the ground predicate is set to true
      if (!burningIn && newAssignment) numTrue_[i]++;
    }
      // If keeping track of true clause groundings
    if (!burningIn && trackClauseTrueCnts_)
      state_->getNumClauseGndings(clauseTrueCnts_, true);

    if (mcmcdebug) cout << "End of Gibbs step" << endl;
  }

  /**
   * Updates the weights of affected ground predicates. These are the ground
   * predicates which are in clauses of predicates which have had their truth
   * value changed.
   * 
   * @param gndPreds Ground predicates whose weights should be updated.
   * @param chainIdx Index of chain where updating occurs.
   */
  void updateWtsForGndPredsH(GroundPredicateHashArray& gndPreds,
	  Array<int>& gndPredIndices, const int& chainIdx)
  {
	  if (mcmcdebug) cout << "Entering updateWtsForGndPreds" << endl;
	  // for each ground predicate whose MB has changed
	  for (int g = 0; g < gndPreds.size(); g++)
	  {
		  double wtIfNoChange = 0, wtIfInverted = 0, wt;
		  // Ground clauses in which this pred occurs
		  Array<int>& negGndClauses =
			  hstate_->getNegOccurenceArray(gndPredIndices[g] + 1);
		  Array<int>& posGndClauses =
			  hstate_->getPosOccurenceArray(gndPredIndices[g] + 1);

		  int gndClauseIdx;
		  bool sense;
		  if (mcmcdebug)
		  {
			  cout << "Ground clauses in which pred " << g << " occurs neg.: "
				  << negGndClauses.size() << endl;
			  cout << "Ground clauses in which pred " << g << " occurs pos.: "
				  << posGndClauses.size() << endl;
		  }

		  for (int i = 0; i < negGndClauses.size() + posGndClauses.size(); i++)
		  {
			  if (i < negGndClauses.size())
			  {
				  gndClauseIdx = negGndClauses[i];
				  if (mcmcdebug) cout << "Neg. in clause " << gndClauseIdx << endl;
				  sense = false;
			  }
			  else
			  {
				  gndClauseIdx = posGndClauses[i - negGndClauses.size()];
				  if (mcmcdebug) cout << "Pos. in clause " << gndClauseIdx << endl;
				  sense = true;
			  }

			  GroundClause* gndClause = hstate_->getGndClause(gndClauseIdx);
			  if (gndClause->isHardClause())
				  wt = hstate_->getClauseCost(gndClauseIdx);
			  else
				  wt = gndClause->getWt();
			  // NumTrueLits are stored differently for multi-chain
			  int numSatLiterals;
			  if (numChains_ > 1)
				  numSatLiterals = numTrueLits_[gndClauseIdx][chainIdx];
			  else
				  numSatLiterals = hstate_->getNumTrueLits(gndClauseIdx);
			  if (numSatLiterals > 1)
			  {
				  // Some other literal is making it sat, so it doesn't matter
				  // if pos. clause. If neg., nothing can be done to unsatisfy it.
				  if (wt > 0)
				  {
					  wtIfNoChange += wt;
					  wtIfInverted += wt;
				  }
			  }
			  else 
				  if (numSatLiterals == 1) 
				  {
					  if (wt > 0) wtIfNoChange += wt;
					  // Truth values are stored differently for multi-chain
					  bool truthValue;
					  if (numChains_ > 1)
						  truthValue = truthValues_[gndPredIndices[g]][chainIdx];
					  else
						  truthValue = gndPreds[g]->getTruthValue();
					  // If the current truth value is the same as its sense in gndClause
					  if (truthValue == sense) 
					  {
						  // This gndPred is the only one making this function satisfied
						  if (wt < 0) wtIfInverted += abs(wt);
					  }
					  else 
					  {
						  // Some other literal is making it satisfied
						  if (wt > 0) wtIfInverted += wt;
					  }
				  }
				  else
					  if (numSatLiterals == 0) 
					  {
						  // None satisfy, so when gndPred switch to its negative, it'll satisfy
						  if (wt > 0) wtIfInverted += wt;
						  else if (wt < 0) wtIfNoChange += abs(wt);
					  }
		  } // for each ground clause that gndPred appears in

		  if (mcmcdebug)
		  {
			  cout << "wtIfNoChange of pred " << g << ": "
				  << wtIfNoChange << endl;
			  cout << "wtIfInverted of pred " << g << ": "
				  << wtIfInverted << endl;
		  }

		  // Clause info is stored differently for multi-chain
		  if (numChains_ > 1)
		  {
			  if (truthValues_[gndPredIndices[g]][chainIdx]) 
			  {
				  wtsWhenTrue_[gndPredIndices[g]][chainIdx] = wtIfNoChange;
				  wtsWhenFalse_[gndPredIndices[g]][chainIdx] = wtIfInverted;
			  }
			  else 
			  {
				  wtsWhenFalse_[gndPredIndices[g]][chainIdx] = wtIfNoChange;
				  wtsWhenTrue_[gndPredIndices[g]][chainIdx] = wtIfInverted;
			  }
		  }
		  else
		  { // Single chain
			  if (gndPreds[g]->getTruthValue())
			  {
				  gndPreds[g]->setWtWhenTrue(wtIfNoChange);
				  gndPreds[g]->setWtWhenFalse(wtIfInverted);
			  }
			  else
			  {
				  gndPreds[g]->setWtWhenFalse(wtIfNoChange);
				  gndPreds[g]->setWtWhenTrue(wtIfInverted);            
			  }
		  }
	  } // for each ground predicate whose MB has changed
	  if (mcmcdebug) cout << "Leaving updateWtsForGndPreds" << endl;
  }

  


  /**
   * Updates the weights of affected ground predicates. These are the ground
   * predicates which are in clauses of predicates which have had their truth
   * value changed.
   * 
   * @param gndPreds Ground predicates whose weights should be updated.
   * @param chainIdx Index of chain where updating occurs.
   */
  void updateWtsForGndPreds(GroundPredicateHashArray& gndPreds,
                            Array<int>& gndPredIndices,
                            const int& chainIdx)
  {
    if (mcmcdebug) cout << "Entering updateWtsForGndPreds" << endl;
      // for each ground predicate whose MB has changed
    for (int g = 0; g < gndPreds.size(); g++)
    {
      double wtIfNoChange = 0, wtIfInverted = 0, wt;
        // Ground clauses in which this pred occurs
      Array<int>& negGndClauses =
        state_->getNegOccurenceArray(gndPredIndices[g] + 1);
      Array<int>& posGndClauses =
        state_->getPosOccurenceArray(gndPredIndices[g] + 1);
      int gndClauseIdx;
      bool sense;
      
      if (mcmcdebug)
      {
        cout << "Ground clauses in which pred " << g << " occurs neg.: "
             << negGndClauses.size() << endl;
        cout << "Ground clauses in which pred " << g << " occurs pos.: "
             << posGndClauses.size() << endl;
      }
      
      for (int i = 0; i < negGndClauses.size() + posGndClauses.size(); i++)
      {
        if (i < negGndClauses.size())
        {
          gndClauseIdx = negGndClauses[i];
          if (mcmcdebug) cout << "Neg. in clause " << gndClauseIdx << endl;
          sense = false;
        }
        else
        {
          gndClauseIdx = posGndClauses[i - negGndClauses.size()];
          if (mcmcdebug) cout << "Pos. in clause " << gndClauseIdx << endl;
          sense = true;
        }
        
        GroundClause* gndClause = state_->getGndClause(gndClauseIdx);
        if (gndClause->isHardClause())
          wt = state_->getClauseCost(gndClauseIdx);
        else
          wt = gndClause->getWt();
          // NumTrueLits are stored differently for multi-chain
        int numSatLiterals;
        if (numChains_ > 1)
          numSatLiterals = numTrueLits_[gndClauseIdx][chainIdx];
        else
          numSatLiterals = state_->getNumTrueLits(gndClauseIdx);

        if (numSatLiterals > 1)
        {
            // Some other literal is making it sat, so it doesn't matter
            // if pos. clause. If neg., nothing can be done to unsatisfy it.
          if (wt > 0)
          {
            wtIfNoChange += wt;
            wtIfInverted += wt;
          }
        }
        else 
        if (numSatLiterals == 1) 
        {
          if (wt > 0) wtIfNoChange += wt;
            // Truth values are stored differently for multi-chain
          bool truthValue;
          if (numChains_ > 1)
            truthValue = truthValues_[gndPredIndices[g]][chainIdx];
          else
            truthValue = gndPreds[g]->getTruthValue();
            // If the current truth value is the same as its sense in gndClause
          if (truthValue == sense) 
          {
            // This gndPred is the only one making this function satisfied
            if (wt < 0) wtIfInverted += abs(wt);
          }
          else 
          {
              // Some other literal is making it satisfied
            if (wt > 0) wtIfInverted += wt;
          }
        }
        else
        if (numSatLiterals == 0) 
        {
          // None satisfy, so when gndPred switch to its negative, it'll satisfy
          if (wt > 0) wtIfInverted += wt;
          else if (wt < 0) wtIfNoChange += abs(wt);
        }
      } // for each ground clause that gndPred appears in

      if (mcmcdebug)
      {
        cout << "wtIfNoChange of pred " << g << ": "
             << wtIfNoChange << endl;
        cout << "wtIfInverted of pred " << g << ": "
             << wtIfInverted << endl;
      }

        // Clause info is stored differently for multi-chain
      if (numChains_ > 1)
      {
        if (truthValues_[gndPredIndices[g]][chainIdx]) 
        {
          wtsWhenTrue_[gndPredIndices[g]][chainIdx] = wtIfNoChange;
          wtsWhenFalse_[gndPredIndices[g]][chainIdx] = wtIfInverted;
        }
        else 
        {
          wtsWhenFalse_[gndPredIndices[g]][chainIdx] = wtIfNoChange;
          wtsWhenTrue_[gndPredIndices[g]][chainIdx] = wtIfInverted;
        }
      }
      else
      { // Single chain
        if (gndPreds[g]->getTruthValue())
        {
          gndPreds[g]->setWtWhenTrue(wtIfNoChange);
          gndPreds[g]->setWtWhenFalse(wtIfInverted);
        }
        else
        {
          gndPreds[g]->setWtWhenFalse(wtIfNoChange);
          gndPreds[g]->setWtWhenTrue(wtIfInverted);            
        }
      }
    } // for each ground predicate whose MB has changed
    if (mcmcdebug) cout << "Leaving updateWtsForGndPreds" << endl;
  }

  /**
   * Chooses an atom from a block according to their probabilities.
   * 
   * @param chainIdx Index of chain.
   * @param block Block of predicate indices from which one is chosen.
   * @param invTemp InvTemp used in simulated tempering.
   * 
   * @return Index of chosen atom in the block.
   */
  int gibbsSampleFromBlock(const int& chainIdx, const int& blockIndex,
                           const double& invTemp)
  {
    Array<double> numerators;
    double denominator = 0;

    int blockSize = state_->getDomain()->getBlockSize(blockIndex);
    for (int i = 0; i < blockSize; i++)
    {
      const Predicate* pred =
        state_->getDomain()->getPredInBlock(i, blockIndex);
      GroundPredicate* gndPred = new GroundPredicate((Predicate*)pred);
      int idx = state_->getIndexOfGroundPredicate(gndPred);
      
      delete gndPred;
      delete pred;

        // Prob is 0 if atom not in state
      double prob = 0.0;
        // Pred is in the state; otherwise, prob is zero
      if (idx >= 0)
        prob = getProbabilityOfPred(idx, chainIdx, invTemp);

      numerators.append(prob);
      denominator += prob;
    }

    double r = random();
    double numSum = 0.0;
    for (int i = 0; i < blockSize; i++)
    {
      numSum += numerators[i];
      if (r < ((numSum / denominator) * RAND_MAX))
      {
        return i;
      }
    }
    return blockSize - 1;
  }

  /**
   * Updates information when a ground predicate is flipped and retrieves the
   * Markov blanket of the ground predicate.
   * 
   * @param gndPredIdx Index of ground pred which was flipped.
   * @param chainIdx Index of chain in which the flipping occured.
   * @param affectedGndPreds Holds the Markov blanket of the ground predicate.
   */
  void gndPredFlippedUpdates(const int& gndPredIdx, const int& chainIdx,
                             GroundPredicateHashArray& affectedGndPreds,
                             Array<int>& affectedGndPredIndices)
  {
    if (mcmcdebug) cout << "Entering gndPredFlippedUpdates" << endl;
    int numAtoms = state_->getNumAtoms();
    GroundPredicate* gndPred = state_->getGndPred(gndPredIdx);
    int idx = affectedGndPreds.append(gndPred, numAtoms);
    if (idx >= 0) affectedGndPredIndices.append(gndPredIdx);
    assert(affectedGndPreds.size() <= numAtoms);

    Array<int>& negGndClauses =
      state_->getNegOccurenceArray(gndPredIdx + 1);
    Array<int>& posGndClauses =
      state_->getPosOccurenceArray(gndPredIdx + 1);
    int gndClauseIdx;
    GroundClause* gndClause; 
    bool sense;

      // Find the Markov blanket of this ground predicate
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
      gndClause = state_->getGndClause(gndClauseIdx);

        // Different for multi-chain
      if (numChains_ > 1)
      {
        if (truthValues_[gndPredIdx][chainIdx] == sense)
          numTrueLits_[gndClauseIdx][chainIdx]++;
        else
          numTrueLits_[gndClauseIdx][chainIdx]--;
      }
      else
      { // Single chain
        if (gndPred->getTruthValue() == sense)
          state_->incrementNumTrueLits(gndClauseIdx);
        else
          state_->decrementNumTrueLits(gndClauseIdx);
      }
      
      for (int j = 0; j < gndClause->getNumGroundPredicates(); j++)
      {
        const GroundPredicateHashArray* gpha = state_->getGndPredHashArrayPtr();
        GroundPredicate* pred = 
          (GroundPredicate*)gndClause->getGroundPredicate(j,
            (GroundPredicateHashArray*)gpha);
        int idx = affectedGndPreds.append(pred, numAtoms);
        if (idx >= 0) affectedGndPredIndices.append(
                               abs(gndClause->getGroundPredicateIndex(j)) - 1);
        assert(affectedGndPreds.size() <= numAtoms);
      }
    }
    if (mcmcdebug) cout << "Leaving gndPredFlippedUpdates" << endl;
  }

  double getProbTrue(const int& predIdx) const { return numTrue_[predIdx]; }
  
  void setProbTrue(const int& predIdx, const double& p)
  { 
    assert(p >= 0);
    numTrue_[predIdx] = p;
  }

  /**
   * The atom assignment in the best state is saved to a chain in the
   * ground predicates.
   * 
   * @param chainIdx Index of chain to which the atom assigment is saved
   */
  void saveLowStateToChain(const int& chainIdx)
  {
    for (int i = 0; i < state_->getNumAtoms(); i++)
      truthValues_[i][chainIdx] = state_->getValueOfLowAtom(i + 1);
  }

  /**
   * Sets the user-set parameters for this MCMC algorithm.
   * 
   * @param params MCMC parameters for this algorithm.
   */
  void setMCMCParameters(MCMCParams* params)
  {
      // User-set parameters
    numChains_ = params->numChains;
    burnMinSteps_ = params->burnMinSteps;
    burnMaxSteps_ = params->burnMaxSteps;
    minSteps_ = params->minSteps;
    maxSteps_ = params->maxSteps;
    maxSeconds_ = params->maxSeconds;    
  }

  void scaleSamples(double factor)
  {
    minSteps_ = (int)(minSteps_ * factor);
    maxSteps_ = (int)(maxSteps_ * factor);
  }
  
 public: 
  // truth value assignment for continuous variables for each chain
  Array<Array<double> > truthValuesCont_;
  // storing values for the discrete parts of hybrid clauses for each chain
  Array<Array<bool> > hybridClauseDisPartValueMCMC_; 
  // storing values for the continuous parts of hybrid clauses for each chain
  Array<Array<double> > hybridClauseContPartValueMCMC_; 
  
 protected:
 
    ////////// BEGIN: User parameters ///////////
    // No. of chains which MCMC will use
  int numChains_;
    // Min. no. of burn-in steps MCMC will take per chain
  int burnMinSteps_;
    // Max. no. of burn-in steps MCMC will take per chain
  int burnMaxSteps_;
    // Min. no. of sampling steps MCMC will take per chain
  int minSteps_;
    // Max. no. of sampling steps MCMC will take per chain
  int maxSteps_;
    // Max. no. of seconds MCMC should run
  int maxSeconds_;
    ////////// END: User parameters ///////////

    // Truth values in each chain for each ground predicate (truthValues_[p][c])
  Array<Array<bool> > truthValues_;
    // Wts when false in each chain for each ground predicate
  Array<Array<double> > wtsWhenFalse_;
    // Wts when true in each chain for each groud predicate
  Array<Array<double> > wtsWhenTrue_;

    // Number of times each ground predicate is set to true
    // overloaded to hold probability that ground predicate is true
  Array<double> numTrue_; // numTrue_[p]

    // Num. of satisfying literals in each chain for each groud predicate
    // numTrueLits_[clause][chain]
  Array<Array<int> > numTrueLits_;
};

#endif /*MCMC_H_*/
