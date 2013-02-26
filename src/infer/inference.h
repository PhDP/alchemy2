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
#ifndef INFERENCE_H_
#define INFERENCE_H_

#include "variablestate.h"
#include "hvariablestate.h"

  // Default seed when none specified (some random number)
const long int DEFAULT_SEED = 2350877;

/**
 * Abstract class from which all inference algorithms are derived.
 * At least one function is pure virtual making this an abstract class
 * (it can not be instantiated).
 */
class Inference
{
 public:

  /**
   * Constructor: Every inference algorithm is required to have a VariableState
   * representing the state of variables and clauses and a seed for any
   * randomization in the algorithm. If there is no randomization, seed is not
   * used.
   * 
   * @param state State of the variables and clauses of the inference.
   * @param seed Seed used to initialize randomization in the algorithm.
   * @param trackClauseTrueCnts Indicates if true counts for each first-order
   * clause are being kept
   */
  Inference(VariableState* state, long int seed,
            const bool& trackClauseTrueCnts,
            Array<Array<Predicate* >* >* queryFormulas = NULL)
      : seed_(seed), state_(state), saveAllCounts_(false),
        clauseTrueCnts_(NULL), clauseTrueSqCnts_(NULL),
        numSamples_(0),
        allClauseTrueCnts_(NULL), oldClauseTrueCnts_(NULL),
        oldAllClauseTrueCnts_(NULL), queryFormulas_(queryFormulas)
  {
      // If seed not specified, then init always to same random number
    if (seed_ == -1) seed_ = DEFAULT_SEED;
    srandom(seed_);

    trackClauseTrueCnts_ = trackClauseTrueCnts;
    if (trackClauseTrueCnts_ && state_)
    {
      int numClauses = state_->getMLN()->getNumClauses();

        // clauseTrueCnts_ and clauseTrueSqCnts_ will hold the true 
        // counts (and squared true counts) for each first-order clause
      clauseTrueCnts_ = new Array<double>(numClauses, 0);
      clauseTrueSqCnts_ = new Array<double>(numClauses, 0);
    }
    
    if (queryFormulas_)
      qfProbs_ = new Array<double>(queryFormulas_->size(), 0);
    else
      qfProbs_ = NULL;
  }
  
  Inference(HVariableState* state, long int seed,
            const bool& trackClauseTrueCnts)
      : seed_(seed), hstate_(state), saveAllCounts_(false),
        clauseTrueCnts_(NULL), clauseTrueCntsCont_(NULL), 
        clauseTrueSqCnts_(NULL), numSamples_(0),
        allClauseTrueCnts_(NULL), oldClauseTrueCnts_(NULL),
        oldAllClauseTrueCnts_(NULL), queryFormulas_(NULL),
        qfProbs_(NULL)
  {
        // If seed not specified, then init always to same random number
      if (seed_ == -1) seed_ = DEFAULT_SEED;
      srandom(seed_);
	  
	  trackClauseTrueCnts_ = trackClauseTrueCnts;
	  if (trackClauseTrueCnts_ && hstate_)
	  {
		  // clauseTrueCnts_ will hold the true counts for each first-order
          // clause
		  clauseTrueCnts_ = new Array<double>;
		  clauseTrueCnts_->growToSize(hstate_->getMLN()->getNumClauses(), 0);
		  clauseTrueCntsCont_ = new Array<double>;
		  clauseTrueCntsCont_->growToSize(hstate_->getNumContFormulas(), 0);
	  }
  }
  
  /**
   * Virtual destructor.
   */
  virtual ~Inference()
  {
    delete clauseTrueCnts_;
    delete clauseTrueSqCnts_;
    delete allClauseTrueCnts_;

    delete oldAllClauseTrueCnts_;
    delete oldClauseTrueCnts_;
    
    if (qfProbs_) delete qfProbs_;
  }


  void saveAllCounts(bool saveCounts=true)
  {
    if (saveAllCounts_ == saveCounts)
      return;

    saveAllCounts_ = saveCounts;
    if (saveCounts)
    {
      allClauseTrueCnts_ = new Array<Array<double> >;
      oldAllClauseTrueCnts_ = new Array<Array<double> >;
    }
    else
    {
      delete allClauseTrueCnts_;
      delete oldAllClauseTrueCnts_;
      allClauseTrueCnts_ = NULL;
      oldAllClauseTrueCnts_ = NULL;
    }
  }


  /**
   * Initializes the inference algorithm.
   */
  virtual void init() = 0;

  /**
   * Performs the inference algorithm.
   */
  virtual void infer() = 0;

  /**
   * Prints out the network. Currently only available in BP (uses factor graph).
   */
  virtual void printNetwork(ostream& out) = 0; 

  /**
   * Prints the probabilities of each predicate to a stream.
   */
  virtual void printProbabilities(ostream& out) = 0;
  
  /**
   * Puts the predicates whose probability has changed (more than probDelta in
   * the case of prob. inference) with respect to the reference vector oldProbs
   * in string form and the corresponding probabilities of each predicate in two
   * vectors.
   */
  virtual void getChangedPreds(vector<string>& changedPreds,
                               vector<float>& probs,
                               vector<float>& oldProbs,
                               const float& probDelta) = 0;

  
  /**
   * Prints the predicates inferred to be true to a stream.
   */
  virtual void printTruePreds(ostream& out) = 0;
  virtual void printTruePredsH(ostream& out) = 0;

  /**
   * Gets the probability of a ground predicate.
   */
  virtual double getProbability(GroundPredicate* const& gndPred) = 0;
  virtual double getProbabilityH(GroundPredicate* const& gndPred) = 0;

  /**
   * Print probabilities of the query formulas.
   */
  void printQFProbs(ostream& out, Domain* domain)
  {
    if (qfProbs_)
    {
      for (int i = 0; i < queryFormulas_->size(); i++)
      {
        Array<Predicate* >* formula = (*queryFormulas_)[i];
        for (int j = 0; j < formula->size(); j++)
        {
          (*formula)[j]->printWithStrVar(out, domain);
          if (j != formula->size() - 1) out << " ^ ";
        }
        out << " " << (*qfProbs_)[i] << endl;
      }
    }
  }

  long int getSeed() { return seed_; }
  void setSeed(long int s) { seed_ = s; }
  
  VariableState* getState() { return state_; }
  void setState(VariableState* s) { state_ = s; }

  HVariableState* getHState() { return hstate_; }
  void setHState(HVariableState* s) { hstate_ = s; }
  
  /**
   * Increase or decrease the number of MCMC samples by a multiplicative
   * factor.  (For MaxWalkSAT, this could perhaps change the number
   * of flips or something that trades off speed and accuracy.)
   */
  virtual void scaleSamples(double factor) { /* Override this... */ }


  
  const Array<double>* getClauseTrueCnts()   { return clauseTrueCnts_; }
  const Array<double>* getClauseTrueSqCnts() { return clauseTrueSqCnts_; }
  int getNumSamples() const { return numSamples_; }

  // Compute the full Hessian matrix from the stored inferred counts.
  // Only works when saveAllCounts(true) has been called before inference.
  // Caller is responsible for deleting Hessian.
  //
  // WARNING: The size of the Hessian is equal to the number of weights
  // squared.  Therefore, this is not practical in a model with thousands
  // of weights!
  const Array<Array<double> >* getHessian()
  {
    int numClauses = state_->getMLN()->getNumClauses();
    int numSamples = allClauseTrueCnts_->size();

    // Allocate Hessian
    Array<Array<double> >* hessian = new Array<Array<double> >(numClauses);
    for (int i = 0; i < numClauses; i++)
      (*hessian)[i].growToSize(numClauses);

    // The i jth element is:
    // E[n_i] E[n_j] - E[n_i * n_j]
    // where n is the vector of all clause counts

    for (int i = 0; i < numClauses; i++)
    {
      for (int j = 0; j < numClauses; j++)
      {
        double ni = 0.0;
        double nj = 0.0;
        double ninj = 0.0;
        for (int s = 0; s < numSamples; s++)
        {
          ni += (*allClauseTrueCnts_)[s][i];
          nj += (*allClauseTrueCnts_)[s][j];
          ninj += (*allClauseTrueCnts_)[s][i] 
                * (*allClauseTrueCnts_)[s][j]; 
        }
        double n = numSamples;
        (*hessian)[i][j] = ni/n * nj/n - ninj/n;
      }
    }

    return hessian;
  }


  // Alternate way to compute product of Hessian with vector
  // WARNING: much less efficient!!
  const Array<double>* getHessianVectorProduct2(Array<double>& v)
  {
    int numClauses = state_->getMLN()->getNumClauses();
    const Array<Array<double> >* hessian = getHessian();
    Array<double>* product = new Array<double>(numClauses,0);

    for (int clauseno = 0; clauseno < numClauses; clauseno++)
    {
      (*product)[clauseno] = 0.0;
      for (int i = 0; i < numClauses; i++)
        (*product)[clauseno] += (*hessian)[clauseno][i] * v[i];
    }

    delete hessian;
    return product;
  }


  const Array<double>* getHessianVectorProduct(const Array<double>& v)
  {
    int numClauses = state_->getMLN()->getNumClauses();
    int numSamples = allClauseTrueCnts_->size();

      // For minimizing the negative log likelihood, 
      // the ith element of H v is:
      //   E[n_i * vn] - E[n_i] E[vn]
      // where n is the vector of all clause counts
      // and vn is the dot product of v and n.

    double sumVN = 0;
    Array<double> sumN(numClauses, 0);
    Array<double> sumNiVN(numClauses, 0);

    // Get sufficient statistics from each sample, 
    // so we can compute expectations
    for (int s = 0; s < numSamples; s++) 
    {
      Array<double>& n = (*allClauseTrueCnts_)[s];

      // Compute v * n
      double vn = 0;
      for (int i = 0; i < numClauses; i++)
        vn += v[i] * n[i];
      
      // Tally v*n, n_i, and n_i v*n
      sumVN += vn;
      for (int i = 0; i < numClauses; i++)
      {
        sumN[i]    += n[i];
        sumNiVN[i] += n[i] * vn;
      }
    }

    // Compute actual product from the sufficient stats
    Array<double>* product = new Array<double>(numClauses,0);
    for (int clauseno = 0; clauseno < numClauses; clauseno++)
    {
      double E_vn = sumVN/numSamples;
      double E_ni = sumN[clauseno]/numSamples;
      double E_nivn = sumNiVN[clauseno]/numSamples;
      (*product)[clauseno] = E_nivn - E_ni * E_vn;
    }

    return product;
  }


  void resetCnts() 
  {
    for (int clauseno = 0; clauseno < clauseTrueCnts_->size(); clauseno++)
    {
      (*clauseTrueCnts_)[clauseno]   = 0;
      (*clauseTrueSqCnts_)[clauseno] = 0;
    }
    numSamples_ = 0;

    if (saveAllCounts_)
    {
      delete allClauseTrueCnts_;
      allClauseTrueCnts_ = new Array<Array<double> >;
    }
  }


  void saveCnts()
  {
    if (!saveAllCounts_)
      return;

    // DEBUG
    //cout << "Saving counts.  numSamples_ = " << numSamples_ << endl;

    delete oldAllClauseTrueCnts_;
    oldAllClauseTrueCnts_ = new Array<Array<double> > (*allClauseTrueCnts_);

    /* DEBUG
    cout << "old counts size: " << oldAllClauseTrueCnts_->size() << endl;
    cout << "marker1: " << (*allClauseTrueCnts_)[0][0] << endl;
    cout << "marker2: " << (*clauseTrueCnts_)[0] << endl;
    */
  }


  void restoreCnts()
  {
    if (!saveAllCounts_)
      return;

    resetCnts();

    *allClauseTrueCnts_ = *oldAllClauseTrueCnts_;
    for (int i = 0; i < allClauseTrueCnts_->size(); i++) 
    {
      int numcounts = (*allClauseTrueCnts_)[i].size();
      for (int j = 0; j < numcounts; j++)
      {
        double currcount = (*allClauseTrueCnts_)[i][j];
        (*clauseTrueCnts_)[j]   += currcount;
        (*clauseTrueSqCnts_)[j] += currcount * currcount;
      }
      numSamples_++;
    }

    /* DEBUG
    cout << "marker1: " << (*allClauseTrueCnts_)[0][0] << endl;
    cout << "marker2: " << (*clauseTrueCnts_)[0] << endl;
    cout << "numSamples_ = " << numSamples_ << endl;
    */
  }


  void tallyCntsFromState()
  {
    int numcounts = clauseTrueCnts_->size();
    Array<double> currCounts(numcounts, 0.0);
    state_->getNumClauseGndings(&currCounts, true);

    if (saveAllCounts_)
    {
      allClauseTrueCnts_->append(Array<double>());
      (*allClauseTrueCnts_)[numSamples_].growToSize(numcounts);
    }

    for (int i = 0; i < numcounts; i++)
    {
      if (saveAllCounts_)
        (*allClauseTrueCnts_)[numSamples_][i] = currCounts[i];
      
      (*clauseTrueCnts_)[i]   += currCounts[i];
      (*clauseTrueSqCnts_)[i] += currCounts[i] * currCounts[i];
    }
    numSamples_++;
  }

 protected:
  
    // Seed for randomizer. If not set, then date + time is used
  long int seed_;
    // State of atoms and clauses used during inference
    // Does not belong to inference (do not delete)
  VariableState* state_;
  HVariableState* hstate_;
  
    // Save all counts for all iterations
  bool saveAllCounts_;
    // Holds the average number of true groundings of each first-order clause
    // in the mln associated with this inference
  Array<double>* clauseTrueCnts_;
  Array<double>* clauseTrueCntsCont_;
  
    // Holds the average number of true groundings squared of each 
    // first-order clause in the mln associated with this inference
  Array<double>* clauseTrueSqCnts_;
    // Number of samples taken of the true counts
  int numSamples_;
    // Indicates if true counts for each first-order clause are being kept
  bool trackClauseTrueCnts_;
    // Where these counts are stored: (*allClauseTrueCnts_)[i][j] has the true
    // counts for clause j in sample i
  Array<Array<double> >* allClauseTrueCnts_;

  Array<double>* oldClauseTrueCnts_;
  Array<Array<double> >* oldAllClauseTrueCnts_;
  
  Array<Array<Predicate* >* >* queryFormulas_;
  Array<double>* qfProbs_;
};

#endif /*INFERENCE_H_*/
