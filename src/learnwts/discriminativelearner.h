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
#ifndef DISCRIMINATIVE_LEARNER_H_OCT_30_2005
#define DISCRIMINATIVE_LEARNER_H_OCT_30_2005

#include "infer.h"
#include "clause.h"
#include "timer.h"
#include "indextranslator.h"
#include "maxwalksat.h"
#include "bp.h"

const bool dldebug = false;
const double EPSILON = .00001;

#ifndef PI
#define PI 3.14159265
#endif

double dotprod(const double* v1, const double* v2, int dim)
{
  double total = 0.0;
  for (int i = 0; i < dim; i++)
    total += v1[i] * v2[i];
  return total;
}

double dotprod(const double* v1, const Array<double>& v2, int dim)
{
  double total = 0.0;
  for (int i = 0; i < dim; i++)
    total += v1[i] * v2[i];
  return total;
}

double dotprod(const Array<double>& v1, const Array<double>& v2, int dim)
{
  double total = 0.0;
  for (int i = 0; i < dim; i++)
  {
    total += v1[i] * v2[i];
  }
  return total;
}

double vlen(const double* v1, int dim)
{ return sqrt(dotprod(v1,v1,dim)); }
double vlen(const Array<double>& v1, int dim)
{ return sqrt(dotprod(v1,v1,dim)); }

bool isOutputIter(int iter)
{
  while (iter % 10 == 0)
    iter /= 10;

  if (iter == 1 || iter == 2 || iter == 5)
    return true;
  else
    return false;
}



/**
 * Discriminative learning algorithms (see "Discriminative Training of Markov
 * Logic Networks", Singla and Domingos, 2005 and "Efficient Weight Learning for
 * Markov Logic Networks", Lowd and Domingos, 2007).
 */
class DiscriminativeLearner
{
 public:

  /**
   * Constructor. Various variables are initialized, relevant clauses are
   * determined and weights and inference procedures are initialized.
   *
   * @param inferences Array of inference procedures to be used for inference
   * in each domain.
   * @param nonEvidPredNames Names of non-evidence predicates. This is used to
   * determine the relevant clauses.
   * @param idxTrans IndexTranslator needed when multiple dbs are used and they
   * don't line up.
   * @param lazyInference If true, lazy inference is used.
   * @param withEM If true, EM is used to fill in missing values.
   * @param rescaleGradient If true, use per-weight learning rates
   * @param method Determines how direction and step size are chosen
   * @param lambda Initial value of lambda for SMD or CG
   * @param preconditionCG Whether or not to use a preconditioner with
   * scaled conjugate gradient
   * @param maxLambda Maximum value of lambda for CG
   * @param min_ll_change Minimum change in likelihood neede to continue
   * iterations
   */
  DiscriminativeLearner(const Array<Inference*>& inferences,
                        const StringHashArray& nonEvidPredNames,
                        IndexTranslator* const & idxTrans,
                        const bool& lazyInference, const bool& withEM,
                        const Array<Inference*>& emInferences,
                        const bool& rescaleGradient, const int& method,
                        const double& lambda, const bool& preconditionCG,
                        const double& maxLambda, const double& min_ll_change)
    : domainCnt_(inferences.size()), idxTrans_(idxTrans),
      lazyInference_(lazyInference), rescaleGradient_(rescaleGradient),
      method_(method),
      // HACK: for now, we'll use the SMD lambda value for CG, even
      // though the two represent *very* different things!
      cg_lambda_(lambda), preconditionCG_(preconditionCG),
      maxBacktracks_(1000), backtrackCount_(0),
      cg_max_lambda_(maxLambda), min_ll_change_(min_ll_change), withEM_(withEM)
  {
    cout << endl << "Constructing discriminative learner..." << endl << endl;

    inferences_.append(inferences);
    emInferences_.append(emInferences);
    logOddsPerDomain_.growToSize(domainCnt_);
    clauseCntPerDomain_.growToSize(domainCnt_);

    for (int i = 0; i < domainCnt_; i++)
    {
        // Check if using BP
      if ((bp_ = dynamic_cast<BP*>(inferences_[i])))
      {
        clauseCntPerDomain_[i] =
          bp_->getFactorGraph()->getMLN()->getNumClauses();
      }
      else
      {
        clauseCntPerDomain_[i] =
          inferences_[i]->getState()->getMLN()->getNumClauses();
      }
      logOddsPerDomain_[i].growToSize(clauseCntPerDomain_[i], 0);
    }

    totalTrueCnts_.growToSize(domainCnt_);
    totalFalseCnts_.growToSize(domainCnt_);
    defaultTrueCnts_.growToSize(domainCnt_);
    defaultFalseCnts_.growToSize(domainCnt_);
    relevantClausesPerDomain_.growToSize(domainCnt_);
    //relevantClausesFormulas_ is set in findRelevantClausesFormulas()

    findRelevantClauses(nonEvidPredNames);
    findRelevantClausesFormulas();

      // Initialize the clause wts
    initializeWts(nonEvidPredNames);

      // Initialize the inference / state
    for (int i = 0; i < inferences_.size(); i++)
    {
        // Check if using BP
      if ((bp_ = dynamic_cast<BP*>(inferences_[i])))
      {
        continue;
      }
      else
      {
        inferences_[i]->init();
      }
    }

    if (withEM_)
    {
        // Initialize the inference / state
      for (int i = 0; i < emInferences_.size(); i++)
      {
        VariableState* state = inferences_[i]->getState();
        const Domain* domain = state->getDomain();
        const GroundPredicateHashArray* knePreds = state->getKnePreds();
        const Array<TruthValue>* knePredValues = state->getKnePredValues();

          // Mark known non-evidence preds as evidence
        domain->getDB()->setValuesToGivenValues(knePreds, knePredValues);

        emInferences_[i]->init();

          // Set non-evidence values to unknown
        Array<TruthValue> tmpValues;
        tmpValues.growToSize(knePreds->size());
        domain->getDB()->setValuesToUnknown(knePreds, &tmpValues);
      }
    }
  }

  ~DiscriminativeLearner()
  {
    for (int i = 0; i < trainTrueCnts_.size(); i++)
    {
      delete[] trainTrueCnts_[i];
      delete[] trainFalseCnts_[i];
    }
    if (withEM_)
    {
      for (int i = 0; i < emDiffCnts_.size(); i++)
      {
        delete[] emDiffCnts_[i];
      }
    }
  }

    // set the prior means and std devs.
  void setMeansStdDevs(const int& arrSize, const double* const & priorMeans,
                       const double* const & priorStdDevs)
  {
    if (arrSize < 0)
    {
      usePrior_ = false;
      priorMeans_ = NULL;
      priorStdDevs_ = NULL;
    }
    else
    {
      //cout << "arr size = " << arrSize<<", clause count = "<<clauseCnt_<<endl;
      usePrior_ = true;
      priorMeans_ = priorMeans;
      priorStdDevs_ = priorStdDevs;

      //cout << "\t\t Mean \t\t Std Deviation" << endl;
      //for (int i = 0; i < arrSize; i++)
      //  cout << i << "\t\t" << priorMeans_[i]<<"\t\t"<<priorStdDevs_[i]<<endl;
    }
  }


  void setMLNWeights(double* const& weights)
  {
      // If there is one db or the clauses for multiple databases line up
    if (idxTrans_ == NULL)
    {
      int clauseCnt = clauseCntPerDomain_[0];
      for (int i = 0; i < domainCnt_; i++)
      {
        assert(clauseCntPerDomain_[i] == clauseCnt);
        const MLN* mln;

          // Check if using BP
        if ((bp_ = dynamic_cast<BP*>(inferences_[i])))
        {
          mln = bp_->getFactorGraph()->getMLN();
        }
        else
        {
          mln = inferences_[i]->getState()->getMLN();
        }

        // Clauses are typically shared among MLNs for different
        // domains, so many of these calls to setWt() are redundant.
        // However, it doesn't hurt to set the same weight twice.
        for (int j = 0; j < clauseCnt; j++)
        {
          Clause* c = (Clause*) mln->getClause(j);
          c->setWt(weights[j]);
        }
      }
    }
    else
    {   // The clauses for multiple databases do not line up
      Array<Array<double> >* wtsPerDomain = idxTrans_->getWtsPerDomain();
      const Array<Array<Array<IdxDiv>*> >* cIdxToCFIdxsPerDomain
        = idxTrans_->getClauseIdxToClauseFormulaIdxsPerDomain();

      for (int i = 0; i < domainCnt_; i++)
      {
        Array<double>& wts = (*wtsPerDomain)[i];
        memset((double*)wts.getItems(), 0, wts.size()*sizeof(double));

          //map clause/formula weights to clause weights
        for (int j = 0; j < wts.size(); j++)
        {
          Array<IdxDiv>* idxDivs = (*cIdxToCFIdxsPerDomain)[i][j];
          for (int k = 0; k < idxDivs->size(); k++)
            wts[j] += weights[ (*idxDivs)[k].idx ] / (*idxDivs)[k].div;
        }
      }

      for (int i = 0; i < domainCnt_; i++)
      {
        Array<bool>& relevantClauses = relevantClausesPerDomain_[i];
        int clauseCnt = clauseCntPerDomain_[i];
        Array<double>& wts = (*wtsPerDomain)[i];
        assert(wts.size() == clauseCnt);
        const MLN* mln;

          // Check if using BP
        if ((bp_ = dynamic_cast<BP*>(inferences_[i])))
        {
          mln = bp_->getFactorGraph()->getMLN();
        }
        else
        {
          mln = inferences_[i]->getState()->getMLN();
        }

        for (int j = 0; j < clauseCnt; j++)
        {
          Clause* c = (Clause*) mln->getClause(j);
          if (relevantClauses[j]) c->setWt(wts[j]);
          else                    c->setWt(0);
        }
      }
    }
  }


  void setLogOddsWeights(double* weights, int numWeights)
  {
      // If there is one db or the clauses for multiple databases line up
    if (idxTrans_ == NULL)
    {
      for (int i = 0; i < domainCnt_; i++)
      {
        Array<double>& logOdds = logOddsPerDomain_[i];
        assert(numWeights == logOdds.size());
        for (int j = 0; j < logOdds.size(); j++) weights[j] += logOdds[j];
      }
    }
    else
    { //the clauses for multiple databases do not line up
      const Array<Array<Array<IdxDiv>*> >* cIdxToCFIdxsPerDomain
        = idxTrans_->getClauseIdxToClauseFormulaIdxsPerDomain();

      Array<int> numLogOdds;
      Array<double> wtsForDomain;
      numLogOdds.growToSize(numWeights);
      wtsForDomain.growToSize(numWeights);

      for (int i = 0; i < domainCnt_; i++)
      {
        memset((int*)numLogOdds.getItems(), 0, numLogOdds.size()*sizeof(int));
        memset((double*)wtsForDomain.getItems(), 0,
               wtsForDomain.size()*sizeof(double));

        Array<double>& logOdds = logOddsPerDomain_[i];

          // Map the each log odds of a clause to the weight of a
          // clause/formula
        for (int j = 0; j < logOdds.size(); j++)
        {
          Array<IdxDiv>* idxDivs =(*cIdxToCFIdxsPerDomain)[i][j];
          for (int k = 0; k < idxDivs->size(); k++)
          {
            wtsForDomain[ (*idxDivs)[k].idx ] += logOdds[j];
            numLogOdds[ (*idxDivs)[k].idx ]++;
          }
        }

        for (int j = 0; j < numWeights; j++)
          if (numLogOdds[j] > 0) weights[j] += wtsForDomain[j]/numLogOdds[j];
      }
    }

    for (int i = 0; i < numWeights; i++)
      weights[i] /= domainCnt_;
  }


  void setRandGaussWeights(double* weights, int numWeights, const double* means,
                           const double* stdevs)
  {
    if (dldebug) cout << "Setting random weights: " << endl;
    for (int i = 0; i < numWeights; i++)
    {
      weights[i] = ExtRandom::gaussRandom(means[i], stdevs[i]);
      cout << weights[i] << endl;
    }
  }

  void setRandomWeights(double* weights, int numWeights, const double& min,
                        const double& max)
  {
    if (dldebug) cout << "Setting random weights: " << endl;
    for (int i = 0; i < numWeights; i++)
    {
      weights[i] = (max - min)*(double)rand()/(double)RAND_MAX + min;
      cout << weights[i] << endl;
    }
  }


  void getVariance(Array<double>& variance, int numWeights)
  {
    if (dldebug)
      cout << "Variances:" << endl;

      // Compute variance and gradient for each clause
    for (int clauseno = 0; clauseno < numWeights; clauseno++)
    {
      variance[clauseno] = 0.0;

        // Takes into account the prior
      if (usePrior_)
      {
         double sd = priorStdDevs_[clauseno];
         variance[clauseno] = 1.0/(sd * sd);
      }

        // Sum variance over all domains
      for (int i = 0; i < domainCnt_; i++)
      {
        const Array<double>* trueCnts = inferences_[i]->getClauseTrueCnts();
        const Array<double>* trueSqCnts = inferences_[i]->getClauseTrueSqCnts();
        int totalSamples = inferences_[i]->getNumSamples();

        double x   = (*trueCnts)[clauseno];
        double xsq = (*trueSqCnts)[clauseno];
        double n   = totalSamples;

          // Add variance for this domain
        variance[clauseno] += xsq/n - (x/n)*(x/n);
      }

      if (dldebug)
        cout << clauseno << ": " << variance[clauseno] << endl;
    }
  }


  const Array<double>* getHessianVectorProduct(const Array<double>& d)
  {
      // Compute Hessian vector product using counts from all inferences
    const Array<double>* Hd = inferences_[0]->getHessianVectorProduct(d);
    for (int i = 1; i < domainCnt_; i++)
    {
      const Array<double>* Hd_i = inferences_[i]->getHessianVectorProduct(d);
      for (int j = 0; j < Hd->size(); j++)
        (*Hd)[j] += (*Hd_i)[j];
      delete Hd_i;
    }

    return Hd;
  }


  double computeQuadraticStepLength(double* gradient, const Array<double>& d,
          const Array<double>* Hd, double lambda)
  {
    int numWeights = d.size();

      // Compute step length using trust region approach
    double dHd = dotprod(d, (*Hd), numWeights);
    double dd = dotprod(d, d, numWeights);
    double dg = dotprod(gradient, d, numWeights);
    double alpha = -dg/(dHd + lambda * dd);

      // Debug messages; can be helpful for diagnosing what is going on
    if (dldebug)
    {
      cout << "dHd = " << dHd << endl;
      cout << "dd = " << dd << endl;
      cout << "dg = " << dg << endl;
      cout << "alpha = " << alpha << endl;
    }

      // Because the problem is convex, the Hessian should always
      // be positive definite, and so alpha should always be non-negative.
      // Since our Hessian is approximate, we may end up with a negative
      // alpha.  In these cases, we used to make alpha zero (i.e., take a
      // step of size zero), but now we return the negative alpha and
      // let the caller deal with it.
    if (alpha < 0.0)
    {
        //alpha = 0.0;
      if (dldebug)
      {
        cout << "Alpha < 0!  Bad direction or Hessian." << endl;
        //cout << "Flipped alpha!  Setting step length to zero." << endl;
      }
    }

    return alpha;
  }


    // learn the weights
  void learnWeights(double* const & weights, const int& numWeights,
                    const int& maxIter, const double& maxSec,
                    const double& learningRate, const double& momentum,
                    bool initWithLogOdds, const int& mwsMaxSubsequentSteps,
                    bool periodicMLNs)
  {
    Timer timer;
    double begSec = timer.time();

    //cout << "Learning weights discriminatively... " << endl;
    memset(weights, 0, numWeights*sizeof(double));

    double* averageWeights = new double[numWeights];
    double* gradient = new double[numWeights];

      // For predicting and assessing changes to the gradient
    const Array<double>* delta_pred = NULL;
    const Array<double>* Hd = NULL;

      // For all methods
    double alpha = 1;
    Array<double> d(numWeights, 0.0);
    double oldgradient[numWeights];
    Array<double> oldd(numWeights, 0.0);

    //double pre_backtrack_lambda = cg_lambda_;

      // Set the initial weight to the average log odds across domains/databases
    if (initWithLogOdds)
      setLogOddsWeights(weights, numWeights);
    else if (withEM_)
    {
        // If using EM, init weights according to prior, or randomly
      if (usePrior_)
        setRandGaussWeights(weights, numWeights, priorMeans_, priorStdDevs_);
      else
        setRandomWeights(weights, numWeights, -1.0, 1.0);
    }

    // hard clause
    const MLN* mln = inferences_[0]->getState()->getMLN();
	for (int i=0; i<numWeights; i++) {
		const Clause* clause=mln->getClause(i);
		if (clause->isHardClause()) {
			weights[i]=10;
		}
	}

	  // Initialize averageWeights
	memcpy(averageWeights, weights, numWeights*sizeof(double));


      // Make sure to save all counts, for computing gradient
    for (int i = 0; i < domainCnt_; i++)
      inferences_[i]->saveAllCounts(true);

      // For DN method
    Array<double> currVariance(numWeights, 0);

      // Used in CG
    Array<double> oldweights(numWeights, 0.0);

      // True if we're backtracking -- prevents us from making a forward
      // step if we've already committed to reverting our weights.
    bool backtracked = false;

    //int itersSinceSignificant = 0;

    for (int iter = 1; iter <= maxIter; iter++)
    {
      cout << endl << "Iteration " << iter << " : " << endl << endl;

      double totalsec = timer.time() - begSec;
      int hour = (int)totalsec/3600;
      int min = (int)(totalsec - 3600*hour)/60;
      int sec = (int)(totalsec - 3600*hour - 60*min);

        // Print time so far
      cout << "Elapsed time: " << totalsec << "s";
      if (hour > 0)
        cout << " (" << hour << "h " << min << "m " << sec << "s)";
      else if (min > 0)
        cout << " (" << min << "m " << sec << "s)";
      cout << endl;

      if (maxSec > 0 && totalsec > maxSec)
      {
        cout << "Time limit exceeded.  Stopping learning." << endl;
        break;
      }

        // In 3rd iteration, we want to tell MWS to perform subsequentSteps
        // (Iter. 1 is random assigment if initial weights are 0)
      if (iter == 3)
      {
        for (int i = 0; i < inferences_.size(); i++)
        {
            // Check if using MWS
          if (MaxWalkSat* mws = dynamic_cast<MaxWalkSat*>(inferences_[i]))
          {
            mws->setMaxSteps(mwsMaxSubsequentSteps);
          }
        }
      }

#if 0
      int tcounts[numWeights];
      int tcountsTotal[numWeights];
      int icounts[numWeights];
      int icountsTotal[numWeights];
#endif

      // Set the weights and run inference
      setMLNWeights(weights);

      if (withEM_) fillInMissingValues();
      cout << "Running inference ..." << endl;
      infer();
      cout << "Done with inference" << endl;

      cout << "Getting gradient..." << endl;
      //getGradient(weights, gradient, numWeights,
      //        tcounts, tcountsTotal, icounts, icountsTotal);
      getGradient(weights, gradient, numWeights);

      double realdist = 1.0;
      double preddist = 1.0;
      if (iter > 1 && delta_pred != NULL && !backtracked)
      {
        Array<double> dist(numWeights, 0.0);
        for (int i = 0; i < numWeights; i++)
          dist[i] = weights[i] - oldweights[i];

        // Below are four different methods for estimating
        // the real and predicted change in the function value.
        // These values are used to size of the trust region
        // by updating lambda.
#if 1
        // Predicted change is quadratic approximation
        Array<double> avgPred(numWeights);
        avgPred.growToSize(numWeights);
        for (int i = 0; i < numWeights; i++)
          avgPred[i] = oldgradient[i] + (*delta_pred)[i]/2.0;
        preddist = dotprod(avgPred, dist, numWeights);

        // Real change is lower bound on actual change
        realdist = dotprod(gradient, dist, numWeights);
#elif 0
        // Actual distance is a lower bound on the change
        // in our function: new gradient * dist
        // Predicted change adopts the same form, but uses
        // the predicted gradient instead: pred gradient * dist
        Array<double> predgradient(numWeights, 0.0);
        for (int i = 0; i < numWeights; i++)
          predgradient[i] = oldgradient[i] + (*delta_pred)[i];

        preddist = dotprod(predgradient, dist, numWeights);
        realdist = dotprod(gradient, dist, numWeights);
#elif 0
        // Act as if we're optimizing the dot product of the
        // gradient with the search direction instead of
        // the function value itself.  In practice, this doesn't work
        // well.
        Array<double> delta_actual(numWeights, 0.0);
        for (int i = 0; i < numWeights; i++)
          delta_actual[i] = gradient[i] - oldgradient[i];

        preddist = dotprod(*delta_pred, dist, numWeights);
        realdist = dotprod(delta_actual, dist, numWeights);
#else
        // Assume quadratic behavior of function
        // If this test passes, then either our quadratic approx is
        // good, or we were strangely lucky.  In practice, this isn't
        // very stable, since our function doesn't act quadratically.
        Array<double> avgPred(numWeights, 0.0);
        Array<double> avgActual(numWeights, 0.0);
        for (int i = 0; i < numWeights; i++)
        {
          avgPred[i] = oldgradient[i] + (*delta_pred)[i]/2.0;
          avgActual[i] = (oldgradient[i] + gradient[i])/2.0;
        }

        preddist = dotprod(avgPred, dist, numWeights);
        realdist = dotprod(avgActual, dist, numWeights);
#endif

          // LOG: Print out the trust-region diagnostics.
        cout << "Pred*dist = " << preddist << endl;
        cout << "Real*dist = " << realdist << endl;
        if (preddist > 0)
          cout << "Distratio = " << realdist/preddist << endl;
      }

      if (iter > 1 && (method_ == CG || method_ == DN))
      {

        // Adjust lambda using technique of (Fletcher, 1987)
        double delta = realdist/preddist;

        if (!backtracked && preddist == 0)
          cg_lambda_ /= 4;

#if 0
        // In general, overshooting the predicted gain in
        // function value is a good thing.  However, since we're
        // doing everything approximately, we may wish to adjust
        // lambda so that our quadratic approximation is always
        // close.  These update formulas penalize overshooting.
        if (delta > 0.75 && delta < 1.333)
          cg_lambda_ /= 2;
        if (delta < 0.25 || delta > 4.0)
        {
          if (cg_lambda_ * 4 > cg_max_lambda_)
            cg_lambda_ = cg_max_lambda_;
          else
            cg_lambda_ *= 4;
        }
#else
        if (!backtracked && preddist != 0 && delta > 0.75)
          cg_lambda_ /= 2;
        // if (delta < 0.25)   // Update schedule from (Fletcher, 1987)
        if (delta < 0.0)       // Gentler update schedule, to handle noise
        {
          if (cg_lambda_ * 4 > cg_max_lambda_)
            cg_lambda_ = cg_max_lambda_;
          else
            cg_lambda_ *= 4;
        }
#endif
        cout << "Lambda = " << cg_lambda_ << endl;

        if (delta < 0.0 && backtrackCount_ < maxBacktracks_)
        {
          cout << "Backtracking!" << endl;
          iter--;
#if 0
          double oldalpha = alpha;
          alpha = computeQuadraticStepLength(oldgradient, d, Hd, cg_lambda_);
          for (int i = 0; i < numWeights; i++)
          {
            weights[i] = oldweights[i] + alpha * d[i];
            (*delta_pred)[i] *= alpha/oldalpha;
          }
#else
          for (int i = 0; i < numWeights; i++)
            weights[i] = oldweights[i];

#endif
          for (int i = 0; i < domainCnt_; i++)
            inferences_[i]->restoreCnts();

          backtracked = true;
          backtrackCount_++;
        }
        else
        {
          backtracked = false;
          backtrackCount_ = 0;
        }
      }

      if (!backtracked)
      {
        for (int i = 0; i < domainCnt_; i++)
          inferences_[i]->saveCnts();
      }


#if 0
      // Check for convergence.  If we're backtracking, skip these tests.
      if (!backtracked)
      {
       for (int i = 0; i < numWeights; i++)
       {
        // Smooth with Laplace prior
        double priorCounts = 1.0;
        tcounts[i]      += priorCounts;
        tcountsTotal[i] += priorCounts * 2.0;
        icounts[i]      += priorCounts;
        icountsTotal[i] += priorCounts * 2.0;

        // Adjust training counts by the Gaussian prior
        double sd = priorStdDevs_[i];
        double effectiveCounts = tcounts[i] -
            (weights[i] - priorMeans_[i])/(sd*sd);

        if (tcountsTotal[i] == 0)
          continue;

        double p = effectiveCounts/tcountsTotal[i];
        if (p < 0.0) { p = 0.0; }
        if (p > 1.0) { p = 1.0; }
        double cdf = // compute CDF for binomial here...

        // TODO: make the significance level a command-line parameter
        if (cdf < 0.05 || cdf > 0.95)
        {
          // DEBUG
          cout << "i = " << i << endl;
          cout << "tc' = " << effectiveCounts << endl;
          cout << "tc = " << tcounts[i] << "; tct = " << tcountsTotal[i] << endl;
          cout << "ic = " << icounts[i] << "; ict = " << icountsTotal[i] << endl;
          cout << "Base prob: " << p << endl;
          cout << "Inferred prob: " << (double)icounts[i]/icountsTotal[i] << endl;
          cout << "cdf = " << cdf << endl;

          itersSinceSignificant = 0;
          break;
        }
       }

       // If no difference has been significant for 10 iterations, stop.
       // TODO: make the constant 10 a command-line parameter.
       itersSinceSignificant++;
       if (itersSinceSignificant > 10)
       {
          cout << "No significant difference in any dimension after "
              << "10 iterations; halting!";
          break;
       }
      }
#endif


        // Scaled conjugate gradient: compute conjugate gradient direction
        //     using Polak-Ribiere and step size using Hessian matrix.
        //     Optional preconditioner (inverse diagonalized Hessian)
      if (method_ == CG && !backtracked)
      {
        // Precond stores the diagonal entries of preconditioning matrix
        Array<double> precond(numWeights, 1.0);

        // If preconditioning, use the inverse diagonalized Hessian
        if (preconditionCG_)
        {
          Array<double> variance(numWeights, 0.0);
          getVariance(variance, numWeights);

          for (int clauseno = 0; clauseno < numWeights; clauseno++)
            precond[clauseno] = 1.0/variance[clauseno];
            //precond[clauseno] = 1.0/(variance[clauseno] + 1.0);
        }

        double beta = 0.0;

          // Compute beta using Polak-Ribiere form:
          //   beta = g_j+1 (g_j+1 - g_j) / (g_j g_j)
          // Preconditioned:
          //   beta = g_j+1 M-1 (g_j+1 - g_j) / (g_j M-1 g_j)
        double beta_num = 0.0;
        double beta_denom = 0.0;

        if (iter > 1)
        {
          for (int i = 0; i < numWeights; i++)
          {
            beta_num   += gradient[i] * precond[i]
                * (gradient[i] - oldgradient[i]);
            beta_denom += oldgradient[i] * precond[i] * oldgradient[i];
          }
          beta = beta_num/beta_denom;
        }
        else
          beta = 0.0;

        if (dldebug)
          cout << "Beta = " << beta << endl;

          // Compute new direction
        for (int w = 0; w < numWeights; w++)
          d[w] = -precond[w]*gradient[w] + beta*oldd[w];

        delete Hd;
        Hd = getHessianVectorProduct(d);
        alpha = computeQuadraticStepLength(gradient, d, Hd, cg_lambda_);

        // HACK DEBUG -- if this isn't working, ignore the conjugacy
        // criteria since it sometimes messes us up.
        if (alpha < 0.0)
        {
          for (int w = 0; w < numWeights; w++)
            d[w] = -precond[w]*gradient[w];

          delete Hd;
          Hd = getHessianVectorProduct(d);
          alpha = computeQuadraticStepLength(gradient, d, Hd, cg_lambda_);
        }
      }

        // Diagonal Newton: divide gradient by clause variance,
        //     and use Hessian to pick step length
      if (method_ == DN && !backtracked)
      {
        Array<double> variance(numWeights, 0.0);
        getVariance(variance, numWeights);

        for (int w = 0; w < numWeights; w++)
        {
          if (fabs(gradient[w]) < 0.00000001)
            gradient[w] = 0;
          d[w] = -gradient[w]/variance[w];
        }

        delete Hd;
        Hd = getHessianVectorProduct(d);
        alpha = computeQuadraticStepLength(gradient, d, Hd, cg_lambda_);
        // alpha = learningRate; (Alternate method for choosing step length.)
      }

      if (method_ == SIMPLE)
      {
          // Search direction is just the gradient
        for (int w = 0; w < numWeights; w++)
          d[w] = -gradient[w];

          // Step length scalar is the learning rate
        alpha = learningRate;
      }

      if (!backtracked && alpha <= 0.0)
      {
        // If alpha is negative, then either the direction or the
        // Hessian is in error.  We call this a backtrack so that
        // we can gather more samples while keeping the old samples.
        backtracked = true;
        iter--;
      }

      if (!backtracked)
      {
          // Compute total weight change
        Array<double> wchange(numWeights, 0.0);
        for (int w = 0; w < numWeights; w++) {
			const Clause* clause=mln->getClause(w);
			if (clause->isHardClause()) continue;
			wchange[w] = d[w] * alpha + (weights[w] - oldweights[w]) * momentum;
		}

        // Convergence criteria for 2nd order methods:
        // Stop when the maximum predicted improvement in log likelihood
        // is very small.
        double maxchange = -dotprod(gradient, wchange, numWeights);
        cout << "Maximum estimated improvement: " << maxchange << endl;
        if ((method_ == CG || method_ == DN) &&
            maxchange < min_ll_change_)
        {
          cout << "Upper bound is less than " << min_ll_change_
               << "; halting learning." << endl;
          break;
        }

        // Save weights, gradient, and direction and adjust the weights
        for (int w = 0; w < numWeights; w++)
        {
          oldweights[w] = weights[w];
          oldd[w] = d[w];
          oldgradient[w] = gradient[w];

          weights[w] += wchange[w];
          averageWeights[w] = ((iter-1) * averageWeights[w] + weights[w])/iter;
        }

        // Predict the next change to the gradient using the Hessian
        delete delta_pred;
        delta_pred = getHessianVectorProduct(wchange);

        // Reset counts.
        // (If we're backtracking, we DO NOT want to
        // reset counts, since we're restoring the old counts.
        // That's why this is handled here.)
        for (int i = 0; i < domainCnt_; i++)
          inferences_[i]->resetCnts();
      }

      bool print = !backtracked && (!periodicMLNs || isOutputIter(iter));
      if (print)
      {
        cout << "Weights:\tNew\t\tAverage\n";
        for (int w = 0; w < numWeights; w++)
          cout << w << ":\t\t" << weights[w] << "\t\t"
               << averageWeights[w] << endl;
      }

      // done with an iteration
    }

    cout << endl << "Learned Weights : " << endl;
    for (int w = 0; w < numWeights; w++)
    {
        // If using MWS, then assign average weights to MLN
      if (dynamic_cast<MaxWalkSat*>(inferences_[0]))
      {
        weights[w] = averageWeights[w];
      }
      cout << w << ":" << weights[w] << endl;
    }

    delete [] averageWeights;
    delete [] gradient;

    resetDBs();
  }


 private:

  /**
   * Resets the values of non-evidence predicates as they were before learning.
   */
  void resetDBs()
  {
    if (!lazyInference_)
    {
      for (int i = 0; i < domainCnt_; i++)
      {
          // Check if using BP
        if ((bp_ = dynamic_cast<BP*>(inferences_[i])))
        {
          continue;
        }

        VariableState* state = inferences_[i]->getState();
        Database* db = state->getDomain()->getDB();
          // Change known NE to original values
        const GroundPredicateHashArray* knePreds = state->getKnePreds();
        const Array<TruthValue>* knePredValues = state->getKnePredValues();
        db->setValuesToGivenValues(knePreds, knePredValues);
          // Set unknown NE back to UKNOWN
        const GroundPredicateHashArray* unePreds = state->getUnePreds();
        for (int predno = 0; predno < unePreds->size(); predno++)
          db->setValue((*unePreds)[predno], UNKNOWN);
      }
    }
  }

  /**
   * Assign true to the elements in the relevantClauses_ bool array
   * corresponding to indices of clauses which would be relevant for list of
   * non-evidence predicates.
   */
  void findRelevantClauses(const StringHashArray& nonEvidPredNames)
  {
    for (int d = 0; d < domainCnt_; d++)
    {
      int clauseCnt = clauseCntPerDomain_[d];
      Array<bool>& relevantClauses = relevantClausesPerDomain_[d];
      relevantClauses.growToSize(clauseCnt);
      memset((bool*)relevantClauses.getItems(), false,
             relevantClauses.size()*sizeof(bool));
      const Domain* domain;
      const MLN* mln;

        // Check if using BP
      if ((bp_ = dynamic_cast<BP*>(inferences_[d])))
      {
        domain = bp_->getFactorGraph()->getDomain();
        mln = bp_->getFactorGraph()->getMLN();
      }
      else
      {
        domain = inferences_[d]->getState()->getDomain();
        mln = inferences_[d]->getState()->getMLN();
      }

      const Array<IndexClause*>* indclauses;
      const Clause* clause;
      int predid, clauseid;
      for (int i = 0; i < nonEvidPredNames.size(); i++)
      {
        predid = domain->getPredicateId(nonEvidPredNames[i].c_str());
        //cout << "finding the relevant clauses for predid = " << predid
        //     << " in domain " << d << endl;
        indclauses = mln->getClausesContainingPred(predid);
        if (indclauses)
        {
          for (int j = 0; j < indclauses->size(); j++)
          {
            clause = (*indclauses)[j]->clause;
            clauseid = mln->findClauseIdx(clause);

              // If clause is external to this mln, then it stays irrelevant
            if (!mln->isExternalClause(clauseid))
              relevantClauses[clauseid] = true;
          }
        }
      }
    }
  }


  void findRelevantClausesFormulas()
  {
    if (idxTrans_ == NULL)
    {
      Array<bool>& relevantClauses = relevantClausesPerDomain_[0];
      relevantClausesFormulas_.growToSize(relevantClauses.size());
      for (int i = 0; i < relevantClauses.size(); i++)
        relevantClausesFormulas_[i] = relevantClauses[i];
    }
    else
    {
      idxTrans_->setRelevantClausesFormulas(relevantClausesFormulas_,
                                            relevantClausesPerDomain_[0]);
      cout << "Relevant clauses/formulas:" << endl;
      idxTrans_->printRelevantClausesFormulas(cout, relevantClausesFormulas_);
      cout << endl;
    }
  }


  /**
   * Calculate true/false/unknown counts for all clauses for the given domain.
   *
   * @param trueCnt Number of true groundings for each clause is stored here.
   * @param falseCnt Number of false groundings for each clause is stored here.
   * @param domainIdx Index of domain where the groundings are counted.
   * @param hasUnknownPreds If true, the domain has predicates with unknown
   * truth values. Otherwise it contains only predicates with known values.
   */
  void calculateCounts(Array<double>& trueCnt, Array<double>& falseCnt,
                       const int& domainIdx, const bool& hasUnknownPreds)
  {
    Clause* clause;
    double tmpUnknownCnt;
    int clauseCnt = clauseCntPerDomain_[domainIdx];
    Array<bool>& relevantClauses = relevantClausesPerDomain_[domainIdx];
    const Domain* domain;
    const MLN* mln;

      // Check if using BP
    if ((bp_ = dynamic_cast<BP*>(inferences_[domainIdx])))
    {
      domain = bp_->getFactorGraph()->getDomain();
      mln = bp_->getFactorGraph()->getMLN();
    }
    else
    {
      domain = inferences_[domainIdx]->getState()->getDomain();
      mln = inferences_[domainIdx]->getState()->getMLN();
    }


    for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
    {
      if (!relevantClauses[clauseno])
      {
        continue;
        //cout << "\n\nthis is an irrelevant clause.." << endl;
      }
      clause = (Clause*) mln->getClause(clauseno);
      clause->getNumTrueFalseUnknownGroundings(domain, domain->getDB(),
                                               hasUnknownPreds,
                                               trueCnt[clauseno],
                                               falseCnt[clauseno],
                                               tmpUnknownCnt);
      assert(hasUnknownPreds || (tmpUnknownCnt==0));
    }
  }


  /**
   * Initializes the weights of the clauses. Assumption is that the inferences_
   * (and their states) have been constructed.
   *
   * @param nonEvidPredNames Names of the non-evidence predicate names in the
   * domain.
   */
  void initializeWts(const StringHashArray& nonEvidPredNames)
  {
    cout << "Initializing weights ..." << endl;

    bool hasUnknownPreds;

    Array<Predicate*> gpreds;
    Array<Predicate*> ppreds;
    Array<TruthValue> gpredValues;
    Array<TruthValue> tmpValues;

    if (!lazyInference_)
    { // Eager inference
      trainTrueCnts_.growToSize(domainCnt_);
      trainFalseCnts_.growToSize(domainCnt_);
      if (withEM_) emDiffCnts_.growToSize(domainCnt_);
    }

    for (int i = 0; i < domainCnt_; i++)
    {
      if (dldebug) cout << "Domain " << i << endl;
      int clauseCnt = clauseCntPerDomain_[i];
      if (dldebug) cout << "Clause count: " << clauseCnt << endl;
      Domain* domain;
        // Check if using BP
      if ((bp_ = dynamic_cast<BP*>(inferences_[i])))
      {
        FactorGraph* fg = bp_->getFactorGraph();
        domain = (Domain*)fg->getDomain();
      }
      else
      {
        VariableState* state = inferences_[i]->getState();
        domain = (Domain*)state->getDomain();
      }

      if (lazyInference_ || (bp_ = dynamic_cast<BP*>(inferences_[i])))
      {
        domain->getDB()->setPerformingInference(false);

        //cout << endl << "Getting the counts for the domain " << i << endl;
        gpreds.clear();
        gpredValues.clear();
        tmpValues.clear();
        for (int predno = 0; predno < nonEvidPredNames.size(); predno++)
        {
          ppreds.clear();
          int predid = domain->getPredicateId(nonEvidPredNames[predno].c_str());
          Predicate::createAllGroundings(predid, domain, ppreds);
          gpreds.append(ppreds);
        }

        domain->getDB()->alterTruthValue(&gpreds, UNKNOWN, FALSE, &gpredValues);

        hasUnknownPreds = false;

        Array<double>& trueCnt = totalTrueCnts_[i];
        Array<double>& falseCnt = totalFalseCnts_[i];
        trueCnt.growToSize(clauseCnt);
        falseCnt.growToSize(clauseCnt);
        calculateCounts(trueCnt, falseCnt, i, hasUnknownPreds);

        //cout << "got the total counts..\n\n\n" << endl;

        hasUnknownPreds = true;

        domain->getDB()->setValuesToUnknown(&gpreds, &tmpValues);

        Array<double>& dTrueCnt = defaultTrueCnts_[i];
        Array<double>& dFalseCnt = defaultFalseCnts_[i];
        dTrueCnt.growToSize(clauseCnt);
        dFalseCnt.growToSize(clauseCnt);
        calculateCounts(dTrueCnt, dFalseCnt, i, hasUnknownPreds);

        for (int predno = 0; predno < gpreds.size(); predno++)
          delete gpreds[predno];

        domain->getDB()->setPerformingInference(true);
      }
      else
      { // Eager inference
        VariableState* state = inferences_[i]->getState();
        if (withEM_)
        {
          VariableState* emState = emInferences_[i]->getState();
          const GroundPredicateHashArray* emUnePreds = emState->getUnePreds();
          const GroundPredicateHashArray* emKnePreds = emState->getKnePreds();

          emDiffCnts_[i] = new double[clauseCnt];

          if (dldebug)
          {
            if (emUnePreds)
            {
              cout << "Unknown non-evid preds (E): " << emUnePreds->size()
                   << endl;
            }
            if (emKnePreds)
            {
              cout << "Known non-evid preds (E): " << emKnePreds->size()
                   << endl;
            }
          }
          int emUnePredsSize = 0;
          int emKnePredsSize = 0;
          if (emUnePreds) emUnePredsSize = emUnePreds->size();
          if (emKnePreds) emKnePredsSize = emKnePreds->size();

          int emTotalPreds = emUnePredsSize + emKnePredsSize;
          // Used to store gnd preds to be ignored in the count because they
          // are UNKNOWN
          Array<bool>* emUnknownPred = new Array<bool>;
          emUnknownPred->growToSize(emTotalPreds, false);
          for (int predno = 0; predno < emTotalPreds; predno++)
          {
            GroundPredicate* p;
            if (predno < emUnePredsSize)
              p = (*emUnePreds)[predno];
            else
              p = (*emKnePreds)[predno - emUnePredsSize];
            TruthValue tv = emState->getDomain()->getDB()->getValue(p);

            //assert(tv != UNKNOWN);
            if (tv == TRUE)
            {
              emState->setValueOfAtom(predno + 1, true, false, -1);
              p->setTruthValue(true);
            }
            else
            {
              emState->setValueOfAtom(predno + 1, false, false, -1);
              p->setTruthValue(false);
                // Can have unknown truth values when using EM. We want to
                // ignore these when performing the counts
              if (tv == UNKNOWN)
              {
                (*emUnknownPred)[predno] = true;
              }
            }
          }

          emState->initMakeBreakCostWatch();
          //cout<<"getting true cnts => "<<endl;
          emState->getNumClauseGndingsWithUnknown(emDiffCnts_[i], clauseCnt,
                                                  true, emUnknownPred);
          delete emUnknownPred;
          if (dldebug)
          {
            for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
            {
              cout << clauseno << " : emDiffCnt = " << emDiffCnts_[i][clauseno]
                   << endl;
            }
          }
        }

        const GroundPredicateHashArray* unePreds = state->getUnePreds();
        const GroundPredicateHashArray* knePreds = state->getKnePreds();

        trainTrueCnts_[i] = new double[clauseCnt];
        trainFalseCnts_[i] = new double[clauseCnt];

        if (dldebug)
        {
          cout << "Unknown non-evid preds: " << unePreds->size() << endl;
          cout << "Known non-evid preds: " << knePreds->size() << endl;
        }
        int totalPreds = unePreds->size() + knePreds->size();
          // Used to store gnd preds to be ignored in the count because they are
          // UNKNOWN
        Array<bool>* unknownPred = new Array<bool>;
        unknownPred->growToSize(totalPreds, false);
        for (int predno = 0; predno < totalPreds; predno++)
        {
          GroundPredicate* p;
          if (predno < unePreds->size())
            p = (*unePreds)[predno];
          else
            p = (*knePreds)[predno - unePreds->size()];
          TruthValue tv = state->getDomain()->getDB()->getValue(p);

          //assert(tv != UNKNOWN);
          //bool activate = true;
          bool activate = false;
          if (tv == TRUE)
          {
            state->setValueOfAtom(predno + 1, true, activate, -1);
            p->setTruthValue(true);
          }
          else
          {
            state->setValueOfAtom(predno + 1, false, activate, -1);
            p->setTruthValue(false);
              // Can have unknown truth values when using EM. We want to ignore
              // these when performing the counts
            if (tv == UNKNOWN)
            {
              (*unknownPred)[predno] = true;
            }
          }
        }

        state->initMakeBreakCostWatch();
        //cout<<"getting true cnts => "<<endl;
        state->getNumClauseGndingsWithUnknown(trainTrueCnts_[i], clauseCnt,
                                              true, unknownPred);
        if (withEM_)
        {
          for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
          {
            emDiffCnts_[i][clauseno] = trainTrueCnts_[i][clauseno] -
                                       emDiffCnts_[i][clauseno];
          }
        }
        //cout<<endl;
        //cout<<"getting false cnts => "<<endl;
        state->getNumClauseGndingsWithUnknown(trainFalseCnts_[i], clauseCnt,
                                              false, unknownPred);
        delete unknownPred;
        if (dldebug)
        {
          for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
          {
            cout << clauseno << " : tc = " << trainTrueCnts_[i][clauseno]
                 << " ** fc = " << trainFalseCnts_[i][clauseno] << endl;
          }
        }
      }
    }

    double tc,fc;
    int nonEvidPreds = 0;
    cout << "List of CNF Clauses : " << endl;
    for (int clauseno = 0; clauseno < clauseCntPerDomain_[0]; clauseno++)
    {
      tc = 0.0; fc = 0.0; nonEvidPreds = 0;
      Domain* domain0;
        // Check if using BP
      if ((bp_ = dynamic_cast<BP*>(inferences_[0])))
      {
        domain0 = (Domain*)bp_->getFactorGraph()->getDomain();
      }
      else
      {
        domain0 = (Domain*)inferences_[0]->getState()->getDomain();
      }

      for (int i = 0; i < domainCnt_; i++)
      {
        const MLN* mln;
        Domain* domain;
          // Check if using BP
        if ((bp_ = dynamic_cast<BP*>(inferences_[i])))
        {
          FactorGraph* fg = bp_->getFactorGraph();
          domain = (Domain*)fg->getDomain();
          mln = fg->getMLN();
        }
        else
        {
          VariableState* state = inferences_[i]->getState();
          domain = (Domain*)state->getDomain();
          mln = state->getMLN();
        }

        Array<bool>& relevantClauses = relevantClausesPerDomain_[i];
        Array<double>& logOdds = logOddsPerDomain_[i];

        if (!relevantClauses[clauseno])
        {
          domain->setNumTrueNonEvidGndings(clauseno, 0);
          domain->setNumFalseNonEvidGndings(clauseno, 0);
          logOdds[clauseno] = 0.0;
          continue;
        }

        cout << clauseno << ":";
        const Clause* clause = mln->getClause(clauseno);
        clause->print(cout, domain0);
        cout << endl;

        double domainTrueCnts = 0.0;
        double domainFalseCnts = 0.0;
        if (lazyInference_ || (bp_ = dynamic_cast<BP*>(inferences_[i])))
        {
          domainTrueCnts = totalTrueCnts_[i][clauseno] -
                           defaultTrueCnts_[i][clauseno];
          domainFalseCnts = totalFalseCnts_[i][clauseno] -
                            defaultFalseCnts_[i][clauseno];
        }
        else
        {
          domainTrueCnts = trainTrueCnts_[i][clauseno];
          domainFalseCnts = trainFalseCnts_[i][clauseno];
        }

        tc += domainTrueCnts;
        fc += domainFalseCnts;
        domain->setNumTrueNonEvidGndings(clauseno, domainTrueCnts);
        domain->setNumFalseNonEvidGndings(clauseno, domainFalseCnts);

        if (dldebug)
          cout << clauseno << " : tc = " << tc << " ** fc = "<< fc <<endl;

        for (int i = 0; i < clause->getNumPredicates(); i++)
        {
          const char* predName =
            clause->getPredicate(i)->getTemplate()->getName();
          if (nonEvidPredNames.contains(predName))
            nonEvidPreds++;
        }
      }

      double weight = 0.0;
      double totalCnt = tc + fc;

      if (totalCnt == 0 || nonEvidPreds == 0)
      {
        //cout << "NOTE: Total count is 0 for clause " << clauseno << endl;
        weight = EPSILON;
      }
      else
      {
        if (fc == 0.0)
          fc = 0.00001;
        if (tc == 0.0)
          tc = 0.00001;

        double priorodds = (pow(2.0, nonEvidPreds) - 1.0)/1.0;
        weight = log(tc/fc) - log(priorodds);
      }

      for (int i = 0; i < domainCnt_; i++)
      {
      	Array<double>& logOdds = logOddsPerDomain_[i];
        logOdds[clauseno] = weight;
      }
    }
  }


  /**
   * Runs inference using the current set of parameters.
   */
  void infer()
  {
    for (int i = 0; i < domainCnt_; i++)
    {
      if (dldebug) cout << "in domain " << i << endl;

        // Check if using BP
      if ((bp_ = dynamic_cast<BP*>(inferences_[i])))
      {
        inferences_[i]->infer();
      }
      else
      {
        VariableState* state = inferences_[i]->getState();
        state->setGndClausesWtsToSumOfParentWts();

        // MWS: Search is started from state at end of last iteration
        state->init();
        inferences_[i]->infer();
        state->saveLowStateToGndPreds();
      }

      if (dldebug)
      {
        cout << "Inferred following values (M): " << endl;
        inferences_[i]->printProbabilities(cout);
      }
    }
  }

  /**
   * Infers values for predicates with unknown truth values and uses these
   * values to compute the training counts.
   */
  void fillInMissingValues()
  {
    assert(withEM_);
    cout << "Filling in missing data ..." << endl;
      // Get values of initial unknown preds by producing MAP state of
      // unknown preds given known evidence and non-evidence preds (VPEM)
    Array<Array<TruthValue> > ueValues;
    ueValues.growToSize(domainCnt_);
    for (int i = 0; i < domainCnt_; i++)
    {
      VariableState* state = inferences_[i]->getState();
      const Domain* domain = state->getDomain();
      const GroundPredicateHashArray* knePreds = state->getKnePreds();
      const Array<TruthValue>* knePredValues = state->getKnePredValues();

        // Mark known non-evidence preds as evidence
      domain->getDB()->setValuesToGivenValues(knePreds, knePredValues);

      VariableState* emState = emInferences_[i]->getState();
        // Infer missing values
      emState->setGndClausesWtsToSumOfParentWts();
        // MWS: Search is started from state at end of last iteration
      emState->init();
      emInferences_[i]->infer();
      emState->saveLowStateToGndPreds();

      if (dldebug)
      {
        cout << "Inferred following values (E): " << endl;
        emInferences_[i]->printProbabilities(cout);
      }

        // Compute counts
      int clauseCnt = clauseCntPerDomain_[i];
      emState->initMakeBreakCostWatch();
      //cout<<"getting true cnts => "<<endl;
      const Array<double>* clauseTrueCnts =
        emInferences_[i]->getClauseTrueCnts();
      assert(clauseTrueCnts->size() == clauseCnt);
      int numSamples = emInferences_[i]->getNumSamples();
      for (int j = 0; j < clauseCnt; j++)
      {
        //trainTrueCnts_[i][j] = (*clauseTrueCnts)[j]/numSamples;
        trainTrueCnts_[i][j] = (*clauseTrueCnts)[j]/numSamples +
                               emDiffCnts_[i][j];
      }

        // Set evidence values back
      //assert(uePreds.size() == ueValues[i].size());
      //domain->getDB()->setValuesToGivenValues(&uePreds, &ueValues[i]);
        // Set non-evidence values to unknown
      Array<TruthValue> tmpValues;
      tmpValues.growToSize(knePreds->size());
      domain->getDB()->setValuesToUnknown(knePreds, &tmpValues);
    }
    cout << "Done filling in missing data" << endl;
  }

#if 0
  // This was used by the statistical significance convergence test,
  // to see if our counts could be drawn from the same distribution
  // as the true counts.

  void getCountsForDomain(int* const & clauseTrainCnts,
          int* const & clauseTrainTotal, int* const & clauseInferredCnts,
          int* const & clauseInferredTotal, int domainIdx)
  {
    Array<bool>& relevantClauses = relevantClausesPerDomain_[domainIdx];
    int clauseCnt = clauseCntPerDomain_[domainIdx];
    double* trainCnts = NULL;
    double* inferredCnts = NULL;
    Array<double>& totalTrueCnts = totalTrueCnts_[domainIdx];
    //Array<double>& totalFalseCnts = totalFalseCnts_[domainIdx];
    Array<double>& defaultTrueCnts = defaultTrueCnts_[domainIdx];
    const MLN* mln = inferences_[domainIdx]->getState()->getMLN();
    const Domain* domain = inferences_[domainIdx]->getState()->getDomain();

    memset(clauseTrainCnts, 0, clauseCnt*sizeof(double));
    memset(clauseInferredCnts, 0, clauseCnt*sizeof(double));
    memset(clauseTrainTotal, 0, clauseCnt*sizeof(double));
    memset(clauseInferredTotal, 0, clauseCnt*sizeof(double));

    if (!(lazyInference_ || (bp_ = dynamic_cast<BP*>(inferences_[domainIdx]))))
    {
      if (!inferredCnts) inferredCnts = new double[clauseCnt];

      const Array<double>* clauseTrueCnts =
        inferences_[domainIdx]->getClauseTrueCnts();
      assert(clauseTrueCnts->size() == clauseCnt);
      for (int i = 0; i < clauseCnt; i++)
        inferredCnts[i] = (*clauseTrueCnts)[i];

      trainCnts = trainTrueCnts_[domainIdx];
    }
      //loop over all the training examples
    //cout << "\t\ttrain count\t\t\t\tinferred count" << endl << endl;
    for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
    {
      if (!relevantClauses[clauseno]) continue;

      // Compute total groundings
      int totalGndings = (int)(trainFalseCnts_[domainIdx][clauseno]
              + trainTrueCnts_[domainIdx][clauseno]);
      clauseTrainTotal[clauseno] += totalGndings;
      clauseInferredTotal[clauseno] += totalGndings *
          inferences_[domainIdx]->getNumSamples();

      if (lazyInference_ || (bp_ = dynamic_cast<BP*>(inferences_[domainIdx])))
      {
      	Clause* clause = (Clause*) mln->getClause(clauseno);

      	double trainCnt = totalTrueCnts[clauseno];
      	double inferredCnt =
          clause->getNumTrueGroundings(domain, domain->getDB(), false);
      	trainCnt -= defaultTrueCnts[clauseno];
      	inferredCnt -= defaultTrueCnts[clauseno];
      	clauseTrainCnts[clauseno] += (int)trainCnt;
      	clauseInferredCnts[clauseno] += (int)inferredCnt;
      }
      else
      {
        clauseTrainCnts[clauseno] += (int)trainCnts[clauseno];
        clauseInferredCnts[clauseno] += (int)inferredCnts[clauseno];
      }
      //cout << clauseno << ":\t\t" <<trainCnt<<"\t\t\t\t"<<inferredCnt<<endl;
    }

    delete[] inferredCnts;
  }
#endif

  // Get the gradient of the negative log likelihood for one domain
  void getGradientForDomain(double* const & gradient, const int& domainIdx)
  {
    Array<bool>& relevantClauses = relevantClausesPerDomain_[domainIdx];
    int clauseCnt = clauseCntPerDomain_[domainIdx];
    double* trainCnts = NULL;
    double* inferredCnts = NULL;
    double* clauseTrainCnts = new double[clauseCnt];
    double* clauseInferredCnts = new double[clauseCnt];
    //double trainCnt, inferredCnt;
    //Array<double>& totalTrueCnts = totalTrueCnts_[domainIdx];
    //Array<double>& defaultTrueCnts = defaultTrueCnts_[domainIdx];
    //const MLN* mln = inferences_[domainIdx]->getState()->getMLN();
    //const Domain* domain = inferences_[domainIdx]->getState()->getDomain();

    memset(clauseTrainCnts, 0, clauseCnt*sizeof(double));
    memset(clauseInferredCnts, 0, clauseCnt*sizeof(double));

    if (!inferredCnts) inferredCnts = new double[clauseCnt];

    const Array<double>* clauseTrueCnts =
      inferences_[domainIdx]->getClauseTrueCnts();
    cout<<clauseTrueCnts->size()<<endl;
    for(int i=0;i<clauseTrueCnts->size();i++)
    {
    	cout<<(*clauseTrueCnts)[i]<<" ";
    }
    cout<<endl;
    assert(clauseTrueCnts->size() == clauseCnt);
    int numSamples = inferences_[domainIdx]->getNumSamples();
    for (int i = 0; i < clauseCnt; i++)
    {
      if ((bp_ = dynamic_cast<BP*>(inferences_[domainIdx])))
        inferredCnts[i] = (*clauseTrueCnts)[i];
      else
        inferredCnts[i] = (*clauseTrueCnts)[i] / numSamples;
    }
    if (!(lazyInference_ || (bp_ = dynamic_cast<BP*>(inferences_[domainIdx]))))
      trainCnts = trainTrueCnts_[domainIdx];

    if (dldebug)
      cout << "numSamples = " << numSamples << endl;

      //loop over all the training examples
    //cout << "\t\ttrain count\t\t\t\tinferred count" << endl << endl;
    for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
    {
      if (!relevantClauses[clauseno]) continue;
      if (lazyInference_ || (bp_ = dynamic_cast<BP*>(inferences_[domainIdx])))
      {
        clauseTrainCnts[clauseno] = clauseTrainCnts[clauseno] +
                                    totalTrueCnts_[domainIdx][clauseno] -
                                    defaultTrueCnts_[domainIdx][clauseno];
        clauseInferredCnts[clauseno] += inferredCnts[clauseno];
      }
      else
      {
        clauseTrainCnts[clauseno] += trainCnts[clauseno];
        clauseInferredCnts[clauseno] += inferredCnts[clauseno];
      }
      //cout << clauseno << ":\t\t" <<trainCnt<<"\t\t\t\t"<<inferredCnt<<endl;
    }

    if (dldebug)
    {
      cout << "Domain " << domainIdx << " net counts : " << endl;
      cout << "\t\ttrain count\t\t\t\tinferred count" << endl << endl;
    }

    for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
    {
      if (!relevantClauses[clauseno]) continue;

      if (dldebug)
        cout << clauseno << ":\t\t" << clauseTrainCnts[clauseno] << "\t\t\t\t"
             << clauseInferredCnts[clauseno] << endl;
        // NOTE: rescaling has been implemented to affect the gradient,
        // but not the Hessian.
      if (rescaleGradient_ && clauseTrainCnts[clauseno] > 0)
      {
        gradient[clauseno] -=
          (clauseTrainCnts[clauseno] - clauseInferredCnts[clauseno])
            / clauseTrainCnts[clauseno];
      }
      else
      {
        gradient[clauseno] -= clauseTrainCnts[clauseno] -
                            clauseInferredCnts[clauseno];
      }
    }

    delete[] inferredCnts;
    delete[] clauseTrainCnts;
    delete[] clauseInferredCnts;
  }


    // Get the gradient of the negative log likelihood
  void getGradient(double* const & weights, double* const & gradient,
                   const int numWts)
  {
      // Compute the gradient
    memset(gradient, 0, numWts*sizeof(double));

      // There is one DB or the clauses of multiple DBs line up
    if (idxTrans_ == NULL)
    {
      for (int i = 0; i < domainCnt_; i++)
      {
        //cout << "For domain number " << i << endl << endl;
        getGradientForDomain(gradient, i);
      }
    }
    else
    {
        // The clauses for multiple databases do not line up
      Array<Array<double> >* gradsPerDomain = idxTrans_->getGradsPerDomain();
      const Array<Array<Array<IdxDiv>*> >* cIdxToCFIdxsPerDomain
        = idxTrans_->getClauseIdxToClauseFormulaIdxsPerDomain();

      for (int i = 0; i < domainCnt_; i++)
      {
        //cout << "For domain number " << i << endl << endl;

        Array<double>& grads = (*gradsPerDomain)[i];
        memset((double*)grads.getItems(), 0, grads.size()*sizeof(double));

        getGradientForDomain((double*)grads.getItems(), i);

          // map clause gradient to clause/formula gradients
        assert(grads.size() == clauseCntPerDomain_[i]);
        for (int j = 0; j < grads.size(); j++)
        {
          Array<IdxDiv>* idxDivs = (*cIdxToCFIdxsPerDomain)[i][j];
          for (int k = 0; k < idxDivs->size(); k++)
            gradient[ (*idxDivs)[k].idx ] += grads[j] / (*idxDivs)[k].div;
        }
      }
    }

      // Add the deriative of the prior
    if (usePrior_)
    {
	  for (int i = 0; i < numWts; i++)
      {
        //if (!relevantClausesFormulas_[i]) continue;
        double sd = priorStdDevs_[i];
        double priorDerivative = (weights[i]-priorMeans_[i])/(sd*sd);
        //cout << i << " : " << "gradient : " << gradient[i]
        //     << "  prior gradient : " << priorDerivative;
        gradient[i] += priorDerivative;
	    //cout << "  net gradient : " << gradient[i] << endl;
      }
    }
  }

 public:
    // Different gradient descent methods
    // SIMPLE is Voted Perceptron
  const static int SIMPLE = 0;
  const static int DN = 2;
  const static int CG = 3;


 private:
  int domainCnt_;
  Array<Array<double> > logOddsPerDomain_;
  Array<int> clauseCntPerDomain_;

	// Used in lazy version
  Array<Array<double> > totalTrueCnts_;
  Array<Array<double> > totalFalseCnts_;
  Array<Array<double> > defaultTrueCnts_;
  Array<Array<double> > defaultFalseCnts_;

  Array<Array<bool> > relevantClausesPerDomain_;
  Array<bool> relevantClausesFormulas_;

	// Used to compute cnts from mrf
  Array<double*> trainTrueCnts_;
  Array<double*> trainFalseCnts_;

  bool usePrior_;
  const double* priorMeans_, * priorStdDevs_;

  IndexTranslator* idxTrans_; //not owned by object; don't delete

  bool lazyInference_;
  bool isQueryEvidence_;
  bool rescaleGradient_;
  int  method_;
  double cg_lambda_;
  bool preconditionCG_;

  int maxBacktracks_;
  int backtrackCount_;
  double cg_max_lambda_;
  double min_ll_change_;

  Array<Inference*> inferences_;

    // Using EM to fill in missing values?
  bool withEM_;
  Array<Inference*> emInferences_;
    // True counts from clauses satisfied by
  Array<double*> emDiffCnts_;
    // Used to check if BP is being used (handled differently, because it works
    // on a factor graph
  BP* bp_;
};


#endif
