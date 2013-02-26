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
#ifndef STRUCTLEARN_H_NOV_2_2005
#define STRUCTLEARN_H_NOV_2_2005

#include "clausefactory.h"
#include "mln.h"
//#include "pseudologlikelihood.h"
#include "lbfgsb.h"
#include "discriminativelearner.h"

//class ExistFormula;
/////////////////////// code to handle existential formulas //////////////////
struct ExistFormula
{
  ExistFormula(const string& fformula) 
    : formula(fformula), gain(0), wt(0), newScore(0), numPreds(0) {}
  ~ExistFormula() 
  { 
    for (int i = 0; i < cnfClausesForDomains.size(); i++)
      cnfClausesForDomains[i].deleteItemsAndClear();
  }

  string formula;
    //cnfClausesForDomains[d] are the CNF clauses formed with the constants 
    //in domain d
  Array<Array<Clause*> > cnfClausesForDomains;
  double gain; 
  double wt;
  Array<double> wts;
  double newScore;
  int numPreds;
};

struct IndexCountDomainIdx { IndexAndCount* iac; int domainIdx; };
struct ClauseAndICDArray {Clause* clause;Array<IndexCountDomainIdx> icdArray;};


class StructLearn
{
 public:
    //If there are more than one domain, there must be a one-to-one 
    //correspondence among the clauses in mlns.
  StructLearn(Array<MLN*>* const & mlns, const bool& startFromEmptyMLN,
              const string& outMLNFileName, Array<Domain*>* const & domains,
              const Array<string>* const & nonEvidPredNames,
              const int& maxVars, const int& maxNumPredicates,
              const bool& cacheClauses, const double& maxCacheSizeMB,
              const bool& tryAllFlips, const bool& sampleClauses, 
              const double& delta, const double& epsilon, 
              const int& minClauseSamples, const int& maxClauseSamples,
              const bool& hasPrior, const double& priorMean, 
              const double& priorStdDev, const bool& wtPredsEqually, 
              const int& lbMaxIter, const double& lbConvThresh,
              const int& looseMaxIter, const double& looseConvThresh, 
              const int& beamSize, const int& bestGainUnchangedLimit,
              const int& numEstBestClauses, 
              const double& minWt, const double& penalty,
              const bool& sampleGndPreds, const double& fraction, 
              const int& minGndPredSamples, const int& maxGndPredSamples,
              const bool& reEvaluateBestCandidatesWithTightParams,
              const bool& structGradDescent, const bool& withEM)
    : mln0_((*mlns)[0]), mlns_(mlns), startFromEmptyMLN_(startFromEmptyMLN), 
      outMLNFileName_(outMLNFileName), domains_(domains), 
      preds_(new Array<Predicate*>), areNonEvidPreds_(new Array<bool>),
      clauseFactory_(new ClauseFactory(maxVars, maxNumPredicates,
                     (*domains_)[0])), 
      cacheClauses_(cacheClauses), origCacheClauses_(cacheClauses), 
      cachedClauses_((cacheClauses) ? (new ClauseHashArray) : NULL), 
      cacheSizeMB_(0), maxCacheSizeMB_(maxCacheSizeMB), 
      tryAllFlips_(tryAllFlips), 
      sampleClauses_(sampleClauses), origSampleClauses_(sampleClauses),
      pll_(NULL), hasPrior_(hasPrior), priorMean_(priorMean), 
      priorStdDev_(priorStdDev), wtPredsEqually_(wtPredsEqually), 
      lbfgsb_(NULL), lbMaxIter_(lbMaxIter), lbConvThresh_(lbConvThresh), 
      looseMaxIter_(looseMaxIter), looseConvThresh_(looseConvThresh),
      beamSize_(beamSize), bestGainUnchangedLimit_(bestGainUnchangedLimit),
      numEstBestClauses_(numEstBestClauses), 
      minGain_(0), minWt_(minWt), penalty_(penalty), 
      sampleGndPreds_(sampleGndPreds), fraction_(fraction), 
      minGndPredSamples_(minGndPredSamples), 
      maxGndPredSamples_(maxGndPredSamples),
      reEvalBestCandsWithTightParams_(reEvaluateBestCandidatesWithTightParams),
      candCnt_(0), iter_(-1), bsiter_(-1), startSec_(-1), indexTrans_(NULL),
      structGradDescent_(structGradDescent), withEM_(withEM)
  { 
    assert(minWt_ >= 0);
    assert(domains_->size() == mlns_->size());
    
    areNonEvidPreds_->growToSize((*domains_)[0]->getNumPredicates(), false);
    for (int i = 0; i < nonEvidPredNames->size(); i++)
    {
      int predId=(*domains_)[0]->getPredicateId((*nonEvidPredNames)[i].c_str());
      if (predId < 0)
      {
        cout << "ERROR: in StructLearn::StructLearn(). Predicate " 
             << (*nonEvidPredNames)[i] << " undefined." << endl;
        exit(-1);
      }
      (*areNonEvidPreds_)[predId] = true;
    }
    
    (*domains_)[0]->createPredicates(preds_, true);

    if (origSampleClauses_)
    {
      ClauseSampler* cs = new ClauseSampler(delta, epsilon, minClauseSamples, 
                                            maxClauseSamples);
      Clause::setClauseSampler(cs);
      for (int i = 0; i < domains_->size(); i++)
        (*domains_)[i]->newTrueFalseGroundingsStore();
    }    
  }

  
  ~StructLearn() 
  {
    if (pll_) delete pll_;
    if (lbfgsb_) delete lbfgsb_;
    preds_->deleteItemsAndClear();
    delete preds_;
    delete areNonEvidPreds_;
    delete clauseFactory_;
    if (cachedClauses_)
    {
      cachedClauses_->deleteItemsAndClear();
      delete cachedClauses_;
    }
    if (origSampleClauses_) delete Clause::getClauseSampler();
    if (indexTrans_) delete indexTrans_;
  }

  void run()
  {
    startSec_ = timer_.time();

    bool needIndexTrans = IndexTranslator::needIndexTranslator(*mlns_,*domains_);

      //if we are starting from an empty MLN, and including MLN clauses among 
      //the candidates in first step of beam search
    Array<Clause*> initialMLNClauses; 
    Array<ExistFormula*> existFormulas;
    if (startFromEmptyMLN_)
    {
      getMLNClauses(initialMLNClauses, existFormulas);
      removeClausesFromMLNs();
      for (int i = 0; i < initialMLNClauses.size(); i++) 
      {
        Clause* c = initialMLNClauses[i];
        c->newAuxClauseData(); 
        c->trackConstants();
        c->getAuxClauseData()->gain = 0;
        c->getAuxClauseData()->op = OP_ADD;
        c->getAuxClauseData()->removedClauseIdx = -1;
        c->getAuxClauseData()->hasBeenExpanded = false;
        c->getAuxClauseData()->lastStepExpanded = -1;
        c->getAuxClauseData()->lastStepOverMinWeight = -1;
      }
    }

      //add unit clauses to the MLN
    cout << "adding unit clauses to MLN..." << endl << endl;
    addUnitClausesToMLNs();

      //create auxiliary data for each clause in MLN
    for (int i = 0; i < mln0_->getNumClauses(); i++) 
    {
      Clause* c = (Clause*) mln0_->getClause(i);
      c->newAuxClauseData(); 
        //only clauses that are not in any exist. quant. formula's CNF may need 
        //to be translated across domains
      if (isModifiableClause(i)) c->trackConstants();
    }

    indexTrans_ = (needIndexTrans)? new IndexTranslator(mlns_, domains_) : NULL;
    if (indexTrans_) 
      cout << "The weights of clauses in the CNFs of existential formulas wiil "
           << "be tied" << endl;
        
    //if (structGradDescent_) runStructGradDescent(initialMLNClauses, existFormulas);
    runStructLearning(initialMLNClauses, existFormulas);
  }


  void runStructLearning(Array<Clause*> initialMLNClauses,
                         Array<ExistFormula*> existFormulas)
  {
/*
    startSec_ = timer_.time();

    bool needIndexTrans = IndexTranslator::needIndexTranslator(*mlns_,*domains_);

      //if we are starting from an empty MLN, and including MLN clauses among 
      //the candidates in first step of beam search
    Array<Clause*> initialMLNClauses; 
    Array<ExistFormula*> existFormulas;
    if (startFromEmptyMLN_)
    {
      getMLNClauses(initialMLNClauses, existFormulas);
      removeClausesFromMLNs();
      for (int i = 0; i < initialMLNClauses.size(); i++) 
      {
        Clause* c = initialMLNClauses[i];
        c->newAuxClauseData(); 
        c->trackConstants();
        c->getAuxClauseData()->gain = 0;
        c->getAuxClauseData()->op = OP_ADD;
        c->getAuxClauseData()->removedClauseIdx = -1;
      }
    }

      //add unit clauses to the MLN
    cout << "adding unit clauses to MLN..." << endl << endl;
    addUnitClausesToMLNs();

      //create auxiliary data for each clause in MLN
    for (int i = 0; i < mln0_->getNumClauses(); i++) 
    {
      Clause* c = (Clause*) mln0_->getClause(i);
      c->newAuxClauseData(); 
        //only clauses that are not in any exist. quant. formula's CNF may need 
        //to be translated across domains
      if (isModifiableClause(i)) c->trackConstants();
    }

    indexTrans_ = (needIndexTrans)? new IndexTranslator(mlns_, domains_) : NULL;
    if (indexTrans_) 
      cout << "The weights of clauses in the CNFs of existential formulas wiil "
           << "be tied" << endl;
*/
      //create PseudoLogLikelihood and LBFGSB
    pll_ = new PseudoLogLikelihood(areNonEvidPreds_, domains_, wtPredsEqually_, 
                                   sampleGndPreds_, fraction_, 
                                   minGndPredSamples_, maxGndPredSamples_);
    pll_->setIndexTranslator(indexTrans_);

    int numClausesFormulas = getNumClausesFormulas();
    lbfgsb_ = new LBFGSB(-1, -1, pll_, numClausesFormulas);

    useTightParams();

      //compute initial counts for clauses in all MLNs
    cout << "computing counts for initial MLN clauses..." << endl;
    double begSec = timer_.time();
    pllComputeCountsForInitialMLN();
    cout << "computing counts for initial MLN clauses took "; 
    timer_.printTime(cout, timer_.time()-begSec); cout << endl << endl;


      //learn initial wt/score of MLN, and set MLN clauses to learned weights
    cout << "learning the initial weights and score of MLN..." << endl << endl;
    begSec = timer_.time();
    double score;
    Array<double> wts;
    if (!learnAndSetMLNWeights(score)) return; 
    printMLNClausesWithWeightsAndScore(score, -1);
    cout << "learning the initial weights and score of MLN took ";
    timer_.printTime(cout, timer_.time()-begSec); cout << endl << endl;

    
      //add unit clauses with diff combinations of variables
    cout <<"trying to add unit clause with diff variable combinations to MLN..."
         << endl << endl;
    begSec = timer_.time();
    appendUnitClausesWithDiffCombOfVar(score); // score is updated
    printMLNClausesWithWeightsAndScore(score, -1);
    cout <<"adding unit clause with diff variable combinations to MLN took ";
    timer_.printTime(cout, timer_.time()-begSec); cout << endl << endl;


      //for each clause in MLN, try all sign flips and keep best one in MLN
    if (tryAllFlips_)
    {
      cout << "trying to flip the senses of MLN clauses..." << endl << endl;
      begSec = timer_.time();
      flipMLNClauses(score); 
      printMLNClausesWithWeightsAndScore(score, -1); // score is updated
      cout << "trying to flip the senses of MLN clauses took ";;
      timer_.printTime(cout, timer_.time()-begSec); cout << endl << endl;
    }

    iter_ = -1;
    //printMLNToFile(NULL, iter_);

      //start searching for best candidates
    double bsec;
    Array<Clause*> initialClauses;
    Array<Clause*> bestCandidates;
    while (true)
    {
      begSec = timer_.time();
      iter_++;

      //if (iter_ == 2) break; // for testing purposes only
      cout << "Iteration " << iter_ << endl << endl;

      minGain_ = 0;

      Array<ExistFormula*> highGainWtExistFormulas;
      if (startFromEmptyMLN_ && !existFormulas.empty()) 
      {
        useTightParams();
        cout << "evaluating the gains of existential formulas..." << endl<<endl;
        bsec = timer_.time(); 
          //score not updated
        minGain_ = evaluateExistFormulas(existFormulas, highGainWtExistFormulas,
                                         score);
        cout << "evaluating the gains of existential formulas took ";
        timer_.printTime(cout, timer_.time()-bsec); cout << endl << endl;
        cout << "setting minGain to min gain among existential formulas. "
             << "minGain = " << minGain_ << endl;
      }

      useLooseParams();
 
        //add clauses in MLN to clauses array
      initialClauses.clear();
      for (int i = 0; i < mln0_->getNumClauses(); i++) 
      {
        Clause* c = (Clause*) mln0_->getClause(i);
        c->getAuxClauseData()->reset();
        if (isModifiableClause(i)) 
        {
          c->trackConstants();
          initialClauses.append(new Clause(*c));
        }
      }

      bestCandidates.clear();
      beamSearch(initialClauses, initialMLNClauses, score, bestCandidates);
      bool noCand = (startFromEmptyMLN_ && !existFormulas.empty()) ?
          (bestCandidates.empty() && highGainWtExistFormulas.empty()) 
        : bestCandidates.empty();
      if (noCand)
      { 
        cout << "Beam is empty. Ending search for MLN clauses." << endl; 
        printIterAndTimeElapsed(begSec);
        break;
      }

      useTightParams();

      if (reEvalBestCandsWithTightParams_) 
      {
        cout << "reevaluating top candidates... " << endl << endl;
        bsec = timer_.time();
        reEvaluateCandidates(bestCandidates, score);
        cout << "reevaluating top candidates took ";
        timer_.printTime(cout, timer_.time()-bsec); cout << endl << endl;
      }

        // effect the single best candidate on the MLN
      bsec = timer_.time();
      bool ok;
      if (startFromEmptyMLN_ && !existFormulas.empty())
        ok = effectBestCandidateOnMLNs(bestCandidates, existFormulas,
                                       highGainWtExistFormulas, score); 
      else
        ok = effectBestCandidateOnMLNs(bestCandidates, score);
      //score was updated when best candidate was effected on MLN
      cout << "effecting best candidates took ";
      timer_.printTime(cout, timer_.time()-bsec); cout << endl << endl;

      if (!ok) 
      { 
        cout << "failed to effect any of the best candidates on MLN" << endl
             << "stopping search for MLN clauses..." << endl;
        break;
      }

      printIterAndTimeElapsed(begSec);
    } // while (true)

    cout << "done searching for MLN clauses" << endl << endl;
    int numIterTaken = iter_+1;
    iter_= -1;

    useTightParams();
                                 
      //Prune each non-unit clause in MLN if it does not decrease score    
    cout << "pruning clauses from MLN..." << endl << endl;
    begSec = timer_.time();
    pruneMLN(score);
    cout << "pruning clauses from MLN took ";
    timer_.printTime(cout, timer_.time()-begSec); cout << endl << endl;

    printMLNClausesWithWeightsAndScore(score, -1);
    printMLNToFile(NULL, -2);
    cout << "num of iterations taken = " << numIterTaken << endl;

    cout << "time taken for structure learning = ";
    timer_.printTime(cout, timer_.time()-startSec_); cout << endl << endl;

    initialMLNClauses.deleteItemsAndClear();
    deleteExistFormulas(existFormulas);
  }//run()
  

  bool learnAndSetMLNWeights(double& score)
  {
    Array<double> priorMeans, priorStdDevs;
    double tmpScore = score;
    pllSetPriorMeansStdDevs(priorMeans, priorStdDevs, 0, NULL);
    int numClausesFormulas = getNumClausesFormulas();
    Array<double> wts;
    wts.growToSize(numClausesFormulas + 1);
      //indexTrans_ is already created
    int iter; bool error; double elapsedSec;
    tmpScore = maximizeScore(numClausesFormulas, 0, &wts, NULL, NULL,
                             iter, error, elapsedSec);
    if (error) 
    { 
      cout << "LBFGSB failed to find wts" << endl; 
      return false;
    }
    else 
    { 
      score = tmpScore;
      printNewScore((Clause*)NULL, NULL, iter, elapsedSec, score, 0, 0);
    }
    
    updateWts(wts, NULL, NULL);

    return true;
  }

  //////////////////////////// beam search ////////////////////////////////
 private:
    // priorMeans and priorStdDevs should be one larger than the current
    // num of clauses in the MLN. 
    // The clauses in initClauses must have their AuxClauseData created.
  void beamSearch(const Array<Clause*>& initClauses, 
                  const Array<Clause*>& initMLNClauses, const double& prevScore,
                  Array<Clause*>& bestClauses)
  {
    int iterBestGainUnchanged = 0; bsiter_ = -1;
    Array<double> priorMeans, priorStdDevs;
    ClauseOpHashArray* beam = new ClauseOpHashArray;
    for (int i = 0; i < initClauses.size(); i++) beam->append(initClauses[i]);

    int numClausesFormulas = getNumClausesFormulas();
    Array<double> wts;
    wts.growToSize(numClausesFormulas + 2); // first slot is unused
    Array<double> origWts(numClausesFormulas); // current clause/formula wts
    if (indexTrans_) indexTrans_->getClauseFormulaWts(origWts);
    else             mln0_->getClauseWts(origWts);

      //set the prior means and std deviations
    setPriorMeansStdDevs(priorMeans, priorStdDevs, 1, NULL);
    if (indexTrans_) indexTrans_->appendClauseIdxToClauseFormulaIdxs(1, 1);

    bool error;
    double begIterSec, bsec;
    Array<Clause*> candidates;
    while (!beam->empty() && iterBestGainUnchanged < bestGainUnchangedLimit_)
    {
      begIterSec = timer_.time();
      bsiter_++; 

      //if (bsiter_ == 1) break; // for testing purposes only
      cout << endl << "BEAM SEARCH ITERATION " << bsiter_ << endl << endl;

      cout << "creating candidate clauses..." << endl;
      candidates.clear();
      bsec = timer_.time();
      
      if (bsiter_==0) createCandidateClauses(beam, candidates, &initMLNClauses);
      else            createCandidateClauses(beam, candidates, NULL);
      
      cout << "num of candidates created = " << candidates.size() << "; took ";
      timer_.printTime(cout, timer_.time()-bsec); cout << endl << endl;

      cout << "evaluating gain of candidates..." << endl << endl;
      bsec = timer_.time();
      for (int i = 0; i < candidates.size(); i++)
        countAndMaxScoreEffectCandidate(candidates[i], &wts, &origWts,prevScore,
                                        true, priorMeans, priorStdDevs, error,
                                        NULL, NULL);
      cout << "evaluating gain of candidates took ";
      timer_.printTime(cout, timer_.time()-bsec); cout << endl << endl;


      cout << "finding best candidates..." << endl;
      bsec = timer_.time();
      ClauseOpHashArray* newBeam = new ClauseOpHashArray(beamSize_);
      bool newBestClause = findBestClausesAmongCandidates(candidates, newBeam, 
                                                          bestClauses);
      cout << "finding best candidates took ";
      timer_.printTime(cout, timer_.time()-bsec); cout << endl << endl;

      beam->deleteItemsAndClear();
      delete beam;
      beam = newBeam;

      if (newBestClause)
      {
        iterBestGainUnchanged = 0;
        cout << "found new best clause in beam search iter " << bsiter_ << endl;
      }
      else
        iterBestGainUnchanged++;


        //print the best clauses
      cout << "best clauses found in beam search iter " << bsiter_ << endl;
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
      for (int i = 0; i < bestClauses.size(); i++)
      {
        cout << i << "\t";
        bestClauses[i]->printWithoutWtWithStrVar(cout,(*domains_)[0]);
        cout << endl
             << "\tgain = " << bestClauses[i]->getAuxClauseData()->gain
             << ",  op = " 
             << Clause::getOpAsString(bestClauses[i]->getAuxClauseData()->op);
        if (bestClauses[i]->getAuxClauseData()->op != OP_REMOVE)
          cout << ",  wt = " << bestClauses[i]->getWt();
        cout << endl;
      }
      cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;


      cout << "BEAM SEARCH ITERATION " << bsiter_ << " took ";
      timer_.printTime(cout, timer_.time()-begIterSec); cout << endl << endl;
      cout << "Time elapsed = ";
      timer_.printTime(cout, timer_.time()-startSec_); cout << endl << endl;
    }

    if (beam->empty()) cout << "Beam search ended because beam is empty" 
                            << endl << endl;
    else cout << "Beam search ended because best clause did not change in "
              << iterBestGainUnchanged << " iterations" << endl << endl;
    
    beam->deleteItemsAndClear();
    delete beam;
    bsiter_ = -1;

    if (indexTrans_) indexTrans_->removeClauseIdxToClauseFormulaIdxs(1, 1);
  }//beamSearch()


  /////////////// functions dealing with counts and maximizing score ///////////
 private:
    //size of wts array must be numClauses+1, its first slot is unused
    //if indexTrans_ is used, the sizes of its instance variables must be 
    //correctly set before calling this function
  double maximizeScore(const int& numClausesFormulas, const int& numExtraWts,
                       Array<double>* const & wts, 
                       const Array<double>* const & origWts,
                       const Array<int>* const & removedClauseFormulaIdxs, 
                       int& iter, bool& error, double& elapsedSec)
  {
    if (origWts) 
    { for (int i=1 ; i<=numClausesFormulas; i++) (*wts)[i] = (*origWts)[i-1];}
    else         
    { for (int i=1; i<=numClausesFormulas; i++) (*wts)[i] = 0; }

    for (int i = 1; i <= numExtraWts; i++) (*wts)[numClausesFormulas+i] = 0;

    if (removedClauseFormulaIdxs) 
      for (int i = 0; i < removedClauseFormulaIdxs->size(); i++) 
        (*wts)[ (*removedClauseFormulaIdxs)[i]+1 ] = 0;

    //commented out: this is done before this function
    //if (indexTrans_) 
    //  indexTrans_->appendClauseIdxToClauseFormulaIdxs(numExtraWts, 1);
  
    double* wwts = (double*) wts->getItems();
    double begSec = timer_.time();
    double newScore 
      = lbfgsb_->minimize(numClausesFormulas + numExtraWts, wwts, iter, error);
    newScore = -newScore;
    elapsedSec = timer_.time() - begSec;

    //commented out: this is done after this function
    //if (indexTrans_) 
    //  indexTrans_->removeClauseIdxToClauseFormulaIdxs(numExtraWts, 1);

    return newScore;
  }


    //Compute/remove counts and maximize score when all the candidates are 
    //effected upon MLN. It computes/removes the counts in PseudoLogLikelihood 
    //and then finds the optimal pseudoLogLikelihood and weights, but it  
    //does not add/remove/replace clauses in the MLN.
    //if resetPriors is false, then priorMeans and priorStdDevs must 
    //already contain the mean/std dev of the existing MLN clauses, and be
    //of size #MLN clauses + #candidates
  double countAndMaxScoreEffectAllCandidates(const Array<Clause*>& candidates, 
                                             Array<double>* const & wts,
                                             const Array<double>*const& origWts,
                                             const double& prevScore,
                                             const bool& resetPriors,
                                             Array<double>& priorMeans,
                                             Array<double>& priorStdDevs,
                                             bool& error,
                                            Array<UndoInfo*>* const& uundoInfos,
                         Array<ClauseAndICDArray*>* const & appendedClauseInfos)
  {
    Array<UndoInfo*>* undoInfos = (uundoInfos)? uundoInfos:new Array<UndoInfo*>;
    Array<int> remClauseFormulaIdxs; //for domain 0 only
    Array<int*> idxPtrs; //for clean up
    Array<Clause*> toBeRemovedClauses;
    Array<Clause*> toBeAppendedClauses;
    int numClausesFormulas = getNumClausesFormulas();
      
    for (int i = 0; i < candidates.size(); i++)
    {
      Clause* cand = candidates[i];
      AuxClauseData* acd = cand->getAuxClauseData();
      int op = acd->op;

        //remove counts for candidate
      if (op == OP_REMOVE || op == OP_REPLACE || op == OP_REPLACE_ADDPRED || 
          op == OP_REPLACE_REMPRED)
      {
        Clause* remClause = (mlns_->size() > 1) ? 
          (Clause*) mln0_->getClause(acd->removedClauseIdx) : NULL;
        Array<int*> idxs; idxs.growToSize(mlns_->size());
        Array<Clause*> remClauses; remClauses.growToSize(mlns_->size());
        for (int d = 0; d < mlns_->size(); d++)
        {
          int remIdx;
          if (d == 0)
          {
            remIdx = acd->removedClauseIdx;
            remClauseFormulaIdxs.append(remIdx);
          }
          else
          {
            if (remClause->containsConstants())
              remClause->translateConstants((*domains_)[d-1], (*domains_)[d]);
            remIdx = (*mlns_)[d]->findClauseIdx(remClause);
          }
          idxs[d] = (*mlns_)[d]->getMLNClauseInfoIndexPtr(remIdx);
          remClauses[d] = (Clause*) (*mlns_)[d]->getClause(remIdx);
        }
        if (mlns_->size() > 1 && remClause->containsConstants())
          remClause->translateConstants((*domains_)[mlns_->size()-1], 
                                        (*domains_)[0]);
        pllRemoveCountsForClause(remClauses, idxs, undoInfos);

        if (op == OP_REMOVE) toBeRemovedClauses.append(cand);
      }

        //compute counts for candidate
      if (op == OP_ADD || op == OP_REPLACE || op == OP_REPLACE_ADDPRED || 
          op == OP_REPLACE_REMPRED)
      {
        Array<int*> tmpIdxs; tmpIdxs.growToSize(mlns_->size());
        for (int d = 0; d < mlns_->size(); d++)
        {
          int* tmpClauseIdxInMLN = new int( (*mlns_)[d]->getNumClauses() + 
                                            toBeAppendedClauses.size() );
          tmpIdxs[d] = tmpClauseIdxInMLN;
          idxPtrs.append(tmpClauseIdxInMLN);
        }

        Array<UndoInfo*>* tmpInfos = (appendedClauseInfos)?  
                                     new Array<UndoInfo*> : undoInfos;
        pllComputeCountsForClause(cand, tmpIdxs, tmpInfos); 
        toBeAppendedClauses.append(cand);
        
        if (appendedClauseInfos) 
        {
          undoInfos->append(*tmpInfos);
          ClauseAndICDArray* ca = new ClauseAndICDArray;
          ca->clause = cand;
          appendedClauseInfos->append(ca);          
          for (int j = 0; j < tmpInfos->size(); j++)
          {
            ca->icdArray.append(IndexCountDomainIdx());
            ca->icdArray.lastItem().iac=(*tmpInfos)[j]->affectedArr->lastItem();
            ca->icdArray.lastItem().domainIdx = (*tmpInfos)[j]->domainIdx;
          }
          delete tmpInfos;
        }
      }      
    }

    //find optimal weights and score

      // find the clause/formula indexes for indexes in remClauseFormulaIdxs 
    if (indexTrans_) 
      indexTrans_->getClauseFormulaIndexes(remClauseFormulaIdxs, 0);
    
    Array<double> removedValues;
    if (resetPriors) 
    {
      setPriorMeansStdDevs(priorMeans, priorStdDevs, toBeAppendedClauses.size(),
                           &remClauseFormulaIdxs);
      if (indexTrans_) indexTrans_->appendClauseIdxToClauseFormulaIdxs(
                                                 toBeAppendedClauses.size(), 1);
    }
    else
    {
      assert(priorMeans.size() == numClausesFormulas + candidates.size());
      assert(priorStdDevs.size() == numClausesFormulas + candidates.size());
      if (indexTrans_) 
        assert(indexTrans_->checkCIdxWtsGradsSize(candidates.size()));

      setRemoveAppendedPriorMeans(numClausesFormulas, priorMeans, 
                                  remClauseFormulaIdxs, 
                                  toBeAppendedClauses.size(), removedValues);
    }

    if (hasPrior_) 
      pll_->setMeansStdDevs(priorMeans.size(), priorMeans.getItems(), 
                            priorStdDevs.getItems());
    else
      pll_->setMeansStdDevs(-1, NULL, NULL);

    wts->growToSize(numClausesFormulas + toBeAppendedClauses.size() + 1);

    int iter; double elapsedSec; 
    double newScore = maximizeScore(numClausesFormulas, 
                                    toBeAppendedClauses.size(), 
                                    wts, origWts, &remClauseFormulaIdxs,
                                    iter, error, elapsedSec);

      //set the weights of the candidates
    for (int i = 0; i < toBeAppendedClauses.size(); i++)
      toBeAppendedClauses[i]->setWt( (*wts)[numClausesFormulas+i+1] );
      //set to large weight so that they are filtered out because of small wt
    for (int i = 0; i < toBeRemovedClauses.size(); i++)
      toBeRemovedClauses[i]->setWt(111111111); 
    
    Array<double> penalties;
    for (int i = 0; i < candidates.size(); i++)
      penalties.append(getPenalty(candidates[i]));
                     
    if (error) { newScore=prevScore; cout<<"LBFGSB failed to find wts"<<endl; }
    printNewScore(candidates, (*domains_)[0], iter, elapsedSec, 
                  newScore, newScore-prevScore, penalties);

    for (int i = 0; i < candidates.size(); i++)
      candidates[i]->getAuxClauseData()->gain = newScore-prevScore-penalties[i];

    if (uundoInfos == NULL) 
    { 
      pll_->undoAppendRemoveCounts(undoInfos); 
      delete undoInfos;
    }

    if (resetPriors) 
    {
      if (indexTrans_) indexTrans_->removeClauseIdxToClauseFormulaIdxs(
                                                 toBeAppendedClauses.size(), 1);
    }
    else
    {   //restore the prior means that were set to zero
      for (int i = 0; i < remClauseFormulaIdxs.size(); i++)
        priorMeans[remClauseFormulaIdxs[i]] = removedValues[i];
    }

    idxPtrs.deleteItemsAndClear();

    return newScore;
  } //countAndMaxScoreEffectAllCandidates()


  double countAndMaxScoreEffectCandidate(Clause* const & candidate, 
                                         Array<double>* const & wts,
                                         const Array<double>* const & origWts, 
                                         const double& prevScore,
                                         const bool& resetPriors,
                                         Array<double>& priorMeans,
                                         Array<double>& priorStdDevs,
                                         bool& error,
                                         Array<UndoInfo*>* const& undoInfos,
                         Array<ClauseAndICDArray*>* const & appendedClauseInfos)
  {
    Array<Clause*> candidates; candidates.append(candidate);
    return countAndMaxScoreEffectAllCandidates(candidates, wts, origWts, 
                                               prevScore, resetPriors,
                                               priorMeans, priorStdDevs, error,
                                               undoInfos, appendedClauseInfos);
  }


    //returns true if the candidates are all effected on MLN resulting in a
    //better score
  bool appendToAndRemoveFromMLNs(const Array<Clause*>& candidates,double& score,
                                const bool& makeChangeForEqualScore=false)
  {
    Array<double> wts;
    Array<UndoInfo*> undoInfos;
    Array<ClauseAndICDArray*> appendedClauseInfos;
    Array<double> priorMeans, priorStdDevs;
    bool error;

    double newScore
      = countAndMaxScoreEffectAllCandidates(candidates, &wts, NULL, score,
                                            true, priorMeans,priorStdDevs,error,
                                            &undoInfos, &appendedClauseInfos);

    bool improve = (makeChangeForEqualScore) ? (newScore >= score) 
                                             : (newScore > score);

    if (!error && improve)
    {
      score = newScore;

        //set MLN clauses/formulas to new weights
        //weights of appended clauses are set when they are added to the MLN 
      updateWts(wts, NULL, NULL);

      int numClausesFormulas = getNumClausesFormulas();

      for (int i = 0; i < appendedClauseInfos.size(); i++)
        appendedClauseInfos[i]->clause->setWt(wts[numClausesFormulas+i+1]);

        //effect candidates on MLN
      int appClauseIdx = 0;
      for (int i = 0; i < candidates.size(); i++)
      {
        Clause* cand = candidates[i];
        AuxClauseData* acd = cand->getAuxClauseData();
        int op = acd->op;

          // NOTE: cand was not canonicalized so as to easily compare it to the
          //       mln clause it is to replace. Canonicalize it now before 
          //       appending to MLN.
        if (op == OP_REPLACE) cand->canonicalize();

        if (op == OP_REMOVE || op == OP_REPLACE || op == OP_REPLACE_ADDPRED || 
            op == OP_REPLACE_REMPRED)
        {
          Clause* remClause = (Clause*) mln0_->getClause(acd->removedClauseIdx);

          for (int d = 0; d < mlns_->size(); d++)
          {
            if (d == 0)
            {
              Clause* r = removeClauseFromMLN(acd->removedClauseIdx, d);
              cout << "Modified MLN: Removed clause from MLN: "; 
              r->printWithoutWtWithStrVar(cout,(*domains_)[0]); cout << endl;
              if (op == OP_REMOVE && cand != r) delete cand;
              //delete r; //this is remClause which has to be used below
            }
            else
            {
              if (remClause->containsConstants())
                remClause->translateConstants((*domains_)[d-1], (*domains_)[d]);
              int remIdx = (*mlns_)[d]->findClauseIdx(remClause);
              delete removeClauseFromMLN(remIdx, d);
            }
          }
          delete remClause;
        }
        
        if (op == OP_ADD || op == OP_REPLACE || op == OP_REPLACE_ADDPRED || 
            op == OP_REPLACE_REMPRED)
        {
          Array<int*> idxPtrs; idxPtrs.growToSize(mlns_->size());
          for (int d = 0; d < mlns_->size(); d++)
          {
            Clause* c = cand;
            if (d > 0)
            {
              if (cand->containsConstants()) 
                cand->translateConstants((*domains_)[d-1], (*domains_)[d]);
              c = new Clause(*cand);
            }
            int idx = appendClauseToMLN(c, d);
            idxPtrs[d] = (*mlns_)[d]->getMLNClauseInfoIndexPtr(idx);
          }

          Array<IndexCountDomainIdx>& icds 
            = appendedClauseInfos[appClauseIdx++]->icdArray;
          for (int j = 0; j < icds.size(); j++)
            icds[j].iac->index = idxPtrs[ icds[j].domainIdx ];
            
          cout << "Modified MLN: Appended clause to MLN: "; 
          cand->printWithoutWtWithStrVar(cout,(*domains_)[0]); 
          cout << endl;
        }
      }

      assert(pll_->checkNoRepeatedIndex());
      assert(appClauseIdx == appendedClauseInfos.size());

      undoInfos.deleteItemsAndClear();

        //MLNs has changed, the index translation must be recreated
      if (indexTrans_) indexTrans_->createClauseIdxToClauseFormulaIdxsMap();
    }
    else
    {
      cout << "undoing candidates because score did not improve..."<<endl<<endl;
      pll_->undoAppendRemoveCounts(&undoInfos);
    }
    
    appendedClauseInfos.deleteItemsAndClear();
    return improve;
  }//appendToAndRemoveFromMLNs()


    //returns true if candidate is effected on MLN resulting in a better score
  bool appendToAndRemoveFromMLNs(Clause* const & candidate, double& score,
                                 const bool& makeChangeForEqualScore=false)
  {
    Array<Clause*> candidates; candidates.append(candidate);
    return appendToAndRemoveFromMLNs(candidates, score,makeChangeForEqualScore);
  }


    //Clause c enters & exits  function with its AuxClauseData's cache == NULL
  void pllComputeCountsForClause(Clause* const & c,
                                 const Array<int*>& clauseIdxInMLNs,
                                 Array<UndoInfo*>* const & undoInfos)
  {
    assert(c->getAuxClauseData()->cache == NULL);
    assert(clauseIdxInMLNs.size() == domains_->size());
    double begSec = timer_.time();
    
    if (cacheClauses_)    
    {
      int i;
      if ((i = cachedClauses_->find(c)) >= 0) // if clause and counts are cached
      {
        //cout << "using cached counts ";
        //cout << "for "; c->printWithoutWtWithStrVar(cout, (*domains_)[0]); 
        //cout << endl;
        pll_->insertCounts(clauseIdxInMLNs, undoInfos, 
                           (*cachedClauses_)[i]->getAuxClauseData()->cache);
        cout << "using cached counts took ";
        timer_.printTime(cout, timer_.time()-begSec); cout << endl;
        return;
      }
      else
      {
        assert(c->getAuxClauseData()->cache == NULL);
        if (cacheSizeMB_ <  maxCacheSizeMB_)
          c->newCache(domains_->size(), (*domains_)[0]->getNumPredicates());
        else
        {
          static bool printCacheFull = true;
          if (printCacheFull)
          {
            cout << "Cache is full, approximate size = " << cacheSizeMB_ 
                 << " MB" << endl;
            printCacheFull = false;
          }
        }
        //cout << "cache size (MB) = " << cacheSizeMB_ << endl;
      }
    }
    
    //we do not have the counts

      //if we don't need to translate constant ids from one domain to another
    if (!c->containsConstants())
    {
      for (int i = 0; i < domains_->size(); i++)
      {
        int* clauseIdxInMLN = clauseIdxInMLNs[i];
        pll_->computeCountsForNewAppendedClause(c, clauseIdxInMLN, i,
                                                undoInfos, sampleClauses_, 
                                                c->getAuxClauseData()->cache,
                                                NULL);
      }
    }
    else
    {  //we need to translate constant ids from one domain to another
      int i; 
      for (i = 0; i < domains_->size(); i++)
      {
        if (i > 0) c->translateConstants((*domains_)[i-1], (*domains_)[i]);
        int* clauseIdxInMLN = clauseIdxInMLNs[i];
        pll_->computeCountsForNewAppendedClause(c, clauseIdxInMLN, i,
                                                undoInfos, sampleClauses_, 
                                                c->getAuxClauseData()->cache,
                                                NULL);
      }
      if (i > 1) c->translateConstants((*domains_)[i-1], (*domains_)[0]);
    }


    if (c->getAuxClauseData()->cache)
    {
      if (cacheSizeMB_ < maxCacheSizeMB_)
      {
        cacheSizeMB_ += c->sizeMB();
        Array<Array<Array<CacheCount*>*>*>* cache =c->getAuxClauseData()->cache;
        c->getAuxClauseData()->cache = NULL;
        Clause* copyClause = new Clause(*c);
        copyClause->getAuxClauseData()->cache = cache;
        copyClause->compress();
        cachedClauses_->append(copyClause);
      }
      else
      {
        c->getAuxClauseData()->deleteCache();
        c->getAuxClauseData()->cache = NULL; 
      }
    }

    cout << "Computing counts took ";
    timer_.printTime(cout, timer_.time()-begSec); 
    //cout << " for " << endl ;
    //cout << "\t"; c->printWithoutWtWithStrVar(cout, (*domains_)[0]); 
    cout << endl;
  }


  void pllRemoveCountsForClause(const Array<Clause*>& remClauses, 
                                const Array<int*>& clauseIdxInMLNs,
                                Array<UndoInfo*>* const & undoInfos)
  { 
    assert(clauseIdxInMLNs.size() == domains_->size());
    double begSec = timer_.time();        
    for (int i = 0; i < domains_->size(); i++)
      pll_->removeCountsForClause(remClauses[i],clauseIdxInMLNs[i],i,undoInfos);
    cout << "Removing counts took ";
    timer_.printTime(cout, timer_.time()-begSec); 
    cout << endl;
  }



  void pllComputeCountsForInitialMLN()
  {
    for (int i = 0; i < mlns_->size(); i++)
    {
      cout << "computing counts for clauses in domain " << i << "..." << endl;
      MLN* mln = (*mlns_)[i];
      for (int j = 0; j < mln->getNumClauses(); j++)
      {
        Clause* c = (Clause*)mln->getClause(j);
        cout << "Clause " << j << ": ";
        c->printWithoutWtWithStrVar(cout, (*domains_)[i]); cout << endl;
        int* clauseIdxInMLN = mln->getMLNClauseInfoIndexPtr(j);
        pll_->computeCountsForNewAppendedClause(c, clauseIdxInMLN, i, NULL,
                                                sampleClauses_, NULL, NULL);
      }
    }
  }


  //////////////////// adding/removing/replacing MLN clauses /////////////////
 private:
  void addUnitClausesToMLNs()
  {
    Array<Predicate*> nonEvidPreds;
    for (int i = 0; i < preds_->size(); i++)
      if ((*areNonEvidPreds_)[(*preds_)[i]->getId()])
        nonEvidPreds.append((*preds_)[i]);

    Array<Clause*> unitClauses;
    bool allowEqualPreds = true;
    ClauseFactory::createUnitClauses(unitClauses,nonEvidPreds,allowEqualPreds);

    for (int i = 0; i < unitClauses.size(); i++)
    {
      if (mln0_->containsClause(unitClauses[i]))
      { delete unitClauses[i]; continue; }
      ostringstream oss; int idx; 
      unitClauses[i]->printWithoutWtWithStrVar(oss, (*domains_)[0]);
      
      for (int j = 0; j < mlns_->size(); j++)
      {
        Clause* c = (j == 0) ? unitClauses[i] : new Clause(*unitClauses[i]);
        (*mlns_)[j]->appendClause(oss.str(), false, c, priorMean_, false, idx,
                                  false, false, false, 0.0);
        ((MLNClauseInfo*)(*mlns_)[j]->getMLNClauseInfo(idx))->priorMean 
          = priorMean_;
      }
    }
  }


  void appendUnitClausesWithDiffCombOfVar(double& score)
  {
    bool allowEqualPreds = false;
    for (int i = 0; i < preds_->size(); i++)
    {
      if (!(*areNonEvidPreds_)[(*preds_)[i]->getId()]) continue;

      Clause* origClause = ClauseFactory::createUnitClause((*preds_)[i], 
                                                           allowEqualPreds);
      if (origClause == NULL) continue;
      assert(origClause->getAuxClauseData() == NULL);
      origClause->setAuxClauseData(new AuxClauseData);
      origClause->trackConstants();

      ClauseOpHashArray newUnitClauses;
      clauseFactory_->createUnitClausesWithDiffCombOfVar((*preds_)[i], OP_ADD,
                                                         -1, newUnitClauses);
      for (int j = 0; j < newUnitClauses.size(); j++)
      {
        Clause* newClause = newUnitClauses[j];

        if (origClause->same(newClause) || 
            !clauseFactory_->validClause(newClause) ||
            mln0_->containsClause(newClause))
        {
          newUnitClauses.removeItemFastDisorder(j);
          delete newClause; 
          j--;
          continue;
        }
          //score is updated
        if (!appendToAndRemoveFromMLNs(newClause, score)) delete newClause;
      }
      delete origClause;
    } // for each predicate
  }


  void flipMLNClauses(double& score)
  {
    Array<Clause*> mlnClauses;
    for (int i = 0; i < mln0_->getNumClauses(); i++)
    {
        //do not flip unit clauses or those derived from existential formulas
      if (mln0_->getClause(i)->getNumPredicates()==1 || !isModifiableClause(i))
        continue;
      mlnClauses.append((Clause*)mln0_->getClause(i));
    }

    for (int i = 0; i < mlnClauses.size(); i++)
    {
      Clause* origClause = mlnClauses[i];
      int origIdx = mln0_->findClauseIdx(origClause);

        // NOTE: new clauses are not canonicalized so that they can be easily 
        //       compared with the original to determine its penalty.
      bool canonicalizeNewClauses = false;
      ClauseOpHashArray newClauses;
      clauseFactory_->flipSensesInClause(origClause, OP_REPLACE, origIdx,
                                         newClauses, canonicalizeNewClauses);

      Array<double> priorMeans, priorStdDevs;
      setPriorMeansStdDevs(priorMeans, priorStdDevs, 1, NULL);
      Array<double> wts;
      wts.growToSize(getNumClausesFormulas()+2);
      Clause* bestClause = NULL;
      double bestScore = score;
      bool error;

      for (int j = 0; j < newClauses.size(); j++) // for each 'flipped' clause
      {
        Clause* newClause = newClauses[j];
        if (origClause->same(newClause) || newClause->hasRedundantPredicates() 
            || mln0_->containsClause(newClause)) { delete newClause; continue; }

        double newScore 
          = countAndMaxScoreEffectCandidate(newClause, &wts, NULL, score, true,
                                            priorMeans, priorStdDevs, error, 
                                            NULL, NULL);
        if (newScore > bestScore)
        {
          bestScore = newScore;
          if (bestClause) delete bestClause;
          bestClause = newClause;          
        }
        else
          delete newClause;
      }

      if (bestClause && !appendToAndRemoveFromMLNs(bestClause, score)) 
        delete bestClause;
    }// for each MLN clause
  }//flipMLNClauses()


    //NOTE: formulas with existentially quantified variables, or variables
    //      with mutually exclusive and exhaustive values are not pruned.
  void pruneMLN(double& score)
  {
    Array<Clause*> mlnClauses;
    for (int i = 0; i < mln0_->getNumClauses(); i++)
    {
        //do not prune unit clauses or clauses derived from existential formulas
      if (mln0_->getClause(i)->getNumPredicates() == 1 || 
          !isModifiableClause(i)) continue;
      mlnClauses.append((Clause*)mln0_->getClause(i));
    }
    
    for (int i = 0; i < mlnClauses.size(); i++)
    {
      Clause* origClause = mlnClauses[i];
      int origIdx = mln0_->findClauseIdx(origClause);
      Clause* copy = new Clause(*origClause);
      copy->setAuxClauseData(new AuxClauseData(0, OP_REMOVE, origIdx));
      copy->trackConstants();
      if (!appendToAndRemoveFromMLNs(copy, score, true)) delete copy;
    }
  }


  int appendClauseToMLN(Clause* const c, const int& domainIdx)
  {
    ostringstream oss; int idx;
    c->printWithoutWtWithStrVar(oss, (*domains_)[domainIdx]);
    MLN* mln = (*mlns_)[domainIdx];
    bool ok = mln->appendClause(oss.str(), false, c, c->getWt(), false, idx,
                                false, false, false, 0.0);
    if (!ok) { cout << "ERROR: failed to insert " << oss.str() <<" into MLN"
                    << endl; exit(-1); }
    ((MLNClauseInfo*)mln->getMLNClauseInfo(idx))->priorMean = priorMean_;
    return idx;
  }


  Clause* removeClauseFromMLN(const int& remClauseIdx, const int& domainIdx)
  {
    Clause* remClause = (*mlns_)[domainIdx]->removeClause(remClauseIdx);
    if (remClause == NULL)
    { 
      cout << "ERROR: failed to remove "; 
      remClause->printWithoutWtWithStrVar(cout, (*domains_)[0]);
      cout << " from MLN" << endl;
      exit(-1);
    }
    return remClause;
  }

  
  ////////////////////// functions for clause creation ////////////////////////
 private:
  void addPredicateToClause(Clause* const & beamClause, const int& op, 
                            const int& removeClauseIdx,
                            ClauseOpHashArray& newClauses)
  {
      // create new clauses by adding predicates to beamClause
    clauseFactory_->addPredicateToClause(*preds_,beamClause,op,removeClauseIdx,
                                         true,newClauses,false);
      // if this is a unit clause, flip its sense and create new clauses by
      // adding preds to it
    if (beamClause->getNumPredicates() == 1)
    {
      beamClause->getPredicate(0)->invertSense();        
      clauseFactory_->addPredicateToClause(*preds_, beamClause, op,
                                           removeClauseIdx, true, newClauses,
                                           false);
      beamClause->getPredicate(0)->invertSense();
    }    
  }


    // create new clauses by removing a predicate from beamClause
  void removePredicateFromClause(Clause* const & beamClause, const int& op,
                                 const int& removeClauseIdx,
                                 ClauseOpHashArray& newClauses)
  {
    if (beamClause->getNumPredicates() > 2)
      clauseFactory_->removePredicateFromClause(beamClause, op, removeClauseIdx,
                                                newClauses);
  }


  void addNewClauseToCandidates(Clause* const & newClause,
                                Array<Clause*>& candidates, 
                                Array<Clause*>* const & thrownOut)
  {
      // remove any new clauses that are already in MLN
    if (thrownOut && mln0_->containsClause(newClause)) 
    { thrownOut->append(newClause); return;}
    newClause->setWt(0);
    candidates.append(newClause);
  }


  void createCandidateClauses(const ClauseOpHashArray* const & beam, 
                              Array<Clause*>& candidates,
                              const Array<Clause*>* const & initMLNClauses)
  {
      // new clauses that are thrown out because they are in MLN
    Array<Clause*> thrownOutClauses; 
    ClauseOpHashArray newClauses;
    for (int i = 0; i < beam->size(); i++)
    {
      Clause* beamClause = (*beam)[i];
      AuxClauseData* beamacd = beamClause->getAuxClauseData();      
      int op = beamacd->op;
      int newClausesBegIdx = newClauses.size();

      if (op == OP_ADD)
      {
        int remIdx = beamacd->removedClauseIdx;
        addPredicateToClause(beamClause, op, remIdx, newClauses);
        removePredicateFromClause(beamClause, op, remIdx, newClauses);
        for (int j = newClausesBegIdx; j < newClauses.size(); j++)
          addNewClauseToCandidates(newClauses[j], candidates,&thrownOutClauses);
      }
      else
      if (op == OP_REPLACE_ADDPRED)
      {
        addPredicateToClause(beamClause, op, beamacd->removedClauseIdx,
                             newClauses);
        for (int j = newClausesBegIdx; j < newClauses.size(); j++)
          addNewClauseToCandidates(newClauses[j], candidates,&thrownOutClauses);
      }
      else
      if (op == OP_REPLACE_REMPRED)
      {
        removePredicateFromClause(beamClause, op, beamacd->removedClauseIdx,
                                  newClauses);
        for (int j = newClausesBegIdx; j < newClauses.size(); j++)
          addNewClauseToCandidates(newClauses[j], candidates,&thrownOutClauses);
      }
      else
      if (op == OP_NONE)
      {
        int idx = mln0_->findClauseIdx(beamClause);
        bool beamClauseInMLN = (idx >= 0);
        bool removeBeamClause = (beamClauseInMLN && 
                                 beamClause->getNumPredicates() > 1);
        bool isModClause = (!beamClauseInMLN || isModifiableClause(idx));

        if (isModClause)
        {
          addPredicateToClause(beamClause, OP_ADD, -1, newClauses);
          for (int j = newClausesBegIdx; j < newClauses.size(); j++)
            addNewClauseToCandidates(newClauses[j], candidates,
                                     &thrownOutClauses);

          if (removeBeamClause)
          {
            newClausesBegIdx = newClauses.size();
            addPredicateToClause(beamClause, OP_REPLACE_ADDPRED, idx, 
                                 newClauses);
            for (int j = newClausesBegIdx; j < newClauses.size(); j++)
              addNewClauseToCandidates(newClauses[j], candidates,
                                       &thrownOutClauses);            
          }

          newClausesBegIdx = newClauses.size();
          removePredicateFromClause(beamClause, OP_ADD, -1, newClauses);
          for (int j = newClausesBegIdx; j < newClauses.size(); j++)
            addNewClauseToCandidates(newClauses[j], candidates,
                                     &thrownOutClauses);          

          if (removeBeamClause)
          {
            newClausesBegIdx = newClauses.size();          
            removePredicateFromClause(beamClause, OP_REPLACE_REMPRED, idx,
                                      newClauses);
            for (int j = newClausesBegIdx; j < newClauses.size(); j++)
              addNewClauseToCandidates(newClauses[j], candidates,
                                       &thrownOutClauses);          
          }

          if (removeBeamClause)
          {
            Clause* c = new Clause(*beamClause);
            c->getAuxClauseData()->op = OP_REMOVE;
            c->getAuxClauseData()->removedClauseIdx = idx;
            addNewClauseToCandidates(c, candidates, NULL);
          }
        }        
      }
      else
        assert(op == OP_REMOVE || op == OP_REPLACE);
    } // for each item in beam


    if (initMLNClauses)
    {
        // add the MLN clauses to the candidates in the first step of beam 
        // search when we start from an empty MLN
      int newClausesBegIdx = newClauses.size();
      for (int i = 0; i < initMLNClauses->size(); i++)
      {
        Clause* c = new Clause(*((*initMLNClauses)[i]));
        if (newClauses.append(c) < 0) delete c;
      }
      for (int i = newClausesBegIdx; i < newClauses.size(); i++)
        addNewClauseToCandidates(newClauses[i], candidates, &thrownOutClauses);
    }

    for (int i = 0; i < thrownOutClauses.size(); i++) 
      delete thrownOutClauses[i];
  }


  //////////////////////// helper functions ////////////////////////////////

 private:
  void useTightParams()
  {
    sampleClauses_ = false;
    pll_->setSampleGndPreds(false);
    lbfgsb_->setMaxIter(lbMaxIter_);
    lbfgsb_->setFtol(lbConvThresh_);
    cacheClauses_ = false;
  }


  void useLooseParams()
  {
    sampleClauses_ = origSampleClauses_;
    pll_->setSampleGndPreds(sampleGndPreds_);
    if (looseMaxIter_ >= 0) lbfgsb_->setMaxIter(looseMaxIter_);
    if (looseConvThresh_ >= 0) lbfgsb_->setFtol(looseConvThresh_);
    cacheClauses_ = origCacheClauses_;
  }

  
    //i is the clause's index in mln_[0]
  bool isModifiableClause(const int& i) const
  { return (!mln0_->isExistClause(i) && !mln0_->isExistUniqueClause(i)); }


  bool isNonModifiableFormula(const FormulaAndClauses* const & fnc) const
  { return (fnc->hasExist || fnc->isExistUnique); }


  void printIterAndTimeElapsed(const double& begSec)
  {
    cout << "Iteration " << iter_ << " took ";
    timer_.printTime(cout, timer_.time()-begSec); cout << endl << endl;
    cout << "Time elapsed = ";
    timer_.printTime(cout, timer_.time()-startSec_); cout << endl << endl;
  }


  void removeClausesFromMLNs()
  {
    for (int i = 0; i < mlns_->size(); i++)
      (*mlns_)[i]->removeAllClauses(NULL);
  }


  void reEvaluateCandidates(Array<Clause*>& candidates,
                            const double& prevScore)
  {
    int numClausesFormulas = getNumClausesFormulas();
    Array<double> wts; 
    wts.growToSize(numClausesFormulas + 2);//slot 0 is unused

    Array<double> priorMeans, priorStdDevs;
    setPriorMeansStdDevs(priorMeans, priorStdDevs, 1, NULL);
    if (indexTrans_) indexTrans_->appendClauseIdxToClauseFormulaIdxs(1, 1);

    bool error;
    for (int i = 0; i < candidates.size(); i++)
    {
      candidates[i]->getAuxClauseData()->gain = 0;
      countAndMaxScoreEffectCandidate(candidates[i], &wts, NULL, prevScore,
                                      true, priorMeans, priorStdDevs, error,
                                      NULL, NULL);
    }

    if (indexTrans_) indexTrans_->removeClauseIdxToClauseFormulaIdxs(1, 1);

    Array<Clause*> tmpCand(candidates);
    rquicksort(tmpCand);
    candidates.clear();
    for (int i = 0; i < tmpCand.size(); i++)
    {
      if (tmpCand[i]->getAuxClauseData()->gain > minGain_ && 
          fabs(tmpCand[i]->getWt()) >= minWt_)
        candidates.append(tmpCand[i]);
      else
        delete tmpCand[i];
    }

    cout << "reevaluated top " << candidates.size() 
         << " candidates with tight params:" << endl;
    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    for (int i = 0; i < candidates.size(); i++)
    {
      cout << i << "\t";
      candidates[i]->printWithoutWtWithStrVar(cout,(*domains_)[0]);
      cout << endl 
           << "\tgain = " << candidates[i]->getAuxClauseData()->gain 
           << ",  wt = " << candidates[i]->getWt()
           << ",  op = "
           << Clause::getOpAsString(candidates[i]->getAuxClauseData()->op) 
           << endl;
    }
    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl << endl;;
  }


    //Copy the best beamSize_ candidates into beam, and place the best 
    //numEstBestClauses_ clauses among those in candidates and bestClauses into
    //bestClauses. Returns true if the best clause changed.
  bool findBestClausesAmongCandidates(Array<Clause*>& candidates,
                                      ClauseOpHashArray* const & beam,
                                      Array<Clause*>& bestClauses)
  {
    assert(beam->size() == 0);

    // fill the beam with best beamSize_ candidates

      //sort candidates in order of decreasing gain
    rquicksort(candidates); 

    for (int i = 0; i < candidates.size(); i++)
    {
      if (beam->size() >= beamSize_) break;
      double candGain = candidates[i]->getAuxClauseData()->gain;
      double candAbsWt = fabs(candidates[i]->getWt());
      if (candGain > minGain_ && candAbsWt >= minWt_)
      {
        int a = beam->append(new Clause(*(candidates[i])));
        assert(a >= 0); a = 0; //avoid compilation warning
      }
    }
    
    if (beam->size() == 0)
    {
      for (int i = 0; i < candidates.size(); i++) delete candidates[i];
      return false;
    }


    Clause* prevBestClause = (bestClauses.size() > 0) ? 
                             new Clause(*(bestClauses[0])) : NULL;

    // pick the best numEstBestClauses_ clauses

    ClauseOpSet cset;
    ClauseOpSet::iterator it;
    for (int i = 0; i < candidates.size(); i++) cset.insert(candidates[i]);

    for (int i = 0; i < bestClauses.size(); i++) 
    {
      if ((it=cset.find(bestClauses[i])) == cset.end())
      {
        candidates.append(bestClauses[i]);
        cset.insert(bestClauses[i]);
      }
      else
      {
        assert((*it)->getAuxClauseData()->gain == 
               bestClauses[i]->getAuxClauseData()->gain);
        delete bestClauses[i];
      }
    }
    
    rquicksort(candidates);

    bestClauses.clear();
    for (int i = 0; i < candidates.size(); i++)
    {
      if (bestClauses.size() < numEstBestClauses_)
      {
        double candGain = candidates[i]->getAuxClauseData()->gain;
        double candAbsWt = fabs(candidates[i]->getWt());
        if (candGain > minGain_ && candAbsWt >= minWt_) 
          bestClauses.append(candidates[i]);
        else
        {
          //cout << "\tDiscarding candidate because of low gain or abs weight:";
          //cout << " gain = " << candGain << ", wt = " 
          //     << candidates[i]->getWt() << endl;
          //cout << "\t"; 
          //candidates[i]->printWithoutWtWithStrVar(cout, (*domains_)[0]);
          //cout << endl;
          delete candidates[i];
        }
      }
      else
        delete candidates[i];
    }

      //check whether the best clause has changed
    bool bestClauseChanged;
    if (bestClauses.size() > 0)
    {
      if (prevBestClause == NULL) bestClauseChanged = true;
      else
      {
          // if the prev and cur best clauses are the same 
        if (bestClauses[0]->same(prevBestClause) && 
            bestClauses[0]->getAuxClauseData()->op == 
            prevBestClause->getAuxClauseData()->op) 
          bestClauseChanged =  false;
        else
        { // prev and cur best clauses are different
            //if they have the same gain
          if (bestClauses[0]->getAuxClauseData()->gain >
              prevBestClause->getAuxClauseData()->gain) 
            bestClauseChanged = true;
          else
            bestClauseChanged = false;
        }
      }
    }
    else
      bestClauseChanged = false;

    if (prevBestClause) delete prevBestClause;
    return bestClauseChanged;
  }


  bool effectBestCandidateOnMLNs(Clause* const & cand, double& score)
  {
    if (appendToAndRemoveFromMLNs(cand, score))
    { 
      printMLNClausesWithWeightsAndScore(score, iter_); 
      //printMLNToFile(NULL, iter_);
      return true;
    }
    return false;
  }


  bool effectBestCandidateOnMLNs(Array<Clause*>& bestCandidates, double& score)
  {
    cout << "effecting best candidate on MLN..." << endl << endl;
    bool ok = false;
    int i;
    for (i = 0; i < bestCandidates.size(); i++)
    {
      cout << "effecting best candidate " << i << " on MLN..." << endl << endl;
      if ((ok = effectBestCandidateOnMLNs(bestCandidates[i], score))) break;
      cout << "failed to effect candidate on MLN." << endl;
      delete bestCandidates[i];
    }
    for (int j = i+1; j < bestCandidates.size();j++) delete bestCandidates[j];
    return ok;
  }


  double getPenalty(const Clause* const & cand)
  {
    int op = cand->getAuxClauseData()->op;

    if (op == OP_ADD) return cand->getNumPredicates() * penalty_;

    if (op == OP_REPLACE_ADDPRED) 
    {
      int remIdx = cand->getAuxClauseData()->removedClauseIdx;
      int origlen = mln0_->getClause(remIdx)->getNumPredicates();        
      int diff = cand->getNumPredicates() - origlen;
      assert(diff > 0);
      return diff * penalty_;
    }

    if (op == OP_REPLACE_REMPRED) 
    {
      int remIdx = cand->getAuxClauseData()->removedClauseIdx;
      int origlen = mln0_->getClause(remIdx)->getNumPredicates();        
      int diff = origlen - cand->getNumPredicates();
      assert(diff > 0);
      return diff * penalty_;
    }

    if (op == OP_REPLACE)
    {
        //NOTE: op is REPLACE only when an MLN clause is to be replaced with
        //      one that is identical except for its predicates' senses
      int remIdx = cand->getAuxClauseData()->removedClauseIdx;
      const Clause* mlnClause = mln0_->getClause(remIdx);
      assert(cand->getNumPredicates() == mlnClause->getNumPredicates());
      int diff = 0;
      for (int i = 0; i < cand->getNumPredicates(); i++)
      {
        Predicate* cpred = (Predicate*) cand->getPredicate(i);
        Predicate* mpred = (Predicate*) mlnClause->getPredicate(i);
        assert(cpred->same(mpred));
        if (cpred->getSense() != mpred->getSense()) diff++;
      }
      return diff * penalty_;      
    }

    if (op == OP_REMOVE) return cand->getNumPredicates() * penalty_;

    assert(false);
    return 88888;
  }


  void printNewScore(const Array<Clause*>& carray, const Domain* const & domain,
                     const int& lbfgsbIter, const double& lbfgsbSec,
                     const double& newScore, const double& gain,
                     const Array<double>& penalties)
  {
    cout << "*************************** " << candCnt_++ << ", iter " << iter_ 
         << ", beam search iter " << bsiter_ << endl;
    for (int i = 0; i < carray.size(); i++)
    {
      if (carray[i]) 
      { 
        cout << "candidate     : ";
        carray[i]->printWithoutWtWithStrVar(cout,domain); cout << endl;
        cout << "op            : ";
        cout << Clause::getOpAsString(carray[i]->getAuxClauseData()->op) <<endl;
        cout << "removed clause: ";
        int remIdx = carray[i]->getAuxClauseData()->removedClauseIdx;
        if (remIdx < 0) cout << "NULL";
        else { mln0_->getClause(remIdx)->printWithoutWtWithStrVar(cout,domain);}
        cout << endl;
        if (carray[i]->getAuxClauseData()->prevClauseStr.length() > 0)
        {
          cout << "prevClause    : ";
          cout << carray[i]->getAuxClauseData()->prevClauseStr << endl;
        }
        if (carray[i]->getAuxClauseData()->addedPredStr.length() > 0)
        {
          cout << "addedPred     : ";
          cout << carray[i]->getAuxClauseData()->addedPredStr << endl;
        }
        if (carray[i]->getAuxClauseData()->removedPredIdx >= 0)
        {
          cout << "removedPredIdx: ";
          cout << carray[i]->getAuxClauseData()->removedPredIdx << endl;
        }

        cout << "score    : " << newScore << endl;
        cout << "gain     : " << gain << endl;
        cout << "penalty  : " << penalties[i] << endl;
        cout << "net gain : " << gain - penalties[i] << endl;
        if (carray[i]->getAuxClauseData()->op != OP_REMOVE)
          cout << "wt       : " << carray[i]->getWt() << endl;
      }
    }
    cout << "num LBFGSB iter      = " << lbfgsbIter << endl;
    cout << "time taken by LBFGSB = "; Timer::printTime(cout, lbfgsbSec);
    cout << endl;
    cout << "*************************** " << endl << endl;;
  }


  void printNewScore(const Clause* const & c, const Domain* const & domain, 
                     const int& lbfgsbIter, const double& lbfgsbSec,
                     const double& newScore, const double& gain,
                     const double& penalty)
  {
    Array<Clause*> carray;
    Array<double> parray;
    if (c) { carray.append((Clause*)c); parray.append(penalty); }
    printNewScore(carray, domain, lbfgsbIter, lbfgsbSec, newScore, gain,parray);
  }


  void printMLNClausesWithWeightsAndScore(const double& score, const int& iter)
  {    
    if (iter >= 0) cout << "MLN in iteration " << iter << ":" << endl;
    else           cout << "MLN:" << endl;
    cout << "------------------------------------" << endl;
    if (indexTrans_) indexTrans_->printClauseFormulaWts(cout, true);
    else             mln0_->printMLNClausesFormulas(cout, (*domains_)[0], true);
    cout << "------------------------------------" << endl;
    cout << "score = "<< score << endl << endl;
  }


  void printMLNToFile(const char* const & appendStr, const int& iter)
  {
    string fname = outMLNFileName_;

    if (appendStr) fname.append(".").append(appendStr);

    if (iter >= -1) 
    { 
      char buf[100]; 
      sprintf(buf, "%d", iter);
      fname.append(".iter").append(buf);
    }

    if (appendStr || iter >= -1) fname.append(".mln");
    
    ofstream out(fname.c_str());
    if (!out.good()) { cout << "ERROR: failed to open " <<fname<<endl;exit(-1);}

      // output the predicate declaration
    out << "//predicate declarations" << endl;
    (*domains_)[0]->printPredicateTemplates(out);
    out << endl;

      // output the function declarations
    if ((*domains_)[0]->getNumFunctions() > 0) 
    {
      out << "//function declarations" << endl;
      (*domains_)[0]->printFunctionTemplates(out);
      out << endl;
    }

    if (indexTrans_) indexTrans_->printClauseFormulaWts(out, false);
    else             mln0_->printMLNClausesFormulas(out, (*domains_)[0], false);

    out << endl;
    out.close();
  } 
  

  void printClausesInBeam(const ClauseOpHashArray* const & beam)
  {
    cout.setf(ios_base::left, ios_base::adjustfield);
    cout << "^^^^^^^^^^^^^^^^^^^ beam ^^^^^^^^^^^^^^^^^^^" << endl;
    for (int i = 0; i < beam->size(); i++)
    {
      cout << i << ":  ";
      cout.width(10); cout << (*beam)[i]->getWt(); cout.width(0); cout << " ";
      (*beam)[i]->printWithoutWtWithStrVar(cout, (*domains_)[0]); cout << endl;
    }
    cout.width(0);
    cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
    
  }

  void setPriorMeansStdDevs(Array<double>& priorMeans, 
                            Array<double>& priorStdDevs, const int& addSlots,
                            const Array<int>* const& removedSlotIndexes)
  {
    priorMeans.clear(); priorStdDevs.clear();

    if (indexTrans_)
    {
      indexTrans_->setPriorMeans(priorMeans);
      priorStdDevs.growToSize(priorMeans.size());
      for (int i = 0; i < priorMeans.size(); i++)
        priorStdDevs[i] = priorStdDev_;
    }
    else
    {
      for (int i = 0; i < mln0_->getNumClauses(); i++)
      {
        priorMeans.append(mln0_->getMLNClauseInfo(i)->priorMean);
        priorStdDevs.append(priorStdDev_);
      }
    }

    if (removedSlotIndexes)
    {
      for (int i = 0; i < removedSlotIndexes->size(); i++)
        priorMeans[ (*removedSlotIndexes)[i] ] = 0;
    }

    for (int i = 0; i < addSlots; i++)
    {
      priorMeans.append(priorMean_);
      priorStdDevs.append(priorStdDev_);      
    }
  }


  void pllSetPriorMeansStdDevs(Array<double>& priorMeans, 
                               Array<double>& priorStdDevs, const int& addSlots,
                               const Array<int>* const & removedSlotIndexes)
  {
    if (hasPrior_)
    {
      setPriorMeansStdDevs(priorMeans, priorStdDevs, 
                           addSlots, removedSlotIndexes);
      pll_->setMeansStdDevs(priorMeans.size(), priorMeans.getItems(), 
                            priorStdDevs.getItems());
    }
    else
      pll_->setMeansStdDevs(-1, NULL, NULL);
  }


  void setRemoveAppendedPriorMeans(const int& numClausesFormulas,
                                   Array<double>& priorMeans, 
                                   const Array<int>& removedSlotIndexes, 
                                   const int& addSlots,
                                   Array<double>& removedValues)
  {
    for (int i = 0; i < removedSlotIndexes.size(); i++)
    {
      removedValues.append(priorMeans[ removedSlotIndexes[i] ]);
      priorMeans[ removedSlotIndexes[i] ] = 0;
    }

    for (int i = 0; i < addSlots; i++)
      priorMeans[numClausesFormulas+i] = priorMean_;
  }

  
  int getNumClausesFormulas()
  {
    if (indexTrans_) return indexTrans_->getNumClausesAndExistFormulas();
    return mln0_->getNumClauses();
  }


  void updateWts(const Array<double>& wts,
                 const Array<Clause*>* const & appendedClauses,
                 const Array<string>* const & appendedFormulas)
  {
    if (indexTrans_)
    {
      Array<double> tmpWts; 
      tmpWts.growToSize(wts.size()-1);
      for (int i = 1; i < wts.size(); i++) tmpWts[i-1] = wts[i];
      indexTrans_->updateClauseFormulaWtsInMLNs(tmpWts, appendedClauses,
                                                appendedFormulas);
    }
    else
    {
      for (int i = 0; i < mlns_->size(); i++)
        for (int j = 0; j < (*mlns_)[i]->getNumClauses(); j++)          
          ((Clause*)(*mlns_)[i]->getClause(j))->setWt(wts[j+1]);
    }
  }


  /////////////////////////// quicksort ////////////////////////////////

 private:
    //sort clauses in decreasing order of gain
  void rquicksort(Array<Clause*>& clauses, const int& l, const int& r)  
  {
    Clause** items = (Clause**) clauses.getItems();
    if (l >= r) return;
    Clause* tmp = items[l];
    items[l] = items[(l+r)/2];
    items[(l+r)/2] = tmp;

    int last = l;
    for (int i = l+1; i <= r; i++)
      if (items[i]->getAuxClauseData()->gain > 
          items[l]->getAuxClauseData()->gain)
      {
        ++last;
        Clause* tmp = items[last];
        items[last] = items[i];
        items[i] = tmp;
      }
    
    tmp = items[l];
    items[l] = items[last];
    items[last] = tmp;
    rquicksort(clauses, l, last-1);
    rquicksort(clauses, last+1, r); 
  }

  void rquicksort(Array<Clause*>& ca) { rquicksort(ca, 0, ca.size()-1); }



  ////////////////// functions to handle existential formulas ///////////////
 private:

  void deleteExistFormulas(Array<ExistFormula*>& existFormulas)
  { existFormulas.deleteItemsAndClear(); }


  void getMLNClauses(Array<Clause*>& initialMLNClauses,
                     Array<ExistFormula*>& existFormulas)
  {
    for (int i = 0; i < mln0_->getNumClauses(); i++)
      if (isModifiableClause(i)) 
      {
        Clause* c = (Clause*) mln0_->getClause(i);
        initialMLNClauses.append(new Clause(*c));
      }

    const FormulaAndClausesArray* fnca = mln0_->getFormulaAndClausesArray();
    for (int i = 0; i < fnca->size(); i++)
      if (isNonModifiableFormula((*fnca)[i]))
        existFormulas.append(new ExistFormula((*fnca)[i]->formula));

    for (int i = 0; i < existFormulas.size(); i++)
    {
      string formula = existFormulas[i]->formula;      
      FormulaAndClauses tmp(formula, 0, false, false, false);
      const FormulaAndClausesArray* fnca
        = (*mlns_)[0]->getFormulaAndClausesArray();
      int a = fnca->find(&tmp);            
      existFormulas[i]->numPreds = (*fnca)[a]->numPreds;

      Array<Array<Clause*> >& cnfClausesForDomains 
        = existFormulas[i]->cnfClausesForDomains;
      cnfClausesForDomains.growToSize(mlns_->size());
      for (int j = 0; j < mlns_->size(); j++)
      {
        Array<Clause*>& cnfClauses = cnfClausesForDomains[j];
        fnca = (*mlns_)[j]->getFormulaAndClausesArray();
        a = fnca->find(&tmp);
        assert(a >= 0);
        IndexClauseHashArray* indexClauses = (*fnca)[a]->indexClauses;
        for (int k = 0; k < indexClauses->size(); k++)
        {
          Clause* c = new Clause(*((*indexClauses)[k]->clause));
          c->newAuxClauseData();
          cnfClauses.append(c);
        }
        cnfClauses.compress();        
      }
    }
  }


    //cnfClauses enter & exit function with its AuxClauseData's cache == NULL
  inline void pllCountsForExistFormula(Clause* cnfClause, 
                                       const int& domainIdx,
                                       int* clauseIdxInMln,
                                       Array<UndoInfo*>* const & undoInfos)
  {
    assert(cnfClause->getAuxClauseData()->cache == NULL);
    bool inCache = false; 
    bool hasDomainCounts = false; 
    double prevCNFClauseSize = 0;

    if (cacheClauses_)
    {
      int i; 
      if ((i = cachedClauses_->find(cnfClause)) >= 0) //if clause is in cache
      {
        inCache = true;
        Array<Array<Array<CacheCount*>*>*>* cache =
          (*cachedClauses_)[i]->getAuxClauseData()->cache;
        Array<Array<CacheCount*>*>*  domainCache = (*cache)[domainIdx];
        for (int j = 0; j < domainCache->size(); j++)
          if ((*domainCache)[j] != NULL) { hasDomainCounts = true; break; }
      }
     
      if (hasDomainCounts) // if clause and counts for domain are cached
      {
        pll_->insertCounts(clauseIdxInMln, undoInfos,
                           (*cachedClauses_)[i]->getAuxClauseData()->cache,
                           domainIdx);
        return;
      }
      else
      {
        if (cacheSizeMB_ <  maxCacheSizeMB_)
        {
          if (inCache) //if clause is in cache but not counts for domain
          {
            cnfClause->getAuxClauseData()->cache 
              = (*cachedClauses_)[i]->getAuxClauseData()->cache;
            prevCNFClauseSize = cnfClause->sizeMB();
          }
          else
            cnfClause->newCache(domains_->size(),
                                (*domains_)[0]->getNumPredicates());
        }
        else
        {
          static bool printCacheFull = true;
          if (printCacheFull)
          {
            cout << "Cache is full, approximate size = " << cacheSizeMB_ 
                 <<" MB" << endl;
            printCacheFull = false;
          }
        }
      }
    }
    
      //we do not have the counts
    pll_->computeCountsForNewAppendedClause(cnfClause,clauseIdxInMln,domainIdx,
                                            undoInfos, sampleClauses_, 
                                            cnfClause->getAuxClauseData()->cache,
                                            NULL);

    if (cnfClause->getAuxClauseData()->cache)
    {
      if (inCache)
      {
        cacheSizeMB_ += cnfClause->sizeMB() - prevCNFClauseSize;
        cnfClause->getAuxClauseData()->cache = NULL; 
      }
      else
      {
        if (cacheSizeMB_ < maxCacheSizeMB_)
        {
          cacheSizeMB_ += cnfClause->sizeMB();
          Array<Array<Array<CacheCount*>*>*>* cache 
            = cnfClause->getAuxClauseData()->cache;
          cnfClause->getAuxClauseData()->cache = NULL;
          Clause* copyClause = new Clause(*cnfClause);
          copyClause->getAuxClauseData()->cache = cache;
          copyClause->compress();
          cachedClauses_->append(copyClause);
        }
        else
        {
          cnfClause->getAuxClauseData()->deleteCache();
          cnfClause->getAuxClauseData()->cache = NULL; 
        }
      }
    }
  }


  inline void printNewScore(const string existFormStr,const int& lbfgsbIter,
                            const double& lbfgsbSec, const double& newScore,
                            const double& gain, const double& penalty, 
                            const double& wt)
  {
    cout << "*************************** " << candCnt_++ << ", iter " << iter_ 
         << ", beam search iter " << bsiter_ << endl;

    cout << "exist formula : " << existFormStr << endl;
    cout << "op            : OP_ADD" << endl; 
    cout << "score    : " << newScore << endl;
    cout << "gain     : " << gain << endl;
    cout << "penalty  : " << penalty << endl;
    cout << "net gain : " << gain - penalty << endl;
    cout << "wt       : " << wt << endl;
    cout << "num LBFGSB iter      = " << lbfgsbIter << endl;
    cout << "time taken by LBFGSB = "; Timer::printTime(cout, lbfgsbSec);
    cout << endl;
    cout << "*************************** " << endl << endl;;
  }


  inline void rquicksort(Array<ExistFormula*>& existFormulas, const int& l, 
                         const int& r)  
  {
    ExistFormula** items = (ExistFormula**) existFormulas.getItems();
    if (l >= r) return;
    ExistFormula* tmp = items[l];
    items[l] = items[(l+r)/2];
    items[(l+r)/2] = tmp;

    int last = l;
    for (int i = l+1; i <= r; i++)
      if (items[i]->gain > items[l]->gain)
      {
        ++last;
        ExistFormula* tmp = items[last];
        items[last] = items[i];
        items[i] = tmp;
      }
    
    tmp = items[l];
    items[l] = items[last];
    items[last] = tmp;
    rquicksort(existFormulas, l, last-1);
    rquicksort(existFormulas, last+1, r); 
  }


  inline void rquicksort(Array<ExistFormula*>& ef) 
  { rquicksort(ef, 0, ef.size()-1); }


    // Compute counts for formula, and learn optimal weights and gain.
    // If we are computing counts only, then the counts are added permanently into
    // PseudoLogLikelihood.
  inline void evaluateExistFormula(ExistFormula* const& ef, 
                                   const bool& computeCountsOnly,
              Array<Array<Array<IndexAndCount*> > >* const & iacArraysPerDomain,
                                   const double& prevScore)
  {
    bool undo = !computeCountsOnly;
    bool evalGainLearnWts = !computeCountsOnly;

    Array<int*> idxPtrs; //for clean up
    Array<UndoInfo*>* undoInfos = (undo) ? new Array<UndoInfo*> : NULL;

      //compute counts for the CNF clauses

    Array<Array<Clause*> >& cnfClausesForDomains = ef->cnfClausesForDomains;
    for (int d = 0; d < cnfClausesForDomains.size(); d++)
    {
      Array<Clause*>& cnfClauses = cnfClausesForDomains[d];    
      MLN* mln = (*mlns_)[d];
      for (int k = 0; k < cnfClauses.size(); k++)
      {
        Clause* c = cnfClauses[k];
          //if we adding counts permanently, then don't add counts that were added
          //previously; otherwise we can add the duplicate counts because they 
          //will be undone later.
        if (!undo && mln->containsClause(c)) continue;
        int* tmpClauseIdxInMLN = new int(mln->getNumClauses() + k);
        idxPtrs.append(tmpClauseIdxInMLN);
        if (iacArraysPerDomain)
        {
          Array<Array<IndexAndCount*> >& iacArrays = (*iacArraysPerDomain)[d];
          assert(iacArrays.size() == cnfClauses.size());
          Array<IndexAndCount*>& iacArray = iacArrays[k];

          Array<UndoInfo*> tmpUndoInfos;
          pllCountsForExistFormula(c, d, tmpClauseIdxInMLN, &tmpUndoInfos);
          for (int i = 0; i < tmpUndoInfos.size(); i++)
            iacArray.append(tmpUndoInfos[i]->affectedArr->lastItem());
          if (undoInfos) undoInfos->append(tmpUndoInfos);
          else           tmpUndoInfos.deleteItemsAndClear();
        }
        else
          pllCountsForExistFormula(c, d, tmpClauseIdxInMLN, undoInfos);
      }
    }

      //find optimal weights and gain of formula 

    if (evalGainLearnWts)
    {
      Array<double> priorMeans, priorStdDevs;    
          //either add one existentially quantified formula (when clauses do not
          //line up perfectly across databases) or add all CNF clauses (when
          //they do)
      int numAdded = (indexTrans_) ? 1 : cnfClausesForDomains[0].size();
      setPriorMeansStdDevs(priorMeans, priorStdDevs, numAdded, NULL);

      if (hasPrior_) 
        pll_->setMeansStdDevs(priorMeans.size(), priorMeans.getItems(), 
                              priorStdDevs.getItems());
      else           
        pll_->setMeansStdDevs(-1, NULL, NULL);

      if (indexTrans_)
      {
        for (int d = 0; d < cnfClausesForDomains.size(); d++)
        {
          int numCNFClauses = cnfClausesForDomains[d].size();
          indexTrans_->appendClauseIdxToClauseFormulaIdxs(1, numCNFClauses);
        }
      }
      
      int iter; bool error; double elapsedSec; 
      int numClausesFormulas = getNumClausesFormulas();
      Array<double>* wts = &(ef->wts);
      wts->growToSize(numClausesFormulas + numAdded + 1);

      double newScore = maximizeScore(numClausesFormulas, numAdded, wts, 
                                      NULL, NULL, iter, error, elapsedSec);

      if (indexTrans_)
      {
        for (int d = 0; d < cnfClausesForDomains.size(); d++)
        {
          int numCNFClauses = cnfClausesForDomains[d].size();
          indexTrans_->removeClauseIdxToClauseFormulaIdxs(1, numCNFClauses);
        }
      }

      double totalWt = 0;
      for (int i = 0; i < numAdded; i++)
        totalWt += (*wts)[numClausesFormulas+i+1];

      double penalty = ef->numPreds * penalty_;
      
      ef->gain = newScore-prevScore-penalty;
      ef->wt = totalWt;
      ef->newScore = newScore;

      if (error) {newScore=prevScore;cout<<"LBFGSB failed to find wts"<<endl;}
      printNewScore(ef->formula, iter, elapsedSec, newScore, 
                    newScore-prevScore, penalty, totalWt);
    }
    
    if (undo) { pll_->undoAppendRemoveCounts(undoInfos); delete undoInfos; }
    idxPtrs.deleteItemsAndClear();
  }//evaluateExistFormula()



  double evaluateExistFormulas(Array<ExistFormula*>& existFormulas, 
                               Array<ExistFormula*>& highGainWtFormulas,
                               const double& prevScore)
  {
    if (existFormulas.size() == 0) return 0;
    for (int i = 0; i < existFormulas.size(); i++)
      evaluateExistFormula(existFormulas[i], false, NULL, prevScore);
    
    Array<ExistFormula*> tmp(existFormulas);
    rquicksort(tmp);
    double minGain = DBL_MAX;
    cout << "evaluated existential formulas " << endl;
    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;  
    for (int i = 0; i < tmp.size(); i++)
    {
      existFormulas[i] = tmp[i];
      double gain = existFormulas[i]->gain;
      double wt = existFormulas[i]->wt;
      cout << i << "\t" << existFormulas[i]->formula << endl
           << "\tgain = " << gain << ",  wt = " << wt << ",  op = OP_ADD" 
           << endl;
      if (gain > minGain_ && wt >= minWt_) 
      {
        highGainWtFormulas.append(tmp[i]);
        if (gain < minGain)  minGain = gain;
      }
    }
    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl << endl;

    if (minGain == DBL_MAX) minGain = 0;
    return minGain;
  } 


  inline void appendExistFormulaToMLNs(ExistFormula* const & ef)
  {
    Array<Array<Array<IndexAndCount*> > > iacsPerDomain;
    Array<Array<Clause*> >& cnfClausesForDomains = ef->cnfClausesForDomains;
    iacsPerDomain.growToSize(cnfClausesForDomains.size());
    for (int d = 0; d < cnfClausesForDomains.size(); d++)
    {
      Array<Array<IndexAndCount*> >& iacArrays = iacsPerDomain[d];
      iacArrays.growToSize( cnfClausesForDomains[d].size() );
    }
    evaluateExistFormula(ef, true, &iacsPerDomain, 0);

    Array<double>& wts = ef->wts;
    int numClausesFormulas = getNumClausesFormulas();

      //update weights before CNF clauses are added as duplicate clause wts are
      //accumulated as they are added.
    if (indexTrans_ == NULL) updateWts(wts, NULL, NULL);

      //append CNF clauses to MLN
    for (int d = 0; d < cnfClausesForDomains.size(); d++)
    {
      MLN* mln = (*mlns_)[d];
      Array<Array<IndexAndCount*> >& iacArrays = iacsPerDomain[d];
      Array<Clause*>& cnfClauses = cnfClausesForDomains[d];

      for (int i = 0; i < cnfClauses.size(); i++)
      {
        int idx;
          //when we are learning the weight of a formula (i.e. indexsTrans !=NULL)
          //its CNF clause weight don't matter, so set them to 0.
        double wt = (indexTrans_) ?  0 : wts[numClausesFormulas+i+1];
          //if cnfClauses[i] is already in MLN, its weights will be correctly set
          //in updateWts() later
        mln->appendClause(ef->formula, true, new Clause(*cnfClauses[i]),
                          wt, false, idx, false, false, false, 0.0);
        mln->setFormulaPriorMean(ef->formula, priorMean_);
        ((MLNClauseInfo*)mln->getMLNClauseInfo(idx))->priorMean 
          += priorMean_/cnfClauses.size();

        int* idxPtr = mln->getMLNClauseInfoIndexPtr(idx);
        Array<IndexAndCount*>& iacs = iacArrays[i]; 
        for (int j = 0; j < iacs.size(); j++)  iacs[j]->index = idxPtr;
      }
    }
    assert(pll_->checkNoRepeatedIndex());


    if (indexTrans_)
    {
         //update weights after the formula has been added to the MLN
      Array<string> appendedFormula; 
      appendedFormula.append(ef->formula);
      updateWts(wts, NULL, &appendedFormula);
        //MLNs has changed, the index translation must be recreated
      indexTrans_->createClauseIdxToClauseFormulaIdxsMap();
    }

    cout << "Modified MLN: Appended formula to MLN: " << ef->formula << endl;
  }


  inline bool effectExistFormulaOnMLNs(ExistFormula* ef, 
                                       Array<ExistFormula*>& existFormulas, 
                                       double& score)
  {
    cout << "effecting existentially quantified formula " << ef->formula 
         << " on MLN..." << endl;
    appendExistFormulaToMLNs(ef);
    score = ef->newScore;
    printMLNClausesWithWeightsAndScore(score, iter_); 
    printMLNToFile(NULL, iter_);

    int r = existFormulas.find(ef);
    assert(r >= 0);
    ExistFormula* rf = existFormulas.removeItemFastDisorder(r);
    assert(rf == ef);
    delete rf;
    return true;
  }


  bool effectBestCandidateOnMLNs(Array<Clause*>& bestCandidates, 
                                 Array<ExistFormula*>& existFormulas,
                                 Array<ExistFormula*>& highGainWtFormulas,
                                 double& score)
  {
    cout << "effecting best candidate among existential formulas and "
         << "best candidates on MLN..." << endl << endl;

    int a = 0, b = 0;
    bool ok = false;
    int numCands = bestCandidates.size() + highGainWtFormulas.size();
    for (int i = 0; i < numCands; i++)
    {
      if (a >= bestCandidates.size())
      {
        if ((ok = effectExistFormulaOnMLNs(highGainWtFormulas[b++],
                                           existFormulas, score))) break;
      }
      else
      if (b >= highGainWtFormulas.size())
      {
        cout << "effecting best candidate " << a << " on MLN..." << endl;
        if ((ok = effectBestCandidateOnMLNs(bestCandidates[a++], score))) break;
        cout << "failed to effect candidate on MLN." << endl;
        delete bestCandidates[a-1];
      }
      else
      if (highGainWtFormulas[b]->gain > 
          bestCandidates[a]->getAuxClauseData()->gain)
      {
        if ((ok = effectExistFormulaOnMLNs(highGainWtFormulas[b++], 
                                           existFormulas, score))) break;
      }
      else
      {
        cout << "effecting best candidate " << a << " on MLN..." << endl;
        if ((ok = effectBestCandidateOnMLNs(bestCandidates[a++], score))) break;
        cout << "failed to effect candidate on MLN." << endl;
        delete bestCandidates[a-1];
      }
    }
  
    for (int i = a; i < bestCandidates.size(); i++) delete bestCandidates[i];
    return ok;
  }

 /*
  void deleteExistFormulas(Array<ExistFormula*>& existFormulas);

  void getMLNClauses(Array<Clause*>& initialMLNClauses,
                     Array<ExistFormula*>& existFormulas);

    //cnfClauses enter & exit function with its AuxClauseData's cache == NULL
  inline void pllCountsForExistFormula(Clause* cnfClause, const int& domainIdx,
                                       int* clauseIdxInMLN,
                                       Array<UndoInfo*>* const & undoInfos);

  inline void printNewScore(const string existFormStr, const int& lbfgsbIter, 
                            const double& lbfgsbSec, const double& newScore, 
                            const double& gain, const double& penalty, 
                            const double& wt);

  inline void rquicksort(Array<ExistFormula*>& existFormulas, const int& l, 
                         const int& r);

  inline void rquicksort(Array<ExistFormula*>& ef);
 
  inline void evaluateExistFormula(ExistFormula* const& ef,
                                   const bool& computeCountsOnly,
              Array<Array<Array<IndexAndCount*> > >* const & iacArraysPerDomain,
                                   const double& prevScore);

  double evaluateExistFormulas(Array<ExistFormula*>& existFormulas, 
                               Array<ExistFormula*>& highGainWtFormulas,
                               const double& prevScore);

  inline void appendExistFormulaToMLNs(ExistFormula* const & ef);


  inline bool effectExistFormulaOnMLNs(ExistFormula* ef, 
                                       Array<ExistFormula*>& existFormulas, 
                                       double& score);

  bool effectBestCandidateOnMLNs(Array<Clause*>& bestCandidates, 
                                 Array<ExistFormula*>& existFormulas,
                                 Array<ExistFormula*>& highGainWtFormulas,
                                 double& score);

  void updateWeightsSGD(int step);
*/
  /////////////////////////////////////////////////////////////////////////
 private:
    // mln0_ and mlns_ are not owned by StructLearn; do not delete;
    // if there are more than one domains, mln0_ corresponds to domains_[0]
  MLN* mln0_; 
  Array<MLN*>* mlns_;
  bool startFromEmptyMLN_;
  string outMLNFileName_;
  Array<Domain*>* domains_; // not owned by StructLearn; do not delete;
  Array<Predicate*>* preds_;
  Array<bool>* areNonEvidPreds_;
  ClauseFactory* clauseFactory_;

  bool cacheClauses_;
  bool origCacheClauses_;
  ClauseHashArray* cachedClauses_;
  double cacheSizeMB_;
  double maxCacheSizeMB_;

  bool tryAllFlips_;
  bool sampleClauses_;
  bool origSampleClauses_;

  PseudoLogLikelihood* pll_;
  bool hasPrior_;
  double priorMean_;
  double priorStdDev_;
  bool wtPredsEqually_;

  LBFGSB* lbfgsb_;
  int lbMaxIter_; 
  double lbConvThresh_;
  int looseMaxIter_; 
  double looseConvThresh_;

  int beamSize_;
  int bestGainUnchangedLimit_;
  int numEstBestClauses_;
  double minGain_;
  double minWt_;
  double penalty_;

  bool sampleGndPreds_;
  double fraction_;
  int minGndPredSamples_;
  int maxGndPredSamples_;

  bool reEvalBestCandsWithTightParams_;

  Timer timer_;
  int candCnt_;

  int iter_;
  int bsiter_;
  double startSec_;

    //Used to map clause indexes in MLNs to the weights indexes that are 
    //presented to LBFGSB. This is required if the clauses do not line up
    //across databases (as when an existentially quantified formula has a 
    //different number of formulas in its CNF for different databases).
  IndexTranslator* indexTrans_;

  bool structGradDescent_;
  bool withEM_;
};


#endif
