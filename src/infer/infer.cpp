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

#include <unistd.h>
#include <fstream>
#include <climits>
#include <sys/times.h>
#include "fol.h"
#include "arguments.h"
#include "util.h"
#include "infer.h"

extern const char* ZZ_TMP_FILE_POSTFIX; //defined in fol.y

  // TODO: List the arguments common to learnwts and inference in
  // inferenceargs.h. This can't be done with a static array.
ARGS ARGS::Args[] = 
{
    // BEGIN: Common arguments
  ARGS("i", ARGS::Req, ainMLNFiles, 
       "Comma-separated input .mln files."),

  ARGS("cw", ARGS::Opt, aClosedWorldPredsStr,
       "Specified non-evidence atoms (comma-separated with no space) are "
       "closed world, otherwise, all non-evidence atoms are open world. Atoms "
       "appearing here cannot be query atoms and cannot appear in the -o "
       "option."),

  ARGS("ow", ARGS::Opt, aOpenWorldPredsStr,
       "Specified evidence atoms (comma-separated with no space) are open "
       "world, while other evidence atoms are closed-world. "
       "Atoms appearing here cannot appear in the -c option."),
    // END: Common arguments

  ARGS("queryEvidence", ARGS::Tog, aisQueryEvidence, 
       "If this flag is set, then all the groundings of query preds not in db "
       "are assumed false evidence."),

    // BEGIN: Common inference arguments
  ARGS("m", ARGS::Tog, amapPos, 
       "Run MAP inference and return only positive query atoms."),

  ARGS("a", ARGS::Tog, amapAll, 
       "Run MAP inference and show 0/1 results for all query atoms."),

  ARGS("p", ARGS::Tog, agibbsInfer, 
       "Run inference using MCMC (Gibbs sampling) and return probabilities "
       "for all query atoms."),
  
  ARGS("ms", ARGS::Tog, amcsatInfer,
       "Run inference using MC-SAT and return probabilities "
       "for all query atoms"),

  ARGS("bp", ARGS::Tog, abpInfer,
       "Run inference using belief propagation and return probabilities "
       "for all query atoms"),

  ARGS("efbp", ARGS::Tog, aefbpInfer,
       "Run inference using expanding frontier belief propagation and return "
       "probabilities for all query atoms"),

  ARGS("decision", ARGS::Tog, adecisionInfer,
       "Run decision inference using BP (or EFBP) and return "
       "max. utility assignment of action atoms"),

  ARGS("simtp", ARGS::Tog, asimtpInfer,
       "Run inference using simulated tempering and return probabilities "
       "for all query atoms"),

  ARGS("outputNetwork", ARGS::Tog, aoutputNetwork,
       "Build the network and output to results file, instead of "
       "running inference"),

  ARGS("counts", ARGS::Tog, aclauseCounts,
       "Write clause counts, not atom values or probabilities"),

  ARGS("seed", ARGS::Opt, aSeed,
       "[2350877] Seed used to initialize the randomizer in the inference "
       "algorithm. If not set, seed is initialized from a fixed random number."),

  ARGS("lazy", ARGS::Opt, aLazy, 
       "[false] Run lazy version of inference if this flag is set."),
  
  ARGS("lazyNoApprox", ARGS::Opt, aLazyNoApprox, 
       "[false] Lazy version of inference will not approximate by deactivating "
       "atoms to save memory. This flag is ignored if -lazy is not set."),
  
  ARGS("memLimit", ARGS::Opt, aMemLimit, 
       "[-1] Maximum limit in kbytes which should be used for inference. "
       "-1 means main memory available on system is used."),

  ARGS("PrintSamplePerIteration", ARGS::Opt, aPrintSamplePerIteration,
       "Whether to print out variable values at each HMCS sample round."),
	
  ARGS("SAInterval", ARGS::Opt, saInterval, "SA interval"),

  ARGS("MaxSeconds", ARGS::Opt, aMaxSeconds, "Max seconds for HMWS and SA."),

  ARGS("SATempDownRatio", ARGS::Opt, aSATempDownRatio, "Simulated annealing "
       "temperature degrade ratio."),

  ARGS("SA", ARGS::Opt, aSA, "simulated annealing inference."),

  ARGS("hybrid", ARGS::Opt, aHybrid, 
	   "Flag for HMLN inference."),

  ARGS("propstdev", ARGS::Opt, aProposalStdev, 
	   "[1.0]Proposal stdev for SA step in HybridSAT."),

  ARGS("contSamples", ARGS::Opt, aContSamples, 
	   "output file for continuous variable samples."),

    // END: Common inference arguments

    // BEGIN: MaxWalkSat args
  ARGS("mwsMaxSteps", ARGS::Opt, amwsMaxSteps,
       "[100000] (MaxWalkSat) The max number of steps taken."),

  ARGS("tries", ARGS::Opt, amwsTries, 
       "[1] (MaxWalkSat) The max number of attempts taken to find a solution."),

  ARGS("targetWt", ARGS::Opt, amwsTargetWt,
       "[the best possible] (MaxWalkSat) MaxWalkSat tries to find a solution "
       "with weight <= specified weight."),

  ARGS("breakHardClauses", ARGS::Tog, amwsHard, 
       "[false] (MaxWalkSat) If true, MaxWalkSat can break a hard clause in "
       "order to satisfy a soft one."),
  
  ARGS("heuristic", ARGS::Opt, amwsHeuristic,
       "[2] (MaxWalkSat) Heuristic used in MaxWalkSat (0 = RANDOM, 1 = BEST, "
       "2 = TABU, 3 = SAMPLESAT)."),
  
  ARGS("tabuLength", ARGS::Opt, amwsTabuLength,
       "[5] (MaxWalkSat) Minimum number of flips between flipping the same "
       "atom when using the tabu heuristic in MaxWalkSat." ),

  ARGS("lazyLowState", ARGS::Opt, amwsLazyLowState, 
       "[false] (MaxWalkSat) If false, the naive way of saving low states "
       "(each time a low state is found, the whole state is saved) is used; "
       "otherwise, a list of variables flipped since the last low state is "
       "kept and the low state is reconstructed. This can be much faster for "
       "very large data sets."),  
    // END: MaxWalkSat args

    // BEGIN: MCMC args
  ARGS("burnMinSteps", ARGS::Opt, amcmcBurnMinSteps,
       "[100] (MCMC) Minimun number of burn in steps (-1: no minimum)."),

  ARGS("burnMaxSteps", ARGS::Opt, amcmcBurnMaxSteps,
       "[100] (MCMC) Maximum number of burn-in steps (-1: no maximum)."),

  ARGS("minSteps", ARGS::Opt, amcmcMinSteps, 
       "[-1] (MCMC) Minimum number of MCMC sampling steps."),

  ARGS("maxSteps", ARGS::Opt, amcmcMaxSteps, 
       "[1000] (MCMC) Maximum number of MCMC sampling steps."),

  ARGS("maxSeconds", ARGS::Opt, amcmcMaxSeconds, 
       "[-1] (MCMC) Max number of seconds to run MCMC (-1: no maximum)."),
    // END: MCMC args
  
    // BEGIN: Simulated tempering args
  ARGS("subInterval", ARGS::Opt, asimtpSubInterval,
        "[2] (Simulated Tempering) Selection interval between swap attempts"),

  ARGS("numRuns", ARGS::Opt, asimtpNumST,
        "[3] (Simulated Tempering) Number of simulated tempering runs"),

  ARGS("numSwap", ARGS::Opt, asimtpNumSwap,
        "[10] (Simulated Tempering) Number of swapping chains"),
    // END: Simulated tempering args

    // BEGIN: BP args
  ARGS("lifted", ARGS::Tog, aliftedInfer, 
       "[false] If true, lifted inference is run"),
  
  ARGS("useHC", ARGS::Tog, auseHC, 
       "[false] If true (and lifted is also true), use the hypercube representation"),
  
  ARGS("useCT", ARGS::Tog, auseCT, 
       "[false] If true (and lifted and hc is also true), use the constraints for hypercube " 
	   "representation"),

  ARGS("convThresh", ARGS::Opt, abpConvergenceThresh,
        "[1e-4] (BP) Max difference in probabilities to determine convergence"),

  ARGS("convIterations", ARGS::Opt, abpConvergeRequiredItrCnt,
        "[20] (BP) Number of converging iterations to determine convergence"),

  ARGS("explicitRep", ARGS::Tog, aexplicitRep, 
       "[false] If true, explicit representation type is used in lifted "
       "inference; otherwise, implicit representation type is used"), 
  
  ARGS("hcCreateType", ARGS::Opt, ahcCreateType, 
       "[Basic] Type of method used for creating hypercubes. DT/BAsic"),
  
  ARGS("hcCreateNoise", ARGS::Opt, ahcCreateNoise, 
       "[0.0] Amount of noise while creating hypercubes"),
  
  ARGS("lncIter", ARGS::Opt, alncIter, 
       "[0] Number of LNC Iterations (0 means run till end)"),
  
  ARGS("noHC", ARGS::Opt, anoHCPredsStr,
       "Comma separated list of predicates, for which hyperCubes should not be created i.e. "
	   " HyperCubes are not created for the specified predicates i.e. each ground tuple "
	   " is a hypercube"),
  
  // END: BP args

    // BEGIN: SampleSat args
  ARGS("numSolutions", ARGS::Opt, amwsNumSolutions,
       "[10] (MC-SAT) Return nth SAT solution in SampleSat"),

  ARGS("saRatio", ARGS::Opt, assSaRatio,
       "[0] (MC-SAT) Ratio of sim. annealing steps mixed with WalkSAT in "
       "MC-SAT"),

  ARGS("saTemperature", ARGS::Opt, assSaTemp,
        "[10] (MC-SAT) Temperature (/100) for sim. annealing step in "
        "SampleSat"),

  ARGS("lateSa", ARGS::Tog, assLateSa,
       "[true] Run simulated annealing from the start in SampleSat"),
    // END: SampleSat args

    // BEGIN: Gibbs sampling args
  ARGS("numChains", ARGS::Opt, amcmcNumChains, 
       "[10] (Gibbs) Number of MCMC chains for Gibbs sampling (there must be "
       "at least 2)."),

  ARGS("delta", ARGS::Opt, agibbsDelta,
       "[0.05] (Gibbs) During Gibbs sampling, probabilty that epsilon error is "
       "exceeded is less than this value."),

  ARGS("epsilonError", ARGS::Opt, agibbsEpsilonError,
       "[0.01] (Gibbs) Fractional error from true probability."),

  ARGS("fracConverged", ARGS::Opt, agibbsFracConverged, 
       "[0.95] (Gibbs) Fraction of ground atoms with probabilities that "
       "have converged."),

  ARGS("walksatType", ARGS::Opt, agibbsWalksatType, 
       "[1] (Gibbs) Use Max Walksat to initialize ground atoms' truth values "
       "in Gibbs sampling (1: use Max Walksat, 0: random initialization)."),

  ARGS("testConvergence", ARGS::Opt, agibbsTestConvergence, 
       "[false] Perform convergence test for Gibbs sampling."),

  ARGS("samplesPerTest", ARGS::Opt, agibbsSamplesPerTest, 
       "[100] Perform convergence test once after this many number of samples "
       "per chain."),
    // END: Gibbs sampling args

    // BEGIN: Args specific to stand-alone inference
  ARGS("e", ARGS::Req, aevidenceFiles, 
       "Comma-separated .db files containing known ground atoms (evidence), "
       "including function definitions."),

  ARGS("r", ARGS::Req, aresultsFile,
       "The probability estimates are written to this file."),
    
  ARGS("q", ARGS::Opt, aqueryPredsStr, 
       "Query atoms (comma-separated with no space)  "
       ",e.g., cancer,smokes(x),friends(Stan,x). Query atoms are always "
       "open world."),

  ARGS("f", ARGS::Opt, aqueryFile,
       "A .db file containing ground query atoms, "
       "which are are always open world."),
    // END: Args specific to stand-alone inference

  ARGS()
};


void printResults(const string& queryFile, const string& queryPredsStr,
				  Domain *domain, ostream& out, 
				  GroundPredicateHashArray* const &queries,
				  Inference* const &inference, HVariableState* const &state)
{
    // Lazy version: Have to generate the queries from the file or query string.
    // This involves calling createQueryFilePreds / createComLineQueryPreds
  if (aLazy)
  {
    const GroundPredicateHashArray* gndPredHashArray = NULL;
    Array<double>* gndPredProbs = NULL;
      // Inference algorithms with probs: have to retrieve this info from state.
      // These are the ground preds which have been brought into memory. All
      // others have always been false throughout sampling.
    if (!(amapPos || amapAll))
    {
      gndPredHashArray = state->getGndPredHashArrayPtr();
      gndPredProbs = new Array<double>;
      gndPredProbs->growToSize(gndPredHashArray->size());
      for (int i = 0; i < gndPredProbs->size(); i++)
        (*gndPredProbs)[i] = inference->getProbabilityH((*gndPredHashArray)[i]);
    }

    if (queryFile.length() > 0)
    {
      cout << "Writing query predicates that are specified in query file..."
           << endl;
      bool ok = createQueryFilePreds(queryFile, domain, domain->getDB(),
                                     NULL, NULL, true, out, amapPos,
                                     gndPredHashArray, gndPredProbs, NULL);
      if (!ok)
      {
        cout <<"Failed to create query predicates."<< endl;
        exit(-1);
      }
    }

    Array<int> allPredGndingsAreQueries;
    allPredGndingsAreQueries.growToSize(domain->getNumPredicates(), false);
    if (queryPredsStr.length() > 0)
    {
      cout << "Writing query predicates that are specified on command line..." 
           << endl;
      bool ok = createComLineQueryPreds(queryPredsStr, domain,
                                        domain->getDB(), NULL, NULL,
                                        &allPredGndingsAreQueries, true,
                                        out, amapPos, gndPredHashArray,
                                        gndPredProbs, NULL);
      if (!ok)
      {
        cout <<"Failed to create query predicates."<< endl; exit(-1);
      }
    }

    if (!(amapPos || amapAll))
      delete gndPredProbs;
  }
    // Eager version: Queries have already been generated and we can get the
    // information directly from the state
  else
  {
    if (amapPos)
      inference->printTruePredsH(out);
    else
    {
      for (int i = 0; i < queries->size(); i++)
      {
          // Prob is smoothed in inference->getProbability
        double prob = inference->getProbabilityH((*queries)[i]);
        (*queries)[i]->print(out, domain); out << " " << prob << endl;
      }
    }
  }
}


/**
 * Prints the results of inference to a stream.
 * 
 * @param queryFile File name containing the query predicates. This is only
 * used with lazy inference. If empty, this is not used.
 * @param query String of query predicates separated by commas without spaces.
 * This is only used with lazy inference. If empty, this is not used.
 * @param domain Domain in which the predicates exist.
 * @param out Stream to which the results are printed.
 * @param queries Queries already built from query file or string. This is only
 * used with eager inference.
 * @param inference Inference algorithm used which contains the results.
 * @param state VariableState used by the inference algorithm which contains
 * the results.
 */
void printResults(const string& queryFile, const string& queryPredsStr,
                  Domain *domain, ostream& out, 
                  GroundPredicateHashArray* const &queries,
                  Inference* const &inference, VariableState* const &state,
                  Array<Predicate*> const &queryPreds,
                  Array<TruthValue> const &queryPredValues)
{
    // Lazy version: Have to generate the queries from the file or query string.
    // This involves calling createQueryFilePreds / createComLineQueryPreds
  if (aLazy)
  {
      // Inference algorithms with probs: have to retrieve this info from state.
      // These are the ground preds which have been brought into memory. All
      // others have always been false throughout sampling.

    for (int i = 0; i < queryPreds.size(); i++)
    {
      int val = (queryPredValues[i] == TRUE) ? 1 : 0;

        // Prob is smoothed in inference->getProbability
      double prob = inference->getProbability((GroundPredicate*)queryPreds[i]);

      queryPreds[i]->print(out, domain);
      out << " " << prob << " " << val << endl;
    }

	/*
	if (!(amapPos || amapAll))
    {
      gndPredHashArray = state->getGndPredHashArrayPtr();
      gndPredProbs = new Array<double>;
      gndPredProbs->growToSize(gndPredHashArray->size());
      for (int i = 0; i < gndPredProbs->size(); i++)
        (*gndPredProbs)[i] =
          inference->getProbability((*gndPredHashArray)[i]);
    }
    
    if (queryFile.length() > 0)
    {
      cout << "Writing query predicates that are specified in query file..."
           << endl;
      bool ok = createQueryFilePreds(queryFile, domain, domain->getDB(), NULL,
                                     NULL, true, out, amapPos,
                                     gndPredHashArray, gndPredProbs);
      if (!ok) { cout <<"Failed to create query predicates."<< endl; exit(-1); }
    }

    Array<int> allPredGndingsAreQueries;
    allPredGndingsAreQueries.growToSize(domain->getNumPredicates(), false);
    if (queryPredsStr.length() > 0)
    {
      cout << "Writing query predicates that are specified on command line..." 
           << endl;
      bool ok = createComLineQueryPreds(queryPredsStr, domain, domain->getDB(), 
                                        NULL, NULL, &allPredGndingsAreQueries,
                                        true, out, amapPos, gndPredHashArray,
                                        gndPredProbs);
      if (!ok) { cout <<"Failed to create query predicates."<< endl; exit(-1); }
    }
    
    if (!(amapPos || amapAll))
      delete gndPredProbs;
	  */
  }
    // Eager version: Queries have already been generated and we can get the
    // information directly from the state
  else
  {
    if (amapPos)
      inference->printTruePreds(out);
    else
    {
      for (int i = 0; i < queryPreds.size(); i++)
      {
        int val = (queryPredValues[i] == TRUE) ? 1 : 0;

          // Prob is smoothed in inference->getProbability
        double prob = inference->getProbability((GroundPredicate*)queryPreds[i]);

        queryPreds[i]->print(out, domain);
        out << " " << prob << " " << val << endl;
      }
    }
  }
}

void printResults(const string& queryFile, const string& queryPredsStr,
                  Domain *domain, ostream& out,
                  GroundPredicateHashArray* const &queries,
                  Inference* const &inference, VariableState* const &state)
{
	// Lazy version: Have to generate the queries from the file or query string.
	// This involves calling createQueryFilePreds / createComLineQueryPreds
  if (aLazy)
  {
    const GroundPredicateHashArray* gndPredHashArray = NULL;
    Array<double>* gndPredProbs = NULL;
      // Inference algorithms with probs: have to retrieve this info from state.
      // These are the ground preds which have been brought into memory. All
      // others have always been false throughout sampling.
    if (!(amapPos || amapAll))
    {
      gndPredHashArray = state->getGndPredHashArrayPtr();
      gndPredProbs = new Array<double>;
      gndPredProbs->growToSize(gndPredHashArray->size());
      for (int i = 0; i < gndPredProbs->size(); i++)
        (*gndPredProbs)[i] =
          inference->getProbability((*gndPredHashArray)[i]);
    }

    if (queryFile.length() > 0)
    {
      cout << "Writing query predicates that are specified in query file..."
           << endl;
      bool ok = createQueryFilePreds(queryFile, domain, domain->getDB(), NULL,
                                     NULL, true, out, amapPos,
                                     gndPredHashArray, gndPredProbs, NULL);
      if (!ok) { cout <<"Failed to create query predicates."<< endl; exit(-1); }
    }

    Array<int> allPredGndingsAreQueries;
    allPredGndingsAreQueries.growToSize(domain->getNumPredicates(), false);
    if (queryPredsStr.length() > 0)
    {
      cout << "Writing query predicates that are specified on command line..." 
           << endl;
      bool ok = createComLineQueryPreds(queryPredsStr, domain, domain->getDB(), 
                                        NULL, NULL, &allPredGndingsAreQueries,
                                        true, out, amapPos, gndPredHashArray,
                                        gndPredProbs, NULL);
      if (!ok) { cout <<"Failed to create query predicates."<< endl; exit(-1); }
    }

    if (!(amapPos || amapAll))
      delete gndPredProbs;
  }
    // Eager version: Queries have already been generated and we can get the
    // information directly from the state
  else
  {
    if (amapPos)
      inference->printTruePreds(out);
    else
    {
      inference->printQFProbs(out, domain);
      if (abpInfer || aefbpInfer)
      {
        inference->printProbabilities(out);
      }
      else
      {
        for (int i = 0; i < queries->size(); i++)
        {
            // Prob is smoothed in inference->getProbability
          double prob = inference->getProbability((*queries)[i]);
          (*queries)[i]->print(out, domain); out << " " << prob << endl;
        }
      }
    }
  }
}


/**
 * The specified inference algorithm is run. First, the MLN and evidence files
 * are parsed and the database is filled. All evidence predicates are
 * closed-world by default (this can be changed with the -o option) and all
 * non-evidence predicates (query and hidden predicates) are open-world by
 * default (this can be changed with the -c option, however query atoms are
 * always open-world).
 */
int main(int argc, char* argv[])
{
  ///////////////////////////  read & prepare parameters ///////////////////////
  ARGS::parse(argc, argv, &cout);
  Timer timer;
  double begSec = timer.time(); 

  Array<Predicate *> queryPreds;
  Array<TruthValue> queryPredValues;

  ofstream resultsOut(aresultsFile);
  if (!resultsOut.good())
  { cout << "ERROR: unable to open " << aresultsFile << endl; return -1; }

  Domain* domain = NULL;
  Inference* inference = NULL;
  if (buildInference(inference, domain, aisQueryEvidence, queryPreds,
                     queryPredValues) > -1)
  {
	  //for(int i=0;i<queryPreds.size();i++)
	  //{
		//  queryPreds[i]->print(cout,domain);
		 // cout<<endl;
	  //}
    double initTime, runTime;
	Timer timer1;
	
	timer1.reset();
	inference->init();
    initTime = timer1.time();
	
	timer1.reset();
	
	  // No inference, just output network
    if (aoutputNetwork)
    {
      cout << "Writing network to file ..." << endl;
      inference->printNetwork(resultsOut);
    }
      // Perform inference
    else
    {
      if (adecisionInfer)
      {
        BP* bp = dynamic_cast<BP*>(inference);
        if (bp) bp->runDecisionBP();
      }
      else
      {
        inference->infer();
      }

      runTime = timer1.time();
      cout<<"Time-Results: Init "<<initTime<<" Run "<<runTime<<" Total "<<(initTime+runTime)<<endl;
	
	
      if (aHybrid)
	  { 
        printResults(queryFile, queryPredsStr, domain, resultsOut, &queries,
                     inference, inference->getHState());	
	  }
      else
      {
        if (adecisionInfer)
        {
          BP* bp = dynamic_cast<BP*>(inference);
          if (bp) bp->printDecisionResults(resultsOut);          
        }
        else
        {
	      printResults(queryFile, queryPredsStr, domain, resultsOut, &queries,
                       inference, inference->getState());
        }
      }
    }
  }

  resultsOut.close();
  if (domain) delete domain;
  for (int i = 0; i < knownQueries.size(); i++)
    if (knownQueries[i]) delete knownQueries[i];
  if (inference) delete inference;
  
  cout << "total time taken = "; Timer::printTime(cout, timer.time()-begSec);
  cout << endl;
}

