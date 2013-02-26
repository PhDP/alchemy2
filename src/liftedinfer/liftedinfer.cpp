#include <unistd.h>
#include <fstream>
#include <climits>
#include <sys/times.h>

#include "infer.h"
#include "arguments.h"
#include "liftedalgsconvertor.h"
#include "fileutils.h"
#include "paramconstants.h"
#include "queryupdater.h"

extern const char* ZZ_TMP_FILE_POSTFIX; //defined in fol.y

char* aInputClusterFile = NULL;
char* aOutputClusterFile = NULL;


bool aisUseDumpFiles = false;

bool aisLiftedGibbsWithCluster = false;
bool aisLiftedGibbs = false;
bool aisClustering = false;
bool aisRandomizedClustering = false;
bool aisSequentialClustering = false;

bool aisWMPTPEZ = false;
bool aisWMPTPAZ = false;

bool aisLIS = false;

double adaptiveLearningRate = -1;
int proposalUpdateInterval = -1;

int aimd = 1;
int acmode = 0;
int ptpcacheSize = 0;
double abmf = -1;

  // TODO: List the arguments common to learnwts and inference in
  // inferenceargs.h. This can't be done with a static array.
ARGS ARGS::Args[] = 
{
    // BEGIN: Common arguments
  ARGS("i", ARGS::Req, ainMLNFiles,
       "Comma-separated input .mln files."),

    // BEGIN: Common inference arguments
    // BEGIN: MCMC args
  ARGS("burnMaxSteps", ARGS::Opt, amcmcBurnMaxSteps,
       "[50] (MCMC) Maximum number of burn-in steps."),

  ARGS("maxSteps", ARGS::Opt, amcmcMaxSteps, 
       "[1000] (MCMC) Maximum number of MCMC sampling steps."),

  ARGS("maxSeconds", ARGS::Opt, amcmcMaxSeconds, 
       "[10000] (MCMC) Max number of seconds to run MCMC."),
    // END: MCMC args
  

    // BEGIN: Args specific to stand-alone inference
  ARGS("e", ARGS::Opt, aevidenceFiles,
       "Comma-separated .db files containing known ground atoms (evidence), "
       "including function definitions."),

  ARGS("r", ARGS::Opt, aresultsFile,
       "[result.dat] The probability estimates are written to this file."),
    
  ARGS("q", ARGS::Opt, aqueryPredsStr, 
       "Query atoms (comma-separated with no space)  "
       ",e.g., cancer,smokes(x),friends(Stan,x). Query atoms are always "
       "open world."),

  ARGS("ic", ARGS::Opt, aInputClusterFile,
	   "[cluster.cls] Input cluster file containing the clustering,"
	   "to be used with -lvg option"),

  ARGS("oc", ARGS::Opt, aOutputClusterFile,
       "[cluster.cls] Output cluster file containing the clustering,"
	   "used with the -cls option"),

   //algorithm options
   ARGS("lvgc", ARGS::Tog, aisLiftedGibbsWithCluster,
      "Run Lifted Blocked Gibbs sampling and clustering in parallel for marginals,"
	  "Note:The two algorithms are run as concurrent processes"),

   ARGS("lvg", ARGS::Tog, aisLiftedGibbs,
      "Compute marginals using Lifted Blocked Gibbs sampling with a given clustering,"
	  "specified using -ic option"),
  
  ARGS("cls", ARGS::Tog, aisClustering,
     "Generate lifted clusters using the clustering algorithm"),

  ARGS("cmode", ARGS::Opt, acmode,
     "[0] Clustering mode 0:sequential 1:randomized"),

  ARGS("ptpe", ARGS::Tog, aisWMPTPEZ,
       "Compute exact partition function,marginals using the WM PTP algorithm"),

  //ARGS("ptpz", ARGS::Tog, aisWMPTPAZ,
	//	"Compute approx partition function using the WM PTP algorithm"),

  ARGS("lis", ARGS::Tog, aisLIS,
        "Compute partition function, marginals using the LIS algorithm"),

 ARGS("imd", ARGS::Opt, aimd,
   	  "[0] Type of importance distribution 0:Informed 1:Binomial 2:Uniform"),

  ARGS("lr", ARGS::Opt, adaptiveLearningRate,
	  "[0.001] Learning rate for the adaptive LIS algorithm"),

  ARGS("ui",ARGS::Opt,proposalUpdateInterval,
	   "[50] update interval when using LIS with informed distribution"),
  ARGS("cache",ARGS::Opt,ptpcacheSize,"[0] cache size to be used with PTP Exact Inference"),
       // END: Args specific to stand-alone inference
  ARGS("bmf",ARGS::Opt,abmf,"[1.25] baseline multiplicative factor for clustering;"
		  "The ptp cost of cluster graph <= bmf*propositionalgibbscost"),
  ARGS()
};


int processInput(Domain* domain, MLN* mln, GroundPredicateHashArray& queries)
{
	  string inMLNFile, wkMLNFile, evidenceFile;
	  StringHashArray queryPredNames;
	  StringHashArray owPredNames;
	  StringHashArray cwPredNames;
	  Array<string> constFilesArr;
	  Array<string> evidenceFilesArr;

	  
	  //the second .mln file to the last one in ainMLNFiles _may_ be used 
	  //to hold constants, so they are held in constFilesArr. They will be
	  //included into the first .mln file.

		//extract .mln, .db file names
	  extractFileNames(ainMLNFiles, constFilesArr);
	  assert(constFilesArr.size() >= 1);
	  inMLNFile.append(constFilesArr[0]);
	  constFilesArr.removeItem(0);
	  extractFileNames(aevidenceFiles, evidenceFilesArr);
	  
	  if (aqueryPredsStr) queryPredsStr.append(aqueryPredsStr);
	  if (aqueryFile) queryFile.append(aqueryFile);

	  if (queryPredsStr.length() == 0 && queryFile.length() == 0 &&
			  (aisLiftedGibbs || aisLiftedGibbsWithCluster))
	  { cout << "No query predicates specified" << endl; return -1; }


		//extract names of all query predicates
	  if (queryPredsStr.length() > 0 || queryFile.length() > 0)
	  {
		if (!extractPredNames(queryPredsStr, &queryFile, queryPredNames)) return -1;
	  }

		//extract names of open-world evidence predicates
	  if (aOpenWorldPredsStr)
	  {
		if (!extractPredNames(string(aOpenWorldPredsStr), NULL, owPredNames)) 
		  return -1;
		assert(owPredNames.size() > 0);
	  }

		//extract names of closed-world non-evidence predicates
	  if (aClosedWorldPredsStr)
	  {
		if (!extractPredNames(string(aClosedWorldPredsStr), NULL, cwPredNames)) 
		  return -1;
		assert(cwPredNames.size() > 0);
		if (!checkQueryPredsNotInClosedWorldPreds(queryPredNames, cwPredNames))
		  return -1;
	  }

	  //////////////////// read in clauses & evidence predicates //////////////////

	  cout << "Reading formulas and evidence predicates..." << endl;

		// Copy inMLNFile to workingMLNFile & app '#include "evid.db"'
	  string::size_type bslash = inMLNFile.rfind("/");
	  string tmp = (bslash == string::npos) ? 
				   inMLNFile:inMLNFile.substr(bslash+1,inMLNFile.length()-bslash-1);
	  char buf[100];
	  sprintf(buf, "%s%d%s", tmp.c_str(), getpid(), ZZ_TMP_FILE_POSTFIX);
	  wkMLNFile = buf;
	  copyFileAndAppendDbFile(inMLNFile, wkMLNFile,
							  evidenceFilesArr, constFilesArr);

		// Parse wkMLNFile, and create the domain, MLN, database
	  bool addUnitClauses = false;
	  bool mustHaveWtOrFullStop = true;
	  bool warnAboutDupGndPreds = true;
	  bool flipWtsOfFlippedClause = true;
	  bool allPredsExceptQueriesAreCW = false;
	  //bool allPredsExceptQueriesAreCW = owPredNames.empty();
	  Domain* forCheckingPlusTypes = NULL;

		// Parse as if lazy inference is set to true to set evidence atoms in DB
		// If lazy is not used, this is removed from DB
	  cout<<wkMLNFile.c_str()<<endl;
	  if (!runYYParser(mln, domain, wkMLNFile.c_str(), allPredsExceptQueriesAreCW, 
					   &owPredNames, &cwPredNames, &queryPredNames, addUnitClauses, 
					   warnAboutDupGndPreds, 0, mustHaveWtOrFullStop, 
					   forCheckingPlusTypes, true, flipWtsOfFlippedClause))
	  {
		unlink(wkMLNFile.c_str());
		//return -1;
		exit(-1);
	  }

	  unlink(wkMLNFile.c_str());
	  const FormulaAndClausesArray* fca = mln->getFormulaAndClausesArray();
	  for (int i = 0; i < fca->size(); i++)
	  {
		IndexClauseHashArray* indexClauses = (*fca)[i]->indexClauses;
		for (int j = 0; j < indexClauses->size(); j++)
		{
		  int idx = (*indexClauses)[j]->index;
		  Clause* c = (*indexClauses)[j]->clause;
		  cout << "idx " << idx << ": ";
		  c->printWithWtAndStrVar(cout, domain);
		  cout << endl;
		}
	  }
		//////////////////////////// run inference /////////////////////////////////

		///////////////////////// read & create query predicates ///////////////////
	  Array<int> allPredGndingsAreQueries;
	  Array<Array<Predicate* >* >* queryFormulas =  new Array<Array<Predicate*> *>;
		if (queryFile.length() > 0)
		{
		  cout << "Reading query predicates that are specified in query file..."
			   << endl;
		  bool ok = createQueryFilePreds(queryFile, domain, domain->getDB(),
										 &queries, &knownQueries, queryFormulas);
		  if (!ok) { cout<<"Failed to create query predicates."<<endl; exit(-1); }
		}

		allPredGndingsAreQueries.growToSize(domain->getNumPredicates(), false);
		if (queryPredsStr.length() > 0)
		{
		  // unePreds = unknown non-evidence predicates
		  // nePreds  = known non-evidence predicates
		  GroundPredicateHashArray unePreds;
		  GroundPredicateHashArray knePreds;
		  bool ok = createComLineQueryPreds(queryPredsStr, domain, 
									  domain->getDB(), &unePreds, &knePreds,
									  &allPredGndingsAreQueries, queryFormulas);
		  if (!ok) { cout<<"Failed to create query predicates."<<endl; exit(-1); }
		  //evidence groundings
		  queries = unePreds;
		}
		return 0;
}

void runLIS(LiftedAlgsConvertor* liftedAlgsConvertor,bool computeZ = true)
{
	//Lifted Importance Sampling
	LvrParams* params = new LvrParams;
	params->maxSteps     = amcmcMaxSteps;
	params->maxSeconds   = amcmcMaxSeconds;
	params->useDumpFiles = aisUseDumpFiles;
	if(aresultsFile)
		params->resultFile = aresultsFile;
	else
		params->resultFile = RESULTFILE;
	if(aimd == 0)
		params->samplingMode = EINFORMED;
	else if(aimd == 1)
		params->samplingMode = EBINOMIAL;
	else
		params->samplingMode = EUNIFORM;
	params->learningRate = adaptiveLearningRate;
	params->proposalUpdateInterval = proposalUpdateInterval;
	params->isWeightLearning = false;
	if(computeZ)
		liftedAlgsConvertor->doLISApproxZ(params);
	else
		liftedAlgsConvertor->doLISApproxMar(params);
	delete params;
}

void runLBGibbs(LiftedAlgsConvertor* liftedAlgsConvertor)
{
	LvrParams* params = new LvrParams;
	//Gibbs parameters
	params->burnMaxSteps = amcmcBurnMaxSteps;
	params->burnMinSteps = amcmcBurnMinSteps;
	params->maxSteps     = amcmcMaxSteps;
	params->maxSeconds   = amcmcMaxSeconds;
	params->useDumpFiles = aisUseDumpFiles;
	params->baselineCostMultiplicativeFactor = abmf;
	if(aresultsFile)
		params->resultFile = aresultsFile;
	else
		params->resultFile = RESULTFILE;
	if(aInputClusterFile)
		params->inClusterFile = aInputClusterFile;
	else
		cout<<"Input cluster file is not specified..using default "<<CLUSTERFILE<<endl;
	if(aOutputClusterFile)
		params->outClusterFile = aOutputClusterFile;

	if(aisLiftedGibbs)
		params->operationMode = ESAMPLER;
	else if(aisClustering)
		params->operationMode = ECLUSTERING;
	else
		params->operationMode = EPARALLEL;

	if(acmode==0)
		params->clusteringMode = ESEQUENTIAL;
	else
		params->clusteringMode = ERANDOM;
	params->isWeightLearning = false;
	liftedAlgsConvertor->doLBGMar(params);
	delete params;

}
void runPTP(LiftedAlgsConvertor* liftedAlgsConvertor)
{
	//Lifted WMC
	LvrParams* params = new LvrParams;
	params->maxSteps     = amcmcMaxSteps;
	params->maxSeconds   = amcmcMaxSeconds;
	if(aresultsFile)
		params->resultFile = aresultsFile;
	else
		params->resultFile = RESULTFILE;

	if(aisWMPTPAZ)
		liftedAlgsConvertor->doWMPTPApproxZ(params);
	else
	{
		params->ptpCacheSize = ptpcacheSize;
		liftedAlgsConvertor->doWMPTPExactZ(params);
	}
	delete params;
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
  amcmcBurnMaxSteps = -1;
  amcmcBurnMinSteps = -1;
  amcmcMaxSteps = -1;
  amcmcMaxSeconds = -1;
  ARGS::parse(argc, argv, &cout);
  Timer timer;
  double begSec = timer.time();
  Array<Predicate *> queryPreds;
  Array<TruthValue> queryPredValues;
  GroundPredicateHashArray queries;
  Domain* domain = new Domain;
  MLN* mln = new MLN();
  if(processInput(domain,mln,queries)==-1)
  {
	  cout<<"Error processing input..exiting"<<endl;
	  return -1;
  }
  LiftedAlgsConvertor* liftedAlgsConvertor = new LiftedAlgsConvertor(mln,domain,queries);
    if(aisLIS)
	{
    	bool computeZ=true;
    	//check if marginals are to be computed
    	if(queries.size() > 0)
    		computeZ=false;
    	runLIS(liftedAlgsConvertor,computeZ);
	}
	else if(aisLiftedGibbsWithCluster || aisLiftedGibbs || aisClustering)
	{
		runLBGibbs(liftedAlgsConvertor);
	}
	else if(aisWMPTPEZ || aisWMPTPAZ)
	{
		runPTP(liftedAlgsConvertor);
	}
	else
	{
		if(queries.size() == 0)
		{
			cout<<"No algorithm or query is specified to be run. Running LIS to compute approxZ.."<<endl;
			aisLIS = true;
			runLIS(liftedAlgsConvertor);
		}
		else
		{
			cout<<"No algorithm is specified to be run. Running Lifted Gibbs Sampling to compute Marginals.."<<endl;
			aisLiftedGibbsWithCluster=true;
			runLBGibbs(liftedAlgsConvertor);
		}
	}
  if (domain) delete domain;
  if(mln) delete mln;
  delete liftedAlgsConvertor;
  cout << "total time taken = "; Timer::printTime(cout, timer.time()-begSec);
  cout << endl;
  return 0;
}

