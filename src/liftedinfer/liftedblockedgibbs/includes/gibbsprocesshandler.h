#ifndef GIBBS_PROCESS_HANDLER_H
#define GIBBS_PROCESS_HANDLER_H

#include <sys/types.h>
#ifndef _MSC_VER 
#include <unistd.h>
#include <cstdio>
#endif
#include "clustercreator.h"
#include "lbgsampler.h"
#include "lvrmln.h"
#include "lvparams.h"
using namespace std;

class GibbsProcessHandler
{
 public:
	GibbsProcessHandler(LvrMLN& mln_):mln(mln_)
	{
		lSampler = new LBGSampler(mln);
		lCCreator = new LClusterCreator(&mln);
	}

	~GibbsProcessHandler()
	{
		delete lCCreator;
		delete lSampler;
	}

	void startLBGForMar(LvrParams* params)
	{
		if(params->operationMode == ECLUSTERING)
		{
			cout<<"Starting Lifted Blocked Gibbs in clustering mode.."<<endl;
			if(params->outClusterFile.size()==0)
				params->outClusterFile = CLUSTERFILE;
			//remove the existing cluster file
	#ifndef _MSC_VER
				remove(params->outClusterFile.c_str());
	#endif
			if(params->clusteringMode == ERANDOM)
				lCCreator->startRandomClustering(params);
			else
				lCCreator->startClustering(params);
		}
		else if(params->operationMode == ESAMPLER)
		{
			cout<<"Starting Lifted Blocked Gibbs in sampling mode.."<<endl;
			if(params->inClusterFile.size()==0)
				params->inClusterFile = CLUSTERFILE;
			lSampler->startLVBGibbs(params);
		}
		else //parallel
		{
#ifndef _MSC_VER 
			cout<<"Starting Lifted Blocked Gibbs in parallel mode.."<<endl;
			params->inClusterFile = CLUSTERFILE;
			params->outClusterFile = CLUSTERFILE;
			remove(CLUSTERFILE);
			cout<<"Forking LBG Sampling and Clustering Processes"<<endl;
			pid_t pID = fork();
			if (pID == 0)
			{
				//start a child process to do the clustering
				if(params->clusteringMode == ERANDOM)
					lCCreator->startRandomClustering(params);
				else
					lCCreator->startClustering(params);
				//can exit child once clustering is complete
				exit(0);
			}
			else
			{
				lSampler->startLVBGibbs(params);
			}
#endif
		}
	}


	void initLBGForWeightLearning(MCMCParams* params)
	{
	}

	void stepLBGForWeightLearning(LvrParams* params)
	{
		//no clustering simply sample
		lSampler->startLVBGibbs(params);
	}

 private:
  LClusterCreator* lCCreator;
  LvrMLN& mln;
  vector<Atom*> evidenceAtoms;
  vector<Atom*> queryAtoms;
  LBGSampler* lSampler;
};

#endif
