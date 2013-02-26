#ifndef __LBGSAMPLER_H
#define __LBGSAMPLER_H
#include "lvrmln.h"
#include "ptpsampler.h"
#include "lvparams.h"
#include "lvrcluster.h"

struct LBGSampler
{
	LBGSampler(LvrMLN& mln_):mln(mln_)
	{
		clusterLevel = -1;
		originalNumPredicates = -1;
		numTimesSampled = 0;
		iterationsSinceLastMerge = 0;
		iterationsAcrossCalls = 0;
		lp = new LPTPSampler(mln);
	}
	~LBGSampler()
	{
		cleanup(LVRClusterList);
		delete lp;
	}
	void startLVBGibbs(LvrParams* params);
private:
	int clusterLevel;
	int originalNumPredicates;
	int numTimesSampled;
	int windowSize;
	int iterationsSinceLastMerge;
	//used in weight learning only, since called multiple times
	int iterationsAcrossCalls;
	LPTPSampler* lp;
	vector<LVRCluster*> LVRClusterList;
	LvrMLN& mln;
	bool useEvidenceCaching;
	void mergeCluster(int clusterId1,int clusterId2);
	void computeClusterLiftedSharings();
	int getClusterIndex(LVRCluster* cluster);
	void initializeClusters();
	void computeClusterMB();	
	void addPTPEvidence(int clusterIndex,vector<WClause*>& newClauses, bool storeEvidence=false);
	void getAtomsMatchingTerm(Atom* atom, int termIndex, vector<Atom*>& matchingMLNAtoms);
	void groundClusterWithMB(int clusterIndex,vector<WClause*>& newClauses);
};


#endif
