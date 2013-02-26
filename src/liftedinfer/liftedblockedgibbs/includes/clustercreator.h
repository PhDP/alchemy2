#ifndef __LCLUSTER_CREATOR__
#define __LCLUSTER_CREATOR__
#include "lvrmln.h"
#include "ptpsampler.h"
#include "normalizer.h"
#include "mcmcparams.h"
#include "lvparams.h"

struct LVRCluster;
struct LTuple
{
	vector<int> clusterIndexes;
	int CF;
};
struct LClusterCreator
{
	LClusterCreator(LvrMLN* mln_)
	{
		//make a copy of the MLN
		mln = new LvrMLN();
		for(int i=0;i<mln_->clauses.size();i++)
		{
			mln->clauses.push_back(LvrMLN::create_new_clause(mln_->clauses[i]));
		}
		for(int i=0;i<mln_->symbols.size();i++)
		{
			mln->symbols.push_back(LvrMLN::create_new_symbol(mln_->symbols[i]));
		}
		mln->maxDegree = mln_->maxDegree;
		mln->max_predicate_id = mln_->max_predicate_id;
		clusterLevel = -1;
		originalNumPredicates = -1;
		lp = new LPTPSampler(*mln);
	}

	~LClusterCreator()
	{
		delete mln;
		delete lp;
		cleanup(LVRClusterList);
	}
	void startClustering(LvrParams* params);
	void startRandomClustering(LvrParams* params);
private:
	int clusterLevel;
	int originalNumPredicates;
	int iterations;
	int windowSize;
	LPTPSampler* lp;
	vector<LVRCluster*> LVRClusterList;
	LvrMLN* mln;
	bool useEvidenceCaching;
	void mergeCluster(int clusterId1,int clusterId2);
	void computeClusterLiftedSharings();
	int getClusterIndex(LVRCluster* cluster);
	void initializeClusters();
	void computeClusterMB();
	
	void makeClusterPairs(vector<LTuple*>& clusterCombinations);
	void addPTPEvidence(int clusterIndex,vector<WClause*>& newClauses, bool storeEvidence=false);
	
	int computeCouplingFactor(int clauseIndex1,int clauseIndex2);
	int numCoOccurrences(Atom* a1,Atom* a2);
	void getAtomsMatchingTerm(Atom* atom, int termIndex, vector<Atom*>& matchingMLNAtoms);
	void groundClusterWithMB(int clusterIndex,vector<WClause*>& newClauses);

};


#endif
