#ifndef LBLOCK_INFERENCE_H_
#define LBLOCK_INFERENCE_H_

#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;
#include "lvrmln.h"
#include "decomposer.h"
#include "splitter.h"
#include "logdouble.h"
#include "heuristics.h"
#include "normalizer.h"
#include "clusterutil.h"
#include "ptptree.h"
#include "lvrnormpropagation.h"
#include "lvrsingletonnormpropagation.h"
#include "lvrptptreesampling.h"

struct LBlockExactInference
{
	LBlockExactInference(LvrMLN& mln);
	~LBlockExactInference()
	{
		delete decomposer;
		delete heuristics;
		delete lvrPTPTreeSampling;
	}
	LogDouble weightedModelCount(vector<WClause*>& CNF,int parentId = -1,int branchNo = -1);
	int doExactInferenceOnCluster(vector<WClause*>& CNF);
	int doMockInferenceOnCluster(vector<WClause*>& CNF,int& lvpCost,int remainingCost);
	LPTPTree& getPTPTreeRef()
	{
		return lvrPTPTreeSampling->getTreeRef();
	}
private:
	bool decomposeCNF(vector<WClause*>& CNF,int& powerFactor,vector<map<int,int> >& decomposerMappings,
		vector<int>& domainPowerFactors);
	LogDouble CNFWeight(vector<WClause*>& CNF);
	LvrMLN& mln;
	LDecomposer* decomposer;
	LHeuristics* heuristics;
	void seperateDisjointCNF(vector<WClause*> CNF, vector<vector<WClause*> >& disjointCNFList);
	int mockWeightedModelCount(vector<WClause*>& CNF);
	int id;
	bool nonPolyNodeFound;
	bool costExceeded;
	int pCost;
	int baselineCost;
	int status;
	LvrSingletonNormPropagation* lvrSingletonNormPropagation;
	time_t startTime;
	LvrPTPTreeSampling* lvrPTPTreeSampling;
};


#endif
