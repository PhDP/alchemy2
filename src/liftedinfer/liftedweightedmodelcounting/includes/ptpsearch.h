#ifndef LPTP_SEARCH_H_
#define LPTP_SEARCH_H_

#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;
#include "lvrmln.h"
#include "decomposer.h"
#include "splitter.h"
#include "logdouble.h"
#include "extensions.h"
#include "heuristics.h"
#include "normalizer.h"
#include "clusterutil.h"
#include "lvparams.h"
#include "lvrptptreesampling.h"
#include "lvrptpsearchtree.h"

#define SPLITCHUNKSIZE 25
struct LPTPSearch
{
	LvrMLN& mln;
	LSplitter* splitter;
	LDecomposer* decomposer;
	LExtensions* extensions;
	LHeuristics* heuristics;
	LNormalizer* normalizer;
	time_t starttime;
	bool dounitpropagate;

	LogDouble startApproxWeightedModelCounting(LvrParams* params);
	LogDouble startExactWeightedModelCounting(LvrParams* params);
	void seperateDisjointCNF(vector<WClause*> CNF, vector<vector<WClause*> >& disjointCNFList);
	LPTPSearch(LvrMLN& mln);
	~LPTPSearch()
	{
		delete splitter;
		delete decomposer;
		delete extensions;
		delete heuristics;
		delete normalizer;
		delete lvrPTPTreeSampling;
		delete lvrptpSearchTree;
	}

private:
	bool decomposeCNF(vector<WClause*>& CNF,int& powerFactor);
	LogDouble CNFWeight(vector<WClause*>& CNF,set<int>* dontCareQueries=NULL);
	int CheckCNFSatisfied(vector<WClause*>& CNF,double& satPercent);
	//LogDouble weightedModelCount(vector<WClause*>& CNF,int parentId,int parentBranchNo);
	LogDouble weightedModelCount(vector<WClause*>& CNF,LvrPTPNode* parent=NULL,int parentBranchNo=-1);
	LogDouble weightedModelCountApprox(vector<WClause*>& CNF);
	int id;
	LvrPTPTreeSampling* lvrPTPTreeSampling;
	LvrPTPSearchTree* lvrptpSearchTree;
	int CheckCNFSatisfiedWithQueries(vector<WClause*>& CNF);
};


#endif /* SEARCH_H_ */
