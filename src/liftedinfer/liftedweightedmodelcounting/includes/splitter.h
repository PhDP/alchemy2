#ifndef LSPLITTER_H_
#define LSPLITTER_H_
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;

#include "lvrmln.h"
#include "sampler.h"
#include "resolver.h"

struct LSplitter
{
	LSampler* sampler;
	LResolver* resolver;
	LSplitter();
	bool doApproxSplit(vector<WClause*> CNF,Atom* atom, int numGroundings, vector<bool> completedGroups, int completedCount,
		LogDouble& probOfSample, int& numTrue, vector<vector<WClause*> >& resolvedList);
	void doFullGrounding(vector<WClause*> CNF,Atom* atom,vector<WClause*>& groundCNF,vector<vector<int> >& allGroundings,int& numGroundings);
	void doCompleteSplitInChunks(vector<WClause*> groundCNF,int atomId, vector<vector<int> > allGroundings, int numGroundings,
		vector<vector<WClause*> >& resolvedCNFList,int startnum, int chunkSize,bool& done);
	void doLiftedSplit(vector<WClause*> CNF, Atom* atom, vector<vector<WClause*> >& resolvedCNFList,int singletonIndex);
	void doLiftedSplitSelfJoins(vector<WClause*> CNF, Atom* atom, vector<vector<WClause*> >& resolvedCNFList,int singletonIndex);
	~LSplitter();

};
#endif