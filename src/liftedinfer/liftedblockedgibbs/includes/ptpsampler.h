#ifndef __PTPSAMPLER__
#define __PTPSAMPLER__
#include "lvrmln.h"
#include "blockexactinference.h"
#include "lvrcluster.h"

struct LPTPSampler
{
	LPTPSampler(LvrMLN& mln);
	~LPTPSampler();
	int runPTP(vector<WClause*>& CNF,LVRCluster* cluster);
	int runMockPTP(vector<WClause*>& CNF,LVRCluster* lvrCluster,int& lvpCost,int remainingCost);
	void setBurnIn(bool val)
	{
		burnIn = val;
	}
	bool getBurnIn()
	{
		return burnIn;
	}
private:
	LBlockExactInference* ls;
	bool burnIn;
	vector<Atom*> queryAtoms;
	vector<vector<int> > queryValues;
	vector<vector<int> > currentSampleValues;
	//0-not updated,1-updated, -1 - evidence atom
	vector<vector<int> > queryUpdatedInIteration;
	void doLightWeightSampling(LVRCluster* lvrCluster,vector<WClause*> clauses);
};

#endif
