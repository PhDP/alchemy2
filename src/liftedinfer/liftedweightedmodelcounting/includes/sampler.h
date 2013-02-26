#ifndef __LSAMPLER__
#define __LSAMPLER__
#include <vector>
#include <iostream>
using namespace std;
#include "logdouble.h"
#define DTABLESIZE 100

struct LSampler
{
	void doUniformDistSampling(int domainSize, vector<bool> completedGroups, int completedCount, 
	LogDouble& probOfSample, int& numTrue, vector<int>& truthValue);
	LSampler();
	~LSampler();
private:
	vector<vector<LogDouble> > cStoredDistribution;
	void computeCumulativeDistribution(int domainSize,vector<LogDouble>& cdf);
	void computeDistributionTable();
};
#endif