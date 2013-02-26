#ifndef __LSamplingAlgs
#define __LSamplingAlgs
#include "logdouble.h"
#include "proposalstructure.h"

#include<vector>
using namespace std;
#include "lvparams.h"

struct LSamplingAlgs
{
	static LSamplingAlgs* Instance();
	void sample(int domainSize,int numValuesToSample,LogDouble& sampleWeight,vector<int>& sampledValues,
		ESamplingMode samplingMode = EUNIFORM, LProposalDistributionElement* lpe = NULL);
private:
	LSamplingAlgs();
	static LSamplingAlgs* m_pInstance;

};

#endif