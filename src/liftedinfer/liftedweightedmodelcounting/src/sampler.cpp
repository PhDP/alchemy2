#include "sampler.h"
#include "randomgenutil.h"
#include <algorithm>
using namespace std;

LSampler::LSampler()
{
	//compute the static distribution table
	computeDistributionTable();
}

LSampler::~LSampler()
{
	cStoredDistribution.clear();
}


void LSampler::computeCumulativeDistribution(int domainSize,vector<LogDouble>& cdf)
{
	LogDouble norm_const(0,false);
	cdf.resize(domainSize+1);
	vector<LogDouble> prob(domainSize+1);
	for(unsigned int j=0;j<=domainSize;j++)
	{
		cdf[j] = LogDouble::Binomial(domainSize,j);
		norm_const += cdf[j];
	}
	cdf[0] = cdf[0]/norm_const;
	for(unsigned int j=1;j<cdf.size()-1;j++)
	{
		cdf[j] = cdf[j]/norm_const;
		cdf[j] = cdf[j-1] + cdf[j];
	}
	cdf[cdf.size()-1] = LogDouble(1,false);

}

void LSampler::computeDistributionTable()
{
	cStoredDistribution.resize(DTABLESIZE+1);
	for(unsigned int i=1;i<=DTABLESIZE;i++)
	{
		LogDouble cumulativeDist;
		LogDouble norm_const(0,false);
		computeCumulativeDistribution(i, cStoredDistribution[i]);
	}
#ifdef __DEBUG_PRINT1__
	for(unsigned int i=0;i<cStoredDistribution.size();i++)
	{
		for(unsigned int j=0;j<cStoredDistribution[i].size();j++)
			cout<<cStoredDistribution[i].at(j)<<" ";
		cout<<endl;
	}
#endif
}


void LSampler::doUniformDistSampling(int domainSize, vector<bool> completedGroups, int completedCount, 
	LogDouble& probOfSample, int& numTrue, vector<int>& truthValue)
{
	vector<LogDouble> cDistribution;
	if(domainSize <= DTABLESIZE)
	{
		//use stored distribution
		cDistribution = cStoredDistribution[domainSize];
	}
	else
	{
		//compute the distribution
		computeCumulativeDistribution(domainSize,cDistribution);
	}

	int groupNumber = 0;
	//check if we need to change the distribution due to zero values
	if(completedCount > 0)
	{
		LogDouble new_norm;
		vector<LogDouble> new_distribution;
		vector<int> numOnes;
		for(unsigned int t=0;t<completedGroups.size();t++)
		{
			if(!completedGroups[t])
			{
				LogDouble pVal;
				if(t==0)
					pVal = cDistribution[0];
				else
					pVal = LogDouble(exp(cDistribution[t].value) - exp(cDistribution[t-1].value),false);
				new_distribution.push_back(pVal);
				numOnes.push_back(t);
				new_norm += pVal;
			}
		}
		new_distribution[0] = new_distribution[0]/new_norm;
		for(unsigned int t=1;t<new_distribution.size()-1;t++)
		{
			new_distribution[t] = new_distribution[t]/new_norm;
			new_distribution[t] = new_distribution[t] + new_distribution[t-1];
		}
		new_distribution[new_distribution.size()-1] = LogDouble(1,false);
		double u = LvRandomGenUtil::Instance()->getNormRand();
		LogDouble lu(u,false);
		int index = lower_bound(new_distribution.begin(),new_distribution.end(),lu) - new_distribution.begin();
		groupNumber = numOnes[index];
		if(index == 0)
			probOfSample = new_distribution[0];
		else
			probOfSample = LogDouble(exp(new_distribution[index].value)-exp(new_distribution[index-1].value),false);
	}
	else
	{
		//sample from original distribution
		double u = LvRandomGenUtil::Instance()->getNormRand();
		LogDouble lu(u,false);
		groupNumber = lower_bound(cDistribution.begin(),cDistribution.end(),lu) - cDistribution.begin();
		if(groupNumber == 0)
			probOfSample = cDistribution[0];
		else
			probOfSample = LogDouble(exp(cDistribution[groupNumber].value)-exp(cDistribution[groupNumber-1].value),false);
	}
#ifdef __DEBUG_PRINT__
	cout<<"Probability of group = "<<probOfSample<<endl;
#endif
	
	truthValue.resize(domainSize);
	//Generate a random 0-1 vector with appropriate number of 1's
	int oneCount=0;
	while(oneCount != groupNumber)
	{
		int positionOnes = LvRandomGenUtil::Instance()->getRandomPosition(domainSize);
		if(truthValue[positionOnes]!=1)
		{
			truthValue[positionOnes]=1;
			oneCount++;
		}
	}
	numTrue=groupNumber;
#ifdef __DEBUG_PRINT__
	cout<<"Generated Truth Vector = [";
	for(unsigned int i=0;i<truthValue.size();i++)
		cout<<truthValue[i]<<" ";
	cout<<"]"<<endl;
#endif
}
