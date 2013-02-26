#include "mathutils.h"
#include "randomgenutil.h"
#include <algorithm>
using namespace std;

LMathUtils* LMathUtils::m_pInstance = NULL;

LMathUtils* LMathUtils::Instance()
{
   if (!m_pInstance)
      m_pInstance = new LMathUtils;
   return m_pInstance;
}

LMathUtils::LMathUtils()
{
	constructMatrix();
	computeDistributionTable();
}
/*
double LMathUtils::factorial(int n)
{
	if(n==0)
		return 1.0;
	double fact = n;
	for(int i=n-1;i>1;i--)
	{
		fact *= i;
	}
	return fact;
}

double LMathUtils::binomVal(int n,int r)
{
	double fact = factorial(n)/(factorial(n-r)*factorial(r));	
	return fact;
}
*/
void LMathUtils::constructMatrix()
{
	binaryCoeffMatrix.resize(MAXCOEFFTABLESIZE + 1);
	for(unsigned int n=0;n<=MAXCOEFFTABLESIZE;n++)
	{
		binaryCoeffMatrix[n].resize(n+1);
		for(unsigned int i=0;i<=n;i++)
		{
			binaryCoeffMatrix[n].at(i) = LogDouble::Binomial(n,i);
		}
	}
}

LogDouble LMathUtils::getCoeff(int n,int i)
{
	if(i > n)
		return LogDouble(0,false);
	if(n > MAXCOEFFTABLESIZE)
	{
		//compute the value
		return LogDouble::Binomial(n,i);
	}
	else
		return binaryCoeffMatrix[n].at(i);
}

void LMathUtils::computeCumulativeDistribution(int domainSize,vector<LogDouble>& cdf)
{
	LogDouble norm_const(0,false);
	cdf.resize(domainSize+1);
	vector<LogDouble> prob(domainSize+1);
	for(unsigned int j=0;j<=domainSize;j++)
	{
		cdf[j] = getCoeff(domainSize,j);
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

void LMathUtils::computeDistributionTable()
{
	cStoredDistribution.resize(MAXCOEFFTABLESIZE+1);
	for(unsigned int i=1;i<=MAXCOEFFTABLESIZE;i++)
	{
		LogDouble cumulativeDist;
		LogDouble norm_const(0,false);
		computeCumulativeDistribution(i, cStoredDistribution[i]);
	}
}

void LMathUtils::sample(int domainSize,LogDouble& probOfSample, int& value)
{
	//sample from stored distribution
	vector<LogDouble> cDistribution;
	if(domainSize <= MAXCOEFFTABLESIZE)
	{
		//use stored distribution
		cDistribution = cStoredDistribution[domainSize];
	}
	else
	{
		//construct distribution
		computeCumulativeDistribution(domainSize,cDistribution);
	}
	double u = LvRandomGenUtil::Instance()->getNormRand();
	LogDouble lu(u,false);
	value = lower_bound(cDistribution.begin(),cDistribution.end(),lu) - cDistribution.begin();
	if(value == 0)
		probOfSample = cDistribution[0];
	else
		probOfSample = LogDouble(exp(cDistribution[value].value)-exp(cDistribution[value-1].value),false);
}

void LMathUtils::toBinary(int num,vector<int>& binVector,int numGroundings)
{
	int mask = 1;
	int iter=numGroundings-1;
	binVector.resize(numGroundings);
	while(num)
	{
		int x = num & mask;
		binVector[iter--] = x;
		num = num >> 1;
	}

}