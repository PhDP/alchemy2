#include "samplingalgs.h"
#include "mathutils.h"
#include "randomgenutil.h"
#include <algorithm>
using namespace std;

LSamplingAlgs* LSamplingAlgs::m_pInstance = NULL;

LSamplingAlgs* LSamplingAlgs::Instance()
{
   if (!m_pInstance)
      m_pInstance = new LSamplingAlgs;
   return m_pInstance;
}

LSamplingAlgs::LSamplingAlgs()
{
}

void LSamplingAlgs::sample(int domainSize,int numValuesToSample,LogDouble& sampleWeight,
	vector<int>& sampledValues,ESamplingMode samplingMode, LProposalDistributionElement* lpe)
{
	switch(samplingMode) 
	{
		case EUNIFORM:
		{
				LogDouble binCoeff(1,false);
				for(unsigned int i=0;i<numValuesToSample;i++)
				{
					sampledValues.push_back(LvRandomGenUtil::Instance()->getRandomPosition(domainSize + 1));
					LogDouble b = LMathUtils::Instance()->getCoeff(domainSize,sampledValues[i]);
					binCoeff *= b;
				}
				double p = (double)numValuesToSample/(double)(domainSize + 1);
				sampleWeight = binCoeff/LogDouble(p,false);
				break;
		}
		case EBINOMIAL:
		{
			LogDouble prob(1,false);
			LogDouble binCoeff(1,false);
			for(int i=0;i<numValuesToSample;i++)
			{
				LogDouble p(1,false);
				int val;
				LMathUtils::Instance()->sample(domainSize,p,val);
				prob *= p;
				if(val<0)
					val = 0;
				sampledValues.push_back(val);
				LogDouble b = LMathUtils::Instance()->getCoeff(domainSize,sampledValues[i]);
				binCoeff *= b;
			}
			sampleWeight = binCoeff/prob;
			break;
		}
		default:
		{
			if(lpe!=NULL)
			{
				LogDouble binCoeff(1,false);
				LogDouble prob(1,false);
				if(samplingMode == EINFORMED)
					lpe->sample(prob,sampledValues,domainSize);
				else
					lpe->sample(prob,sampledValues);
				if(sampledValues.size() < numValuesToSample)
				{
					int currSize = sampledValues.size();
					//sample binomially
					for(unsigned int i=0;i<numValuesToSample - currSize;i++)
					{
						LogDouble p(1,false);
						int val;
						LMathUtils::Instance()->sample(domainSize,p,val);
						prob *= p;
						sampledValues.push_back(val);
					}
				}
				for(int i=0;i<numValuesToSample;i++)
				{
					LogDouble b = LMathUtils::Instance()->getCoeff(domainSize,sampledValues[i]);
					binCoeff *= b;
				}
				sampleWeight = binCoeff/prob;
			}
			else
			{
				//use uniform
				LogDouble binCoeff(1,false);
				for(unsigned int i=0;i<numValuesToSample;i++)
				{
					sampledValues.push_back(LvRandomGenUtil::Instance()->getRandomPosition(domainSize + 1));
					LogDouble b = LMathUtils::Instance()->getCoeff(domainSize,sampledValues[i]);
					binCoeff *= b;
				}
				double p = (double)numValuesToSample/(double)(domainSize + 1);
				sampleWeight = binCoeff/LogDouble(p,false);
			}
			break;
		}
	}
}
