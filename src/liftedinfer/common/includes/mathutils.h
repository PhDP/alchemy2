#ifndef __LMATHUTILS
#define __LMATHUTILS
#include <vector>
using namespace std;
#include "logdouble.h"

#define MAXCOEFFTABLESIZE 50
//singleton class
struct LMathUtils
{
	static LMathUtils* Instance();
	LogDouble getCoeff(int n,int r);
	void sample(int domainSize,LogDouble& probOfSample, int& value);
	void toBinary(int num,vector<int>& binVector,int numGroundings);
private:
	vector<vector<LogDouble> > binaryCoeffMatrix;
	vector<vector<LogDouble> > cStoredDistribution;
	LMathUtils();
	void constructMatrix();
	void computeCumulativeDistribution(int domainSize,vector<LogDouble>& cdf);
	void computeDistributionTable();
	//double factorial(int n);
	//double binomVal(int n,int r);
	static LMathUtils* m_pInstance;
};
#endif