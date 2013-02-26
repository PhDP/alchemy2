#ifndef __LVRIMPORTANCE_DISTRIBUTION
#define __LVRIMPORTANCE_DISTRIBUTION
#include <map>
#include <vector>
using namespace std;
struct LvrDistributionElement
{
	map<int,vector<double> > distributionElement;
};
struct LvrImportanceDistribution
{
	map<int,map<int,vector<LvrDistributionElement*> > distributions;
	LvrImportanceDistribution(){}
	~LvrImportanceDistribution(){}


};
#endif
