#ifndef __LVRPERMUTATIONS
#define __LVRPERMUTATIONS
#include "lvrmln.h"
struct LvrPermutations
{
static void permuteTerms(vector<LvrTerm*> terms, vector<vector<int> >& permutedDomainList)
{
	if(terms.size()==0)
		return;
	map<int,int> domainIndex;
	int numberOfGroundings=1;
	for(unsigned int i=0;i<terms.size();i++)
	{
		//Initialize to first position
		domainIndex.insert(pair<int,int>(i,0));
		numberOfGroundings *= terms[i]->domain.size();
	}	
	int iter = 0;
	while(1)
	{
		vector<int> domainValues;
		for(unsigned int i=0;i<terms.size();i++)
		{
			//get current index
			int currDomainIndex = domainIndex[i];
			domainValues.push_back(terms[i]->domain[currDomainIndex]);
		}
		permutedDomainList.push_back(domainValues);
		//Update pointers
		int ind=terms.size()-1;
		domainIndex[ind]++;
		while(domainIndex[ind] == terms[ind]->domain.size())
		{
			domainIndex[ind]=0;
			if(ind!=0)
				domainIndex[ind-1]++;
			else
				break;
			ind--;
		}
		iter++;
		//check if we finished all domains
		if(iter >= numberOfGroundings)
			break;
	}
}
};
#endif