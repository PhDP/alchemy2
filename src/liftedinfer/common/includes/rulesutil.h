#ifndef LRULESUTIL_H_
#define LRULESUTIL_H_
#include "lvrmln.h"
#include "cleanuputils.h"
#include "normalizer.h"
#include "lvrpermutations.h"

struct LRulesUtil
{
	static bool computeIsolatedTerms(Atom* atom, vector<WClause*> clauses,vector<bool>& isolatedTerms)
	{
		isolatedTerms.clear();
		isolatedTerms.resize(atom->terms.size());
		int falseCnt = 0;
		for(unsigned int i=0;i<isolatedTerms.size();i++)
		{
			//a constant is never isolated
			if(atom->terms[i]->domain.size()==1)
			{
				falseCnt++;
				isolatedTerms[i] = false;
			}
			else
			{
				isolatedTerms[i] = true;
			}
		}
		if(falseCnt == isolatedTerms.size())
			return false;

		for(unsigned int i=0;i<clauses.size();i++)
		{
			for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
			{
				if(clauses[i]->atoms[j]->symbol->id == atom->symbol->id)
				{
					for(unsigned int k=0;k<clauses[i]->atoms[j]->terms.size();k++)
					{
						if(isolatedTerms[k])
						{
							if(clauses[i]->atoms[j]->terms[k]->domain.size()>1)
							{
								bool shared = checkOccurrence(clauses[i]->atoms[j]->terms[k],clauses[i],j);
								if(shared)
								{
									isolatedTerms[k]=false;
									falseCnt++;
									if(falseCnt == isolatedTerms.size())
										return false;
								}
							}
						}
					}
				}
			}
		}
		return true;
	}

	static bool isBlocked(Atom* atom, int singletonIndex,vector<WClause*> clauses)
	{
		for(unsigned int i=0;i<clauses.size();i++)
		{
			vector<LvrTerm*> terms;
			vector<int> occurences;
			for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
			{
				if(clauses[i]->atoms[j]->symbol->id == atom->symbol->id)
				{
					terms.push_back(clauses[i]->atoms[j]->terms[singletonIndex]);
					occurences.push_back(j);
				}
			}
			if(terms.size()==0)
				continue;
			for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
			{
				if(clauses[i]->atoms[j]->isConstant())
					continue;
				if(find(occurences.begin(),occurences.end(),j)!=occurences.end())
					continue;
				int foundcount=0;
				for(unsigned int k=0;k<clauses[i]->atoms[j]->terms.size();k++)
				{
					if(clauses[i]->atoms[j]->terms[k]->domain.size()==1)
						continue;
					if(find(terms.begin(),terms.end(),clauses[i]->atoms[j]->terms[k])!=terms.end())
						foundcount++;
				}
				if(foundcount!=0 && foundcount!=terms.size())
					return false;
			}
		}
		return true;
	}

private:
	static bool checkOccurrence(LvrTerm* term, WClause* clause, int exceptIndex)
	{
		for(unsigned int i=0;i<clause->atoms.size();i++)
		{
			if(i == exceptIndex)
				continue;
			if(clause->atoms[i]->isConstant())
				continue;
			for(unsigned int j=0;j<clause->atoms[i]->terms.size();j++)
			{
				if(clause->atoms[i]->terms[j]->domain.size() != 1
					&& clause->atoms[i]->terms[j] == term)
				{
					return true;
				}
			}
		}
		return false;
	}

};
#endif
