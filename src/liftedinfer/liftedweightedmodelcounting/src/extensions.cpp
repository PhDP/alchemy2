#include "extensions.h"
#include "cleanuputils.h"
#include "normalizer.h"
#include "queryupdater.h"

LExtensions::LExtensions(LvrMLN &mln_):mln(mln_),cache(NULL)
{
}

LExtensions::~LExtensions()
{
	if(cache)
		delete cache;
}


LogDouble LExtensions::getUnitClauseWeightWMPTP(WClause* clause)
{
	int groundings = 1;
	for(unsigned int i=0;i<clause->atoms[0]->terms.size();i++)
	{
		groundings = groundings * clause->atoms[0]->terms[i]->domain.size();
	}
	LogDouble val(1,false);
	if(clause->sign[0])
		LogDouble::LDPower(clause->atoms[0]->symbol->nweight,groundings,val);
	else
		LogDouble::LDPower(clause->atoms[0]->symbol->pweight,groundings,val);
	return val;
}


LogDouble LExtensions::unitPropagation(vector<WClause*>& CNF,bool ptpexactmar)
{
	LogDouble val(1,false);
	while(1)
	{
		//maintain a set of current unit clause predicate id's
		vector<int> unitPredicateIds;
		vector<bool> unitPredicateIdSign;
		for(unsigned int i=0;i<CNF.size();i++)
		{
			if(CNF[i]->atoms.size()==1 && !CNF[i]->satisfied)
			{
				if(ptpexactmar)
				{
					//do not unit propagate on query predicates
					if(LvrQueryUpdater::isInstanceCreated())
					{
						if(LvrQueryUpdater::Instance()->isNormIdInQuery(CNF[i]->atoms[0]->symbol->normParentId))
							continue;
					}
				}
				int matchingId = -1;
				for(unsigned int j=0;j<unitPredicateIds.size();j++)
				{
					if(unitPredicateIds[j]==CNF[i]->atoms[0]->symbol->id)
					{
						matchingId=j;
						break;
					}
				}
				if(matchingId != -1)
				{
					//if signs are the same, delete one of them since equivalent
					if(CNF[i]->sign[0] == unitPredicateIdSign[matchingId])
					{
						removeItem(CNF,i);
						i--;
					}
					else
					{
						//conflicting unit clauses can never be satisfied
						return LogDouble(0,false);
					}
				}
				else
				{
					//set unit clause to satisfied
					CNF[i]->satisfied = true;
					//compute weight
					val=val*getUnitClauseWeightWMPTP(CNF[i]);
					unitPredicateIds.push_back(CNF[i]->atoms[0]->symbol->id);
					unitPredicateIdSign.push_back(CNF[i]->sign[0]);
					//remove atom from unit clause and count its weight
					CNF[i]->removeAtom(0);
					i--;
				}
			}
		}
		bool isCNFChanged = false;
		for(unsigned int i=0;i<unitPredicateIds.size();i++)
		{
			//delete all clauses with conflicting sign for the unit atom
			for(unsigned int j=0;j<CNF.size();j++)
			{
				for(unsigned int k=0;k<CNF[j]->atoms.size();k++)
				{
					if(CNF[j]->atoms[k]->symbol->id == unitPredicateIds[i])
					{
						if(CNF[j]->sign[k] != unitPredicateIdSign[i])
						{
							//conflicting atom in clause, remove it
							CNF[j]->removeAtom(k);
							isCNFChanged = true;
							k--;
							//check if we get an empty clause
							if(CNF[j]->atoms.size()==0 && !CNF[j]->satisfied)
								return LogDouble(0,false);
						}
						else
						{
							//set clause to satisfied
							CNF[j]->satisfied = true;
							//satisfied atom in clause, remove it
							CNF[j]->removeAtom(k);
							k--;
						}
					}
				}
			}
		}
		if(!isCNFChanged)
			break;
	}
	return val;
}


bool LExtensions::getFromCache(vector<WClause*> CNF, LogDouble& val)
{
	if(cache)
	{
		vector<WClause*> clauses;
		mln.copyAllClauses(CNF,clauses);
		bool ret = cache->findInCache(clauses,val);
		cleanup(clauses);
		return ret;
	}
	else
		return false;
}

void LExtensions::storeToCache(vector<WClause*> CNF, LogDouble& val)
{
	if(cache)
		cache->addToCache(CNF,val);
}
