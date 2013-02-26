#include "resolver.h"

void LResolver::getMatchingTerms(WClause* clause,Atom* atom, vector<LvrTerm*>& terms)
{
	for(unsigned int i=0;i<clause->atoms.size();i++)
	{
		if(clause->atoms[i]->symbol->id==atom->symbol->id)
		{
			//collect terms if not already collected
			for(unsigned int j=0;j<clause->atoms[i]->terms.size();j++)
			{
				bool found=false;
				for(unsigned int k=0;k<terms.size();k++)
				{
					if(terms[k]==clause->atoms[i]->terms[j])
					{
						found=true;
						break;
					}
				}
				if(!found)
					terms.push_back(clause->atoms[i]->terms[j]);
			}
		}
	}
}

void LResolver::permute(vector<LvrTerm*> terms, vector<vector<int> >& permutedDomainList)
{
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

int LResolver::getIndex(vector<vector<int> > allGroundings, Atom* groundedAtom)
{
	for(unsigned int i=0;i<allGroundings.size();i++)
	{
		bool found=true;
		for(unsigned int j=0;j<(allGroundings[i]).size();j++)
		{
			if((allGroundings[i])[j]!=groundedAtom->terms[j]->domain[0])
			{
				found=false;
				break;
			}
		}
		if(found)
			return i;
	}
	return -1;
}

void LResolver::resolveGroundClause(WClause* clause,int predicateId, vector<int> truthVector, vector<vector<int> > allGroundings, WClause& resolvedClause)
{
	int sizeReduction=0;
	for(unsigned int k=0;k<clause->atoms.size();k++)
	{
		int newClauseIndex=k-sizeReduction;
		if(predicateId==resolvedClause.atoms[newClauseIndex]->symbol->id)
		{
			//get the Index into the correct truthvector value corresponding to the domain
			int index=getIndex(allGroundings,resolvedClause.atoms[newClauseIndex]);

			if((truthVector[index] == 1 && !resolvedClause.sign[newClauseIndex] )||(truthVector[index] == 0 && resolvedClause.sign[newClauseIndex]))
				resolvedClause.satisfied=true;
			resolvedClause.removeAtom(k-sizeReduction);
			sizeReduction++;
		}
	}
}


void LResolver::doFullGrounding(vector<WClause*> CNF,Atom* atom,vector<WClause*>& groundCNF,vector<vector<int> >& allGroundings,int& numGroundings)
{
	for(unsigned int i=0;i<CNF.size();i++)
	{
		vector<LvrTerm*> matchingTerms;
		getMatchingTerms(CNF[i],atom,matchingTerms);
		if(matchingTerms.size()==0)
		{
			WClause* newClause = LvrMLN::create_new_clause(CNF[i]);
			groundCNF.push_back(newClause);
			continue;
		}
		vector<vector<int> > permutedDomainList;
		permute(matchingTerms,permutedDomainList); 
		for(unsigned int j=0;j<permutedDomainList.size();j++)
		{
			vector<int> permutedDomain=permutedDomainList[j];
			WClause* newClause = LvrMLN::create_new_clause(CNF[i]);
			for(unsigned int k=0;k<matchingTerms.size();k++)
			{
				for(unsigned int m=0;m<CNF[i]->atoms.size();m++)
				{
					for(unsigned int n=0;n<CNF[i]->atoms[m]->terms.size();n++)
					{
						//ignore constants
						if(CNF[i]->atoms[m]->terms[n]->domain.size()==0)
							continue;
						if(CNF[i]->atoms[m]->terms[n]==matchingTerms[k])
						{
							//substitute term in new clause with the appropriate domain value
							newClause->atoms[m]->terms[n]->domain.clear();
							newClause->atoms[m]->terms[n]->domain.push_back(permutedDomain[k]);

						}
					}
				}
			}
			groundCNF.push_back(newClause);
		}
	}
	numGroundings=1;
	for(unsigned int i=0;i<atom->terms.size();i++)
	{
		numGroundings *= atom->terms[i]->domain.size();
	}
	
	permute(atom->terms,allGroundings);
}


bool isSingletonShared(WClause* clause, Atom* atom, int singletonIndex)
{
	for(unsigned int i=0;i<clause->atoms.size();i++)
	{
		if(clause->atoms[i]->symbol->id == atom->symbol->id)
			continue;
		for(unsigned int j=0;j<clause->atoms[i]->terms.size();j++)
		{
			if(clause->atoms[i]->terms[j]->domain.size()==1)
				continue;
			if(clause->atoms[i]->terms[j] == atom->terms[singletonIndex])
				return true;
		}
	}
	return false;
}

void LResolver::reduceDomains(WClause* clause, Atom* atom, vector<int> truthValues, vector<WClause*>& domainModifiedClauses, int singletonIndex)
{
	int numTrue=0;
	int numFalse=0;
	for(unsigned int i=0;i<truthValues.size();i++)
	{
		if(truthValues[i]==1)
			numTrue++;
		else
			numFalse++;
	}
	if(numTrue==0 || numFalse==0)
	{
		WClause* nClause = LvrMLN::create_new_clause(clause);
		//simply remove the matching atoms
		int rmcnt=0;
		for(unsigned int i=0;i<nClause->atoms.size();i++)
		{
			if(nClause->atoms[i]->symbol->id == atom->symbol->id)
			{
				if(nClause->sign[i] && numTrue == 0)
					nClause->satisfied=true;
				else if(!nClause->sign[i] && numFalse==0)
					nClause->satisfied=true;
				nClause->removeAtom(i);
				i--;
			}
		}
		domainModifiedClauses.push_back(nClause);
		return;
	}
	bool changed=false;
	domainModifiedClauses.push_back(LvrMLN::create_new_clause(clause));
	while(1)
	{
		changed=false;
		//find the singleton
		for(unsigned int i=0;i<domainModifiedClauses.size();i++)
		{
			for(unsigned int j=0;j<domainModifiedClauses[i]->atoms.size();j++)
			{
				if(domainModifiedClauses[i]->atoms[j]->symbol->id == atom->symbol->id)
				{
					//reduce its domain
					WClause* nClausep = LvrMLN::create_new_clause(domainModifiedClauses[i]);
					WClause* nClausen = LvrMLN::create_new_clause(domainModifiedClauses[i]);
					nClausep->atoms[j]->terms[singletonIndex]->domain.erase(
						nClausep->atoms[j]->terms[singletonIndex]->domain.begin()+numTrue,
						nClausep->atoms[j]->terms[singletonIndex]->domain.end());
					nClausen->atoms[j]->terms[singletonIndex]->domain.erase(
						nClausen->atoms[j]->terms[singletonIndex]->domain.begin(),
						nClausen->atoms[j]->terms[singletonIndex]->domain.begin()+numTrue);
					if(domainModifiedClauses[i]->sign[j])
						nClausen->satisfied=true;
					else
						nClausep->satisfied=true;
					nClausep->removeAtom(j);
					nClausen->removeAtom(j);
					delete domainModifiedClauses[i];
					domainModifiedClauses[i] = nClausep;
					domainModifiedClauses.push_back(nClausen);
					changed=true;
				}
			}
		}
		if(!changed)
			break;
	}
}



void LResolver::reduceDomainsSelfJoins(WClause* clause, Atom* atom, vector<int> truthValues, vector<WClause*>& domainModifiedClauses, int singletonIndex)
{
	vector<int> matchingIds;
	bool posSignFound = false;
	bool negSignFound = false;
	for(unsigned int i=0;i<clause->atoms.size();i++)
	{
		if(clause->atoms[i]->symbol->id==atom->symbol->id)
		{
			matchingIds.push_back(i);
			if(clause->sign[i])
				negSignFound = true;
			else
				posSignFound = true;
		}
	}
	if(matchingIds.size()==0)
	{
		//splitting atom not found, return original clause
		WClause* nClause= LvrMLN::create_new_clause(clause);
		domainModifiedClauses.push_back(nClause);
		return;
	}

	int numOne=0;
	int numZero=0;
	//for a singleton atom, #truthValues = #domainvalues(or size of domain)
	for(unsigned int i=0;i<truthValues.size();i++)
	{
		if(truthValues[i]==1)
		{
			numOne++;
		}
		else
		{
			numZero++;
		}
	}
	if((numOne == 0 && negSignFound) ||
		(numZero == 0 && posSignFound))
	{
		//all are satisfied
		WClause* nClause = LvrMLN::create_new_clause(clause);
		for(unsigned int i=0;i<matchingIds.size();i++)
		{
			nClause->removeAtom(matchingIds[i]-i);
		}
		nClause->satisfied=true;
		domainModifiedClauses.push_back(nClause);
		return;
	}
	int rmCnt = 0;
	WClause* remainingClause = LvrMLN::create_new_clause(clause);
	vector<WClause*> satClauses;
	for(unsigned int i=0;i<matchingIds.size();i++)
	{
		int ind = matchingIds[i]-rmCnt;
		if(!remainingClause->sign[ind])
		{
			if(numOne > 0)
			{
				WClause* posClause = LvrMLN::create_new_clause(remainingClause);
				//+ve atom
				posClause->atoms[ind]->terms[singletonIndex]->domain.erase(posClause->atoms[ind]->terms[singletonIndex]->domain.begin()
					+ numOne,posClause->atoms[ind]->terms[singletonIndex]->domain.end());
				posClause->satisfied = true;
				posClause->removeAtom(ind);
				satClauses.push_back(posClause);
				//posClause->print();
			}
			remainingClause->atoms[ind]->terms[singletonIndex]->domain.erase(remainingClause->atoms[ind]->terms[singletonIndex]->domain.begin(),
				remainingClause->atoms[ind]->terms[singletonIndex]->domain.begin()+numOne);
			remainingClause->removeAtom(ind);
			rmCnt++;
		}
		else if(remainingClause->sign[ind])
		{
			if(numZero > 0)
			{
				WClause* posClause = LvrMLN::create_new_clause(remainingClause);
				//+ve atom
				//posClause->atoms[ind]->terms[singletonIndex]->domain.erase(posClause->atoms[ind]->terms[singletonIndex]->domain.begin()
					//+ numZero,posClause->atoms[ind]->terms[singletonIndex]->domain.end());
				posClause->atoms[ind]->terms[singletonIndex]->domain.erase(posClause->atoms[ind]->terms[singletonIndex]->domain.begin()
					,posClause->atoms[ind]->terms[singletonIndex]->domain.end() - numZero);
				posClause->satisfied = true;
				posClause->removeAtom(ind);
				satClauses.push_back(posClause);
				//posClause->print();
			}
			//remainingClause->atoms[ind]->terms[singletonIndex]->domain.erase(remainingClause->atoms[ind]->terms[singletonIndex]->domain.begin(),
				//remainingClause->atoms[ind]->terms[singletonIndex]->domain.begin()+numZero);
			remainingClause->atoms[ind]->terms[singletonIndex]->domain.erase(remainingClause->atoms[ind]->terms[singletonIndex]->domain.end()-numZero,
			remainingClause->atoms[ind]->terms[singletonIndex]->domain.end());
			remainingClause->removeAtom(ind);
			rmCnt++;
		}
	}
	for(unsigned int t=0;t<satClauses.size();t++)
	{
		for(unsigned int i=0;i<satClauses[t]->atoms.size();i++)
		{
			if(satClauses[t]->atoms[i]->symbol->id == atom->symbol->id)
			{
				satClauses[t]->removeAtom(i);
				i--;
			}
		}
		domainModifiedClauses.push_back(satClauses[t]);
	}
	domainModifiedClauses.push_back(remainingClause);
	//for(unsigned int i=0;i<domainModifiedClauses.size();i++)
	//{
		//domainModifiedClauses[i]->print();
	//}
	//remainingClause->print();
}
