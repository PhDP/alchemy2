#ifndef __LNORMPRESERVING_GROUNDER
#define __LNORMPRESERVING_GROUNDER

#include "lvrmln.h"
#include "hashalgorithm.h"
#include "lvrpermutations.h"
#include "cleanuputils.h"
#include "randomgenutil.h"

struct LvrNormPreservingGrounder
{
	static void propagateNormalizedCNF(vector<WClause*>& CNF, Atom* atomToSplit,
			vector<bool> isolatedTerms)
	{
		doPropagateNormalizedCNF(CNF,atomToSplit,isolatedTerms);
	}

private:
	static void insertOrUpdate( map<int,vector<bool> >&  sharedTermsMap, Atom* atom, int indexToSet)
	{
		map<int,vector<bool> >::iterator it = sharedTermsMap.find(atom->symbol->id);
		if(it==sharedTermsMap.end())
		{
			vector<bool> value(atom->terms.size());
			value[indexToSet] = true;
			sharedTermsMap.insert(pair<int,vector<bool> >(atom->symbol->id,value));
		}
		else
		{
			sharedTermsMap[atom->symbol->id].at(indexToSet) = true;
		}
	}

	static void updateSharedTerms(WClause* clause, int atomId, map<int,vector<bool> >&  sharedTermsMap,vector<bool>& dirtyAtoms)
	{
		dirtyAtoms.resize(clause->atoms.size());
		vector<int> atomInd;
		vector<vector <LvrTerm*> > terms;
		for(unsigned int j=0;j<clause->atoms.size();j++)
		{
			//if(clause->atoms[j]->isConstant())
				//continue;
			bool allconstants = true;
			if(atomId == clause->atoms[j]->symbol->id)
			{
				vector<LvrTerm*> tmpTerms;
				for(unsigned int k=0;k<clause->atoms[j]->terms.size();k++)
				{
					if(clause->atoms[j]->terms[k]->domain.size() > 1)
					{
						allconstants = false;
						tmpTerms.push_back(clause->atoms[j]->terms[k]);
					}
				}
				if(tmpTerms.size() > 0)
					terms.push_back(tmpTerms);
			}
			if(!allconstants)
				atomInd.push_back(j);
		}

		if(atomInd.size() > 0)
		{
			for(unsigned int j=0;j<clause->atoms.size();j++)
			{
				map<int,vector<bool> >::iterator it = sharedTermsMap.find(clause->atoms[j]->symbol->id);
				if(it==sharedTermsMap.end())
				{
					vector<bool> shterms(clause->atoms[j]->terms.size());
					sharedTermsMap.insert(pair<int,vector<bool> >(clause->atoms[j]->symbol->id,shterms));
					it = sharedTermsMap.find(clause->atoms[j]->symbol->id);
				}
				if(find(atomInd.begin(),atomInd.end(),j)!=atomInd.end())
					continue;
				for(unsigned int k=0;k<clause->atoms[j]->terms.size();k++)
				{
					if(it->second[k])
						continue;
					for(unsigned int m=0;m<atomInd.size();m++)
					{
						if(find(terms[m].begin(),terms[m].end(),clause->atoms[j]->terms[k])!=terms[m].end())
						{
							//k is a shared term
							insertOrUpdate(sharedTermsMap,clause->atoms[j],k);
							dirtyAtoms[j] = true;
							break;
						}
					}
				}
			}
		}
	}


	static void doPropagateNormalizedCNF(vector<WClause*>& CNF, Atom* atomToSplit,vector<bool> isolatedTerms)
	{
		map<int,vector<bool> > sharedTermsMap;
		//shared terms inverse of isolated terms
		vector<bool> shTerms(isolatedTerms.size());
		for(unsigned int i=0;i<isolatedTerms.size();i++)
			shTerms[i] = (!isolatedTerms[i]);
		sharedTermsMap.insert(pair<int,vector<bool> >(atomToSplit->symbol->id,shTerms));
		vector<int> activeAtoms;
		activeAtoms.push_back(atomToSplit->symbol->id);
		while(!activeAtoms.empty())
		{
			int id = activeAtoms.back();
			activeAtoms.pop_back();
			for(unsigned int i=0;i<CNF.size();i++)
			{
				vector<bool> dirtyAtoms;
				updateSharedTerms(CNF[i],id,sharedTermsMap,dirtyAtoms);
				//add all dirty atoms
				for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
				{
					if(dirtyAtoms[j])
						activeAtoms.push_back(CNF[i]->atoms[j]->symbol->id);
				}
			}
		}
		groundCNF(CNF,sharedTermsMap,atomToSplit->symbol->id);
	}

	static bool isGrounded(Atom* atom,vector<bool> sharedTerms)
	{
		for(unsigned int i=0;i<sharedTerms.size();i++)
		{
			if(sharedTerms[i] && atom->terms[i]->domain.size() > 1)
				return false;
		}
		return true;
	}

	static void groundClause(WClause* clause, int atomInd, vector<bool> sharedTerms,vector<WClause*>& newClauses,int atomIdToRemove)
	{
		vector<LvrTerm*> terms;
		for(unsigned int i=0;i<sharedTerms.size();i++)
		{
			if(sharedTerms[i])
				terms.push_back(clause->atoms[atomInd]->terms[i]);
		}
		vector<vector<int> > permutedList;
		LvrPermutations::permuteTerms(terms,permutedList);
		for(unsigned int i=0;i<permutedList.size();i++)
		{
			WClause* nClause = LvrMLN::create_new_clause(clause);
			int iter=0;
			for(unsigned int j=0;j<nClause->atoms[atomInd]->terms.size();j++)
			{
				if(sharedTerms[j])
				{
					nClause->atoms[atomInd]->terms[j]->domain.clear();
					nClause->atoms[atomInd]->terms[j]->domain.push_back(permutedList[i].at(iter++));
				}
			}
			//change the id of the predicate
			nClause->atoms[atomInd]->symbol->id = LvrHashAlgorithm::convertToHash(nClause->atoms[atomInd]);
			newClauses.push_back(nClause);
		}
	}
	
	static void groundCNF(vector<WClause*>& CNF, map<int,vector<bool> > sharedTermsMap,int atomIdToRemove)
	{
		bool changed=false;
		while(1)
		{
			changed=false;
			for(unsigned int i=0;i<CNF.size();i++)
			{
				for(int j=0;j<CNF[i]->atoms.size();j++)
				{
					map<int,vector<bool> >::iterator it = sharedTermsMap.find(CNF[i]->atoms[j]->symbol->id);
					if(it==sharedTermsMap.end())
						continue;
					//check if grounded
					if(!isGrounded(CNF[i]->atoms[j],it->second))
					{
						//ground and add to CNF
						vector<WClause*> newClauses;
						groundClause(CNF[i],j,it->second,newClauses,atomIdToRemove);
						if(newClauses.size()==0)
							continue;
						//removeItem(CNF,i);
						changed = true;
						delete CNF[i];
						CNF[i] = newClauses[0];
						for(unsigned int jj=1;jj<newClauses.size();jj++)
							CNF.push_back(newClauses[jj]);
						break;
					}
				}
				if(changed)
					break;
			}
			if(!changed)
				break;
		}
		for(unsigned int i=0;i<CNF.size();i++)
		{
			for(int j=0;j<CNF[i]->atoms.size();j++)
			{
				CNF[i]->atoms[j]->symbol->id = LvrHashAlgorithm::convertToHash(CNF[i]->atoms[j]);

			}
		}
	}
};

#endif
