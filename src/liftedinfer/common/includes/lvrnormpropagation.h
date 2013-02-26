#ifndef __LVRNORM_PROPAGATION
#define __LVRNORM_PROPAGATION
#include "lvrmln.h"
#include "hashalgorithm.h"
#include "lvrpermutations.h"
#include "cleanuputils.h"
#include "randomgenutil.h"

struct LvrNormPropagation
{
	LvrMLN& mln;
	set<int> hashCodesToSetTrue;
	map<int,int> isolatedHashCodeMap;
	enum ERules{AGB,IV,NONE};
	ERules ruleToApply;
	vector<bool>isolatedTerms;
	LvrNormPropagation(LvrMLN& mln_):mln(mln_){}

	void propagateNormalizedCNF(vector<WClause*>& CNF, Atom* atomToSplit,vector<bool> isolatedTerms_)
	{
		if(atomToSplit->isConstant())
		{
			//efficient to split
			processConstants(CNF,atomToSplit,false);
		}
		else
		{
			isolatedTerms = isolatedTerms_;
			ruleToApply = NONE;
			doPropagateNormalizedCNF(CNF,atomToSplit);
		}
	}
	/*
	void processConstants(vector<WClause*>& CNF, Atom* atom, bool setTrue)
	{
		unsigned int hash = LvrHashAlgorithm::convertToHash(atom);
		for(int i=0;i<CNF.size();i++)
		{
			for(int j=0;j<CNF[i]->atoms.size();j++)
			{
				unsigned int hash1 = LvrHashAlgorithm::convertToHash(CNF[i]->atoms[j]);
				if(hash == hash1)
				{
					if(CNF[i]->sign[j] && !setTrue)
						CNF[i]->satisfied= true;
					else if(!CNF[i]->sign[j] && setTrue)
						CNF[i]->satisfied= true;
					CNF[i]->removeAtom(j);
					j--;
				}
			}
			if(CNF[i]->atoms.size()==0 && !CNF[i]->satisfied)
			{
				removeItem(CNF,i);
				i--;
			}
		}
		
	}
	*/
	void processConstants(vector<WClause*>& CNF, Atom* atom, bool setTrue)
	{
		for(int i=0;i<CNF.size();i++)
		{
			for(int j=0;j<CNF[i]->atoms.size();j++)
			{
				if(CNF[i]->atoms[j]->symbol->id == atom->symbol->id)
				{
					if(CNF[i]->sign[j] && !setTrue)
						CNF[i]->satisfied= true;
					else if(!CNF[i]->sign[j] && setTrue)
						CNF[i]->satisfied= true;
					CNF[i]->removeAtom(j);
					j--;
				}
			}
			if(CNF[i]->atoms.size()==0 && !CNF[i]->satisfied)
			{
				removeItem(CNF,i);
				i--;
			}
		}
		
	}

	void propagateNormalizedCNF(vector<WClause*>& CNF, Atom* atomToSplit,vector<bool> isolatedTerms_,vector<int> numGroundingsToSetTrue)
	{
		isolatedTerms = isolatedTerms_;
		hashCodesToSetTrue.clear();
		isolatedHashCodeMap.clear();
		//determine the hashCodes of groundings that are to be set to true
		if(numGroundingsToSetTrue.size()==1)
		{
			//check for constants that can be removed efficiently
			if(atomToSplit->isConstant())
			{
				if(numGroundingsToSetTrue[0]==0)
				{
					processConstants(CNF,atomToSplit,false);
					return;
				}
				else if(numGroundingsToSetTrue[0]==atomToSplit->getNumberOfGroundings())
				{
					processConstants(CNF,atomToSplit,true);
					return;
				}
			}
			ruleToApply = AGB;
			//no isolated variables, find hash codes of random assignments
			vector<int> indexesToSetTrue;
			int countInd = 0;
			while(countInd!=numGroundingsToSetTrue[0])
			{
				int ind = LvRandomGenUtil::Instance()->getRandomPosition(atomToSplit->getNumberOfGroundings());
				if(find(indexesToSetTrue.begin(),indexesToSetTrue.end(),ind)==indexesToSetTrue.end())
				{
					indexesToSetTrue.push_back(ind);
					countInd++;
				}
			}
			/*
			vector<vector<int> > permutedList;
			LvrPermutations::permuteTerms(atomToSplit->terms,permutedList);

			for(unsigned int i=0;i<indexesToSetTrue.size();i++)
			{
				vector<int> randGrounding;
				randGrounding.push_back(atomToSplit->symbol->parentId);
				for(unsigned int j=0;j<permutedList[indexesToSetTrue[i]].size();j++)
					randGrounding.push_back(permutedList[indexesToSetTrue[i]].at(j));
				hashCodesToSetTrue.insert(LvrHashAlgorithm::DJBHash(randGrounding));
			}
			*/
			//find the addresses of the indexes to set to true
			vector<int> domSizes(atomToSplit->terms.size());
			for(unsigned int i=0;i<atomToSplit->terms.size();i++)
			{
				domSizes[i] = atomToSplit->terms[i]->domain.size();
			}
			for(unsigned int i=0;i<indexesToSetTrue.size();i++)
			{
				int num = indexesToSetTrue[i];
				vector<int> randGrounding(domSizes.size()+1);
				randGrounding[0] = atomToSplit->symbol->normParentId;
				for(int j=domSizes.size()-1;j>= 0;j--)
				{
					if(domSizes[j]!=1)
					{
					int val = num % domSizes[j];
					num = num/domSizes[j];
					randGrounding[j+1] = atomToSplit->terms[j]->domain[val];
					}
					else
						randGrounding[j+1] = atomToSplit->terms[j]->domain[0];
				}
				hashCodesToSetTrue.insert(LvrHashAlgorithm::DJBHash(randGrounding));
			}

		}
		else
		{
			ruleToApply = IV;
			//contains isolated terms, hence for each grounding to isolated terms, set appropriate numtrue
			vector<LvrTerm*> nonIsolatedTerms;
			for(unsigned int i=0;i<isolatedTerms.size();i++)
			{
				if(!isolatedTerms[i])
					nonIsolatedTerms.push_back(atomToSplit->terms[i]);
			}
			vector<vector<int> > permutedList;
			LvrPermutations::permuteTerms(nonIsolatedTerms,permutedList);
			vector<vector<int> > signature;
			vector<int> id(1);id[0] = atomToSplit->symbol->normParentId;
			signature.push_back(id);
			//signature.push_back(atomToSplit->symbol->normParentId);
			for(unsigned int i=0;i<permutedList.size();i++)
			{
				int iter=0;
				for(unsigned int j=0;j<isolatedTerms.size();j++)
				{
					if(isolatedTerms[j])
					{
						//signature.push_back(atomToSplit->terms[j]->domain.size());
						vector<int> doms;
						for(unsigned int k=0;k<atomToSplit->terms[j]->domain.size();k++)
						{
							doms.push_back(atomToSplit->terms[j]->domain[k]);
						}
						signature.push_back(doms);
					}
					else
					{
						vector<int> dm(1);dm[0] = permutedList[i].at(iter++);
						signature.push_back(dm);
						//signature.push_back(permutedList[i].at(iter++));
					}
				}
				isolatedHashCodeMap.insert(pair<int,int>(LvrHashAlgorithm::DJBLHash(signature),numGroundingsToSetTrue[i]));
			}
		}
		doPropagateNormalizedCNF(CNF,atomToSplit);
	}


void insertOrUpdate( map<int,vector<bool> >&  sharedTermsMap, Atom* atom, int indexToSet)
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

void updateSharedTerms(WClause* clause, int atomId, map<int,vector<bool> >&  sharedTermsMap,vector<bool>& dirtyAtoms)
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


	void doPropagateNormalizedCNF(vector<WClause*>& CNF, Atom* atomToSplit)
	{
		map<int,vector<bool> > sharedTermsMap;
		//sharedTermsMap.insert(pair<int,vector<bool> >(atomToSplit->symbol->id,isolatedTerms));
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

	bool isGrounded(Atom* atom,vector<bool> sharedTerms)
	{
		for(unsigned int i=0;i<sharedTerms.size();i++)
		{
			if(sharedTerms[i] && atom->terms[i]->domain.size() > 1)
				return false;
		}
		return true;
	}

	void groundClause(WClause* clause, int atomInd, vector<bool> sharedTerms,vector<WClause*>& newClauses,int atomIdToRemove)
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
			if(nClause->atoms[atomInd]->symbol->id == atomIdToRemove)
			{
				if(ruleToApply == AGB)
				{
					//check if we need to set some groundings to true
					int hash = LvrHashAlgorithm::convertToHash(nClause->atoms[atomInd]);
					bool found = (hashCodesToSetTrue.find(hash) != hashCodesToSetTrue.end());
					if( (found && !nClause->sign[atomInd])
						|| (!found && nClause->sign[atomInd]))
					{
						nClause->satisfied = true;
					}
					nClause->removeAtom(atomInd);
					if(nClause->satisfied || nClause->atoms.size()!=0)
						newClauses.push_back(nClause);
					else
						delete nClause;
				}
				else if(ruleToApply == IV)
				{
					int hash = LvrHashAlgorithm::convertToHash(nClause->atoms[atomInd]);
					map<int,int>::iterator it = isolatedHashCodeMap.find(hash);
					if(it!=isolatedHashCodeMap.end())
					{
						WClause* nClause_copy = LvrMLN::create_new_clause(nClause);
						LogDouble pwtChange(it->second,false);
						LogDouble nwtChange(nClause->atoms[atomInd]->getNumberOfGroundings() - it->second,false);
						//update the weight
						if(!pwtChange.is_zero)
						{
							nClause->weight = pwtChange*nClause->weight;
							if(!nClause->sign[atomInd])
								nClause->satisfied=true;
							nClause->removeAtom(atomInd);
							if(nClause->satisfied || nClause->atoms.size()!=0)
								newClauses.push_back(nClause);
							else
								delete nClause;

						}
						else
						{
							delete nClause;
						}
						if(!nwtChange.is_zero)
						{
							nClause_copy->weight = nwtChange*nClause_copy->weight;
							if(nClause_copy->sign[atomInd])
								nClause_copy->satisfied = true;
							nClause_copy->removeAtom(atomInd);
							if(nClause->satisfied || nClause->atoms.size()!=0)
								newClauses.push_back(nClause_copy);
							else
								delete nClause;
						}
						else
						{
							delete nClause_copy;
						}
					}
				}
				else
				{
					nClause->removeAtom(atomInd);
					newClauses.push_back(nClause);
				}
			}
			else
			{
				//change the id of the predicate
				nClause->atoms[atomInd]->symbol->id = LvrHashAlgorithm::convertToHash(nClause->atoms[atomInd]);
				newClauses.push_back(nClause);
			}
			
		}
	}

	void groundCNF(vector<WClause*>& CNF, map<int,vector<bool> > sharedTermsMap,int atomIdToRemove)
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
					else if(CNF[i]->atoms[j]->symbol->id == atomIdToRemove)
					{
						CNF[i]->removeAtom(j);
						j--;
					}
				}
				if(changed)
					break;
			}
			if(!changed)
				break;
		}
		//for(unsigned int i=0;i<CNF.size();i++)
			//CNF[i]->print();
		for(unsigned int i=0;i<CNF.size();i++)
		{
			for(int j=0;j<CNF[i]->atoms.size();j++)
			{
				if(CNF[i]->atoms[j]->symbol->id == atomIdToRemove)
				{
					CNF[i]->removeAtom(j);
					j--;
				}
				else
				{
					CNF[i]->atoms[j]->symbol->id = LvrHashAlgorithm::convertToHash(CNF[i]->atoms[j]);
				}
			}
		}
	}

};
#endif
