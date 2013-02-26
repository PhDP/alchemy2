#ifndef __LVRSINGLETON_NORM_PROPAGATION
#define __LVRSINGLETON_NORM_PROPAGATION
#include "lvrmln.h"
#include "hashalgorithm.h"

struct LvrSingletonNormPropagation
{

	static void propagateNormalizedCNF(vector<WClause*>& CNF, Atom* atomToSplit,int singletonIndex, int numTrue)
	{
		int numFalse = atomToSplit->getNumberOfGroundings() - numTrue;
		if(numFalse==0)
		{
			//simply remove all the matching predicates and set to satisfied
			for(unsigned int i=0;i<CNF.size();i++)
			{
				for(int j=0;j<CNF[i]->atoms.size();j++)
				{
					if(CNF[i]->atoms[j]->symbol->id == atomToSplit->symbol->id)
					{
						if(!CNF[i]->sign[j])
						{
							CNF[i]->satisfied=true;
						}
						CNF[i]->removeAtomWithWeightUpdation(j);
						j--;
					}
				}
			}
			return;
		}
		else if(numTrue==0)
		{
			//simply remove all the matching predicates and set to satisfied
			for(unsigned int i=0;i<CNF.size();i++)
			{
				for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
				{
					if(CNF[i]->atoms[j]->symbol->id == atomToSplit->symbol->id)
					{
						if(CNF[i]->sign[j])
						{
							CNF[i]->satisfied=true;
						}
						CNF[i]->removeAtomWithWeightUpdation(j);
						j--;
					}
				}
			}
			return;
		}

		vector<int> tmp;
		tmp.push_back(singletonIndex);
		map<int,vector<int> > domainReductionTerms;
		bool changed = false;
		domainReductionTerms.insert(pair<int,vector<int> >(atomToSplit->symbol->id,tmp));
		while(1)
		{
			changed = false;
			//find all domain reduce positions
			for(map<int,vector<int> >::iterator it=domainReductionTerms.begin();it!=domainReductionTerms.end();it++)
			{
				for(unsigned int j=0;j<CNF.size();j++)
				{
					vector<int> foundIds;
					vector<vector<LvrTerm*> > termsToCheck;
					for(unsigned int k=0;k<CNF[j]->atoms.size();k++)
					{
						if(it->first == CNF[j]->atoms[k]->symbol->id)
						{
							//collect the terms
							foundIds.push_back(k);
							vector<LvrTerm*> terms;
							for(unsigned int m=0;m<it->second.size();m++)
							{
								terms.push_back(CNF[j]->atoms[k]->terms[it->second.at(m)]);
							}
							termsToCheck.push_back(terms);
						}
					}
					for(unsigned int k=0;k<foundIds.size();k++)
					{
						for(unsigned int m=0;m<CNF[j]->atoms.size();m++)
						{
							if(m==foundIds[k])
								continue;
							for(unsigned int jj=0;jj<CNF[j]->atoms[m]->terms.size();jj++)
							{
								if(find(termsToCheck[k].begin(),termsToCheck[k].end(),CNF[j]->atoms[m]->terms[jj])!=
									termsToCheck[k].end())
								{
									//need to add to domain reduction terms
									map<int,vector<int> >::iterator it1 = domainReductionTerms.find(CNF[j]->atoms[m]->symbol->id);
									if(it1!=domainReductionTerms.end())
									{
										if(find(it1->second.begin(),it1->second.end(),jj)==it1->second.end())
										{
											changed = true;
											it1->second.push_back(jj);
										}
									}
									else
									{
										changed = true;
										vector<int> tmp(1);
										tmp[0]=jj;
										domainReductionTerms.insert(pair<int,vector<int> >(CNF[j]->atoms[m]->symbol->id,tmp));
									}

								}
							}
						}
					}
				}
			}
			if(!changed)
				break;
		}
		vector<WClause*> newCNF;
		for(unsigned int i=0;i<CNF.size();i++)
		{
			if(CNF[i]->isPropositional())
			{
				//just change the ids of the atoms
				for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
					CNF[i]->atoms[j]->symbol->id = LvrHashAlgorithm::convertToHash(CNF[i]->atoms[j]);
				newCNF.push_back(CNF[i]);
				continue;
			}
			vector<WClause*> outClauses;
			splitClause(CNF[i],atomToSplit->symbol->id,domainReductionTerms,numTrue,numFalse,outClauses);
			for(unsigned int j=0;j<outClauses.size();j++)
				newCNF.push_back(outClauses[j]);
		}
		CNF.clear();
		CNF=newCNF;
	}

	static void splitClause(WClause* clause,int atomToRemove,map<int,vector<int> > domainReduceMap,
		int numTrue,int numFalse, vector<WClause*>& outClauses)
	{
		outClauses.push_back(clause);
		bool changed=false;
		while(1)
		{
			changed=false;
			for(unsigned int t=0;t<outClauses.size();t++)
			{
				WClause* clause1 = outClauses[t];
				if(clause1->isPropositional())
					continue;
				vector<int> positionsToRemove;
				vector<int> otherpositions;
				for(unsigned int i=0;i<clause1->atoms.size();i++)
				{
					if(clause1->atoms[i]->symbol->id==atomToRemove)
					{
						positionsToRemove.push_back(i);
					}
					else
						otherpositions.push_back(i);
				}
				while(!positionsToRemove.empty() || !otherpositions.empty())
				{
					int i = -1;
					if(!positionsToRemove.empty())
					{
						i=positionsToRemove.back();
						positionsToRemove.pop_back();
					}
					else
					{
						i=otherpositions.back();
						otherpositions.pop_back();
					}
					if(clause1->atoms[i]->isConstant())
						continue;
					//choose the atom to be removed, if it exists before others

					map<int,vector<int> >::iterator it = domainReduceMap.find(clause1->atoms[i]->symbol->id);
					if(it==domainReduceMap.end() || it->second.size()==0)
						continue;
					//check if there exists a candidate for domain reduction
					int domainReduceIndex = -1;

					for(unsigned int j=0;j<it->second.size();j++)
					{
						if(clause1->atoms[i]->terms[it->second[j]]->domain.size() != numTrue
							&& clause1->atoms[i]->terms[it->second[j]]->domain.size() != numFalse)
						{
							domainReduceIndex=it->second[j];
							break;
						}
					}
					if(domainReduceIndex==-1)
						continue;
					if(clause1->atoms[i]->terms[domainReduceIndex]->domain.size()!=(numTrue+numFalse))
					{
						//maybe a collision the hash algorithm
						continue;
					}
					//can simply reduce domain
					WClause* nClausep = LvrMLN::create_new_clause(clause1);
					WClause* nClausen = LvrMLN::create_new_clause(clause1);
					nClausep->atoms[i]->terms[domainReduceIndex]->domain.erase
						(nClausep->atoms[i]->terms[domainReduceIndex]->domain.begin()+numTrue,
						nClausep->atoms[i]->terms[domainReduceIndex]->domain.end());
					nClausen->atoms[i]->terms[domainReduceIndex]->domain.erase(
						nClausen->atoms[i]->terms[domainReduceIndex]->domain.begin(),
						nClausen->atoms[i]->terms[domainReduceIndex]->domain.begin()+numTrue);
					if(clause1->atoms[i]->symbol->id == atomToRemove)
					{
						if(clause1->sign[i])
							nClausen->satisfied=true;
						else
							nClausep->satisfied=true;
						nClausep->removeAtomWithWeightUpdation(i);
						nClausen->removeAtomWithWeightUpdation(i);
					}
					delete clause1;
					outClauses[t]=nClausep;
					outClauses.push_back(nClausen);
					changed=true;
					break;
				}
			}
			if(!changed)
				break;
		}

		for(unsigned int i=0;i<outClauses.size();i++)
		{
			for(int j=0;j<outClauses[i]->atoms.size();j++)
			{
				if(outClauses[i]->atoms[j]->symbol->id == atomToRemove)
				{
					outClauses[i]->removeAtomWithWeightUpdation(j);
					j--;
				}
				else
				{
					//change id
					outClauses[i]->atoms[j]->symbol->id = LvrHashAlgorithm::convertToHash(outClauses[i]->atoms[j]);
				}
			}
		}
	}
};
#endif
