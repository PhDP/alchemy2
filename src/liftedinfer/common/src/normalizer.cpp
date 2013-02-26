#include "normalizer.h"
#include "cleanuputils.h"
#include "setutils.h"
#include "lvrpermutations.h"
#include "hashalgorithm.h"
#include "rulesutil.h"
#include "normpreservinggrounder.h"
#include <assert.h>
using namespace std;

bool LNormalizer::checkEquivalency(Atom* at1,Atom* at2)
{
	if(at1->symbol->symbol.compare(at2->symbol->symbol)!=0)
		return false;
	if(at1->terms.size() != at2->terms.size())
		return false;
	//#terms in at1 == #terms in at2
	for(unsigned int i=0;i<at1->terms.size();i++)
	{
		if(at1->terms[i]->domain.size()!=at2->terms[i]->domain.size())
			return false;
		else if(at1->terms[i]->domain.size() == 1 && at2->terms[i]->domain.size() == 1)
		{
			if(at1->terms[i]->domain[0]!=at2->terms[i]->domain[0])
				return false;
		}
		else if(!includes(at1->terms[i]->domain.begin(),at1->terms[i]->domain.end(),
			at2->terms[i]->domain.begin(),at2->terms[i]->domain.end()))
		{
			return false;
		}
	}

	return true;
}

void LNormalizer::mergeIdenticalClauses(vector<WClause*>& CNF)
{
	for(unsigned int i=0;i<CNF.size();i++)
	{
		int matchingIndex=-1;
		for(unsigned int j=i+1;j<CNF.size();j++)
		{
			if(CNF[i]->atoms.size()!=CNF[j]->atoms.size())
				continue;
			if((CNF[i]->atoms.size() == 0 && CNF[j]->atoms.size()==0) ||
				compareClauses(CNF[i],CNF[j]))
			{
				matchingIndex = j;
				break;
			}
		}
		if(matchingIndex!=-1)
		{
			//check if any of the 2 clauses are satisfied
			if(CNF[i]->satisfied == CNF[matchingIndex]->satisfied)
			{
				LogDouble tmpWt = CNF[i]->weight;
				double val = exp(tmpWt.value);
				LogDouble valL(val,true);
				double val1 = exp( CNF[matchingIndex]->weight.value);
				LogDouble valL1(val1,true);
				

				WClause* temp=CNF[i];
				CNF[i]=CNF[matchingIndex];
				CNF[matchingIndex]=temp;
				if(!tmpWt.is_zero)
				{
					CNF[matchingIndex]->weight = LogDouble(val+val1,false);
				}
				//delete ith clause
				removeItem(CNF,i);
				i--;
			}
		}
	}
}

void LNormalizer::normalizeClausesWithoutJoins(vector<WClause*>& clauses, bool weightsOnAtoms,bool resetIds)
{
	vector<PredicateSymbol*> constDomainSymbols;
	vector<bool> completedIds(mln.getMaxPredicateId());
	for(int j=0;j<clauses.size();j++)
	{
		for(int k=0;k<clauses[j]->atoms.size();k++)
		{
			//constDomainSymbols.insert(MLN::create_new_symbol(clauses[j]->atoms[k]->symbol));
			if(!completedIds[clauses[j]->atoms[k]->symbol->id])
			{
				completedIds[clauses[j]->atoms[k]->symbol->id] = true;
				constDomainSymbols.push_back(LvrMLN::create_new_symbol(clauses[j]->atoms[k]->symbol));
			}
		}
	}

	while(1)
	{
		int clauseCount;
		int normClauseCount;
		clauseCount=clauses.size();
		//for(set<PredicateSymbol*>::iterator it=constDomainSymbols.begin();it!=constDomainSymbols.end();it++)
		for(vector<PredicateSymbol*>::iterator it=constDomainSymbols.begin();it!=constDomainSymbols.end();it++)
		{
			toNormalForm(*it,clauses);
		}
		normClauseCount=clauses.size();
		if(clauseCount==normClauseCount)
			break;
	}
	cleanup(constDomainSymbols);
#ifdef __DEBUG_PRINT1__
	LvrMLN::print(clauses,"NORM-STEP");
#endif
	renamePredicates(clauses,resetIds);
#ifdef __DEBUG_PRINT1__
	LvrMLN::print(clauses,"RENAME-STEP");
#endif
	if(weightsOnAtoms)
		mergeIdenticalClausesPTP(clauses);
	else
		mergeIdenticalClauses(clauses);
#ifdef __DEBUG_PRINT1__
	LvrMLN::print(clauses,"MERGE-STEP");
#endif

}


void LNormalizer::normalizeClauses(vector<WClause*>& clauses,bool weightsOnAtoms,bool resetIds)
{
	convertToNormalForm(clauses);
#ifdef __DEBUG_PRINT1__
	LvrMLN::print(clauses,"NORM-STEP");
#endif
	renamePredicates(clauses,resetIds);
#ifdef __DEBUG_PRINT1__
	LvrMLN::print(clauses,"RENAME-STEP");
#endif
	
	if(weightsOnAtoms)
		mergeIdenticalClausesPTP(clauses);
	else
		mergeIdenticalClauses(clauses);
		
#ifdef __DEBUG_PRINT1__
	LvrMLN::print(clauses,"MERGE-STEP");
#endif
}

bool LNormalizer::compareClauses(WClause* c1,WClause* c2)
{
	assert(c1->atoms.size()==c2->atoms.size());
	bool c1prop = c1->isPropositional();
	bool c2prop = c2->isPropositional();
	if(c1prop!=c2prop)
		return false;
	if(!c1prop)
	{
		return false;
	}
	else
	{
		vector<bool> found(c1->atoms.size());
		for(unsigned int i=0;i<c1->atoms.size();i++)
		{
			bool equFound=false;
			for(unsigned int j=0;j<c2->atoms.size();j++)
			{
				if(found[j])
					continue;
				if(c1->sign[i]==c2->sign[j] && checkEquivalency(c1->atoms[i],c2->atoms[j]))
				{
					found[j]=true;
					equFound = true;
					break;
				}
			}
			if(!equFound)
				return false;
		}
	}
	return true;
}

void LNormalizer::mergeIdenticalClausesPTP(vector<WClause*>& CNF)
{
	for(int i=0;i<CNF.size();i++)
	{
		int matchingIndex=-1;
		for(int j=i+1;j<CNF.size();j++)
		{
			if(CNF[i]->atoms.size()!=CNF[j]->atoms.size())
				continue;
			if((CNF[i]->atoms.size() == 0 && CNF[j]->atoms.size()==0) ||
				compareClauses(CNF[i],CNF[j]))
			{
				matchingIndex = j;
				break;
			}
		}
		if(matchingIndex!=-1)
		{
			//check if any of the 2 clauses are satisfied
			//CNF[i]->print();
			//CNF[matchingIndex]->print();
			if(CNF[i]->satisfied && CNF[matchingIndex]->satisfied)
			{
				//delete ith clause
				removeItem(CNF,i);
			}
			else if(CNF[i]->satisfied && !CNF[matchingIndex]->satisfied)
			{
				//delete ith clause
				removeItem(CNF,i);
			}
			else if(!CNF[i]->satisfied && CNF[matchingIndex]->satisfied)
			{
				//swap clauses i and matchingIndex
				WClause* temp=CNF[i];
				CNF[i]=CNF[matchingIndex];
				CNF[matchingIndex]=temp;
				//delete ith clause
				removeItem(CNF,i);
			}
			else
			{
				//delete clause i
				removeItem(CNF,i);
			}
			i--;
		}
	}
}

bool LNormalizer::compareClauses(vector<WClause*> c1,vector<WClause*> c2)
{
	if(c1.size()!=c2.size())
		return false;
	for(unsigned int i=0;i<c1.size();i++)
	{
		c1[i]->print();
		c2[i]->print();
		if(c1[i]->atoms.size()!=c2[i]->atoms.size())
			return false;
		for(unsigned int j=0;j<c1[i]->atoms.size();j++)
		{
			if(c1[i]->atoms[j]->symbol->id!=c2[i]->atoms[j]->symbol->id)
			{
				return false;
			}
			if(c1[i]->atoms[j]->terms.size()!=c2[i]->atoms[j]->terms.size())
			{
				return false;
			}
			for(unsigned int k=0;k<c1[i]->atoms[j]->terms.size();k++)
			{
				if(c1[i]->atoms[j]->terms[k]->domain.size()!=c2[i]->atoms[j]->terms[k]->domain.size())
					return false;
				if(c1[i]->atoms[j]->terms[k]->domain.size()==1)
				{
					if(c1[i]->atoms[j]->terms[k]->domain[0]!=c2[i]->atoms[j]->terms[k]->domain[0])
						return false;
				}
			}
		}
	}
	return true;
}



void LNormalizer::renamePredicates(vector<WClause*>& CNF,bool resetIds)
{
	int curr_predicate_id=mln.symbols.size();
	vector<vector<Atom*> > equivAtoms;
	if(resetIds)
		curr_predicate_id = 0;
    for(unsigned int i=0;i<CNF.size();i++)
    {
		set<int> completedAtoms;
        for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
        {
			bool bucketFound=false;
			
			
			for(unsigned int k=0;k<equivAtoms.size();k++)
			{
				vector<Atom*> atomList=equivAtoms[k];
				if(atomList.size()>0)
				{
					if(checkEquivalency(atomList[0],CNF[i]->atoms[j]))
					{
						equivAtoms[k].push_back(CNF[i]->atoms[j]);
						bucketFound=true;
						break;
					}
				}
			}
			
			if(!bucketFound)
			{
				//push into new bucket
				vector<Atom*> newValue;
				newValue.push_back(CNF[i]->atoms[j]);
				equivAtoms.push_back(newValue);
			}
		}
	}
	//for each bucket choose a unique predicate symbol

	for(unsigned int k=0;k<equivAtoms.size();k++)
	{
		//create a new symbol based on the symbols in the bucket
		PredicateSymbol* ps = new PredicateSymbol(curr_predicate_id++,
			(equivAtoms[k])[0]->symbol->symbol, (equivAtoms[k])[0]->symbol->variable_types,
			(equivAtoms[k])[0]->symbol->pweight,
			(equivAtoms[k])[0]->symbol->nweight,(equivAtoms[k])[0]->symbol->normParentId);
		
		if(!resetIds)
			ps->parentId = (equivAtoms[k])[0]->symbol->parentId;
		else
			ps->parentId = ps->id;

		//replace with new symbol
		for(unsigned int m=0;m<equivAtoms[k].size();m++)
		{
			delete (equivAtoms[k])[m]->symbol;
			(equivAtoms[k])[m]->symbol=LvrMLN::create_new_symbol(ps);
		}
		delete ps;
		
	}
	mln.setMaxPredicateId(curr_predicate_id);
}




void LNormalizer::findNPSelfJoinedAtomToGround(vector<WClause*>& CNF,Atom*& atom,vector<bool>& isolatedTerms)
{
	atom = NULL;
	//for all non-poly self joins ground
	for(unsigned int i=0;i<CNF.size();i++)
	{
		//find a candidate to ground self-join (i.e. a non poly self join)
		map<int,vector<int> > selfJoins;
		for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
		{
			map<int,vector<int> >::iterator it = selfJoins.find(CNF[i]->atoms[j]->symbol->id);
			if(it==selfJoins.end())
			{
				vector<int> tmp(1);
				tmp[0] = j;
				selfJoins.insert(pair<int,vector<int> >(CNF[i]->atoms[j]->symbol->id,tmp));
			}
			else
			{
				it->second.push_back(j);
			}
		}
		WClause* clause  = CNF[i];
		int rmcnt = 0;
		for(map<int,vector<int> >::iterator it = selfJoins.begin();it!=selfJoins.end();it++)
		{
			if(it->second.size() <= 1)
				continue;
			int singletonIndex;
			bool singleton = clause->atoms[it->second[0]-rmcnt]->isSingletonAtom(singletonIndex);
			//check if domain is a constant
			if(clause->atoms[it->second[0]-rmcnt]->isConstant())
			{
				//process self joined constants trivially
				//remove all but one
				int totalSize = it->second.size();
				int signCnt = clause->sign[it->second[0]]?0:1;
				int rmcnt1 = 0;
				for(unsigned int j=1;j<it->second.size();j++)
				{
					signCnt += clause->sign[it->second[j]-rmcnt]?0:1;
					clause->removeAtom(it->second[j]-rmcnt);
					rmcnt1++;
				}
				//check if there are atleast 2 opposite signs
				if(signCnt!=totalSize && signCnt!=0)
				{
					clause->removeAtom(it->second[0]-rmcnt);
					clause->satisfied = true;
					rmcnt++;
				}
				rmcnt = rmcnt1;
			}
			else
			{
				if(singleton)
				{
					//if blocked
					bool blocked = LRulesUtil::isBlocked(clause->atoms[it->second[0]],singletonIndex,CNF);
					//bool blocked = false;
					if(!blocked)
					{
						//ground
						isolatedTerms.resize(clause->atoms[it->second[0]-rmcnt]->terms.size());
						atom = clause->atoms[it->second[0]-rmcnt];
						return;
					}
				}
				else
				{
					//ground
						
					LRulesUtil::computeIsolatedTerms(clause->atoms[it->second[0]-rmcnt],CNF,isolatedTerms);
					atom = clause->atoms[it->second[0]-rmcnt];
					return;
				}
			}
		}
	}
}

void LNormalizer::removeNPSelfJoins(vector<WClause*>& CNF)
{
	while(1)
	{
		Atom* atom=NULL;
		vector<bool> isolatedTerms;
		findNPSelfJoinedAtomToGround(CNF,atom,isolatedTerms);
		if(atom==NULL)
			break;
		LvrNormPreservingGrounder::propagateNormalizedCNF(CNF,atom,isolatedTerms);
	}
	for(unsigned int i=0;i<CNF.size();i++)
	{
		for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
		{
			CNF[i]->atoms[j]->symbol->id = LvrHashAlgorithm::convertToHash(CNF[i]->atoms[j]);
		}
	}
}




void LNormalizer::convertToNormalForm(vector<WClause*>& clauses) 
{
	while(1)
	{
		bool changed = false;
		for(unsigned int i=0;i<clauses.size();i++)
		{
			for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
			{
				for(unsigned int k=0;k<clauses.size();k++)
				{
					for(unsigned int m=0;m<clauses[k]->atoms.size();m++)
					{
						if(k==i && j==m)
							continue;
						if(clauses[i]->atoms[j]->symbol->id != clauses[k]->atoms[m]->symbol->id)
							continue;
						if(LvrHashAlgorithm::convertToHash(clauses[i]->atoms[j]) == LvrHashAlgorithm::convertToHash(clauses[k]->atoms[m]))
							continue;
						//need to normalize one or more terms
						vector<int> commonElems;
						int termId = -1;
						for(unsigned int t=0;t<clauses[k]->atoms[m]->terms.size();t++)
						{
							//check for the terms intersection, if no intersection does not need normalization
							vector<int> common(clauses[k]->atoms[m]->terms[t]->domain.size());
							vector<int>::iterator it = set_intersection(clauses[i]->atoms[j]->terms[t]->domain.begin(),
								clauses[i]->atoms[j]->terms[t]->domain.end(),
								clauses[k]->atoms[m]->terms[t]->domain.begin(),
								clauses[k]->atoms[m]->terms[t]->domain.end(),common.begin());
							int numElems = int(it - common.begin());
							if(numElems==0)
								continue;
							if(numElems == clauses[k]->atoms[m]->terms[t]->domain.size()
								&& numElems == clauses[i]->atoms[j]->terms[t]->domain.size())
								continue;
							common.erase(it,common.end());
							commonElems = common;
							termId = t;
							break;
						}
						if(termId==-1)
							continue;
						//compute the remainder in each atom
						vector<int> remainder1;
						vector<int> remainder2;
						if(commonElems.size()!=clauses[i]->atoms[j]->terms[termId]->domain.size())
						{
							vector<int> remaining1(clauses[i]->atoms[j]->terms[termId]->domain.size()-commonElems.size());
							set_difference(clauses[i]->atoms[j]->terms[termId]->domain.begin(),
								clauses[i]->atoms[j]->terms[termId]->domain.end()
								,commonElems.begin(),commonElems.end(),remaining1.begin());
							remainder1 = remaining1;
						}
						if(commonElems.size()!=clauses[k]->atoms[m]->terms[termId]->domain.size())
						{
							vector<int> remaining2(clauses[k]->atoms[m]->terms[termId]->domain.size()-commonElems.size());
							set_difference(clauses[k]->atoms[m]->terms[termId]->domain.begin(),
								clauses[k]->atoms[m]->terms[termId]->domain.end()
								,commonElems.begin(),commonElems.end(),remaining2.begin());
							remainder2 = remaining2;
						}
						//two atoms in distinct clauses need to be normalized at termId
						if(remainder1.size()>0)
						{
							changed=true;
							clauses[i]->atoms[j]->terms[termId]->domain.clear();
							clauses[i]->atoms[j]->terms[termId]->domain = remainder1;
						}
						if(remainder2.size()>0)
						{
							changed=true;
							clauses[k]->atoms[m]->terms[termId]->domain.clear();
							clauses[k]->atoms[m]->terms[termId]->domain = remainder2;
						}
						bool usedcommonvec=false;
						if(commonElems.size()>0 && remainder1.size() > 0)
						{
							changed=true;
							WClause* c1 = LvrMLN::create_new_clause(clauses[i]);
							c1->atoms[j]->terms[termId]->domain.clear();
							c1->atoms[j]->terms[termId]->domain = commonElems;
							usedcommonvec = true;
							clauses.push_back(c1);
						}
						if(commonElems.size()>0 && remainder2.size() > 0)
						{
							changed=true;
							WClause* c2 = LvrMLN::create_new_clause(clauses[k]);
							c2->atoms[m]->terms[termId]->domain.clear();
							if(!usedcommonvec)
								c2->atoms[m]->terms[termId]->domain = commonElems;
							else
							{
								//make a seperate copy
								c2->atoms[m]->terms[termId]->domain.resize(commonElems.size());
								for(unsigned int jj=0;jj<commonElems.size();jj++)
								{
									c2->atoms[m]->terms[termId]->domain[jj] = commonElems[jj];
								}
							}
							clauses.push_back(c2);
						}
					}
				}
			}
		}
		if(!changed)
			break;
	}
}

void LNormalizer::toNormalForm(PredicateSymbol* symbol,vector<WClause*>& clauses) {

	vector<WClause*> new_clauses;
	vector<WClause*> symbol_clauses;
	vector<int> symbol_loc;
	for (size_t i = 0; i < clauses.size(); i++) {
		bool relevant = false;
		for (size_t j = 0; j < clauses[i]->atoms.size(); j++) {
			if (clauses[i]->atoms[j]->symbol->id == symbol->id)
			{
				symbol_clauses.push_back(LvrMLN::create_new_clause(clauses[i]));
				symbol_loc.push_back(j);
				relevant = true;
				break;
			}
		}
		if (!relevant)
			new_clauses.push_back(LvrMLN::create_new_clause(clauses[i]));
	}
	for (size_t i = 0; i < symbol->variable_types.size(); i++) {
		bool changed = true;
		while (changed) {
			changed = false;
			for (size_t a = 0; a < symbol_clauses.size(); a++) {
				if (!symbol_clauses[a]->valid()) continue;
				for (size_t b = a + 1; b < symbol_clauses.size(); b++) {
					if (!symbol_clauses[b]->valid()) continue;
					int j = symbol_loc[a];
					int k = symbol_loc[b];
					vector<int>& vec1 =
							symbol_clauses[a]->atoms[j]->terms[i]->domain;
					vector<int>& vec2 =
							symbol_clauses[b]->atoms[k]->terms[i]->domain;
					// if the domains are disjoint

					if (vec1 == vec2 || is_disjoint(vec1, vec2)) {
						continue;
					} else {
						changed = true;
						WClause* clause1 = LvrMLN::create_new_clause(symbol_clauses[a]);
						WClause* clause2 = LvrMLN::create_new_clause(symbol_clauses[b]);
						do_set_difference(vec1, vec2,
								clause1->atoms[j]->terms[i]->domain);
						do_set_difference(vec2, vec1,
								clause2->atoms[k]->terms[i]->domain);
						do_set_intersection(vec1, vec2, vec1);
						vec2 = vec1;
						symbol_clauses.push_back(clause1);
						symbol_loc.push_back(j);
						symbol_clauses.push_back(clause2);
						symbol_loc.push_back(k);
						break;
					}
				}
				if (changed)
					break;
			}
		}
	}
	for (int i = 0; i < symbol_clauses.size(); i++){
		if(symbol_clauses[i]->valid())
			new_clauses.push_back(LvrMLN::create_new_clause(symbol_clauses[i]));
	}

	//cleanup old clauses that are no longer used
	//removeUnusedItems(clauses,new_clauses);
	cleanup(clauses);
	cleanup(symbol_clauses);
	clauses = new_clauses;
}
