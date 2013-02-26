#include "heuristics.h"
#include "cleanuputils.h"
#include "queryupdater.h"

Atom* LHeuristics::getAtomToSplit(vector<WClause*> CNF,bool ptpexactmar)
{
	//get the atom which is singleton and has minimum domain
	Atom* minDomainSingleton = NULL;
	int minNumDomainsST = 100000;
	Atom* minDomainNonSingleton = NULL;
	int minNumDomainsNST = 100000;
	Atom* minDomainUCAtom = NULL;
	int minNumDomainsUCAtom = 100000;
	int maxDegreeSingleton = 0;
	int maxDegreeNST = 0;
	int maxDegreeUC = 0;
	int maxDegreeConstant = 0;
	Atom* atomMaxDegree_S = NULL;
	Atom* atomMaxDegree_NS = NULL;
	Atom* atomMaxDegree_UC = NULL;
	Atom* atomMaxDegree_Const = NULL;
	//compute the number of unsatisfied clauses that each predicate participates in
	map<int,int> degree;
	for(unsigned int i=0;i<CNF.size();i++)
	{
		bool checkforsatclause=false;
		if(ptpexactmar)
		{
			if(CNF[i]->satisfied)
			{
				if(LvrQueryUpdater::isInstanceCreated())
				{
					for(unsigned int jj=0;jj<CNF[i]->atoms.size();jj++)
					{
						if(LvrQueryUpdater::Instance()->isNormIdInQuery(CNF[i]->atoms[jj]->symbol->normParentId))
						{
							checkforsatclause=true;
							break;
						}
					}
				}
			}
		}
		
		if(!CNF[i]->satisfied || checkforsatclause)
		//if(!CNF[i]->satisfied)
		{
			for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
			{
				map<int,int>::iterator it = degree.find(CNF[i]->atoms[j]->symbol->id);
				if(it!=degree.end())
				{
					//increment counter
					it->second++;
				}
				else
				{
					//insert new element to map
					degree.insert(pair<int,int>(CNF[i]->atoms[j]->symbol->id,1));
				}
			}
		}
	}
	for(unsigned int i=0;i<CNF.size();i++)
	{
		
		bool checkforsatclause=false;
		if(ptpexactmar)
		{
			if(CNF[i]->satisfied)
			{
				if(LvrQueryUpdater::isInstanceCreated())
				{
					for(unsigned int jj=0;jj<CNF[i]->atoms.size();jj++)
					{
						if(LvrQueryUpdater::Instance()->isNormIdInQuery(CNF[i]->atoms[jj]->symbol->normParentId))
						{
							checkforsatclause=true;
							break;
						}
					}
				}
			}
		}
		
		//do not choose atoms in satisfied clauses
		if(CNF[i]->satisfied && !checkforsatclause)
		//if(CNF[i]->satisfied)
			continue;
		for(unsigned int k=0;k<CNF[i]->atoms.size();k++)
		{
			
			int index;
			int atomDeg = 0;
			map<int,int>::iterator it = degree.find(CNF[i]->atoms[k]->symbol->id);
			if(it!=degree.end())
			{
				atomDeg = it->second;
			}
			if(CNF[i]->atoms[k]->isSingletonAtom(index))
			{
				//if singleton, update the min domains and max degree
				if(CNF[i]->atoms[k]->terms[index]->domain.size() < minNumDomainsST)
				{
					minDomainSingleton=CNF[i]->atoms[k];
					minNumDomainsST=CNF[i]->atoms[k]->terms[index]->domain.size();
				}
				if(atomDeg > maxDegreeSingleton)
				{
					maxDegreeSingleton = atomDeg;
					atomMaxDegree_S = CNF[i]->atoms[k];
				}
			}
			else if(CNF[i]->atoms[k]->terms.size()>1)
			{
				int totalSize=1;
				for(unsigned int m=0;m<CNF[i]->atoms[k]->terms.size();m++)
				{
					totalSize*=CNF[i]->atoms[k]->terms[m]->domain.size();
				}
				if(totalSize < minNumDomainsNST )
				{
					minDomainNonSingleton=CNF[i]->atoms[k];
					minNumDomainsNST=totalSize;
				}
				if(atomDeg > maxDegreeNST)
				{
					maxDegreeNST= atomDeg;
					atomMaxDegree_NS = CNF[i]->atoms[k];
				}

				if(CNF[i]->atoms.size()==1)
				{
					//unit clause
					if(totalSize < minNumDomainsUCAtom)
					{
						minNumDomainsUCAtom = totalSize;
						minDomainUCAtom = CNF[i]->atoms[k];

					}
					if(atomDeg > maxDegreeUC)
					{
						maxDegreeUC= atomDeg;
						atomMaxDegree_UC = CNF[i]->atoms[k];
					}

				}
			}
			else if(CNF[i]->atoms[k]->isConstant())
			{
				if(atomDeg > maxDegreeConstant)
				{
					maxDegreeConstant= atomDeg;
					atomMaxDegree_Const = CNF[i]->atoms[k];
				}
			}
		
		}
	}
	if(atomMaxDegree_Const)
		return atomMaxDegree_Const;
	if(minDomainSingleton)
		return minDomainSingleton;
	else
	{
		//find best non singleton using lookahead
		/*Atom* atom = decomposition_lookahead(CNF);
		if(atom!=NULL)
			return atom;
		else
		*/
			return minDomainNonSingleton;
	}
}


Atom* LHeuristics::decomposition_lookahead(vector<WClause*> clauses)
{
	int maxSize = 0;
	Atom* atom = NULL;
	map<int,Atom*> predicates;
	//collect atoms
	for(unsigned int i=0;i<clauses.size();i++)
	{
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			if(!clauses[i]->atoms[j]->isConstant() && clauses[i]->atoms[j]->symbol->symbol.find(NEWPREDICATEPREFIX)==string::npos)
				predicates.insert(pair<int,Atom*>(clauses[i]->atoms[j]->symbol->id,clauses[i]->atoms[j]));
		}
	}
	//for each atom split, try to find a decomposer
	for(map<int,Atom*>::iterator it=predicates.begin();it!=predicates.end();it++)
	{
		//work with a CNF copy
		vector<WClause*> CNF;
		LvrMLN::copyAllClauses(clauses,CNF);
		for(unsigned int i=0;i<CNF.size();i++)
		{
			for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
			{
				if(CNF[i]->atoms[j]->symbol->id==it->first)
				{
					//found atom
					for(unsigned int k=0;k<CNF[i]->atoms[j]->terms.size();k++)
					{
						//make that term any constant
						CNF[i]->atoms[j]->terms[k]->domain.clear();
						CNF[i]->atoms[j]->terms[k]->domain.push_back(-1);
					}
					//remove atom
					CNF[i]->removeAtom(j);
					j--;
				}
			}
		}
		//Find decomposer in CNF
		vector<vector<WClause*> > decomposedList;
		vector<Decomposer*> decomp;
		decomposer.find_decomposer(CNF,decomp);
		if(decomp.size() > maxSize)
		{
			maxSize = decomp.size();
			atom = it->second;
		}
		cleanup(CNF);
	}
	return atom;
}


Atom* LHeuristics::getAtomToSplitH2(vector<WClause*> CNF)
{
	Atom* minDomSingleton = NULL;
	int minDomSingletonSize = RAND_MAX;
	Atom* constant = NULL;
	Atom* minDomNonSingleton = NULL;
	int minDomNonSingletonSize = RAND_MAX;

	for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->satisfied)
			continue;
		for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
		{
			int singletonind;
			if(CNF[i]->atoms[j]->isSingletonAtom(singletonind))
			{
				if(CNF[i]->atoms[j]->getNumberOfGroundings() < minDomSingletonSize)
					minDomSingleton = CNF[i]->atoms[j]; 
			}
			else if(CNF[i]->atoms[j]->isConstant())
			{
				constant = CNF[i]->atoms[j];
			}
			else
			{
				if(CNF[i]->atoms[j]->getNumberOfGroundings() < minDomNonSingletonSize)
				{
					minDomNonSingleton = CNF[i]->atoms[j];
				}
			}
		}
	}
	if(minDomSingleton)
		return minDomSingleton;
	else if(constant)
		return constant;
	else
		return minDomNonSingleton;
}

int LHeuristics::proposalCouplingScore(vector<WClause*> origClauses,Atom* a1, Atom* a2)
{
	int totalScore = 0;
	for(unsigned int i=0;i<origClauses.size();i++)
	{
		//check if the two atoms occur together in clauses
		int a1Ind=-1;
		int a2Ind=-1;
		for(unsigned int j=0;j<origClauses[i]->atoms.size();j++)
		{
			if(origClauses[i]->atoms[j]->symbol->parentId == a1->symbol->parentId)
				a1Ind = j;
			else if(origClauses[i]->atoms[j]->symbol->parentId == a2->symbol->parentId)
				a2Ind = j;
		}
		if(a1Ind==-1 || a2Ind ==-1)
			continue;
		//check if the groundings contain a1, a2
		WClause* tmpClause = LvrMLN::create_new_clause(origClauses[i]);
		//ground to make it equivalent to a1
		for(unsigned int j=0;j<a1->terms.size();j++)
		{
			if(a1->terms[j]->domain.size() != tmpClause->atoms[a1Ind]->terms[j]->domain.size())
			{
				tmpClause->atoms[a1Ind]->terms[j]->domain.clear();
				for(unsigned int k=0;k<a1->terms[j]->domain.size();k++)
					tmpClause->atoms[a1Ind]->terms[j]->domain.push_back(a1->terms[j]->domain[k]);
			}
		}
		
		bool mismatch = false;
		//ground to make it equivalent to a2
		for(unsigned int j=0;j<a2->terms.size();j++)
		{
			if(tmpClause->atoms[a2Ind]->terms[j]->domain.size() == 1)
			{
				if(a2->terms[j]->domain[0]!=tmpClause->atoms[a2Ind]->terms[j]->domain[0])
				{
					mismatch = true;
					break;
				}
				else
				{
					continue;
				}
			}
			if(!includes(tmpClause->atoms[a2Ind]->terms[j]->domain.begin(),tmpClause->atoms[a2Ind]->terms[j]->domain.end(),
				a2->terms[j]->domain.begin(),a2->terms[j]->domain.end()))
			{
				mismatch=true;
				break;
			}
			tmpClause->atoms[a2Ind]->terms[j]->domain.clear();
			for(unsigned int k=0;k<a2->terms[j]->domain.size();k++)
				tmpClause->atoms[a2Ind]->terms[j]->domain.push_back(a2->terms[j]->domain[k]);

		}
		if(mismatch)
			continue;
		totalScore += a1->getNumberOfGroundings()*a2->getNumberOfGroundings();
		delete tmpClause;
	}
	return totalScore;
}


