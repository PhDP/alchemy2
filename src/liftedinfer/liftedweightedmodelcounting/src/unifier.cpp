#include "unifier.h"
#include <algorithm>
using namespace std;
#include "cleanuputils.h"

bool compareFunc(Atom* a1, Atom* a2)
{
	return a1->symbol->symbol.compare(a2->symbol->symbol)<0;
}


void LUnifier::getTermRelations(WClause* clause,vector<LvrTerm*>& clause_terms, vector<vector<int> >& clause_positions)
{
	
	for(unsigned int i=0;i<clause->atoms.size();i++){
		if(clause->atoms[i]->isConstant())
			continue;
		for(unsigned int j=0;j<clause->atoms[i]->terms.size();j++){
			bool found = false;
			for(unsigned int k=0;k<clause_terms.size();k++)
			{
				if(clause_terms[k]==clause->atoms[i]->terms[j])
				{
					found=true;
					break;
				}
				
			}
			if(!found)
				clause_terms.push_back(clause->atoms[i]->terms[j]);
		}
	}
	clause_positions.resize(clause_terms.size());
	for(unsigned int i=0;i<clause->atoms.size();i++){
		for(unsigned int j=0;j<clause->atoms[i]->terms.size();j++){
			int matchingIndex = -1;
			for(unsigned int k=0;k<clause_terms.size();k++)
			{
				if(clause_terms[k]==clause->atoms[i]->terms[j])
				{
					matchingIndex = k;
					break;
				}
			}
			if(matchingIndex!=-1)
				clause_positions[matchingIndex].push_back(j);
		}
	}
}

bool LUnifier::clauseUnify(WClause* c1, WClause* c2,map<int,vector<int> >& originalIdConsts)
{ 
	if(c1->atoms.size() != c2->atoms.size() || c1->satisfied != c2->satisfied)
		return false;
	
	//c1->print();
	sort(c1->atoms.begin(),c1->atoms.end(),compareFunc);
	//c1->print();
	sort(c2->atoms.begin(),c2->atoms.end(),compareFunc);
	//c2->print();
	for(unsigned int i=0;i<c1->atoms.size();i++)
	{
		if(c1->sign[i]!=c2->sign[i] || c1->atoms[i]->symbol->symbol.compare(c2->atoms[i]->symbol->symbol)!=0)
			return false;
		//check if term sizes are the same
		if(c1->atoms[i]->terms.size()!=c2->atoms[i]->terms.size())
			return false;
		//check if constants are the same
		
		for(unsigned int j=0;j<c1->atoms[i]->terms.size();j++)
		{
			if(c1->atoms[i]->terms[j]->domain.size() != c2->atoms[i]->terms[j]->domain.size())
				return false;
			else
			{
				if(c1->atoms[i]->terms[j]->domain.size() == 1)
				{
					if(c1->atoms[i]->terms[j]->domain[0]!=c2->atoms[i]->terms[j]->domain[0])
					{
						//if already mapped return false else map to new value
						if(originalIdConsts.find(c1->atoms[i]->symbol->parentId)!=originalIdConsts.end())
						{
							if(originalIdConsts[c1->atoms[i]->symbol->parentId].at(j)==-1)
							{
								originalIdConsts[c1->atoms[i]->symbol->parentId].at(j) = c2->atoms[i]->terms[j]->domain[0];
							}
							else if(originalIdConsts[c1->atoms[i]->symbol->parentId].at(j)!=c2->atoms[i]->terms[j]->domain[0])
								return false;
						}
					}
				}
			}
		}
		
	}
	vector<LvrTerm*> terms1;
	vector<LvrTerm*> terms2;
	vector<vector<int> > pos1;
	vector<vector<int> > pos2;
	getTermRelations(c1,terms1,pos1);
	getTermRelations(c2,terms2,pos2);
	if(terms1.size()!=terms2.size())
		return false;
	int size = terms1.size();
	for(unsigned int i=0;i<size;i++)
	{
		if(terms1[i]->domain.size()!=terms2[i]->domain.size())
			return false;
		else
		{
			
			//else if the domain sizes are identical we need to check equvalency of content
			//if(!includes(terms1[i]->domain.begin(),terms1[i]->domain.end(),terms2[i]->domain.begin(),terms2[i]->domain.end()))
			//	return false;
			if(pos1[i].size()!=pos2[i].size())
				return false;
			else
			{
				for(unsigned int j=0;j<pos1.size();j++)
				{
					if(pos1[j]!=pos2[j])
						return false;
				}
			}
		}
	}
	return true;
}

bool LUnifier::CNFUnify(vector<WClause*> CNF1, vector<WClause*> CNF2)
{
	vector<bool> completedIndex(CNF2.size());
	//vector<vector<int> > originalIdConsts(mln.symbols.size());
	map<int,vector<int> > originalIdConsts;
	int maxDegree = mln.getMaxDegree();
	for(unsigned int i=0;i<mln.symbols.size();i++)
	{
		vector<int> temp(maxDegree);
		for(unsigned int j=0;j<temp.size();j++)
			temp[j] = -1;
		originalIdConsts[mln.symbols[i]->parentId] = temp;
	}

	for(unsigned int i=0;i<CNF1.size();i++)
	{
		bool unified = false;
		for(unsigned int j=0;j<CNF2.size();j++)
		{
			if(completedIndex[j])
				continue;
			if(clauseUnify(CNF1[i],CNF2[j],originalIdConsts))
			{
				completedIndex[j]=true;
				unified = true;
				break;
			}
		}
		if(!unified)
			return false;
	}
	
	return true;
}
