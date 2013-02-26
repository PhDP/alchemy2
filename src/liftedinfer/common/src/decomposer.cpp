#include "decomposer.h"
#include "cleanuputils.h"

Decomposer::~Decomposer()
{
	decomposer_terms.clear();
	predicate_positions.clear();
	atomCounter.clear();
}

//creates a new decomposer item
Decomposer* LDecomposer::create_new_decomposer(LvrTerm* term,vector<int> atom_list,vector<int> position_list,vector<int> norm_atom_list)
{
	vector<LvrTerm*> decomposer_terms;
	map<int,int> predicate_positions;
	map<int,int> atomCounter;
	map<int,int> norm_predicate_positions;
	decomposer_terms.push_back(term);
	for(unsigned int k=0;k<atom_list.size();k++)
	{
		predicate_positions[atom_list[k]]=position_list[k];
		norm_predicate_positions[norm_atom_list[k]] = position_list[k];
		atomCounter[atom_list[k]]=1;
	}
	Decomposer* d = new Decomposer(decomposer_terms,predicate_positions,atomCounter,norm_predicate_positions);
	return d;
}
//0 - did not unify, 1- unified, -1 - disjoint
int LDecomposer::unify(Decomposer* d1, Decomposer* d2)
{
	bool found=false;
	for(map<int,int>::iterator it1=d1->predicate_positions.begin();it1!=d1->predicate_positions.end();it1++)
	{
		map<int,int>::iterator tmp = d2->predicate_positions.find(it1->first);
		if(tmp!=d2->predicate_positions.end())
		{
			found=true;
			//found d1's predicate in d2, match the positions
			if(tmp->second != it1->second)
				return 0;
		}
	}
	if(found)
		return 1;
	else
		return -1;
}

Decomposer* LDecomposer::merge(Decomposer* d1, Decomposer* d2)
{
	vector<LvrTerm*> decomposerTerms;
	map<int,int> predicate_positions;
	map<int,int> atomCounter;
	map<int,int> norm_predicate_positions;
	for(unsigned int i=0;i<d1->decomposer_terms.size();i++)
	{
		decomposerTerms.push_back(d1->decomposer_terms[i]);
	}
	for(unsigned int i=0;i<d2->decomposer_terms.size();i++)
	{
		decomposerTerms.push_back(d2->decomposer_terms[i]);
	}
	for(map<int,int>::iterator it1=d1->predicate_positions.begin();it1!=d1->predicate_positions.end();it1++)
	{
		predicate_positions.insert(pair<int,int>(it1->first,it1->second));
	}
	for(map<int,int>::iterator it1=d2->predicate_positions.begin();it1!=d2->predicate_positions.end();it1++)
	{
		predicate_positions.insert(pair<int,int>(it1->first,it1->second));
	}
	for(map<int,int>::iterator it1=d1->atomCounter.begin();it1!=d1->atomCounter.end();it1++)
	{
		atomCounter.insert(pair<int,int>(it1->first,it1->second));
	}
	for(map<int,int>::iterator it1=d2->atomCounter.begin();it1!=d2->atomCounter.end();it1++)
	{
		map<int,int>::iterator tmp = atomCounter.find(it1->first);
		if(tmp!=atomCounter.end())
		{
			tmp->second = it1->second+tmp->second;
		}
		else
			atomCounter.insert(pair<int,int>(it1->first,it1->second));
	}
	for(map<int,int>::iterator it1=d1->norm_predicate_positions.begin();it1!=d1->norm_predicate_positions.end();it1++)
	{
		norm_predicate_positions.insert(pair<int,int>(it1->first,it1->second));
	}
	for(map<int,int>::iterator it1=d2->norm_predicate_positions.begin();it1!=d2->norm_predicate_positions.end();it1++)
	{
		norm_predicate_positions.insert(pair<int,int>(it1->first,it1->second));
	}

	d1->deletionMarker = true;
	d2->deletionMarker = true;
	Decomposer* newDecomp = new Decomposer(decomposerTerms,predicate_positions,atomCounter,norm_predicate_positions);
	return newDecomp;
}

//appends new decomposers to the decomposer list
bool LDecomposer::append_to_decomposer (vector<int> atom_list,vector<LvrTerm*> terms,vector<vector<int> >positions,
	vector<Decomposer*>& decomposer_list,vector<int> norm_atom_list)
{
	bool unified = true;
	vector<Decomposer*> potentialDecomposers;
	//Make new decomposer and Append to potential decomposers list
	for(unsigned int t1=0;t1<positions.size();t1++)
	{
		vector<int> position_list = positions[t1];
		LvrTerm* term=terms[t1];
		potentialDecomposers.push_back(create_new_decomposer(term,atom_list,position_list,norm_atom_list));
	}
	if(decomposer_list.size()==0)//decomposer is empty
	{
		//return the potential decomposer list
		decomposer_list = potentialDecomposers;
		return unified;
	}

	//match every potential decomposer with existing list of decomposers
	vector<Decomposer*> decomposer_add;
	int nonUnifiedTermCount = 0;
	for(unsigned int i=0;i<potentialDecomposers.size();i++)
	{
		int notFoundCount=0;
		vector<int> unifiedIndex;
		int numNonUnify = 0;
		for(unsigned int j=0;j<decomposer_list.size();j++)
		{
			//Unify new potential decomposer with existing decomposers
			int status = unify(potentialDecomposers[i],decomposer_list[j]);
			if(status==0)
			{
				numNonUnify++;
				decomposer_list[j]->deletionMarker = true;
			}
			else if(status == 1)
			{
				//merge the decomposers into 1
				unifiedIndex.push_back(j);
			}
			else
				notFoundCount++;
		}
		if(numNonUnify == decomposer_list.size())
		{
			//did not unify with anything
			nonUnifiedTermCount++;
		}
		if(notFoundCount == decomposer_list.size())
		{
			//disjoint with all existing sets,append
			decomposer_add.push_back(potentialDecomposers[i]);
		}
		else if(unifiedIndex.size() > 0)
		{
			//merge all unified rows into 1
			Decomposer* tempD = potentialDecomposers[i];
			for(unsigned int k=0;k<unifiedIndex.size();k++)
			{
				Decomposer* newD = merge(tempD,decomposer_list[unifiedIndex[k]]);
				tempD = newD;
			}

			decomposer_add.push_back(tempD);
		}
	}
	
	if(nonUnifiedTermCount == potentialDecomposers.size())
	{
		//no terms unified with any decomposer row
		unified = false;
	}
	//add the new decomposer rows
	for(unsigned int i=0;i<decomposer_add.size();i++)
	{
		decomposer_list.push_back(decomposer_add[i]);
	}

	for(unsigned int i=0;i<decomposer_list.size();i++)
	{
		if(decomposer_list.size() == 0)
			break;
		if(decomposer_list[i]->deletionMarker)
		{
			removeItem(decomposer_list,i);
			i--;
		}
	}
	return unified;
}


void LDecomposer::removeRowsFromDecomposer(WClause* clause, vector<Decomposer*>& decomposer_list)
{
	//remove from the potential decomposers list all predicates intersecting with clause i's predicates
	vector<int> rowsToDelete;
	for(unsigned int j=0;j<clause->atoms.size();j++)
	{
		for(unsigned int k=0;k<decomposer_list.size();k++)
		{
			if(decomposer_list[k]->predicate_positions.count(clause->atoms[j]->symbol->id)> 0)
			{
				//remove entire decomposer row
				rowsToDelete.push_back(k);
			}
		}
	}
	int removedRows=0;
	for(unsigned int k=0;k<rowsToDelete.size();k++)
	{
		if(decomposer_list.empty())
			break;
		int updatedRow=rowsToDelete[k]-removedRows-1;
		if(updatedRow<0)
			updatedRow=0;
		if(updatedRow>decomposer_list.size()-1)
			updatedRow=decomposer_list.size()-1;
		removeItem(decomposer_list,updatedRow);
		removedRows++;
	}

}

void LDecomposer::find_decomposer(vector<WClause*>& CNF,vector<Decomposer*>& decomposer_list)
{
	//vector<bool> nonUnifiedList(mln.getMaxPredicateId());
	set<int> nonUnifiedList;
	//map<int,int> predicateCounter;
	map<int,int>::iterator it;
	//vector<int> predicateCounter(mln.getMaxPredicateId());
	map<int,int> predicateCounter;
	for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->atoms.size()==0)
			continue;
		if(CNF[i]->isPropositional())
			continue;
		vector<LvrTerm*> terms;
		vector<vector<int> > positions;
		vector<int> atom_list;
		vector<int> norm_atom_list;
		//vector<bool> completedAtoms(mln.getMaxPredicateId());
		set<int> completedAtoms;
		for(unsigned int k=0;k<CNF[i]->atoms.size();k++)
		{
			int id = CNF[i]->atoms[k]->symbol->id;
			atom_list.push_back(id);
			norm_atom_list.push_back(CNF[i]->atoms[k]->symbol->normParentId);
			if(completedAtoms.count(id)==0)
			{
				completedAtoms.insert(id);
				predicateCounter[id]++;
			}

		}
		//check if any atom occurs in failed list
		bool nonunified = false;
		for(unsigned int k=0;k<atom_list.size();k++)
		{
			if(nonUnifiedList.count(atom_list[k]) != 0)
			{
				nonunified = true;
				break;
			}
		}
		if(!nonunified)
			CNF[i]->find_decomposer(terms,positions);
		if(terms.size()==0)
		{
			removeRowsFromDecomposer(CNF[i],decomposer_list);
			//add all atoms to non unified list
			for(unsigned int k=0;k<atom_list.size();k++)
				nonUnifiedList.insert(atom_list[k]);
			continue;
		}
		bool unified = append_to_decomposer(atom_list,terms,positions,decomposer_list,norm_atom_list);
		if(!unified)
		{
			removeRowsFromDecomposer(CNF[i],decomposer_list);
			//add all atoms to non unified list
			for(unsigned int k=0;k<atom_list.size();k++)
				nonUnifiedList.insert(atom_list[k]);
		}
#ifdef __DEBUG_PRINT2__
		for(unsigned int x=0;x<decomposer_list.size();x++)
		{
			cout<<"[";
			for(unsigned int y=0;y<decomposer_list[x]->decomposer_terms.size();y++)
				cout<<decomposer_list[x]->decomposer_terms[y]<<" ";
			cout<<"]"<<endl;
			cout<<"{";
			for(map<int,int>::iterator it1=decomposer_list[x]->predicate_positions.begin();it1!=decomposer_list[x]->predicate_positions.end();it1++)
			{
				cout<<it1->first<<"::"<<it1->second<<",";
			}
			cout<<"}"<<endl;
			cout<<"(";
			for(map<int,int>::iterator it1=decomposer_list[x]->atomCounter.begin();it1!=decomposer_list[x]->atomCounter.end();it1++)
			{
				cout<<it1->first<<"::"<<it1->second<<",";
			}
			cout<<")"<<endl;
		}
		cout<<"--------------------------------"<<endl;
		
#endif		
	}
	
		for(int x=0;x<decomposer_list.size();x++)
		{
			if(decomposer_list.size()==0)
				break;
			bool isDecomposer = true;
			for(map<int,int>::iterator it1=decomposer_list[x]->atomCounter.begin();it1!=decomposer_list[x]->atomCounter.end();it1++)
			{
				if(predicateCounter[it1->first]!=it1->second)
				{
					isDecomposer = false;
					break;
				}
			}
			if(!isDecomposer)
			{
				removeItem(decomposer_list,x);
				x--;
			}
		}
}

bool LDecomposer::decomposeCNF(vector<WClause*>& CNF,int& powerFactor)
{
	vector<Decomposer*> decomposer_list;
	find_decomposer(CNF,decomposer_list);
	
	if(decomposer_list.size()==0)
		return false;
	else
	{
		//replace all decomposer terms with a constant from the domain
		powerFactor = 1;
		for(unsigned int i=0;i<decomposer_list.size();i++)
		{
			powerFactor*=decomposer_list[i]->decomposer_terms[0]->domain.size();
			for(unsigned int j=0;j<decomposer_list[i]->decomposer_terms.size();j++)
			{
				//Erase all but one grounding
				decomposer_list[i]->decomposer_terms[j]->domain.erase(decomposer_list[i]->decomposer_terms[j]->domain.begin(),decomposer_list[i]->decomposer_terms[j]->domain.begin()+decomposer_list[i]->decomposer_terms[j]->domain.size()-1);

			}
		}
	}
	cleanup(decomposer_list);
	return true;
}


