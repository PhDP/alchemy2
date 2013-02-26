#ifndef LDECOMPOSER_H_
#define LDECOMPOSER_H_
#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;

#include "lvrmln.h"

struct Decomposer
{
	vector<LvrTerm*> decomposer_terms;
	//key is predicatesymbol id, value is position the term occupies
	map<int,int> predicate_positions;
	map<int,int> norm_predicate_positions;
	//key is predicatesymbol id, value is number of predicates sharing a common position
	map<int,int> atomCounter;
	bool deletionMarker;
	Decomposer(){}
	Decomposer(vector<LvrTerm*> decomposer_terms_, map<int,int> predicate_positions_,map<int,int> atomCounter_,
		map<int,int> norm_predicate_positions_):decomposer_terms(decomposer_terms_),
		predicate_positions(predicate_positions_),atomCounter(atomCounter_),deletionMarker(false),norm_predicate_positions(norm_predicate_positions_)
	{}
	void print()
	{
		cout<<"[";
		for(unsigned int j=0;j<decomposer_terms.size();j++)
		{
			cout<<decomposer_terms[j];
			if(j!=decomposer_terms.size()-1)
				cout<<" , ";
		}
		cout<<"]";
	}

	~Decomposer();
};

struct LDecomposer
{
	LvrMLN& mln;
	bool decomposeCNF(vector<WClause*>& CNF,int& powerFactor);
	void find_decomposer(vector<WClause*>& CNF, vector<Decomposer*>& decomposer_list);
	//void find_decomposer_new(vector<WClause*>& CNF,vector<Decomposer*>& decomposer_list, int curr_predicate_count);
	LDecomposer(LvrMLN& mln_):mln(mln_){}
private:
	bool append_to_decomposer (vector<int> atom_list,vector<LvrTerm*> terms,vector<vector<int> >positions,vector<Decomposer*>& decomposer_list,
		vector<int> norm_atom_list);
	Decomposer* create_new_decomposer(LvrTerm* term,vector<int> atom_list,vector<int> position_list,vector<int> norm_atom_list);
	void removeRowsFromDecomposer(WClause* clause, vector<Decomposer*>& decomposer_list);
	Decomposer* merge(Decomposer* d1, Decomposer* d2);
	int unify(Decomposer* d1, Decomposer* d2);
};

#endif