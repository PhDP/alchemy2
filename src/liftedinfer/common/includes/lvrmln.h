#ifndef LVRMLN_H_
#define LVRMLN_H_

#include <vector>
#include <iostream>
#include <set>
#include <map>
#include <algorithm>
#include "logdouble.h"
using namespace std;
#include "paramconstants.h"

struct LvrTerm
{
	// There is a one to one mapping between constants and the domain
	int type;
	vector<int> domain;
	//pre-deomposition domain
	int origDomainSize;
	vector<int> origDomain;

	int currentSampledPos;
	LvrTerm():origDomainSize(-1),currentSampledPos(0){}
	LvrTerm(int type_, vector<int>& domain_): type(type_),domain(domain_),origDomainSize(-1),currentSampledPos(0){}
	LvrTerm(int type_,int value):type(type),domain(vector<int>(1)),origDomainSize(-1),currentSampledPos(0){ domain[0]=value;}
	LvrTerm(const LvrTerm& term)
	{
		type = term.type;
		for(unsigned int k=0;k<term.domain.size();k++)
			domain.push_back(term.domain[k]);
		origDomainSize = term.origDomainSize;
	}
	~LvrTerm()
	{
		domain.clear();
	}
};
struct PredicateSymbol
{
	int id;
	//USED IN CACHING
	int parentId;
	string symbol;
	vector<int> variable_types;
	LogDouble pweight;
	LogDouble nweight;
	bool isOriginalSymbol;
	//USED IN NORMALIZATION
	int normParentId;
	PredicateSymbol(){}
	PredicateSymbol(int id_, string symbol_, vector<int> var_types,LogDouble pweight_,LogDouble nweight_,
		int normParentId_=0, bool isOriginalSymbol_ = true): id(id_),symbol(symbol_),variable_types(var_types),
		pweight(pweight_),nweight(nweight_),isOriginalSymbol(isOriginalSymbol_),normParentId(normParentId_)
	{
		parentId = id;
	}
	~PredicateSymbol()
	{
		variable_types.clear();
	}
};
struct Atom
{
	PredicateSymbol* symbol;
	//terms may be shared across atoms
	vector<LvrTerm*> terms;

	// Conditions: terms.size()=symbol->variables.size() and terms[i].type=symbol->variables[i]
	Atom(){}
	Atom(PredicateSymbol* symbol_,vector<LvrTerm*> terms_): symbol(symbol_),terms(terms_){}
	Atom(const Atom& atom)
	{
		symbol = atom.symbol;
		for(unsigned int i=0;i<atom.terms.size();i++)
		{
			LvrTerm* tm(atom.terms[i]);
			terms.push_back(tm);
		}
	}
	bool isConstant();
	bool isSingletonAtom(int& variableIndex);
	void print(bool verbose=true);
	int getNumberOfGroundings();
	~Atom();
};

struct ClauseDecomposer
{
	vector<LvrTerm*> terms;
	vector<vector<int> > positions;
};

struct WClause
{
	vector<Atom*> atoms;
	vector<bool> sign;
	LogDouble weight;
	bool satisfied;
	int originalClauseIndex;
	WClause():satisfied(false){}
	bool valid(){
		for(unsigned int i=0;i<atoms.size();i++){
			for(unsigned int j=0;j<atoms[i]->terms.size();j++){
							if(atoms[i]->terms[j]->domain.empty()) return false;
			}
		}
		
		return true;
	}
	void print(){
		if(satisfied)
			cout<<"satisfied v ";
		if(atoms.size()==0)
			cout<<"{ }";
		for(unsigned int i=0;i<atoms.size();i++){
			if(sign[i])
				cout<<"!";
			cout<<atoms[i]->symbol->symbol.c_str()<<"["<< atoms[i]->symbol->id<<","<<atoms[i]->symbol->parentId<<"{"<<atoms[i]->symbol->pweight<<","<<atoms[i]->symbol->nweight<<"}]"<<" (";
			for(unsigned int j=0;j<atoms[i]->terms.size();j++){
				if((int)atoms[i]->terms[j]->domain.size()==1){
					cout<<atoms[i]->terms[j]->domain[0];
				}
				else{
					cout<<atoms[i]->terms[j]<<"[ #"<<atoms[i]->terms[j]->domain.size()<<" ]";
					if(atoms[i]->terms[j]->origDomainSize!=-1)
					{
						cout<<"D#"<<atoms[i]->terms[j]->origDomainSize;
					}
				}
				if(j!=atoms[i]->terms.size()-1)
					cout<<",";
			}
			cout<<")";
			if(i!=atoms.size()-1)
				cout<<" V ";

		}
		cout<<"::"<<weight<<endl;
	}
	void find_decomposer(vector<LvrTerm*>& terms, vector<vector<int> >& positions);
	~WClause();
	void removeAtom(int index);
	void removeAtomWithWeightUpdation(int index);
	void findSelfJoins(map<int,vector<int> >& selfJoinedAtoms);
	bool isSelfJoinedOnAtom(Atom* atom);
	bool isSelfJoined()
	{
		for(unsigned int i=0;i<atoms.size();i++)
		{
			if(isSelfJoinedOnAtom(atoms[i]))
			{
				return true;
			}
		}
		return false;
	}
	bool isPropositional();
	int getClauseSize()
	{
		set<LvrTerm*> terms;
		for(unsigned int i=0;i<atoms.size();i++)
		{
			for(unsigned int j=0;j<atoms[i]->terms.size();j++)
				terms.insert(atoms[i]->terms[j]);
		}
		int sz =1;
		for(set<LvrTerm*>::iterator it=terms.begin();it!=terms.end();it++)
			sz *= (*it)->domain.size();
		return sz;
	}
	bool isAtomInClause(Atom* atom)
	{
		for(unsigned int i=0;i<atoms.size();i++)
		{
			if(atom->symbol->id == atoms[i]->symbol->id)
				return true;
		}
		return false;
	}
	int getNumberOfGroundedClauses();
};
struct Formula
{
	int MLNClauseStartIndex;
	int MLNClauseEndIndex;
	bool isEvidence;
	LogDouble weight;
	Formula(int MLNClauseStartIndex_,int MLNClauseEndIndex_,LogDouble weight_, bool isEvidence_ = false):MLNClauseStartIndex(MLNClauseStartIndex_),
		MLNClauseEndIndex(MLNClauseEndIndex_),weight(weight_),isEvidence(isEvidence_){}
};

struct LvrMLN
{
	
	vector<PredicateSymbol*> symbols;
	vector<WClause*> clauses;
	vector<Formula*> formulas;
	static WClause* create_new_clause(WClause* clause);
	static Atom* create_new_atom(Atom* atom);
	static LvrTerm* create_new_term(LvrTerm* term);
	static PredicateSymbol* create_new_symbol(PredicateSymbol* symbol);
	static void copyAllClauses(vector<WClause*> origClauses, vector<WClause*>& newClauses);
	static void print(vector<WClause*> clauses,string banner="CNF")
	{
		string tmp("###########");
		tmp.append(banner);
		tmp.append("-START");
		tmp.append("###########");
		cout<<tmp.c_str()<<endl;
		for(unsigned int i=0;i<clauses.size();i++)
			clauses[i]->print();
		string tmp1("###########");
		tmp1.append(banner);
		tmp1.append("-END");
		tmp1.append("###########");
		cout<<tmp1.c_str()<<endl;
	}
	void setMaxPredicateId(int Id)
	{
		max_predicate_id = max(max_predicate_id,Id);
	}
	int getMaxPredicateId()
	{
		/*if(max_predicate_id == 0)
			setMaxPredicateId(symbols.size());
		return max_predicate_id;*/
		max_predicate_id = -1;
		for(unsigned int i=0;i<symbols.size();i++)
		{
			if(symbols[i]->id > max_predicate_id)
				max_predicate_id = symbols[i]->id;
		}
		return max_predicate_id + 1;
	}
	int getMaxNormId()
	{
		/*if(max_predicate_id == 0)
			setMaxPredicateId(symbols.size());
		return max_predicate_id;*/
		int max_norm_id = -1;
		for(unsigned int i=0;i<symbols.size();i++)
		{
			if(symbols[i]->normParentId > max_norm_id)
				max_norm_id = symbols[i]->normParentId;
		}
		return max_norm_id + 1;
	}
	void preprocessEvidenceWMPTP();
	void putWeightsOnClauses();
	void clearData()
	{
		symbols.clear();
		clauses.clear();
		formulas.clear();
		maxDegree = -1;
		max_predicate_id = 0;
	}
	void setMLNPoperties()
	{
		//max degree of LvrMLN
		int maxterms = 0;
		for(unsigned int i=0;i<symbols.size();i++)
		{
			if(symbols[i]->variable_types.size() > maxterms)
				maxterms = symbols[i]->variable_types.size();
		}
		maxDegree = maxterms;
	}
	int getMaxDegree()
	{
		if(maxDegree == -1)
		{
			setMLNPoperties();
		}
		return maxDegree;
	}
	~LvrMLN();
	LvrMLN():max_predicate_id(0),maxDegree(-1){}
	void preprocessEvidence(vector<vector<int> > evidenceIntRep);
	void preprocessEvidence(vector<Atom*> evidenceAtoms);
	void doPreprocessEvidence(map<int,int> evidenceHashMap);
	void reduceEvidence();
	int max_predicate_id;
	int maxDegree;
	LogDouble conversionFactor;
};



#endif /* MLN_H_ */
