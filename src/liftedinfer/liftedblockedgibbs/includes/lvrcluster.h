#ifndef __LVRCLUSTER__H
#define __LVRCLUSTER__H
#include "lvrmln.h"
#include "cleanuputils.h"
#include <algorithm>
using namespace std;

struct LVAssignment
{
	vector<int> grounding;
	int nTrueValues;
	int nFalseValues;
	LVAssignment(Atom* gAtom, int val)
	{
		grounding.resize(gAtom->terms.size());
		int totSize = 1;
		for(unsigned int i=0;i<gAtom->terms.size();i++)
		{
			if(gAtom->terms[i]->domain.size() == 1)
				grounding[i] = gAtom->terms[i]->domain[0];
			else
				//for the lifted terms, store -1*domainsize
				grounding[i] = gAtom->terms[i]->domain.size()*-1;
			totSize *= gAtom->terms[i]->domain.size();
		}
		nTrueValues = val;
		nFalseValues = totSize - nTrueValues;
	}
	void print()
	{
		cout<<"[(";
		for(unsigned int i=0;i<grounding.size();i++)
		{
			cout<<grounding[i]<<",";
		}
		cout<<") "<<nTrueValues<<"] ";
	}
};

struct LTermShare
{
	vector<bool> shared;
	//do not delete pointers
	vector<vector<Atom*> > sharedAtoms;
	LTermShare(int sz)
	{
		shared.resize(sz);
		sharedAtoms.resize(sz);
	}
	void print()
	{
		for(int i=0;i<shared.size();i++)
		{
			if(shared[i])
			{
				for(int j=0;j<sharedAtoms[i].size();j++)
					sharedAtoms[i].at(j)->print(false);
			}
			cout<<" :: ";
		}
		cout<<endl;
	}
	~LTermShare()
	{
		for(unsigned int i=0;i<sharedAtoms.size();i++)
		{
			sharedAtoms[i].clear();
		}
		sharedAtoms.clear();
	}
};

struct LVRCluster
{
	int clusterId;
	vector<Atom*> elements;
	vector<LVRCluster*> MB;
	vector<LTermShare*> sharedTerms;
	LvrMLN& mln;
	vector<vector<LVAssignment*> > lAssignments;
	vector<vector<int> > sharedLiftedPositions;
	vector<vector<int> > sharedGroundedPositions;
	vector<int> totalTrues;
	int numSingletonSplits;
	int numGroundedAtoms;
	int approxPTPCost;
	int evidenceReductionCost;
	
	LVRCluster(int clusterId_,LvrMLN& mln_):clusterId(clusterId_),mln(mln_){}
	LVRCluster& operator=(const LVRCluster& other);
	LVRCluster(const LVRCluster& other);
	int inCluster(Atom* atom);
	int getEquivalentElementIndex(Atom* atom);
	void initializeAssignments();
	int createNewPredicateId(int elementId,Atom* groundedAtom);
	void renamePredicates(vector<WClause*>& clauses);
	~LVRCluster()
	{
		cleanup(elements);
		while(!lAssignments.empty())
		{
			cleanup(lAssignments.back());
			lAssignments.pop_back();
		}
		cleanup(sharedTerms);
		MB.clear();
	}
	void print()
	{
		for(unsigned int i=0;i<elements.size();i++)
			elements[i]->print(false);
		cout<<endl;
	}
	void printAssignment()
	{
		for(unsigned int i=0;i<lAssignments.size();i++)
		{
			for(unsigned int j=0;j<lAssignments[i].size();j++)
			{
				lAssignments[i].at(j)->print();
			}
			cout<<endl;
		}
		cout<<endl;
	}
	void printSharings()
	{
		for(unsigned int i=0;i<sharedTerms.size();i++)
		{
			cout<<elements[i]->symbol->symbol.c_str()<<"(";
			for(unsigned int j=0;j<sharedTerms[i]->shared.size();j++)
			{
				if(sharedTerms[i]->shared[j])
				{
					cout<<"Shared,";
				}
				else
				{
					cout<<"UnShared,";
				}
			}
			cout<<")";
		}
		cout<<endl;
	}
	void addElement(Atom* atom)
	{
		elements.push_back(LvrMLN::create_new_atom(atom));
		//initializeAssignments();
	}
	void removeElements(vector<Atom*> atoms);
	Atom* createGroundedAtom(int elementIndex,vector<bool> sharedPos, vector<int> grounding);
	WClause* groundClusterInClause(WClause* clause, vector<vector<int> > grounding);
	LVRCluster* cloneClusterElements();
	int totalSpace()
	{
		//return assignments.size();
		int sz = 0;
		for(unsigned int i=0;i<lAssignments.size();i++)
			sz += lAssignments[i].size();
		return sz;
	}
	
	bool isInClause(WClause* clause)
	{
		for(unsigned int i=0;i<clause->atoms.size();i++)
		{
			int ind = inCluster(clause->atoms[i]);
			if(ind!=-1)
				return true;
		}
		return false;
	}
	//void groundSharedTermsInCNF(vector<WClause*> Clauses,vector<WClause*>& groundedClauses);
	void groundSharedTermsInCNF(vector<WClause*> &clauses);
	void reduceCNFByEvidence(vector<WClause*>& clauses);
	void getNumTrueAssignments(int elementIndex, Atom* atom,int& nLVTrue,int& nLVFalse);
	void doReduceCNFByEvidence(int elementIndex,vector<WClause*>& clauses,vector<WClause*>& groundedClauses);
	void groundAllSharedTermsInCNF(vector<WClause*> clauses,vector<WClause*>& groundedClauses);
	bool doExactInference();
	void partitionPositions();
	void resetAllAssignments();
	void setAllDontCareAssignments();
	void doReduceClauseWithSelfJoinsByEvidence(int elementIndex,WClause* clause,vector<WClause*>& groundedClauses);
	bool isUniConstantCluster()
	{
		if(elements.size()==1 && elements[0]->isConstant())
			return true;
		else
			return false;
	}
};

#endif
