#ifndef __LPTPTREE__
#define __LPTPTREE__
#include "lvrmln.h"
#include "cleanuputils.h"

enum NodeType
{
	ESPLIT,
	EDECOMPOSER,
	EDISJOINT
};

struct LAssignment
{
	vector<Atom*> atoms;
	vector<int> assignment;
	vector<LogDouble> probabilities;
	~LAssignment()
	{
		//cleanup(atoms);
	}
	void cleanAssignments()
	{
		cleanup(atoms);
		assignment.clear();
		probabilities.clear();
	}
};
struct LPTPNode
{
	int type;
	Atom* atom;
	int powerFactor;
	vector<LogDouble> zValues;
	vector<LogDouble> binCoeff;
	LogDouble norm;
	vector<LPTPNode*> children;
	int id;
	int parentId;
	int parentBranchNo;
	int assignment;
	vector<LogDouble> cdf;
	LogDouble downValue;
	vector<LogDouble> branchWeights;
	LogDouble nodeValue;
	//key = hash(normid+termid),value = domainvalue
	vector<map<int,int> > decomposerMappings;
	vector<int> domainPowerFactors;
	//vector<vector<WClause*> > decomposedList;
	LPTPNode():downValue(LogDouble(1,false))
	{
	}
	~LPTPNode()
	{
		if(type==ESPLIT)
			delete atom;
	}
	
	LogDouble getProbability(int index)
	{
		if(type==ESPLIT)
		{
			return (zValues[index]*binCoeff[index])/norm;
		}
		return LogDouble(0,false);
	}

	void constructCDF()
	{
		if(type==ESPLIT)
		{
			cdf.clear();
			cdf.resize(zValues.size());
			cdf[0] = (zValues[0]*binCoeff[0])/norm;
			for(unsigned int i=1;i<zValues.size()-1;i++)
			{
				cdf[i]=cdf[i-1] + (zValues[i]*binCoeff[i])/norm;
			}
			cdf[zValues.size()-1]=LogDouble(1,false);			
		}
		
}


	void print()
	{
		cout<<"Id="<<id<<endl;
		cout<<"ParentId="<<parentId<<endl;
		cout<<"ParentBranchNo="<<parentBranchNo<<endl;
		
		if(type == 0)
		{
			atom->print(false);
			for(unsigned int kk=0;kk<zValues.size();kk++)
				cout<<binCoeff[kk]*zValues[kk]<<" ";
			cout<<endl;
			cout<<"Norm="<<norm<<endl;
		}
		else if(type == 1)
		{
			cout<<"PF="<<powerFactor<<endl;
		}
		else
		{
			cout<<"Decomp Node"<<endl;
		}
		cout<<"Children=";
		for(unsigned int i=0;i<children.size();i++)
			cout<<children[i]->id<<" ";
		cout<<"["<<nodeValue<<","<<downValue<<"]"<<endl;
		cout<<endl<<"*************"<<endl;
	}
};
struct LPTPTree
{
	vector<LPTPNode*> treeNodes;
	LAssignment sampledAssignments;
	~LPTPTree()
	{
	}
	void cleanTree()
	{
		cleanup(treeNodes);
		sampledAssignments.cleanAssignments();
	}
	void printAssignments()
	{
		for(unsigned int j=0;j<sampledAssignments.atoms.size();j++)
		{
			sampledAssignments.atoms[j]->print(false);
			cout<<"$$"<<sampledAssignments.assignment[j]<<endl;
		}
	}
	void makeCDF()
	{
		for(unsigned int i=0;i<treeNodes.size();i++)
			treeNodes[i]->constructCDF();
	}
	void printTree()
	{
		for(unsigned int i=0;i<treeNodes.size();i++)
			treeNodes[i]->print();
	}
};
#endif