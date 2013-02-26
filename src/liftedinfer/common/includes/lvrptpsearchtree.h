#ifndef __LVRPTP_SEARCH_TREE
#define __LVRPTP_SEARCH_TREE
#include "lvrmln.h"
#include "cleanuputils.h"
#include "queryupdater.h"

enum LvrPTPNodeType
{
	ESPLITTER,
	EPOWER,
	EAND,
	ELEAF
};

struct LvrPTPNode
{
	LvrPTPNodeType type;
	Atom* atom;
	int powerFactor;
	vector<LvrPTPNode*> children;
	int parentBranchNo;
	LogDouble downValue;
	vector<LogDouble> branchWeights;
	LogDouble nodeValue;
	set<int>* dontCareQueries;
	LvrPTPNode(LvrPTPNodeType type_):type(type_),downValue(LogDouble(1,false)),dontCareQueries(NULL)
	{
	}
	~LvrPTPNode()
	{
		if(type==ESPLITTER)
			delete atom;
		if(dontCareQueries)
			delete dontCareQueries;
	}

	void print()
	{
		cout<<"ParentBranchNo="<<parentBranchNo<<endl;
		cout<<"Pointer value="<<this<<endl;
		if(type == ESPLITTER)
		{
			atom->print(false);
			for(unsigned int i=0;i<atom->terms.size();i++)
				cout<<"#"<<atom->terms[i]->origDomain.size();
		}
		else if(type == EPOWER)
		{
			cout<<"PowerNode, PF="<<powerFactor<<endl;
		}
		else if(type==EAND)
		{
			cout<<"AND Node"<<endl;
		}
		else if(type==ELEAF)
		{
			cout<<"Leaf node ";
			if(dontCareQueries && dontCareQueries->size() > 0)
				cout<<"#"<<dontCareQueries->size();
			cout<<endl;
		}
		cout<<"Children pointers=";
		for(unsigned int i=0;i<children.size();i++)
			cout<<children[i]<<" ";
		cout<<"["<<nodeValue<<","<<downValue<<"]"<<endl;
		cout<<endl;
	}
	/*
	void print()
	{
		if(type == ESPLITTER)
		{
			cout<<"SN=";
			atom->print(false);
			for(unsigned int i=0;i<atom->terms.size();i++)
				cout<<"#"<<atom->terms[i]->origDomain.size();
		}
		else if(type == EPOWER)
		{
			cout<<"PF="<<powerFactor;
		}
		else if(type==EAND)
		{
			cout<<"AND Node";
		}
		else if(type==ELEAF)
		{
			cout<<"Leaf";
			if(dontCareQueries && dontCareQueries->size() > 0)
				cout<<"#"<<dontCareQueries->size();
			//cout<<endl;
		}
		cout<<"["<<nodeValue<<","<<downValue<<"]";
	}
*/
	static LvrPTPNode* create_new_node(LvrPTPNodeType type_)
	{
		LvrPTPNode* temp = new LvrPTPNode(type_);
		return temp;
	}
};

struct LvrPTPSearchTree
{
	LvrPTPSearchTree();
	~LvrPTPSearchTree();
	void cleanTree();
	void printTree(LvrPTPNode* node);
	LvrPTPNode* createNewNode(LvrPTPNodeType type,LvrPTPNode* parent,int parentbranchno);
	void setRoot(LvrPTPNode* root_);
	bool isRootSet();
	LvrPTPNode* getRoot();
	void getHashValues(Atom* atom,set<int>& hashes);
	void getHashValuesNPAtom(Atom* atom,vector<int>& hashes);
	void downPropagationMAR(LvrPTPNode* node);
	void updateQueries(LvrPTPNode* node);
private:
	LvrPTPNode* root;
};

#endif
