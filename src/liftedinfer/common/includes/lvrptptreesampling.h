#ifndef __LVRPTPTREE_SAMPLING__
#define __LVRPTPTREE_SAMPLING__
#include "ptptree.h"

struct LvrPTPTreeSampling
{
	LvrPTPTreeSampling();
	~LvrPTPTreeSampling();
	void addTreeNode();
	void sampleTree();
	LPTPTree& getTreeRef();
	void cleanTree();
	void traverseTree(int id);
	void sampleTree(LPTPNode* node,vector<map<int,int>* >& addedDecompMappings);
	void addTreeNode(Atom* atom,vector<LogDouble> weights,vector<LogDouble> binCoeffs,
		LogDouble norm,NodeType nodeType,int nodeId, int parentId, int parentBranchNo,
		vector<LogDouble>* branchWeights = NULL,LogDouble nodeValue = LogDouble(0,false));
	LPTPNode* addTreeNode(NodeType nodeType,int nodeId, int parentId, int parentBranchNo,int powerFactor = -1, LogDouble nodeValue=LogDouble(1,false));
	int startNewSampling();
	LPTPNode* getRoot();
	void print();
private:
	LPTPTree* ptpTree;
};
#endif
