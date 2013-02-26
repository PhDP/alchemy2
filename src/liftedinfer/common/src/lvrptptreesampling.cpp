#include "lvrptptreesampling.h"
#include "randomgenutil.h"
#include <assert.h>
#include "hashalgorithm.h"

bool compareChildren(LPTPNode* l1, LPTPNode* l2)
{
	if(l1->parentBranchNo < l2->parentBranchNo)
		return true;
	else
		return false;
}

LvrPTPTreeSampling::LvrPTPTreeSampling()
{
	ptpTree = new LPTPTree();
}

LvrPTPTreeSampling::~LvrPTPTreeSampling()
{
	delete ptpTree;
}

void LvrPTPTreeSampling::addTreeNode(Atom* atom,vector<LogDouble> weights,vector<LogDouble> binCoeffs,
	LogDouble norm,NodeType nodeType,int nodeId, int parentId, 
	int parentBranchNo,vector<LogDouble>* branchWeights,LogDouble nodeValue)
{
	assert(nodeType == ESPLIT);
	Atom* newAtom = LvrMLN::create_new_atom(atom);
	LPTPNode* nd = new LPTPNode();
	nd->atom = newAtom;
	nd->zValues = weights;
	nd->norm = norm;
	nd->binCoeff = binCoeffs;
	nd->type = nodeType;
	nd->id = nodeId;
	nd->parentId = parentId;
	nd->parentBranchNo = parentBranchNo;
	if(branchWeights!=NULL)
		nd->branchWeights = (*branchWeights);
	ptpTree->treeNodes.push_back(nd);
	nd->nodeValue = nodeValue;
}

LPTPNode* LvrPTPTreeSampling::addTreeNode(NodeType nodeType,int nodeId, int parentId, int parentBranchNo,int powerFactor,LogDouble nodeValue)
{
	LPTPNode* nd = new LPTPNode();
	nd->type = nodeType;
	nd->id = nodeId;
	nd->parentId = parentId;
	nd->parentBranchNo = parentBranchNo;
	nd->powerFactor = powerFactor;
	nd->nodeValue = nodeValue;
	ptpTree->treeNodes.push_back(nd);
	return nd;
}

void LvrPTPTreeSampling::sampleTree(LPTPNode* node,vector<map<int,int>* >& addedDecompMappings)
{
	if(node->type == ESPLIT)
	{
		//sample atom to obtain assignment
		//store assignment
		int nextInd = -1;
		//apply all added decomposer mappings
		if(addedDecompMappings.size()!=0)
		{
			double u = LvRandomGenUtil::Instance()->getNormRand();
			LogDouble lu(u,false);
			//cout<<"lu="<<lu<<endl;
			int index = lower_bound(node->cdf.begin(),node->cdf.end(),lu) - node->cdf.begin();
			int dum=0;
			if(dum==1)
			{
				index=0;
			}
			Atom* newAtom = LvrMLN::create_new_atom(node->atom);
			vector<int> val(2);
			val[0] = newAtom->symbol->normParentId;
			for(unsigned int i=0;i<addedDecompMappings.size();i++)
			{
				map<int,int>* dMap = addedDecompMappings[i];
				//check if a term has this mapping
				for(unsigned int j=0;j<newAtom->terms.size();j++)
				{
					val[1]=j;
					unsigned int key = LvrHashAlgorithm::DJBHashUS(val);
					map<int,int>::iterator it1 = (*dMap).find(key);
					if(it1==(*dMap).end())
						continue;
					if(it1->second >= newAtom->terms[j]->origDomain.size())
						continue;
					//change the domain to the one specified in the mapping
					newAtom->terms[j]->domain.clear();
					newAtom->terms[j]->domain.push_back(newAtom->terms[j]->origDomain[it1->second]);
				}

			}
			ptpTree->sampledAssignments.atoms.push_back(newAtom);
			ptpTree->sampledAssignments.assignment.push_back(index);
			//newAtom->print();cout<<endl;
			nextInd = index;
		}
		else
		{
			//no decomposed position, sample as for given value
			double u = LvRandomGenUtil::Instance()->getNormRand();
			LogDouble lu(u,false);
			int index = lower_bound(node->cdf.begin(),node->cdf.end(),lu) - node->cdf.begin();
			ptpTree->sampledAssignments.atoms.push_back(LvrMLN::create_new_atom(node->atom));
			ptpTree->sampledAssignments.assignment.push_back(index);
			//node->atom->print();cout<<endl;
			nextInd = index;
		}
		//take sampled branch child node
		if(nextInd < node->children.size())
			sampleTree(node->children[nextInd],addedDecompMappings);
		else
			return;
	}
	else if(node->type == EDECOMPOSER)
	{
		if(node->children.size()==0)
			return;

		for(unsigned int i=0;i<node->domainPowerFactors.size();i++)
		{
			int pf = node->domainPowerFactors[i];
			for(unsigned int j=0;j<pf;j++)
			{
				//increment the mapping for all downstream nodes
				for(map<int,int>::iterator it = node->decomposerMappings[i].begin();
					it!=node->decomposerMappings[i].end();it++)
				{
					it->second = j;
				}
				//check if the mapping already exists
				if(find(addedDecompMappings.begin(),addedDecompMappings.end(),&node->decomposerMappings[i])
					== addedDecompMappings.end())
					addedDecompMappings.push_back(&node->decomposerMappings[i]);
				//has exactly 1 child
				for(unsigned int jj=0;jj<node->children.size();jj++)
					sampleTree(node->children[jj],addedDecompMappings);
			}
		}
	}
	else if(node->type == EDISJOINT)
	{
		//sample all children
		for(unsigned int i=0;i<node->children.size();i++)
		{
			sampleTree(node->children[i],addedDecompMappings);
		}
	}

}

LPTPTree& LvrPTPTreeSampling::getTreeRef()
{
	return (*ptpTree);
}

void LvrPTPTreeSampling::cleanTree()
{
	ptpTree->cleanTree();
}

void LvrPTPTreeSampling::traverseTree(int id)
{
	//search for the id
	int index = -1;
	for(unsigned int i=0;i<ptpTree->treeNodes.size();i++)
	{
		if(ptpTree->treeNodes[i]->id == id)
		{
			index = i;
			break;
		}
	}
	if(index!=-1)
	{
		//search for all nodes that have this node as its parent
		vector<LPTPNode*> children;
		for(unsigned int i=0;i<ptpTree->treeNodes.size();i++)
		{
			if(ptpTree->treeNodes[i]->parentId == ptpTree->treeNodes[index]->id)
			{
				children.push_back(ptpTree->treeNodes[i]);
			}
		}
		//sort children according to parent branch numbers
		sort(children.begin(),children.end(),compareChildren);
		ptpTree->treeNodes[index]->children = children;
		for(unsigned int i=0;i<children.size();i++)
		{
			traverseTree(children[i]->id);
		}
	}
}

LPTPNode* LvrPTPTreeSampling::getRoot()
{
	for(unsigned int i=0;i<ptpTree->treeNodes.size();i++)
	{
		if(ptpTree->treeNodes[i]->id == 0)
		{
			return ptpTree->treeNodes[i];
		}
	}
	return NULL;
}

int LvrPTPTreeSampling::startNewSampling()
{
	traverseTree(0);
	LPTPNode* root;
	root = getRoot();
	if(root==NULL)
		return 0;
	ptpTree->makeCDF();
	ptpTree->sampledAssignments.assignment.clear();
	ptpTree->sampledAssignments.atoms.clear();
	vector<map<int,int>* > addedDecompMappings;
	sampleTree(root,addedDecompMappings);
	return 1;
}

void LvrPTPTreeSampling::print()
{
	ptpTree->printTree();
}
