#include "lvrptpsearchtree.h"
#include "mathutils.h"

LvrPTPSearchTree::LvrPTPSearchTree():root(NULL)
{
}

LvrPTPSearchTree::~LvrPTPSearchTree()
{
}

void LvrPTPSearchTree::cleanTree()
{
	//TODO!! traverse and clean
	root=NULL;
}

void LvrPTPSearchTree::printTree(LvrPTPNode* node)
{
	if(node==NULL)
		return;
	node->print();
	for(unsigned int i=0;i<node->children.size();i++)
	{
		printTree(node->children[i]);
	}
}

LvrPTPNode* LvrPTPSearchTree::createNewNode(LvrPTPNodeType type,LvrPTPNode* parent,int parentbranchno)
{
	LvrPTPNode* newnode = LvrPTPNode::create_new_node(type);
	if(parent==NULL)
	{
		root=newnode;
	}
	else
	{
		parent->children.push_back(newnode);
	}
	newnode->parentBranchNo = parentbranchno;
	return newnode;
}

void LvrPTPSearchTree::setRoot(LvrPTPNode* root_)
{
	if(root!=NULL)
		root = root_;
}
bool LvrPTPSearchTree::isRootSet()
{
	if(!root)
		return false;
	else
		return true;
}
LvrPTPNode* LvrPTPSearchTree::getRoot()
{
	return root;
}

void LvrPTPSearchTree::getHashValues(Atom* atom,set<int>& hashes)
{
	if(!LvrQueryUpdater::Instance()->isNormIdInQuery(atom->symbol->normParentId))
		return;
	Atom* tmpAtom = LvrMLN::create_new_atom(atom);
	for(unsigned int i=0;i<tmpAtom->terms.size();i++)
	{
		//replace domain by original domain
		if(tmpAtom->terms[i]->origDomain.size() > 0)
		{
			tmpAtom->terms[i]->domain.clear();
			tmpAtom->terms[i]->domain = tmpAtom->terms[i]->origDomain;
		}
	}
	if(tmpAtom->isConstant())
	{
		hashes.insert(LvrHashAlgorithm::convertToHash(tmpAtom));
		delete tmpAtom;
		return;
	}
	vector<vector<int> > permutedList;
	LvrPermutations::permuteTerms(tmpAtom->terms,permutedList);
	vector<int> hvalue(tmpAtom->terms.size()+1);
	hvalue[0]=tmpAtom->symbol->normParentId;
	for(unsigned int i=0;i<permutedList.size();i++)
	{
		for(unsigned int j=0;j<permutedList[i].size();j++)
		{
			hvalue[j+1] = permutedList[i].at(j);
		}
		hashes.insert(LvrHashAlgorithm::DJBHash(hvalue));
	}
	delete tmpAtom;
}

void LvrPTPSearchTree::getHashValuesNPAtom(Atom* atom,vector<int>& hashes)
{
	if(!LvrQueryUpdater::Instance()->isNormIdInQuery(atom->symbol->normParentId))
		return;
	Atom* tmpAtom = LvrMLN::create_new_atom(atom);
	for(unsigned int i=0;i<tmpAtom->terms.size();i++)
	{
		//replace domain by original domain
		if(tmpAtom->terms[i]->origDomain.size() > 0)
		{
			tmpAtom->terms[i]->domain.clear();
			tmpAtom->terms[i]->domain = tmpAtom->terms[i]->origDomain;
		}
	}
	if(tmpAtom->isConstant())
	{
		hashes.push_back(LvrHashAlgorithm::convertToHash(tmpAtom));
		delete tmpAtom;
		return;
	}
	vector<vector<int> > permutedList;
	LvrPermutations::permuteTerms(tmpAtom->terms,permutedList);
	vector<int> hvalue(tmpAtom->terms.size()+1);
	hvalue[0]=tmpAtom->symbol->normParentId;
	for(unsigned int i=0;i<permutedList.size();i++)
	{
		for(unsigned int j=0;j<permutedList[i].size();j++)
		{
			hvalue[j+1] = permutedList[i].at(j);
		}
		hashes.push_back(LvrHashAlgorithm::DJBHash(hvalue));
	}
	delete tmpAtom;
}


void LvrPTPSearchTree::downPropagationMAR(LvrPTPNode* node)
{
	if(node==NULL)
		return;
	if(node->type == ESPLITTER)
	{
		//for each child dv' = dv*wt
		for(unsigned int i=0;i<node->children.size();i++)
		{
			if(node->children.size()-1 == node->atom->getNumberOfGroundings())
			{
				//add weight for binomial split
				LogDouble bin;
				bin = LogDouble::Binomial(node->children.size()-1,i);
				node->children[i]->downValue = node->downValue*node->branchWeights[i]*bin;
			}
			else
			{
				//full grounding do not add coeffecient
				node->children[i]->downValue = node->downValue*node->branchWeights[i];
			}
			downPropagationMAR(node->children[i]);
		}
	}
	else if(node->type == EPOWER)
	{
		//for each child dv' = dv^(|x|-1)
		for(unsigned int i=0;i<node->children.size();i++)
		{
			node->children[i]->downValue = LogDouble(1,false);
			LogDouble::LDPower(node->children[i]->nodeValue,node->powerFactor-1,node->children[i]->downValue);
			node->children[i]->downValue = node->children[i]->downValue*node->downValue;
			downPropagationMAR(node->children[i]);
		}
	}
	else if(node->type == EAND)
	{
		for(unsigned int i=0;i<node->children.size();i++)
		{
			node->children[i]->downValue = LogDouble(1,false);
			for(unsigned int j=0;j<node->children.size();j++)
			{
				if(i==j)
					continue;
				node->children[i]->downValue = node->children[i]->downValue*node->children[j]->nodeValue;
			}
			node->children[i]->downValue = node->children[i]->downValue*node->downValue;
			downPropagationMAR(node->children[i]);
		}	
	}
	else if(node->type == ELEAF)
	{
		if(node->dontCareQueries)
		{
			if(node->dontCareQueries->size() > 0)
			{
				LogDouble dontcareval = node->nodeValue*node->downValue/LogDouble(2,false);
				LvrQueryUpdater::Instance()->updateExactQueryWeights((*node->dontCareQueries),dontcareval);
			}
		}
	}
}
	
void LvrPTPSearchTree::updateQueries(LvrPTPNode* node)
{
	if(node==NULL)
		return;
	if(node->type == ESPLITTER)
	{
		if(node->atom->getNumberOfGroundings()!=node->children.size()-1)
		{
			//full grounding
			vector<int> hashes;
			getHashValuesNPAtom(node->atom,hashes);
			if(hashes.size() > 0)
			{
				int numGoundings = node->atom->getNumberOfGroundings();
				for(unsigned int i=1;i<node->children.size();i++)
				{
					//convert i to binary sequence
					vector<int> binvector;
					LMathUtils::Instance()->toBinary(i,binvector,numGoundings);
					for(unsigned int k=0;k<binvector.size();k++)
					{
						//update only if the atom has a true assignment
						if(binvector[k]==1)
						{
							set<int> tmphash;
							tmphash.insert(hashes[k]);
							LogDouble value = node->children[i]->nodeValue*node->children[i]->downValue;
							LvrQueryUpdater::Instance()->updateExactQueryWeights(tmphash,value);
						}
					}
				}
			}
		}
		else
		{
			//update all the query predicates;
			set<int> queryHashes;
			getHashValues(node->atom,queryHashes);
			if(queryHashes.size() > 0)
			{
				//binomial split
				LogDouble combinedVal(0,false);
				//compute the value for this node and set in query
				int domainSize = node->children.size()-1;
				LogDouble lvnorm(domainSize,false);
				for(unsigned int i=1;i<node->children.size();i++)
				{
					if(node->children[i]->nodeValue.is_zero || node->children[i]->downValue.is_zero)
						continue;
					combinedVal += (node->children[i]->nodeValue*node->children[i]->downValue)*LogDouble(i,false)/lvnorm;
				}
				LvrQueryUpdater::Instance()->updateExactQueryWeights(queryHashes,combinedVal);
			}
		}
	}
	for(unsigned int j=0;j<node->children.size();j++)
	{
		updateQueries(node->children[j]);	
	}
}
