#include "blockexactinference.h"
#include "cleanuputils.h"
#include "propositionalresolution.h"
#include "randomgenutil.h"
#ifndef _MSC_VER 
#include <unistd.h>
#endif

LogDouble LBlockExactInference::CNFWeight(vector<WClause*>& CNF)
{
	vector<Atom*> removedAtoms;
	set<int> completedIds;
	LogDouble totalWeight(1,false);
	vector<WClause*> remainingCNF;
	for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->satisfied)
		{
			for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
			{
				if(completedIds.count(CNF[i]->atoms[j]->symbol->id) == 0)
				{
					removedAtoms.push_back(CNF[i]->atoms[j]);
					completedIds.insert(CNF[i]->atoms[j]->symbol->id);
				}
				
			}
			//get weight of formula
			set<LvrTerm*> terms;
			for(unsigned int kk=0;kk<CNF[i]->atoms.size();kk++)
			{
				for(unsigned int mm=0;mm<CNF[i]->atoms[kk]->terms.size();mm++)
				{
					terms.insert(CNF[i]->atoms[kk]->terms[mm]);
				}
			}
			int sz=1;
			for(set<LvrTerm*>::iterator it=terms.begin();it!=terms.end();it++)
			{
				sz *= (*it)->domain.size();
			}
			//if(!CNF[i]->weight.is_zero)
			{
				LogDouble tmpL = CNF[i]->weight*LogDouble(sz,false);
				double tmp = exp(tmpL.value);
				LogDouble newWt(tmp,true);
				totalWeight *= newWt;
			}
		}
		else
		{
			remainingCNF.push_back(LvrMLN::create_new_clause(CNF[i]));
		}
		
	}
	int totalDontcare = 0;
	for(vector<Atom*>::iterator it=removedAtoms.begin();it!=removedAtoms.end();it++)
	{
		//find the number of groundings
		int sz = 1;
		for(unsigned int i=0;i<(*it)->terms.size();i++)
		{
			sz *= (*it)->terms[i]->domain.size();
		}
		totalDontcare += sz;
	}
	LogDouble d(totalDontcare,false);
	LogDouble out(1,false);
	LogDouble c(2,false);
	LogDouble::LDPower(c,totalDontcare,out);
	LogDouble tVal = totalWeight*out;
	cleanup(CNF);
	CNF = remainingCNF;

	return tVal;
}


bool LBlockExactInference::decomposeCNF(vector<WClause*>& CNF,int& powerFactor,vector<map<int,int> >& decomposerMappings,
	vector<int>& domainPowerFactors)
{
	vector<Decomposer*> decomposer_list;
	decomposer->find_decomposer(CNF,decomposer_list);
	
#ifdef __DEBUG_PRINT__
	cout<<"Decomposer={ ";
	for(unsigned int i=0;i<decomposer_list.size();i++)
	{
		decomposer_list[i]->print();
	}
	cout<<"}"<<endl;
#endif
	if(decomposer_list.size()==0)
		return false;
	else
	{
		//replace all decomposer terms with a constant from the domain
		powerFactor = 1;
		for(unsigned int i=0;i<decomposer_list.size();i++)
		{
			powerFactor*=decomposer_list[i]->decomposer_terms[0]->domain.size();
			int domSize = decomposer_list[i]->decomposer_terms[0]->domain.size();
			domainPowerFactors.push_back(domSize);
			for(unsigned int j=0;j<decomposer_list[i]->decomposer_terms.size();j++)
			{
				//make a copy of the pre-decomp domain
				vector<int> origDomain;
				for(unsigned int jj=0;jj<decomposer_list[i]->decomposer_terms[j]->domain.size();jj++)
					origDomain.push_back(decomposer_list[i]->decomposer_terms[j]->domain[jj]);
				decomposer_list[i]->decomposer_terms[j]->origDomain = origDomain;

				//Erase all but one grounding
				decomposer_list[i]->decomposer_terms[j]->domain.erase(decomposer_list[i]->decomposer_terms[j]->domain.begin(),decomposer_list[i]->decomposer_terms[j]->domain.begin()+decomposer_list[i]->decomposer_terms[j]->domain.size()-1);
				//insert the original size into the term
				decomposer_list[i]->decomposer_terms[j]->origDomainSize = domSize;
				
				map<int,int> decompMapping;
				for(map<int,int>::iterator it=decomposer_list[i]->norm_predicate_positions.begin();
					it!=decomposer_list[i]->norm_predicate_positions.end();it++)
				{
					vector<int> vals(2);
					vals[0] = it->first;vals[1]=it->second;
					unsigned int hash = LvrHashAlgorithm::DJBHashUS(vals);
					decompMapping.insert(pair<int,int>(hash,0));
				}
				decomposerMappings.push_back(decompMapping);
			}
		}
	}
	cleanup(decomposer_list);
	return true;
}

LogDouble LBlockExactInference::weightedModelCount(vector<WClause*>& CNF,int parentId,int parentBranchNo)
{
	if(nonPolyNodeFound)
	{
		cleanup(CNF);
		return LogDouble(0,false);
	}
	int completedCount=0;
	for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->atoms.size()==0 || CNF[i]->satisfied)
		{
			completedCount++;
		}
	}
	if(completedCount == CNF.size())
	{
		LogDouble weight = CNFWeight(CNF);
		cleanup(CNF);
		return weight;
	}
	vector<int> powerFactors;
	vector<bool> isDecomposedCNF;
	vector<vector<WClause*> > decomposedList;
	LClusterUtil::seperateDisjointCNF(CNF,decomposedList);
	//CNF can be cleaned up since split into clusters
	cleanup(CNF);
	LogDouble totalVal(1,false);

	bool isParentChanged=false;
	int changedParent;
	if(decomposedList.size() > 1)
	{
		lvrPTPTreeSampling->addTreeNode(EDISJOINT,id,parentId,parentBranchNo);
		id++;
		isParentChanged = true;
		changedParent = id-1;
	}


	for(unsigned int t=0;t<decomposedList.size();t++)
	{
		if(isParentChanged)
		{
		parentId = changedParent;
		parentBranchNo = t;
		}

		int powerFactor;
		vector<map<int,int> > decomposerMappings;
		vector<int> domainPowerFactors;
		bool isDecomposed = decomposeCNF(decomposedList[t],powerFactor,decomposerMappings,domainPowerFactors);
		if(isDecomposed)
		{
			int idToUse = id++;
			//use power rule
			LogDouble mcnt = weightedModelCount(decomposedList[t],idToUse,0);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			LPTPNode* node = lvrPTPTreeSampling->addTreeNode(EDECOMPOSER,idToUse,parentId,parentBranchNo,powerFactor);
			node->decomposerMappings = decomposerMappings;
			node->domainPowerFactors = domainPowerFactors;
			continue;
		}
		else
		{
			vector<vector<WClause*> > tempCNFList;
			Atom* tmpatom = heuristics->getAtomToSplit(decomposedList[t]);
			if(tmpatom==NULL)
			{	
				LogDouble wt1 = CNFWeight(decomposedList[t]);
				cleanup(decomposedList[t]);
				totalVal *= wt1;
				continue;
			}
			Atom* atom = LvrMLN::create_new_atom(tmpatom);
			int singletonIndex;
			bool singleton = atom->isSingletonAtom(singletonIndex);			
			//check for self joins
			bool selfJoined = false;
			for(unsigned int m=0;m<decomposedList[t].size();m++)
			{
				if(decomposedList[t].at(m)->isSelfJoinedOnAtom(atom))
				{
					selfJoined = true;
					break;
				}
			}
			if(atom->isConstant())
			{
				int idToUse = id++;
				//propositional atom
				if(!selfJoined)
				{
					PropositionalResolution::doPropositionalSplit(decomposedList[t],atom,tempCNFList);
					LogDouble w1 = weightedModelCount(tempCNFList[0],idToUse,0);
					LogDouble w2 = weightedModelCount(tempCNFList[1],idToUse,1);
					totalVal *= (w1+w2);
					vector<LogDouble> wts;
					wts.push_back(w1);
					wts.push_back(w2);
					LogDouble norm = w1+w2;
					vector<LogDouble> bincoefs;
					bincoefs.push_back(LogDouble(1,false));
					bincoefs.push_back(LogDouble(1,false));
					lvrPTPTreeSampling->addTreeNode(atom,wts,bincoefs,norm,ESPLIT,idToUse,parentId,parentBranchNo);
					delete atom;
					pCost += 2;
					cleanup(tempCNFList[0]);
					cleanup(tempCNFList[1]);
					cleanup(decomposedList[t]);

				}
				else
				{
					PropositionalResolution::doPropositionalSplitSelfJoins(decomposedList[t],atom,tempCNFList);
					LogDouble w1 = weightedModelCount(tempCNFList[0],idToUse,0);
					LogDouble w2 = weightedModelCount(tempCNFList[1],idToUse,1);
					totalVal *= (w1+w2);
					vector<LogDouble> wts;
					wts.push_back(w1);
					wts.push_back(w2);
					LogDouble norm = w1+w2;
					vector<LogDouble> bincoefs;
					bincoefs.push_back(LogDouble(1,false));
					bincoefs.push_back(LogDouble(1,false));
					lvrPTPTreeSampling->addTreeNode(atom,wts,bincoefs,norm,ESPLIT,idToUse,parentId,parentBranchNo);

					cleanup(tempCNFList[0]);
					cleanup(tempCNFList[1]);
					tempCNFList.clear();
					delete atom;
					cleanup(decomposedList[t]);
					pCost += 2;
				}
			}
			else if(singleton)//first order atom for which GBR applies
			{
				int idToUse = id++;
				int totalSize = atom->getNumberOfGroundings();
				LogDouble v;
				vector<LogDouble> zvals;
				vector<LogDouble> binceoffs;
				LogDouble norm;
				for(unsigned int i=0;i<=totalSize;i++)
				{
					vector<WClause*> clausesCopy;
					LvrMLN::copyAllClauses(decomposedList[t],clausesCopy);
					LvrSingletonNormPropagation::propagateNormalizedCNF(clausesCopy,
						atom,singletonIndex,i);
					LogDouble binCoeff = LogDouble::Binomial(totalSize,i);
					LogDouble mcnt = weightedModelCount(clausesCopy,idToUse,i);
					if(nonPolyNodeFound)
					{
						cleanup(clausesCopy);
						delete atom;
						cleanupItems(decomposedList,t);
						return LogDouble(0,false);
					}
					cleanup(clausesCopy);
					v = v + binCoeff*mcnt;
					zvals.push_back(mcnt);
					binceoffs.push_back(binCoeff);
					norm = norm + binCoeff*mcnt;
				}

				totalVal *= v;
				lvrPTPTreeSampling->addTreeNode(atom,zvals,binceoffs,norm,ESPLIT,idToUse,parentId,parentBranchNo);
				delete atom;
				cleanup(decomposedList[t]);
			}
			else
			{
				cleanupItems(decomposedList,t);
				nonPolyNodeFound = true;
				atom->print();
				delete atom;
				return LogDouble(0,false);
			}
		}
	}
	return totalVal;
}


int LBlockExactInference::mockWeightedModelCount(vector<WClause*>& CNF)
{
#ifndef _MSC_VER 
		//check if parent process has exited
		if(getppid()==1)
		{
			cout<<"Polled parent has exited.. hence terminating"<<getppid()<<endl;
			exit(0);
		}
#endif
	
	time_t curTime;
	time(&curTime);
	if(difftime(curTime,startTime) > MAXTIMEOUTCLUSTERING)
	{
		cout<<"Timeout limit reached"<<endl;
		cleanup(CNF);
		status = -2;
		return 0;
	}
	int completedCount=0;
	for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->atoms.size()==0 || CNF[i]->satisfied)
		{
			completedCount++;
		}
	}
	if(completedCount == CNF.size())
	{
		cleanup(CNF);
		return 1;
	}
	vector<int> powerFactors;
	vector<bool> isDecomposedCNF;
	vector<vector<WClause*> > decomposedList;
	LClusterUtil::seperateDisjointCNF(CNF,decomposedList);
	//CNF can be cleaned up since split into clusters
	cleanup(CNF);
	int cost = 1;
	for(unsigned int t=0;t<decomposedList.size();t++)
	{
		int powerFactor;
		vector<map<int,int> > decomposerMappings;
		vector<int> domainPowerFactors;
		bool isDecomposed = decomposeCNF(decomposedList[t],powerFactor,decomposerMappings,domainPowerFactors);
		if(isDecomposed)
		{
			//use power rule
			int ret = mockWeightedModelCount(decomposedList[t]);
			if(status!=0)
			{
				cleanupItems(decomposedList,t);
				return ret;
			}
			cost += ret*1;
			continue;
		}
		else
		{
			vector<vector<WClause*> > tempCNFList;
			Atom* tmpatom = heuristics->getAtomToSplit(decomposedList[t]);
			if(tmpatom==NULL)
			{	
				cleanup(decomposedList[t]);
				continue;
			}
			Atom* atom = LvrMLN::create_new_atom(tmpatom);
			int singletonIndex;
			bool singleton = atom->isSingletonAtom(singletonIndex);			
			//check for self joins
			bool selfJoined = false;
			
			for(unsigned int m=0;m<decomposedList[t].size();m++)
			{
				if(decomposedList[t].at(m)->isSelfJoinedOnAtom(atom))
				{
					selfJoined = true;
					break;
				}
			}
			if(atom->isConstant())
			{
				if(!selfJoined)
					PropositionalResolution::doPropositionalSplit(decomposedList[t],atom,tempCNFList);
				else
					PropositionalResolution::doPropositionalSplitSelfJoins(decomposedList[t],atom,tempCNFList);
				double r = LvRandomGenUtil::Instance()->getNormRand();
				int ret;
				if(r<0.5)
				{
					ret = mockWeightedModelCount(tempCNFList[0]);	
				}
				else
				{
					ret = mockWeightedModelCount(tempCNFList[1]);
				}

				cleanup(tempCNFList[0]);
				cleanup(tempCNFList[1]);
				tempCNFList.clear();
				cleanup(decomposedList[t]);
				cost += ret*2;
				if(status!=0)
				{
					delete atom;
					cleanupItems(decomposedList,t+1);
					return ret;
				}
				if(cost > baselineCost)
				{
					delete atom;
					cleanupItems(decomposedList,t+1);
					status = -1;
					return 0;
				}
				delete atom;
			}
			else if(singleton)//first order atom for which GBR applies
			{
				int sz = atom->getNumberOfGroundings()+1;
				int r = LvRandomGenUtil::Instance()->getRandomPosition(sz);
				LvrSingletonNormPropagation::propagateNormalizedCNF(decomposedList[t],atom,singletonIndex,r);
				int ret = mockWeightedModelCount(decomposedList[t]);
				cost += ret*sz;
				if(status!=0)
				{
					delete atom;
					cleanupItems(decomposedList,t);
					return ret;
				}
				if(cost > baselineCost)
				{
					status = -1;
					delete atom;
					cleanupItems(decomposedList,t);
					return -1;
				}
				delete atom;
				cleanup(decomposedList[t]);
			}				
			else
			{
				status = -2;
				delete atom;
				cleanupItems(decomposedList,t);
				return -2;
			}
		}
	}
	return cost;
}

LBlockExactInference::LBlockExactInference(LvrMLN& mln_): mln(mln_)
{
	decomposer = new LDecomposer(mln);
	heuristics = new LHeuristics(*decomposer,mln);
	lvrPTPTreeSampling = new LvrPTPTreeSampling();
	nonPolyNodeFound = false;
}

int LBlockExactInference::doExactInferenceOnCluster(vector<WClause*>& CNF)
{
	lvrPTPTreeSampling->cleanTree();
	id = 0;

	nonPolyNodeFound = false;
	weightedModelCount(CNF,-1,-1);
	if(nonPolyNodeFound)
	{
		return -1;
	}	
	return lvrPTPTreeSampling->startNewSampling();
}

int LBlockExactInference::doMockInferenceOnCluster(vector<WClause*>& CNF,int& lvpCost,int remainingCost)
{
	status = 0;
	time(&startTime);
	baselineCost = remainingCost;
	int ret = mockWeightedModelCount(CNF);
	lvpCost = ret;
	return status;
}
