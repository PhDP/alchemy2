#include <time.h>
#include "clustercreator.h"	
#include "normalizer.h"
#include "randomgenutil.h"
#include <fstream>
#include <sstream>
using namespace std;
#include <sys/types.h>
#ifndef _MSC_VER 
#include <unistd.h>
#endif

bool sortFunc (LTuple* l1,LTuple* l2) { return (l1->CF > l2->CF); }


int LClusterCreator::getClusterIndex(LVRCluster* cluster)
{
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		if(LVRClusterList[i]->clusterId == cluster->clusterId)
		{
			return i;
		}
	}
	return -1;
}

void LClusterCreator::mergeCluster(int clusterInd1,int clusterInd2)
{
	for(unsigned int i=0;i<LVRClusterList[clusterInd2]->elements.size();i++)
		LVRClusterList[clusterInd1]->addElement(LVRClusterList[clusterInd2]->elements[i]);
	LVRClusterList[clusterInd1]->initializeAssignments();
	removeItem(LVRClusterList,clusterInd2);
}

void LClusterCreator::getAtomsMatchingTerm(Atom* atom, int termIndex, vector<Atom*>& matchingMLNAtoms)
{
	set<int> completedIds;
	for(unsigned int i=0;i<mln->clauses.size();i++)
	{
		LvrTerm* termToMatch = NULL;
		int atomInd = -1;
		for(unsigned int j=0;j<mln->clauses[i]->atoms.size();j++)
		{
			if(completedIds.count(mln->clauses[i]->atoms[j]->symbol->id)>0)
				continue;
			if(mln->clauses[i]->atoms[j]->symbol->id == atom->symbol->id)
			{
				//capture the term
				termToMatch = mln->clauses[i]->atoms[j]->terms[termIndex];
				atomInd = j;
				break;
			}
		}
		if(termToMatch == NULL)
			continue;
		for(unsigned int j=0;j<mln->clauses[i]->atoms.size();j++)
		{
			if(j==atomInd)
				continue;
			for(unsigned int k=0;k<mln->clauses[i]->atoms[j]->terms.size();k++)
			{
				if(mln->clauses[i]->atoms[j]->terms[k] == termToMatch &&
					completedIds.count(mln->clauses[i]->atoms[j]->symbol->id)==0)
				{
					completedIds.insert(mln->clauses[i]->atoms[j]->symbol->id);
					matchingMLNAtoms.push_back(mln->clauses[i]->atoms[j]);
				}
			}
		}
	}
}

void LClusterCreator::computeClusterLiftedSharings()
{
	//for each cluster set which clusters it shares a term with
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		for(unsigned int j=0;j<LVRClusterList[i]->sharedTerms.size();j++)
		{
			for(unsigned int k=0;k<LVRClusterList[i]->sharedTerms[j]->shared.size();k++)
			{
				if(LVRClusterList[i]->sharedTerms[j]->shared[k])
				{
					vector<Atom*> matchingAtoms;
					getAtomsMatchingTerm(LVRClusterList[i]->elements[j],k,matchingAtoms);
					LVRClusterList[i]->sharedTerms[j]->sharedAtoms[k] = matchingAtoms;
				}
			}
		}
	}
	/*
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		for(unsigned int j=0;j<LVRClusterList[i]->sharedTerms.size();j++)
		{
			LVRClusterList[i]->sharedTerms[j]->print();
		}
	}
	*/
}
void LClusterCreator::computeClusterMB()
{
	vector<vector<int> > clauseIndexes(LVRClusterList.size());
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		vector<int> clauseIndex;
		for(unsigned int j=0;j<mln->clauses.size();j++)
		{
			for(unsigned int k=0;k<mln->clauses[j]->atoms.size();k++)
			{
				int ind = LVRClusterList[i]->inCluster(mln->clauses[j]->atoms[k]);
				if(ind!=-1)
				{
					clauseIndex.push_back(j);
					break;
				}
			}
		}
		clauseIndexes[i] = clauseIndex;
	}
	for(unsigned int i=0;i<clauseIndexes.size();i++)
	{
		LVRClusterList[i]->MB.clear();
		for(unsigned int j=0;j<clauseIndexes.size();j++)
		{
			if(i==j)
				continue;
			vector<int> commonSet(clauseIndexes[i].size());
			for(unsigned int k=0;k<commonSet.size();k++)
				commonSet[k] = -1;
			set_intersection(clauseIndexes[i].begin(),clauseIndexes[i].end(),clauseIndexes[j].begin(),clauseIndexes[j].end(),commonSet.begin());
			bool emptyIntersect = true;
			for(unsigned int k=0;k<commonSet.size();k++)
			{
				if(commonSet[k] != -1)
				{
					emptyIntersect = false;
					break;
				}
			}
			if(!emptyIntersect)
			{
				LVRClusterList[i]->MB.push_back(LVRClusterList[j]);
			}
		}
	}
	//for(unsigned int i=0;i<LVRClusterList.size();i++)
	//{
		//for(unsigned int j=0;j<LVRClusterList[i]->MB.size();j++)
			//LVRClusterList[i]->MB[j]->print();
		//cout<<endl;
	//}
	computeClusterLiftedSharings();
}
//void LClusterCreator::makeClusterPairs(vector<vector<int> >& pairedClusterIndexes,vector<LTuple>& clusterCombinations)
void LClusterCreator::makeClusterPairs(vector<LTuple*>& clusterCombinations)
{
	vector<vector<int> > completedPairs;

	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		int index_i = getClusterIndex(LVRClusterList[i]);
		//pair with every cluster in its MB
		for(unsigned int j=0;j<LVRClusterList[i]->MB.size();j++)
		{
			int index_j = getClusterIndex(LVRClusterList[i]->MB[j]);
			vector<int> np(2);
			np[0] = index_i;
			np[1] = index_j;
			bool completed = false;
			for(unsigned int jj=0;jj<completedPairs.size();jj++)
			{
				if((completedPairs[jj].at(0) == np[0] && completedPairs[jj].at(1) == np[1])
				|| (completedPairs[jj].at(0) == np[1] && completedPairs[jj].at(1) == np[0]))
				{
					completed = true;
					break;
				}
			}
			if(!completed)
			{
				//pairedClusterIndexes.push_back(np);
				LTuple* lt = new LTuple();
				lt->clusterIndexes = np;
				clusterCombinations.push_back(lt);
				completedPairs.push_back(np);
			}
		}
	}
	//each pair has a weight = CF
	for(unsigned int i=0;i<clusterCombinations.size();i++)
	{
		int cf = computeCouplingFactor(clusterCombinations[i]->clusterIndexes[0],clusterCombinations[i]->clusterIndexes[1]);
		clusterCombinations[i]->CF = cf;
	}

	//for(unsigned int i=0;i<pairedClusterIndexes.size();i++)
		//cout<<"("<<pairedClusterIndexes[i].at(0)<<","<<pairedClusterIndexes[i].at(1)<<")"<<endl;
}
int LClusterCreator::numCoOccurrences(Atom* a1,Atom* a2)
{
	int sz = 0;
	for(unsigned int i=0;i<mln->clauses.size();i++)
	{
		if(mln->clauses[i]->isAtomInClause(a1) &&
			mln->clauses[i]->isAtomInClause(a2))
			sz += mln->clauses[i]->getNumberOfGroundedClauses();
	}
	return sz;
}

int LClusterCreator::computeCouplingFactor(int clusterIndex1,int clusterIndex2)
{
	//check how many times do elements of the 2 clusters occur together in formulas
	int cf = 0;
	for(unsigned int i=0;i<LVRClusterList[clusterIndex1]->elements.size();i++)
	{
		for(unsigned int j=0;j<LVRClusterList[clusterIndex2]->elements.size();j++)
		{
			cf += numCoOccurrences(LVRClusterList[clusterIndex1]->elements[i],
				LVRClusterList[clusterIndex2]->elements[j]);
		}
	}
	return cf;
}


void LClusterCreator::startRandomClustering(LvrParams* params)
{
	int level = 0;
	initializeClusters();
	int baselinePTPCost = 0;
	int baselineSpaceCost = 0;
	if(params->baselineCostMultiplicativeFactor < 1)
		params->baselineCostMultiplicativeFactor = 1.25;
	cout<<"Clustering Process::Setting baseline.. level "<<level<<endl;
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		vector<WClause*> newClauses;
		addPTPEvidence(i,newClauses);
		int lvpCost;
		if(lp->runMockPTP(newClauses,LVRClusterList[i],lvpCost,RAND_MAX) < 0)
		{
			cout<<"Cluster Creator::Non Poly node, terminate current clustering"<<endl;
			return;
		}
		baselinePTPCost += (double)lvpCost*params->baselineCostMultiplicativeFactor;
		baselineSpaceCost += LVRClusterList[i]->totalSpace();
		newClauses.clear();
		//LVRClusterList[0]->printAssignment();
	}
	int bestCouplingFactor = 0;
	int iter=0;
	while(1)
	{

#ifndef _MSC_VER 
		if(getppid()==1)
		{
			cout<<"Clustering Process:: Parent process has exited, hence terminating"<<endl;
			exit(0);
		}
		//IF RUNNING in parallel mode sleep during each iteration, give priority to sampler
		if(params->operationMode == EPARALLEL)
			sleep(SLEEPTIME);
#endif

		//cout<<iter++<<endl;
		vector<LTuple*> clusterCombinations;
		makeClusterPairs(clusterCombinations);
		if(clusterCombinations.size() == 0)
			break;
		vector<LVRCluster*> LVRClusterList_orig(LVRClusterList.size());
		int ind = LvRandomGenUtil::Instance()->getRandomPosition(clusterCombinations.size());

		for(unsigned int jj=0;jj<LVRClusterList.size();jj++)
		{
			LVRClusterList_orig[jj] = LVRClusterList[jj]->cloneClusterElements();
		}

		mergeCluster(clusterCombinations[ind]->clusterIndexes.at(0),clusterCombinations[ind]->clusterIndexes.at(1));
		computeClusterMB();
		int currPTPCost = 0;
		int currSpaceCost = 0;
		int couplingFactor = 0;
		bool clusteringNP = false;
		bool clusterCostExceeded = false;

		for(unsigned int jj=0;jj<LVRClusterList.size();jj++)
		{
				vector<WClause*> newClauses;
				addPTPEvidence(jj,newClauses);
				if(currPTPCost > baselinePTPCost)
				{
					clusterCostExceeded = true;
					cleanup(newClauses);
					break;
				}
				int lvpCost;
				int remainingCost = baselinePTPCost-currPTPCost;
				int retVal = lp->runMockPTP(newClauses,LVRClusterList[jj],lvpCost,remainingCost);
				//cout<<remainingCost<<" "<<lvpCost<<endl;
				if(retVal == -1)
				{
					clusteringNP = true;
					newClauses.clear();
					break;
				}
				else if(retVal == -2)
				{
					clusterCostExceeded = true;
					newClauses.clear();
					break;
				}

				currPTPCost += lvpCost;
				currSpaceCost += LVRClusterList[jj]->totalSpace();
				newClauses.clear();
		}
		//if(clusteringNP)
			//cout<<"Clustering Process::Clustering was infeasible, trying next clustering in level "<<level<<endl;
		//else if(clusterCostExceeded || currSpaceCost > baselineSpaceCost)
			//cout<<"Clustering Process::Cost exceeds baseline, trying next clustering in level "<<level<<endl;
		if(!clusteringNP && !clusterCostExceeded && currSpaceCost <= baselineSpaceCost)
		{
			level++;
			int id1 =clusterCombinations[ind]->clusterIndexes.at(0);
			int id2 = clusterCombinations[ind]->clusterIndexes.at(1);
			ofstream outfile;
			outfile.open(params->outClusterFile.c_str(),ios::app);
			outfile<<level<<" "<<id1<<" "<<id2<<endl;
			outfile.close();
			cout<<"Clustering Process::merged ";
			LVRClusterList_orig[id1]->print();
			LVRClusterList_orig[id2]->print();
		}
		else
		{
			//revert back to previous state
			mln->max_predicate_id = originalNumPredicates;
			cleanup(LVRClusterList);
			LVRClusterList.resize(LVRClusterList_orig.size());
				
			for(unsigned int i=0;i<LVRClusterList_orig.size();i++)
			{
				LVRClusterList[i] = new LVRCluster(LVRClusterList_orig[i]->clusterId,*mln);
				for(unsigned int j=0;j<LVRClusterList_orig[i]->elements.size();j++)
					LVRClusterList[i]->addElement(LVRClusterList_orig[i]->elements[j]);
			}
			for(unsigned int i=0;i<LVRClusterList.size();i++)
				LVRClusterList[i]->initializeAssignments();
			computeClusterMB();
		}
		cleanup(clusterCombinations);
		cleanup(LVRClusterList_orig);
	}
	cout<<"Clustering Process:: completed random clustering, stable level = "<<level<<endl;
}

void LClusterCreator::startClustering(LvrParams* params)
{
	int level = 0;
	initializeClusters();
	LogDouble PTPCost;
	int spaceCost = 0;
	int CF = 0;
	LNormalizer ln(*mln);
	int baselinePTPCost=0;
	int baselineSpaceCost = 0;
	int currPTPCost = 0;
	int currSpaceCost = 0;
	int originalNumPredicates = mln->getMaxPredicateId();
	//cout<<"Sampling from cluster level "<<level<<endl;
	if(params->baselineCostMultiplicativeFactor < 1)
		params->baselineCostMultiplicativeFactor = 1.25;

	cout<<"Clustering Process::Setting baseline.. level "<<level<<endl;
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		vector<WClause*> newClauses;
		addPTPEvidence(i,newClauses);
		int lvpCost;
		if(lp->runMockPTP(newClauses,LVRClusterList[i],lvpCost,RAND_MAX) < 0)
		{
			cout<<"Cluster Creator::Non Poly node, terminate current clustering"<<endl;
			return;
		}
		currPTPCost += lvpCost;
		currSpaceCost += LVRClusterList[i]->totalSpace();
		newClauses.clear();
		//LVRClusterList[0]->printAssignment();
	}
	//cout<<currPTPCost<<endl;
	baselinePTPCost = (double)currPTPCost*params->baselineCostMultiplicativeFactor;
	baselineSpaceCost = currSpaceCost;
	//for(unsigned int i=0;i<LVRClusterList.size();i++)
		//LVRClusterList[i]->printAssignment();
	bool changed = false;
	mln->max_predicate_id = originalNumPredicates;
	while(1)
	{
		//check if parent process has exited, then exit
		//if parent id is 1 then parent has exited
#ifndef _MSC_VER 
		if(getppid()==1)
		{
			cout<<"Clustering Process:: Parent has exited, hence terminating";
			exit(0);
		}
		//IF RUNNING in parallel mode sleep during each iteration, give priority to sampler
		if(params->operationMode == EPARALLEL)
			sleep(SLEEPTIME);
#endif
		level++;
		//cout<<"Sampling from cluster level "<<level<<endl;
		cout<<"Cluster Creator::Clustering level "<<level<<endl;
		mln->max_predicate_id = originalNumPredicates;
		vector<LTuple*> clusterCombinations;
		makeClusterPairs(clusterCombinations);
		changed = false;
		if(clusterCombinations.size() == 0)
			break;
		sort(clusterCombinations.begin(),clusterCombinations.end(),sortFunc);
		for(unsigned int i=0;i<clusterCombinations.size();i++)
		{
			vector<LVRCluster*> LVRClusterList_orig(LVRClusterList.size());
			for(unsigned int jj=0;jj<LVRClusterList.size();jj++)
			{
				LVRClusterList_orig[jj] = LVRClusterList[jj]->cloneClusterElements();
			}
			mergeCluster(clusterCombinations[i]->clusterIndexes.at(0),clusterCombinations[i]->clusterIndexes.at(1));
			computeClusterMB();
			LogDouble totalPTPCost;
			int totalSpace = 0;
			int couplingFactor = 0;
			bool clusteringNP = false;
			bool clusterCostExceeded = false;
			currPTPCost = 0;
			currSpaceCost = 0;
			for(unsigned int jj=0;jj<LVRClusterList.size();jj++)
			{
				vector<WClause*> newClauses;
				addPTPEvidence(jj,newClauses);
				if(currPTPCost > baselinePTPCost)
				{
					clusterCostExceeded = true;
					cleanup(newClauses);
					break;
				}
				int lvpCost;
				int remainingCost = baselinePTPCost-currPTPCost;
				int retVal = lp->runMockPTP(newClauses,LVRClusterList[jj],lvpCost,remainingCost);
				if(retVal == -1)
				{
					clusteringNP = true;
					newClauses.clear();
					break;
				}
				else if(retVal == -2)
				{
					clusterCostExceeded = true;
					newClauses.clear();
					break;
				}
				currPTPCost += lvpCost;
				currSpaceCost += LVRClusterList[jj]->totalSpace();
				newClauses.clear();
			}
			//if(clusteringNP)
				//cout<<"Clustering Process::Clustering was infeasible, trying next clustering in level "<<level<<endl;
			//else if(clusterCostExceeded || currSpaceCost > baselineSpaceCost)
				//cout<<"Clustering Process::Cost exceeds baseline, trying next clustering in level "<<level<<endl;
			if(!clusteringNP && !clusterCostExceeded && currSpaceCost <= baselineSpaceCost)
			{
				changed = true;
				int id1 =clusterCombinations[i]->clusterIndexes.at(0);
				int id2 = clusterCombinations[i]->clusterIndexes.at(1);
				ofstream outfile;
				outfile.open(params->outClusterFile.c_str(),ios::app);
				outfile<<level<<" "<<id1<<" "<<id2<<endl;
				outfile.close();
				cout<<"Cluster Creator::merged ";
				LVRClusterList_orig[id1]->print();
				LVRClusterList_orig[id2]->print();
				cleanup(LVRClusterList_orig);
				break;
			}
			else
			{
				//revert back to previous state
				mln->max_predicate_id = originalNumPredicates;
				cleanup(LVRClusterList);
				LVRClusterList.resize(LVRClusterList_orig.size());
				
				for(unsigned int i=0;i<LVRClusterList_orig.size();i++)
				{
					LVRClusterList[i] = new LVRCluster(LVRClusterList_orig[i]->clusterId,*mln);
					for(unsigned int j=0;j<LVRClusterList_orig[i]->elements.size();j++)
						LVRClusterList[i]->addElement(LVRClusterList_orig[i]->elements[j]);
				}
				for(unsigned int i=0;i<LVRClusterList.size();i++)
					LVRClusterList[i]->initializeAssignments();
				computeClusterMB();
			}
			cleanup(LVRClusterList_orig);
		}
		if(!changed)
		{
			//Level x+1 clustering did not work
			//stay at level x
			cout<<"Clustering Process::Cluster level "<<level<<" is infeasible, reverting to level"<<level-1<<endl;
			//cout<<"Clustering complete. Stable clustering reached. Sampling from stable clusters..."<<endl;
			cout<<"Clustering Process::Clustering complete. Stable clustering reached."<<endl;
			break;
		}
		cleanup(clusterCombinations);
	}
}


void LClusterCreator::groundClusterWithMB(int clusterIndex,vector<WClause*>& newClauses)
{
	LVRClusterList[clusterIndex]->groundAllSharedTermsInCNF(mln->clauses,newClauses);
	//ground everything in its MB
	//for(unsigned int i=0;i<newClauses.size();i++)
		//newClauses[i]->print();
	for(unsigned int i=0;i<LVRClusterList[clusterIndex]->MB.size();i++)
	{
		LVRClusterList[clusterIndex]->MB[i]->groundSharedTermsInCNF(newClauses);
	}
}

void LClusterCreator::addPTPEvidence(int clusterIndex,vector<WClause*>& newClauses,bool storeEvidence)
{
	if(LVRClusterList[clusterIndex]->doExactInference())
	{
		for(unsigned int i=0;i<mln->clauses.size();i++)
		{
			if(LVRClusterList[clusterIndex]->isInClause(mln->clauses[i]))
			{
				newClauses.push_back(LvrMLN::create_new_clause(mln->clauses[i]));
			}
		}
	}
	else
	{
		groundClusterWithMB(clusterIndex,newClauses);
	}

	LVRClusterList[clusterIndex]->evidenceReductionCost = 0;
	for(unsigned int i=0;i<LVRClusterList[clusterIndex]->MB.size();i++)
	{
		LVRClusterList[clusterIndex]->MB[i]->reduceCNFByEvidence(newClauses);
		LVRClusterList[clusterIndex]->evidenceReductionCost += LVRClusterList[clusterIndex]->MB[i]->evidenceReductionCost;
	}
	LVRClusterList[clusterIndex]->renamePredicates(newClauses);	
	for(unsigned int i=0;i<newClauses.size();i++)
		newClauses[i]->satisfied=false;
}

void LClusterCreator::initializeClusters()
{
	cleanup(LVRClusterList);
	int clusterId = 0;
	vector<Atom*> completedAtoms;
	//vector<bool> completed(mln->getMaxPredicateId());
	set<int> completed;

	for(unsigned int i=0;i<mln->clauses.size();i++)
	{
		if(mln->clauses[i]->satisfied)
			continue;
		for(unsigned int j=0;j<mln->clauses[i]->atoms.size();j++)
		{
			if(completed.count(mln->clauses[i]->atoms[j]->symbol->id)>0)
			{
				continue;
			}
			LVRCluster* c =  new LVRCluster(clusterId++,*mln);
			c->addElement(mln->clauses[i]->atoms[j]);
			LVRClusterList.push_back(c);
			completed.insert(mln->clauses[i]->atoms[j]->symbol->id);
		}
	}
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		LVRClusterList[i]->initializeAssignments();
	}
	//for(unsigned int i=0;i<LVRClusterList.size();i++)
		//LVRClusterList[i]->printAssignment();
	computeClusterMB();
}

