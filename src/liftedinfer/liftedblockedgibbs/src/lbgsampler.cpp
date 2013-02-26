#include <time.h>
#include "lbgsampler.h"
#include "normalizer.h"
#include <fstream>
#include <sstream>
using namespace std;
#include "queryupdater.h"
#include "stringconversionutils.h"


int LBGSampler::getClusterIndex(LVRCluster* cluster)
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

void LBGSampler::mergeCluster(int clusterInd1,int clusterInd2)
{
	for(unsigned int i=0;i<LVRClusterList[clusterInd2]->elements.size();i++)
		LVRClusterList[clusterInd1]->addElement(LVRClusterList[clusterInd2]->elements[i]);
	LVRClusterList[clusterInd1]->initializeAssignments();
	removeItem(LVRClusterList,clusterInd2);
}

void LBGSampler::getAtomsMatchingTerm(Atom* atom, int termIndex, vector<Atom*>& matchingMLNAtoms)
{
	//vector<bool> completedIds(mln.getMaxPredicateId());
	set<int> completedIds;
	for(unsigned int i=0;i<mln.clauses.size();i++)
	{
		LvrTerm* termToMatch = NULL;
		int atomInd = -1;
		for(unsigned int j=0;j<mln.clauses[i]->atoms.size();j++)
		{
			if(completedIds.count(mln.clauses[i]->atoms[j]->symbol->id)>0)
				continue;
			if(mln.clauses[i]->atoms[j]->symbol->id == atom->symbol->id)
			{
				//capture the term
				termToMatch = mln.clauses[i]->atoms[j]->terms[termIndex];
				atomInd = j;
				break;
			}
		}
		if(termToMatch == NULL)
			continue;
		for(unsigned int j=0;j<mln.clauses[i]->atoms.size();j++)
		{
			if(j==atomInd)
				continue;
			for(unsigned int k=0;k<mln.clauses[i]->atoms[j]->terms.size();k++)
			{
				if(mln.clauses[i]->atoms[j]->terms[k] == termToMatch &&
					completedIds.count(mln.clauses[i]->atoms[j]->symbol->id)==0)
				{
					completedIds.insert(mln.clauses[i]->atoms[j]->symbol->id);
					matchingMLNAtoms.push_back(mln.clauses[i]->atoms[j]);
				}
			}
		}
	}
}

void LBGSampler::computeClusterLiftedSharings()
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
}

void LBGSampler::computeClusterMB()
{
	vector<vector<int> > clauseIndexes(LVRClusterList.size());
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		vector<int> clauseIndex;
		for(unsigned int j=0;j<mln.clauses.size();j++)
		{
			for(unsigned int k=0;k<mln.clauses[j]->atoms.size();k++)
			{
				int ind = LVRClusterList[i]->inCluster(mln.clauses[j]->atoms[k]);
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
	computeClusterLiftedSharings();
}

void LBGSampler::startLVBGibbs(LvrParams* params)
{
	int iterationsInCurrCall = 0;
	//if sampling mode no need to poll, make window size very large
	if(params->operationMode == ESAMPLER)
		windowSize = RAND_MAX;
	else
		windowSize = CLUSTERINGWINDOWSIZE;
	bool printResults = false;
	time_t start,end;
	time (&start);
	//set defaults
	int printInterval = PRINTRESULTSINTERVAL;
	if(params->maxSteps <= 0)
	{
		if(params->isWeightLearning)
			params->maxSteps = MAXSTEPSWL;
		else
			params->maxSteps = MAXSTEPSINFER;
	}
	if(params->maxSeconds <= 0)
	{
		if(params->isWeightLearning)
			params->maxSeconds = MAXSECONDSWL;
		else
			params->maxSeconds = MAXSECONDSINFER;
	}
	if(params->burnMaxSteps <= 0)
		params->burnMaxSteps = BURNINSTEPS;
	bool clustersRead = false;
	while(1)
	{
		if(clusterLevel < 0)
		{
			//initialize to cluster level 0
			initializeClusters();
			clusterLevel = 0;
			originalNumPredicates = mln.getMaxPredicateId();
			iterationsSinceLastMerge = 0;
		}
		//decide on the cluster level
		//read the levels, clusterings from cluster file
		fstream filestr(params->inClusterFile.c_str());
		vector<vector<int> > mergeIndexList;
		vector<int> levels;
		int maxLevel = 0;
		if(filestr != NULL)
		{
			cout<<"Reading current cluster file = "<<params->inClusterFile.c_str()<<endl;
			clustersRead = true;
			char* buf = new char[1024];
			while(filestr)
			{
				filestr.getline(buf,1024);
				string line(buf);
				if(line.size()==0)
					continue;
				string s(line);
				vector<string> vals;
				LStringConversionUtils::tokenize(s,vals," ");
				stringstream st(vals[0]);
				int level;
				st >> level;
				levels.push_back(level);
				if(level > maxLevel)
				{
					maxLevel = level;
				}
				vector<int> mergeIndexes;
				for(unsigned int jj=1;jj<vals.size();jj++)
				{
					int lv;
					stringstream ss(vals[jj]);
					ss >> lv;
					mergeIndexes.push_back(lv);
				}
				mergeIndexList.push_back(mergeIndexes);
			}
			filestr.close();
			delete[] buf;
		}
		if(maxLevel!=clusterLevel)
		{
			//level has changed
			//do all the merges after the current level in order
			cout<<"Clustering has changed, restructuring cluster graph..."<<endl;
			for(unsigned int i=0;i<levels.size();i++)
			{
				if(levels[i] > clusterLevel)
				{
					mergeCluster(mergeIndexList[i].at(0),mergeIndexList[i].at(1));
					clusterLevel = levels[i];
				}
			}
			computeClusterMB();
			iterationsSinceLastMerge = 0;
			mln.max_predicate_id = originalNumPredicates;
		}

		for(unsigned int t=0;t<windowSize;t++)
		{
			time_t t1;
			int currIterationSize;
			if(params->isWeightLearning)
				currIterationSize = iterationsAcrossCalls;
			else
				currIterationSize = iterationsSinceLastMerge;
			if(currIterationSize + 1 > params->burnMaxSteps)
			{
				//burn in is complete
				if(lp->getBurnIn())
				{
					cout<<"Burn In complete."<<endl;
					lp->setBurnIn(false);
				}
			}
			else
			{
				if(!lp->getBurnIn())
				{
					cout<<"Burning samples..."<<endl;
					lp->setBurnIn(true);
				}
				//cout<<"LBG Sampling Process::Burning in level "<<clusterLevel<<"...Remaining to burn ="<<params->burnMaxSteps - currIterationSize<<endl;
			}
			for(unsigned int i=0;i<LVRClusterList.size();i++)
			{
				vector<WClause*> newClauses;
				addPTPEvidence(i,newClauses,false);
				int lvpCost;
				if(lp->runPTP(newClauses,LVRClusterList[i]) == -1)
				{
					cout<<"Lifted Gibbs Sampler::Error,trying to reinitialize!!"<<endl;
					windowSize = RAND_MAX;
					initializeClusters();
					break;
				}
				newClauses.clear();
			}

			iterationsSinceLastMerge++;
			iterationsAcrossCalls++;
			if(!lp->getBurnIn())
			{
				iterationsInCurrCall++;
				//cout<<"LBG Sampling Process::Sampling from cluster level "<<clusterLevel<<endl;
				numTimesSampled++;
				LvrQueryUpdater::Instance()->updateDontCare();
				if(!params->isWeightLearning)
				{
					if(numTimesSampled%printInterval == 0)
					{
						LvrQueryUpdater::Instance()->writeToFile(numTimesSampled);
						cout<<"LBG Sampling Process::Sampling from cluster level "<<clusterLevel<<",iter="<<numTimesSampled<<endl;
					}
				}
				
			}
			time (&end);
			double dif = difftime (end,start);
			if(iterationsInCurrCall > params->maxSteps
					|| dif > params->maxSeconds)
			{
				if(!params->isWeightLearning)
				{
					LvrQueryUpdater::Instance()->writeToFile(numTimesSampled);
					cout<<"Sampling Process exiting"<<endl;
				}
				return;
			}
		}
	}
}

void LBGSampler::groundClusterWithMB(int clusterIndex,vector<WClause*>& newClauses)
{
	LVRClusterList[clusterIndex]->groundAllSharedTermsInCNF(mln.clauses,newClauses);
	//ground everything in its MB
	for(unsigned int i=0;i<LVRClusterList[clusterIndex]->MB.size();i++)
	{
		LVRClusterList[clusterIndex]->MB[i]->groundSharedTermsInCNF(newClauses);
	}
}

void LBGSampler::addPTPEvidence(int clusterIndex,vector<WClause*>& newClauses,bool storeEvidence)
{
	//ground this cluster
	if(LVRClusterList[clusterIndex]->doExactInference())
	{
		for(unsigned int i=0;i<mln.clauses.size();i++)
		{
			if(LVRClusterList[clusterIndex]->isInClause(mln.clauses[i]))
			{
				newClauses.push_back(LvrMLN::create_new_clause(mln.clauses[i]));
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
}

void LBGSampler::initializeClusters()
{
	cout<<"Initializing the cluster graph..."<<endl;
	cleanup(LVRClusterList);
	int clusterId = 0;
	vector<Atom*> completedAtoms;
	set<int> completed;
	for(unsigned int i=0;i<mln.clauses.size();i++)
	{
		for(unsigned int j=0;j<mln.clauses[i]->atoms.size();j++)
		{
			if(completed.count(mln.clauses[i]->atoms[j]->symbol->id) > 0)
			{
				continue;
			}
			LVRCluster* c =  new LVRCluster(clusterId++,mln);
			c->addElement(mln.clauses[i]->atoms[j]);
			LVRClusterList.push_back(c);
			completed.insert(mln.clauses[i]->atoms[j]->symbol->id);
		}
	}
	for(unsigned int i=0;i<LVRClusterList.size();i++)
	{
		LVRClusterList[i]->initializeAssignments();
	}
	computeClusterMB();
}

