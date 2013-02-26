#include "ptpsearch.h"
#include "cleanuputils.h"
#include "propositionalresolution.h"
#include "rulesutil.h"
#include "queryupdater.h"
#include "lvrsingletonnormpropagation.h"


//1-satisfied, 0-false, -1 - neither
int LPTPSearch::CheckCNFSatisfied(vector<WClause*>& CNF, double& satPercent)
{
	int satisfiedcount = 0;
	for(int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->satisfied)
			satisfiedcount++;
		else if(CNF[i]->atoms.size()==0)
		{
			return 0;
		}
	}
	if(satisfiedcount == CNF.size())
	{
		satPercent = 1;
		return 1;
	}
	else
	{
		satPercent = satisfiedcount/(double)CNF.size();
		return -1;
	}
}


LogDouble LPTPSearch::CNFWeight(vector<WClause*>& CNF,set<int>* dontCareQueries)
{
	//compute the total size of the domains of atoms in CNF
	set<int> completedAtomIds;
	map<int,LogDouble> weights;
	map<int,int> numGroundings;
	for(int i=0;i<mln.symbols.size();i++)
	{
		numGroundings[mln.symbols[i]->parentId]=0;
		weights[mln.symbols[i]->parentId] = mln.symbols[i]->pweight+mln.symbols[i]->nweight;
	}
	for(int i=0;i<CNF.size();i++)
	{
		for(int j=0;j<CNF[i]->atoms.size();j++)
		{
			if(completedAtomIds.count(CNF[i]->atoms[j]->symbol->id) == 0)
			{
				completedAtomIds.insert(CNF[i]->atoms[j]->symbol->id);
				if(dontCareQueries)
				{
					//check the query table and add hashes of groundings present in the query table
					lvrptpSearchTree->getHashValues(CNF[i]->atoms[j],(*dontCareQueries));
				}
			}
			else
			{
				//identical atom has already been considered do not recount atom
				continue;
			}
			int totalSize=1;
			for(int k=0;k<CNF[i]->atoms[j]->terms.size();k++)
			{
				totalSize *= CNF[i]->atoms[j]->terms[k]->domain.size();
			}
			map<int,int>::iterator it = numGroundings.find(CNF[i]->atoms[j]->symbol->parentId);
			if(it!=numGroundings.end())
			{
				it->second += totalSize;
			}
		}
	}
	LogDouble totalWeight=LogDouble(1,false);
	for(map<int,int>::iterator it = numGroundings.begin();it!=numGroundings.end();it++)
	{
		if(it->second != 0)
		{
			map<int,LogDouble>::iterator it1 = weights.find(it->first);
			LogDouble atomWt=LogDouble(1,false);
			LogDouble::LDPower(it1->second,it->second,atomWt);
			totalWeight *= atomWt;
		}
	}

	if(totalWeight.is_zero)
		return LogDouble(1,false);
	else
		return totalWeight;
}

//1-satisfied, 0-false, -1 - neither
int LPTPSearch::CheckCNFSatisfiedWithQueries(vector<WClause*>& CNF)
{
	int satisfiedcount = 0;
	int queryCount = 0;
	for(int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->satisfied)
			satisfiedcount++;
		else if(CNF[i]->atoms.size()==0)
		{
			return 0;
		}
		//check if the CNF contains query atoms
		if(LvrQueryUpdater::isInstanceCreated())
		{
			for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
			{
				if(LvrQueryUpdater::Instance()->isNormIdInQuery(CNF[i]->atoms[j]->symbol->normParentId))
					queryCount++;
			}
		}
	}
	if(satisfiedcount == CNF.size())
	{
		//check if any queries were present
		if(queryCount!=0)
		{
			//should continue splitting on the query terms
			return -1;
		}
		return 1;
	}
	else
	{
		return -1;
	}
}

bool LPTPSearch::decomposeCNF(vector<WClause*>& CNF,int& powerFactor)
{
	vector<Decomposer*> decomposer_list;
	decomposer->find_decomposer(CNF,decomposer_list);
	
#ifdef __DEBUG_PRINT__
	cout<<"Decomposer={ ";
	for(int i=0;i<decomposer_list.size();i++)
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
		for(int i=0;i<decomposer_list.size();i++)
		{
			powerFactor*=decomposer_list[i]->decomposer_terms[0]->domain.size();
			for(int j=0;j<decomposer_list[i]->decomposer_terms.size();j++)
			{
				//store pre-decomp domain
				decomposer_list[i]->decomposer_terms[j]->origDomain.clear();
				for(unsigned int jj=0;jj<decomposer_list[i]->decomposer_terms[j]->domain.size();jj++)
					decomposer_list[i]->decomposer_terms[j]->origDomain.push_back(decomposer_list[i]->decomposer_terms[j]->domain[jj]);
				//Erase all but one grounding
				decomposer_list[i]->decomposer_terms[j]->domain.erase(decomposer_list[i]->decomposer_terms[j]->domain.begin(),decomposer_list[i]->decomposer_terms[j]->domain.begin()+decomposer_list[i]->decomposer_terms[j]->domain.size()-1);

			}
		}
	}
	cleanup(decomposer_list);
	return true;
}

LogDouble LPTPSearch::weightedModelCount(vector<WClause*>& CNF,LvrPTPNode* parent,int parentBranchNo)
{
	double satPercent1 = 0;
	int cnfStatus=CheckCNFSatisfied(CNF,satPercent1);	
	if( cnfStatus==0)
	{
		if(LvrQueryUpdater::isInstanceCreated())
			lvrptpSearchTree->createNewNode(ELEAF,parent,parentBranchNo);
		cleanup(CNF);
		return LogDouble(0,false);
	}
	vector<int> powerFactors;
	vector<bool> isDecomposedCNF;
	vector<vector<WClause*> > decomposedList;
	LClusterUtil::seperateDisjointCNF(CNF,decomposedList);
	//CNF can be cleaned up since split into clusters
	cleanup(CNF);
	LogDouble totalVal(1,false);
	LvrPTPNode* djnode = NULL;
	if(LvrQueryUpdater::isInstanceCreated())
		djnode = lvrptpSearchTree->createNewNode(EAND,parent,parentBranchNo);

	for(int t=0;t<decomposedList.size();t++)
	{
		int dum=1;
		bool isDecomposed=false;
		int powerFactor;
		isDecomposed = decomposeCNF(decomposedList[t],powerFactor);
		if(isDecomposed)
		{
#ifdef __DEBUG_PRINT__		
			LvrMLN::print(decomposedList[t],"DECOMPOSED");
#endif
			//use power rule
			LvrPTPNode* decompNode = NULL;
			if(LvrQueryUpdater::isInstanceCreated())
			{
				decompNode = lvrptpSearchTree->createNewNode(EPOWER,djnode,t);
				decompNode->powerFactor = powerFactor;
			}
			LogDouble mcnt = weightedModelCount(decomposedList[t],decompNode,0);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			if(LvrQueryUpdater::isInstanceCreated())
				decompNode->nodeValue = val;
			continue;
		}

		else
		{
			//normalizer->normalizeClauses(decomposedList[t]);
			//do unit propagation
#ifdef __DEBUG_PRINT1__
	LvrMLN::print(decomposedList[t],"NORMALIZED");
#endif
			if(dounitpropagate)
			{
				LogDouble unitPropVal;
				if(LvrQueryUpdater::isInstanceCreated())
					unitPropVal = extensions->unitPropagation(decomposedList[t],true);
				else
					unitPropVal = extensions->unitPropagation(decomposedList[t]);
				if(unitPropVal.is_zero)
				{
					cleanupItems(decomposedList,t);
	#ifdef __DEBUG_PRINT__
				cout<<"***********UnitPROP***********"<<endl;
				cout<<"{}"<<endl;
				cout<<"************UNITProp***************"<<endl;
	#endif
					if(LvrQueryUpdater::isInstanceCreated())
						lvrptpSearchTree->createNewNode(ELEAF,djnode,t);
					return unitPropVal;
				}
				else
				{
					totalVal = totalVal*unitPropVal;
				}
			}
	
#ifdef __DEBUG_PRINT__
LvrMLN::print(decomposedList[t],"UNITPROP");	
#endif
			double satPercent = 0;
			int status=CheckCNFSatisfied(decomposedList[t],satPercent);
			//int status=CheckCNFSatisfiedWithQueries(decomposedList[t]);
			LogDouble wt = LogDouble();
			if(status == 1)
			{
				if(LvrQueryUpdater::isInstanceCreated())
				{
					set<int>* dontcarequeries = new set<int>();
					//if the leaf contains any query atoms (they become dont care), store the hash representations of these
					LogDouble wt = CNFWeight(decomposedList[t],dontcarequeries);
					LvrPTPNode* leaf = lvrptpSearchTree->createNewNode(ELEAF,djnode,t);
					leaf->dontCareQueries = dontcarequeries;
					leaf->nodeValue = wt;
					totalVal = totalVal * wt;
				}
				else
				{
					totalVal = totalVal * CNFWeight(decomposedList[t]);
				}
				cleanup(decomposedList[t]);
				continue;
			}

			
			vector<vector<WClause*> > tempCNFList;
			Atom* atom = heuristics->getAtomToSplit(decomposedList[t]);
#ifdef __DEBUG_PRINT__
				cout<<"Splitting on Atom ";
				atom->print();
				cout<<endl;
#endif
			bool cached = false;
			if(satPercent1 < 0.75)
			{
				cached = extensions->getFromCache(decomposedList[t],wt);
			}
			if(!cached)
			{
				//check for self joins
				bool selfJoined = false;
				for(int m=0;m<decomposedList[t].size();m++)
				{
					if(decomposedList[t].at(m)->isSelfJoinedOnAtom(atom))
					{
						selfJoined = true;
						break;
					}
				}
				vector<LogDouble> branchWeights;
				if(atom->isConstant())
				{
					int idToUse = id++;
					//propositional atom
					if(!selfJoined)
						PropositionalResolution::doPropositionalSplit(decomposedList[t],atom,tempCNFList);
					else
						PropositionalResolution::doPropositionalSplitSelfJoins(decomposedList[t],atom,tempCNFList);
					LogDouble symPosWt = atom->symbol->pweight;
					LogDouble symNegWt = atom->symbol->nweight;
					LvrPTPNode* splitNode = NULL;
					if(LvrQueryUpdater::isInstanceCreated())
					{
						splitNode = lvrptpSearchTree->createNewNode(ESPLITTER,djnode,t);
						splitNode->atom = LvrMLN::create_new_atom(atom);
					}
					LogDouble w1 = symNegWt*weightedModelCount(tempCNFList[0],splitNode,0);
					LogDouble w2 = symPosWt*weightedModelCount(tempCNFList[1],splitNode,1);
					//cleanup(decomposedList[t]);
					wt += w1 + w2;
					if(LvrQueryUpdater::isInstanceCreated())
					{
						splitNode->branchWeights.push_back(symNegWt);
						splitNode->branchWeights.push_back(symPosWt);
						splitNode->nodeValue = w1+w2;
					}
					cleanup(tempCNFList[0]);
					cleanup(tempCNFList[1]);
					tempCNFList.clear();
				}
				else //first order atom
				{
					//check in cache
					int singletonIndex;
					bool singleton = atom->isSingletonAtom(singletonIndex);
					bool nonpoly = false;
					if(!singleton)
						nonpoly = true;
					else if(selfJoined)
					{
						//check for blocking rule
						bool blocked = LRulesUtil::isBlocked(atom,singletonIndex,decomposedList[t]);
						if(!blocked)
							nonpoly=true;
					}
					//lifted domain reduction
					if(!nonpoly)
					{
						/*
						splitter->doLiftedSplit(decomposedList[t],atom,tempCNFList,singletonIndex);
						LogDouble symPosWt = atom->symbol->pweight;
						LogDouble symNegWt = atom->symbol->nweight;
						int domSize = atom->terms[singletonIndex]->domain.size();
						//int i = tempCNFList.size()-1;
						LvrPTPNode* splitNode = NULL;
						if(LvrQueryUpdater::isInstanceCreated())
						{
							splitNode = lvrptpSearchTree->createNewNode(ESPLITTER,djnode,t);
							splitNode->atom = LvrMLN::create_new_atom(atom);
							branchWeights.clear();
							branchWeights.resize(tempCNFList.size());
						}
						LogDouble norm;
						for(unsigned int jj=0;jj<tempCNFList.size();jj++)
						{
							LogDouble pos = LogDouble(1,false);
							LogDouble neg = LogDouble(1,false);
							LogDouble::LDPower(symPosWt,jj,pos);
							LogDouble::LDPower(symNegWt,domSize-jj,neg);
							LogDouble binCoeff = LogDouble::Binomial(domSize,jj);
							normalizer->normalizeClauses(tempCNFList[jj],true,false);
							LogDouble mcnt = weightedModelCount(tempCNFList[jj],splitNode,jj);
							wt+=binCoeff*pos*neg*mcnt;
							if(LvrQueryUpdater::isInstanceCreated())
								branchWeights[jj] = pos*neg;
							norm += binCoeff*pos*neg*mcnt; 
						}
						while(!tempCNFList.empty())
						{
							cleanup(tempCNFList.back());
							tempCNFList.pop_back();
						}
						*/
						Atom* tmpAtom = LvrMLN::create_new_atom(atom);
						LogDouble norm;
						LogDouble symPosWt = atom->symbol->pweight;
						LogDouble symNegWt = atom->symbol->nweight;
						int domSize = atom->terms[singletonIndex]->domain.size();

						LvrPTPNode* splitNode = NULL;
						if(LvrQueryUpdater::isInstanceCreated())
						{
							splitNode = lvrptpSearchTree->createNewNode(ESPLITTER,djnode,t);
							splitNode->atom = LvrMLN::create_new_atom(atom);
							branchWeights.clear();
							branchWeights.resize(domSize+1);
						}

						for(int jj=0;jj<=domSize;jj++)
						{
							LogDouble pos = LogDouble(1,false);
							LogDouble neg = LogDouble(1,false);
							LogDouble::LDPower(symPosWt,jj,pos);
							LogDouble::LDPower(symNegWt,domSize-jj,neg);
							LogDouble binCoeff = LogDouble::Binomial(domSize,jj);
							vector<WClause*> lvClausesCopy;
							LvrMLN::copyAllClauses(decomposedList[t],lvClausesCopy);
							LvrSingletonNormPropagation::propagateNormalizedCNF(lvClausesCopy,tmpAtom,singletonIndex,jj);
							LogDouble mcnt = weightedModelCount(lvClausesCopy,splitNode,jj);
							cleanup(lvClausesCopy);
							wt+=binCoeff*pos*neg*mcnt;
							if(LvrQueryUpdater::isInstanceCreated())
								branchWeights[jj] = pos*neg;
							norm += binCoeff*pos*neg*mcnt;
						}
						delete tmpAtom;
						
						if(LvrQueryUpdater::isInstanceCreated())
						{
							splitNode->nodeValue = norm;
							splitNode->branchWeights = branchWeights;
						}
					}
					else //NON-Singleton or selfjoined atom,non polynomial grounding
					{
						
						cout<<"Non polynomial splitting atom encountered...may take a long time to compute"<<endl;
						vector<WClause*> groundCNF;
						vector<vector<int> > allGroundings;
						int numGroundings;
						splitter->doFullGrounding(decomposedList[t],atom,groundCNF,allGroundings,numGroundings);
						/*
						//Do resolution in chunks to avoid memory allocation problems
						int startnum = 0;
						int chunkSize = SPLITCHUNKSIZE;
						bool done = false;
						LogDouble symPosWt = atom->symbol->pweight;
						LogDouble symNegWt = atom->symbol->nweight;
						int atomId = atom->symbol->id;
						while(!done)
						{
							vector<vector<WClause*> > resolvedCNFList;
							//resolve exactly chunksize truth assignments
							splitter->doCompleteSplitInChunks(groundCNF,atomId, allGroundings, numGroundings,
							resolvedCNFList,startnum, chunkSize,done);
							int endnum = startnum + resolvedCNFList.size()-1;
							while(!resolvedCNFList.empty())
							{
								//get the number of 1's in the truth assignment
								int n = endnum;
								int count=0;
								while(n)
								{
									count++;
									n &= (n-1);
								}

								LogDouble pos = LogDouble(1,false);
								LogDouble neg = LogDouble(1,false);
								LogDouble::LDPower(symPosWt,count,pos);
								LogDouble::LDPower(symNegWt,numGroundings-count,neg);
								normalizer->normalizeClauses(resolvedCNFList.back(),true,false);
								LogDouble mcnt = weightedModelCount(resolvedCNFList.back(),0,endnum);
								cleanup(resolvedCNFList.back());
								resolvedCNFList.pop_back();
								endnum--;
								wt=wt+pos*neg*mcnt;
							}
							//increment the start position for next chunk
							startnum += chunkSize;
						}
						*/
						vector<vector<WClause*> > resolvedCNFList;
						bool done = false;
						LogDouble symPosWt = atom->symbol->pweight;
						LogDouble symNegWt = atom->symbol->nweight;
						int atomId = atom->symbol->id;
						splitter->doCompleteSplitInChunks(groundCNF,atomId, allGroundings, numGroundings,
						resolvedCNFList,0, pow(2.0,numGroundings),done);
						LvrPTPNode* splitNode = NULL;
						if(LvrQueryUpdater::isInstanceCreated())
						{
							splitNode = lvrptpSearchTree->createNewNode(ESPLITTER,djnode,t);
							splitNode->atom = LvrMLN::create_new_atom(atom);
							branchWeights.clear();
							branchWeights.resize(resolvedCNFList.size());
						}
						LogDouble norm;
						for(unsigned int jj=0;jj<resolvedCNFList.size();jj++)
						{
							//get the number of 1's in the truth assignment
							int n = jj;
							int count=0;
							while(n)
							{
								count++;
								n &= (n-1);
							}

							LogDouble pos = LogDouble(1,false);
							LogDouble neg = LogDouble(1,false);
							LogDouble::LDPower(symPosWt,count,pos);
							LogDouble::LDPower(symNegWt,numGroundings-count,neg);
							normalizer->normalizeClauses(resolvedCNFList[jj],true,false);
							LogDouble mcnt = weightedModelCount(resolvedCNFList[jj],splitNode,jj);
							wt=wt+pos*neg*mcnt;
							if(LvrQueryUpdater::isInstanceCreated())
								branchWeights[jj] = pos*neg;
							norm += pos*neg*mcnt; 

						}
						while(!resolvedCNFList.empty())
						{
							cleanup(resolvedCNFList.back());
							resolvedCNFList.pop_back();
						}
						if(LvrQueryUpdater::isInstanceCreated())
						{
							splitNode->nodeValue = norm;
							splitNode->branchWeights = branchWeights;
						}
						cleanup(groundCNF);				
						
					}//non singleton
				}//first order atom
				if(satPercent < 0.75)
					extensions->storeToCache(decomposedList[t],wt);
			}//not cached
#ifdef __DEBUG_PRINT__
		cout<<"Wt="<<wt<<endl;
#endif
			cleanup(decomposedList[t]);
			if(wt.is_zero)
			{
				if(LvrQueryUpdater::isInstanceCreated())
					lvrptpSearchTree->createNewNode(ELEAF,djnode,t);
				cleanupItems(decomposedList,t+1);
				return wt;
			}
			totalVal = totalVal*wt;
		}//end split step
	}//end for
	if(LvrQueryUpdater::isInstanceCreated())
		djnode->nodeValue = totalVal;
	return totalVal;
}

/*
LogDouble LPTPSearch::weightedModelCount(vector<WClause*>& CNF,int parentId,int parentBranchNo)
{
	double satPercent1 = 0;
	int cnfStatus=CheckCNFSatisfied(CNF,satPercent1);	
	if( cnfStatus==0)
	{
		cleanup(CNF);
		return LogDouble(0,false);
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
	LPTPNode* lnode = NULL;
	if(decomposedList.size() > 1)
	{
		lnode = lvrPTPTreeSampling->addTreeNode(EDISJOINT,id,parentId,parentBranchNo);
		id++;
		isParentChanged = true;
		changedParent = id-1;
	}

	for(int t=0;t<decomposedList.size();t++)
	{
		if(isParentChanged)
		{
		parentId = changedParent;
		parentBranchNo = t;
		}

		int powerFactor;
		bool isDecomposed = decomposeCNF(decomposedList[t],powerFactor);
		if(isDecomposed)
		{
#ifdef __DEBUG_PRINT__		
			LvrMLN::print(decomposedList[t],"DECOMPOSED");
#endif
			int idToUse = id++;
			//use power rule
			LogDouble mcnt = weightedModelCount(decomposedList[t],idToUse,0);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			lvrPTPTreeSampling->addTreeNode(EDECOMPOSER,idToUse,parentId,parentBranchNo,powerFactor,totalVal);
			continue;
		}
		else
		{
			//normalizer->normalizeClauses(decomposedList[t]);
			//do unit propagation
#ifdef __DEBUG_PRINT1__
	LvrMLN::print(decomposedList[t],"NORMALIZED");
#endif

			LogDouble unitPropVal = extensions->unitPropagation(decomposedList[t]);
			if(unitPropVal.is_zero)
			{
				cleanupItems(decomposedList,t);
#ifdef __DEBUG_PRINT__
			cout<<"***********UnitPROP***********"<<endl;
			cout<<"{}"<<endl;
			cout<<"************UNITProp***************"<<endl;
#endif
				return unitPropVal;
			}
			else
			{
				totalVal = totalVal*unitPropVal;
			}
		
#ifdef __DEBUG_PRINT__
LvrMLN::print(decomposedList[t],"UNITPROP");	
#endif
			double satPercent = 0;
			int status=CheckCNFSatisfied(decomposedList[t],satPercent);
			LogDouble wt = LogDouble();
			if(status == 1)
			{
				totalVal = totalVal * CNFWeight(decomposedList[t]);
				cleanup(decomposedList[t]);
				continue;
			}

			
			vector<vector<WClause*> > tempCNFList;
			Atom* atom = heuristics->getAtomToSplit(decomposedList[t]);
#ifdef __DEBUG_PRINT__
				cout<<"Splitting on Atom ";
				atom->print();
				cout<<endl;
#endif
			bool cached = false;
			if(satPercent1 < 0.75)
			{
				cached = extensions->getFromCache(decomposedList[t],wt);
			}
			if(!cached)
			{
				//check for self joins
				bool selfJoined = false;
				for(int m=0;m<decomposedList[t].size();m++)
				{
					if(decomposedList[t].at(m)->isSelfJoinedOnAtom(atom))
					{
						selfJoined = true;
						break;
					}
				}
				vector<LogDouble>* branchWeights = new vector<LogDouble>();
				if(atom->isConstant())
				{
					int idToUse = id++;
					//propositional atom
					if(!selfJoined)
						PropositionalResolution::doPropositionalSplit(decomposedList[t],atom,tempCNFList);
					else
						PropositionalResolution::doPropositionalSplitSelfJoins(decomposedList[t],atom,tempCNFList);
					LogDouble symPosWt = atom->symbol->pweight;
					LogDouble symNegWt = atom->symbol->nweight;
					LogDouble w1 = symNegWt*weightedModelCount(tempCNFList[0],idToUse,0);
					LogDouble w2 = symPosWt*weightedModelCount(tempCNFList[1],idToUse,1);
					//cleanup(decomposedList[t]);
					wt += w1 + w2;
					vector<LogDouble> wts;
					wts.push_back(w1);
					wts.push_back(w2);
					LogDouble norm = w1+w2;
					vector<LogDouble> bincoefs;
					bincoefs.push_back(LogDouble(1,false));
					bincoefs.push_back(LogDouble(1,false));
					(*branchWeights).push_back(symNegWt);
					(*branchWeights).push_back(symPosWt);
					lvrPTPTreeSampling->addTreeNode(atom,wts,bincoefs,norm,ESPLIT,idToUse,parentId,parentBranchNo,branchWeights,wt);
					cleanup(tempCNFList[0]);
					cleanup(tempCNFList[1]);
					tempCNFList.clear();
				}
				else //first order atom
				{
					int idToUse = id++;
					//check in cache
					int singletonIndex;
					bool singleton = atom->isSingletonAtom(singletonIndex);
					bool nonpoly = false;
					if(!singleton)
						nonpoly = true;
					else if(selfJoined)
					{
						//check for blocking rule
						bool blocked = LRulesUtil::isBlocked(atom,singletonIndex,decomposedList[t]);
						if(!blocked)
							nonpoly=true;
					}
					vector<LogDouble> zvals;
					vector<LogDouble> binceoffs;
					LogDouble norm;
					//lifted domain reduction
					if(!nonpoly)
					{
						splitter->doLiftedSplit(decomposedList[t],atom,tempCNFList,singletonIndex);
						LogDouble symPosWt = atom->symbol->pweight;
						LogDouble symNegWt = atom->symbol->nweight;
						int domSize = atom->terms[singletonIndex]->domain.size();
						int i = tempCNFList.size()-1;
						while(!tempCNFList.empty())
						{
							LogDouble pos = LogDouble(1,false);
							LogDouble neg = LogDouble(1,false);
							LogDouble::LDPower(symPosWt,i,pos);
							LogDouble::LDPower(symNegWt,domSize-i,neg);
							LogDouble binCoeff = LogDouble::Binomial(domSize,i);
							normalizer->normalizeClauses(tempCNFList.back(),true,false);
							//for(int jj=0;jj<tempCNFList.back().size();jj++)
								//tempCNFList.back().at(jj)->print();
							LogDouble mcnt = weightedModelCount(tempCNFList.back(),idToUse,i);
							//cout<<binCoeff<<" "<<pos<<" "<<neg<<" "<<mcnt<<endl;
							cleanup(tempCNFList.back());
							tempCNFList.pop_back();
							wt+=binCoeff*pos*neg*mcnt;
							i--;
							zvals.push_back(pos*neg*mcnt);
							binceoffs.push_back(binCoeff);
							norm = norm + binCoeff*mcnt;
							(*branchWeights).push_back(pos*neg);
						}
						lvrPTPTreeSampling->addTreeNode(atom,zvals,binceoffs,norm,ESPLIT,idToUse,parentId,parentBranchNo,branchWeights,wt);
					}
					else //NON-Singleton or selfjoined atom,non polynomial grounding
					{
						cout<<"Non polynomial splitting atom encountered...may take a long time to compute"<<endl;
						vector<WClause*> groundCNF;
						vector<vector<int> > allGroundings;
						int numGroundings;
						splitter->doFullGrounding(decomposedList[t],atom,groundCNF,allGroundings,numGroundings);
						//Do resolution in chunks to avoid memory allocation problems
						int startnum = 0;
						int chunkSize = SPLITCHUNKSIZE;
						bool done = false;
						LogDouble symPosWt = atom->symbol->pweight;
						LogDouble symNegWt = atom->symbol->nweight;
						int atomId = atom->symbol->id;
						while(!done)
						{
							vector<vector<WClause*> > resolvedCNFList;
							//resolve exactly chunksize truth assignments
							splitter->doCompleteSplitInChunks(groundCNF,atomId, allGroundings, numGroundings,
							resolvedCNFList,startnum, chunkSize,done);
							int endnum = startnum + resolvedCNFList.size()-1;
							while(!resolvedCNFList.empty())
							{
								//get the number of 1's in the truth assignment
								int n = endnum;
								int count=0;
								while(n)
								{
									count++;
									n &= (n-1);
								}

								LogDouble pos = LogDouble(1,false);
								LogDouble neg = LogDouble(1,false);
								LogDouble::LDPower(symPosWt,count,pos);
								LogDouble::LDPower(symNegWt,numGroundings-count,neg);
								normalizer->normalizeClauses(resolvedCNFList.back(),true,false);
								LogDouble mcnt = weightedModelCount(resolvedCNFList.back(),idToUse,endnum);
								cleanup(resolvedCNFList.back());
								resolvedCNFList.pop_back();
								endnum--;
								wt=wt+pos*neg*mcnt;
								zvals.push_back(pos*neg*mcnt);
								binceoffs.push_back(LogDouble(1,false));
								norm = norm + mcnt;
								(*branchWeights).push_back(pos*neg);
							}
							//increment the start position for next chunk
							startnum += chunkSize;
							lvrPTPTreeSampling->addTreeNode(atom,zvals,binceoffs,norm,ESPLIT,idToUse,parentId,parentBranchNo,branchWeights,wt);
						}
						cleanup(groundCNF);				
					}//non singleton
				}//first order atom
			if(satPercent < 0.75)
				extensions->storeToCache(decomposedList[t],wt);
			}//not cached
#ifdef __DEBUG_PRINT__
		cout<<"Wt="<<wt<<endl;
#endif
			cleanup(decomposedList[t]);
			if(wt.is_zero)
			{
				cleanupItems(decomposedList,t+1);
				return wt;
			}
			totalVal = totalVal*wt;
		}//end split step
	}//end for
	if(lnode)
		lnode->nodeValue = totalVal;
	return totalVal;
}
*/

LogDouble LPTPSearch::weightedModelCountApprox(vector<WClause*>& CNF)
{
	double satPercent1 = 0;
	int cnfStatus=CheckCNFSatisfied(CNF,satPercent1);	
	if( cnfStatus==0)
	{
		cleanup(CNF);
		return LogDouble(0,false);
	}

	vector<int> powerFactors;
	vector<bool> isDecomposedCNF;
	vector<vector<WClause*> > decomposedList;
	LClusterUtil::seperateDisjointCNF(CNF,decomposedList);
	//CNF can be cleaned up since split into clusters
	cleanup(CNF);
	LogDouble totalVal(1,false);
	for(int t=0;t<decomposedList.size();t++)
	{		
		int powerFactor;
		bool isDecomposed = decomposeCNF(decomposedList[t],powerFactor);
		if(isDecomposed)
		{
#ifdef __DEBUG_PRINT__		
	LvrMLN::print(decomposedList[t],"DECOMPOSED");
#endif

			//use power rule
			LogDouble mcnt = weightedModelCountApprox(decomposedList[t]);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			continue;
		}
		else
		{
			//normalizer->normalizeClauses(decomposedList[t]);
#ifdef __DEBUG_PRINT__
	LvrMLN::print(decomposedList[t],"NORMALIZED");
#endif
			LogDouble unitPropVal = extensions->unitPropagation(decomposedList[t]);
			if(unitPropVal.is_zero)
			{
				cleanupItems(decomposedList,t);
#ifdef __DEBUG_PRINT__
	cout<<"***********UnitPROP***********"<<endl;
	cout<<"{}"<<endl;
	cout<<"************UNITProp***************"<<endl;
#endif
				return unitPropVal;
			}
			else
			{
				totalVal = totalVal*unitPropVal;
			}
#ifdef __DEBUG_PRINT__
	LvrMLN::print(decomposedList[t],"UNITPROP");	
#endif
			double satPercent = 0;
			int status=CheckCNFSatisfied(decomposedList[t],satPercent);
			LogDouble wt = LogDouble();

			if(status == 1)
			{
				totalVal = totalVal * CNFWeight(decomposedList[t]);
				cleanup(decomposedList[t]);
				continue;
			}
			else if( status==0)
			{
				cleanup(decomposedList[t]);
				return LogDouble(0,false);
			}
			
			vector<vector<WClause*> > tempCNFList;
			Atom* atom = heuristics->getAtomToSplit(decomposedList[t]);
#ifdef __DEBUG_PRINT__
	cout<<"Splitting on Atom ";
	atom->print();
	cout<<")"<<endl;
#endif

#ifdef __DEBUG_PRINT__
	LvrMLN::print(decomposedList[t],"SPLIT-INPUT");	
#endif

				
			int completedCount = 0;
			vector<int> reductionFactor;
			int numGroundings = atom->getNumberOfGroundings();
			vector<bool> completedGroups(numGroundings + 1);
			while(1)
			{
				int numTrue;
				LogDouble probOfSample;

				splitter->doApproxSplit(decomposedList[t], atom, numGroundings, completedGroups, completedCount, probOfSample, numTrue, tempCNFList);
				//cout<<"here"<<endl;
				int numFalse = numGroundings - numTrue;
				//just 1 sample collected
				LogDouble pos = LogDouble(1,false);
				LogDouble neg = LogDouble(1,false);
				LogDouble::LDPower(atom->symbol->pweight,numTrue,pos);
				LogDouble::LDPower(atom->symbol->nweight,numFalse,neg);
				LogDouble binCoeff = LogDouble::Binomial(numTrue+numFalse,numTrue);
				normalizer->normalizeClauses(tempCNFList.back(),true,false);
				LogDouble mcnt = weightedModelCountApprox(tempCNFList.back());
				cleanup(tempCNFList.back());
				tempCNFList.pop_back();
				if(mcnt.is_zero)
				{
					//cout<<"backtrack!"<<endl;
					completedGroups[numTrue] = true;
					completedCount++;
					if(completedCount == completedGroups.size())
					{
						wt = 0;
						break;
					}
					else
						continue;
				}
				wt+=binCoeff*pos*neg*mcnt/probOfSample;
				break;
			}
#ifdef __DEBUG_PRINT__
cout<<"Wt="<<wt<<endl;
#endif
			cleanup(decomposedList[t]);
			if(wt.is_zero)
			{
				cleanupItems(decomposedList,t+1);
				return wt;
			}
			totalVal = totalVal*wt;
		}
	}
	return totalVal;
}


LogDouble LPTPSearch::startApproxWeightedModelCounting(LvrParams* params)
{
	//normalizer->normalizeClauses(mln.clauses,true,false);
	cout<<"Starting Lifted Weighted Model Counting..."<<endl;
	time_t start;
	time(&start);
	cout<<"Time ="<<ctime(&start)<<endl;
	LogDouble ZApprox;
	int iterations = 0;
	if(params->maxSteps <= 0)
		params->maxSteps = MAXSTEPSINFER;
	if(params->maxSeconds <= 0)
		params->maxSeconds = MAXSECONDSINFER;

	while(1)
	{
		//make a copy of normalized clauses
		vector<WClause*> clauses;
		LvrMLN::copyAllClauses(mln.clauses,clauses);
		LogDouble currZ = weightedModelCountApprox(clauses)*mln.conversionFactor;
		//update approximation
		ZApprox = ZApprox*LogDouble(iterations,false) + currZ;
		iterations++;
		//take average
		ZApprox = ZApprox / LogDouble(iterations,false);
		//cout<<"iteration="<<iterations<<", currZ="<<currZ<<"currZ (Log Space)="<<currZ.value<<", ZApprox="<<ZApprox<<", ZApprox(Log Space)="<<ZApprox.value<<endl;
		cout<<"iteration="<<iterations<<", currZ : ";
		currZ.printValue();
		cout<<", ZApprox : ";
		ZApprox.printValue();
		cout<<endl;
		if(iterations >= params->maxSteps)
			return ZApprox; 
		time_t currTime;
		time(&currTime);
		if(difftime(currTime,start) > params->maxSeconds)
			break;
	}
	return 0;
}

LogDouble LPTPSearch::startExactWeightedModelCounting(LvrParams* params)
{
	dounitpropagate = true;
	cout<<"Starting Exact Lifted Model Counting"<<endl;
	id = 0;
	normalizer->normalizeClauses(mln.clauses,true,false);
	//make a copy of normalized clauses
	vector<WClause*> clauses;
	LvrMLN::copyAllClauses(mln.clauses,clauses);
	if(!LvrQueryUpdater::isInstanceCreated())
		extensions->setCaching(params->ptpCacheSize);
	else
	{
		//turn off unit propagation
		dounitpropagate = false;
		//try to see if we can decompose, if yes only consider part where queries occur
		vector<WClause*> clausesToConsider;
		vector<vector<WClause*> > tmpClauses;
		LClusterUtil::seperateDisjointCNF(mln.clauses,tmpClauses);
		for(unsigned int i=0;i<tmpClauses.size();i++)
		{
			bool queryInClauses=false;
			for(unsigned int j=0;j<tmpClauses[i].size();j++)
			{
				for(unsigned int k=0;k<tmpClauses[i].at(j)->atoms.size();k++)
				{
					if(LvrQueryUpdater::Instance()->isNormIdInQuery(tmpClauses[i].at(j)->atoms[k]->symbol->normParentId))
					{
						queryInClauses=true;
						break;
					}
				}
				if(queryInClauses)
					break;
			}
			if(queryInClauses)
			{
				for(unsigned int j=0;j<tmpClauses[i].size();j++)
				{
					clausesToConsider.push_back(tmpClauses[i].at(j));
				}
			}
		}
		cleanup(clauses);
		clauses=clausesToConsider;
	}
	//takes ownership of clauses
	LogDouble wt1 = weightedModelCount(clauses);
	//cout<<wt1<<endl;
	//LogDouble wt1 = weightedModelCount(clauses,NULL,-1);
	LogDouble wt = wt1*mln.conversionFactor;
	//lvrPTPTreeSampling->startSamplingMAR();
	if(LvrQueryUpdater::isInstanceCreated())
	{
		lvrptpSearchTree->downPropagationMAR(lvrptpSearchTree->getRoot());
		//lvrptpSearchTree->printTree(lvrptpSearchTree->getRoot());
		//lvrptpSearchTree->bfsprintTree(lvrptpSearchTree->getRoot());
		lvrptpSearchTree->updateQueries(lvrptpSearchTree->getRoot());
		LvrQueryUpdater::Instance()->writeToFile(wt1);
	}
	return wt;
}

LPTPSearch::LPTPSearch(LvrMLN& mln_): mln(mln_)
{
	splitter = new LSplitter();
	decomposer = new LDecomposer(mln);
	extensions = new LExtensions(mln);
	heuristics = new LHeuristics(*decomposer,mln);
	normalizer = new LNormalizer(mln);
	lvrPTPTreeSampling = new LvrPTPTreeSampling();
	lvrptpSearchTree = new LvrPTPSearchTree();
}

