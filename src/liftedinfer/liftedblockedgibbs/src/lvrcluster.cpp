#include "lvrcluster.h"
#include "clusterutil.h"
#include "decomposer.h"
#include "splitter.h"
#include "normalizer.h"
#include "lvrpermutations.h"
#include "hashalgorithm.h"
#include "randomgenutil.h"

LVRCluster& LVRCluster::operator=(const LVRCluster& other)
{
	clusterId = other.clusterId;
	//elements.clear();
	cleanup(elements);
	for(unsigned int i=0;i<other.elements.size();i++)
		addElement(other.elements[i]);
	return *this;
}
LVRCluster* LVRCluster::cloneClusterElements()
{
	LVRCluster* cluster = new LVRCluster(clusterId,mln);
	for(unsigned int i=0;i<elements.size();i++)
	{
		cluster->elements.push_back(LvrMLN::create_new_atom(elements[i]));
	}
	return cluster;
}

LVRCluster::LVRCluster(const LVRCluster& other):clusterId(other.clusterId),mln(other.mln)
{
	for(unsigned int i=0;i<other.elements.size();i++)
		addElement(other.elements[i]);
}
void LVRCluster::removeElements(vector<Atom*> atoms)
{
	for(unsigned int i=0;i<atoms.size();i++)
	{
		int ind = inCluster(atoms[i]);
		if(ind != -1)
			removeItem(elements,ind);
	}

	//cleanup(assignments);
	MB.clear();
	initializeAssignments();
}

int LVRCluster::inCluster(Atom* atom)
{
	for(unsigned int i=0;i<elements.size();i++)
	{
		if(atom->symbol->id == elements[i]->symbol->id)
		//if(atom->symbol->parentId == elements[i]->symbol->id)
			return i;
	}
	return -1;
}

int LVRCluster::getEquivalentElementIndex(Atom* atom)
{
	for(unsigned int i=0;i<elements.size();i++)
	{
		if(atom->symbol->parentId == elements[i]->symbol->id)
		//if(atom->symbol->symbol.compare(elements[i]->symbol->symbol) == 0)
			return i;

	}
	return -1;
}

Atom* LVRCluster::createGroundedAtom(int elementIndex,vector<bool> sharedPos, vector<int> grounding)
{
	Atom* atom = LvrMLN::create_new_atom(elements[elementIndex]);
	int iter=0;
	for(unsigned int i=0;i<sharedPos.size();i++)
	{
		if(sharedPos[i])
		{
			//int val = atom->terms[i]->domain.at(grounding[iter++]);
			atom->terms[i]->domain.clear();
			atom->terms[i]->domain.push_back(grounding[iter++]);
			//atom->terms[i]->domain.push_back(grounding[iter++]);
		}
	}
	return atom;
}

void LVRCluster::initializeAssignments()
{

	//sharedPos.clear();
	//sharedPos.resize(elements.size());
	vector<vector<bool> > sharedPos(elements.size());
	for(unsigned int i=0;i<sharedPos.size();i++)
	{
		sharedPos[i].resize(elements[i]->terms.size());
	}
	for(unsigned int i=0;i<mln.clauses.size();i++)
	{
		vector<vector<LvrTerm*> > terms(elements.size());
		for(unsigned int j=0;j<mln.clauses[i]->atoms.size();j++)
		{
			int ind = inCluster(mln.clauses[i]->atoms[j]);
			if(ind!=-1)
			{
				//collect its terms
				for(unsigned int k=0;k<mln.clauses[i]->atoms[j]->terms.size();k++)
				{
					terms[ind].push_back(mln.clauses[i]->atoms[j]->terms[k]);
				}
				for(unsigned int k=0;k<mln.clauses[i]->atoms.size();k++)
				{
					if(inCluster(mln.clauses[i]->atoms[k])!=-1)
						continue;
					for(unsigned int m=0;m<mln.clauses[i]->atoms[k]->terms.size();m++)
					{
						int matchPos = -1;
						for(unsigned int n=0;n<terms[ind].size();n++)
						{
							if(terms[ind].at(n)->domain.size()> 1 && terms[ind].at(n) == mln.clauses[i]->atoms[k]->terms[m])
							{
								matchPos = n;
								break;
							}
						}
						if(matchPos!=-1)
						{
							//shared position
							//get the shared position
							int sharePosition = 0;
							for(unsigned int jj=0;jj<mln.clauses[i]->atoms.size();jj++)
							{
								if(elements[ind]->symbol->id != mln.clauses[i]->atoms[jj]->symbol->id)
									continue;
								int pos = -1;
								for(unsigned int kk=0;kk<mln.clauses[i]->atoms[jj]->terms.size();kk++)
								{
									if(mln.clauses[i]->atoms[jj]->terms[kk]==
										terms[ind].at(matchPos))
									{
										pos = kk;
										break;
									}
								}
								if(pos!=-1)
								{
									sharePosition = pos;
									break;
								}
							}
							//sharedPos[ind].at(matchPos) = true;
							sharedPos[ind].at(sharePosition) = true;
						}
					}
				}
			}
		}
	}

	//find which clusters the shared positions are shared with
	cleanup(sharedTerms);
	sharedTerms.resize(sharedPos.size());
	for(unsigned int i=0;i<sharedPos.size();i++)
	{
		sharedTerms[i] = new LTermShare(sharedPos[i].size());
		sharedTerms[i]->shared = sharedPos[i];
	}

	
	//cleanup(assignments);
	while(!lAssignments.empty())
	{
		cleanup(lAssignments.back());
		lAssignments.pop_back();
	}
	lAssignments.resize(elements.size());
	totalTrues.clear();
	totalTrues.resize(elements.size());
	//permute shared positions and store groundings
	for(unsigned int i=0;i<elements.size();i++)
	{
		vector<LvrTerm*> sTerms;
		int unsharedSize = 1;
		for(unsigned int j=0;j<sharedPos[i].size();j++)
		{
			if(sharedPos[i].at(j))
			{
				sTerms.push_back(elements[i]->terms[j]);
			}
			else
				unsharedSize *= elements[i]->terms[j]->domain.size();
		}
		
		if(sTerms.size()>0)
		{
			vector<vector<int> > groundings;
			LvrPermutations::permuteTerms(sTerms,groundings);
			for(unsigned int j=0;j<groundings.size();j++)
			{
				Atom* tmpAtom = createGroundedAtom(i,sharedPos[i],groundings[j]);
				int val = LvRandomGenUtil::Instance()->getRandomPosition(unsharedSize+1);
				LVAssignment* lAssg = new LVAssignment(tmpAtom,val);
				lAssignments[i].push_back(lAssg);
				totalTrues[i] += val;
				delete tmpAtom;
			}
		}
		else
		{
			int val = LvRandomGenUtil::Instance()->getRandomPosition(unsharedSize+1);
			LVAssignment* lAssg = new LVAssignment(elements[i],val);
			lAssignments[i].push_back(lAssg);
			totalTrues[i] += val;
		}
	}
	partitionPositions();	
}

void LVRCluster::groundSharedTermsInCNF(vector<WClause*> &clauses)
{
	for(unsigned int i=0;i<clauses.size();i++)
	{
		//ground every term part of the cluster (not already grounded) that is shared within the clause
		vector<LvrTerm*> termsToGround;
		vector<int> atomPos;
		vector<int> termPos;
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			int ind = inCluster(clauses[i]->atoms[j]);
			if(ind!=-1)
			{
				for(unsigned int m=0;m<clauses[i]->atoms.size();m++)
				{
					int ind = inCluster(clauses[i]->atoms[m]);
					if(ind!=-1)
						continue;
					for(unsigned int k=0;k<clauses[i]->atoms[j]->terms.size();k++)
					{
						if(clauses[i]->atoms[j]->terms[k]->domain.size()==1)
							continue;
						for(unsigned int k1=0;k1<clauses[i]->atoms[m]->terms.size();k1++)
						{
							if(clauses[i]->atoms[j]->terms[k] == clauses[i]->atoms[m]->terms[k1])
							{
								//term k needs to be grounded
								if(find(termsToGround.begin(),termsToGround.end(),clauses[i]->atoms[j]->terms[k])
									== termsToGround.end())
								{
									termsToGround.push_back(clauses[i]->atoms[j]->terms[k]);
									termPos.push_back(k);
									atomPos.push_back(j);
								}
							}
						}
					}
				}
				//ground all permutations of terms
				if(termsToGround.size()!=0)
				{
					WClause* clauseToGround =  LvrMLN::create_new_clause(clauses[i]);
					vector<vector<int> > permutedDoms;
					LvrPermutations::permuteTerms(termsToGround,permutedDoms);
					for(unsigned int jj=0;jj<permutedDoms.size();jj++)
					{
						WClause* nClause =  LvrMLN::create_new_clause(clauseToGround);
						for(unsigned int kk=0;kk<permutedDoms[jj].size();kk++)
						{
							nClause->atoms[atomPos[kk]]->terms[termPos[kk]]->domain.clear();
							nClause->atoms[atomPos[kk]]->terms[termPos[kk]]->domain.push_back(permutedDoms[jj].at(kk));
						}
						//replace original clause with first grounding
						if(jj==0)
						{
							delete clauses[i];
							clauses[i] = nClause;
						}
						//push rest of the groundings back
						else
						{
							clauses.push_back(nClause);
						}
					}
				}
			}
		}
	}
}

bool LVRCluster::doExactInference()
{
	for(unsigned int i=0;i<sharedTerms.size();i++)
	{
		for(unsigned int j=0;j<sharedTerms[i]->shared.size();j++)
		{
			if(sharedTerms[i]->shared[j])
				return false;
		}
	}
	return true;
}

int LVRCluster::createNewPredicateId(int elementInd,Atom* groundedAtom)
{
	vector<int> sharedTermValues;
	//get the domain sizes of the shared terms
	vector<int> sharedDomainSizes;
	for(unsigned int i=0;i<sharedTerms[elementInd]->shared.size();i++)
	{
		if(sharedTerms[elementInd]->shared[i])
		{
			sharedDomainSizes.push_back(elements[elementInd]->terms[i]->domain.size());
			sharedTermValues.push_back(groundedAtom->terms[i]->domain[0]);
		}
	}
	if(sharedTermValues.size()==0)
		return 0;
	int val = sharedTermValues[sharedTermValues.size()-1];
	for(int j=sharedDomainSizes.size()-2;j>=0;j--)
	{
		for(unsigned int k=j+1;k<sharedDomainSizes.size();k++)
			val += sharedTermValues[j]*sharedDomainSizes[k];
	}

	return val;
}
/*
void LVRCluster::renamePredicates(vector<WClause*>& clauses)
{
	//int upperBoundSize = 0;
	//for(unsigned int i=0;i<clauses.size();i++)
	//{
		//upperBoundSize += clauses[i]->getNumberOfGroundedClauses();
	//}
	int upperBoundSize = clauses.size();
	//reset mln predicae ids
	mln.max_predicate_id = 0;
	//make new predicateids depending on signature of shared terms in predicate
	for(unsigned int i=0;i<clauses.size();i++)
	{
		for(unsigned int kk=0;kk<clauses[i]->atoms.size();kk++)
		{
			int ind = inCluster(clauses[i]->atoms[kk]);
			if(ind!=-1)
			{
				int pId = createNewPredicateId(ind,clauses[i]->atoms[kk]);
				clauses[i]->atoms[kk]->symbol->id = pId + ind*upperBoundSize;
				mln.setMaxPredicateId(clauses[i]->atoms[kk]->symbol->id+1);
			}
		}
	}
}
*/
void LVRCluster::renamePredicates(vector<WClause*>& clauses)
{
	vector<map<int,int> > idGroundingMap(elements.size());
	vector<int> currIds(elements.size());
	//assuming no self joins max number of occurences/predicate =  size of CNF
	for(unsigned int i=0;i<currIds.size();i++)
	{
		currIds[i] = i*clauses.size();
	}
	//reset mln predicae ids
	mln.max_predicate_id = 0;
	//make new predicateids depending on signature of shared terms in predicate
	for(unsigned int i=0;i<clauses.size();i++)
	{
		for(unsigned int kk=0;kk<clauses[i]->atoms.size();kk++)
		{
			int ind = inCluster(clauses[i]->atoms[kk]);
			if(ind!=-1)
			{
				int pId = createNewPredicateId(ind,clauses[i]->atoms[kk]);
				map<int,int>::iterator it = idGroundingMap[ind].find(pId);
				if(it!=idGroundingMap[ind].end())
				{
					clauses[i]->atoms[kk]->symbol->id = it->second;
				}
				else
				{
					clauses[i]->atoms[kk]->symbol->id = currIds[ind]++;
					idGroundingMap[ind].insert(pair<int,int>(pId,clauses[i]->atoms[kk]->symbol->id));
					mln.setMaxPredicateId(clauses[i]->atoms[kk]->symbol->id+1);
				}
			}
		}
	}
}


void LVRCluster::groundAllSharedTermsInCNF(vector<WClause*> clauses,vector<WClause*>& groundedClauses)
{
	for(unsigned int i=0;i<clauses.size();i++)
	{
		vector<LvrTerm*> termsToGround;
		vector<int> atomPos;
		vector<int> termPos;
		bool found = false;
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			int ind = inCluster(clauses[i]->atoms[j]);
			if(ind!=-1)
			{
				for(unsigned int k=0;k<clauses[i]->atoms[j]->terms.size();k++)
				{
					if(clauses[i]->atoms[j]->terms[k]->domain.size()==1)
						continue;
					if(!sharedTerms[ind]->shared[k])
						continue;
					if(find(termsToGround.begin(),termsToGround.end(),clauses[i]->atoms[j]->terms[k])
						== termsToGround.end())
					{
						termsToGround.push_back(clauses[i]->atoms[j]->terms[k]);
						termPos.push_back(k);
						atomPos.push_back(j);
					}
				}
				found=true;
			}
		}
		if(!found)
			continue;
		//ground all permutations of terms
		if(termsToGround.size() == 0)
		{
			WClause* nClause =  LvrMLN::create_new_clause(clauses[i]);
			groundedClauses.push_back(nClause);
			continue;
		}
		vector<vector<int> > permutedDoms;
		LvrPermutations::permuteTerms(termsToGround,permutedDoms);
		for(unsigned int jj=0;jj<permutedDoms.size();jj++)
		{
			WClause* nClause =  LvrMLN::create_new_clause(clauses[i]);
			for(unsigned int kk=0;kk<permutedDoms[jj].size();kk++)
			{
				nClause->atoms[atomPos[kk]]->terms[termPos[kk]]->domain.clear();
				nClause->atoms[atomPos[kk]]->terms[termPos[kk]]->domain.push_back(permutedDoms[jj].at(kk));
			}
			groundedClauses.push_back(nClause);
		}
	}
}

void LVRCluster::getNumTrueAssignments(int elementIndex, Atom* atom,int& nLVTrue,int& nLVFalse)
{
	nLVTrue = 0;
	nLVFalse = 0;
	/*
	if(atom->isConstant() && elements[elementIndex]->terms.size()==1)
	{
		nLVTrue = lAssignments[elementIndex].at(atom->terms[0]->domain[0])->nTrueValues;
		nLVFalse = lAssignments[elementIndex].at(atom->terms[0]->domain[0])->nFalseValues;
		return;
	}
	*/
	if(LvrHashAlgorithm::convertToHash(atom)==LvrHashAlgorithm::convertToHash(elements[elementIndex]))
	{
		nLVTrue = totalTrues[elementIndex];
		nLVFalse = atom->getNumberOfGroundings() - nLVTrue;
		return;
	}
	//find out how many values would match with atom's partial grounding
	int sz = 1;
	for(unsigned int i=0;i<sharedTerms[elementIndex]->shared.size();i++)
	{
		//stored as grounded, but need lifted
		if(sharedTerms[elementIndex]->shared[i] && atom->terms[i]->domain.size() != 1)
		{
			sz *= atom->terms[i]->domain.size();
		}
	}

	int matches = 0;
	int iterationcnt=0;
	for(unsigned int i=0;i<lAssignments[elementIndex].size();i++)
	{
		bool found = true;
		iterationcnt++;
		for(unsigned int j=0;j<lAssignments[elementIndex].at(i)->grounding.size();j++)
		{
			if(atom->terms[j]->domain.size() == 1)
			{
				if(lAssignments[elementIndex].at(i)->grounding[j] != atom->terms[j]->domain[0])
				{
					found = false;
					break;
				}
			}
		}
		if(found)
		{
			matches++;
			nLVTrue += lAssignments[elementIndex].at(i)->nTrueValues;
			nLVFalse += lAssignments[elementIndex].at(i)->nFalseValues;
			
			if(matches ==  sz)
				break;
		}
	}
	//cout<<iterationcnt<<"$";
	//atom->print();cout<<"$";elements[elementIndex]->print();
	//cout<<endl;
	//if(iterationcnt > 100)
		//cout<<"$$"<<iterationcnt<<endl;
}
void LVRCluster::reduceCNFByEvidence(vector<WClause*>& clauses)
{
	//cout<<"Sz="<<elements.size()*clauses.size()<<endl;
	evidenceReductionCost = 0;
	for(unsigned int i=0;i<elements.size();i++)
	{
		evidenceReductionCost += clauses.size();
		vector<WClause*> groundedClauses;
		doReduceCNFByEvidence(i,clauses,groundedClauses);
		clauses.clear();
		clauses = groundedClauses;
	}
}




void LVRCluster::doReduceClauseWithSelfJoinsByEvidence(int elementIndex,WClause* clause,vector<WClause*>& groundedClauses)
{
	//partition into equivalent groundings all the elements that are in cluster
	vector<vector<int> > groundings;
	vector<vector<int> > atomIndexes;
	vector<int> sameSign;
	vector<int> groundingSizes;
	for(unsigned int j=0;j<clause->atoms.size();j++)
	{
		int ind = inCluster(clause->atoms[j]);
		if(ind == elementIndex)
		{
			int mIndex = -1;
			for(unsigned int k=0;k<groundings.size();k++)
			{
				bool found=true;
				for(unsigned int m=0;m<groundings[k].size();m++)
				{
					if(groundings[k].at(m) != clause->atoms[j]->terms[m]->domain[0])
					{
						found = false;
						break;
					}
				}
				if(found)
				{
					mIndex = k;
					break;
				}
			}
			if(mIndex==-1)
			{
				vector<int> tmpDom;
				int sz = 1;
				for(unsigned int jj=0;jj<clause->atoms[j]->terms.size();jj++)
				{
					sz *= clause->atoms[j]->terms[jj]->domain.size();
					if(clause->atoms[j]->terms[jj]->domain.size()>1)
						continue;
					tmpDom.push_back(clause->atoms[j]->terms[jj]->domain[0]);
				}
				groundingSizes.push_back(sz);
				groundings.push_back(tmpDom);
				vector<int> tmpInd;
				tmpInd.push_back(j);
				atomIndexes.push_back(tmpInd);
				if(clause->sign[j])
					sameSign.push_back(-1);
				else
					sameSign.push_back(1);
			}
			else
			{
				atomIndexes[mIndex].push_back(j);
				if(sameSign[mIndex]!=0)
				{
					if(clause->sign[j] && sameSign[mIndex]==1)
					{
						sameSign[mIndex]=0;
					}
					else if(!clause->sign[j] && sameSign[mIndex]==-1)
					{
						sameSign[mIndex]=0;
					}
				}
			}
		}
	}
	vector<int> nLVTrueList;
	vector<int> nLVFalseList;
	for(unsigned int i=0;i<atomIndexes.size();i++)
	{
		int nLVTrue;
		int nLVFalse;
		getNumTrueAssignments(elementIndex,clause->atoms[atomIndexes[i].at(0)],nLVTrue,nLVFalse);
		nLVTrueList.push_back(nLVTrue);
		nLVFalseList.push_back(nLVFalse);
	}
	int totalSatisfied = 0;
	int totalUnSatisfied = 0;
	for(unsigned int i=0;i<nLVTrueList.size();i++)
	{
		//+ve and -ve grounding in self join
		if(sameSign[i]==0)
		{
			//all satisfied
			int mFac=1;
			if(atomIndexes[i].size() > 1)
			{
				mFac = groundingSizes[i]*atomIndexes[i].size()-1;
			}
			totalSatisfied += nLVTrueList[i]* mFac;
		}
		//all +ve groundings
		else if(sameSign[i]==1)
		{
			totalSatisfied += nLVTrueList[i];
			totalUnSatisfied += nLVFalseList[i];
		}
		//all -ve groundings
		else if(sameSign[i]==-1)
		{
			totalSatisfied += nLVFalseList[i];
			totalUnSatisfied += nLVTrueList[i];
		}
	}
	WClause* satClause = LvrMLN::create_new_clause(clause);
	WClause* unsatClause = LvrMLN::create_new_clause(clause);
	int rmCnt = 0;
	for(unsigned int i=0;i<atomIndexes.size();i++)
	{
		for(unsigned int j=0;j<atomIndexes[i].size();j++)
		{
			satClause->removeAtom(atomIndexes[i].at(j) - rmCnt);
			unsatClause->removeAtom(atomIndexes[i].at(j) - rmCnt);
			rmCnt++;
		}
	}
	if(totalSatisfied > 0)
	{
		satClause->satisfied = true;
		satClause->weight = satClause->weight*LogDouble(totalSatisfied,false);
		groundedClauses.push_back(satClause);
	}
	else
	{
		delete satClause;
	}
	if(totalUnSatisfied > 0)
	{
		unsatClause->weight = unsatClause->weight*LogDouble(totalUnSatisfied,false);
		groundedClauses.push_back(unsatClause);
	}
	else
	{
		delete unsatClause;
	}
}

void LVRCluster::doReduceCNFByEvidence(int elementIndex,vector<WClause*>& clauses,vector<WClause*>& groundedClauses)
{
	int cnt = 0;
	for(unsigned int i=0;i<clauses.size();i++)
	{
		int atomInd = -1;
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			int ind = inCluster(clauses[i]->atoms[j]);
			if(ind==elementIndex)
			{
				atomInd = j;
				break;
			}
		}
		if(atomInd == -1)
		{
			//WClause* nClause = LvrMLN::create_new_clause(clauses[i]);
			//groundedClauses.push_back(nClause);
			groundedClauses.push_back(clauses[i]);
			//delete clauses[i];
			continue;
		}
		//set evidence
		if(clauses[i]->isSelfJoinedOnAtom(clauses[i]->atoms[atomInd]))
		{
			doReduceClauseWithSelfJoinsByEvidence(elementIndex,clauses[i],groundedClauses);
			delete clauses[i];
			continue;
		}

		int nLVTrue;
		int nLVFalse;
		getNumTrueAssignments(elementIndex,clauses[i]->atoms[atomInd],nLVTrue,nLVFalse);
		if(nLVTrue!=0 && nLVFalse!=0)
		{
			clauses[i]->weight = clauses[i]->weight*LogDouble(nLVTrue,false);
			bool nsat = false;
			if(!clauses[i]->sign[atomInd])
				clauses[i]->satisfied = true;
			else
				nsat = true;
			clauses[i]->removeAtom(atomInd);
			WClause* nClause = LvrMLN::create_new_clause(clauses[i]);
			nClause->weight = nClause->weight*LogDouble(nLVFalse,false);
			nClause->satisfied = nsat;
			groundedClauses.push_back(clauses[i]);
			groundedClauses.push_back(nClause);
		}
		else
		{
			//clauses[i]->weight = clauses[i]->weight*(LogDouble(nLVTrue,false)+LogDouble(nLVFalse,false));
			if((clauses[i]->sign[atomInd] && nLVTrue==0) ||
					(!clauses[i]->sign[atomInd] && nLVFalse==0))
				clauses[i]->satisfied = true;
			clauses[i]->removeAtomWithWeightUpdation(atomInd);
			groundedClauses.push_back(clauses[i]);
		}

/*		if(nLVTrue!=0)
		{
			WClause* nClause = LvrMLN::create_new_clause(clauses[i]);
			//update the weight
			if(!clauses[i]->sign[atomInd])
				nClause->satisfied = true;
			nClause->weight = nClause->weight*LogDouble(nLVTrue,false);
			nClause->removeAtom(atomInd);
			groundedClauses.push_back(nClause);
		}
		if(nLVFalse!=0)
		{
			WClause* nClause = LvrMLN::create_new_clause(clauses[i]);
			if(clauses[i]->sign[atomInd])
				nClause->satisfied = true;
			nClause->weight = nClause->weight*LogDouble(nLVFalse,false);
			nClause->removeAtom(atomInd);
			groundedClauses.push_back(nClause);
		}
		delete clauses[i];
*/
	}
	clauses.clear();
	//cout<<"NPO="<<cnt<<endl;
}

void LVRCluster::partitionPositions()
{
	sharedLiftedPositions.clear();
	sharedGroundedPositions.clear();
	sharedLiftedPositions.resize(elements.size());
	sharedGroundedPositions.resize(elements.size());
	for(unsigned int i=0;i<elements.size();i++)
	{
		for(unsigned int j=0;j<sharedTerms[i]->shared.size();j++)
		{
			if(sharedTerms[i]->shared[j])
				sharedGroundedPositions[i].push_back(j);
			else
				sharedLiftedPositions[i].push_back(j);
		}
	}
}

void LVRCluster::resetAllAssignments()
{
	for(unsigned int i=0;i<lAssignments.size();i++)
	{
		for(unsigned int j=0;j<lAssignments[i].size();j++)
		{
			int tot = lAssignments[i].at(j)->nTrueValues + lAssignments[i].at(j)->nFalseValues; 
			lAssignments[i].at(j)->nTrueValues = -1;
			lAssignments[i].at(j)->nFalseValues = tot;
		}
		totalTrues[i] = 0;
	}
}

void LVRCluster::setAllDontCareAssignments()
{
	for(unsigned int i=0;i<lAssignments.size();i++)
	{
		for(unsigned int j=0;j<lAssignments[i].size();j++)
		{
			if(lAssignments[i].at(j)->nTrueValues == -1)
			{
				int val = LvRandomGenUtil::Instance()->getRandomPosition(lAssignments[i].at(j)->nFalseValues+1);
				lAssignments[i].at(j)->nTrueValues = val;
				lAssignments[i].at(j)->nFalseValues -= val;
				totalTrues[i] += val;
			}
		}
	}
}
