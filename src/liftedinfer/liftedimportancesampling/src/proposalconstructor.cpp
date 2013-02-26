#include "proposalconstructor.h"
#include "heuristics.h"
#include "lvratomhashtemplate.h"


LProposalConstructor::LProposalConstructor(LvrMLN& mln_):mln(mln_)
{
	lpd = NULL;
	lc = new LClusterUtil();
	//ptpsearch = new LPTPSearch(mln);
	lsearch = new LISApproxInference(mln);
	lpg = new LvrNormPropagation(mln);
	lvrSingletonNormPropagation = new LvrSingletonNormPropagation();
}

LProposalConstructor::~LProposalConstructor()
{
	if(lpd)
		delete lpd;
	delete lc;
	delete lsearch;
	delete lpg;
	delete lvrSingletonNormPropagation;
}

void LProposalConstructor::getParents(vector<WClause*> clauses,Atom* atom, vector<Atom*> potentialParents,vector<Atom*>& parents)
{
	//select k among the potential parents
	vector<int> counts(potentialParents.size());
	/*vector<int> parentHashes(potentialParents.size());
	for(unsigned int i=0;i<potentialParents.size();i++)
		parentHashes[i] = LvrHashAlgorithm::convertToHash(potentialParents[i]);
	unsigned int atomHash = LvrHashAlgorithm::convertToHash(atom);
	for(unsigned int i=0;i<clauses.size();i++)
	{
		bool found=false;
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			//check if atom appears in clause
			
			if(LvrHashAlgorithm::convertToHash(clauses[i]->atoms[j]) != atomHash)
			{
				found=true;
				break;
			}
		}
		if(!found)
			continue;
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			vector<int>::iterator it = find(parentHashes.begin(),parentHashes.end(),LvrHashAlgorithm::convertToHash(clauses[i]->atoms[j]));
			int pos = it - parentHashes.begin();
			if(it!=parentHashes.end())
			{
				counts[pos]++;
			}
		}
	}*/
	for(unsigned int i=0;i<potentialParents.size();i++)
	{
		counts[i] = LHeuristics::proposalCouplingScore(mln.clauses,atom,potentialParents[i]);
	}
	//cout<<"counts="<<counts.size()<<endl;
	//get the k max count potential parents
	vector<int> index(counts.size());
	for(unsigned int i=0;i<index.size();i++)
		index[i] = i;
	//sort in descending order
	for(unsigned int i=0;i<counts.size();i++)
	{
		for(unsigned int j=0;j<counts.size()-i-1;j++)
		{
			if(counts[j] < counts[j+1])
			{
				int tmp=counts[j];
				counts[j]=counts[j+1];
				counts[j+1]=tmp;
				tmp=index[j];
				index[j]=index[j+1];
				index[j+1]=tmp;			
			}
		}
	}
	int maxVal = PROPOSALNUMPARENTS;
	if(counts.size() < maxVal)
		maxVal = counts.size();
	for(unsigned int i=0;i<maxVal;i++)
	{
		parents.push_back(potentialParents[index[i]]);
	}
}

void LProposalConstructor::constructProposal(vector<WClause*>& CNF1,vector<Atom*>& potentialParents)
{
	if(CNF1.size()==0)
		return;
	vector<vector<WClause*> > dCNFs;
	LClusterUtil::seperateDisjointCNF(CNF1,dCNFs);
	cleanup(CNF1);
	for(unsigned int t=0;t<dCNFs.size();t++)
	{
		vector<WClause*> CNF = dCNFs[t];

		int powerFactor;
		bool isDecomposed = lsearch->decomposeCNF(CNF,powerFactor);
		if(isDecomposed)
		{
			//use power rule
			constructProposal(CNF,potentialParents);
			continue;
		}

		Atom* tmpSelectedAtom = lsearch->selectAtomToCondition(CNF);
		if(tmpSelectedAtom == NULL)
		{
			cleanup(CNF);
			continue;
		}
		Atom* selectedAtom = LvrMLN::create_new_atom(tmpSelectedAtom);
		vector<Atom*> parents;		
		if(potentialParents.size() > 0)
		{
			getParents(CNF,selectedAtom,potentialParents,parents);	
		}
		//see if atom has isolated terms
		vector<bool> isolatedTerms;
		bool isIsolated = LRulesUtil::computeIsolatedTerms(selectedAtom,CNF,isolatedTerms);
		lpd->insertElement(selectedAtom,parents,isolatedTerms);
		potentialParents.push_back(selectedAtom);
		//remove atom from all clauses
		for(unsigned int i=0;i<CNF.size();i++)
		{
			for(int j=0;j<CNF[i]->atoms.size();j++)
			{
				if(CNF[i]->atoms[j]->symbol->id == selectedAtom->symbol->id)
				{
					CNF[i]->removeAtom(j);
					j--;
				}
			}
		}
		constructProposal(CNF,potentialParents);
		cleanup(CNF);
	}
}

/*
void LProposalConstructor::constructProposal(vector<WClause*>& CNF1,vector<Atom*>& potentialParents)
{
	if(CNF1.size()==0)
		return;
	vector<vector<WClause*> > dCNFs;
	LClusterUtil::seperateDisjointCNF(CNF1,dCNFs);
	cleanup(CNF1);
	for(unsigned int t=0;t<dCNFs.size();t++)
	{
		vector<WClause*> CNF = dCNFs[t];

		int powerFactor;
		bool isDecomposed = lsearch->decomposeCNF(CNF,powerFactor);
		if(isDecomposed)
		{
			//use power rule
			constructProposal(CNF,potentialParents);
			continue;
		}
		Atom* tmpSelectedAtom = lsearch->selectAtomToCondition(CNF);
		if(tmpSelectedAtom == NULL)
		{
			cleanup(CNF);
			continue;
		}
		Atom* selectedAtom = LvrMLN::create_new_atom(tmpSelectedAtom);
		vector<Atom*> parents;		
		if(potentialParents.size() > 0)
		{
			getParents(CNF,selectedAtom,potentialParents,parents);	
		}
		//see if atom has isolated terms
		vector<bool> isolatedTerms;
		bool isIsolated = LRulesUtil::computeIsolatedTerms(selectedAtom,CNF,isolatedTerms);
		lpd->insertElement(selectedAtom,parents,isolatedTerms);
		potentialParents.push_back(selectedAtom);
		//lpg->propagateNormalizedCNF(CNF,selectedAtom,isolatedTerms);
		//check for self joins
		bool selfJoined = false;
		for(int m=0;m<CNF.size();m++)
		{
			if(CNF.at(m)->isSelfJoinedOnAtom(selectedAtom))
			{
				selfJoined = true;
				break;
			}
		}

		int singletonIndex;
		bool singleton = selectedAtom->isSingletonAtom(singletonIndex);
		bool nonpoly = false;
		if(!singleton)
			nonpoly = true;
		else if(selfJoined)
		{
			//check for blocking rule
			bool blocked = LRulesUtil::isBlocked(selectedAtom,singletonIndex,CNF);
			if(!blocked)
				nonpoly=true;
		}
		if(!nonpoly)
		{
			selectedAtom->print();
			int domSize = selectedAtom->terms[singletonIndex]->domain.size();
			for(int jj=0;jj<=domSize;jj++)
			{
				vector<WClause*> tmpClauses;
				LvrMLN::copyAllClauses(CNF,tmpClauses);
				LvrSingletonNormPropagation::propagateNormalizedCNF(tmpClauses,selectedAtom,singletonIndex,jj);
				constructProposal(tmpClauses,potentialParents);
				cleanup(tmpClauses);
			}
		}
		else
		{
			lpg->propagateNormalizedCNF(CNF,selectedAtom,isolatedTerms);
			constructProposal(CNF,potentialParents);
		}
		//delete selectedAtom;
		cleanup(CNF);
	}
}
*/
////IMPLEMENTATION1: NO DISJOINT SEPERATION MORE GROUNDING, faster in some cases
void LProposalConstructor::constructProposalV1(vector<WClause*>& CNF,vector<Atom*>& potentialParents)
{
	if(CNF.size()==0)
		return;
	int powerFactor;
	bool isDecomposed = lsearch->decomposeCNF(CNF,powerFactor);
	if(isDecomposed)
	{
		//use power rule
		constructProposalV1(CNF,potentialParents);
		return;
	}
	Atom* tmpSelectedAtom = lsearch->selectAtomToCondition(CNF);
	if(tmpSelectedAtom == NULL)
	{
		cleanup(CNF);
		return;
	}
	Atom* selectedAtom = LvrMLN::create_new_atom(tmpSelectedAtom);
	vector<Atom*> parents;		
	if(potentialParents.size() > 0)
	{
		getParents(CNF,selectedAtom,potentialParents,parents);	
	}
	//see if atom has isolated terms
	vector<bool> isolatedTerms;
	bool isIsolated = LRulesUtil::computeIsolatedTerms(selectedAtom,CNF,isolatedTerms);
	lpd->insertElement(selectedAtom,parents,isolatedTerms);
	potentialParents.push_back(selectedAtom);
	lpg->propagateNormalizedCNF(CNF,selectedAtom,isolatedTerms);
	constructProposalV1(CNF,potentialParents);
	cleanup(CNF);
}

void LProposalConstructor::startConstruction(LvrParams* params)
{
	cout<<"Constructing proposal distribution structure..."<<endl;
	vector<Atom*> potentialParents;
	vector<WClause*> nClauses;
	LvrMLN::copyAllClauses(mln.clauses,nClauses);
	lpd = new LProposalDistribution(params->samplingMode);
	if(params->samplingMode == EINFORMED)
	{
		constructProposal(nClauses,potentialParents);	
	}
	else
	{
		constructProposalV1(nClauses,potentialParents);	
	}
}

void LProposalConstructor::startMARInference(LvrParams* params)
{
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
	if(params->samplingMode == EINFORMED || params->samplingMode == EINFORMEDV1)
	{
		startConstruction(params);
		lsearch->estimateApproxMarginals(params,lpd);
	}
	else
		lsearch->estimateApproxMarginals(params);
}

void LProposalConstructor::startPartitionFunction(LvrParams* params)
{
	//set defaults
	int printInterval = PRINTRESULTSINTERVAL;
	if(params->maxSteps <= 0)
	{
		params->maxSteps = MAXSTEPSINFER;
	}
	if(params->maxSeconds <= 0)
	{
		params->maxSeconds = MAXSECONDSINFER;
	}
	if(params->samplingMode == EINFORMED || params->samplingMode == EINFORMEDV1)
	{
		startConstruction(params);
		lsearch->estimatePartitionFunction(params,lpd);
	}
	else
		lsearch->estimatePartitionFunction(params);
}

void LProposalConstructor::readDistributionFromDumpFile()
{
	lpd->readDistributionFromDumpFile();
}

void LProposalConstructor::dumpDistributionToDumpFile()
{
	lpd->dumpDistributionToFile();
}
