#include "lisapproxinference.h"
#include "cleanuputils.h"
#include "rulesutil.h"
#include <time.h>
#include "mathutils.h"
#include "samplingalgs.h"
#include "queryupdater.h"


LogDouble LISApproxInference::CNFWeight(vector<WClause*>& CNF)
{
	/*for(unsigned int i=0;i<CNF.size();i++)
	{
		if(CNF[i]->satisfied)
			cout<<"Satisfied v";
		CNF[i]->print();
	}
	*/
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
			//totalWeight += CNF[i]->weight*LogDouble(sz,false);
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


bool LISApproxInference::decomposeCNF(vector<WClause*>& CNF,int& powerFactor)
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
			for(unsigned int j=0;j<decomposer_list[i]->decomposer_terms.size();j++)
			{
				//in all clauses that this term occurs change weight
				/*
				for(unsigned int jj=0;jj<CNF.size();jj++)
				{
					for(unsigned int kk=0;kk<CNF[jj]->atoms.size();kk++)
					{
						bool done = false;
						for(unsigned int mm=0;mm<CNF[jj]->atoms[kk]->terms.size();mm++)
						{
							if(CNF[jj]->atoms[kk]->terms[mm] == decomposer_list[i]->decomposer_terms[j])
							{
								//change weight
								CNF[jj]->weight = CNF[jj]->weight*LogDouble(domSize,false);
								done=true;
							}
						}
						if(done)
							break;
					}
				}
				*/
				//store pre-decomp domain
				decomposer_list[i]->decomposer_terms[j]->origDomain.clear();
				for(unsigned int jj=0;jj<decomposer_list[i]->decomposer_terms[j]->domain.size();jj++)
					decomposer_list[i]->decomposer_terms[j]->origDomain.push_back(decomposer_list[i]->decomposer_terms[j]->domain[jj]);

				//Erase all but one grounding
				decomposer_list[i]->decomposer_terms[j]->domain.erase(decomposer_list[i]->decomposer_terms[j]->domain.begin(),decomposer_list[i]->decomposer_terms[j]->domain.begin()+decomposer_list[i]->decomposer_terms[j]->domain.size()-1);
				//insert the original size into the term
				decomposer_list[i]->decomposer_terms[j]->origDomainSize = domSize;
			}
		}
	}
	cleanup(decomposer_list);
	return true;
}
////IMPLEMENTATION1: NO DISJOINT SEPERATION MORE GROUNDING, faster in some cases
LogDouble LISApproxInference::doLvApproxPartitionInformedV1(vector<WClause*>& CNF)
{
	LogDouble totalVal(1,false);
	if(CNF.size()==0)
		return totalVal;
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
	
	int powerFactor;
	bool isDecomposed = decomposeCNF(CNF,powerFactor);
	if(isDecomposed)
	{
		LogDouble mcnt = doLvApproxPartitionInformedV1(CNF);
		LogDouble val = LogDouble(1,false);
		LogDouble::LDPower(mcnt,powerFactor,val);
		totalVal = totalVal*val;
		return totalVal;
	}
	
	Atom* tmpAtom = heuristics->getAtomToSplit(CNF);	
	if(tmpAtom==NULL)
	{	
		LogDouble wt1 = CNFWeight(CNF);
		cleanup(CNF);
		return totalVal*wt1;
			
	}
	Atom* atom = LvrMLN::create_new_atom(tmpAtom);
	//sample according to mode
	LProposalDistributionElement* lpe = lpe = distribution->findElement(atom);
	//use Isolated terms rule
	vector<bool> isolatedTerms;
	bool isIsolated = LRulesUtil::computeIsolatedTerms(atom,CNF,isolatedTerms);
	double binCoeff = 1;
	LogDouble probOfSample(1,false);
	LogDouble sampleWeight(1,false);
	vector<int> sampledValues;
	if(isIsolated)
	{
		int isolatedSize=1;
		int nonisolatedSize = 1;
		for(unsigned int jj=0;jj<isolatedTerms.size();jj++)
		{
			if(isolatedTerms[jj])
			{
				isolatedSize *= atom->terms[jj]->domain.size();
			}
			else
			{
				nonisolatedSize *=  atom->terms[jj]->domain.size();
			}
		}
		LSamplingAlgs::Instance()->sample(isolatedSize,nonisolatedSize,sampleWeight,sampledValues,samplingMode,lpe);
		if(LvrQueryUpdater::isInstanceCreated())
		{
			LvrQueryUpdater::Instance()->updateQueryValues(atom,isolatedTerms,sampledValues);
		}
	}
	else
	{
		int numGroundings = atom->getNumberOfGroundings();
		LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues,samplingMode,lpe);
		if(LvrQueryUpdater::isInstanceCreated())
		{
			LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
		}
		
	}
	lvrNormPropagate->propagateNormalizedCNF(CNF,atom,isolatedTerms,sampledValues);
	LogDouble mcnt = doLvApproxPartitionInformedV1(CNF);
	totalVal *= sampleWeight*mcnt;
	delete atom;
	cleanup(CNF);
	return totalVal;
}


LogDouble LISApproxInference::doLvApproxPartitionInformed(vector<WClause*>& CNF1)
{
	LogDouble totalVal(1,false);
	if(CNF1.size()==0)
		return totalVal;
	vector<vector<WClause*> > dCNFs;
	LClusterUtil::seperateDisjointCNF(CNF1,dCNFs);
	cleanup(CNF1);
	for(unsigned int t=0;t<dCNFs.size();t++)
	{
		vector<WClause*> CNF = dCNFs[t];
		int powerFactor;
		bool isDecomposed = decomposeCNF(CNF,powerFactor);
		if(isDecomposed)
		{
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			continue;
		}		
		Atom* tmpAtom = heuristics->getAtomToSplit(CNF);	
		if(tmpAtom==NULL)
		{	
			LogDouble wt1 = CNFWeight(CNF);
			cleanup(CNF);
		    totalVal = totalVal*wt1;
			continue;
			
		}
		Atom* atom = LvrMLN::create_new_atom(tmpAtom);
		LProposalDistributionElement* lpe = lpe = distribution->findElement(atom);
		//check for nonpoly cases
		//check for self joins
		bool selfJoined = false;
		for(int m=0;m<CNF.size();m++)
		{
			if(CNF.at(m)->isSelfJoinedOnAtom(atom))
			{
				selfJoined = true;
				break;
			}
		}
		int singletonIndex;
		bool singleton = atom->isSingletonAtom(singletonIndex);
		bool nonpoly = false;
		if(!singleton)
			nonpoly = true;
		else if(selfJoined)
		{
			//check for blocking rule
			bool blocked = LRulesUtil::isBlocked(atom,singletonIndex,CNF);
			if(!blocked)
				nonpoly=true;
		}
		if(!nonpoly)
		{
			LogDouble probOfSample(1,false);
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			int domSize = atom->terms[singletonIndex]->domain.size();
			LSamplingAlgs::Instance()->sample(domSize,1,sampleWeight,sampledValues,samplingMode,lpe);
			if(LvrQueryUpdater::isInstanceCreated())
			{
				LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
			}
			LvrSingletonNormPropagation::propagateNormalizedCNF(CNF,atom,singletonIndex,sampledValues[0]);
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			totalVal *= sampleWeight*mcnt;

		}
		else
		{
			//use Isolated terms rule
			vector<bool> isolatedTerms;
			bool isIsolated = LRulesUtil::computeIsolatedTerms(atom,CNF,isolatedTerms);
			double binCoeff = 1;
			LogDouble probOfSample(1,false);
			LogDouble sampleWeight(1,false);
			vector<int> sampledValues;
			if(isIsolated)
			{
				int isolatedSize=1;
				int nonisolatedSize = 1;
				for(unsigned int jj=0;jj<isolatedTerms.size();jj++)
				{
					if(isolatedTerms[jj])
					{
						isolatedSize *= atom->terms[jj]->domain.size();
					}
					else
					{
						nonisolatedSize *=  atom->terms[jj]->domain.size();
					}
				}
				LSamplingAlgs::Instance()->sample(isolatedSize,nonisolatedSize,sampleWeight,sampledValues,samplingMode,lpe);
				if(LvrQueryUpdater::isInstanceCreated())
				{
					LvrQueryUpdater::Instance()->updateQueryValues(atom,isolatedTerms,sampledValues);
				}
			}
			else
			{
				int numGroundings = atom->getNumberOfGroundings();
				LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues,samplingMode,lpe);
				if(LvrQueryUpdater::isInstanceCreated())
				{
					LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
				}
		
			}
			lvrNormPropagate->propagateNormalizedCNF(CNF,atom,isolatedTerms,sampledValues);
			LogDouble mcnt = doLvApproxPartitionInformed(CNF);
			totalVal *= sampleWeight*mcnt;
		}
		
		delete atom;
		cleanup(CNF);
	}
	return totalVal;
}

LogDouble LISApproxInference::doLvApproxPartitionBinomial(vector<WClause*>& CNF1)
{
	LogDouble totalVal(1,false);
	if(CNF1.size()==0)
		return totalVal;
	int completedCount=0;
	for(unsigned int i=0;i<CNF1.size();i++)
	{
		if(CNF1[i]->atoms.size()==0 || CNF1[i]->satisfied)
		{
			completedCount++;
		}
	}
	if(completedCount == CNF1.size())
	{
		LogDouble weight = CNFWeight(CNF1);
		cleanup(CNF1);
		return weight;
	}
	vector<vector<WClause*> > dCNFs;
	LClusterUtil::seperateDisjointCNF(CNF1,dCNFs);
	cleanup(CNF1);
	for(unsigned int t=0;t<dCNFs.size();t++)
	{
		vector<WClause*> CNF = dCNFs[t];
		int powerFactor;
		bool isDecomposed = decomposeCNF(CNF,powerFactor);
		if(isDecomposed)
		{
			LogDouble mcnt = doLvApproxPartitionBinomial(CNF);
			LogDouble val = LogDouble(1,false);
			LogDouble::LDPower(mcnt,powerFactor,val);
			totalVal = totalVal*val;
			continue;
		}		
		Atom* tmpAtom = heuristics->getAtomToSplit(CNF);	
		if(tmpAtom==NULL)
		{	
			LogDouble wt1 = CNFWeight(CNF);
			cleanup(CNF);
		    totalVal*wt1;
			continue;
			
		}
		Atom* atom = LvrMLN::create_new_atom(tmpAtom);
		LogDouble probOfSample(1,false);
		LogDouble sampleWeight(1,false);
		int numGroundings = atom->getNumberOfGroundings();
		//select the largest domain grounding
		int maxDomainTerm = 0;
		int maxDomain = -1;
		for(unsigned int jj=0;jj<atom->terms.size();jj++)
		{
			int up = atom->terms[jj]->domain.size();
			if(up > maxDomain)
			{
				maxDomain = atom->terms[jj]->domain.size();
				maxDomainTerm = jj;
			}
		}
		numGroundings = atom->terms[maxDomainTerm]->domain.size();
		vector<int> sampledValues;
		if(samplingMode==EBINOMIAL)
		{
			LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues,EBINOMIAL);
		}
		else
			LSamplingAlgs::Instance()->sample(numGroundings,1,sampleWeight,sampledValues);
		if(LvrQueryUpdater::isInstanceCreated())
		{
			LvrQueryUpdater::Instance()->updateQueryValues(atom,sampledValues[0]);
		}
		LvrSingletonNormPropagation::propagateNormalizedCNF(CNF,atom,maxDomainTerm,sampledValues[0]);
		LogDouble mcnt = doLvApproxPartitionBinomial(CNF);
		totalVal *= sampleWeight*mcnt;
		delete atom;
		cleanup(CNF);
	}
	return totalVal;
}

LogDouble LISApproxInference::estimatePartitionFunction(LvrParams* params,LProposalDistribution* distribution)
{
	cout<<"Lifted Importance Sampling for estimating Z..."<<endl;
	time_t start;
	time(&start);
	cout<<"Time ="<<ctime(&start)<<endl;
	LogDouble ZApprox;
	int iterations = 0;
	if(distribution)
		setProposalDistribution(distribution);
	setSamplingMode(params->samplingMode);
	if(params->learningRate <=0 )
		params->learningRate = LEARNINGRATE;
	if(params->proposalUpdateInterval <= 0)
		params->proposalUpdateInterval = PROPOSALUPDATEINTERVAL;

	while(1)
	{
		//make a copy of normalized clauses
		vector<WClause*> clauses;
		LvrMLN::copyAllClauses(mln.clauses,clauses);
		LogDouble currZ;
		//currZ = weightedModelCountApprox(clauses);
		if(params->samplingMode == EINFORMED)
			currZ = doLvApproxPartitionInformed(clauses);
		else if(params->samplingMode == EINFORMEDV1)
			currZ = doLvApproxPartitionInformedV1(clauses);
		else
			currZ = doLvApproxPartitionBinomial(clauses);
		//update approximation
		ZApprox = ZApprox*LogDouble(iterations,false) + currZ;
		iterations++;
		if(params->samplingMode == EINFORMED || params->samplingMode == EINFORMEDV1)
		{
			if(iterations % params->proposalUpdateInterval == 0)
				distribution->updateDistributions(params->learningRate);
		}
		//take average
		ZApprox = ZApprox / LogDouble(iterations,false);
		time_t curr;
		time(&curr);
		int seconds = difftime(curr,start);
		if(iterations%PRINTRESULTSINTERVAL == 0 ||
				seconds > params->maxSeconds || iterations >= params->maxSteps)
		{
		cout<<"iteration="<<iterations<<", currZ : ";
		currZ.printValue();
		cout<<", ZApprox : ";
		ZApprox.printValue();
		cout<<endl;
		if(seconds > params->maxSeconds || iterations >= params->maxSteps)
			return ZApprox;
		}
	}
}

void LISApproxInference::estimateApproxMarginals(LvrParams* params,LProposalDistribution* distribution)
{
	setProposalDistribution(distribution);
	cout<<"Estimating Marginals using Lifted Importance Sampling..."<<endl;
	time_t start;
	time(&start);
	cout<<"Time ="<<ctime(&start)<<endl;
	LogDouble ZApprox;
	int iterations = 0;
	int printInterval = PRINTRESULTSINTERVAL;
	setSamplingMode(params->samplingMode);
	LogDouble totalWeight;
	if(params->learningRate <=0 )
		params->learningRate = LEARNINGRATE;
	if(params->proposalUpdateInterval <= 0)
		params->proposalUpdateInterval = PROPOSALUPDATEINTERVAL;
	while(1)
	{
		//make a copy of normalized clauses
		vector<WClause*> clauses;
		LvrMLN::copyAllClauses(mln.clauses,clauses);
		LogDouble currWt;
		if(params->samplingMode == EINFORMED)
			currWt = doLvApproxPartitionInformed(clauses);
		else if(params->samplingMode == EINFORMEDV1)
			currWt = doLvApproxPartitionInformedV1(clauses);
		else
			currWt = doLvApproxPartitionBinomial(clauses);
		LvrQueryUpdater::Instance()->updateDontCare();
		clauses.clear();
		//update the weights for query
		LvrQueryUpdater::Instance()->updateAllImportanceWeights(currWt);

		//update the cumulative weight
		totalWeight += currWt;
		//update approximation
		iterations++;
		if(params->samplingMode == EINFORMED || params->samplingMode == EINFORMEDV1)
		{
			if(iterations % params->proposalUpdateInterval == 0)
			{
				distribution->updateDistributions(params->learningRate);
			}
		}
		time_t curr;
		time(&curr);
		int seconds = difftime(curr,start);
		if(iterations % PRINTRESULTSINTERVAL == 0 ||
				seconds > params->maxSeconds || iterations >= params->maxSteps)
		{
			LvrQueryUpdater::Instance()->writeToFile(totalWeight);
			cout<<"iteration "<< iterations<<endl;
			cout<<"Z-curr = ";
			currWt.printValue();
			cout<<endl;
			cout<<"cumulative-Z = ";
			totalWeight.printValue();
			cout<<endl;
			if(seconds > params->maxSeconds || iterations >= params->maxSteps)
				break;
		}
	}
}

LISApproxInference::LISApproxInference(LvrMLN& mln_): mln(mln_)
{
	decomposer = new LDecomposer(mln);
	heuristics = new LHeuristics(*decomposer,mln);
	lvrNormPropagate = new LvrNormPropagation(mln);
}
