#include "proposalstructure.h"
#include "filedump.h"
#include "cleanuputils.h"
#include "randomgenutil.h"

	LProposalTable::LProposalTable(int parentsSize,int maxTruthValue_):totalSamplesCollected(0),maxTruthValue(maxTruthValue_)
	{
		distribution.clear();
		adaptiveMLECounter.clear();
		totalSamplesCollected.clear();
		normConstants.clear();
		distribution.resize(parentsSize);
		adaptiveMLECounter.resize(parentsSize);
		totalSamplesCollected.resize(parentsSize);
		normConstants.resize(parentsSize);
		for(unsigned int i=0;i<distribution.size();i++)
		{
			distribution[i].resize(maxTruthValue+1);
			adaptiveMLECounter[i].resize(maxTruthValue+1);
			for(unsigned int j=0;j<distribution[i].size();j++)
				distribution[i].at(j) = LogDouble(1.0,false)/LogDouble(maxTruthValue+1,false);
			//set the norm constant based on domain size
			normConstants[i] = distribution[i].at(0) + distribution[i].at(distribution[i].size()-1);
			LogDouble factor;
			for(unsigned int j=1;j<distribution[i].size()-1;j++)
				factor = factor + distribution[i].at(j);
			if(!factor.is_zero)
				factor = factor*LogDouble(maxTruthValue,false);
			else
				factor = LogDouble(maxTruthValue,false);
			normConstants[i] = normConstants[i] + factor;
		}	
	}
	void LProposalTable::updateDistribution(double learningRate)
	{
		for(unsigned int i=0;i<distribution.size();i++)
		{
			if(totalSamplesCollected[i] == 0)
				continue;
			double norm = 0;
			vector<double> newProbabilities(distribution[i].size());
			for(unsigned int j=0;j<distribution[i].size();j++)
			{
				newProbabilities[j] = (double)adaptiveMLECounter[i].at(j)/(double)totalSamplesCollected[i];
				//distribution[i].at(j) = (double)adaptiveMLECounter[i].at(j)/(double)totalSamplesCollected[i];
				//norm += distribution[i].at(j);
				norm += newProbabilities[j];
				//reset the mle counter for the next updation
				adaptiveMLECounter[i].at(j) = 0;
			}
			//Use adaptive IS using the gradient descent rule
			LogDouble newNorm;
			for(unsigned int j=0;j<distribution[i].size();j++)
			{
				distribution[i].at(j) = distribution[i].at(j) + LogDouble((newProbabilities[j]-exp(distribution[i].at(j).value))*learningRate,false);
				newNorm += distribution[i].at(j);
			}
			for(unsigned int j=0;j<distribution[i].size();j++)
			{
				distribution[i].at(j) = distribution[i].at(j)/newNorm;
			}
			//update the norm constants
			normConstants[i] = distribution[i].at(0) + distribution[i].at(distribution[i].size()-1);
			LogDouble factor;
			for(unsigned int j=1;j<distribution[i].size()-1;j++)
				factor = factor + distribution[i].at(j);
			if(!factor.is_zero)
				factor = factor*LogDouble(maxTruthValue,false);
			else
				factor = LogDouble(maxTruthValue,false);
			
			normConstants[i] = normConstants[i] + factor;

			//RESET the sample counter
			totalSamplesCollected[i] = 0;
		}
	}
	void LProposalTable::sampleDistribution(int index)
	{
		if(index < 0 || index > distribution.size())
		{
			return;
		}
		vector<LogDouble> cdf(distribution[index].size());
		cdf[0] = distribution[index].at(0);
		for(unsigned int i=1;i<distribution[index].size()-1;i++)
		{
			cdf[i] = cdf[i-1]+distribution[index].at(i);
		}
		cdf[cdf.size()-1] = LogDouble(1,false);
		double u = LvRandomGenUtil::Instance()->getNormRand();
		LogDouble lu(u,false);
		int nTrues = lower_bound(cdf.begin(),cdf.end(),lu) - cdf.begin();
		sampledValue = nTrues;
		sampledProbability = distribution[index].at(nTrues);
		//update the appropriate counter
		adaptiveMLECounter[index].at(nTrues)++;
		totalSamplesCollected[index]++;
	}

	void LProposalTable::printDistribution()
	{
		for(unsigned int i=0;i<distribution.size();i++)
		{
			cout<<i<<"th distribution = ";
			for(unsigned int j=0;j<distribution[i].size();j++)
			{
				cout<<distribution[i].at(j)<<" ";
			}
			cout<<endl;
		}
		cout<<endl;
	}

	int LProposalTable::getTableSize(int index)
	{
		return distribution[index].size();
	}

	LogDouble LProposalTable::getDistributionValue(int index,int numTrue)
	{
		return distribution[index].at(numTrue);
	}

	LogDouble LProposalTable::getNormConstant(int index)
	{
		return normConstants[index];
	}

	LProposalDistributionElement::LProposalDistributionElement(Atom* atom_,vector<LProposalDistributionElement*> parents_,
		vector<bool> isolatedTerms_,int id_):parents(parents_),isolatedTerms(isolatedTerms_),id(id_),nTrueValues(0)
	{
		atom = LvrMLN::create_new_atom(atom_);
		initializeProposalTables();
	}

	LProposalDistributionElement::~LProposalDistributionElement()
	{
		delete atom;
		cleanup(proposalTables);
		parents.clear();
	}
	void LProposalDistributionElement::initializeProposalTables()
	{
		int nonisolatedSize = 1;
		int isolatedSize = 1;
		int totalSize = 1;
		for(unsigned int i=0;i<atom->terms.size();i++)
		{
			if(!isolatedTerms[i])
			{
				nonisolatedSize *= atom->terms[i]->domain.size();
			}
			else
			{
				isolatedSize *= atom->terms[i]->domain.size();
			}
			totalSize *= atom->terms[i]->domain.size();
		}
		int parentsSize = 1;
		for(unsigned int i=0;i<parents.size();i++)
		{
			parentsSize *= (parents[i]->atom->getNumberOfGroundings()+1);
		}
		if(isolatedSize==1)
		{
			//exactly one table
			LProposalTable* table = new LProposalTable(parentsSize,atom->getNumberOfGroundings());
			proposalTables.push_back(table);
		}
		else
		{
			//multiple tables, one for each grounding to nonisolated terms
			for(unsigned int i=0;i<nonisolatedSize;i++)
			{
				LProposalTable* table = new LProposalTable(parentsSize,isolatedSize);
				proposalTables.push_back(table);
			}
		}
		
	}

	void LProposalDistributionElement::sample(LogDouble& prob,vector<int>& sampledValues, int domainSize)
	{
		//print();
		prob = LogDouble(1,false);
		sampledValues.clear();
		sampledValues.resize(proposalTables.size());
		
		int index = 0;
		//from the current parent's value get the index ofdistribution to sample
		if(parents.size()==0)
		{
			index=0;
		}
		else
		{
			vector<int> parentValues;
			vector<int> parentSizes;
			for(int i=0;i<parents.size();i++)
			{
				parentValues.push_back(parents[i]->nTrueValues);
				parentSizes.push_back(parents[i]->atom->getNumberOfGroundings());
			}
			
			//decode
			index = parentValues[parentValues.size()-1];
			for(int j=parentSizes.size()-2;j>=0;j--)
			{
				for(unsigned int k=j+1;k<parentSizes.size();k++)
					index += parentValues[j]*parentSizes[k];
			}
			
		}
		int totalValue = 0;
		LogDouble probToUpdate(1,false);
		for(unsigned int i=0;i<proposalTables.size();i++)
		{
			proposalTables[i]->sampleDistribution(index);
			//check for mismatch
			double tableSize = proposalTables[i]->getTableSize(index);
			if(domainSize+1 != tableSize && domainSize > 0)
			{
				if(domainSize+1 < tableSize)
				{
				//sampled for a variant of the atom that has a larger domain than atom being currently sampled
				//need to adjust sample value and prob
				//pick x unique numbers and see how many fall within domainSize
				set<int> nums;
				int nlvTrue = 0;
				while(nums.size()!=proposalTables[i]->sampledValue)
				{
					pair<set<int>::iterator,bool> ret;
					int r = LvRandomGenUtil::Instance()->getRandomPosition(tableSize);
					ret = nums.insert(r);
					if(ret.second)
					{
						if(r < domainSize)
							nlvTrue++;
					}
				}
				//update the probability
				int diff = (tableSize - 1) - domainSize;
				LogDouble nprob = proposalTables[i]->getDistributionValue(index,nlvTrue) + proposalTables[i]->getDistributionValue(index,nlvTrue + diff);
				LogDouble factor;
				for(unsigned int jj=nlvTrue+1;jj<nlvTrue + diff-1;jj++)
					factor = factor + proposalTables[i]->getDistributionValue(index,jj);
				if(!factor.is_zero)
					factor = factor*LogDouble(diff,false);
				else
					factor = LogDouble(diff,false);

				LogDouble nprob_norm =  (factor*nprob)/proposalTables[i]->getNormConstant(index);
				prob *= nprob_norm;
				sampledValues[i] = nlvTrue;
				}
			}
			else
			{
				prob *= proposalTables[i]->sampledProbability;
				sampledValues[i] = proposalTables[i]->sampledValue;
			}
			probToUpdate *= proposalTables[i]->sampledProbability;
			totalValue += proposalTables[i]->sampledValue;
		}
		
		this->probability = probToUpdate;
		this->nTrueValues = totalValue;
	}

	void LProposalDistributionElement::print()
	{
		atom->print(false);
		cout<<"[";
		for(unsigned int i=0;i<parents.size();i++)
		{
			parents[i]->atom->print(false);
			cout<<",";
		}
		cout<<"]"<<endl;
		for(unsigned int i=0;i<proposalTables.size();i++)
		{
			cout<<"Proposal Table "<<i<<endl;
			proposalTables[i]->printDistribution();
		}
	}
	
	void LProposalDistributionElement::updateProposals(double learningRate)
	{
		for(unsigned int i=0;i<proposalTables.size();i++)
			proposalTables[i]->updateDistribution(learningRate);
	}



	LProposalDistribution::LProposalDistribution(ESamplingMode samplingMode_):counter(0),samplingMode(samplingMode_)
	{
		hashedDistributions = new LvrAtomHashTemplate<LProposalDistributionElement*>();
	}
	LProposalDistribution::~LProposalDistribution()
	{
		delete hashedDistributions;
	}
	LProposalDistributionElement* LProposalDistribution::insertElement(Atom* atom, vector<Atom*> parents,vector<bool> isolatedTerms)
	{
		//find all parents
		
		vector<LProposalDistributionElement*> lvParents;
		for(unsigned int i=0;i<parents.size();i++)
		{
			LProposalDistributionElement* elem = hashedDistributions->getValue(parents[i]);
			if(elem!=NULL)
				lvParents.push_back(elem);
		}
		//create a proposal element
		LProposalDistributionElement* nElement = new LProposalDistributionElement(atom,lvParents,isolatedTerms,hashedDistributions->size());
		if(samplingMode==EINFORMED)
		{
			//useparent id as key
			hashedDistributions->insert(atom->symbol->parentId,nElement);
		}
		else
		{
			//use atom as key
			hashedDistributions->insert(atom,nElement);
		}
		return nElement;
	}
	
	void LProposalDistribution::updateDistributions(double learningRate)
	{
		cout<<"Updating Distributions..."<<endl;
		vector<LProposalDistributionElement*> allElems = hashedDistributions->getAllData();
		for(vector<LProposalDistributionElement*>::iterator it=allElems.begin();it!=allElems.end();it++)
		{
			(*it)->updateProposals(learningRate);
		}
	}
	void LProposalDistribution::print()
	{
		cout<<"Proposal Structure"<<endl;
		vector<LProposalDistributionElement*> allElems = hashedDistributions->getAllData();
		for(vector<LProposalDistributionElement*>::iterator it=allElems.begin();it!=allElems.end();it++)
		{
			(*it)->print();
			cout<<" ";
		}
		cout<<endl;
	}

	LProposalDistributionElement* LProposalDistribution::findElement(Atom* atom)
	{
		if(samplingMode==EINFORMED)
		{
			LProposalDistributionElement* dist;
			hashedDistributions->getValue(atom->symbol->parentId,dist);
			return dist;
		}
		else
		{
			return hashedDistributions->getValue(atom);
		}
	}

	void LProposalDistribution::dumpDistributionToFile()
	{
		vector<LProposalDistributionElement*> allElems = hashedDistributions->getAllData();
		vector<int> ids;
		vector<Atom*> atoms;
		vector<vector<Atom*> > parents;
		vector<vector<bool> > isolatedTerms;
		for(vector<LProposalDistributionElement*>::iterator it=allElems.begin();it!=allElems.end();it++)
		{
			atoms.push_back((*it)->atom);
			vector<Atom*> tmpParents;
			vector<bool> tmpISTerms;
			for(unsigned int i=0;i<(*it)->parents.size();i++)
				tmpParents.push_back((*it)->parents[i]->atom);
			for(unsigned int i=0;i<(*it)->isolatedTerms.size();i++)
				tmpISTerms.push_back((*it)->isolatedTerms[i]);
			parents.push_back(tmpParents);
			isolatedTerms.push_back(tmpISTerms);
			ids.push_back((*it)->id);
		}

		LFileDump::dumpProposalToFile(atoms,parents,ids);
		LFileDump::dumpProposalITTermsFile(isolatedTerms);
	}


	void LProposalDistribution::readDistributionFromDumpFile()
	{
		vector<Atom*> atoms;
		vector<vector<Atom*> > parents;
		vector<vector<bool> >  isolatedTerms;
		vector<int> ids;
		LFileDump::readDumpedProposalFile(atoms,parents,ids);
		LFileDump::readDumpedProposalIsolatedTermsFile(isolatedTerms);
		//set it in the right order
		vector<Atom*> atoms_n(atoms.size());
		vector<vector<Atom*> > parents_n(parents.size());
		vector<vector<bool> >  isolatedTerms_n(isolatedTerms.size());
		for(unsigned int i=0;i<ids.size();i++)
		{
			atoms_n[ids[i]] = atoms[i];
			parents_n[ids[i]]=parents[i];
			isolatedTerms_n[ids[i]] = isolatedTerms[i];
		}
		for(unsigned int i=0;i<atoms_n.size();i++)
		{
			insertElement(atoms_n[i],parents_n[i],isolatedTerms_n[i]);
		}
	}
