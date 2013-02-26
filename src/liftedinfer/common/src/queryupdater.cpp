#include "queryupdater.h"
#include "randomgenutil.h"
#include <assert.h>
#include "fileutils.h"
LvrQueryUpdater* LvrQueryUpdater::m_pInstance = NULL;

LvrQueryUpdater* LvrQueryUpdater::Instance()
{
	//the instance must be created beforehand
	assert(m_pInstance!=NULL);
   return m_pInstance;
}

void LvrQueryUpdater::createInstance(vector<vector<int> > queryIntRep,vector<string> queryStrings, string resultfilename,
	vector<vector<int> >* evidenceIntRep)
{
	m_pInstance = new LvrQueryUpdater(queryIntRep,queryStrings,evidenceIntRep);
	LFileUtils::createOutFileInstance(resultfilename);
}

bool LvrQueryUpdater::isInstanceCreated()
{
	if(m_pInstance!=NULL)
		return true;
	else
		return false;
}

LvrQueryUpdater::LvrQueryUpdater(vector<vector<int> > queryIntRep,vector<string> queryStrings,vector<vector<int> >* evidenceIntRep)
{
	lvrAtomHashTemplate = new LvrAtomHashTemplate<int>();
	lvrAtomHashUpdateFlags = new LvrAtomHashTemplate<bool>();
	queryHashTemplate = new LvrAtomHashTemplate<string>();
	currentSampledValue = new LvrAtomHashTemplate<int>();
	cumulativeWeight = new LvrAtomHashTemplate<LogDouble*>();
	for(unsigned i=0;i<queryIntRep.size();i++)
	{
		LogDouble* ld = new LogDouble(0,false);
		cumulativeWeight->insert(queryIntRep[i],ld);
		lvrAtomHashTemplate->insert(queryIntRep[i],0);
		lvrAtomHashUpdateFlags->insert(queryIntRep[i],false);
		currentSampledValue->insert(queryIntRep[i],0);
		if(queryStrings.size() > 0)
			queryHashTemplate->insert(queryIntRep[i],queryStrings[i]);
		queryNormIds.insert(pair<int,bool>(queryIntRep[i].at(0),0));
	}
	if(evidenceIntRep)
	{
		for(unsigned i=0;i<evidenceIntRep->size();i++)
		{
			lvrAtomHashTemplate->deleteEntry((*evidenceIntRep)[i]);
			lvrAtomHashUpdateFlags->deleteEntry((*evidenceIntRep)[i]);
			currentSampledValue->deleteEntry((*evidenceIntRep)[i]);
		}
	}
	currentProbabilities.resize(lvrAtomHashTemplate->size());
}

void LvrQueryUpdater::updateQueryValues(Atom* atom,int sampledTrueVal)
{
	if(!isNormIdInQuery(atom->symbol->normParentId))
		return;
	Atom* origDomainAtom = NULL;
	vector<bool> changedDomain(atom->terms.size());
	vector<LvrTerm*> termsToPermute;
	//check if this atom represents a decomposed atom
	for(unsigned int i=0;i<atom->terms.size();i++)
	{
		if(atom->terms[i]->origDomain.size() > atom->terms[i]->domain.size())
		{
			//decomposed atom
			if(!origDomainAtom)
				origDomainAtom = LvrMLN::create_new_atom(atom);
			origDomainAtom->terms[i]->domain.clear();
			origDomainAtom->terms[i]->domain = origDomainAtom->terms[i]->origDomain;
			termsToPermute.push_back(origDomainAtom->terms[i]);
			changedDomain[i] = true;
		}
	}
	if(origDomainAtom == NULL)
	{
		doUpdateQueryValues(atom,sampledTrueVal);
		return;
	}
	vector<vector<int> > permutedList;
	LvrPermutations::permuteTerms(termsToPermute,permutedList);
	for(unsigned int i=0;i<permutedList.size();i++)
	{
		int iter=0;
		for(unsigned int j=0;j<changedDomain.size();j++)
		{
			if(!changedDomain[j])
				continue;
			origDomainAtom->terms[j]->domain.clear();
			origDomainAtom->terms[j]->domain.push_back(permutedList[i].at(iter++));
		}
		doUpdateQueryValues(origDomainAtom,sampledTrueVal);
	}
	delete origDomainAtom;
}

void LvrQueryUpdater::updateQueryValues(Atom* atom,vector<bool> isolatedTerms,vector<int> sampledTrueVals)
{
	if(!isNormIdInQuery(atom->symbol->normParentId))
		return;
	Atom* origDomainAtom = NULL;
	vector<bool> changedDomain(atom->terms.size());
	vector<LvrTerm*> termsToPermute;
	//check if this atom represents a decomposed atom
	for(unsigned int i=0;i<atom->terms.size();i++)
	{
		if(isolatedTerms[i])
			continue;
		if(atom->terms[i]->origDomain.size() > atom->terms[i]->domain.size())
		{
			//decomposed atom
			if(!origDomainAtom)
				origDomainAtom = LvrMLN::create_new_atom(atom);
			origDomainAtom->terms[i]->domain.clear();
			origDomainAtom->terms[i]->domain = origDomainAtom->terms[i]->origDomain;
			termsToPermute.push_back(origDomainAtom->terms[i]);
			changedDomain[i] = true;
		}
	}
	if(origDomainAtom == NULL)
	{
		doUpdateQueryValues(atom,isolatedTerms,sampledTrueVals);
		return;
	}
	vector<vector<int> > permutedList;
	LvrPermutations::permuteTerms(termsToPermute,permutedList);
	for(unsigned int i=0;i<permutedList.size();i++)
	{
		int iter=0;
		for(unsigned int j=0;j<changedDomain.size();j++)
		{
			if(!changedDomain[j])
				continue;
			origDomainAtom->terms[j]->domain.clear();
			origDomainAtom->terms[j]->domain.push_back(permutedList[i].at(iter++));
		}
		doUpdateQueryValues(origDomainAtom,isolatedTerms,sampledTrueVals);
	}
	delete origDomainAtom;
}

void LvrQueryUpdater::doUpdateQueryValues(Atom* atom,int sampledTrueVal)
{
	//if(!isNormIdInQuery(atom->symbol->normParentId))
		//return;
	//ground the atom, if not grounded
	if(atom->isConstant())
	{
		vector<int> intRep(atom->terms.size()+1);
		intRep[0] = atom->symbol->normParentId;
		for(unsigned int i=0;i<atom->terms.size();i++)
		intRep[i+1] = atom->terms[i]->domain[0];
		bool value;
		bool found = lvrAtomHashUpdateFlags->getValue(intRep,value);
		if(!found)
			return;
		if(!value)
		{
			bool assignment = false;
			if(sampledTrueVal!=0)
			{
				assignment = true;
				lvrAtomHashTemplate->incrementValue(atom,sampledTrueVal);
			}
			//currentSampledValue->update(atom,assignment);
			currentSampledValue->update(atom,sampledTrueVal);
			lvrAtomHashUpdateFlags->update(atom,true);
		}
	}
	else
	{
		vector<vector<int> > permutedList;
		LvrPermutations::permuteTerms(atom->terms,permutedList);
		//increment random "sampledTrueVal" groundings
		set<int> indexes;
		while(indexes.size()!=sampledTrueVal)
		{
			int ind = LvRandomGenUtil::Instance()->getRandomPosition(permutedList.size());
			indexes.insert(ind);
		}
		for(unsigned int i=0;i<permutedList.size();i++)
		{
			vector<int> intRep(permutedList[i].size()+1);
			intRep[0] = atom->symbol->normParentId;
			for(unsigned int jj=0;jj<permutedList[i].size();jj++)
				intRep[jj+1] = permutedList[i].at(jj);
			bool value;
			bool found = lvrAtomHashUpdateFlags->getValue(intRep,value);
			if(!found)
				continue;
			if(!value)
			{
				//bool assignment = false;
				int assignment = 0;
				if(indexes.find(i)!=indexes.end())
				{
					lvrAtomHashTemplate->incrementValue(intRep,1);
					//assignment = true;
					assignment = 1;
				}
				///currentSampledValue->update(intRep,assignment);
				currentSampledValue->update(intRep,assignment);
				lvrAtomHashUpdateFlags->update(intRep,true);
			}
		}
	}
}

void LvrQueryUpdater::doUpdateQueryValues(Atom* atom,vector<bool> isolatedTerms,vector<int> sampledTrueVals)
{
	//if(!isNormIdInQuery(atom->symbol->normParentId))
		//return;
	//ground the atom, if not grounded
	vector<LvrTerm*> iterms;
	vector<LvrTerm*> niterms;
	for(unsigned int i=0;i<atom->terms.size();i++)
	{
		if(!isolatedTerms[i])
		{
			niterms.push_back(atom->terms[i]);
		}
		else
		{
			iterms.push_back(atom->terms[i]);
		}
	}
	vector<vector<int> > permutedList;
	LvrPermutations::permuteTerms(niterms,permutedList);
	vector<vector<int> > permutedList1;
	LvrPermutations::permuteTerms(iterms,permutedList1);
	vector<int> grounding(atom->terms.size());
	for(unsigned int i=0;i<permutedList.size();i++)
	{
		int iter=0;
		for(unsigned int jj=0;jj<grounding.size();jj++)
		{
			if(!isolatedTerms[jj])
				grounding[jj] = permutedList[i].at(iter);
		}
		//set random groundings to true
		//set<int> indexes;
		vector<int> indexesToSetTrue;
		int countInd = 0;
		while(countInd!=sampledTrueVals[i])
		{
			int ind = LvRandomGenUtil::Instance()->getRandomPosition(permutedList1.size());
			if(find(indexesToSetTrue.begin(),indexesToSetTrue.end(),ind)==indexesToSetTrue.end())
			{
				//indexes.insert(ind);
				indexesToSetTrue.push_back(ind);
				countInd++;
			}
		}
		//for(set<int>::iterator it = indexes.begin();it!=indexes.end();it++)
		for(unsigned int kk=0;kk<permutedList1.size();kk++)
		{
			int iter = 0;
			for(unsigned int jj=0;jj<grounding.size();jj++)
			{
				if(isolatedTerms[jj])
					grounding[jj] = permutedList1[kk].at(iter);
			}
			vector<int> intRep(grounding.size()+1);
			//intRep[0]=atom->symbol->parentId;
			intRep[0]=atom->symbol->normParentId;
			for(unsigned int jj=0;jj<grounding.size();jj++)
				intRep[jj+1] = grounding[jj];
			bool value;
			bool found = lvrAtomHashUpdateFlags->getValue(intRep,value);
			if(!found)
				continue;
			if(!value)
			{
				//bool assignment=false;
				int assignment=0;
				if(find(indexesToSetTrue.begin(),indexesToSetTrue.end(),kk)!=indexesToSetTrue.end())
				{
					lvrAtomHashTemplate->incrementValue(intRep,1);
					//assignment=true;
					assignment = 1;
				}
				currentSampledValue->update(intRep,assignment);
				lvrAtomHashUpdateFlags->update(intRep,true);
			}
		}

		/*
		vector<int> intRep1(grounding.size()+1);
		//intRep[0]=atom->symbol->parentId;
		intRep1[0] = atom->symbol->normParentId;
		for(unsigned int jj=0;jj<grounding.size();jj++)
			intRep1[jj+1] = grounding[jj];
		lvrAtomHashUpdateFlags->update(intRep1,true);
		*/
	}
}

void LvrQueryUpdater::updateQueryValuesLVGibbs(Atom* atom,int sampledTrueVal)
{
	if(!isNormIdInQuery(atom->symbol->normParentId))
		return;
	//ground the atom, if not grounded
	if(atom->isConstant())
	{
		vector<int> intRep(atom->terms.size()+1);
		intRep[0] = atom->symbol->normParentId;
		for(unsigned int i=0;i<atom->terms.size();i++)
		intRep[i+1] = atom->terms[i]->domain[0];
		bool value;
		bool found = lvrAtomHashUpdateFlags->getValue(intRep,value);
		if(!found)
			return;
		if(!value)
		{
			bool assignment = false;
			if(sampledTrueVal!=0)
			{
				assignment = true;
				lvrAtomHashTemplate->incrementValue(atom,sampledTrueVal);
			}
			//currentSampledValue->update(atom,assignment);
			currentSampledValue->update(atom,sampledTrueVal);
			lvrAtomHashUpdateFlags->update(atom,true);
		}
	}
	else
	{
		vector<vector<int> > permutedList;
		LvrPermutations::permuteTerms(atom->terms,permutedList);
		//increment random "sampledTrueVal" groundings
		set<int> indexes;
		while(indexes.size()!=sampledTrueVal)
		{
			int ind = LvRandomGenUtil::Instance()->getRandomPosition(permutedList.size());
			indexes.insert(ind);
		}
		for(unsigned int i=0;i<permutedList.size();i++)
		{
			vector<int> intRep(permutedList[i].size()+1);
			intRep[0] = atom->symbol->normParentId;
			for(unsigned int jj=0;jj<permutedList[i].size();jj++)
				intRep[jj+1] = permutedList[i].at(jj);
			bool value;
			bool found = lvrAtomHashUpdateFlags->getValue(intRep,value);
			if(!found)
				continue;
			if(!value)
			{
				//bool assignment = false;
				int assignment = 0;
				if(indexes.find(i)!=indexes.end())
				{
					lvrAtomHashTemplate->incrementValue(intRep,1);
					//assignment = true;
					assignment = 1;
				}
				///currentSampledValue->update(intRep,assignment);
				currentSampledValue->update(intRep,assignment);
				lvrAtomHashUpdateFlags->update(intRep,true);
			}
		}
	}
}

//called to sample query atoms that were set as don't care in iteration
void LvrQueryUpdater::updateDontCare()
{
	vector<bool> values;
	vector<int> keys;
	lvrAtomHashUpdateFlags->getAllKeyValuePairs(keys,values);
	for(unsigned int i=0;i<values.size();i++)
	{
		if(!values[i])
		{
			bool assignment = false;
			double v = LvRandomGenUtil::Instance()->getNormRand();
			if(v > 0.5)
			{
				assignment = true;
				lvrAtomHashTemplate->incrementIntValue(keys[i]);
			}
			//currentSampledValue->setValue(keys[i],assignment);
			currentSampledValue->setValue(keys[i],-1);
		}
		lvrAtomHashUpdateFlags->setValue(keys[i],false);
	}
}

void LvrQueryUpdater::normalize(int iterations)
{
	vector<int> mleCounts = lvrAtomHashTemplate->getAllData();
	for(unsigned int i=0;i<mleCounts.size();i++)
	{
		currentProbabilities[i] = (double)mleCounts[i]/(double)iterations;
	}
}

void LvrQueryUpdater::writeToFile(int numIterations)
{
	normalize(numIterations);
	vector<string> queryStrings = queryHashTemplate->getAllData();
	LFileUtils::Instance()->updateFile(currentProbabilities,queryStrings);
}

void LvrQueryUpdater::updateAllImportanceWeights(LogDouble currentIterWeight)
{
	vector<LogDouble*> values;
	vector<int> keys;
	cumulativeWeight->getAllKeyValuePairs(keys,values);
	for(unsigned int i=0;i<keys.size();i++)
	{
		//bool sampleValue;
		int sampleValue;
		bool found = currentSampledValue->getValue(keys[i],sampleValue);
		if(!found)
			continue;
		if(sampleValue == 0)
			continue;
		else if(sampleValue==1)
			(*values[i]) = (*values[i]) + currentIterWeight;
		else if(sampleValue == -1)
			(*values[i]) = (*values[i]) + currentIterWeight/LogDouble(2,false);
	}
}

void LvrQueryUpdater::writeToFile(LogDouble totalWeight)
{
	vector<LogDouble*> values;
	vector<int> keys;
	cumulativeWeight->getAllKeyValuePairs(keys,values);
	vector<LogDouble> probabilities(values.size());
	for(unsigned int i=0;i<values.size();i++)
	{
		LogDouble p = (*values[i])/totalWeight;
		/*if(p.is_zero)
			probabilities[i] = 0;
		else
			probabilities[i] = exp(p.value);
			*/
		probabilities[i] = p;
	}

	vector<string> queryStrings = queryHashTemplate->getAllData();
	LFileUtils::Instance()->updateFile(probabilities,queryStrings);
}

int LvrQueryUpdater::getCurrentSampleValue(vector<int> queryintRep)
{
	//bool value;
	int value;
	bool found = currentSampledValue->getValue(queryintRep,value);
	if(!found)
		return -1;
	if(value < 0)
	{
		//dont care atom, sample either 1 or 0 as the current sampled value
		double r=LvRandomGenUtil::Instance()->getNormRand();
		if(r < 0.5)
			value = 0;
		else
			value = 1;
	}
	return value;
}

void LvrQueryUpdater::updateExactQueryWeights(set<int> queryHashes,LogDouble value)
{
	if(value.is_zero)
		return;
	for(set<int>::iterator it = queryHashes.begin();it!=queryHashes.end();it++)
	{
		LogDouble* ld = new LogDouble(0,false);
		cumulativeWeight->getValue((*it),ld);
		(*ld) += value;
		cumulativeWeight->setValue((*it),ld);
	}
}
