#include "ptpsampler.h"
#include <fstream>
using namespace std;
#include "lvrpermutations.h"
#include "queryupdater.h"


LPTPSampler::LPTPSampler(LvrMLN& mln)
{
	ls= new LBlockExactInference(mln);
}

LPTPSampler::~LPTPSampler()
{
	delete ls;
}


int LPTPSampler::runMockPTP(vector<WClause*>& CNF,LVRCluster* lvrCluster,int& lvpCost,int remainingCost)
{
	return ls->doMockInferenceOnCluster(CNF,lvpCost,remainingCost);
}

int LPTPSampler::runPTP(vector<WClause*>& CNF,LVRCluster* lvrCluster)
{
	if(lvrCluster->isUniConstantCluster())
	{
		doLightWeightSampling(lvrCluster,CNF);
		return 0;
	}
	int ret = ls->doExactInferenceOnCluster(CNF);
	if( ret <= 0)
		return ret;
	lvrCluster->resetAllAssignments();	
	for(unsigned int i=0;i<ls->getPTPTreeRef().sampledAssignments.atoms.size();i++)
	{
		int ind = lvrCluster->getEquivalentElementIndex(ls->getPTPTreeRef().sampledAssignments.atoms[i]);
		if(ind==-1)
			continue;
		for(unsigned int j=0;j<lvrCluster->lAssignments[ind].size();j++)
		{
			bool match = true;
			for(unsigned int k=0;k<lvrCluster->sharedGroundedPositions[ind].size();k++)
			{
				if(ls->getPTPTreeRef().sampledAssignments.atoms[i]->terms[lvrCluster->sharedGroundedPositions[ind].at(k)]->domain[0] !=
					lvrCluster->lAssignments[ind].at(j)->grounding[lvrCluster->sharedGroundedPositions[ind].at(k)])
				{
					match = false;
					break;
				}
			}
			if(match)
			{
				if(lvrCluster->lAssignments[ind].at(j)->nTrueValues == -1)
				{
					lvrCluster->lAssignments[ind].at(j)->nTrueValues = ls->getPTPTreeRef().sampledAssignments.assignment[i];
					lvrCluster->lAssignments[ind].at(j)->nFalseValues = lvrCluster->lAssignments[ind].at(j)->nFalseValues
						- ls->getPTPTreeRef().sampledAssignments.assignment[i];
				}
				else
				{
					lvrCluster->lAssignments[ind].at(j)->nTrueValues += ls->getPTPTreeRef().sampledAssignments.assignment[i];
					lvrCluster->lAssignments[ind].at(j)->nFalseValues = lvrCluster->lAssignments[ind].at(j)->nFalseValues
						- ls->getPTPTreeRef().sampledAssignments.assignment[i];
				}
				lvrCluster->totalTrues[ind] += ls->getPTPTreeRef().sampledAssignments.assignment[i];
			}
		}
		//check if not burning in,only then update query
		if(!burnIn)
		{
			LvrQueryUpdater::Instance()->updateQueryValuesLVGibbs(ls->getPTPTreeRef().sampledAssignments.atoms[i],ls->getPTPTreeRef().sampledAssignments.assignment[i]);
		}
	}
	lvrCluster->setAllDontCareAssignments();
	return 0;
}



//if we know our cluster is a single propositional element, we can do fast sampling
void LPTPSampler::doLightWeightSampling(LVRCluster* lvrCluster, vector<WClause*> clauses)
{
	int atomId = lvrCluster->elements[0]->symbol->id;
	LogDouble posWt(1,false);
	LogDouble negWt(1,false);
	bool dontcare=true;
	for(unsigned int i=0;i<clauses.size();i++)
	{
		if(clauses[i]->satisfied)
		{
			double tmp = exp(clauses[i]->weight.value);
			posWt *= LogDouble(tmp,true);
			negWt *= LogDouble(tmp,true);
			delete clauses[i];
			continue;
		}
		bool posSign=false;
		bool negSign = false;
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			if(clauses[i]->sign[j])
				negSign = true;
			else
				posSign = true;
		}
		if(negSign && posSign)
		{
			double tmp = exp(clauses[i]->weight.value);
			posWt *= LogDouble(tmp,true);
			negWt *= LogDouble(tmp,true);
			dontcare = false;
			
		}
		else if(negSign)
		{
			double tmp = exp(clauses[i]->weight.value);
			negWt *= LogDouble(tmp,true);
			dontcare = false;
			
		}
		else if(posSign)
		{
			double tmp = exp(clauses[i]->weight.value);
			posWt *= LogDouble(tmp,true);
			dontcare = false;
		}
		delete clauses[i];
	}
	int assignment = 0;
	if(!dontcare)
	{
		LogDouble posWt1 = posWt/(posWt + negWt);
		LogDouble negWt1 = negWt/(posWt + negWt);
		double r = LvRandomGenUtil::Instance()->getNormRand();
		LogDouble lr(r,false);
		//sample the assignment		
		if(lr > negWt1)
			assignment = 1;
	}
	else
	{
		double r = LvRandomGenUtil::Instance()->getNormRand();
		//sample the assignment
		if(r > 0.5)
			assignment = 1;
	}
	lvrCluster->lAssignments[0].at(0)->nTrueValues = assignment;
	lvrCluster->lAssignments[0].at(0)->nFalseValues = 1 - assignment;
	if(!burnIn)
	{
		LvrQueryUpdater::Instance()->updateQueryValuesLVGibbs(lvrCluster->elements[0],assignment);
	}

	clauses.clear();
}