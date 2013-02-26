#include "propositionalresolution.h"
#include "cleanuputils.h"


void PropositionalResolution::doPropositionalSplit(vector<WClause*> CNF, Atom* atom,vector<vector<WClause*> >& resolvedCNFList)
{
	vector<WClause*> posClauses;
	vector<WClause*> negClauses;
	LvrMLN::copyAllClauses(CNF,posClauses);
	LvrMLN::copyAllClauses(CNF,negClauses);
	resolvedCNFList.resize(2);
	for(unsigned int i=0;i<CNF.size();i++)
	{
		vector<int> matchingIndexes;
		for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
		{
			if(CNF[i]->atoms[j]->symbol->id == atom->symbol->id)
			{
				if(matchingIndexes.size()==0)
				{
					if(CNF[i]->sign[j])
					{
						//negation
						negClauses[i]->satisfied = true;
					}
					else
					{
						posClauses[i]->satisfied = true;
					}
				}
				matchingIndexes.push_back(j);
			}
		}
		if(matchingIndexes.size() > 1)
		{
			cout<<"Error"<<endl;
		}
		int removedCnt = 0;
		for(unsigned int j=0;j<matchingIndexes.size();j++)
		{
			negClauses[i]->removeAtom(matchingIndexes[j]-removedCnt);
			posClauses[i]->removeAtom(matchingIndexes[j]-removedCnt);
			removedCnt++;
		}
	}
	resolvedCNFList[0] = negClauses;
	resolvedCNFList[1] = posClauses;
}

void PropositionalResolution::doPropositionalSplit(vector<WClause*> CNF, Atom* atom,vector<WClause*>& posClauses)
{
	LvrMLN::copyAllClauses(CNF,posClauses);
	for(unsigned int i=0;i<CNF.size();i++)
	{
		bool found=false;
		for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
		{
			if(CNF[i]->atoms[j]->symbol->id == atom->symbol->id)
			{
				if(CNF[i]->sign[j])
				{
					//negation
					CNF[i]->satisfied = true;
				}
				else
				{
					posClauses[i]->satisfied = true;
				}
				CNF[i]->removeAtom(j);
				posClauses[i]->removeAtom(j);
				found=true;
				break;
			}
		}
		if(!found)
			posClauses.push_back(LvrMLN::create_new_clause(CNF[i]));
	}
}


void PropositionalResolution::doPropositionalSplitSelfJoins(vector<WClause*> CNF, Atom* atom,vector<vector<WClause*> >& resolvedCNFList)
{
	vector<WClause*> posClauses;
	vector<WClause*> negClauses;
	//cout<<"**************"<<endl;
	//atom->print(false);
	LvrMLN::copyAllClauses(CNF,posClauses);
	LvrMLN::copyAllClauses(CNF,negClauses);
	resolvedCNFList.resize(2);
	for(unsigned int i=0;i<CNF.size();i++)
	{
		//CNF[i]->print();
		vector<int> matchingIndexes;
		int signCnt = 0;
		for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
		{
			if(CNF[i]->atoms[j]->symbol->id == atom->symbol->id)
			{				
				matchingIndexes.push_back(j);
				if(!CNF[i]->sign[j])
					signCnt++;
			}
		}
		if(matchingIndexes.size()==0)
			continue;
		//if signCnt is 0-- ALL are -ve signs
		//if signCnt = matchingIndexes size, all are +ve signs
		// else opposite signs exist
		if(signCnt == 0)
		{
			negClauses[i]->satisfied = true;
		}
		else if(signCnt == matchingIndexes.size())
		{
			posClauses[i]->satisfied = true;
		}
		else
		{
			negClauses[i]->satisfied = true;
			posClauses[i]->satisfied = true;
		}
		int removedCnt = 0;
		for(unsigned int j=0;j<matchingIndexes.size();j++)
		{
			negClauses[i]->removeAtom(matchingIndexes[j]-removedCnt);
			posClauses[i]->removeAtom(matchingIndexes[j]-removedCnt);
			removedCnt++;
		}

	}
	resolvedCNFList[0] = negClauses;
	resolvedCNFList[1] = posClauses;
	/*cout<<"****"<<endl;
	for(unsigned int i=0;i<negClauses.size();i++)
		negClauses[i]->print();
	cout<<"****"<<endl;
	for(unsigned int i=0;i<posClauses.size();i++)
		posClauses[i]->print();
		*/
}
