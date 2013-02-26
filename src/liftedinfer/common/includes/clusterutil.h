#ifndef __LCLusterUtil__
#define __LCLusterUtil__
#include "lvrmln.h"
#include "clusterutil.h"
#include "cleanuputils.h"

struct LClusterUtil
{
	static void seperateDisjointCNF(vector<WClause*> CNF, vector<vector<WClause*> >& disjointCNFList)
	{
		vector<vector<int> > clauseIndexes;
		vector<vector<int> > atomIds;
		int numEmptyClauses=0;
		vector<LogDouble> satWts;
		vector<bool> lvsat;
		for(unsigned int i=0;i<CNF.size();i++)
		{
			if(CNF[i]->atoms.size()==0)
			{
				numEmptyClauses++;
				satWts.push_back(CNF[i]->weight);
				lvsat.push_back(CNF[i]->satisfied);
			continue;
			}
			vector<int> foundPos;
			for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
			{
				for(unsigned int k=0;k<atomIds.size();k++)
				{
					if(find(atomIds[k].begin(),atomIds[k].end(),CNF[i]->atoms[j]->symbol->id)!=atomIds[k].end())
					{
						if(find(foundPos.begin(),foundPos.end(),k)==foundPos.end())
							foundPos.push_back(k);
						break;
					}
				}
			}
			if(foundPos.size()==0)
			{
				vector<int> tmpIds;
				for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
					tmpIds.push_back(CNF[i]->atoms[j]->symbol->id);
				atomIds.push_back(tmpIds);
				vector<int> tmpCind;
				tmpCind.push_back(i);
				clauseIndexes.push_back(tmpCind);
			}
			else
			{
				sort(foundPos.begin(),foundPos.end());
				int itx=0;
				vector<int> tmpAtomIds;
				vector<int> tmpCInds;
				tmpCInds.push_back(i);
				for(unsigned int j=0;j<CNF[i]->atoms.size();j++)
				{
					tmpAtomIds.push_back(CNF[i]->atoms[j]->symbol->id);
				}
				int rmCnt=0;
				for(vector<int>::iterator it=foundPos.begin();it!=foundPos.end();it++)
				{
					int newInd = *it-rmCnt;
					if(newInd < 0)
						break;
					for(unsigned int j=0;j<atomIds[newInd].size();j++)
						tmpAtomIds.push_back(atomIds[newInd].at(j));
					for(unsigned int j=0;j<clauseIndexes[newInd].size();j++)
						tmpCInds.push_back(clauseIndexes[newInd].at(j));
					atomIds.erase(atomIds.begin()+newInd);
					clauseIndexes.erase(clauseIndexes.begin()+newInd);
					rmCnt++;
				}
				atomIds.push_back(tmpAtomIds);
				clauseIndexes.push_back(tmpCInds);
			}
		}
		for(unsigned int i=0;i<clauseIndexes.size();i++)
		{
			vector<WClause*> tmpClauses;
			for(unsigned int j=0;j<clauseIndexes[i].size();j++)
			{
				tmpClauses.push_back(LvrMLN::create_new_clause(CNF[clauseIndexes[i].at(j)]));
			}
			disjointCNFList.push_back(tmpClauses);
		}
		for(unsigned int i=0;i<numEmptyClauses;i++)
		{
			vector<WClause*> tmpClauses;
			WClause* nClause= new WClause();
			nClause->satisfied = lvsat[i];
			nClause->weight = satWts[i];
			tmpClauses.push_back(nClause);
			disjointCNFList.push_back(tmpClauses);
		}
	}

};
#endif
