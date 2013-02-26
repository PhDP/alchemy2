#include "splitter.h"
#include "cleanuputils.h"

LSplitter::LSplitter()
{
	sampler = new LSampler();
	resolver = new LResolver();
}

LSplitter::~LSplitter()
{
	delete sampler;
	delete resolver;
}

void LSplitter::doFullGrounding(vector<WClause*> CNF,Atom* atom,vector<WClause*>& groundCNF,
			vector<vector<int> >& allGroundings,int& numGroundings)
{
	resolver->doFullGrounding(CNF,atom,groundCNF,allGroundings,numGroundings);
}

void LSplitter::doCompleteSplitInChunks(vector<WClause*> groundCNF,int atomId, vector<vector<int> > allGroundings, int numGroundings,
	vector<vector<WClause*> >& resolvedCNFList,int startnum, int chunkSize,bool& done)
{

	vector<int> binVector(numGroundings);
	for(unsigned int i=0;i<numGroundings;i++)
	{
		binVector[i] = 0;
	}
	int mask = 1;
	int iter=numGroundings-1;
	while(startnum)
	{
		int x = startnum & mask;
		binVector[iter--] = x;
		startnum = startnum >> 1;

	}
	iter = 0;
	
	while(1)
	{
		vector<WClause*> resolvedCNF;
		for(unsigned int i=0;i<groundCNF.size();i++)
		{
			WClause* resolvedClause=LvrMLN::create_new_clause(groundCNF[i]);
			resolver->resolveGroundClause(groundCNF[i],atomId,binVector,allGroundings,*resolvedClause);
			resolvedCNF.push_back(resolvedClause);
		}
		resolvedCNFList.push_back(resolvedCNF);
		int ind = binVector.size()-1;
		binVector[ind]++;
		while(binVector[ind] == 2)
		{
			binVector[ind]=0;
			if((ind-1) >= 0)
				binVector[ind-1]++;
			else
			{
				done=true;
				break;
			}
			ind--;
		}
		iter++;
		//reached chunk size assignments
		if(iter == chunkSize)
			break;
		//reached the end of all assignments
		if(done)
			break;
	}
	
}

void LSplitter::doLiftedSplit(vector<WClause*> CNF, Atom* atom, vector<vector<WClause*> >& resolvedCNFList,int singletonIndex)
{
	int numberOfGroundings=1;
	for(unsigned int i=0;i<atom->terms.size();i++)
	{
		numberOfGroundings *= atom->terms[i]->domain.size();
	}
	vector<int> truthValues(numberOfGroundings);
	for(unsigned int i=0;i<numberOfGroundings;i++)
	{
		truthValues[i] = 0;
	}
	resolvedCNFList.resize(numberOfGroundings+1);
	for(unsigned int t=0;t<=numberOfGroundings;t++)
	{
		for(unsigned int m=0;m<t;m++)
			truthValues[m]=1;
		vector<WClause*> reducedDomainCNF;
		for(unsigned int i=0;i<CNF.size();i++)
		{
			vector<WClause*> newClauses;
			resolver->reduceDomains(CNF[i],atom,truthValues,newClauses,singletonIndex);
			for(unsigned int k=0;k<newClauses.size();k++)
				reducedDomainCNF.push_back(newClauses[k]);
		}
		resolvedCNFList[t] = reducedDomainCNF;
	}
}



bool LSplitter::doApproxSplit(vector<WClause*> CNF,Atom* atom, int numGroundings, vector<bool> completedGroups, int completedCount,
	LogDouble& probOfSample, int& numTrue, vector<vector<WClause*> >& resolvedList)
{
	int singletonIndex;
	bool singleton = atom->isSingletonAtom(singletonIndex);	
	vector<int> truthValues;
	sampler->doUniformDistSampling(numGroundings, completedGroups, completedCount, probOfSample, numTrue, truthValues);
	if(!singleton)
	{
		//Non singleton atom, do the full grounding and sample
		vector<WClause*> groundCNF;
		vector<vector<int> > allGroundings;
		int numGroundings;
		doFullGrounding(CNF,atom,groundCNF,allGroundings,numGroundings);
		vector<WClause*> resolvedCNF;
		for(unsigned int i=0;i<groundCNF.size();i++)
		{
			WClause* resolvedClause=LvrMLN::create_new_clause(groundCNF[i]);
			resolver->resolveGroundClause(groundCNF[i],atom->symbol->id,truthValues,allGroundings,*resolvedClause);
			resolvedCNF.push_back(resolvedClause);
		}
		resolvedList.push_back(resolvedCNF);
		cleanup(groundCNF);
	}
	else
	{
		//check for self joins on the singleton
		bool selfJoined = false;
		for(unsigned int i=0;i<CNF.size();i++)
		{
			if(CNF[i]->isSelfJoinedOnAtom(atom))
			{
				selfJoined = true;
				break;
			}
		}
		
		if(!selfJoined)
		{
			//do an efficient domain reduction for these clauses based on sampled truth vector
			vector<WClause*> reducedDomainCNF;
			for(unsigned int i=0;i<CNF.size();i++)
			{
				vector<WClause*> newClauses;
				resolver->reduceDomains(CNF[i],atom,truthValues,newClauses,singletonIndex);
				for(unsigned int k=0;k<newClauses.size();k++)
					reducedDomainCNF.push_back(newClauses[k]);
			}
			resolvedList.push_back(reducedDomainCNF);
		}
		else
		{
			//Self joined atom, do the full grounding and sample
			vector<WClause*> groundCNF;
			vector<vector<int> > allGroundings;
			int numGroundings;
			doFullGrounding(CNF,atom,groundCNF,allGroundings,numGroundings);
			vector<WClause*> resolvedCNF;
			for(unsigned int i=0;i<groundCNF.size();i++)
			{
				WClause* resolvedClause=LvrMLN::create_new_clause(groundCNF[i]);
				resolver->resolveGroundClause(groundCNF[i],atom->symbol->id,truthValues,allGroundings,*resolvedClause);
				resolvedCNF.push_back(resolvedClause);
			}
			resolvedList.push_back(resolvedCNF);
			cleanup(groundCNF);
		}
	}
	return false;
}

void LSplitter::doLiftedSplitSelfJoins(vector<WClause*> CNF, Atom* atom, vector<vector<WClause*> >& resolvedCNFList,int singletonIndex)
{
	int numberOfGroundings=1;
	for(unsigned int i=0;i<atom->terms.size();i++)
	{
		numberOfGroundings *= atom->terms[i]->domain.size();
	}
	vector<int> truthValues(numberOfGroundings);
	for(unsigned int i=0;i<numberOfGroundings;i++)
	{
		truthValues[i] = 0;
	}
	resolvedCNFList.resize(numberOfGroundings+1);
	for(unsigned int t=0;t<=numberOfGroundings;t++)
	{
		for(unsigned int m=0;m<t;m++)
			truthValues[m]=1;
		vector<WClause*> reducedDomainCNF;
		for(unsigned int i=0;i<CNF.size();i++)
		{
			vector<WClause*> newClauses;
			resolver->reduceDomainsSelfJoins(CNF[i],atom,truthValues,newClauses,singletonIndex);
			for(unsigned int k=0;k<newClauses.size();k++)
				reducedDomainCNF.push_back(newClauses[k]);
		}
		resolvedCNFList[t] = reducedDomainCNF;
	}
}

