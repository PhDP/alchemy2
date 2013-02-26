#include "cache.h"
#include "cleanuputils.h"

LCacheEntry::LCacheEntry(vector<WClause*> CNF_,LogDouble wCount_):wCount(wCount_),deleteFlag(false),hits(0)
{
	//store a copy
	LvrMLN::copyAllClauses(CNF_,CNF);
}

LCacheEntry::~LCacheEntry()
{
	cleanup(CNF);
}


LCache::LCache(LvrMLN& mln_,int cacheSize_):mln(mln_),currentIndex(0),successfulHits(0),misses(0),hits(0),numEntries(0),mostRecentCounter(0)
{
	cout<<"Initializing cache with size = "<<cacheSize_<<endl;
	cacheEntries.resize(cacheSize_);
	for(unsigned int i=0;i<cacheEntries.size();i++)
	{
		cacheEntries[i] = NULL;
	}
	unifier = new LUnifier(mln);
	mostRecentList.resize(MOSTRECENTSIZE);
	cacheIndex.resize(MAXCNFSIZE);
}
LCache::~LCache()
{
	cleanup(cacheEntries);
	cacheIndex.clear();
	delete unifier;
}

bool LCache::addToCache(vector<WClause*> CNF, LogDouble wCount)
{
	if(CNF.size() > MAXCNFSIZE || CNF.size() == 0)
		return false;
	if(!isFull())
	{
	for(unsigned int i=0;i<cacheEntries.size();i++)
	{
		if(cacheEntries[i] == NULL)
		{
			LCacheEntry* lc = new LCacheEntry(CNF,wCount);
			cacheEntries[i] = lc;
			numEntries++;
			if(mostRecentCounter >= MOSTRECENTSIZE)
				mostRecentCounter = 0;
			mostRecentList[mostRecentCounter++] = i;
			if(isFull())
			{
				misses = 0;
				hits = 0;
			}
			cacheIndex[CNF.size()-1].push_back(cacheEntries[i]);
			break;
		}
	}
	}

	return true;
}

bool LCache::findInCache(vector<WClause*> CNF, LogDouble& val)
{
	if(CNF.size() > MAXCNFSIZE || CNF.size() == 0)
		return false;
	//for(unsigned int i=0;i<currentIndex;i++)
	vector<LCacheEntry*> tempCache = cacheIndex[CNF.size() - 1];
	//for(unsigned int i=0;i<cacheEntries.size();i++)
	for(unsigned int i=0;i<tempCache.size();i++)
	{
		if(tempCache[i])
		{
			if(unifier->CNFUnify(tempCache[i]->CNF,CNF))
			{
				val = tempCache[i]->wCount;
				tempCache[i]->hits++;
				successfulHits++;
#ifdef __DEBUG_PRINT__				
				cout<<"Cache Hit "<<i<<endl;
				for(unsigned int j=0;j<tempCache[i]->CNF.size();j++)
				{
					if(tempCache[i]->CNF[j]->satisfied)
						cout<<"Satisfied V ";
					tempCache[i]->CNF[j]->print();
				}
				cout<<"******************"<<endl;
				for(unsigned int j=0;j<CNF.size();j++)
				{
					if(CNF[j]->satisfied)
						cout<<"Satisfied V ";
					CNF[j]->print();
				}
				cout<<"sucesssful hit end"<<endl;
#endif
				hits++;
				return true;
			}
		}
	}
	misses++;
	if(isFull() && (double)misses/(double)(misses+hits) > 0.5)
	{
		cleanCache(0);
		hits = 0;
		misses = 0;
	}
	return false;
}

void LCache::cleanCache(int numItemsToRemove)
{
	//cout<<"clean start"<<numEntries<<endl;
	for(unsigned int i=0;i<cacheEntries.size();i++)
	{
		if(cacheEntries[i])
		{
			if(cacheEntries[i]->hits == 0)
			{
				//set to NULL in the Index
				int index = cacheEntries[i]->CNF.size() - 1;
				for(unsigned int j=0;j<cacheIndex[index].size();j++)
				{
					if(cacheEntries[i] == cacheIndex[index].at(j))
					{
						cacheIndex[index].at(j) = NULL;
						break;
					}
				}
				delete cacheEntries[i];
				cacheEntries[i]=NULL;
				numEntries--;
			}
		}
	}
	if(numEntries == cacheEntries.size())//still full
	{
		//delete everything other than most recent list
		for(unsigned int i=0;i<cacheEntries.size();i++)
		{
			bool found = false;
			for(unsigned int j=0;j<mostRecentList.size();j++)
			{
				if(i==mostRecentList[j])
				{
					found = true;
					break;
				}
			}
			if(!found)
			{
				//set to NULL in the Index
				int index = cacheEntries[i]->CNF.size() - 1;
				for(unsigned int j=0;j<cacheIndex[index].size();j++)
				{
					if(cacheEntries[i] == cacheIndex[index].at(j))
					{
						cacheIndex[index].at(j) = NULL;
						break;
					}
				}				delete cacheEntries[i];
				cacheEntries[i]=NULL;
				numEntries--;
			}
		}
	}
#ifdef __DEBUG_PRINT__
	cout<<"cache cleaned, numentries left="<<numEntries<<endl;
#endif
}

void LCache::printCache()
{
	for(unsigned int i=0;i<cacheEntries.size();i++)
	{
		if(cacheEntries[i])
		{
		for(unsigned int j=0;j<cacheEntries[i]->CNF.size();j++)
		{
			if(cacheEntries[i]->CNF[j]->satisfied)
				cout<<"Satisfied V ";
			cacheEntries[i]->CNF[j]->print();
		}
		cout<<cacheEntries[i]->wCount<<endl;
		cout<<"*************"<<endl;
		}
	}
	
}