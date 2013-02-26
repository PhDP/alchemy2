#ifndef __LEXTENSIONS__
#define __LEXTENSIONS__
#include "lvrmln.h"
#include "cache.h"

struct LExtensions
{
	LogDouble unitPropagation(vector<WClause*>& CNF,bool ptpexactmar=false);
	bool getFromCache(vector<WClause*> CNF, LogDouble& val);
	void storeToCache(vector<WClause*> CNF, LogDouble& val);
	LExtensions(LvrMLN& mln_);
	void setCaching(int cacheSize)
	{
		if(cacheSize > 0)
		{
			if(cacheSize > RAND_MAX)
				cacheSize = CACHESIZE;
			cache = new LCache(mln,cacheSize);
		}
	}
	~LExtensions();
private:
	LCache* cache;
	LvrMLN& mln;
	LogDouble getUnitClauseWeightWMPTP(WClause* clause);
};
#endif
