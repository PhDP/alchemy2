#ifndef __LRESOLVER__
#define __LRESOLVER__
#include "lvrmln.h"
struct LResolver
{
	void getMatchingTerms(WClause* clause,Atom* atom, vector<LvrTerm*>& terms);
	void permute(vector<LvrTerm*> terms, vector<vector<int> >& permutedDomainList);
	int getIndex(vector<vector<int> > allGroundings, Atom* groundedAtom);
	void resolveGroundClause(WClause* clause,int predicateId, vector<int> truthVector, 
		vector<vector<int> > allGroundings, WClause& resolvedClause);
	void doFullGrounding(vector<WClause*> CNF,Atom* atom,vector<WClause*>& groundCNF,
		vector<vector<int> >& allGroundings,int& numGroundings);
	void reduceDomains(WClause* clause, Atom* atom, vector<int> truthValues, 
		vector<WClause*>& domainModifiedClauses, int singletonIndex);
	void reduceDomainsSelfJoins(WClause* clause, Atom* atom, vector<int> truthValues, 
		vector<WClause*>& domainModifiedClauses, int singletonIndex);
};
#endif