#ifndef __LNORMALIZER__
#define __LNORMALIZER__
#include "lvrmln.h"
struct LNormalizer
{
	LvrMLN& mln;
	LNormalizer(LvrMLN& mln_):mln(mln_){}
	bool checkEquivalency(Atom* at1,Atom* at2);
	void normalizeClauses(vector<WClause*>& clauses,bool weightsOnAtoms=false,bool resetIds = false);
	void removeNPSelfJoins(vector<WClause*>& CNF);
	void normalizeClausesWithoutJoins(vector<WClause*>& clauses,bool weightsOnAtoms=false,bool resetIds = false);
	bool compareClauses(WClause* c1, WClause* c2);
	bool compareClauses(vector<WClause*> c1,vector<WClause*> c2);
private:
	void findNPSelfJoinedAtomToGround(vector<WClause*>& CNF,Atom*& atom,vector<bool>& isolatedTerms);
	void mergeIdenticalClauses(vector<WClause*>& CNF);
	void mergeIdenticalClausesPTP(vector<WClause*>& CNF);
	void convertToNormalForm(vector<WClause*>& clauses);
	void toNormalForm(PredicateSymbol* symbol,vector<WClause*>& clauses);
	void renamePredicates(vector<WClause*>& CNF,bool resetIds = false);
};
#endif
