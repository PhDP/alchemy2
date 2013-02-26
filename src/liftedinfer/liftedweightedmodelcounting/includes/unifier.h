#ifndef LUNIFIER__
#define LUNIFIER__
#include "lvrmln.h"
struct LUnifier
{
	bool CNFUnify(vector<WClause*> CNF1, vector<WClause*> CNF2);
	LUnifier(LvrMLN& mln_):mln(mln_){}
	LvrMLN& mln;
private:
	void getTermRelations(WClause* clause,vector<LvrTerm*>& terms, vector<vector<int> >& positions);
	bool clauseUnify(WClause* c1, WClause* c2,map<int,vector<int> >& originalIdConsts);
};
#endif