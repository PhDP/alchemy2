#ifndef __HEURISTICS__
#define __HEURISTICS__
#include "decomposer.h"

struct LHeuristics
{
	LDecomposer& decomposer;
	LvrMLN& mln;
	LHeuristics(LDecomposer& decomposer_,LvrMLN& mln_):decomposer(decomposer_),mln(mln_){}
	Atom* decomposition_lookahead(vector<WClause*> clauses);
	Atom* getAtomToSplit(vector<WClause*> CNF,bool ptpexactmar=false);
	static Atom* getAtomToSplitH2(vector<WClause*> CNF);
	static int proposalCouplingScore(vector<WClause*> origClauses,Atom* a1, Atom* a2);
};
#endif