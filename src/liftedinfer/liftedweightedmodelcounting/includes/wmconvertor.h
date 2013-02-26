#ifndef __LWM_CONVERTOR__
#define __LWM_CONVERTOR__
#include "lvrmln.h"




struct AtomLocation
{
	int clauseIndex;
	int atomIndex;
	int termIndex;
	AtomLocation(int clauseIndex_,int atomIndex_,int termIndex_):clauseIndex(clauseIndex_),atomIndex(atomIndex_),termIndex(termIndex_){}
	bool isEqual(AtomLocation* at)
	{
		if(at->atomIndex == atomIndex && at->clauseIndex == clauseIndex && at->termIndex == termIndex)
			return true;
		else
			return false;
	}
};


//Convert a PKB (weights on Formulae) into one with weights on the literals
struct LWMConvertor
{
	LvrMLN& mln;
	LWMConvertor(LvrMLN& mln_):mln(mln_)
	{
		conversionFactor = LogDouble(1,false);
	}
	void convertMLN();
	void updateMLN(int formulaIndex,int symbolIndex,vector<WClause*>& newClausestoAdd);
	//vector<WClause*> generateAtomCombinations(int index,int startIndex, LvrMLN& mln);
	vector<WClause*> generateAtomCombinations(int index,int startIndex,vector<vector<AtomLocation*> >& locations);
	void getNewAtomTerms(int formulaIndex,vector<LvrTerm*>& terms,vector<vector<AtomLocation*> >& locationIndex);
	void appendNewPredicate1(vector<LvrTerm*> newTerms, vector<vector<AtomLocation*> > termLocs,vector<WClause*>& clauses,PredicateSymbol* ps);
	void appendNewPredicate(vector<LvrTerm*> newTerms, vector<vector<AtomLocation*> > termLocs, 
	vector<vector<AtomLocation*> > atomLocs,vector<WClause*>& clauses,PredicateSymbol* ps);
	int toSymbolIndex(int predId);
	//int getMatchedLocation(vector<AtomLocation*> aLoc, vector<AtomLocation*> termLocs);
	//AtomLocation* getMatchedLocation(vector<AtomLocation*> storedTermLocs, vector<AtomLocation*> newTermLocs);
	AtomLocation* getMatchedLocation(vector<AtomLocation*> atomLocs, vector<AtomLocation*> termLocs,
	WClause* clause);
	LogDouble conversionFactor;
};
#endif