#ifndef __LMPARSER
#define __LMPARSER
#include "lvrmln.h"

#define DOMAINSTART "#domains"
#define PREDICATESTART "#predicates"
#define FORMULASTART "#formulas"
#define LEFTPRNTH "("
#define RIGHTPRNTH ")"
#define NOTOPERATOR "!"
#define ANDOPERATOR "^"
#define OROPERATOR "|"
#define LEFTFLOWER "{"
#define RIGHTFLOWER "}"
#define WEIGHTSEPARATOR "::"
#define COMMASEPARATOR ","
#define EQUALSTO "="




enum ParsingState
{
	Domains,
	Predicates,
	Formulas
};
struct PDomain
{
	string name;
	vector<string> values;
	PDomain(string name_,vector<string> values_):name(name_),values(values_){}
};

struct LMParser
{
	LvrMLN& mln;
	int predicateId;
	//set of domains with values in string format
	vector<PDomain*> domainList;
	//map key:predicate Id, value:for each of its terms, Index into the domainList vector
	map<int,vector<int> > predicateDomainMap;
	void parseInputMLNFile(string filename);
	string convert_atom_to_string(Atom* atom);
	LvrTerm* create_new_term(int domainSize);
	WClause* create_new_clause(vector<int> predicateSymbolIndex,vector<bool> sign,
						   vector<vector<LvrTerm*> > iTermsList);
	void parseCNFString(string formula,fstream& filestr);
	void parsePredicateString(string line);
	void parseDomainString(string line,fstream& filestr);
	void parseDB(string filename);
	void getAtomDomainIndex(int predicateId,vector<int>& domainIndex);
	bool checkTermsValidity(int predicateId,vector<string> terms, vector<int>& matchedIndexList);
	bool isTermConstant(string term);
	LMParser(LvrMLN& mln_):mln(mln_),predicateId(0){}
	~LMParser();
};

#endif
