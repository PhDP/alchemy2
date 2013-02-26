#ifndef PROPOSITIONAL_RESOLUTION__
#define PROPOSITIONAL_RESOLUTION__
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;

#include "lvrmln.h"
struct PropositionalResolution
{	
	static void doPropositionalSplit(vector<WClause*> CNF, Atom* atom,vector<vector<WClause*> >& resolvedCNFList);
	static void doPropositionalSplit(vector<WClause*> CNF, Atom* atom,vector<WClause*>& posClauses);
	static void doPropositionalSplitSelfJoins(vector<WClause*> CNF, Atom* atom,vector<vector<WClause*> >& resolvedCNFList);
};
#endif