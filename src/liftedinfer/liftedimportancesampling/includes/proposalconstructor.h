#ifndef __LPROPOSAL_CONSTRUCTOR
#define __LPROPOSAL_CONSTRUCTOR
#include "lvrmln.h"
#include "proposalstructure.h"
#include "cleanuputils.h"
#include "normalizer.h"
#include "clusterutil.h"
#include "filedump.h"
#include "ptpsearch.h"
#include "heuristics.h"
#include "rulesutil.h"
#include "lisapproxinference.h"
#include "lvrnormpropagation.h"
#include "lvparams.h"


struct LProposalConstructor
{
LProposalConstructor(LvrMLN& mln_);
~LProposalConstructor();
void getParents(vector<WClause*> clauses,Atom* atom, vector<Atom*> potentialParents,vector<Atom*>& parents);
void constructProposal(vector<WClause*>& clauses,vector<Atom*>& potentialParents);
void constructProposalV1(vector<WClause*>& CNF,vector<Atom*>& potentialParents);
void print(){}
void startConstruction(LvrParams* params);
void startMARInference(LvrParams* params);
void startPartitionFunction(LvrParams* params);
void readDistributionFromDumpFile();
void dumpDistributionToDumpFile();
private:
LProposalDistribution* lpd;
LvrMLN& mln;
LClusterUtil* lc;
LISApproxInference* lsearch;
LvrNormPropagation* lpg;
LvrSingletonNormPropagation* lvrSingletonNormPropagation;
};

#endif
