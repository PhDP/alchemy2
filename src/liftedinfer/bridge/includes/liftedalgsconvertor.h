#ifndef __LALGSCONVERTOR
#define __LALGSCONVERTOR
//alchemy formats
#include "mln.h"
#include "variablestate.h"

//lifted algs format
#include "lvrmln.h"
#include "liftedalgshandler.h"
//utility to dump conversions
#include "filedump.h"

#include "lvparams.h"


//conversion utility between the alchemy internal representation and the PTP algorithms internal representation
class LiftedAlgsConvertor
{
public:
	//constructor for weight learning where information needs to be converted back and forth
	//between the alchemy format and the lifted inference format
	LiftedAlgsConvertor(VariableState* state);
	//constructor for marginal estimations where information only needs to be converted
	//from alchemy format to lifted inference format
	LiftedAlgsConvertor(MLN* mln, Domain* domain,GroundPredicateHashArray queries);
	~LiftedAlgsConvertor();

	//functions to interface with the algorithms
	void doLBGMar(LvrParams* params);
	void doLISApproxZ(LvrParams* params);
	void doLISApproxMar(LvrParams* params);
	void doWMPTPApproxZ(LvrParams* params);
	void doWMPTPExactZ(LvrParams* params);

	//Weight learning specific functions
	void updateState();
	void updateMLNWeights(vector<double> weights);
	void processWeightLearningInput(LvrParams* params,bool initialization = false);

private:
	//conversion functions
	//convert from Clause (alchemy format) to WClause (lifted algorithms  format)
	WClause* convertClauseToLifted(Clause* clause, const Domain* domain,
			vector<PredicateSymbol*>& symbolsToAdd,bool& containsGroundedClause);
	//convert from mln (alchemy format) to LvrMLN (lifted algorithms  format)
	void convertInputMLNToLifted(const MLN* mln, const Domain* domain,bool& containsGroundedClause);
	PredicateSymbol* createPredicateSymbol(Predicate* pred);
	void addEvidence(vector<vector<int> >& evidenceIntRep);
	void processSelfJoins();
	void doConversion(vector<vector<int> >& queriesIntRep, vector<vector<int> > & evidenceIntRep,
		vector<string>& queryStrings, bool removeSelfJoins=false);


	LvrMLN* lvrMln;
	LiftedAlgsHandler* liftedAlgsHandler;

	//not owned
	VariableState* state_;
	const MLN* mln_;
	const Domain* domain_;
	GroundPredicateHashArray queries_;
};
#endif

