#ifndef __LIFTEDALGS_HANDLER
#define __LIFTEDALGS_HANDLER
#include "lvrmln.h"
#include "lvparams.h"
#include "clustercreator.h"
#include "gibbsprocesshandler.h"
#include "proposalconstructor.h"

struct LiftedAlgsHandler
{
	LiftedAlgsHandler(LvrMLN& mln_);
	~LiftedAlgsHandler();
	//Lifted Blocked Gibbs Sampler
	void LBGApproxMarginals(LvrParams* params);
	//Lifted Importance Sampler
	void LISApproxPartition(LvrParams* params);
	void LISApproxMarginals(LvrParams* params);
	//Lifted Weighted Model Counter (PTP)
	void WMPTPZApprox(LvrParams* params);
	void WMPTPZExact(LvrParams* params);
	void buildProposal(LvrParams* params);
private:
	LvrMLN& mln;
	LPTPSearch* ptpSearch;
	LClusterCreator* clusterCreator;
	LProposalConstructor* proposalConstructor;
	GibbsProcessHandler* gibbsProcessHandler;
};

#endif
