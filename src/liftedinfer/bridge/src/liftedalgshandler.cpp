#include "liftedalgshandler.h"
#include "wmconvertor.h"


LiftedAlgsHandler::LiftedAlgsHandler(LvrMLN& mln_):mln(mln_)
{
	ptpSearch = new LPTPSearch(mln);
	clusterCreator = new LClusterCreator(&mln);
	proposalConstructor = new LProposalConstructor(mln);
	gibbsProcessHandler = new GibbsProcessHandler(mln);
}

LiftedAlgsHandler::~LiftedAlgsHandler()
{
	delete ptpSearch;
	delete clusterCreator;
	delete proposalConstructor;
	delete gibbsProcessHandler;
}

void LiftedAlgsHandler::LBGApproxMarginals(LvrParams* params)
{
	gibbsProcessHandler->startLBGForMar(params);
}

void LiftedAlgsHandler::LISApproxPartition(LvrParams* params)
{
	/*if(!params->isWeightLearning)
	{
		if(params->samplingMode == EINFORMED)
		{
			buildProposal(params);
		}
	}
	*/
	proposalConstructor->startPartitionFunction(params);
}

void LiftedAlgsHandler::LISApproxMarginals(LvrParams* params)
{
	/*
	if(!params->isWeightLearning)
	{
		if(params->samplingMode == EINFORMED)
		{
			buildProposal(params);
		}
	}
	*/
	proposalConstructor->startMARInference(params);
}

void LiftedAlgsHandler::buildProposal(LvrParams* params)
{
	proposalConstructor->startConstruction(params);
	cout<<"done construction"<<endl;
}

void LiftedAlgsHandler::WMPTPZApprox(LvrParams* params)
{
	//convert to WM PTP
	cout<<"Pre-processing..Converting to Weighted Model Counting..."<<endl;
	LWMConvertor* lw = new LWMConvertor(mln);
	lw->convertMLN();
	delete lw;		
	cout<<"WM Converted MLN"<<endl;
	for(int i=0;i<mln.clauses.size();i++)
		mln.clauses[i]->print();
	cout<<endl;
	LFileDump::dumpMLNToFile(&mln);
	cout<<"done writing dump files ("<<MLNDUMPFILENAME<<","<<SYMBOLSDUMPFILENAME<<")"<<endl;
	ptpSearch->startApproxWeightedModelCounting(params);
}


void LiftedAlgsHandler::WMPTPZExact(LvrParams* params)
{
	cout<<"Pre-processing..Converting to Weighted Model Counting..."<<endl;
	//convert to WM PTP
	LWMConvertor* lw = new LWMConvertor(mln);
	lw->convertMLN();
	delete lw;		
	cout<<"WM Converted MLN"<<endl;
	for(int i=0;i<mln.clauses.size();i++)
		mln.clauses[i]->print();
	cout<<endl;
	LFileDump::dumpMLNToFile(&mln);
	cout<<"done writing dump files ("<<MLNDUMPFILENAME<<","<<SYMBOLSDUMPFILENAME<<")"<<endl;
	cout<<"Exact Z = ";
	ptpSearch->startExactWeightedModelCounting(params).printValue();
	cout<<endl;
}


