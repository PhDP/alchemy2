#ifndef __LPROPOSAL_STRUCTURE
#define __LPROPOSAL_STRUCTURE
#include "lvrmln.h"
#include "lvratomhashtemplate.h"
#include "lvparams.h"

struct LProposalTable
{
	LogDouble sampledProbability;
	int sampledValue;	
	LProposalTable(int parentsSize,int maxTruthValue_);
	void updateDistribution(double learningRate);
	void sampleDistribution(int index);
	void printDistribution();
	int getTableSize(int index);
	LogDouble getDistributionValue(int index,int numTrue);
	LogDouble getNormConstant(int index);
private:
	vector<vector<LogDouble> > distribution;
	vector<vector<int> > adaptiveMLECounter;
	vector<int> totalSamplesCollected;
	int maxTruthValue;
	vector<LogDouble> normConstants;
};

struct LProposalDistributionElement
{
	LProposalDistributionElement(Atom* atom_,vector<LProposalDistributionElement*> parents_,
		vector<bool> isolatedTerms_,int id_);
	~LProposalDistributionElement();
	void initializeProposalTables();
	void sample(LogDouble& prob,vector<int>& sampledValues,int domainSize = -1);
	void print();
	void updateProposals(double learningRate);

	Atom* atom;
	int id;
	vector<LProposalTable*> proposalTables;
	vector<bool> isolatedTerms;
	//do not delete pointers
	vector<LProposalDistributionElement*> parents;
	int nTrueValues;
	LogDouble probability;

};


struct LProposalDistribution
{
	LProposalDistribution(ESamplingMode samplingMode_);
	~LProposalDistribution();
	LProposalDistributionElement* insertElement(Atom* atom, vector<Atom*> parents,vector<bool> isolatedTerms);
	void updateDistributions(double learningRate);
	LProposalDistributionElement* findElement(Atom* atom);
	void dumpDistributionToFile();
	void readDistributionFromDumpFile();
	void print();
private:
	int counter;
	LvrAtomHashTemplate<LProposalDistributionElement*>* hashedDistributions;
	ESamplingMode samplingMode;
};

#endif
