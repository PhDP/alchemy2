#ifndef TESTFORMULAINFER__
#define TESTFORMULAINFER__
#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;
#include "lvrmln.h"
#define NUMTESTS 6
struct TestFormulaInfer
{
	LvrMLN mln;
	int numtests;
	string testfolder;
	TestFormulaInfer(string testfolder_ = "..\\autotestmlnfiles\\",int numtests_ = NUMTESTS):numtests(numtests_),
	testfolder(testfolder_){}
	void runFormulaInferenceTests();
	~TestFormulaInfer(){}
private:
	void doMLNTestExpWeighted(int index,string mlnfile,string dbfile,bool withEvidence=false);
	void readGroundTruth(int id,vector<vector<vector<int> > >& trueIndexes, vector<int>& evidenceIndexes,
		vector<double>& weights,int& size);
	LogDouble manualEstimate(int id, bool withEvidence=false,bool literalWeights=false);
	bool getNextTruthAssignment(vector<int>& truthValues);
	bool isEvidenceSat(vector<int> truthValue,vector<int> evidence);
	void doMLNTestPTPWeighted(int index,string mlnfile,string dbfile,bool withEvidence=false);
};
#endif
