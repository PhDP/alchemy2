#ifndef TESTMODELCOUNT__
#define TESTMODELCOUNT__
#include <vector>
#include <iostream>
#include <set>
#include <map>
using namespace std;
#include "lvrmln.h"
#define NUMAUTOTESTFILES 6
struct TestModelCount
{
	LvrMLN mln;
	TestModelCount(){}
	void runAllWMCTests();
	void createMLNFromId(int testId);
	~TestModelCount(){}
private:
	void createTestMLN1(LvrMLN& mln);
	void createTestMLN2(LvrMLN& mln);
	void createTestMLN3(LvrMLN& mln);
	void createTestMLN4(LvrMLN& mln);
	void createTestMLN5(LvrMLN& mln);
	void createTestMLN6(LvrMLN& mln);
	void createTestMLN7(LvrMLN& mln);
	void createTestMLN8(LvrMLN& mln);
	void createTestMLN9(LvrMLN& mln);
	void createTestMLN10(LvrMLN& mln);
	void createTestMLN11(LvrMLN& mln);
	void createTestMLN12(LvrMLN& mln);
	void createTestMLN13(LvrMLN& mln);
	void createTestMLN15(LvrMLN& mln);
	void createTestMLN16(LvrMLN& mln);
	void createTestMLN17(LvrMLN& mln);
	void createTestMLN18(LvrMLN& mln);
	void createTestMLN19(LvrMLN& mln);
	void createTestMLN20(LvrMLN& mln);
	void createTestMLN21(LvrMLN& mln);
	LogDouble manuallyTestModelCount2();
	LogDouble manuallyTestModelCount3();
	LogDouble manuallyTestModelCount5();
	LogDouble manuallyTestModelCount6();
	LogDouble manuallyTestModelCount8();
	LogDouble manuallyTestModelCount9();
	LogDouble manuallyTestModelCount10();
	LogDouble manuallyTestModelCount11();
	LogDouble manuallyTestModelCount17();
	
	LogDouble manuallyTestModelCountxx();
	void createTestFormula(LvrMLN& mln);
	void manuallyTestModelCount(int n);
	void runMLNTest(string folder,int id);
	bool getNextTruthAssignment(vector<int>& truthValues);
	void createMLNs(int testId);
	void doWMCTests(int testId);
	void manuallyTestWC(int testId);
};
#endif
