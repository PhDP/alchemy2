#include "testmodelcount.h"
#include "ptpsearch.h"
#include "wmconvertor.h"
#include "cleanuputils.h"
#include "parser.h"
#include "blockexactinference.h"

#include <sstream>
using namespace std;

bool TestModelCount::getNextTruthAssignment(vector<int>& truthValues)
{
	int ind = truthValues.size()-1;
	truthValues[ind]++;
	while(truthValues[ind] == 2)
	{
		truthValues[ind]=0;
		if((ind-1) >= 0)
			truthValues[ind-1]++;
		else
		{
			return true;
		}
		ind--;
	}	
	return false;
}

void TestModelCount::manuallyTestModelCount(int n)
{
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		//if((truthValues[0]||truthValues[2])&&(truthValues[0]||truthValues[3])&&(truthValues[1]||truthValues[4])&&(truthValues[1]||truthValues[5]))
		//if(truthValues[0]&&truthValues[1])
		//if(truthValues[0]&&truthValues[1]&&truthValues[2]&&truthValues[3])
		//if((truthValues[0]||truthValues[2])&&(truthValues[1]||truthValues[3]))	
		//if((truthValues[0]||truthValues[4])&&(truthValues[1]||truthValues[5])&&(truthValues[2]||truthValues[6])&&(truthValues[3]||truthValues[7]))
		//if((truthValues[0]||truthValues[4])&&(truthValues[1]||truthValues[5])&&(truthValues[2]||truthValues[6])&&(truthValues[3]||truthValues[7])
		//	&&(truthValues[0]||truthValues[8])&&(truthValues[0]||truthValues[9])&&(truthValues[1]||truthValues[8])&&(truthValues[1]||truthValues[9])
		//	&& (truthValues[2]||truthValues[10])&&(truthValues[2]||truthValues[11])&&(truthValues[3]||truthValues[10])&&(truthValues[3]||truthValues[11]))
		//if((truthValues[0]||truthValues[4])&&(truthValues[8]||truthValues[12])&&(truthValues[8]||truthValues[13])&&(truthValues[9]||truthValues[12])&&(truthValues[9]||truthValues[13]))
		//if((truthValues[0]||truthValues[2])&&(truthValues[1]||truthValues[3])&&(truthValues[0]||truthValues[4])
		//	&&(truthValues[0]||truthValues[5])&&(truthValues[1]||truthValues[4])&&(truthValues[1]||truthValues[5]))
		//if((truthValues[0]||truthValues[1])&&(truthValues[0]||truthValues[3])&&(truthValues[0]||truthValues[4])
		//	&&(truthValues[2]||truthValues[3])&&(truthValues[2]||truthValues[4]))
		//if((truthValues[0]||truthValues[2])&&(truthValues[1]||truthValues[3])&&(truthValues[0]||truthValues[4])
		//	&&(truthValues[0]||truthValues[5])&&(truthValues[3]||truthValues[4])&&(truthValues[3]||truthValues[5]))
		
//		if(truthValues[0]==0 && truthValues[1]==0)
//		{
//		if((truthValues[0]||truthValues[2])&&(truthValues[1]||truthValues[3])&&(truthValues[0]||truthValues[4])
//			&&(truthValues[0]||truthValues[5])&&(truthValues[1]||truthValues[4])&&(truthValues[1]||truthValues[5]))
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			//LDPower(LogDouble(0.1,false),4,t);
			//LDPower(LogDouble(0.1,false),6,t);
			LogDouble::LDPower(LogDouble(0.1,false),6,t);
			//LDPower(LogDouble(0.1,false),2,t);
			//LogDouble t =LogDouble(pow(0.1,6),false);
			//LogDouble d = LogDouble(0.1,false);
			//for(int x=0;x<6;x++)
			//	t*=d;
			wt+=t;
			//wt1+=pow(0.1,8);
			//wt1+=pow(0.1,6);
			//wt1+=pow(0.1,4);
			//wt1+=pow(0.1,2);
		}
		
		getNextTruthAssignment(truthValues);
	}
	

	cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<wt1<<endl;
	cout<<"Sat Count="<<satcnt<<endl;
}

LogDouble TestModelCount::manuallyTestModelCount2()
{
	int n=12;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((truthValues[0]||truthValues[4])&&(truthValues[1]||truthValues[5])&&(truthValues[2]||truthValues[6])&&(truthValues[3]||truthValues[7])
			&&(truthValues[0]||truthValues[8])&&(truthValues[0]||truthValues[9])&&(truthValues[1]||truthValues[8])&&(truthValues[1]||truthValues[9])
			&& (truthValues[2]||truthValues[10])&&(truthValues[2]||truthValues[11])&&(truthValues[3]||truthValues[10])&&(truthValues[3]||truthValues[11]))
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}

LogDouble TestModelCount::manuallyTestModelCount3()
{
	int n=8;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((truthValues[0]||truthValues[4])&&(truthValues[1]||truthValues[5])
		&&(truthValues[2]||truthValues[6])&&(truthValues[3]||truthValues[7]))		
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}

LogDouble TestModelCount::manuallyTestModelCount5()
{
	int n=12;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((truthValues[0])&&(truthValues[1])&&(truthValues[2])&&(truthValues[3])
			&&(truthValues[0]||truthValues[4]||truthValues[8])&&(truthValues[0]||truthValues[4]||truthValues[9])
			&&(truthValues[0]||truthValues[4]||truthValues[10])&&(truthValues[0]||truthValues[4]||truthValues[11])
			&&(truthValues[0]||truthValues[5]||truthValues[8])&&(truthValues[0]||truthValues[5]||truthValues[9])
			&&(truthValues[0]||truthValues[5]||truthValues[10])&&(truthValues[0]||truthValues[5]||truthValues[11])
			&&(truthValues[1]||truthValues[4]||truthValues[8])&&(truthValues[1]||truthValues[4]||truthValues[9])
			&&(truthValues[1]||truthValues[4]||truthValues[10])&&(truthValues[1]||truthValues[4]||truthValues[11])
			&&(truthValues[1]||truthValues[5]||truthValues[8])&&(truthValues[1]||truthValues[5]||truthValues[9])
			&&(truthValues[1]||truthValues[5]||truthValues[10])&&(truthValues[1]||truthValues[5]||truthValues[11])
			&&(truthValues[2]||truthValues[6]||truthValues[8])&&(truthValues[2]||truthValues[6]||truthValues[9])
			&&(truthValues[2]||truthValues[6]||truthValues[10])&&(truthValues[2]||truthValues[6]||truthValues[11])
			&&(truthValues[2]||truthValues[7]||truthValues[8])&&(truthValues[2]||truthValues[7]||truthValues[9])
			&&(truthValues[2]||truthValues[7]||truthValues[10])&&(truthValues[2]||truthValues[7]||truthValues[11])
			&&(truthValues[3]||truthValues[6]||truthValues[8])&&(truthValues[3]||truthValues[6]||truthValues[9])
			&&(truthValues[3]||truthValues[6]||truthValues[10])&&(truthValues[3]||truthValues[6]||truthValues[11])
			&&(truthValues[3]||truthValues[7]||truthValues[8])&&(truthValues[3]||truthValues[7]||truthValues[9])
			&&(truthValues[3]||truthValues[7]||truthValues[10])&&(truthValues[3]||truthValues[7]||truthValues[11]))
	

		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}

LogDouble TestModelCount::manuallyTestModelCount6()
{
	int n=10;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((truthValues[0]||truthValues[2])&&(truthValues[0]||truthValues[3])
		&&(truthValues[1]||truthValues[4])&&(truthValues[1]||truthValues[5])
		&&(truthValues[2]||truthValues[6])&&(truthValues[2]||truthValues[7])
		&&(truthValues[4]||truthValues[6])&&(truthValues[4]||truthValues[7])
		&&(truthValues[3]||truthValues[8])&&(truthValues[3]||truthValues[9])
		&&(truthValues[5]||truthValues[8])&&(truthValues[5]||truthValues[9]))

		//&&truthValues[2]&&truthValues[3]&&truthValues[4]&&truthValues[5])		
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}

LogDouble TestModelCount::manuallyTestModelCount8()
{
	int n=8;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((truthValues[0]||truthValues[4])&&(truthValues[0]||truthValues[5])
		&&(truthValues[2]||truthValues[4])&&(truthValues[2]||truthValues[5])
		&&(truthValues[1]||truthValues[6])&&(truthValues[1]||truthValues[7])
		&&(truthValues[3]||truthValues[6])&&(truthValues[3]||truthValues[7]))
		//&&truthValues[2]&&truthValues[3]&&truthValues[4]&&truthValues[5])		
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}


LogDouble TestModelCount::manuallyTestModelCount9()
{
	int n=8;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((truthValues[0]||truthValues[4])&&(truthValues[0]||truthValues[5])
		&&(truthValues[2]||truthValues[4])&&(truthValues[2]||truthValues[5])
		&&(truthValues[1]||truthValues[6])&&(truthValues[1]||truthValues[7])
		&&(truthValues[3]||truthValues[6])&&(truthValues[3]||truthValues[7])
		&&truthValues[2]&&truthValues[3])		
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}
/*
LogDouble TestModelCount::manuallyTestModelCount9()
{
	int n=8;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((truthValues[0]||truthValues[4])&&(truthValues[0]||truthValues[5])
		&&(truthValues[2]||truthValues[4])&&(truthValues[2]||truthValues[5])
		&&(truthValues[1]||truthValues[6])&&(truthValues[1]||truthValues[7])
		&&(truthValues[3]||truthValues[6])&&(truthValues[3]||truthValues[7])
		&&truthValues[2]&&truthValues[3])		
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}
*/
LogDouble TestModelCount::manuallyTestModelCount10()
{
	int n=4;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if(truthValues[0]&&(truthValues[0]||truthValues[1])&&truthValues[1]&&(truthValues[1]||truthValues[2])
			&&truthValues[2]&&(truthValues[2]||truthValues[3])&&truthValues[3])
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}

LogDouble TestModelCount::manuallyTestModelCount11()
{
	int n=6;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((!truthValues[0]||truthValues[1])&&(!truthValues[1]||truthValues[0])&&(!truthValues[0]||truthValues[2])
			&&(!truthValues[0]||truthValues[4])&&(!truthValues[1]||truthValues[3])&&(!truthValues[1]||truthValues[5])
			&&(!truthValues[0]||!truthValues[3]||truthValues[1])&&(!truthValues[1]||!truthValues[4]||truthValues[0]))
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}


void TestModelCount::createTestMLN1(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}

	vector<LvrTerm*> terms(10);
	vector<int> domain(10);
	for (int i = 0; i < 10; i++)
		domain[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms[i] = new LvrTerm(0, domain);
	}

	mln.clauses = vector<WClause*>(4);
	{
		mln.clauses[0] = new WClause();
		mln.clauses[0]->atoms = vector<Atom*>(3);
		mln.clauses[0]->sign = vector<bool>(3);
		mln.clauses[0]->satisfied=false;
		{
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms[0];
			myterms[1] = terms[1];
			mln.clauses[0]->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms[2];
			myterms[1] = terms[1];
			mln.clauses[0]->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		{
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms[3];
			myterms[1] = terms[1];
			mln.clauses[0]->atoms[2] = new Atom(mln.symbols[2], myterms);
		}
	}

	{
		mln.clauses[1] = new WClause();
		mln.clauses[1]->atoms = vector<Atom*>(2);
		mln.clauses[1]->sign = vector<bool>(2);
		mln.clauses[1]->satisfied=false;
		{
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms[4];
			myterms[1] = terms[5];
			mln.clauses[1]->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms[5];
			myterms[1] = terms[6];
			mln.clauses[1]->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
	}

	{
		mln.clauses[2] = new WClause();
		mln.clauses[2]->atoms = vector<Atom*>(2);
		mln.clauses[2]->sign = vector<bool>(2);
		mln.clauses[2]->satisfied=false;
		{
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms[7];
			myterms[1] = terms[8];
			mln.clauses[2]->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms[9];
			myterms[1] = terms[8];
			mln.clauses[2]->atoms[1] = new Atom(mln.symbols[2], myterms);
		}
	}

	{
		mln.clauses[3] = new WClause();
		mln.clauses[3]->atoms = vector<Atom*>(1);
		mln.clauses[3]->sign = vector<bool>(1);
		mln.clauses[3]->satisfied=false;
		vector<LvrTerm*> myterms1(2);
		myterms1[0] = new LvrTerm(0, 1);
		myterms1[1] = new LvrTerm(0, 2);
		mln.clauses[3]->atoms[0] = new Atom(mln.symbols[0], myterms1);
	}
}

//R(x,y) V S(X,Y);R(X,Y) V T(X,Z)
void TestModelCount::createTestMLN2(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(5);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 5; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{
			vector<LvrTerm*> myterms(2);
			//vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[1] = new LvrTerm(0,25);	
			//myterms[2] = terms[2];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[2];
			myterms[1] = terms1[3];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[2];
			//myterms[1] = terms1[3];
			myterms[1] = terms1[4];
			//myterms[2] = terms[2];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}

}

//R(x,y) V S(X,Y);
void TestModelCount::createTestMLN3(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(5);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 5; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{
			vector<LvrTerm*> myterms(2);
			//vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[1] = new LvrTerm(0,25);	
			//myterms[2] = terms[2];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}

}


void TestModelCount::createTestMLN4(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(5);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 5; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,1);
			myterms[1] = new LvrTerm(0,0);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,1);
			myterms[1] = new LvrTerm(0,1);
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,1);
			myterms[1] = terms1[0];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}

}


void TestModelCount::createTestMLN5(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[2];
			myterms[1] = terms1[3];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[3];
			myterms[1] = terms1[4];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[5];
			myterms[1] = terms1[6];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[2], myterms);
		}

		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		mln.clauses.push_back(clause);
	}
}

void TestModelCount::createTestMLN6(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	

	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_type1(1);
	var_type1[0] = 0;
	vector<int> var_type2(2);
	var_type2[0] = var_type2[1] = 0;

	mln.symbols[0] = new PredicateSymbol(0, symbols[0], var_type1,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[1] = new PredicateSymbol(1, symbols[1], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[2] = new PredicateSymbol(2, symbols[2], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));

	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[0];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
/*		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,25);
			myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[2], myterms);
		}
		*/
		mln.clauses.push_back(clause);
	}
	
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[3];
			myterms[1] = terms1[4];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[4];
			myterms[1] = terms1[5];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}
		
		mln.clauses.push_back(clause);
	}
	
}


void TestModelCount::createTestMLN7(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(10);
	for (int i = 0; i < 10; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,0);
			myterms[1] = terms1[0];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[1];
			myterms[1] = terms1[2];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,1);
			myterms[1] = terms1[3];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[4];
			myterms[1] = terms1[5];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
//		mln.clauses.push_back(clause);
	}

}

void TestModelCount::createTestMLN8(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[3];
			myterms[1] = terms1[4];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[4];
			myterms[1] = terms1[5];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}
		
		mln.clauses.push_back(clause);
	}
	
}

void TestModelCount::createTestMLN9(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = true;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,0);
			myterms[1] = terms1[0];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,1);
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[2];
			myterms[1] = terms1[3];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[3];
			myterms[1] = terms1[4];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}

}

void TestModelCount::createTestMLN10(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[2];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}

		mln.clauses.push_back(clause);
	}
}

void TestModelCount::createTestMLN11(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_type1(1);
	var_type1[0] = 0;
	vector<int> var_type2(2);
	var_type2[0] = var_type2[1] = 0;

	mln.symbols[0] = new PredicateSymbol(0, symbols[0], var_type1,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[1] = new PredicateSymbol(1, symbols[1], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[2] = new PredicateSymbol(2, symbols[2], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[0];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}


		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[0], myterms);
		}

		mln.clauses.push_back(clause);
	}
	
	
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[3];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[4];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[0], myterms);
		}
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[5];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[6];
			myterms[1] = terms1[5];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}
	
}

void TestModelCount::createTestMLN12(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(5);	
	vector<string> symbols(5);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	symbols[3] = "P";
	symbols[4] = "Q";
	vector<int> var_type1(1);
	var_type1[0] = 0;
	vector<int> var_type2(2);
	var_type2[0] = var_type2[1] = 0;
	mln.symbols[0] = new PredicateSymbol(0, symbols[0], var_type1,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[1] = new PredicateSymbol(1, symbols[1], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[2] = new PredicateSymbol(2, symbols[2], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[3] = new PredicateSymbol(3, symbols[3], var_type1,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[4] = new PredicateSymbol(4, symbols[4], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));

	mln.clauses.clear();
	vector<LvrTerm*> terms1(20);
	vector<int> domain1(10);
	for (int i = 0; i < 10; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 20; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}
	
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[0];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[2];
			myterms[1] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[2], myterms);
		}

		mln.clauses.push_back(clause);
	}
	
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[3];
			myterms[1] = terms1[4];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[3];
			clause->atoms[1] = new Atom(mln.symbols[3], myterms);
		}

		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[5];
			clause->atoms[0] = new Atom(mln.symbols[3], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[6];
			myterms[1] = terms1[5];
			clause->atoms[1] = new Atom(mln.symbols[4], myterms);
		}
		mln.clauses.push_back(clause);
	}
	
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[7];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[8];
			myterms[1] = terms1[7];
			clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[4], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[9];
			myterms[1] = terms1[10];
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[10];
			myterms[1] = terms1[11];
			clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[4], myterms);
		}
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[11];
			clause->atoms[2] = new Atom(mln.symbols[0], myterms);
		}

		mln.clauses.push_back(clause);
	}	

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[12];
			clause->atoms[0] = new Atom(mln.symbols[3], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[13];
			myterms[1] = terms1[12];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}	

}


void TestModelCount::createTestMLN13(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(7);	
	vector<string> symbols(7);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	symbols[3] = "P";
	symbols[4] = "Q";
	symbols[5] = "U";
	symbols[6] = "V";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 7; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(30);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 30; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[2];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}


	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[3];
			myterms[1] = terms1[4];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[5];
			myterms[1] = terms1[4];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}

		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[6];
			myterms[1] = terms1[7];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[5], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[8];
			myterms[1] = terms1[9];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[5], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[10];
			myterms[1] = terms1[11];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[6], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[12];
			myterms[1] = terms1[13];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[6], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[14];
			myterms[1] = terms1[15];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[3], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[15];
			myterms[1] = terms1[16];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[4], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[17];
			myterms[1] = terms1[18];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[18];
			myterms[1] = terms1[19];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[3], myterms);
		}

		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = true;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[20];
			myterms[1] = terms1[21];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[21];
			myterms[1] = terms1[22];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[4], myterms);
		}

		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = true;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[23];
			myterms[1] = terms1[24];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = true;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[25];
			myterms[1] = terms1[26];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}

}
void TestModelCount::createMLNFromId(int testId)
{
	createMLNs(testId);
	if(mln.clauses.size()==0)
		return;
	vector<WClause*> tempClauses;
	for(int i=0;i<mln.clauses.size();i++)
	{
		WClause* nClause = LvrMLN::create_new_clause(mln.clauses[i]);
		tempClauses.push_back(nClause);
	}
	mln.clauses.clear();
	mln.clauses=tempClauses;
	for(int i=0;i<mln.clauses.size();i++)
	{
		mln.clauses[i]->print();
	}
}
void TestModelCount::doWMCTests(int testId)
{
	/*createMLNs(testId);
	if(mln.clauses.size()==0)
		return;
	vector<WClause*> tempClauses;
	for(int i=0;i<mln.clauses.size();i++)
	{
		WClause* nClause = LvrMLN::create_new_clause(mln.clauses[i]);
		tempClauses.push_back(nClause);
	}
	mln.clauses.clear();
	mln.clauses=tempClauses;
	for(int i=0;i<mln.clauses.size();i++)
	{
		mln.clauses[i]->print();
	}
	*/
	mln.clearData();
	createMLNFromId(testId);
	if(mln.clauses.size()==0)
		return;

	LPTPSearch ls =LPTPSearch(mln);
	LNormalizer ln(mln);
	ln.normalizeClauses(mln.clauses,true,false);
	LvrParams* params = new LvrParams;
	params->ptpCacheSize = 0;
	cout<<"PTP Value="<<ls.startExactWeightedModelCounting(params)<<endl;
	delete params;
	manuallyTestWC(testId);
	cleanup(mln.clauses);
}

void TestModelCount::createMLNs(int testId)
{
	switch(testId)
	{
	case 1:
		//createTestMLN1(mln);
		break;
	case 2:
		createTestMLN2(mln);
		break;
	case 3:
		createTestMLN3(mln);
		break;
	case 4:
		//createTestMLN4(mln);
		break;
	case 5:
		createTestMLN5(mln);
		break;
	case 6:
		createTestMLN6(mln);
		break;
	case 7:
		//createTestMLN7(mln);
		break;
	case 8:
		createTestMLN8(mln);
		break;
	case 9:
		createTestMLN9(mln);
		break;
	case 10:
		createTestMLN10(mln);
		break;
	case 11:
		createTestMLN11(mln);
		break;
	case 12:
		//createTestMLN12(mln);
		break;
	case 13:
		//createTestMLN13(mln);
		break;
	case 14:
		//createTestMLN14(mln);
		break;
	case 15:
		//createTestMLN15(mln);
		break;
	case 16:
		//createTestMLN16(mln);
		break;
	case 17:
		createTestMLN17(mln);
		break;
	case 18:
		//createTestMLN18(mln);
		break;
	case 19:
		//createTestMLN19(mln);
		break;
	case 20:
		//createTestMLN20(mln);
		break;
	case 21:
		//createTestMLN21(mln);
		break;

	}
	
}

void TestModelCount::manuallyTestWC(int testId)
{
	switch(testId)
	{
	case 1:
		//cout<<"Manual Value="<< manuallyTestModelCount1()<<endl;
		break;
	case 2:
		cout<<"Manual Value="<< manuallyTestModelCount2()<<endl;
		break;
	case 3:
		cout<<"Manual Value="<< manuallyTestModelCount3()<<endl;
		break;
	case 4:
		//cout<<"Manual Value="<< manuallyTestModelCount4()<<endl;
		break;
	case 5:
		cout<<"Manual Value="<< manuallyTestModelCount5()<<endl;
		break;
	case 6:
		cout<<"Manual Value="<< manuallyTestModelCount6()<<endl;
		break;
	case 7:
		//cout<<"Manual Value="<< manuallyTestModelCount7()<<endl;
		break;
	case 8:
		cout<<"Manual Value="<< manuallyTestModelCount8()<<endl;
		break;
	case 9:
		cout<<"Manual Value="<< manuallyTestModelCount9()<<endl;
		break;
	case 10:
		cout<<"Manual Value="<< manuallyTestModelCount10()<<endl;
		break;
	case 11:
		cout<<"Manual Value="<< manuallyTestModelCount11()<<endl;
		break;
	case 12:
		//cout<<"Manual Value="<< manuallyTestModelCount12()<<endl;
		break;
	case 13:
		//cout<<"Manual Value="<< manuallyTestModelCount13()<<endl;
		break;
	case 14:
		//cout<<"Manual Value="<< manuallyTestModelCount14()<<endl;
		break;
	case 15:
		//cout<<"Manual Value="<< manuallyTestModelCount15()<<endl;
		break;
	case 16:
		//cout<<"Manual Value="<< manuallyTestModelCount16()<<endl;
		break;
	case 17:
		cout<<"Manual Value="<< manuallyTestModelCount17()<<endl;
		break;
	case 18:
		//cout<<"Manual Value="<< manuallyTestModelCount18()<<endl;
		break;
	case 19:
		//cout<<"Manual Value="<< manuallyTestModelCount19()<<endl;
		break;
	case 20:
		//cout<<"Manual Value="<< manuallyTestModelCount20()<<endl;
		break;
	case 21:
		//cout<<"Manual Value="<< manuallyTestModelCount21()<<endl;
		break;

	}
	

}
void TestModelCount::runAllWMCTests()
{
	cout<<"*******************Exact Weighted Model Counting(PTP) tests************************"<<endl;
	for(int i=1;i<=21;i++)
	{
		cout<<"***********WMCTEST"<<i<<"********"<<endl;
		doWMCTests(i);
		cout<<"***************************"<<endl;
	}
}

LogDouble TestModelCount::manuallyTestModelCountxx()
{
	int n=10;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((truthValues[0]||truthValues[2]||!truthValues[6])&&(truthValues[0]||truthValues[2]||!truthValues[7])&&
			(truthValues[0]||truthValues[3]||!truthValues[6])&&(truthValues[0]||truthValues[3]||!truthValues[7])&&
			(truthValues[1]||truthValues[4]||!truthValues[8])&&(truthValues[1]||truthValues[4]||!truthValues[9])&&
			(truthValues[1]||truthValues[5]||!truthValues[8])&&(truthValues[1]||truthValues[5]||!truthValues[9])&&
			(truthValues[0]||truthValues[6])&&(truthValues[0]||truthValues[7])&&(truthValues[1]||truthValues[8])&&(truthValues[1]||truthValues[9])&&
			(truthValues[2]||truthValues[6])&&(truthValues[3]||truthValues[7])&&(truthValues[4]||truthValues[8])&&(truthValues[5]||truthValues[9]))
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			if(truthValues[6]==0) 
				t=t*LogDouble(0.5,false);
			if(truthValues[7]==0) 
				t=t*LogDouble(0.5,false);
			if(truthValues[8]==0) 
				t=t*LogDouble(0.5,false);
			if(truthValues[9]==0) 
				t=t*LogDouble(0.5,false);

			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}



void TestModelCount::createTestFormula(LvrMLN& mln)
{
	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_type1(1);
	var_type1[0] = 0;
	vector<int> var_type2(2);
	var_type2[0] = var_type2[1] = 0;

	mln.symbols.push_back( new PredicateSymbol(0, symbols[0], var_type1,LogDouble(0.1,false),LogDouble(0.1,false)));
	mln.symbols.push_back(new PredicateSymbol(1, symbols[1], var_type2,LogDouble(0.1,false),LogDouble(0.1,false)));
	mln.symbols.push_back(new PredicateSymbol(2, symbols[2], var_type2,LogDouble(0.1,false),LogDouble(0.1,false)));
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[0];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}
	
	
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[3];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[3];
			myterms[1] = terms1[4];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}

		mln.clauses.push_back(clause);
	}
	/*
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(1);
			myterms[0] = terms1[5];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[6];
			myterms[1] = terms1[5];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}
	*/
	mln.formulas.push_back(new Formula(0,2,LogDouble(0.5,false)));
	//mln.formulas.push_back(new Formula(2,3,LogDouble(0.75,false)));

}

void TestModelCount::createTestMLN15(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(5);	

	vector<string> symbols(5);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	symbols[3] = "P";
	symbols[4] = "Q";
	
	vector<int> var_type1(1);
	var_type1[0] = 0;
	vector<int> var_type2(2);
	var_type2[0] = var_type2[1] = 0;

	mln.symbols[0] = new PredicateSymbol(0, symbols[0], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[1] = new PredicateSymbol(1, symbols[1], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[2] = new PredicateSymbol(2, symbols[2], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[3] = new PredicateSymbol(3, symbols[3], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[4] = new PredicateSymbol(4, symbols[4], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));

	mln.clauses.clear();
	vector<LvrTerm*> terms1(25);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 25; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[1];
			myterms[1] = terms1[2];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[3];
			//myterms[2] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[2], myterms);
		}
		
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[4];
			myterms[1] = terms1[5];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[3], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[5];
			myterms[1] = terms1[6];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[0], myterms);
		}
		
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[7];
			myterms[1] = terms1[8];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[7];
			myterms[1] = terms1[9];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			//clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[4], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[7];
			myterms[1] = terms1[10];
			//myterms[2] = terms1[1];
			clause->sign[2]=true;
			clause->atoms[2] = new Atom(mln.symbols[3], myterms);
		}
		
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[11];
			myterms[1] = terms1[12];
			//myterms[2] = terms1[1];
			//clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[3], myterms);
		}
		
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[13];
			myterms[1] = terms1[14];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[4], myterms);
		}
		
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[15];
			myterms[1] = terms1[16];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[15];
			myterms[1] = terms1[17];
			clause->sign[1]=true;
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[4], myterms);
		}
		
		mln.clauses.push_back(clause);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[18];
			myterms[1] = terms1[19];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}

		
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[3];
			myterms[1] = terms1[4];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[4];
			myterms[1] = terms1[5];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}
		
		mln.clauses.push_back(clause);
	}
	
}


void TestModelCount::createTestMLN16(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(5);	

	vector<string> symbols(5);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	symbols[3] = "P";
	symbols[4] = "Q";
	
	vector<int> var_type1(1);
	var_type1[0] = 0;
	vector<int> var_type2(2);
	var_type2[0] = var_type2[1] = 0;

	mln.symbols[0] = new PredicateSymbol(0, symbols[0], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[1] = new PredicateSymbol(1, symbols[1], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[2] = new PredicateSymbol(2, symbols[2], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[3] = new PredicateSymbol(3, symbols[3], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));
	mln.symbols[4] = new PredicateSymbol(4, symbols[4], var_type2,LogDouble(0.1,false),LogDouble(0.1,false));

	mln.clauses.clear();
	vector<LvrTerm*> terms1(25);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 25; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[1];
			myterms[1] = terms1[2];
			//myterms[1] = new LvrTerm(0,25);
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[0], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[1];
			myterms[1] = terms1[3];
			//myterms[2] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[1], myterms);
		}
		
		mln.clauses.push_back(clause);
	}
}

void TestModelCount::createTestMLN17(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[1];
			myterms[1] = terms1[2];
			//myterms[2] = terms1[1];
			clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[2];
			//myterms[2] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[0], myterms);
		}


		mln.clauses.push_back(clause);
	}
}

LogDouble TestModelCount::manuallyTestModelCount17()
{
	int n=4;
	bool done =false;
	int satcnt=0;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	double wt1=0;
	//cout<<"Enumerating all combinations..."<<endl;
	while(!done)
	{
		if((!truthValues[0] || !truthValues[0] || truthValues[0])&&
			(!truthValues[0] || !truthValues[1] || truthValues[1])&&
			(!truthValues[1] || !truthValues[2] || truthValues[0])&&
			(!truthValues[1] || !truthValues[3] || truthValues[1])&&
			(!truthValues[2] || !truthValues[0] || truthValues[2])&&
			(!truthValues[2] || !truthValues[1] || truthValues[3])&&
			(!truthValues[3] || !truthValues[2] || truthValues[2])&&
			(!truthValues[3] || !truthValues[3] || truthValues[3]))
		{
			satcnt++;
			LogDouble t= LogDouble(1,false);
			LogDouble::LDPower(LogDouble(0.1,false),n,t);
			wt+=t;
		}
		if(getNextTruthAssignment(truthValues))
			break;
	}
	
	//cout<<"Manually computed Weighted Model Count="<<wt<<endl;
	//cout<<"Sat Count="<<satcnt<<endl;
	return wt;
}

void TestModelCount::createTestMLN18(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(4);
		clause->sign = vector<bool>(4);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[2];
			//myterms[2] = terms1[1];
			clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[0], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[3];
			//myterms[2] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[0], myterms);
		}



		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[4];
			myterms[1] = terms1[0];
			//myterms[2] = terms1[1];
			clause->atoms[3] = new Atom(mln.symbols[1], myterms);
		}

		mln.clauses.push_back(clause);
	}


	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[5];
			myterms[1] = terms1[6];
			//myterms[2] = terms1[1];
			clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[7];
			myterms[1] = terms1[6];
			//myterms[2] = terms1[1];
			clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[6];
			myterms[1] = terms1[8];
			//myterms[2] = terms1[1];
			clause->sign[1]=true;
			clause->atoms[2] = new Atom(mln.symbols[2], myterms);
		}

		mln.clauses.push_back(clause);
	}

}

void TestModelCount::createTestMLN19(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	vector<int> var_types1(2);
	var_types1[0] = 0;
	mln.symbols[0] = new PredicateSymbol(0, symbols[0], var_types1,LogDouble(0.1,false),LogDouble(0.1,false));
	for (int i = 1; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		/*{	
			vector<LvrTerm*> myterms(1);
			//myterms[0] = terms1[0];
			myterms[0] = new LvrTerm(0,1);
			//myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			//clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		*/
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			//myterms[0] = new LvrTerm(0,1);
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			//clause->sign[1]=true;
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[1];
			myterms[1] = terms1[2];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}

		mln.clauses.push_back(clause);
	}

}

void TestModelCount::createTestMLN20(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	vector<int> var_types1(2);
	var_types1[0] = 0;
	mln.symbols[0] = new PredicateSymbol(0, symbols[0], var_types1,LogDouble(0.1,false),LogDouble(0.1,false));
	for (int i = 1; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(2);
	for (int i = 0; i < 2; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(2);
		clause->sign = vector<bool>(2);
		clause->satisfied = false;
		/*{	
			vector<LvrTerm*> myterms(1);
			//myterms[0] = terms1[0];
			myterms[0] = new LvrTerm(0,1);
			//myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			//clause->sign[0]=true;
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		*/
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			//myterms[0] = new LvrTerm(0,1);
			myterms[1] = terms1[1];
			//myterms[2] = terms1[1];
			//clause->sign[1]=true;
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[2];
			myterms[1] = terms1[0];
			//myterms[2] = terms1[1];
			clause->atoms[1] = new Atom(mln.symbols[2], myterms);
		}

		mln.clauses.push_back(clause);
	}

}

void TestModelCount::createTestMLN21(LvrMLN& mln)
{
	mln.symbols = vector<PredicateSymbol*>(3);	
	vector<string> symbols(3);
	symbols[0] = "R";
	symbols[1] = "S";
	symbols[2] = "T";
	vector<int> var_types(2);
	var_types[0] = var_types[1] = 0;
	vector<int> var_types1(2);
	for (int i = 0; i < 3; i++) {
		mln.symbols[i] = new PredicateSymbol(i, symbols[i], var_types,LogDouble(0.1,false),LogDouble(0.1,false));
	}
	mln.clauses.clear();
	vector<LvrTerm*> terms1(10);
	vector<int> domain1(3);
	for (int i = 0; i < 3; i++)
		domain1[i] = i;
	for (size_t i = 0; i < 10; i++) {
		terms1[i] = new LvrTerm(0, domain1);
	}

	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(3);
		clause->sign = vector<bool>(3);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			myterms[1] = terms1[1];
			clause->atoms[0] = new Atom(mln.symbols[0], myterms);
		}
		
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[0];
			//myterms[0] = new LvrTerm(0,1);
			myterms[1] = terms1[2];
			//myterms[2] = terms1[1];
			//clause->sign[1]=true;
			clause->atoms[1] = new Atom(mln.symbols[1], myterms);
		}

		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = terms1[1];
			myterms[1] = terms1[2];
			//myterms[2] = terms1[1];
			clause->atoms[2] = new Atom(mln.symbols[2], myterms);
		}

		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,0);
			myterms[1] = new LvrTerm(0,0);
			clause->atoms[0] = new Atom(mln.symbols[1], myterms);
		}
		mln.clauses.push_back(clause);
	}
	{
		WClause* clause = new WClause();
		clause->atoms = vector<Atom*>(1);
		clause->sign = vector<bool>(1);
		clause->satisfied = false;
		{	
			vector<LvrTerm*> myterms(2);
			myterms[0] = new LvrTerm(0,1);
			myterms[1] = new LvrTerm(0,0);
			clause->atoms[0] = new Atom(mln.symbols[2], myterms);
		}
		mln.clauses.push_back(clause);
	}
}

