#include "testformulainfer.h"
#include "blockexactinference.h"
#include "parser.h"
#include "ptpsearch.h"
#include "wmconvertor.h"
#include "stringconversionutils.h"
#include <sstream>
#include <fstream>
using namespace std;


bool TestFormulaInfer::getNextTruthAssignment(vector<int>& truthValues)
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

bool TestFormulaInfer::isEvidenceSat(vector<int> truthValue,vector<int> evidence)
{
	bool evidenceSat = true;
	for(int j=0;j<evidence.size();j++)
	{
		if(truthValue[evidence[j]]==0)
		{
			evidenceSat = false;
			break;
		}
	}
	return evidenceSat;
}


LogDouble TestFormulaInfer::manualEstimate(int id,bool withEvidence,bool literalWeights)
{
	vector<vector<vector<int> > > trueIndexes;
	vector<int> evidenceIndexes;
	vector<double> weights;
	int size;
	readGroundTruth(id,trueIndexes,evidenceIndexes,weights,size);
	vector<int> evidence;
	if(withEvidence)
	{
		evidence = evidenceIndexes;
	}
	int n=size;
	bool done =false;
	vector<int> truthValues;
	for(int i=0;i<n;i++)
	{
		truthValues.push_back(0);
	}
	LogDouble wt =LogDouble();
	while(!done)
	{
		LogDouble tot(1,false);
		for(unsigned int i=0;i<trueIndexes.size();i++)
		{
			int cnt = 0;
			for(unsigned int j=0;j<trueIndexes[i].size();j++)
			{
				bool trueGrnd=false;
				for(unsigned int k=0;k<trueIndexes[i].at(j).size();k++)
				{
					if(trueIndexes[i].at(j).at(k) < 0)
					{
						if(truthValues[-1*trueIndexes[i].at(j).at(k)-1] == 0)
						{
							trueGrnd=true;
							break;
						}
					}
					else
					{
						if(truthValues[trueIndexes[i].at(j).at(k)-1] == 1)
						{
							trueGrnd=true;
							break;
						}
					}
				}
				if(trueGrnd)
					cnt++;
			}
			if(!literalWeights)
			{
				double coeff = cnt*weights[i];
				tot *= LogDouble(coeff,true);
			}
			else
			{
				double v=pow(weights[i],(double)(trueIndexes[i].size() - cnt));
				tot *= LogDouble(v,false);
			}
		}	
		if(isEvidenceSat(truthValues,evidence))
			wt = wt + tot;
		if(getNextTruthAssignment(truthValues))
			break;
	}	
	return wt;
}

void TestFormulaInfer::doMLNTestExpWeighted(int index,string mlnfile,string dbfile,bool withEvidence)
{
	mln.clearData();
	LBlockExactInference lb(mln);
	LMParser p(mln);
	p.parseInputMLNFile(mlnfile);
	if(withEvidence)
	{
		int origSize=mln.clauses.size();
		p.parseDB(dbfile);
		int newSize=mln.clauses.size();
		vector<Atom*> evAtoms;
		if(origSize!=newSize)
		{
			for(int i=origSize;i<newSize;i++)
				evAtoms.push_back(mln.clauses[i]->atoms[0]);
			mln.preprocessEvidence(evAtoms);			
		}
		else
		{
			cout<<"No evidence file(.db) found for test mln"<<index<<endl;
			return;
		}
	}
	LogDouble mVal;
	LogDouble pVal;
	mVal=lb.weightedModelCount(mln.clauses);
	pVal = manualEstimate(index,withEvidence);
	if(!withEvidence)
	{
		cout<<"LGB Exact inference(formula weights) testmln"<<index<<"(no evidence)="<<mVal<<endl;
		cout<<"Manually computed value testmln"<<index<<"(no evidence)="<<pVal<<endl;
	}
	else
	{
		cout<<"LGB Exact inference(formula weights) testmln"<<index<<"(testdb"<<index<<" evidence)="<<mVal<<endl;
		cout<<"Manually computed value testmln"<<index<<"(testdb"<<index<<" evidence)="<<pVal<<endl;
	}
}

void TestFormulaInfer::doMLNTestPTPWeighted(int index,string mlnfile,string dbfile,bool withEvidence)
{
	mln.clearData();
	LPTPSearch ls =LPTPSearch(mln);
	LMParser p = LMParser(mln);
	p.parseInputMLNFile(mlnfile);
	if(withEvidence)
	{
		int origSize=mln.clauses.size();
		p.parseDB(dbfile);
		int newSize=mln.clauses.size();
		vector<Atom*> evAtoms;
		if(origSize==newSize)
		{
			cout<<"No evidence file(.db) found for test mln"<<index<<endl;
			return;
		}
	}
	LWMConvertor lc = LWMConvertor(mln);
	lc.convertMLN();
	LogDouble mVal;
	LogDouble pVal;
	LvrParams* params = new LvrParams;
	params->ptpCacheSize = 0;
	mVal=ls.startExactWeightedModelCounting(params);
	delete params;
	pVal = manualEstimate(index,withEvidence);
	if(!withEvidence)
	{
		cout<<"PTP Exact inference(Literal weights LWMC) testmln"<<index<<"(no evidence)="<<mVal<<endl;
		cout<<"Manually computed value testmln"<<index<<"(no evidence)="<<pVal<<endl;
	}
	else
	{
		cout<<"PTP Exact inference(Literal weights LWMC) testmln"<<index<<"(testdb"<<index<<" evidence)="<<mVal<<endl;
		cout<<"Manually computed value testmln"<<index<<"(testdb"<<index<<" evidence)="<<pVal<<endl;
	}
}

void TestFormulaInfer::readGroundTruth(int id,vector<vector<vector<int> > >& trueIndexes, vector<int>& evidenceIndexes,
	vector<double>& weights,int& size)
{
	string filename(testfolder);
	//read the ground truth file
	filename.append("groundtruth");
	string st = LStringConversionUtils::toString(id);
	filename.append(st);
	filename.append(".txt");
	ifstream filestr(filename.c_str());
	if(filestr == NULL)
	{
		cout<<"Error reading Ground Truth File"<<endl;
		exit(-1);
	}
	char* buf = new char[1024];
	//line 1=size
	{
		filestr.getline(buf,1024);
		string line(buf);
		size = LStringConversionUtils::toInt(line);
	}

	//line 2=evidence
	{
		filestr.getline(buf,1024);
		string line(buf);
		//check if it has evidence or not
		if(line.find("NOEVIDENCE")==string::npos)
		{
			//collect evidence indexes
			vector<string> evs;
			LStringConversionUtils::tokenize(line,evs," ");
			LStringConversionUtils::toIntArr(evs,evidenceIndexes);
		}
	}
	//line 3=numclauses wt1 wt2 ...
	int numClauses;
	{
		filestr.getline(buf,1024);
		string line(buf);
		vector<string> evs;
		LStringConversionUtils::tokenize(line,evs," ");
		numClauses = LStringConversionUtils::toInt(evs[0]);
		evs.erase(evs.begin());
		LStringConversionUtils::toDoubleArr(evs,weights);
	}
	trueIndexes.resize(numClauses);
	int trueIndexIter=0;
	for(int t=0;t<numClauses;t++)
	{
		filestr.getline(buf,1024);
		string line(buf);
		int gSize = LStringConversionUtils::toInt(line);
		trueIndexes[t].resize(gSize);
		for(int i=0;i<gSize;i++)
		{
			filestr.getline(buf,1024);
			string line(buf);
			vector<string> inds;
			LStringConversionUtils::tokenize(line,inds," ");
			LStringConversionUtils::toIntArr(inds,trueIndexes[t].at(i));
		}
	}
	filestr.close();
}

void TestFormulaInfer::runFormulaInferenceTests()
{
	
	cout<<"***********Exact Inference(literal weights) PTP Tests*************"<<endl;
	for(int i=1;i<=numtests;i++)	
	{
		string mlnFile(testfolder);
		string dbFile(testfolder);
		vector<int> evidence;
		int id=i;
		stringstream st;
		st<<id;
		mlnFile.append("testmln");
		mlnFile.append(st.str());
		mlnFile.append(".txt");
		dbFile.append("testdb");
		dbFile.append(st.str());
		dbFile.append(".txt");
		doMLNTestPTPWeighted(i,mlnFile,dbFile);
		doMLNTestPTPWeighted(i,mlnFile,dbFile,true);
	}
	
	cout<<"***********Exact Inference(clause weights) LBG Tests*************"<<endl;
	for(int i=1;i<=numtests;i++)	
	{
		string mlnFile(testfolder);
		string dbFile(testfolder);
		vector<int> evidence;
		string st=LStringConversionUtils::toString(i);
		mlnFile.append("testmln");
		mlnFile.append(st);
		mlnFile.append(".txt");
		dbFile.append("testdb");
		dbFile.append(st);
		dbFile.append(".txt");
		doMLNTestExpWeighted(i,mlnFile,dbFile);
		doMLNTestExpWeighted(i,mlnFile,dbFile,true);
	}
	
}
