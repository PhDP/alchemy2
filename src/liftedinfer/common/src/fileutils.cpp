
#include <fstream>
#include <vector>
using namespace std;
#include <assert.h>
#include "fileutils.h"
LFileUtils* LFileUtils::m_pInstance = NULL;

LFileUtils* LFileUtils::Instance()
{
	assert(m_pInstance!=NULL);
	return m_pInstance;
}

void LFileUtils::createOutFileInstance(string outfilename)
{
	m_pInstance = new LFileUtils(outfilename);
}

LFileUtils::LFileUtils(string outfilename_):outfilename(outfilename_)
{
	//outfile = new ofstream(outfilename.c_str());
}

void LFileUtils::updateFile(vector<LogDouble> queryValues,vector<string> queryStrings)
{
	ofstream out(outfilename.c_str());
	assert(queryStrings.size() == queryValues.size());
	for(unsigned int i=0;i<queryValues.size();i++)
	{
		//(*outfile)<<queryStrings[i].c_str()<<" "<<queryValues[i]<<",";
		out<<queryStrings[i].c_str()<<" ";
		bool inlogspace=false;
		double val = queryValues[i].getPrintableValue(inlogspace);
		if(inlogspace)
			out<<val<<"(Log space)";
		else
			out<<val;
		out<<endl;
	}
	//(*outfile)<<endl;
	out.close();
}

void LFileUtils::updateFile(vector<double> queryValues,vector<string> queryStrings)
{
	ofstream out(outfilename.c_str());
	assert(queryStrings.size() == queryValues.size());
	for(unsigned int i=0;i<queryValues.size();i++)
	{
		//(*outfile)<<queryStrings[i].c_str()<<" "<<queryValues[i]<<",";
		out<<queryStrings[i].c_str()<<" "<<queryValues[i]<<endl;
	}
	//(*outfile)<<endl;
	out.close();
}

