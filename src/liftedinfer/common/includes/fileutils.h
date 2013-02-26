#ifndef __FILE_UTILS
#define __FILE_UTILS
#include <fstream>
#include <vector>
using namespace std;
#include "logdouble.h"

struct LFileUtils
{
	static LFileUtils* Instance();
	//static void createOutFileInstance(string outfilename_,vector<string> queries_, vector<string> evidence_);
	static void createOutFileInstance(string outfilename);
	void updateFile(vector<LogDouble> queryValues,vector<string> queryStrings);
	void updateFile(vector<double> queryValues,vector<string> queryStrings);
private:
	LFileUtils(string outfilename);
	static LFileUtils* m_pInstance;
	//ofstream* outfile;
	string outfilename;
};
#endif
