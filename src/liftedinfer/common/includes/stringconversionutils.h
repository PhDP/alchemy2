#ifndef __LSTRINGCONVERSION_UTILS
#define __LSTRINGCONVERSION_UTILS
#include <string>
#include <vector>
#include <sstream>
using namespace std;
struct LStringConversionUtils
{
	static void tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
		string substring = str.substr(lastPos, pos - lastPos);
        tokens.push_back(substring);
        // Skip delimiters. 
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

	static int toInt(string s)
	{
		stringstream convert(s);
		int i;
		convert >> i;
		return i;
	}
	static double toDouble(string s)
	{
		stringstream convert(s);
		double i;
		convert >> i;
		return i;
	}
	static void toIntArr(vector<string>& s,vector<int>& intArr)
	{
		for(unsigned int i=0;i<s.size();i++)
			intArr.push_back(toInt(s[i]));
	}
	static void toDoubleArr(vector<string>& s,vector<double>& dArr)
	{
		for(unsigned int i=0;i<s.size();i++)
			dArr.push_back(toDouble(s[i]));
	}
	static string toString(int i)
	{
		stringstream st;
		st<<i;
		return st.str();
	}
};
#endif