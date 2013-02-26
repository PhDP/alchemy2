#include "testmodelcount.h"
#include "arguments.h"
#include "testformulainfer.h"
#include <sstream>
using namespace std;
ARGS ARGS::Args[] =
{
    // BEGIN: Common arguments
  ARGS()
};


int main(int argc, char* argv[])
{
	if(argc != 3)
	{
		cout<<"Usage::./runliftedinfertests autotestfolder(ending with slash) numberoftests, e.g.runliftedinfertests ../exdata/autotestmlnfiles/ 5"<<endl;
		return 0;
	}
	string s(argv[1]);
	string s1(argv[2]);
	stringstream st(s1);
	int numtests;
	st>>numtests;
	cout<<"Running "<<numtests<<" tests..."<<endl;
	TestModelCount tm;
	tm.runAllWMCTests();
	TestFormulaInfer ti(s,numtests);
	ti.runFormulaInferenceTests();
	return 0;
}
