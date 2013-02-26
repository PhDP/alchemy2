#ifndef __LVPARAMS
#define __LVPARAMS
#include "mcmcparams.h"

//Sampling modes for importance sampling
enum ESamplingMode
{
	EUNIFORM,
	EBINOMIAL,
	EINFORMED,
	EINFORMEDV1
};

//Clustering modes for LBG
enum EClusteringMode
{
	ESEQUENTIAL,
	ERANDOM
};

//MODES of operating for LBG
enum ELBGOperationMode
{
	ESAMPLER,
	ECLUSTERING,
	EPARALLEL
};

enum EAlgorithm
{
	LVG,
	LISZ,
	LISM,
	PTPZ,
	PTPE
};

struct LvrParams:public MCMCParams
{
	EAlgorithm algorithm;
	bool isWeightLearning;
	string resultFile;
	string inClusterFile;
	string outClusterFile;
	bool useDumpFiles;
	ESamplingMode samplingMode;
	double learningRate;
	EClusteringMode clusteringMode;
	ELBGOperationMode operationMode;
	int proposalUpdateInterval;
	int ptpCacheSize;
	double baselineCostMultiplicativeFactor;
};

#endif
