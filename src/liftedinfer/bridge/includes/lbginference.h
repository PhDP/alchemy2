#ifndef LBGINFERENCE_H_
#define LBGINFERENCE_H_

#include "mcmc.h"
#include "mcsatparams.h"
#include "liftedalgsconvertor.h"
#include "lvparams.h"


class LBGInference : public MCMC
{
 public:

  LBGInference(VariableState* state, long int seed, const bool& trackClauseTrueCnts,
        LvrParams* lvParams_,
        Array<Array<Predicate* >* >* queryFormulas = NULL)
    : MCMC(state, seed, trackClauseTrueCnts, (MCMCParams*)lvParams_, queryFormulas)
  {
	  algsConvertor = new LiftedAlgsConvertor(state);
	  initialization = true;
	  //mcmcParams = mcmcParams_;
	  lvParams = lvParams_;
  }

  ~LBGInference()
  {
	  delete algsConvertor;
	  delete lvParams;
  }
  
  void init()
  {
    //update the weights
    //algsConvertor->initLiftedGibbsSampler();
  }

  void infer()
  {
	  if(initialization)
	  {
		  algsConvertor->processWeightLearningInput(lvParams,true);
		  initialization = false;
	  }
	  else
	  {
		  algsConvertor->processWeightLearningInput(lvParams,false);
	  }
	  //update the weights of the mln (needed during weight learning)
	  algsConvertor->updateState();
	  Array<double> wts;
	  state_->getMLN()->getClauseWts(wts);
	  vector<double> wtsA;
	  for(int i=0;i<wts.size();i++)
	  {
		  wtsA.push_back(wts[i]);
	  }
	  algsConvertor->updateMLNWeights(wtsA);
	  tallyCntsFromState();
  }
 private:
  LiftedAlgsConvertor* algsConvertor;
  bool initialization;
  //MCMCParams* mcmcParams;
  LvrParams* lvParams;
};

#endif /*MCSAT_H_*/
