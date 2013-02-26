/*
 * All of the documentation and software included in the
 * Alchemy Software is copyrighted by Stanley Kok, Parag
 * Singla, Matthew Richardson, Pedro Domingos, Marc
 * Sumner, Hoifung Poon, Daniel Lowd, and Jue Wang.
 * 
 * Copyright [2004-09] Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang. All rights reserved.
 * 
 * Contact: Pedro Domingos, University of Washington
 * (pedrod@cs.washington.edu).
 * 
 * Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that
 * the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the
 * following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the
 * above copyright notice, this list of conditions and the
 * following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 * 
 * 3. All advertising materials mentioning features or use
 * of this software must display the following
 * acknowledgment: "This product includes software
 * developed by Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang in the Department of
 * Computer Science and Engineering at the University of
 * Washington".
 * 
 * 4. Your publications acknowledge the use or
 * contribution made by the Software to your research
 * using the following citation(s): 
 * Stanley Kok, Parag Singla, Matthew Richardson and
 * Pedro Domingos (2005). "The Alchemy System for
 * Statistical Relational AI", Technical Report,
 * Department of Computer Science and Engineering,
 * University of Washington, Seattle, WA.
 * http://alchemy.cs.washington.edu.
 * 
 * 5. Neither the name of the University of Washington nor
 * the names of its contributors may be used to endorse or
 * promote products derived from this software without
 * specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF WASHINGTON
 * AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OF WASHINGTON OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
#ifndef _SUPERPRED_H_DEC_2007
#define _SUPERPRED_H_DEC_2007

#include "clausecounter.h"
#include "superclause.h"
#include "hypercube.h"
#include "hypercubeoperations.h"
#include "domain.h"

class SuperPred;

//used in the SuperPred class
typedef hash_map<HyperCube*, SuperPred*, HashHyperCube, EqualHyperCube> HyperCubeToSuperPredMap;

//extern functions defined in superpred.cpp

/***************************************************************************************/
extern void initSuperPreds(Array<Array<SuperClause *> *> * const & superClausesArr,
		                  Domain * const & domain);

extern void createSuperPreds(Array<Array<SuperClause *> *> * const & superClausesArr,
		                     Domain * const & domain, 
							 bool isApproximateSuperPreds);

/***************************************************************************************/

class SuperPred {
	 public:
		  SuperPred(int & predId, ClauseCounter * const & clauseCounter, int parentSuperPredId) { 
			   hyperCubes_ = new Array<HyperCube *>;
			   predId_ = predId;
			   
			   //find the id of this superpred - this is simply the
			   //current cnt of the superpreds for this predId
			   Array<SuperPred *> * superPreds = (*(SuperPred::superPredsArr_))[predId];
			   superPredId_ = superPreds->size();
			   superPreds->append(this);
			   clauseCounter_ = clauseCounter;
			   parentSuperPredId_ = parentSuperPredId;
		  }

		  //not responsible for deleting the hyperCubes
		  ~SuperPred() {
			   delete hyperCubes_;
		       
			   //note: though, clauseCounter is allocated memory somewhere outside,
			   //it is delete here
			   delete clauseCounter_;
		  }
          
		  int getPredId() {return predId_;}
	 
		  HyperCube * getHyperCube(int hindex) {return (*hyperCubes_)[hindex];}

		  int getSuperPredId() {return superPredId_;}
		  
		  int getParentSuperPredId() {return parentSuperPredId_;}

		  int getNumHyperCubes(){ return hyperCubes_->size();}
          
          //this is approximate for the case when there are constraints
		  //- there may be some double counting
		  int getApproxNumTuples() {
               int cnt = 0;
			   for(int i=0;i<hyperCubes_->size();i++) {
					cnt += (*hyperCubes_)[i]->getNumTuples();
			   }
			   return cnt;
		  }

		  const ClauseCounter * getClauseCounter() { return clauseCounter_;}
          void clearClauseCounter() { delete clauseCounter_; clauseCounter_ = new ClauseCounter();}

		  void addHyperCube(HyperCube * hyperCube, int predId) {
			   hyperCubes_->append(hyperCube);
               HyperCubeToSuperPredMap * hyperCubeToSuperPred;
			   hyperCubeToSuperPred = (*(SuperPred::hyperCubeToSuperPredArr_))[predId];
               (*hyperCubeToSuperPred)[hyperCube] = this;
			   //int tmp;
               //(*constantsToSuperPred)[constants] = new SuperPred(tmp);
		  }


          void printPredicates(Domain * const & domain, 
							   ostream & out) {
			  Array<IntArray *> * tuples;
			  IntArray *tuple;
			  HyperCube *hyperCube;
			  //const PredicateTemplate *pt;
			  Predicate *gndPred;
			  for(int hindex=0;hindex<hyperCubes_->size();hindex++) {
				  hyperCube = (*hyperCubes_)[hindex];
			      tuples = hyperCube->getTuples();
				  for(int tindex=0;tindex<tuples->size();tindex++) {
                     tuple = (*tuples)[tindex];
				     tuple->removeItem(0);
					 
					 gndPred = domain->getPredicate(tuple,predId_);
					 /*
					  pt = domain->getPredicateTemplate(predId_);
					  gndPred = getGroundPredicate(pt,tuple);
					 */
					 gndPred->printWithStrVar(out,domain);
					 out<<endl;
					 delete gndPred;
					 delete tuple;
				   }
				   tuples->clear();
				   delete tuples;
			   }
		   }


		  //static functions
		  //
		  static int getSuperPredCount(int predId){return (*superPredsArr_)[predId]->size();}
         
		  //get all the superpreds for the given predId
		  static Array<SuperPred *> * getSuperPreds(int predId) { return (*superPredsArr_)[predId];}
		 
		  static void clearClauseCounters(int predCnt) {
			   Array<SuperPred *> *superPreds;
			   for(int predId=0;predId<predCnt;predId++) {
			   superPreds = (*superPredsArr_)[predId];
			   for(int i=0;i<superPreds->size();i++) {
				 (*superPreds)[i]->clearClauseCounter();
			    }
			   }
		  }

		  static void clear(int predCnt) {
			   //if first time access then, need to allocate memory
			   if(isFirstAccessForStaticVars_) {
					superPredsArr_ = new Array<Array<SuperPred *>*>();
					hyperCubeToSuperPredArr_ = new Array<HyperCubeToSuperPredMap *>();
					for(int predId=0;predId<predCnt;predId++) {
					   //cout<<"came here for creaing the hash maps..."<<endl;
					   hyperCubeToSuperPredArr_->append(new HyperCubeToSuperPredMap());
					   superPredsArr_->append(new Array<SuperPred *>());
					}
					isFirstAccessForStaticVars_ = false;
			   } else {
					//this is not the first time access - need to clean up
					Array<SuperPred *> *superPreds;
					HyperCubeToSuperPredMap * hyperCubeToSuperPred;
					HyperCubeToSuperPredMap::iterator hcToSpItr;
					for(int predId=0;predId<predCnt;predId++) {
					   
					   //first delete super preds	 
					   superPreds = (*superPredsArr_)[predId];
					   for(int i=0;i<superPreds->size();i++) {
							delete (*superPreds)[i];
					   }

					   //now delete the constants
					   Array<HyperCube *> keysArr;
					   hyperCubeToSuperPred = (*hyperCubeToSuperPredArr_)[predId];
		               keysArr.clear();
					   for(hcToSpItr = hyperCubeToSuperPred->begin();
			               hcToSpItr != hyperCubeToSuperPred->end();
			               hcToSpItr++) {
							keysArr.append(hcToSpItr->first);
					   }
					   for(int i=0;i<keysArr.size();i++) {
							delete keysArr[i];
					   }
				
					   //reinitialize
					   superPreds->clear();
					   hyperCubeToSuperPred->clear();
			       }
			    }
		  }
		  
	 static SuperPred * getSuperPred(HyperCube * const & hyperCube, int & predId) {
			HyperCubeToSuperPredMap * hyperCubeToSuperPred;
			HyperCubeToSuperPredMap::iterator hcToSpItr;
			if(isFirstAccessForStaticVars_) {
				 return NULL;
			}
			assert(hyperCubeToSuperPredArr_);
			hyperCubeToSuperPred = (*hyperCubeToSuperPredArr_)[predId];
            hcToSpItr = hyperCubeToSuperPred->find(hyperCube); 
			if(hcToSpItr == hyperCubeToSuperPred->end()) {
				 cout<<"Problem here...!!"<<endl;
				 cout<<"Can't find the hypercube "<<endl;
				 hyperCube->print(cout)<<endl;
				 cout<<endl;
			}
			SuperPred *superPred = hcToSpItr->second;
			assert(superPred != NULL);
			return superPred;
     }

	 static int getSuperPredId(HyperCube * const & hyperCube, int & predId) {
			  SuperPred * superPred = getSuperPred(hyperCube,predId);
			  if(superPred) {
			   return superPred->getSuperPredId();
			  } else {
			   return -1;
			  }
	 }

	 static HyperCubeRefinement *getHyperCubeRefinement(int predId) {
		  if(!hcRefinementArr_)
			   return NULL;
		  return (*hcRefinementArr_)[predId];
	 }

	 static Array<HyperCubeRefinement *> *getHyperCubeRefinement() {
		  return hcRefinementArr_;
	 }

     //delete the hypercube refinement array
	 static void deleteHyperCubeRefinement() {
	  if(!hcRefinementArr_)
		   return;
	  
	  for(int i=0;i<hcRefinementArr_->size();i++) {
	   delete (*hcRefinementArr_)[i];
      }
      delete hcRefinementArr_;
      SuperPred::setHyperCubeRefinement(NULL);
	 }

	 static void setHyperCubeRefinement(Array<HyperCubeRefinement *> * const & 
			                            newHCRefinementArr) {
		  hcRefinementArr_ = newHCRefinementArr;
	 }

	 private:
		  int predId_;
		  int superPredId_;
		  int parentSuperPredId_;
		  Array<HyperCube *> *hyperCubes_;
		  ClauseCounter *clauseCounter_;

		  static bool isFirstAccessForStaticVars_;
		  static Array<HyperCubeToSuperPredMap*>* hyperCubeToSuperPredArr_;
		  static Array<Array<SuperPred *>*>* superPredsArr_;
		  static Array<HyperCubeRefinement *> * hcRefinementArr_;
};


class PredHCAuxInfo {
	 public:
	 PredHCAuxInfo(ClauseCounter * const & cc, int superPredId) {
		  cc_ = cc;
		  superPredId_ = superPredId;
	 }
     
	 ~PredHCAuxInfo() {
		  delete cc_;
	 }

	 int getSuperPredId() { return superPredId_;}
	 ClauseCounter * getClauseCounter() {return cc_;}

	 private:
		  int superPredId_;
		  ClauseCounter *cc_;
};

typedef hash_map<ClauseCounter*, SuperPred*, HashClauseCounter, EqualClauseCounter> ClauseCounterToSuperPredMap;
typedef hash_map<HyperCube*, PredHCAuxInfo*, HashHyperCube, EqualHyperCube> HyperCubeToPredHCAuxInfoMap;

#endif

