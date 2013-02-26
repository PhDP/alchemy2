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
#include "util.h"
#include "mrf.h"
#include "superclause.h"
#include "superpred.h"

//this static variable is used to intialize all the static vars appropriately
bool SuperPred::isFirstAccessForStaticVars_ = true;
Array<HyperCubeToSuperPredMap *> * SuperPred::hyperCubeToSuperPredArr_ = NULL;
Array<Array<SuperPred *>*>* SuperPred::superPredsArr_ = NULL;
Array<HyperCubeRefinement *>* SuperPred::hcRefinementArr_ = NULL;

	
//initialize the superPreds - one superpred for each predicate
void initSuperPreds(Array<Array<SuperClause *> *> * const & superClausesArr,
		            Domain * const & domain) {

	  int domainPredCnt = domain->getNumPredicates();
	  Array<HyperCubeHashArray *> * hyperCubesArr = new Array<HyperCubeHashArray *>();
      for(int predId=0;predId<domainPredCnt;predId++) {
		   hyperCubesArr->append(new HyperCubeHashArray());
	  }
      
	  SuperPred::clear(domainPredCnt);
	  
	  HyperCubeRefinement *hcRefinement;
	  Array<SuperClause *> * superClauses;
	  HyperCubeHashArray * hyperCubes;
	  HyperCube *phyperCube, *baseHyperCube;
	  Clause* clause;
	  Predicate *pred;
	  SuperClause *superClause;
	  ClauseCounter  *clauseCounter;
	  int superClauseCnt;
	  
	  Array<SuperPred *> *superPreds;
	  SuperPred *superPred;

      //iterate over all the super clauses
	  for(int arrIndex=0;arrIndex<superClausesArr->size();arrIndex++) {
	   superClauses = (*superClausesArr)[arrIndex];
	   superClauseCnt = superClauses->size();
	  
	   for(int scindex=0;scindex<superClauseCnt;scindex++) {
		  superClause = (*superClauses)[scindex];
		  clause = superClause->getClause();
		  int predCnt = clause->getNumPredicates();

		  //for each hypercube, extract the predicate hypercube
		  int numHyperCubes = superClause->getNumHyperCubes();
		  for(int hindex=0;hindex<numHyperCubes;hindex++) {
			   for(int pindex=0;pindex<predCnt;pindex++) {
					pred = clause->getPredicate(pindex);
					phyperCube = superClause->getPredicateHyperCube(hindex, pred);
	
					int predId = domain->getPredicateId(pred->getName());
					hcRefinement = SuperPred::getHyperCubeRefinement(predId);
					if(hcRefinement) {
					 baseHyperCube = hcRefinement->getBaseHyperCube(phyperCube);
					 //create a new copy
					 baseHyperCube = new HyperCube(baseHyperCube);
					 delete phyperCube;
					} else {
					 baseHyperCube = phyperCube;
					}
                    
					hyperCubes = (*hyperCubesArr)[predId];
					//if this hypercube has already not been added
                    if(hyperCubes->append(baseHyperCube) < 0) {
						 delete baseHyperCube;
						 continue;
					}

					superPreds = SuperPred::getSuperPreds(predId);
					assert((superPreds->size() == 0) || (superPreds->size() == 1));
					
					if(superPreds->size() == 0) {
				       clauseCounter = new ClauseCounter();
					   //this is automatically added to the list of superPreds for predId
					   superPred = new SuperPred(predId, clauseCounter, -1);
					} else {
					   superPred = (*superPreds)[0];
					}
			        //superPred->addHyperCube(new HyperCube(baseHyperCube), predId);
			        superPred->addHyperCube(baseHyperCube, predId);
			   }
		  }
	   }
	  }
	  
	  hyperCubesArr->deleteItemsAndClear();
	  delete hyperCubesArr;
}



//create the super preds given the current set of super features
void createSuperPreds(Array<Array<SuperClause *> *> * const & superClausesArr,
		              Domain * const & domain,
					  bool isApproximateSuperPreds) {

	 //cout<<"*****************************************************"<<endl;
	 //cout<<"Creating Super Preds.........."<<endl;
	 //cout<<"*****************************************************"<<endl;

	 IntHashArray seenPredIds;
	 Array<SuperClause *> * superClauses;

	 int domainPredCnt = domain->getNumPredicates();
	 Array<HyperCubeToPredHCAuxInfoMap* > * hcToAuxInfoArr = NULL;
	 
	 if(!isApproximateSuperPreds) {
      hcToAuxInfoArr = new Array<HyperCubeToPredHCAuxInfoMap*>();
	  hcToAuxInfoArr->growToSize(domainPredCnt);
	  //initialize the mappings from pred constants to their clause counts
	  for(int i=0;i<domainPredCnt;i++) {
		  (*hcToAuxInfoArr)[i] = new HyperCubeToPredHCAuxInfoMap();
	  }
	 }
	 
	 HyperCubeToPredHCAuxInfoMap * hcToAuxInfo;
	 HyperCubeToPredHCAuxInfoMap::iterator hcToAuxInfoItr;
     
	 HyperCubeRefinement *hcRefinement;
     Array<HyperCube *> * refinedHyperCubes = new Array<HyperCube *>();

	 HyperCube *phyperCube, *baseHyperCube, *refinedHyperCube;
	 Clause* clause;
	 Predicate *pred;
	 SuperClause *superClause;
	 PredHCAuxInfo *auxInfo;
	 ClauseCounter  *clauseCounter;
	 int superClauseCnt;
    
	 Array<SuperPred *> *superPreds;
	 SuperPred *superPred;
	 
	 //whether this run needs to consider the refined hypercube creation process
	 bool useRefinedHC;
	 
	 Array<HyperCubeRefinement *> *hcRefinementArr;
	 hcRefinementArr = SuperPred::getHyperCubeRefinement();			
	 if(hcRefinementArr) 
	    useRefinedHC = true;
	 else
		useRefinedHC = false;
	 
	 SuperPred::clearClauseCounters(domainPredCnt);

	 //whether refined hypercube creation process has converged
	 bool refinedHCConverged = true;
    
	 int prevCnt = 0;
	 //iterate over all the super clauses
	 for(int arrIndex=0;arrIndex<superClausesArr->size();arrIndex++) {
	  superClauses = (*superClausesArr)[arrIndex];
	  superClauseCnt = superClauses->size();
	  
	  for(int scindex=0;scindex<superClauseCnt;scindex++) {
		  superClause = (*superClauses)[scindex];
		  clause = superClause->getClause();
		  int predCnt = clause->getNumPredicates();

		  //for each hypercube, extract the predicate hypercube
		  int numHyperCubes = superClause->getNumHyperCubes();
		
		  /*
		  cout<<endl<<endl;
          cout<<"--------------------------------------------------------------"<<endl;
		  cout<<scindex<<" => Clause: ";
		  clause->printWithoutWt(cout,domain);
		  cout<<endl<<"Num tuples is = "<<numTuples<<endl;
		  cout<<endl<<endl;
          */
		  for(int hindex=0;hindex<numHyperCubes;hindex++) {
			   
			   //number of times this hypercube appears in the superclause
			   double hcnt = superClause->getHyperCubeCount(hindex);
               //this will be the count of impilcit constant tuples joining with
			   //a particular predicate
			   
			   /*
			   HyperCube *hc = superClause->getHyperCube(hindex);
			   cout<<"HyperCube is :"<<endl;
			   hc->print(cout);
			   cout<<endl;
			  */

			   /*
		       cout<<"hcnt is = "<<hcnt<<endl; 
			   */

			   int implicitCnt;
			   for(int pindex=0;pindex<predCnt;pindex++) {
					pred = clause->getPredicate(pindex);
					phyperCube = superClause->getPredicateHyperCube(hindex, pred);
	
					//Change from the vanilla (trivial hypercube) case
					//is that we now iterate over the refined subsets of
					//the phypercubes and store the counts for them
					int predId = domain->getPredicateId(pred->getName());
					hcRefinement = SuperPred::getHyperCubeRefinement(predId);

					/*
					cout<<"predicate hypercube is"<<endl;
					phyperCube->print(cout);
					cout<<endl;
					*/
					refinedHyperCubes->clear();
					if(useRefinedHC) {
					 (*refinedHyperCubes) = (*(hcRefinement->getRefinedHyperCubes(phyperCube)));
					 baseHyperCube = hcRefinement->getBaseHyperCube(phyperCube);
					 
					 //hypercube creation has not converged as long there is
					 //some predicate hypercube for which there is more than
					 //one refined hypercube
					 if(refinedHyperCubes->size() > 1) {
						refinedHCConverged = false;
					 }
					 delete phyperCube;
					
					} else {
					 baseHyperCube = phyperCube;
					 refinedHyperCubes->append(phyperCube);
					}
					
					assert(baseHyperCube);
					if(!baseHyperCube) {
						 cout<<"Problem here.!!"<<endl;
						 cout<<"No base hypercube!"<<endl;
					}
			
					implicitCnt = superClause->getNumTuplesJoiningWithHyperCube(hindex, pred);
					seenPredIds.append(predId);
					//hcToAuxInfo = (*hcToAuxInfoArr)[predId];

					for(int rhindex=0;rhindex<refinedHyperCubes->size();rhindex++) {
					  
					  //refinedHyperCube = (*refinedHyperCubes)[rhindex];
					  if(useRefinedHC)
					   refinedHyperCube = new HyperCube((*refinedHyperCubes)[rhindex]);
					  else
					   refinedHyperCube = (*refinedHyperCubes)[rhindex];

					  //Added to do the approximate construction of superPreds
					  if(isApproximateSuperPreds) {
						superPred = SuperPred::getSuperPred(baseHyperCube, predId);
						clauseCounter = (ClauseCounter *)superPred->getClauseCounter();
						int hcNumTuples = refinedHyperCube->getNumTuples();
					    clauseCounter->incrementCount(scindex+prevCnt,pindex,predCnt,hcnt*implicitCnt*hcNumTuples);
					    delete refinedHyperCube;
					  } else {
					   
					   hcToAuxInfo = (*hcToAuxInfoArr)[predId];
					   hcToAuxInfoItr = hcToAuxInfo->find(refinedHyperCube);
			
					   /*
					   cout<<"refined hypercube ("<<pindex<<") => "<<endl;
					   refinedHyperCube->print(cout);
					   cout<<" ** cnt = "<<implicitCnt*hcnt<<endl;
				      */
					   /*
					   if(implicitCnt > 1) {
					    cout<<"In: ";
					    (superClause->getHyperCube(hindex))->print();
					    cout<<"** Pred: ";
					    refinedHyperCube->print(cout);
					    cout<<endl<<" implicit count = "<<implicitCnt<<endl;
					   }*/
	               
					  if(hcToAuxInfoItr == hcToAuxInfo->end()) {
						 //int superPredId = SuperPred::getSuperPredId(refinedHyperCube, predId);
						 int superPredId = SuperPred::getSuperPredId(baseHyperCube, predId);
						 clauseCounter = new ClauseCounter();
                         auxInfo = new PredHCAuxInfo(clauseCounter,superPredId);
						 (*hcToAuxInfo)[refinedHyperCube] = auxInfo;
					   } else {
						 //if(!useRefinedHC)
                         delete refinedHyperCube;
						 auxInfo = hcToAuxInfoItr->second;
						 clauseCounter = auxInfo->getClauseCounter();
					  }
					  clauseCounter->incrementCount(scindex+prevCnt,pindex,predCnt,hcnt*implicitCnt);
				   }
				}
			   }
		  }
	   }
	  prevCnt += superClauseCnt;
	 }

	 refinedHyperCubes->clear();
	 delete refinedHyperCubes;
	 
	 if(isApproximateSuperPreds) {
       for(int predId=0;predId<domainPredCnt;predId++) {
		    if(seenPredIds.find(predId) < 0)
			   continue;
			superPreds = SuperPred::getSuperPreds(predId);
			for(int i=0;i<superPreds->size();i++) {
				 superPred = (*superPreds)[i];
				 clauseCounter = (ClauseCounter *)superPred->getClauseCounter();
				 int numTuples = superPred->getApproxNumTuples();
				 clauseCounter->multiplyCounts(1.0/numTuples);
			}
	   }
	   return;	  
	 }

	 //we do not need the hypercube refinement stuff anymore
	 //- it has converged
	 if(useRefinedHC && refinedHCConverged) {
      cout<<"Refined hypercube creation process has converged...! :-)"<<endl;
	  SuperPred::deleteHyperCubeRefinement();
	  /*
	   for(int i=0;i<hcRefinementArr->size();i++) {
	   delete (*hcRefinementArr)[i];
       }
       delete hcRefinementArr;
       SuperPred::setHyperCubeRefinement(NULL);
	  */
	 }

	 /*now, cluster together all the hypercubes having the same count*/
	 ClauseCounterToSuperPredMap * clauseCounterToSuperPred;
	 ClauseCounterToSuperPredMap::iterator ocToSpItr;
     //SuperPred *superPred;

	 int parentSuperPredId;

	 //clear all the previously stored super preds and their counts
	 SuperPred::clear(domainPredCnt);
     
	 for(int predId=0;predId<domainPredCnt;predId++) {
         //don't need to worry about the preds which did not appear in any of the superclauses
		 if(seenPredIds.find(predId) < 0)
			  continue;

		 //cout<<"**** Doing it for Predicate ***"<<domain->getPredicateName(predId)<<endl;
		 clauseCounterToSuperPred = new ClauseCounterToSuperPredMap(); 
 		 hcToAuxInfo = (*hcToAuxInfoArr)[predId];
		 
		 for(hcToAuxInfoItr = hcToAuxInfo->begin();
			 hcToAuxInfoItr != hcToAuxInfo->end();
			 hcToAuxInfoItr++) {
			 
			 phyperCube = hcToAuxInfoItr->first;
			 auxInfo = hcToAuxInfoItr->second;
			 parentSuperPredId = auxInfo->getSuperPredId();
			 clauseCounter = auxInfo->getClauseCounter();

			 
			/* 
			 cout<<"For the hypercube=>"<<endl;
			 phyperCube->print(cout);
			 cout<<"the clauseCounter object is "<<endl;
			 clauseCounter->print(cout);
			 cout<<endl<<endl;
             */

			 ocToSpItr = clauseCounterToSuperPred->find(clauseCounter);
			 if(ocToSpItr == clauseCounterToSuperPred->end()) {
				  //cout<<"Created a new super pred.."<<endl;
				  superPred = new SuperPred(predId, clauseCounter, parentSuperPredId);
				  (*clauseCounterToSuperPred)[clauseCounter] = superPred;
			 } else {
                  delete clauseCounter;
				  superPred = ocToSpItr->second;
			 }
			 /*
			 cout<<"Hyper cube is=>"<<endl;
			 phyperCube->print(cout);
			 cout<<endl;
			 */
			 superPred->addHyperCube(phyperCube, predId);
		 }
         //cout<<"Number of superpreds for this predicate = "+clauseCounterToSuperPred->size()<<endl;
		 //cout<<"****************************************************"<<endl<<endl;
         clauseCounterToSuperPred->clear();
		 delete clauseCounterToSuperPred;
	 }
	 
	 // clean up
	 for(int predId=0;predId<domainPredCnt;predId++) {
		 hcToAuxInfo = (*hcToAuxInfoArr)[predId];
		 hcToAuxInfo->clear();
		 delete hcToAuxInfo;
	 }

	 delete hcToAuxInfoArr;
}

