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
#ifndef _PARTIAL_GND_HC_MAY_2008
#define _PARTIAL_GND_HC_MAY_2008

#include <ext/hash_map>
#include <ostream>

#include "util.h"
#include "mrf.h"
#include "array.h"
#include "hashint.h"
#include "hypercubeoperations.h"

class PartialGndHCTuple {

	 public:
		  PartialGndHCTuple(Clause * const & clause, bool useCT) {
			   
               clause_ = clause;
               useCT_ = useCT;

			   int varCnt = clause->getNumVariables();
			   int predCnt = clause->getNumPredicates();

			   hyperCube_ = new HyperCube(varCnt);
               
			   isGndPreds_ = new Array<bool>(predCnt,false);
			   joinHyperCubes_ = new Array<HyperCube *>(predCnt,NULL);
		       constrainedHyperCubes_ = NULL;
		  }	   

          PartialGndHCTuple(PartialGndHCTuple * const & pgTuple) {
		       clause_ = pgTuple->getClause();
		       useCT_ = pgTuple->isUseCT();
			   hyperCube_ = new HyperCube(pgTuple->getHyperCube());
			   isGndPreds_ = new Array<bool>(*pgTuple->getIsGndPreds());
			   joinHyperCubes_ = new Array<HyperCube *>(*pgTuple->getJoinHyperCubes());
		       constrainedHyperCubes_ = NULL;
		  }

		  ~PartialGndHCTuple() {
               delete hyperCube_;
			   delete isGndPreds_;
			   delete joinHyperCubes_;
			   if(constrainedHyperCubes_)
					delete constrainedHyperCubes_;
		  }

	  //whether to use the constrained representation
	  bool isUseCT() {
		   return useCT_;
	  }

       // a var is gnded if the corresponding set of var
	   // constants is not null
	   bool isVarGnded(int varId) {
		  return (hyperCube_->getVarConstants(varId) != NULL);
	   } 
	 
	 //is the clause now completely grounded
     bool isGnded() {
		  int predCnt = isGndPreds_->size();
		  for(int i=0;i<predCnt;i++) {
              if(!(*isGndPreds_)[i])
				   return false;
		  }
		  return true;
	 }

	 Clause * getClause() { return clause_;}
     HyperCube * getHyperCube() { return hyperCube_;}
	 Array<HyperCube *> * getJoinHyperCubes() { return joinHyperCubes_;}
	 
	 Array<HyperCube *> *getConstrainedHyperCubes() { return constrainedHyperCubes_;}

	 Array<bool> * getIsGndPreds(){return isGndPreds_;}

	 IntArray * getVarConstants(int varId) {
		  return hyperCube_->getVarConstants(varId);
	 }
     
	 //get the index of first ungrounded predicate
	 int getFirstUnGndPredIndex() {
		  int predCnt = isGndPreds_->size();
		  for(int pindex=0;pindex<predCnt;pindex++) {
			   if(!(*isGndPreds_)[pindex])
					return pindex;
		  }
		  return -1;
	 }

	 //extract the hypercube corresponding to a given predicate
	 HyperCube * getPredicateHyperCube(HyperCube * const & hc, int pindex) {
		  Predicate *pred = clause_->getPredicate(pindex);
		  int predTermCnt = pred->getNumTerms();
          HyperCube *phyperCube = new HyperCube(predTermCnt);
          IntArray * varConstants;

		  if(useCT_) {
		   hc->updateToImplicitTiedRepresentation();
		  }

		  for(int termno=0; termno<predTermCnt;termno++) {
           const Term *term = pred->getTerm(termno);
           int varId = -term->getId();
		   
		   //this variable should be grounded
		   assert(isVarGnded(varId));

		   //copies are created in the end
		   //varConstants = new IntArray(*(hc->getVarConstants(varId)));
		   varConstants = hc->getVarConstants(varId);
		   phyperCube->setVarConstants(varConstants, termno+1);
		  }
		 
		  if(useCT_) {
		   hc->updateToExplicitTiedRepresentation();
		   phyperCube->updateToExplicitTiedRepresentation();
		  }

		  //now create the copies for each of the varConstants
		  for(int termno=0; termno<predTermCnt;termno++) {
		   int termIndex = termno+1;
		   int refTermIndex = phyperCube->getReferenceVarId(termIndex);
		   //this is a reference variable - nothing to do
		   
		   if(!useCT_) {
				assert(refTermIndex == termIndex);
		   }

		   if(refTermIndex != termIndex)
				continue;
		   varConstants = new IntArray(*(phyperCube->getVarConstants(termIndex)));
		   phyperCube->setVarConstants(varConstants,termIndex);
		  }
		  return phyperCube;
	 }

	 //get the id of the variable this variable refers to (or return varId
	 //if it does not refer to anything)
     int getReferenceVarId(int varId) {
		  return hyperCube_->getReferenceVarId(varId);
	 }

	 void join(HyperCube * const & hyperCube,
			   int pindex) {
          Predicate *pred = clause_->getPredicate(pindex);
          IntArray * varConstants, *refVarConstants;
		 
		  if(useCT_) {
		   /*
		   cout<<"Before joining, hypercube is"<<endl;
		   hyperCube_->print(cout);
 		   cout<<endl;
          */
		   hyperCube_->updateToImplicitTiedRepresentation();
		  }
		  
		  /*
		  cout<<"Join hypercube is "<<endl;
		  hyperCube->print(cout);
		  cout<<endl;
          */

		  for(int termno=0; termno<pred->getNumTerms();termno++) {
           const Term *term = pred->getTerm(termno);
           int varId = -term->getId();
		   int termIndex = termno+1;
		   int refTermIndex = hyperCube->getReferenceVarId(termIndex);
		   //cout<<"var Id = "<<varId<<endl;

		   if(!useCT_) {
				assert(refTermIndex == termIndex);
		   } 
		   //If there are variables tied in the input joint hypercube, tie them in 
		   //the partial hypercube representation as well
		   if(refTermIndex != termIndex) {
               int refTermNo = refTermIndex - 1;
			   int refVarId = -(pred->getTerm(refTermNo))->getId();
			   //cout<<"ref var Id = "<<refVarId<<endl;
			   assert(isVarGnded(refVarId));
			   varConstants = hyperCube_->getVarConstants(varId);
			   refVarConstants = hyperCube_->getVarConstants(refVarId);
			   if(varConstants == refVarConstants)
					continue;
			   if(!varConstants) {
                hyperCube_->setVarConstants(refVarConstants,varId);
			   } else {
			    hyperCube_->intersect(varConstants,refVarId);
			    refVarConstants = hyperCube_->getVarConstants(refVarId);
				hyperCube_->setVarConstantsAndDeleteOld(refVarConstants,varId);
			   }
		    }

		   if(isVarGnded(varId)) { 
				//cout<<"joined with existing.."<<endl;
		        varConstants = hyperCube->getVarConstants(refTermIndex);
		       	hyperCube_->intersect(varConstants,varId);
		   } else {
				varConstants = new IntArray(*hyperCube->getVarConstants(refTermIndex));
		       	hyperCube_->setVarConstants(varConstants,varId);
		   }
		  }

          if(useCT_) {
		   /*
		   cout<<"Before updating to explicit rep, hyperCube is"<<endl;
		   hyperCube_->print(cout);
		   cout<<endl;
		   */
		   hyperCube_->updateToExplicitTiedRepresentation();
		  }

		  (*isGndPreds_)[pindex] = true;
		  (*joinHyperCubes_)[pindex] = hyperCube;
		  
		  /*
		  cout<<"After joining, the hypercube is"<<endl;
		  hyperCube_->print(cout);
		  cout<<endl;
		  */
	 }

	 //create the constrained hypercubes - split of the origHyperCube
     void createConstrainedHyperCubes(IntArrayHashArray * const & inpNeqConstraints) {
	  assert(constrainedHyperCubes_ == NULL);
	  
	  HyperCube *origHyperCube = new HyperCube(hyperCube_);

	  /*
	  cout<<"Processing hypercube = "<<endl;
	  origHyperCube->print(cout);
	  cout<<endl;
	  */
	  constrainedHyperCubes_ = new Array<HyperCube *>();
	  if(!useCT_) {
		   constrainedHyperCubes_->append(origHyperCube);
		   return;
	  }
	  
	  IntArrayHashArray * neqConstraints = hyperCube_->getRelevantConstraints(inpNeqConstraints);
	  Array<HyperCube *> *minusHyperCubes;
	  HyperCube *hyperCube, minusHyperCube, *intersectHyperCube;
      IntArray *constraint;
      bool isReadyToAdd;
	  
	  Array<HyperCube *> *hcStack = new Array<HyperCube *>();
	  hcStack->append(origHyperCube);

	  while(hcStack->size() > 0) {
		   hyperCube = hcStack->removeLastItem();
		   isReadyToAdd = true;
		   for(int cno=0;cno<neqConstraints->size();cno++) {
              
			  constraint = (*neqConstraints)[cno];
			  int varId1 = -(*constraint)[0];
			  int varId2 = -(*constraint)[1];
	
			  /*
			  cout<<"Processing hypercube.. "<<endl;
			  hyperCube->print(cout);
			  cout<<endl;
              */
			  
			  if(hyperCube->hasSameVarConstants(varId1,varId2))
				   continue;
              if(hyperCube->hasDisjointVarConstants(varId1,varId2))
				   continue;
              isReadyToAdd = false;
			  
			  intersectHyperCube = new HyperCube(hyperCube);
			 
			  /*
			  cout<<"came till 1 HyperCube = "<<endl;
			  hyperCube->print(cout);
			  cout<<endl;
			  */

			  intersectHyperCube->intersect(varId1,varId2);
			  
			  /*
			  cout<<"Intersect HyperCube = "<<endl;
			  intersectHyperCube->print(cout);
			  cout<<endl;
			  */

			  minusHyperCubes = hyperCube->getMinus(intersectHyperCube);
			  
			  /*
			  cout<<"came till 2 HyperCube = "<<endl;
			  hyperCube->print(cout);
			  cout<<endl;
			  cout<<"Minus hypercubes are :"<<endl;
			  */
			  
			  for(int i=0;i<minusHyperCubes->size();i++) {
			    //(*minusHyperCubes)[i]->print(cout);
			    //cout<<endl;
				hcStack->append((*minusHyperCubes)[i]);
			  }
			  delete minusHyperCubes;
			  hcStack->append(intersectHyperCube);
			 
			  /* 
			  cout<<"came till 3 HyperCube = "<<endl;
			  hyperCube->print(cout);
			  cout<<endl; 

			  cout<<"came till end. HyperCube = "<<endl;
			  hyperCube->print(cout);
			  cout<<endl;
		      */
			  
			  delete hyperCube;
			  break;
		   }
		   //if came here, all the constraints are satisified
           if(isReadyToAdd) {
		    constrainedHyperCubes_->append(hyperCube);
		    /*
			cout<<"Adding hyperCube"<<endl;
			hyperCube->print(cout);
			cout<<endl;
		   */
		  }
	  }
      
	  neqConstraints->deleteItemsAndClear();
	  delete neqConstraints;
	  delete hcStack;
}
				   

	 //update the hypercuberefinement object - add the mapping from
	 //joinHyperCube (the actual hypercube which came in during join) to 
	 //the predHypercube (hypercube obtained after joining)
	 void updateHyperCubeRefinement(Array<HyperCubeRefinement *> * const & hyperCubeRefinementArr,
			                        IntHashArray * const & unknownPredIndices) {
          assert(isGnded());
		  HyperCube *joinHyperCube, *predHyperCube, *baseHyperCube;
		  HyperCubeRefinement *hyperCubeRefinement;
		  Predicate *pred;
		  int predCnt = clause_->getNumPredicates();
		
		  HyperCube *hc;
		  for(int hindex=0;hindex<constrainedHyperCubes_->size();hindex++) {
              hc = (*constrainedHyperCubes_)[hindex];
		     //cout<<"size of unknown pred ids = "<<unknownPredIndices->size()<<endl;
		     for(int pindex=0;pindex<predCnt;pindex++) {
			 
			  //need to build the refinement array only for unknown preds
			  if(unknownPredIndices->find(pindex) < 0)
				  continue;

			  pred = clause_->getPredicate(pindex);
			  int predId =  pred->getId();
			  hyperCubeRefinement = (*hyperCubeRefinementArr)[predId];

			  joinHyperCube = (*joinHyperCubes_)[pindex];
			  //baseHyperCube = joinHyperCube;
			 
			  baseHyperCube = new HyperCube(joinHyperCube);
              predHyperCube = getPredicateHyperCube(hc,pindex);

			  //the method is also responsible for delete the input hypercubes
			  //if they are already present in the refinement
			  hyperCubeRefinement->addSubsetHyperCube(baseHyperCube,predHyperCube);
		  }
		}
	 }

 private:
  Clause *clause_;
  HyperCube *hyperCube_;
  Array<HyperCube *> * joinHyperCubes_;
  Array<HyperCube *> * constrainedHyperCubes_;
  Array<bool> *isGndPreds_;
  //whether to use the constrained representation
  bool useCT_;
};


#endif

