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
#ifndef SUPERCLAUSE_H_NOV_2007
#define SUPERCLAUSE_H_NOV_2007

#include "util.h"
#include "mrf.h"
#include "array.h"
#include "hashint.h"
#include <ext/hash_set>
#include "variable.h"
#include "hypercube.h"
#include "hypercubeoperations.h"
#include "inferutil.h"

using namespace std;
using namespace __gnu_cxx;

class SuperClause
{
  public:
          
    SuperClause(Clause * const & clause, 
                IntArrayHashArray * const & neqConstraints,
                Array<int> * const & varIdToCanonicalVarId,
                bool useCT, double outputWt)
    {
      int parentSuperClauseId = -1;
      outputWt_ = outputWt;
      init(clause,neqConstraints,varIdToCanonicalVarId,useCT,parentSuperClauseId);
    }

          SuperClause(Clause * const & clause, 
                      IntArrayHashArray * const & neqConstraints,
                      Array<int> * const & varIdToCanonicalVarId,
                      bool useCT, 
                      int parentSuperClauseId, double outputWt)
    {
    	outputWt_ = outputWt;
               init(clause,neqConstraints,varIdToCanonicalVarId,useCT,parentSuperClauseId);
          }

          void init(Clause * const & clause, 
                      IntArrayHashArray * const & neqConstraints,
                      Array<int> * const & varIdToCanonicalVarId,
                      bool useCT,
                      int parentSuperClauseId) {
               clause_ = clause;
			   neqConstraints_ = neqConstraints;
               varIdToCanonicalVarId_ = varIdToCanonicalVarId;
               useCT_ = useCT;
               
			   hyperCubes_ = new HyperCubeHashArray();
               hyperCubeCnts_ = new Array<double>();
			   
               superClauseId_ = superClauseIndex__++; 
               parentSuperClauseId_ = parentSuperClauseId;
          }
                    
          ~SuperClause() {
               delete hyperCubes_;
               delete hyperCubeCnts_;
          }
          
          SuperClause * createSuperClauseFromTemplate() {
               //this becomes the parent super clause
               return new SuperClause(clause_, neqConstraints_, varIdToCanonicalVarId_,
                                      useCT_,superClauseId_);
          }

		  IntArrayHashArray * getNeqConstraints() { return neqConstraints_;}

          int getSuperClauseId() { return superClauseId_;}
          
          int getParentSuperClauseId() { return parentSuperClauseId_;}

          Clause * getClause(){ return clause_;}

          double getHyperCubeCount(int hindex) {return (*hyperCubeCnts_)[hindex];}

          int getNumHyperCubes(){ return hyperCubes_->size();}

          int getHyperCubeIndex(HyperCube * const & hyperCube) {return hyperCubes_->find(hyperCube);}

          HyperCubeHashArray *getHyperCubes() {return hyperCubes_;}
          
		  HyperCube * getHyperCube(int hindex) {return (*hyperCubes_)[hindex];}

          Array<int> * getVarIdToCanonicalVarId(){return varIdToCanonicalVarId_;}

          bool isUseCT(){ return useCT_;}

          //increments the count of the given HyperCube. Assumes that
          //hyperCube is already present
          void incrementHyperCubeCount(HyperCube * const & hyperCube, double cnt) {
             int index = hyperCubes_->find(hyperCube);
             (*hyperCubeCnts_)[index] += cnt;
          }

          //add the hyperCube - returns true if the addition was
          //successful i.e. if the hyperCube was not already present
          bool addHyperCube(HyperCube * const & hyperCube) {
               int index = hyperCubes_->find(hyperCube);
               //this tuple does not exist already then add it
               if(index < 0) {
                 hyperCubes_->append(hyperCube);
                 hyperCubeCnts_->append(0.0);
                 return true;
               } else {
                    return false;
               }
          }

          //create a hypercube with the given constants and increment the count
          void addNewConstantsAndIncrementCount(Array<int> * const & constants, double cnt) {
           bool addedNew = false;
           HyperCube *newHyperCube = new HyperCube(constants);
           //add this hyperCube to the superclause
           addedNew = addHyperCube(newHyperCube);   
           incrementHyperCubeCount(newHyperCube,cnt);
           //delete this newly created hypercube if it was already present
           if(!addedNew) {
			//newHyperCube->deleteVarConstants();
            delete newHyperCube;
		   }
         }
          
		  //add a copy of the hypercube and increment the count
          void addNewHyperCubeAndIncrementCount(HyperCube * const & hyperCube, double cnt) {
           bool addedNew = false;
           HyperCube *newHyperCube = new HyperCube(hyperCube);
           //add this hyperCube to the superclause
           addedNew = addHyperCube(newHyperCube);   
           incrementHyperCubeCount(newHyperCube,cnt);
           //delete this newly created hypercube if it was already present
           if(!addedNew) {
			//newHyperCube->deleteVarConstants();	
            delete newHyperCube;
		   }
         }


          //get the HyperCube corresponding to this predicate in the given hypercube id
          HyperCube * getPredicateHyperCube(int hindex, Predicate *pred) {
			   int predTermCnt = pred->getNumTerms();
			   //IntHashArray *varConstants;
			   IntArray *varConstants;
			   HyperCube * phyperCube = new HyperCube(predTermCnt);
               HyperCube * hyperCube = (*hyperCubes_)[hindex];

               if(useCT_) {
			    hyperCube->updateToImplicitTiedRepresentation();
			   }

			   //add the element at index 0
			   //phyperCube->addNewVariable(hyperCube->getVarConstants(0));
               for(int termno=0;termno<predTermCnt;termno++) {
                const Term *term = pred->getTerm(termno);
                int varId = -term->getId();
                assert(varId >= 0);
				
		        //copies are created in the end
				//varConstants = new IntArray(*(hyperCube->getVarConstants(varId)));
				varConstants = hyperCube->getVarConstants(varId);
                phyperCube->setVarConstants(varConstants,termno+1);
			    //phyperCube->addNewVariable(hyperCube->getVarConstants(varId));
               }
		
			   if(useCT_) {
			    hyperCube->updateToExplicitTiedRepresentation();
		        phyperCube->updateToExplicitTiedRepresentation();
			   }

			   //now create the copies for each of the varConstants
		       for(int termno=0; termno<predTermCnt;termno++) {
		        int termIndex = termno+1;
		        int refTermIndex = phyperCube->getReferenceVarId(termIndex);
		        if(!useCT_) {
					assert(refTermIndex == termIndex);
				}
				//this is a reference variable - nothing to do
		        if(refTermIndex != termIndex)
				 continue;
		        varConstants = new IntArray(*(phyperCube->getVarConstants(termIndex)));
		        phyperCube->setVarConstants(varConstants,termIndex);
		      }
			  return phyperCube;
          }
            
          //check whether this hypercube has non zero number of tuples 
		  int hasZeroNumTuples(HyperCube * const & hyperCube) {
			bool needDetailedCheck = false;
			int varCnt = hyperCube->getVarCount();
			for(int varId=1;varId<=varCnt;varId++) {
				  int refVarId = hyperCube->getReferenceVarId(varId);
				  if(varId != refVarId)
					   continue;
                  if(hyperCube->getVarConstantCount(varId) < varCnt) {
					   needDetailedCheck = true;
					   break;
				}
			}

		   if(!needDetailedCheck)
				return false;

		   Array<bool> *fixedIds = new Array<bool>(varCnt+1,false);   
		   int num = getNumTuples(hyperCube,fixedIds);
		   delete fixedIds;
		   return (num <= 0);
		  }

		  //get the number of partial tuples when the variables in the fixedIds
		  //array are fixed
		  int getNumTuples(HyperCube * const & hyperCube, Array<bool> * const & fixedIds) {
			   int num;
			   int varCnt = hyperCube->getVarCount();
			   
			   if(useCT_) 
			   {
				 
				 /*	
			     cout<<"hyperCube is "<<endl;
				 hyperCube->print(cout);
				 cout<<endl;
                 */
				 bool allAreFixed = true;
			     for(int varId=1;varId<=varCnt;varId++) {
					if((*fixedIds)[varId])
						 continue;
					allAreFixed = false;
			    }
			    if(allAreFixed) {
					return 1;
			    }
			    
				//note: it is important to initialize everything to one in the
				//following array of constantCnts
				Array<int> *constantCnts = new IntArray(varCnt+1,1);
				
				for(int varId=1;varId<=varCnt;varId++) {
				  int refVarId = hyperCube->getReferenceVarId(varId);
				  if(varId != refVarId)
					   continue;
				  if((*fixedIds)[varId]) 
					 continue;
                  (*constantCnts)[varId] = hyperCube->getVarConstantCount(varId);
				}

			    //following code takes care of modifying the counts appropriately
				//to incoprorate the constrained representation -
				//NOTE: THE FOLLOWING PIECE OF CODE ASSUMES THAT MAXIMUM SIZE OF AN
				//EQUIVALENCE CLASS OF CONSTRAINTS IS 3. IF THIS CONDITION IS VIOLATED,
				//THE COUNTS MAY NOT BE CORRECT
				
				Array<int> *seenCnts = new IntArray(varCnt+1,0);
				IntArrayHashArray *relevantConstraints = hyperCube->getRelevantConstraints(neqConstraints_);
                IntArray *constraint;
				int affectedVarId, oldCnt;
				for(int i=0;i<relevantConstraints->size();i++) {
					 constraint = (*relevantConstraints)[i];
		             int varId1 = -(*constraint)[0];
					 int varId2 = -(*constraint)[1];
					 affectedVarId = 0;

					 if(hyperCube->hasDisjointVarConstants(varId1,varId2))
						  continue;
                     
					 (*seenCnts)[varId1]++;
					 (*seenCnts)[varId2]++;

					 if((*fixedIds)[varId1] || ((*fixedIds)[varId2])) 
					  if((*fixedIds)[varId1] && ((*fixedIds)[varId2]))
                       continue;
                     
					 //following are various conditions for finding which variable
					 //should be affected for counts
					 //first - fixed ids should not be touched
					 if(!affectedVarId && (*fixedIds)[varId1]) 
					    affectedVarId = varId2;	  
                     if(!affectedVarId && (*fixedIds)[varId2]) 
						  affectedVarId = varId1;
					 
					 //second - change the count of the id which has been seen less number of times
				     if(!affectedVarId && (*seenCnts)[varId1] < (*seenCnts)[varId2])
						affectedVarId = varId1;
				     if(!affectedVarId && (*seenCnts)[varId2] < (*seenCnts)[varId1])
						affectedVarId = varId2;
					 
					 //third - change the count which is lower
					 if(!affectedVarId && (*constantCnts)[varId1] <= (*constantCnts)[varId2])
					   affectedVarId = varId1;
					 if(!affectedVarId && (*constantCnts)[varId1] >= (*constantCnts)[varId2])
					   affectedVarId = varId2;
					 
					 oldCnt = (*constantCnts)[affectedVarId];
				     (*constantCnts)[affectedVarId] = max(0,oldCnt-1);
			   }
				 
			   num = 1;
			   for(int varId=1;varId<=varCnt;varId++) {
			     if((*fixedIds)[varId])
					continue;
				 num = num * (*constantCnts)[varId];
			  }
			  
			  relevantConstraints->deleteItemsAndClear();
			  delete relevantConstraints;
			  delete seenCnts;
			  delete constantCnts;
			   
			 } else {
				 num = 1;
			     for(int varId=1;varId<=varCnt;varId++) {
			     if((*fixedIds)[varId])
					continue;
				 num = num * hyperCube->getVarConstantCount(varId);
			   }
		    }
		     return num;
		  }


          //get the count of partial tuples joining with a predicate in
          //the given hyperCube
          int getNumTuplesJoiningWithHyperCube(int hindex, Predicate *pred) {
               HyperCube *hyperCube = (*hyperCubes_)[hindex];
			   int varCnt = hyperCube->getVarCount();
               Array<bool> *fixedIds = new Array<bool>(varCnt+1,false);

			   //those appearing in this predicate should not be counted
               for(int i=0;i<pred->getNumTerms();i++) {
                const Term *term = pred->getTerm(i);
                int varId = -term->getId();
                assert(varId > 0);
                int refVarId = hyperCube->getReferenceVarId(varId);
				(*fixedIds)[refVarId] = true;
               }

               int num = getNumTuples(hyperCube,fixedIds);
               delete fixedIds;
               return num;
          }
          
          //get the total number of tuples which are represented by the
          //hyperCube at the given index
          int getNumTuples(int hindex) {
               HyperCube *hyperCube = (*hyperCubes_)[hindex];
			   int varCnt = hyperCube->getVarCount();
               Array<bool> *fixedIds = new Array<bool>(varCnt+1,false);
			   
			   int num = getNumTuples(hyperCube,fixedIds);
			  
			   /*
			   cout<<endl;
			   hyperCube->print(cout);
			   cout<<" ** "<<num<<endl;
			   //cout<<"  :  "<<(*hyperCubeCnts_)[hindex]<<" * "<<num<<endl;
			   cout<<endl;
			   */
			   delete fixedIds;
               return num;
		  }

          int getNumTuples(){ 
               int num = 0;
			   //cout<<"***************************************************"<<endl;
			   //cout<<"The Hypercubes and their counts are ="<<endl;
               for(int hindex=0;hindex<hyperCubes_->size();hindex++) {
                 num = num + getNumTuples(hindex);
               }
			   //cout<<"***************************************************"<<endl;
               return num;

          }
          
          double getNumTuplesWithDuplicates(){ 
               double num = 0;
			   //cout<<"***************************************************"<<endl;
			   //cout<<"The Hypercubes and their counts are ="<<endl;
               for(int hindex=0;hindex<hyperCubes_->size();hindex++) {
			     double duplicateCnt = (double)(*hyperCubeCnts_)[hindex];
                 num = num + duplicateCnt*getNumTuples(hindex);
               }
			   //cout<<"***************************************************"<<endl;
               return num;

          }
		  int getNumConstants(int hindex) {
			   return (*hyperCubes_)[hindex]->getNumConstants();
		  }

		  int getNumConstants(){ 
               int num = 0;
               for(int hindex=0;hindex<hyperCubes_->size();hindex++) {
                 num = num + getNumConstants(hindex);
               }
               return num;
          }

  double getOutputWt()
  {
    return outputWt_;
  }

  void addOutputWt(const double& outputWt)
  {
    outputWt_ += outputWt;
  }

		  
		  //merge the given superClause with this one
		  //Assumes that they correspond to the same underlying clause
		  void merge(SuperClause * const & superClause) {
			   HyperCube *hyperCube;	
			   for(int hindex=0;hindex<superClause->getNumHyperCubes();hindex++) {
                      hyperCube = superClause->getHyperCube(hindex);
                      double hcnt = superClause->getHyperCubeCount(hindex);
                      bool addedNew = addHyperCube(hyperCube);
				      incrementHyperCubeCount(hyperCube,hcnt);
					  if(!addedNew) {
						//hyperCube->deleteVarConstants();  
						delete hyperCube;
					  }
			   }
		  }

/*************************** Methods for refining the hypercubes *************************/

		/* get the predicate hypercube reverse index for the predicate
		 * refinements corresponding to this hypercube */
		Array<HyperCubeReverseIndex *> *getRefinedReverseIndex(int hindex, 
				                                               Array<HyperCubeRefinement *> *hcRefinementArr,
															   bool print) 
		{
			 int predCnt = clause_->getNumPredicates();
			 Array<HyperCube *> *refinedPHyperCubes;
			 Array<HyperCubeReverseIndex *> * reverseIndexArr = new Array<HyperCubeReverseIndex *>(predCnt,NULL);
			 HyperCubeRefinement *hcRefinement;
			 HyperCubeReverseIndex *reverseIndex;
			 HyperCube *phyperCube, *refinedPHyperCube;
             
			 Predicate *pred;
             int predId;
			 if(print) {
				  cout<<"Getting reverse index.."<<endl;
			 }

			 for(int pindex=0;pindex<predCnt;pindex++) {
				  pred = clause_->getPredicate(pindex);
				  predId = pred->getId();
				  reverseIndex = new HyperCubeReverseIndex();
				  reverseIndex->createIndex();
				  phyperCube = getPredicateHyperCube(hindex, pred);
				  hcRefinement = (*hcRefinementArr)[predId];
				  
				  /*cout<<"Pred id = "<<predId<<endl;
				  cout<<"Number of subset hypercubes in this refinement = "<<hcRefinement->getNumSubsetHyperCubes()<<endl;
                   */
				  //cout<<"While getting refined reverse index"<<endl;
				  if(print) {
				   cout<<"Pred HyperCube is:"<<endl;
				   phyperCube->print(cout);
				   cout<<endl;
				  }

				  refinedPHyperCubes = hcRefinement->getRefinedHyperCubes(phyperCube);
				  if(print) {
				   cout<<"Number of refined Predicate hypercubes obtained = "<<refinedPHyperCubes->size()<<endl;
				   cout<<"And these are "<<endl;
				  }
				  
				  for(int i=0;i<refinedPHyperCubes->size();i++) {
					   if(print) {
					     cout<<i<<":"<<endl;
					     (*refinedPHyperCubes)[i]->print(cout);
					   }
					    //refinedPHyperCube = new HyperCube((*refinedPHyperCubes)[i]); 
					    refinedPHyperCube = (*refinedPHyperCubes)[i]; 
					    reverseIndex->addHyperCube(refinedPHyperCube);
				 }
				  
				  (*reverseIndexArr)[pindex] = reverseIndex;
				  delete phyperCube;
			 }
			 return reverseIndexArr;
		}
        
		//refine the hypercubes in this superclause according to hypercuberefinement
		//Array. Also update the new hypercuberefinement array
        void refineHyperCubes(Array<HyperCubeRefinement *> *hcRefinementArr,
                              Array<HyperCubeRefinement *> *newHCRefinementArr,
							  Domain * const & domain,
							  bool print) {
		  //print = true;
		  if(print) {
		   cout<<"------------------------------------------------------------"<<endl;
		   cout<<"Ok, came in to get the Refined HyperCubes..."<<endl;
		  }

		  Array<HyperCube *> *newHyperCubes = new Array<HyperCube *>();	 
		  Array<double> *newHyperCubeCnts = new Array<double>();	 
		  
		  Array<HyperCube *> *refinedHyperCubes;
		  Array<HyperCubeReverseIndex *> * reverseIndexArr;
		  IntHashArray *unknownPredIndices = new IntHashArray();
		  for(int predno=0;predno<clause_->getNumPredicates();predno++) {
			   unknownPredIndices->append(predno);
		  }

		  for(int hindex=0;hindex<getNumHyperCubes();hindex++) {

		   //first get the phypercubes and the reverseindex according to
		   //current refinement
		    if(print) {
			 cout<<endl<<"----------------------------------------------"<<endl;
			 cout<<"hindex = "<<hindex<<endl;
		     HyperCube *hyperCube = getHyperCube(hindex);
		    
			 cout<<"Getting the refined hypercubes for "<<endl;
		     hyperCube->print(cout);
			 cout<<endl;
			}
			reverseIndexArr = getRefinedReverseIndex(hindex,hcRefinementArr,print);
		   
		   //now perform the hypercube join (based on reverseindex) and get the new hypercube 
		   //refinement
		   
		   if(print) {
		    cout<<endl;
            cout<<"Now performing the join.."<<endl;
		   }

		    SuperClause *superClause = NULL;
		    refinedHyperCubes = getClauseHyperCubes(clause_,superClause,reverseIndexArr,
					                               newHCRefinementArr,unknownPredIndices,
												   neqConstraints_,useCT_,domain);
		   
		   newHyperCubes->append(refinedHyperCubes);
		   
		   if(print) {
			cout<<"size of refined hypercube = "<<refinedHyperCubes->size()<<endl;
		   }

		   //update the cnts array as well

		   for(int i=0;i<refinedHyperCubes->size();i++) {
		    newHyperCubeCnts->append((*hyperCubeCnts_)[hindex]);
		   }
           
		   delete refinedHyperCubes;

		   //(*hyperCubes_)[hindex]->deleteVarConstants();
		   delete (*hyperCubes_)[hindex];

		   //delete the individual elements
		   reverseIndexArr->deleteItemsAndClear();
           delete reverseIndexArr;
		  }
		  
		  hyperCubes_->clear();
		  hyperCubeCnts_->clear();

		  for(int i=0;i<newHyperCubes->size();i++) { 
			   HyperCube *hc = (*newHyperCubes)[i];
			   double cnt = (*newHyperCubeCnts)[i];
			   int index = hyperCubes_->find(hc);

			   //add this this hypercube if not added earlier
			   if(index < 0) {
                 hyperCubes_->append((*newHyperCubes)[i]);
		         hyperCubeCnts_->append(cnt);
			   } else { //simply update the count
					delete hc;
					(*hyperCubeCnts_)[index] += cnt;
			   }
		  }

		  if(print) {
			   cout<<"old counts were: "<<endl;
			   for(int i=0;i<hyperCubeCnts_->size();i++)
			    cout<<(*hyperCubeCnts_)[i]<<" ";
			   cout<<endl;

			   cout<<"new counts are: "<<endl;
			   for(int i=0;i<newHyperCubeCnts->size();i++)
			    cout<<(*newHyperCubeCnts)[i]<<" ";
			   cout<<endl;
		  }

		  
		  //(*hyperCubeCnts_) = (*newHyperCubeCnts);
		  
		  delete unknownPredIndices;
		  delete newHyperCubes;
		  delete newHyperCubeCnts;
		 
		  if(print) {
		   cout<<"              Done getting the refined hypercubes...       "<<endl;
		   cout<<"------------------------------------------------------------"<<endl;
		  }
		}

/*************************** End Methods for Refining the hypercubes *************************/
        
		
		//print the hyperCubes
        ostream& print(ostream& out) {
           for(int i=0;i<hyperCubes_->size();i++) {
			 (*hyperCubes_)[i]->print(out);
		   }
		   return out;
		}
          
		
		void printClauses(Domain * const & domain, ostream & out) {
			  Array<IntArray *> * tuples;
			  IntArray *tuple;
			  HyperCube *hyperCube;
			  Clause *gndClause;
			  Predicate *pred,*gndPred;
			  IntArray *ptuple;
			  int constant;
			  for(int hindex=0;hindex<hyperCubes_->size();hindex++) {
				  hyperCube = (*hyperCubes_)[hindex];
			      tuples = hyperCube->getTuples();
				  for(int tindex=0;tindex<tuples->size();tindex++) {
                     tuple = (*tuples)[tindex];
				     gndClause = new Clause(); 

					 for(int predno=0;predno<clause_->getNumPredicates();predno++) {
					   ptuple = new IntArray();
                       pred = clause_->getPredicate(predno);
                       for(int termno=0;termno<pred->getNumTerms();termno++) {
                          int varId = -(pred->getTerm(termno)->getId());
						  constant = (*tuple)[varId];
						  ptuple->append(constant);
					   }
					   gndPred = domain->getPredicate(ptuple,pred->getId());
					   //gndPred = getGroundPredicate(pred->getTemplate(),ptuple);
					   gndPred->setSense(pred->getSense());
					   gndClause->appendPredicate(gndPred);
					   delete ptuple;
					 }
					 gndClause->print(out,domain);
				     out<<endl;
					 delete gndClause;
					 delete tuple;
				  }
			      tuples->clear();
				  delete tuples;
			  }
		}


        //static function
         void static resetIndex() { superClauseIndex__ = 0;}

     private:
          Clause* clause_;
          double wtScale_;
          HyperCubeHashArray *hyperCubes_;
          Array<double> * hyperCubeCnts_;
          IntArrayHashArray * neqConstraints_;
          Array<int> * varIdToCanonicalVarId_;
          int superClauseId_;
          int parentSuperClauseId_;
  double outputWt_;
          
		  //whether to use the constrained representation
          bool useCT_;
		  
          
          static int superClauseIndex__;

};

//extern function defined in superpred.cpp
extern void createSuperClauses(Array<Array<SuperClause*>*> * const & superClausesArr,
                               Domain * const & domain);

typedef hash_map<Array<int>*,SuperClause *, HashIntArray, EqualIntArray> IntArrayToSuperClause;

#endif
