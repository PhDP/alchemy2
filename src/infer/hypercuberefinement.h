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
#ifndef _HYPERCUBE_REFINEMENT_H_MAY_2008
#define _HYPERCUBE_REFINEMENT_H_MAY_2008

#include <ext/hash_set>
#include <ext/hash_map>
#include <ostream>

#include "util.h"
#include "array.h"
#include "hashint.h"
#include "mln.h"
#include "domain.h"
#include "hypercube.h"

#define USE_INLINE 0


using namespace std;
using namespace __gnu_cxx;

typedef hash_map<int,Array<int> *, HashInt, EqualInt> IntToIntArrayMap;

class HyperCubeRefinement {
 
public:
 HyperCubeRefinement() {
		subsetHCToBaseHC_ = new HyperCubeToHyperCubeMap();
	    baseHCToSubsetHCs_ = new HyperCubeToHyperCubeHashArrayMap();
		subsetHCToRefinedHCs_ = new HyperCubeToHyperCubeArrayMap();
   }


 ~HyperCubeRefinement() {
	  
	   HyperCube *baseHyperCube;

	   HyperCubeToHyperCubeHashArrayMap::iterator hitr;
	   HyperCubeToHyperCubeArrayMap::iterator itr;

	   HyperCubeHashArray *subsetHyperCubes;
       Array<HyperCube *> *baseHyperCubes;
	  
	   baseHyperCubes = new Array<HyperCube *>(); 

	   //delete the refined hypercubes
	  
	   /*
       Array<HyperCube *> *refinedHyperCubes;
	   for(itr = subsetHCToRefinedHCs_->begin();
		   itr != subsetHCToRefinedHCs_->end();
		   itr++) {
			refinedHyperCubes = itr->second;
            
			//this is deleted somewhere else 
			//refinedHyperCubes->deleteItemsAndClear();
			delete refinedHyperCubes;
	   }*/

	   //deletes each individual refined hypercube as well - 
	   //a cleaner implementation!
	   deleteRefinedHyperCubes();
        
	   //delete subset hypercubes 
		for(hitr = baseHCToSubsetHCs_->begin();
		   hitr != baseHCToSubsetHCs_->end();
		   hitr++) {
			
			baseHyperCube = hitr->first; 
			baseHyperCubes->append(baseHyperCube);
			
			subsetHyperCubes = hitr->second;
            subsetHyperCubes->deleteItemsAndClear();
			delete subsetHyperCubes;
		}
	   
		//delete the base hypercubes
		baseHyperCubes->deleteItemsAndClear();
		delete baseHyperCubes;
  
	   // now, delete the mapping objects
	   delete baseHCToSubsetHCs_;
	   delete subsetHCToRefinedHCs_;	
	   delete subsetHCToBaseHC_;
 }
 
 void deleteRefinedHyperCubes() {
	   HyperCubeToHyperCubeArrayMap::iterator itr;
       Array<HyperCube *> *refinedHyperCubes;
	   HyperCube *refinedHyperCube, *subsetHyperCube;
       HyperCubeHashArray *uniqueRefinedHyperCubes = new HyperCubeHashArray();


	   for(itr = subsetHCToRefinedHCs_->begin();
		   itr != subsetHCToRefinedHCs_->end();
		   itr++) {
			subsetHyperCube = itr->first;
            
			/*
			cout<<endl<<"*********************************"<<endl;
			cout<<"Processing the Subset HyperCube "<<endl;
			subsetHyperCube->print(cout);
			cout<<endl<<"**********************************"<<endl;
            */

			refinedHyperCubes = itr->second;

			for(int i=0;i<refinedHyperCubes->size();i++) {
				 refinedHyperCube = (*refinedHyperCubes)[i];
			     uniqueRefinedHyperCubes->append(refinedHyperCube);
				 //cout<<"hyper cube is: "<<endl;;
				 //refinedHyperCube->print(cout);
				 //delete refinedHyperCube;
			}
			//cout<<endl<<"**********************************"<<endl;
			//refinedHyperCubes->deleteItemsAndClear();
			
			delete refinedHyperCubes;
	   }
	   
	   //cout<<"size being deleted = "<<uniqueRefinedHyperCubes->size()<<endl;
	   uniqueRefinedHyperCubes->deleteItemsAndClear();
	   delete uniqueRefinedHyperCubes;

	   subsetHCToRefinedHCs_->clear();
 }


 int getNumSubsetHyperCubes() {
     return baseHCToSubsetHCs_->size();
 }

 //get the base hypercube corresponding to this subset hypercube
 HyperCube * getBaseHyperCube(HyperCube * const & subsetHyperCube) {
	  HyperCubeToHyperCubeMap::iterator itr;
	  itr = subsetHCToBaseHC_->find(subsetHyperCube);
	  if(itr != subsetHCToBaseHC_->end())
		   return itr->second;
	  else
		   return NULL;
 }
         
 //return the hypercubes corresponding to the given subset hypercube
 Array<HyperCube *> *getRefinedHyperCubes(HyperCube * const & subsetHyperCube) {
	  int i=0;
	  HyperCubeToHyperCubeArrayMap::iterator itr,itr1;
	  //cout<<"number of subset hypercubes = "<<subsetHCToRefinedHCs_->size()<<endl;
	  itr = subsetHCToRefinedHCs_->find(subsetHyperCube);
	  if(itr == subsetHCToRefinedHCs_->end()) {
		   cout<<"Wow, returning Null!!"<<endl;
		   subsetHyperCube->print(cout);
		   cout<<endl;
		  
		   cout<<"Hypercubes are : "<<endl;
		   for(itr1 = subsetHCToRefinedHCs_->begin();
               itr1 != subsetHCToRefinedHCs_->end();
			   itr1++) {
               cout<<i++<<":"<<endl;
			   (itr1->first)->print(cout);
			   cout<<endl;
		   }
		   cout<<"----------------------------------------------------"<<endl;
		   return NULL;
	  }
	  return itr->second;
 }



 //add subset hypercube to the appropriate bucket (corresponding to the
 //given base hypercube) - return true if subsethypercube was succesfully
 //added. return false if it was already present
 void addSubsetHyperCube(HyperCube * const & inpBaseHyperCube,
		                 HyperCube * const & subsetHyperCube) {
	  HyperCube *baseHyperCube;
	  HyperCubeHashArray * subsetHyperCubes;
	  HyperCubeToHyperCubeHashArrayMap::iterator itr;
	  
	  /*
	  cout<<"came to add to mapping:"<<endl;
	  cout<<"Base HyperCube:"<<endl;
	  inpBaseHyperCube->print(cout);
	  cout<<endl;

	  cout<<"Subset HyperCube:"<<endl;
	  subsetHyperCube->print(cout);
	  cout<<endl;
      */
	  
	  //add to the mapping from subsethypercube to basehypercube
	  itr = baseHCToSubsetHCs_->find(inpBaseHyperCube);
	  if(itr == baseHCToSubsetHCs_->end()) {
		   //add to the mapping from basehypercube to subsethypercubes
           //baseHyperCube = new HyperCube(inpBaseHyperCube);
           baseHyperCube = inpBaseHyperCube;
		   subsetHyperCubes = new HyperCubeHashArray();
		   (*baseHCToSubsetHCs_)[baseHyperCube] = subsetHyperCubes;
	  } else {
		   //inpBaseHyperCube->deleteVarConstants();
		   delete inpBaseHyperCube;
		   baseHyperCube = itr->first;
		   subsetHyperCubes = itr->second;
	  }

	  int index = subsetHyperCubes->append(subsetHyperCube);
	  
	  //check if this subset hypercube has already been added to the refinement
	  if(index < 0) {
		   delete subsetHyperCube;
		   //return false;
	  } else {
	     (*subsetHCToBaseHC_)[subsetHyperCube] = baseHyperCube;
		 //return true;
	  }
 }


 //create the set of refined hypercubes and the mapping from subset hypercubes
 //to refined hypercubes
 //Note: this should be done after all the subset hypercubes have been added
 void createRefinedHyperCubeMapping() {
	  //int createdSize = 0;
	  
	  HyperCube *subsetHyperCube, *baseHyperCube; //, *refinedHyperCube;
	  HyperCubeHashArray * subsetHyperCubes;
	  Array<HyperCube *> * refinedHyperCubes, *intersectingHyperCubes; 
      
	  HyperCubeToHyperCubeHashArrayMap::iterator itr;
	  int n=0;
	  //cout<<"Total number of hypercubes to be processed for this Predicate = "<<baseHCToSubsetHCs_->size()<<endl;

	  //while(true) {

	  n = 0;
	  for(itr = baseHCToSubsetHCs_->begin();
		   itr != baseHCToSubsetHCs_->end();
		   itr++) {
		    //cout<<"Processing Base HyperCube "<<n++<<endl;
		    baseHyperCube = itr->first;
			subsetHyperCubes = itr->second;
			
			/*
			cout<<"for the base hypercube "<<endl;
			baseHyperCube->print(cout);
			cout<<endl;
		    cout<<"Processing the following subset hypercubes:"<<endl;
			
			for(int i=0;i<subsetHyperCubes->size();i++) {
				 cout<<i<<":"<<endl;
				 subsetHyperCube = (*subsetHyperCubes)[i];
				 subsetHyperCube->print(cout);
				 cout<<endl;
			}
			 */
            
		    refinedHyperCubes = createRefinedHyperCubes(subsetHyperCubes);
		    
			/*
			refinedHyperCubes->deleteItemsAndClear();
			delete refinedHyperCubes;
			cout<<"going in for next round.."<<endl;
			continue;
            */
	
			/*
			HyperCube *refinedHyperCube;
			cout<<"After refining, the following are refined hypercubes:"<<endl;
			for(int i=0;i<refinedHyperCubes->size();i++) {
				 cout<<i<<":"<<endl;
				 refinedHyperCube = (*refinedHyperCubes)[i];
				 refinedHyperCube->print(cout);
				 cout<<endl;
			}
		   */

			/*
			cout<<"After refining, the following are subset hypercubes:"<<endl;
			for(int i=0;i<subsetHyperCubes->size();i++) {
				 cout<<i<<":"<<endl;
				 subsetHyperCube = (*subsetHyperCubes)[i];
				 subsetHyperCube->print(cout);
				 cout<<endl;
			}
            */
			/*
            cout<<endl;
			cout<<"------------------------------------------------------------ "<<endl;
			cout<<"Now Getting the refined Hypercubes for each subset hypercube "<<endl;
			*/


			//HyperCubeHashArray *allIntersectingHyperCubes = new HyperCubeHashArray();
            //HyperCube *hc;

			//now iterate over each subset hypercube and find
			//the intersecting refined hypercubes
			for(int i=0;i<subsetHyperCubes->size();i++) {
				 subsetHyperCube = (*subsetHyperCubes)[i];
			     
				/*
				 cout<<"-------------------------------------------------------------"<<endl;
			     cout<<"Subset hypercube is = "<<endl;
	             subsetHyperCube->print(cout);
	             cout<<endl;
                */
				 intersectingHyperCubes = getIntersectingHyperCubes(subsetHyperCube,refinedHyperCubes);
				 (*subsetHCToRefinedHCs_)[subsetHyperCube] = intersectingHyperCubes;
				 

				 //test
			     /*
				 for(int j=0;j<intersectingHyperCubes->size();j++) {
				  hc = (*intersectingHyperCubes)[j];
			      allIntersectingHyperCubes->append(hc);
				 }
				 */
				 
				 //cout<<"-------------------------------------------------------------"<<endl;
			}
			//cout<<"------------------------------------------------------------ "<<endl;
		
			/*
            if(allIntersectingHyperCubes->size() != refinedHyperCubes->size()) {
				 cout<<"THERE IS A PROBLEM HERE!!!!!"<<endl;
				 cout<<"refined size = "<<refinedHyperCubes->size()<<endl;
				 cout<<"intersecting size = "<<allIntersectingHyperCubes->size()<<endl;
			}
            
			createdSize += refinedHyperCubes->size();
			delete allIntersectingHyperCubes; 
			*/

			//delete subsetHyperCubes;
			
			delete refinedHyperCubes;
            
			//cout<<endl;
			/*
			cout<<"************************************************************"<<endl;
	        cout<<"         Done procesing this base hypercube                      "<<endl;
            cout<<"************************************************************"<<endl;
	        */
	   }
	  //}
	  //cout<<"created size = "<<createdSize<<endl;
 }


 /************************ Private Methods ************************************************/
private:

 //get the array of refined hypercubes for the given set of hypercubes
 Array<HyperCube *> * createRefinedHyperCubes(HyperCubeHashArray * const & subsetHyperCubes) {
	  
	  //cout<<"=================================================================="<<endl;
	  //cout<<"Creating Refined HyperCubes....."<<endl;
	  //cout<<"=================================================================="<<endl;
	  
	  Array<HyperCube *> *refinedHyperCubes = new Array<HyperCube *>();
	  Array<HyperCube *> *remainingHyperCubes = new Array<HyperCube *>();
	  HyperCube *subsetHyperCube;
	  bool print = true;
	  for(int index=0;index<subsetHyperCubes->size();index++) {
		   subsetHyperCube = (*subsetHyperCubes)[index];
		   
		  
		   /*
		   cout<<endl<<endl<<index<<": subset hypercube is "<<endl;
		   subsetHyperCube->print(cout);
		   cout<<endl;
		   */

		   if(false) {
		   //if(index == 314) {
		   //}
				print = true;
				cout<<"yup! yup! print is true.."<<endl;
		   } else {
				print = false;
		   }
		   
		   HyperCube *hc = new HyperCube(subsetHyperCube);
				
		  
		   /*
		   cout<<"Created a new hypercube: "<<endl;
		   hc->print(cout);
		   cout<<endl;
		   */

		   updateRefinedHyperCubes(hc, refinedHyperCubes,remainingHyperCubes, print, 0,index);
		   refinedHyperCubes->append(remainingHyperCubes);
		   remainingHyperCubes->clear();
	  
	  
        /*
		cout<<"---------------------------------------------------------"<<endl;
		cout<<"Refined HyperCubes after interesting with this subset hypercube are:"<<endl;
			
	    for(int i=0;i<refinedHyperCubes->size();i++) {
				  cout<<i<<":"<<endl;
				  (*refinedHyperCubes)[i]->print(cout);
				  cout<<endl;
	    } 
		cout<<"---------------------------------------------------------"<<endl;
	   */
	  }

	  delete remainingHyperCubes;
	  return refinedHyperCubes;
 }
 
void updateRefinedHyperCubes(HyperCube * const & toIntersectHyperCube, 
		                      Array<HyperCube *> * const & refinedHyperCubes, 
		                      Array<HyperCube *> * const & remainingHyperCubes,
							  bool print, int depth, int index);


//get the array of refined hypercubes which have non null intersection with the
//given hypercube
Array<HyperCube *> *getIntersectingHyperCubes(HyperCube * const & subsetHyperCube,
		                                       Array<HyperCube *> * const & refinedHyperCubes) {
	  Array<HyperCube *> *intersectingHyperCubes = new Array<HyperCube *>();
	  HyperCube *refinedHyperCube;
	  
	  //cout<<"Refined hypercubes are = "<<endl;
	  for(int i=0;i<refinedHyperCubes->size();i++) {
		     refinedHyperCube = (*refinedHyperCubes)[i];
			 if(subsetHyperCube->hasIntersection(refinedHyperCube)) {
				  intersectingHyperCubes->append(refinedHyperCube);
                  /*
				  cout<<i<<":"<<endl;
	              refinedHyperCube->print(cout);
	              cout<<endl;
				  */
			 }
	  }
	  return intersectingHyperCubes; 
 }

/************************ End Private Methods ************************************************/

//variables
private:
 HyperCubeToHyperCubeMap *subsetHCToBaseHC_;
 
 HyperCubeToHyperCubeHashArrayMap *baseHCToSubsetHCs_;
 //Refined hypercubes are disjoint with each other
 HyperCubeToHyperCubeArrayMap *subsetHCToRefinedHCs_;

};

#endif

