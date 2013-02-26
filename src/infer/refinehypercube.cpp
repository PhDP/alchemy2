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
#include <ext/hash_set>
#include <ext/hash_map>

#include "util.h"
#include "hypercube.h"
#include "hypercuberefinement.h"

#define USE_RECURSIVE 0
static bool PRINTED = false;

//update the array of refined hypercubes by taking intersection with the input
//subset hypercube - two implementations below. One iterative and other recursive


#if !USE_RECURSIVE
void HyperCubeRefinement::updateRefinedHyperCubes(HyperCube * const & subsetHyperCube, 
		                      Array<HyperCube *> * const & refinedHyperCubes, 
		                      Array<HyperCube *> * const & remainingHyperCubes,
							  bool print, int depth, int hcIndex) {
	      
 	      if(!PRINTED) {
 	       cout<<"Using Iterative version.."<<endl;
 	       PRINTED=true; 
 	      }
		
		//if there is nothing to intersect, add the subsetHypercube
		//to the refined list and return
	    if(refinedHyperCubes->size() == 0) {
			  refinedHyperCubes->append(subsetHyperCube);
			  return;
		}
		  
	      HyperCube *refinedHyperCube, *toIntersectHyperCube;
	      Array<HyperCube *> * toIntersectHyperCubeComponents, *refinedHyperCubeComponents;
		  Array<HyperCube *> * toIntersectHyperCubes = new Array<HyperCube *>();
          Array<int> * startIndices = new Array<int>();
          
		  toIntersectHyperCubes->append(subsetHyperCube);
          //this array stores an index for each toIntersect hypercube
		  //- where to start intersecting in the refined hypercube array
		  startIndices->append(0);
          
		  int refinedIndex;

		  while(toIntersectHyperCubes->size() > 0) {

			assert(toIntersectHyperCubes->size() == startIndices->size());
			
			toIntersectHyperCube = toIntersectHyperCubes->removeLastItem();
            refinedIndex = startIndices->removeLastItem();     
			
			while(refinedIndex < refinedHyperCubes->size()) {
		      refinedHyperCube = (*refinedHyperCubes)[refinedIndex];

		     //if toInterserct hypercube is same as the refinedHyperCube
		     if(refinedHyperCube->same(toIntersectHyperCube)) {
				 delete toIntersectHyperCube;
		         break;
			 }

		     //if there is no intersection between the toIntersect hypercube and
		     //the first hypercube in the disjoint list
		     bool hasIntersection = refinedHyperCube->hasIntersection(toIntersectHyperCube, print);
		     if(!hasIntersection) {
				  refinedIndex++;
				  continue;
		     } 
			 
		     //all other cases i.e. there is an intersection and it is not trivial
		     toIntersectHyperCubeComponents = toIntersectHyperCube->getMinus(refinedHyperCube);
		     for(int j=0;j<toIntersectHyperCubeComponents->size();j++) {
              toIntersectHyperCubes->append((*toIntersectHyperCubeComponents)[j]);
			  startIndices->append(refinedIndex+1);
			 }
		     delete toIntersectHyperCubeComponents;
		  
		     refinedHyperCubeComponents = refinedHyperCube->getMinus(toIntersectHyperCube);
		     refinedHyperCubes->append(refinedHyperCubeComponents);
		     delete refinedHyperCubeComponents;
		  
		     refinedHyperCube->intersect(toIntersectHyperCube);
		     delete toIntersectHyperCube;
			 break;
			}
			
			//if it came out of the loop because there was nothing to intersect, then append the 
			//toIntersect hypercube at the end of the refined list
			if(refinedIndex >= refinedHyperCubes->size()) {
			 refinedHyperCubes->append(toIntersectHyperCube);
			}
		  }
		  delete toIntersectHyperCubes;
		  delete startIndices;
 } 

#endif

#if USE_RECURSIVE
 void HyperCubeRefinement::updateRefinedHyperCubes(HyperCube * const & toIntersectHyperCube, 
 		                      Array<HyperCube *> * const & refinedHyperCubes, 
 		                      Array<HyperCube *> * const & remainingHyperCubes,
 							  bool print, int depth, int hcIndex) {
 	      static int Max_Depth;
 	 
 		  if(false) {
 		  //if(hcIndex == 400 || hcIndex == 416) {
 		  //} 
 			   print = true;
 		  }
 
 		  //cout<<"Printed is "<<PRINTED<<endl;
 	      if(!PRINTED) {
 	       cout<<"Using Recursive version.."<<endl;
 	       PRINTED=true; 
 	      }
 
 		  if(Max_Depth < depth) {
 		   //cout<<"Max recursion depth is = "<<depth<<endl;
 		   Max_Depth = depth;
 		   }
 	      HyperCube * refinedHyperCube;
 	      Array<HyperCube *> * toIntersectHyperCubeComponents, *refinedHyperCubeComponents;
 
 		  if(print) {
 			   cout<<"yo man! I came here..."<<endl;
 		       cout<<"came 1"<<endl;
 		  }
 
 		  //if there is nothing to intersect, add the toIntersect Hypercube
 		  //to the refined list and return
 	      if(refinedHyperCubes->size() == 0) {
 			   remainingHyperCubes->append(toIntersectHyperCube);
 			   return;
 		  }
 		  if(print) {
 			   cout<<"yo man! I came here..."<<endl;
 		       cout<<"came 2"<<endl;
 		  }
 		  
 		  refinedHyperCube = refinedHyperCubes->removeItemFastDisorder(0);
 
 		  if(print) {
 		   cout<<endl;
 		   cout<<"Processing refined hypercube .."<<endl;
 		   refinedHyperCube->print(cout);
 		   cout<<endl;
 		  }
 
 		  //if toInterserct hypercube is same as the first hypercube in the disjoint list
 		  if(refinedHyperCube->same(toIntersectHyperCube)) {
 				 delete toIntersectHyperCube;
 			     remainingHyperCubes->append(refinedHyperCube);
 		         return;
 		  }
 		  if(print) {
 			   cout<<"yo man! I came here..."<<endl;
 		       cout<<"came 3"<<endl;
 		  }
 
 		  //if there is no intersection between the toIntersect hypercube and
 		  //the first hypercube in the disjoint list
 		  //if(!refinedHyperCube->hasIntersection(toIntersectHyperCube))
 		  
 		  bool hasIntersection = refinedHyperCube->hasIntersection(toIntersectHyperCube, print);
 		  if(!hasIntersection) { 
 		    if(print) {
 			   cout<<"yo man! I came here..."<<endl;
 		       cout<<"came 3.1"<<endl;
 		    } 
 			 updateRefinedHyperCubes(toIntersectHyperCube,refinedHyperCubes,remainingHyperCubes, 
 					                 print, depth+1,hcIndex);
 			 refinedHyperCubes->append(refinedHyperCube);
 		      if(print) {
 			   cout<<"yo man! I came here..."<<endl;
 		       cout<<"came 3.2"<<endl;
 		     }
 		     return;
 		  }
 
 		  
 		  if(print) {
 			   cout<<"yo man! I came here..."<<endl;
 		       cout<<"came 4"<<endl;
 		  }
 
 
 		  //all other cases i.e. there is an intersection and it is not trivial
 		  toIntersectHyperCubeComponents = toIntersectHyperCube->getMinus(refinedHyperCube);
 		  for(int j=0;j<toIntersectHyperCubeComponents->size();j++) {
             updateRefinedHyperCubes((*toIntersectHyperCubeComponents)[j],refinedHyperCubes,remainingHyperCubes, 
 					                print, depth+1, hcIndex);
 		  }
 		  delete toIntersectHyperCubeComponents;
 		  
 		  refinedHyperCubeComponents = refinedHyperCube->getMinus(toIntersectHyperCube);
 		  refinedHyperCubes->append(refinedHyperCubeComponents);
 		  delete refinedHyperCubeComponents;
 		  
 		  refinedHyperCube->intersect(toIntersectHyperCube);
 		  remainingHyperCubes->append(refinedHyperCube);
 		  
 		  //was causing memory leak earlier
 		  delete toIntersectHyperCube;
 		  if(print) {
 			   cout<<"yo man! I came here..."<<endl;
 		       cout<<"came 5"<<endl;
 		  }
  } 
 
#endif
