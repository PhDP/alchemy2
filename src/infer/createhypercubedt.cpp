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

#include "hashint.h"
#include "inferutil.h"
#include "hypercube.h"
#include "hypercuberefinement.h"

using namespace std;
using namespace __gnu_cxx;

//extern function
//merge the input hypercubes into larger hypercubes 
//(make as large as possible)
void mergeHyperCubes(Array<HyperCube *> * const & inpHyperCubes);


IntToIntHashArrayMap * getTupleIndex(IntArrayHashArray * const & tupleConstantsArr,
		                             int varId) {
	 IntToIntHashArrayMap * constantToTupleIndices = new IntToIntHashArrayMap();
	 IntArray *tupleConstants;
	 IntHashArray *tupleIndices;
	 for(int i=0;i<tupleConstantsArr->size();i++) {
		  tupleConstants = (*tupleConstantsArr)[i];
		  int constant = (*tupleConstants)[varId];
		  tupleIndices = (*constantToTupleIndices)[constant];
		  if(tupleIndices == NULL) {
			   tupleIndices = new IntHashArray();
		       (*constantToTupleIndices)[constant] = tupleIndices;
		  }
		  tupleIndices->append(i);
	 }
     return constantToTupleIndices;
}

HyperCube * getBoundingHyperCube(IntArrayHashArray * const & tupleConstantsArr) {
	 if(tupleConstantsArr->size() == 0)
		  return NULL;
	 
	 IntArray * tupleConstants;
	 
	 tupleConstants = (*tupleConstantsArr)[0];
     int varCnt = tupleConstants->size()-1;
	 HyperCube *boundingHyperCube = new HyperCube(varCnt);
	 
	 Array<int> *varConstants;
	 Array<int> *tmpArr = new Array<int>();

	 for(int varId=1;varId<=varCnt;varId++) {
	      varConstants = new Array<int>();
		  for(int i=0;i<tupleConstantsArr->size();i++) {
			  tupleConstants = (*tupleConstantsArr)[i];
			  int constant = (*tupleConstants)[varId];
              tmpArr->clear();
			  tmpArr->append(constant);
			  unionArrayKeepOrder(varConstants,tmpArr);
		  }
		  boundingHyperCube->setVarConstants(varConstants,varId);
	 }
	 delete tmpArr;
	 return boundingHyperCube;
}


void populateMoveOdds(HyperCube * const & hyperCube, 
				      IntToIntHashArrayMap * const & constantToTupleIndices,
	                  Array<double> * const & moveOdds, int splitVarId) {
    
	 Array<int> *varConstants = hyperCube->getVarConstants(splitVarId);
	 IntHashArray *tupleIndices;
     
	 //number of data points that will move by moving one constant from
	 //this (primary) hypercube to the other hypercube
	 int dataMoveCnt = (hyperCube->getNumTuples())/varConstants->size();
     int tupleMoveCnt;

	 for(int i=0;i<varConstants->size();i++) {
      int constant = (*varConstants)[i];
	  tupleIndices = (*constantToTupleIndices)[constant];
      if(tupleIndices == NULL) {
		 tupleMoveCnt = 0;
	  } else {
	     tupleMoveCnt = tupleIndices->size();
	  }
	  moveOdds->append(((double)tupleMoveCnt)/dataMoveCnt);
    }
}


void createHyperCubes(HyperCube * const & hyperCube, 
		              IntArrayHashArray * const & tupleConstantsArr,
					  Array<HyperCube *> * const & hyperCubes,
					  double noise) {
     
	 /*
	 cout<<"******************************************"<<endl;
	 cout<<"Bounding Hypercube is.."<<endl;
	 hyperCube->print(cout);
	 cout<<endl;
	 cout<<endl;
	 */



	 int varCnt = hyperCube->getVarCount();
	 int tupleCnt = tupleConstantsArr->size();
	 int hcNumTuples = hyperCube->getNumTuples();
	 
	 double noiseRatio = noise/hcNumTuples;
	 
	 //cout<<"noise is = "<<noise<<endl;
	 //cout<<"hc numtuples is ="<<hcNumTuples<<endl;
	 double classOdd = (0.0+tupleCnt)/hcNumTuples;
	
	 /*
	 if(classOdd-noiseRatio  <= 0)
		  return;
	 
	 if(classOdd+noiseRatio >= 1) {
      
	  hyperCubes->append(new HyperCube(hyperCube));
	  //cout<<"ok, at least added something..."<<endl;
	  //hyperCube->print(cout);
	  //cout<<endl;
	  return;
	 }*/

	 
	 //cout<<"tupleCnt = "<<tupleCnt<<" hcNumTuples = "<<hcNumTuples<<endl;
	 if((classOdd-noiseRatio  <= 0) ||(classOdd+noiseRatio >= 1)) {
      
	  if(tupleCnt > hcNumTuples/2) {
	   hyperCubes->append(new HyperCube(hyperCube));
	  }
	  
	 
	  //cout<<"ok, at least added something..."<<endl;
	  //hyperCube->print(cout);
	  //cout<<endl;
	 
	  return;
	 }


	 //cout<<"came here to process.."<<endl;
	 
	 IntToIntHashArrayMap *constantToTupleIndices;
	 Array<double> *bestMoveOdds = NULL;
	 Array<double> *moveOdds;
	 
	 double totalOdd = 0;
	 double bestTotalOdd = 0;
	 int bestVarId = -1;

	 for(int varId = 1; varId <= varCnt; varId++) {

	     moveOdds = new Array<double>();
		 //only consider variables for which there are at least two constants
		 if(hyperCube->getVarConstantCount(varId) <= 1)
			  continue;
		 
		 totalOdd = 0;
		 constantToTupleIndices = getTupleIndex(tupleConstantsArr, varId);
		 populateMoveOdds(hyperCube, constantToTupleIndices, moveOdds, varId);
		 
		 for(int i=0;i<moveOdds->size();i++) {
			 if((*moveOdds)[i] > classOdd)
				  totalOdd += (*moveOdds)[i];
		 }

		 if(!bestMoveOdds || bestTotalOdd <= totalOdd) {
            bestTotalOdd = totalOdd;
			bestVarId = varId;
			if(bestMoveOdds)
				 delete bestMoveOdds;
			bestMoveOdds = moveOdds;
		 } else {
			delete moveOdds;
		 }
	     
		 IntToIntHashArrayMap::iterator itr;
		 for(itr = constantToTupleIndices->begin(); itr != constantToTupleIndices->end();itr++) {
			  delete itr->second;
		 }
		 delete constantToTupleIndices;
	 }
         
	 int splitVarId = bestVarId;
	 moveOdds = bestMoveOdds;

	 Array<int> *varConstants = hyperCube->getVarConstants(splitVarId);
     Array<int> *movedVarConstants = new Array<int>();

     for(int i=0;i<moveOdds->size();i++) {
	   if((*moveOdds)[i] > classOdd) {
	        movedVarConstants->append((*varConstants)[i]);
	   }
	 }

	 //special case: when all the odds are same as class odd
	 if(movedVarConstants->size() == 0) {
		  for(int i=0;i<varConstants->size()/2;i++) {
	        movedVarConstants->append((*varConstants)[i]);
		  }
	 }

	 /*
	 cout<<"split var = "<<splitVarId<<endl;
	 cout<<"moved constants are "<<endl;
	 printArray(*movedVarConstants, cout);
	 cout<<endl;
     */

	 HyperCube * childHyperCube;
     IntArrayHashArray *childTupleConstantsArr; 
     IntArray *tupleConstants;
	 
	 /*
	 cout<<"Bounding Hypercube is.."<<endl;
	 hyperCube->print(cout);
	 cout<<endl;
	 cout<<endl;
	 
	 cout<<"Child HyperCubes are: "<<endl;
     */
	 
	 for(int i=0;i<2;i++) {
		  childHyperCube = new HyperCube(hyperCube);
		  childTupleConstantsArr = new IntArrayHashArray();
		  
		  varConstants = childHyperCube->getVarConstants(splitVarId);
          /*
		  cout<<"var constants are "<<endl;
		  printArray(*varConstants, cout);
		  cout<<endl;
		  */
		  if(i == 0) {
			   varConstants->clear();
			   varConstants->append(*movedVarConstants);
		  } else {
			   subtractArrayKeepOrder(varConstants, movedVarConstants);
		  }
		  /*
          cout<<"var constants are "<<endl;
		  printArray(*varConstants, cout);
		  cout<<endl;
		  childHyperCube->print(cout);
		  cout<<endl;
          cout<<"-----------------------------"<<endl;
          */

		  for(int i=0;i<tupleConstantsArr->size();i++) {
			   tupleConstants = (*tupleConstantsArr)[i];
			   int constant = (*tupleConstants)[splitVarId];
			   //cout<<"constant is "<<constant<<endl;
			   if(orderedArrayFind(varConstants,constant) >= 0) {
					//cout<<"found the element!"<<endl;
					childTupleConstantsArr->append(tupleConstants);
			   }
		  }
		  
		  createHyperCubes(childHyperCube,childTupleConstantsArr,hyperCubes, noise);
		  delete childHyperCube;
		  delete childTupleConstantsArr;
	 }
	 
	 //cout<<"******************************************"<<endl;

	 delete movedVarConstants;
	 delete moveOdds;
}


//create the hypercubes by merging individual elements from the given array of tuple constants
Array<HyperCube *> *createHyperCubesDT(IntArrayHashArray * const & tupleConstantsArr,
		                               double noise) {
       
	   Array<HyperCube *> *hyperCubes = new Array<HyperCube *>();
	   HyperCube *boundingHyperCube = getBoundingHyperCube(tupleConstantsArr);
	   

	   //cout<<"*****************************************************"<<endl;
       //cout<<"came before create hypercubes..."<<endl;
	   createHyperCubes(boundingHyperCube, tupleConstantsArr, hyperCubes, noise);            
	   delete boundingHyperCube;

	   //cout<<"size of hypercubes in the end is... "<<hyperCubes->size()<<endl;
	   /*
	   cout<<"HyperCubes Are "<<endl;
	   
	   for(int i=0;i<hyperCubes->size();i++) {
			cout<<i<<":"<<endl;
			(*hyperCubes)[i]->print(cout);
			cout<<endl;
	   }
	   cout<<"*****************************************************"<<endl;
	   */

	   //may need to merge them
	   mergeHyperCubes(hyperCubes);

	   return hyperCubes;
}

