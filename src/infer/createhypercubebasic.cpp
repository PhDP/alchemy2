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

#include "hypercube.h"
#include "hypercuberefinement.h"

using namespace std;
using namespace __gnu_cxx;

//global variable - used only in this file
IntArray * _defaultVarConstants_ = new IntArray();

//manipulate the indices in superclause
void retrieveIndexItr(HyperCube * const & hyperCube,
		  HyperCubeToHyperCubeHashArrayMap * const & indexToHyperCubes,
		  int & varId,
		  HyperCubeToHyperCubeHashArrayMap::iterator & itr) {

	 IntArray * origVarConstants = hyperCube->getVarConstants(varId);
	 //replace the constants at varId by default set of constants
	 hyperCube->setVarConstants(_defaultVarConstants_, varId);
	 itr = indexToHyperCubes->find(hyperCube);
	 //replace back the original set of constants
	 hyperCube->setVarConstants(origVarConstants, varId);
     
	 /*
	 cout<<"**************************************************"<<endl;

	 cout<<"Var Id = "<<varId<<endl;
	 cout<<"Index Constants =>"<<endl;
	 printArray(*origVarConstants,0,cout);
	 cout<<endl;
	 cout<<"Hypercubes in the index are :"<<endl;

	 if(itr != indexToHyperCubes->end()) {
	 HyperCubeHashArray *hyperCubes = itr->second;
	 for(int i=0;i<hyperCubes->size();i++){
		  (*hyperCubes)[i]->print(cout);
		  cout<<endl;
	  }
	 }

	 //All the indices
	 cout<<"All the Indices for Var "<<varId<<" Are "<<endl;
	 HyperCubeToHyperCubeHashArrayMap::iterator myitr;
	 HyperCube *hc;
	 for(myitr = indexToHyperCubes->begin();
	     myitr != indexToHyperCubes->end();
		 myitr++) {
          hc = myitr->first;
		  hc->print(cout);
		  cout<<endl;
	 }
	 */
}

//add the hypercube to the index for the given varId
void addHyperCubeToIndex(HyperCube * const & hyperCube, 
		  HyperCubeToHyperCubeHashArrayMap * const & indexToHyperCubes,
		  HyperCubeHashArray * const & indicesForMerge,
		  int & varId,
		  HyperCubeToHyperCubeHashArrayMap::iterator & itr) {

	 //cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
	 //cout<<"Adding Hypercube......"<<endl;
	 
	 HyperCubeHashArray *hyperCubes;
	 HyperCube *index;
	 retrieveIndexItr(hyperCube, indexToHyperCubes, varId, itr);
	 if(itr == indexToHyperCubes->end()) {
		  index = new HyperCube(hyperCube);
		  //index->setVarConstantsAndDeleteOld(new IntHashArray(*_defaultVarConstants_),varId);
		  index->setVarConstantsAndDeleteOld(new IntArray(*_defaultVarConstants_),varId);
		  hyperCubes = new HyperCubeHashArray();
		  (*indexToHyperCubes)[index] = hyperCubes;
	 } else {
		  index = itr->first;
		  hyperCubes = itr->second;
	 }

	 hyperCubes->append(hyperCube);
	 if(hyperCubes->size() >= 2) {
		  indicesForMerge->append(index);
	 }
	 //cout<<"--------------------------------------------------------"<<endl;
}

//remove the hypercube from the index of the given varId
void removeHyperCubeFromIndex(HyperCube * hyperCube, 
		  HyperCubeToHyperCubeHashArrayMap * const & indexToHyperCubes,
		  HyperCubeHashArray * const & indicesForMerge,
		  int & varId,
		  HyperCubeToHyperCubeHashArrayMap::iterator & itr) {
	
	 //cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
	 //cout<<"Removing Stuff......"<<endl;
	 HyperCubeHashArray *hyperCubes;
	 HyperCube *index;
	 retrieveIndexItr(hyperCube, indexToHyperCubes, varId, itr);
	 assert(itr != indexToHyperCubes->end());
	 index = itr->first;
	 hyperCubes = itr->second;

	 hyperCubes->removeInputItemFastDisorder(hyperCube);
	 
	 //remove this index from the list of mergeable indices
	 //if there is only one hypercube left corresponding to 
	 //this index
	 if(hyperCubes->size() < 2) {
		  indicesForMerge->removeInputItemFastDisorder(index);
	 }
	 
	 //remove the index if there is no hypercube corresponding to it now
	 if(hyperCubes->size() < 1) {
         indexToHyperCubes->erase(itr);
		 delete hyperCubes;
		 delete index;
	 }
		  
    /*
	 cout<<"Removed Hypercube =>"<<endl;
	 hyperCube->print(cout);
	 cout<<endl;
	 cout<<"After removing, hypercubes are .."<<endl;
	 for(int i=0;i<hyperCubes->size();i++){
		  (*hyperCubes)[i]->print(cout);
		  cout<<endl;
	 }
	 cout<<"--------------------------------------------------------"<<endl;
	 */
}


//merge the input hypercubes into larger hypercubes (make as large as possible)
void mergeHyperCubes(Array<HyperCube *> * const & inpHyperCubes) {
		  
     if(inpHyperCubes->size() ==0 )
	      return;
	 
	 HyperCube *hyperCube, *index;
      
     hyperCube = (*inpHyperCubes)[0];
	 int varCnt = hyperCube->getVarCount();

	 Array<HyperCubeToHyperCubeHashArrayMap *> *indexToHyperCubesArr = new Array<HyperCubeToHyperCubeHashArrayMap *>();
	 HyperCubeToHyperCubeHashArrayMap::iterator itr;

	 Array<HyperCubeHashArray *> *indicesForMergeArr = new Array<HyperCubeHashArray*>();

	 indexToHyperCubesArr->growToSize(varCnt+1);
	 indicesForMergeArr->growToSize(varCnt+1);
	 for(int varId=1;varId<=varCnt;varId++) {
		  (*indexToHyperCubesArr)[varId] = new HyperCubeToHyperCubeHashArrayMap();
		  (*indicesForMergeArr)[varId] = new HyperCubeHashArray();
	 }

	 //first populate the indices for each variable
	 // - key is all the subdomains in the hypercube
	 // except for the variable being considered
     for(int i=0;i<inpHyperCubes->size();i++) {
       hyperCube = (*inpHyperCubes)[i];
	   //add this hypercube to appropriate indices
	   for(int varId=1;varId<=varCnt;varId++) {
	      addHyperCubeToIndex(hyperCube,(*indexToHyperCubesArr)[varId],(*indicesForMergeArr)[varId],varId, itr);
	   }
	 }

	 HyperCubeHashArray *hyperCubes;
	 HyperCube *hyperCube1, *hyperCube2;
	 int indexVarId;
	 //cout<<"Done creating the Initial Set of hypercubes.."<<endl;
	 //now, start merging the hypercubes
	 while(true) {
		  indexVarId = -1;
		  for(int varId=1;varId<=varCnt;varId++) {
			   if((*indicesForMergeArr)[varId]->size() > 0) {
					indexVarId = varId;
					break;
			   }
		  }

		  //done merging all the hypercubes
		  if(indexVarId < 0)
			   break;

		  index = (*(*indicesForMergeArr)[indexVarId])[0];
		  itr = ((*indexToHyperCubesArr)[indexVarId])->find(index);
		  //assert(itr != (indexToHyperCubesArr*)[indexVarId]->end());
		  hyperCubes = itr->second;
		  hyperCube1 = (*hyperCubes)[0];
		  hyperCube2 = (*hyperCubes)[1];

		  /*
		  cout<<"indexVarId = "<<indexVarId<<endl;
		  cout<<"Index = "<<endl;
		  index->print(cout);
		  cout<<"H1 = "<<endl;
		  hyperCube1->print(cout);
		  cout<<"H2 = "<<endl;
		  hyperCube2->print(cout);
          */

		  //delete these hypercubes from the indices
		  for(int varId=1;varId<=varCnt;varId++) {
			   removeHyperCubeFromIndex(hyperCube1,(*indexToHyperCubesArr)[varId],(*indicesForMergeArr)[varId],
						 varId, itr);
			   removeHyperCubeFromIndex(hyperCube2,(*indexToHyperCubesArr)[varId],(*indicesForMergeArr)[varId],
						 varId, itr);
		  }

		  hyperCube1->addVarConstants(hyperCube2->getVarConstants(indexVarId), indexVarId);
		  hyperCube2->deleteVarConstants(); 
		  delete hyperCube2;

		  //hypercube1 is now the merge of the two hypercube along the dimenstion
		  //varId
		  for(int varId=1;varId<=varCnt;varId++) {
			   addHyperCubeToIndex(hyperCube1,(*indexToHyperCubesArr)[varId],(*indicesForMergeArr)[varId],varId, itr);
		  }
		  
		  /*
		  cout<<"Merged HyperCube = "<<endl;
		  hyperCube1->print(cout);
		  */
	 }

	 //Now, simply read the hypercubes from any one of the indices
	 //add it to the array of hypercubes
	 Array<HyperCube *> *mergedHyperCubes = new Array<HyperCube *>();
     
	 int hcnt = 0;
	 
	 //cout<<"Final Set of hypercubes obtained is "<<endl;
	 for(itr = ((*indexToHyperCubesArr)[1])->begin();
			   itr != ((*indexToHyperCubesArr)[1])->end();
			   itr++ ) {

		  hyperCubes = itr->second;
		  assert(hyperCubes->size() <= 1);

		  for(int i=0;i<hyperCubes->size();i++) {
			   hcnt++;
			   hyperCube = (*hyperCubes)[i];
			   mergedHyperCubes->append(hyperCube);
			   /*
			   cout<<"HyperCube "<<hcnt<<" =>"<<endl;
			   hyperCube->print(cout);
			   cout<<endl;
			   */
		  }
	 }

	 //now, clear up the stuff

	 //to store the keys of the map. They can not be deleted while iterating
	 //over the map. So, we first store them in an array and delete them
	 //later
	 Array<HyperCube *> keysArr;

	 for(int varId=1;varId<=varCnt;varId++) {
		  keysArr.clear();
		  for(itr = ((*indexToHyperCubesArr)[varId])->begin();
					itr != ((*indexToHyperCubesArr)[varId])->end();
					itr++ ) {
			   keysArr.append(itr->first); 
			   delete itr->second;
		  }
		  for(int i=0;i<keysArr.size();i++) {
			   keysArr[i]->deleteVarConstants();
			   delete keysArr[i];
		  }

		  delete (*indexToHyperCubesArr)[varId];
		  delete (*indicesForMergeArr)[varId];
	 }

	 delete indexToHyperCubesArr;
	 delete indicesForMergeArr;

     inpHyperCubes->clear();
     inpHyperCubes->append(*mergedHyperCubes);
	 /*
	 cout<<"*****************************************************"<<endl;
	 cout<<"Merged hyper cubes are "<<endl;
	 for(int i=0;i<mergedHyperCubes->size();i++) {
		  cout<<i<<":"<<endl;
		  (*mergedHyperCubes)[i]->print(cout);
		  cout<<endl<<endl;
	 }
	 cout<<"*****************************************************"<<endl;
	 */
	 delete mergedHyperCubes;
}


//create the basic ground hypercubes - one hypercube for each tuple */
Array<HyperCube *> *createGroundHyperCubes(IntArrayHashArray * const & tupleConstantsArr) {
	 if(tupleConstantsArr->size() == 0)
		  return new Array<HyperCube *>();
	 Array<int> * tupleConstants;
	 tupleConstants = (*tupleConstantsArr)[0];
     HyperCube *hyperCube;
     Array<HyperCube *> *hyperCubes = new Array<HyperCube *>();

	 for(int i=0;i<tupleConstantsArr->size();i++) {
		  hyperCube = new HyperCube((*tupleConstantsArr)[i]);
		  hyperCubes->append(hyperCube);
	 }
	 return hyperCubes;
}


//create the hypercubes by merging individual elements from the given array of tuple constants
Array<HyperCube *> *createHyperCubesBasic(IntArrayHashArray * const & tupleConstantsArr) {

	 /*
	 if(tupleConstantsArr->size() == 0)
		  return new Array<HyperCube *>();
	 Array<int> * tupleConstants;
	 tupleConstants = (*tupleConstantsArr)[0];
     HyperCube *hyperCube;
     Array<HyperCube *> *hyperCubes = new Array<HyperCube *>();


	 for(int i=0;i<tupleConstantsArr->size();i++) {
		  hyperCube = new HyperCube((*tupleConstantsArr)[i]);
		  hyperCubes->append(hyperCube);
	 }*/
     
	 Array<HyperCube *> *hyperCubes;
	 hyperCubes = createGroundHyperCubes(tupleConstantsArr);
	 mergeHyperCubes(hyperCubes);
	 return hyperCubes;
}

