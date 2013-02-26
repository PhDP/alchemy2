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
#ifndef _HYPERCUBE_REVERSEINDEX_H_MAY_2008
#define _HYPERCUBE_REVERSEINDEX_H_MAY_2008

#include <ext/hash_map>
#include <ostream>

#include "util.h"
#include "mrf.h"
#include "array.h"
#include "hashint.h"
#include "inferutil.h"

//used in this class

class HyperCubeReverseIndex {

 public:

	 HyperCubeReverseIndex() {
	   init();
	 }

	 HyperCubeReverseIndex(Array<HyperCube *> * const & hyperCubes) {
		  init();
		  for(int i=0;i<hyperCubes->size();i++) {
			   hyperCubes_->append((*hyperCubes)[i]);
		  }
	 }

	 void init() {
		  hasIndex_ = false;
		  hyperCubes_ = new Array<HyperCube *>();
          constantToHCIndicesArr_ = new Array<IntToIntHashArrayMap *>();
	 }

     ~HyperCubeReverseIndex() {
          /*
		  HyperCube *hyperCube;
		  for(int i=0;i<hyperCubes_->size();i++) {
			   hyperCube = (*hyperCubes_)[i];
			   hyperCube->deleteVarConstants();
			   delete hyperCube;
		  }*/
		  
		  delete hyperCubes_;
		  deleteIndex();
		  delete constantToHCIndicesArr_;
	 }
     
	 int getNumHyperCubes(){return hyperCubes_->size();}

	 //add a hypercube
	 bool addHyperCube(HyperCube * const & hyperCube) {
		  int index = hyperCubes_->append(hyperCube);
		  if(hasIndex_) {
			   addHyperCubeToIndex(hyperCubes_->size()-1);
		  }
		  
		  //retrun true if the addition was successful
		  if(index >=0)
			   return true;
		  else 
			   return false;
	 }

	 //grow the index to the newsize
	 void growIndex(int newsize) {
            int size = constantToHCIndicesArr_->size();
			for(int varId=size; varId<newsize;varId++) {
                constantToHCIndicesArr_->append(new IntToIntHashArrayMap());
			   }
	 }

	 //add the hypercube at hindex to the reverse index
     void addHyperCubeToIndex(int hindex) {
		  //IntHashArray * varConstants;
		  IntArray * varConstants;
		  IntHashArray *hcIndices;
		  HyperCube *hyperCube;
		   
		  IntToIntHashArrayMap *constantToHCIndices;
		  IntToIntHashArrayMap::iterator itr;
		  int constant;

		  hyperCube = (*hyperCubes_)[hindex];
		  int varCnt = hyperCube->getVarCount();
			   
	      //grow the index to varCnt+1
		  growIndex(varCnt+1);
		  for(int varId=1;varId<=varCnt;varId++) {
		     constantToHCIndices = (*constantToHCIndicesArr_)[varId];
			 varConstants =  hyperCube->getVarConstants(varId);
             
			 //add this hindex to reverse index for each constant in turn
		     for(int i=0;i<varConstants->size();i++) {
			  constant = (*varConstants)[i];
			  itr = constantToHCIndices->find(constant);
			  if(itr == constantToHCIndices->end()) {
				hcIndices = new IntHashArray();
				(*constantToHCIndices)[constant] = hcIndices;
			  } else {
				 hcIndices = itr->second;
			 }
			  hcIndices->append(hindex);
		  }
	   }
	 }

	 void createIndex() {
		  deleteIndex();
		  hasIndex_ = true;
		  for(int hindex=0;hindex<hyperCubes_->size();hindex++) {
			   addHyperCubeToIndex(hindex);
		  }
	 }

	 //delete the index
	 void deleteIndex() {
		  hasIndex_ = false;
		  IntHashArray *hcIndices;
		  IntToIntHashArrayMap *constantToHCIndices;
		  IntToIntHashArrayMap::iterator itr;
          
		  //for(int varId=1;varId<constantToHCIndicesArr_->size();varId++) {
		  for(int varId=0;varId<constantToHCIndicesArr_->size();varId++) {
            constantToHCIndices = (*constantToHCIndicesArr_)[varId];
			
			for(itr = constantToHCIndices->begin();
				itr != constantToHCIndices->end();
				itr++) {
				 hcIndices = itr->second;
				 delete hcIndices;
			}
			constantToHCIndices->clear();
			delete constantToHCIndices;
		  }
		  constantToHCIndicesArr_->clear();
	 }
	 
	 void deleteAllHyperCubes() {
		  HyperCube *hyperCube;
		  Array<HyperCube *> *allHyperCubes = getAllHyperCubes();
		  for(int i =0;i<allHyperCubes->size();i++) {
           hyperCube = (*allHyperCubes)[i];
		   hyperCube->deleteVarConstants();
		   delete hyperCube;
		  }
		  allHyperCubes->clear();
		  delete allHyperCubes;
	 }
	 

	 //get the hypercube indices corresponding to given set of constants
	 IntHashArray * getHyperCubeIndices(IntArray * const & varConstants,
			                            int varId) {
	  IntHashArray *allHCIndices = new IntHashArray();
	  IntToIntHashArrayMap *constantToHCIndices;
	  if(constantToHCIndicesArr_->size() <= 0)
		   return allHCIndices;

	  //cout<<"came at at point 1"<<endl;
	  constantToHCIndices = (*constantToHCIndicesArr_)[varId];
      IntToIntHashArrayMap::iterator itr;
	  IntHashArray *hcIndices; 
      
	  /*
	  cout<<"constants which have index are => "<<endl;
	  for(itr = constantToHCIndices->begin();
	      itr != constantToHCIndices->end();
		  itr++) {
		   cout<<itr->first<<" ";
	  }
	  cout<<endl;
      */

	  for(int i=0;i<varConstants->size();i++) {
		  int constant = (*varConstants)[i];
          itr = constantToHCIndices->find(constant);
          if(itr == constantToHCIndices->end())
			   continue;
	      hcIndices = itr->second;
		  allHCIndices->append(*hcIndices);
	  }
	  return allHCIndices;
	 }

	 //get the hypercubes corresponding to the given set of indices
	 Array<HyperCube *> * getHyperCubes(IntHashArray * const & hcIndices) {
		  Array<HyperCube *> * hyperCubes = new Array<HyperCube *>();
		  for(int i=0;i<hcIndices->size();i++) {
			   int index = (*hcIndices)[i];
			   hyperCubes->append((*hyperCubes_)[index]);
		  }
	      return hyperCubes;
	 }
	 
	 //get the array (newly created) of all the hypercubes 
	 Array<HyperCube *> * getAllHyperCubes() {
          return new Array<HyperCube *>(*hyperCubes_);
	 }

	 //get the id of the variable this variable refers to (or return varId
	 //if it does not refer to anything)
	 int getReferenceVarId(int varId) {
      if(hyperCubes_->size() == 0)
		   return varId;
	  return (*hyperCubes_)[0]->getReferenceVarId(varId);
	 }


 private:
	  Array<HyperCube *> * hyperCubes_;
	  Array<IntToIntHashArrayMap *> *constantToHCIndicesArr_;
	  bool hasIndex_;
};


#endif

