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
#ifndef _HYPERCUBE_H_MAY_2008
#define _HYPERCUBE_H_MAY_2008

#include <ext/hash_set>
#include <ext/hash_map>
#include <ostream>

#include "util.h"
#include "array.h"
#include "hashint.h"
#include "inferutil.h"


using namespace std;
using namespace __gnu_cxx;


class HyperCube {

public:

     HyperCube() {
	  init();
	 }
     
	 HyperCube(int varCnt) {
	  init();
	  varConstantsArr_->growToSize(varCnt+1,NULL);
	 }

	 //input is one constant for each variable
	 HyperCube(Array<int> * const & constants) {
		  init();
		  
		  /*
		  for(int varId=0;varId<constants->size();varId++) {
			   constants_->append((*constants)[varId]);
		  }*/

		  // start comment
		  
		  //IntArray * varConstants;
		  IntArray * varConstants;
	      varConstantsArr_->growToSize(constants->size(),NULL);
		  
		  for(int varId=1;varId<constants->size();varId++) {
			   //varConstants = new IntArray();
			   varConstants = new IntArray();
			   varConstants->append((*constants)[varId]);
			   varConstants->compress();
			   setVarConstants(varConstants,varId);
		  }
		  //end comment
	 }

	 HyperCube(HyperCube * const & hyperCube) {
		  //exit(0);
		  
		  init();
		  //Array<int> * constants = hyperCube->getConstants();
		  //IntArray * constants = hyperCube->getConstants();
		  /*
		  for(int varId=0;varId<constants->size();varId++) {
			   constants_->append((*constants)[varId]);
		  }*/
		  
		  //
		  //start comment
	  
		  int varCnt = hyperCube->getVarCount();
		  varConstantsArr_->growToSize(varCnt+1,NULL);
		  IntArray *varConstants;
		  for(int varId=1;varId<=varCnt;varId++) {
               varConstants = hyperCube->getVarConstants(varId);
			   if(varConstants)
			    setVarConstants(new IntArray(*varConstants),varId);
		  }
		  //end comment
	
	 }

	 //currently, we only allow copying from pointers - need to fix this at some point!!
	 HyperCube(HyperCube const & hyperCube) {
		  cout<<"Operation Not Supported!!"<<endl;
		  assert(false);
	 }

	 //IntArray * getConstants() {return constants_;}
	 //Array<int> * getConstants() {return constants_;}

     ~HyperCube() {
	      deleteVarConstants();
		  //delete constants_;
		  delete varConstantsArr_;
	 }
	
	 void init() {
	   //constants_ = new Array<int>();
	   //constants_ = new IntArray();
	   
	   dirty_ = true;
	   varConstantsArr_ = new Array<IntArray *>();
	 }

	 void deleteVarConstants() {
		  dirty_ = true;
	      for(int i=0;i<varConstantsArr_->size();i++) {
		  	   if(!(*varConstantsArr_)[i])
					continue;
			   (*varConstantsArr_)[i]->clear();
		       delete (*varConstantsArr_)[i];
			   (*varConstantsArr_)[i] = NULL;
		  }
		  varConstantsArr_->clear();
	 }
    
	 //get the number of variables in this hypercube
	 int getVarCount() const {
		  return (varConstantsArr_->size()-1);
	 }
	
	 int getNumConstants() {
		  int numConstants = 0;
          int varCnt = getVarCount();
		  for(int var=1;var<=varCnt;var++) {
			   numConstants += (*varConstantsArr_)[var]->size();
		  }
		  return numConstants;
	 }


	 //get the number of ground tuples represented by this hypercube
	 int getNumTuples() const {
		  int varCnt = getVarCount();
		  int tupleCnt = 1;
          for(int var=1;var<=varCnt;var++) {
            tupleCnt = tupleCnt * (*varConstantsArr_)[var]->size();
		  }
		  return tupleCnt;
	 }

	 IntArray * getVarConstants(int var) const {
		  return (*varConstantsArr_)[var];
	 }
	 
	 int getVarConstantCount(int var) const {
		  return (*varConstantsArr_)[var]->size();
	 }

	 //get the constant tuple represented by the given array of indices
	 Array<int> * getTupleConstantsAtIndex(Array<int> * const & indices) {
          int varCnt = getVarCount();
		  Array<int> * tupleConstants = new Array<int>();
		  //adding the dummy variable
		  tupleConstants->append(0);
		  
		  for(int varId=1;varId<=varCnt;varId++) {
			 int index = (*indices)[varId]; 	
			 tupleConstants->append((*(*varConstantsArr_)[varId])[index]);
		  }
		  return tupleConstants;
	 }

	 bool hasDuplicateConstants(IntArray * const & tupleConstants,
			                    IntArrayHashArray * const & neqConstraints) {
		IntArray *constraint;  
		if(!neqConstraints)
			 return false;
	    
		for(int i=0;i<neqConstraints->size();i++) {
         constraint = (*neqConstraints)[i];
         int varId1 = -(*constraint)[0];
		 int varId2 = -(*constraint)[1];
		 if((*tupleConstants)[varId1] == (*tupleConstants)[varId2]) {
			  if((*tupleConstants)[varId1] >= 0)
				   return true;
		 }
		}
		return false;
	 }


	 void replaceRefVarIdsByConstants(IntArray * const & tupleConstants) {
		  int refVarId;
		  int constId;
		  for(int varId=1;varId<tupleConstants->size();varId++) {
			   constId = (*tupleConstants)[varId];
			   if(constId < 0) {
					refVarId = -constId;
					assert(refVarId <= getVarCount());
					(*tupleConstants)[varId] = (*tupleConstants)[refVarId];
			   }
		  }
	 }

	 //get all the constant tuples represented by this hypercube
	 Array<Array<int> *>* getTuples(IntArrayHashArray * const & neqConstraints) {
         Array<Array<int>*> * tupleConstantsArr = new Array<Array<int> *>;
         IntArray *tupleConstants;
		 bool isDuplicate;
		 int varCnt = getVarCount();
		 Array<int> * stack = new Array<int>();
		 int breakPos;

         //initialize the stack
		 for(int i=0;i<=varCnt;i++) {
			  stack->append(0);
		 }
		 
		 //iterate over each combination one by one
		 while(true) {
		   tupleConstants = getTupleConstantsAtIndex(stack);
		   isDuplicate = hasDuplicateConstants(tupleConstants, neqConstraints);
		   
		   if(isDuplicate) {
				delete tupleConstants;
		   } else {
			 replaceRefVarIdsByConstants(tupleConstants);
			 tupleConstantsArr->append(tupleConstants);
		   }

		   breakPos = varCnt + 1;
		   //now backtrack to get the next combination
		   for(int var=1;var<=varCnt;var++) {
				if((*stack)[var] < getVarConstantCount(var)-1) {
					 breakPos = var;
					 break;
				}
		   }
		   //done with all the elements
		   if(breakPos > varCnt)
				break;
		   
		   int elem = (*stack)[breakPos];
           (*stack)[breakPos] = elem + 1;
		   
		   for(int var=breakPos-1;var>0;var--) {
            (*stack)[var] = 0;
		   }
	    }   
		delete stack;
        return tupleConstantsArr;	 
     }
	 
	 //get all the constant tuples represented by this hypercube
	 Array<Array<int> *>* getTuples() {
		  return getTuples(NULL);
	 }

	 //get the id of the variable this variable refers to (or return varId
	 //if it does not refer to anything)
     int getReferenceVarId(int varId) {
		  IntArray *varConstants = getVarConstants(varId);
		  if(varConstants == NULL)
			   return varId;
		  if(varConstants->size() != 1)
			   return varId;
		  int refVarId = -(*varConstants)[0];
		  if(refVarId > 0) {
			   assert(refVarId != varId);
			   return refVarId;
		  } else {
			   return varId;
		  }
	 }

//get the relevant set of constraints for this hypercube
IntArrayHashArray *getRelevantConstraints(IntArrayHashArray * const & inpConstraints) {	 
	 
     IntArrayHashArray *constraints = new IntArrayHashArray();
	 IntArray *constraint, *inpConstraint;
	 int varId1, varId2, refVarId1, refVarId2;

	 for(int i=0;i<inpConstraints->size();i++) {
        inpConstraint = (*inpConstraints)[i];
        varId1 = -(*inpConstraint)[0];
		varId2 = -(*inpConstraint)[1];
		refVarId1 = getReferenceVarId(varId1);
		refVarId2 = getReferenceVarId(varId2);
		if(refVarId1 == refVarId2)
			 continue;

		if(refVarId1 > refVarId2) {
         int tmp = refVarId1;
		 refVarId1 = refVarId2;
		 refVarId2 = tmp;
		}
		constraint = new IntArray();
		constraint->append(-refVarId1); 
		constraint->append(-refVarId2);
	    /*
		cout<<"Constraint is ="<<endl;
		printArray(*constraint,cout);
		cout<<endl;
		*/
		constraints->append(constraint);
	  }
      return constraints;
}
   
	 //check if the hypercube has a singleton duplicte
     bool hasDuplicateConstants(IntArrayHashArray * const & neqConstraints) {
		  bool foundDuplicate = false;
		  int varCnt = getVarCount();
		  IntArray *constraint = NULL;
		  IntArray *varConstants1, *varConstants2;
		  IntArrayHashArray *relevantConstraints = NULL;
		  for(int varId1=1;varId1<=varCnt;varId1++) {
		    varConstants1 = getVarConstants(varId1);
			if((varConstants1->size() != 1) || ((*varConstants1)[0] < 0))
				 continue;
			
			for(int varId2=varId1+1;varId2<=varCnt;varId2++) {
		     varConstants2 = getVarConstants(varId2);
			 if((varConstants2->size() != 1) || ((*varConstants2)[0] < 0))
				 continue;
			 if((*varConstants1)[0] != (*varConstants2)[0])
				  continue;

			 //ok found the same set of constants. now, need to make sure that the
			 //varId pair belongs the constraint set
			 if(!relevantConstraints) 
			  relevantConstraints = getRelevantConstraints(neqConstraints);
			 
			 if(!constraint) {
			  constraint = new IntArray();
			 } else {
              constraint->clear();
			 }
			constraint->append(-varId1);
			constraint->append(-varId2);
			/*
			cout<<"Superclause constraint is "<<endl;
			printArray(*(*neqConstraints)[0],cout);
			cout<<endl;
			cout<<"varId1 = "<<varId1<<", varId2 = "<<varId2<<endl;
            */
			if(neqConstraints->find(constraint) >= 0) {
				 foundDuplicate = true;
				 break;
			} 
			/*
			else {
				 cout<<"wow! this is the problem!!"<<endl;
				 this->print(cout);
				 cout<<endl;
		    }*/
		   }
		  } 
		  if(constraint)
			   delete constraint;
          if(relevantConstraints) {
			   relevantConstraints->deleteItemsAndClear();
			   delete relevantConstraints;
		  }
		  return foundDuplicate;
	 }


	 bool isReferenceVar(int varId) {
		  IntArray *varConstants = getVarConstants(varId);
		  return (varConstants != NULL && varConstants->size() == 1 && (*varConstants)[0] < 0);
	 }


	 bool isEmpty() {
	  int varCnt = getVarCount();
	  for(int varId=1;varId<=varCnt;varId++) {
		   if(getVarConstantCount(varId) == 0) {
				return true;
		   }
	  }
	  return false;
	 }

	 bool isSubsetOf(IntArray * const & inpVarConstants, int varId) {
		  int constant;
		  IntArray * varConstants = (*varConstantsArr_)[varId];
		  if(varConstants->size() > inpVarConstants->size()) 
			   return false;
		  for(int i=0;i<varConstants->size();i++) {
			   constant = (*varConstants)[i];
			   //if(inpVarConstants->find(constant) < 0)
			   if(orderedArrayFind(inpVarConstants,constant) < 0)
					return false;
		  }
		  return true;
	 }

     bool isSubsetOf(HyperCube * const & hyperCube) {
		  IntArray * inpVarConstants;
		  int varCnt = getVarCount();
		  for(int varId=1;varId<=varCnt;varId++) {
			   inpVarConstants = hyperCube->getVarConstants(varId);
			   if(!isSubsetOf(inpVarConstants, varId))
					return false;
		  }
		  return true;
	 }
	 
	 bool hasIntersection(IntArray * const & inpVarConstants, int varId) {
		  int constant;
		  IntArray * varConstants = (*varConstantsArr_)[varId];
		  //assert(varConstants);
		  for(int i=0;i<varConstants->size();i++) {
			   constant = (*varConstants)[i];
			   //if(inpVarConstants->find(constant) >= 0)
			   if(orderedArrayFind(inpVarConstants,constant) >= 0)
					return true;
		  }
		  return false;
	 }

	 
	 bool hasIntersection(HyperCube * const & hyperCube, bool print=false) {
		  IntArray * inpVarConstants;
		  int varCnt = getVarCount();
		  //assert(varCnt == hyperCube->getVarCount());
		  
		  if(print) {
		   cout<<"-------------------------------------------------"<<endl;
		   cout<<"Getting the intersection of"<<endl;
		   cout<<"A=>"<<endl;
		   this->print(cout);
		   cout<<endl;
		   cout<<"B=>"<<endl;
		   hyperCube->print(cout);
		   cout<<endl;
		   cout<<"-------------------------------------------------"<<endl;
		  }

		  for(int varId=1;varId<=varCnt;varId++) {
			   inpVarConstants = hyperCube->getVarConstants(varId);
			   //assert(inpVarConstants);
			   //each dimension should have a non-zero intersection
			   if(!hasIntersection(inpVarConstants, varId)) {
                    //cout<<"returning false"<<endl;
		            //cout<<"-------------------------------------------------"<<endl;
					return false;
			   }
		  }
          //cout<<"returning true"<<endl;
		  //cout<<"-------------------------------------------------"<<endl;
		  return true;
	 }
	
	 //check whether var Constants at two positions are identical
	 bool hasSameVarConstants(int varId1, int varId2) {
		  IntArray * varConstants1 = getVarConstants(varId1);
		  IntArray * varConstants2 = getVarConstants(varId2);
		  return orderedArrayAreSame(varConstants1, varConstants2);
	 }

	 //check whether var Constants at two positions are disjoint
	 bool hasDisjointVarConstants(int varId1, int varId2) {
		  IntArray * varConstants1 = getVarConstants(varId1);
		  return !hasIntersection(varConstants1, varId2);
	 }


	 //set the ref id of the given variable 
     void setReferenceVarId(int varId, int refId) {
		  dirty_ = true;
		  assert(isReferenceVar(varId));
		  IntArray *varConstants = getVarConstants(varId);
		  (*varConstants)[0] = -refId;
	 }


	 void setToEmpty() {
	  dirty_ = true;
	  int varCnt = getVarCount();
	  for(int varId=1;varId<=varCnt;varId++) {
		  (*varConstantsArr_)[varId]->clear();
	  }
	 }

	 //add a new variable with the given set of constants
     void addNewVariable(IntArray * const & varConstants) {
		  dirty_ = true;
		  varConstantsArr_->append(varConstants);
	 }

	 //Adds the input varConstants to the list of constants already present
     void addVarConstants(IntArray * const & inpVarConstants, int varId) {
		  dirty_ = true;
          IntArray *varConstants = getVarConstants(varId);
		  unionArrayKeepOrder(varConstants,inpVarConstants);
	 }

	 //set the var constants at the given index (varId) to the given array
	 //of constants. Does not delete the previous set of constants - caller
	 //should take care of it
	 void setVarConstants(IntArray * const & varConstants, int varId) {
		  dirty_ = true;
          (*varConstantsArr_)[varId] = varConstants;
	 }

	 
	 //set the var constants at the given index (varId) to the given array
	 //of constants. Deletes the previous set of constants
	 void setVarConstantsAndDeleteOld(IntArray * const & varConstants, int varId) {
		  dirty_ = true;
          delete (*varConstantsArr_)[varId];
		  setVarConstants(varConstants, varId);
	 }
    
	 //canocialize this hypercube according to the given mapping. As we ground away
	 //the variables, we need to maintain a cnt of how many copies of the simplified
	 //hypercube are created
	 int canonicalize(Array<int> * const &varIdToCanonicalVarId) {
        dirty_ = true;
		
		//Need to work on the implicit representation
		updateToImplicitTiedRepresentation();

		int cnt = 1;
		Array<IntArray *> *newVarConstantsArr = new Array<IntArray *>();
        
		//first entry is null
		newVarConstantsArr->growToSize(1,NULL);

		IntArray *varConstants;

		int canonicalVarId;
        int varCnt = getVarCount();
		
		for(int varId=1;varId<=varCnt;varId++) {
          varConstants = (*varConstantsArr_)[varId];
          canonicalVarId = (*varIdToCanonicalVarId)[varId];
          //cout<<"varId = "<<varId<<", canonicalVarId = "<<canonicalVarId<<endl;
		  if(canonicalVarId < 0) {
			   //since hypercube is assumed to be grounded, varConstants should not be null
			   assert(varConstants);
			   cnt = cnt * varConstants->size();
			   delete varConstants;
			   
			   // if(varConstants)
			   // delete varConstants;
               continue;
		  }

          if(newVarConstantsArr->size() < canonicalVarId+1)
               newVarConstantsArr->growToSize(canonicalVarId + 1);
          (*newVarConstantsArr)[canonicalVarId] = varConstants;
     }

	  delete varConstantsArr_;
      varConstantsArr_ = newVarConstantsArr;

	  //back to explicit representation
	  updateToExplicitTiedRepresentation();
	  return cnt;
	 }

	 //the tied variables by represented by storing the same copy of
	 //varConstants at the tied positions
     void updateToImplicitTiedRepresentation() {
		dirty_ = true;
		IntArray *varConstants;
		int varCnt = getVarCount();
		for(int varId = 1; varId <= varCnt; varId++) {
         int refId = getReferenceVarId(varId);
	     if(refId == varId)
			 continue;
		 //cout<<"wow!!, kya chal raha hai 1 1 1 1 !! came here.!!!"<<endl;
         varConstants = getVarConstants(refId);
		 setVarConstantsAndDeleteOld(varConstants,varId);
		}
	 }
     
	 //tied variables are stored by storing an explicit reference
	 void updateToExplicitTiedRepresentation() {
		dirty_ = true;
		IntArray *varConstants, *varConstants1, *varConstants2;
		int varCnt = getVarCount();
		for(int varId1 = 1; varId1 <= varCnt; varId1++) {
			 varConstants1 = getVarConstants(varId1);
			 if(!varConstants1)
				  continue;
			 for(int varId2 = varId1+1;varId2 <= varCnt; varId2++) {
               varConstants2 = getVarConstants(varId2);
			   if(varConstants1 == varConstants2) {
				   //cout<<"wow!!, kya chal raha hai 2 2 2 2 !! came here.!!!"<<endl;
				   //cout<<"varId1 = "<<varId1<<" ** varId2 = "<<varId2<<endl;
				   //cout<<varConstants1<<" ** "<<varConstants2<<endl;
				   varConstants = new IntArray();
				   varConstants->append(-varId1);
				   setVarConstants(varConstants,varId2);
			  }
			 }
		}
	 }
			 
	 //Intersects the constants at varId. Returns true if intersection is non null
	 void intersect(IntArray * const & inpVarConstants,
			        int varId) {
		  dirty_ = true;
		  IntArray *varConstants = (*varConstantsArr_)[varId];
		  //get the intersection of varConstants
		  intersectArrayKeepOrder(varConstants,inpVarConstants);
	 }

     void intersect(HyperCube * const & hyperCube) {
		  dirty_ = true;
		  int varCnt = getVarCount();
		  assert(varCnt == hyperCube->getVarCount());
          for(int varId = 1; varId <= varCnt; varId++) {
			   intersect(hyperCube->getVarConstants(varId), varId);
			   if(getVarConstantCount(varId) <= 0) {
					setToEmpty();
					break;
			   }
		  }
	 }

	 //intersect the constants at two variable positions
     void intersect(int varId1, int varId2) {
		  IntArray *varConstants1, *varConstants2;
		  varConstants1 = getVarConstants(varId1);
		  intersect(varConstants1, varId2);
		  varConstants2 = getVarConstants(varId2);
		  intersect(varConstants2, varId1);
	 }

	 //get the set of hypercubes obtained by set subtraction of the input
	 //hypercube from this hypercube 
     Array<HyperCube *> * getMinus(HyperCube * const & hyperCube) {
		
		 /*
		 cout<<"--------------------------------------------------------"<<endl;
         cout<<"Getting Minus of:"<<endl;
		 cout<<"A =>"<<endl;
		 print(cout);
		 cout<<endl;
		 
		 cout<<"B =>"<<endl;
		 hyperCube->print(cout);
		 cout<<endl;
         */
		 
		 Array<HyperCube *> * minusHyperCubes = new Array<HyperCube *>();
		 HyperCube *primaryMinusHyperCube, *minusHyperCube;
		 IntArray * varConstants;

		 int varCnt = getVarCount();
		 assert(varCnt == hyperCube->getVarCount());
		 
		 //check if this hypercube is same as or a subset of given hypercube
		 if(same(hyperCube) || isSubsetOf(hyperCube)) {
			  /*
			  cout<<"yes, this is a subset.."<<endl;
		      cout<<"--------------------------------------------------------"<<endl;
			  */
			  return minusHyperCubes;
		 }
     
         //check if the two hypercubes do not intersect with each other
		 if(!hasIntersection(hyperCube)) {
			  minusHyperCubes->append(new HyperCube(this));
			  return minusHyperCubes;
		 }


		 //** If it came at this point, then there is a non trivial minus component **/
		 
		 //first get the hypercube obtained by subtracting out each dimension
		 primaryMinusHyperCube = new HyperCube(varCnt);
		 for(int varId=1; varId<=varCnt; varId++) {
              varConstants = new IntArray(*getVarConstants(varId));
              subtractArrayKeepOrder(varConstants,hyperCube->getVarConstants(varId));
			  primaryMinusHyperCube->setVarConstants(varConstants, varId);
		 }
		 
		/* 
		 cout<<"Primary Minus Hypercube is = "<<endl;
         primaryMinusHyperCube->print(cout);
		 cout<<endl;
		 */
         
		 //check if primary minus hypercube is same as this hypercube which 
		 //means the intersection of two hypercubes is null
		 
		 if(this->same(primaryMinusHyperCube)) {
           minusHyperCubes->append(primaryMinusHyperCube);
		   return minusHyperCubes;
		   //cout<<"--------------------------------------------------------"<<endl;
		 }
          
		 //now, obtain all the minus hyperCubes
		 for(int index=1; index<=varCnt; index++) {
			  minusHyperCube = new HyperCube(this);
			  
			  for(int varId=1;varId<index;varId++) {
				   //varConstants = new IntArray(*hyperCube->getVarConstants(varId));
				   varConstants = hyperCube->getVarConstants(varId);
				   minusHyperCube->intersect(varConstants, varId);
			  }
             
			  varConstants = new IntArray(*primaryMinusHyperCube->getVarConstants(index));
			  minusHyperCube->setVarConstantsAndDeleteOld(varConstants, index);
			 
			  /*
			  for(int varId=index+1;varId<=varCnt;varId++) {
				   varConstants = new IntArray(*getVarConstants(varId));
				   minusHyperCube->setVarConstants(varConstants, varId);
			  }
			  */
			  if(!minusHyperCube->isEmpty())
			   minusHyperCubes->append(minusHyperCube);
			  else
			   delete minusHyperCube;
		 }
	     
		 //primaryMinusHyperCube->deleteVarConstants();
		 delete primaryMinusHyperCube;
		 //cout<<"--------------------------------------------------------"<<endl;
		 return minusHyperCubes;
	 }

     bool same(HyperCube * const & hyperCube) {

		 HyperCube *hc1, *hc2;
         hc1 = this;
		 hc2 = hyperCube;
          
		 if(hc1->getHashCode() != hc2->getHashCode())
			   return false;

		 // this exactly a copy of the function in EqualHyperCube class
		 
		 /************ Start Copy *******************************/
		 const int *items1, *items2;
	     int constCnt1, constCnt2;
	     int varCnt1 = hc1->getVarCount();
	     int varCnt2 = hc2->getVarCount();
	  
	     if(varCnt1 != varCnt2)
			return false;

	     for(int var=1;var<=varCnt1;var++) {

		 constCnt1 = hc1->getVarConstantCount(var);
		 constCnt2 = hc2->getVarConstantCount(var);

		 if(constCnt1 != constCnt2)
			  return false;

	     items1 = hc1->getVarConstants(var)->getItems();
		 items2 = hc2->getVarConstants(var)->getItems();
		 
		 bool same = (memcmp(items1, items2,constCnt1*sizeof(int)) == 0);
		 if(!same)
          return false;
	   }

	   return true;
	
	   /************ End Copy *******************************/
  }
	 
	 size_t getHashCode() {	  
		  IntArray *varConstants;
		  if(!dirty_)
			   return hashCode_;
		  int code = 1;
		  
		  //experimental
	
		  /*
		  int cnt = constants_->size();
		  for(int i=0;i<=cnt;i++) {
				code = 31*code + (*constants_)[i];
		  }*/
         

		  //start comment
		  
		  code = 1;
		  int varCnt = getVarCount();
		  for(int var=1;var<=varCnt;var++) {
			   varConstants = (*varConstantsArr_)[var];
			   //if any of the dimensions is zero, this is basically
			   //an empty hypercube
			   if(varConstants->size() == 0) {
					code = 1;
					break;
			   }
			 for(int i=0;i<varConstants->size();i++) { 
				  code = 31*code + (*varConstants)[i];
		     }
		  }
		  
		 //end comment

		  dirty_ = false;
		  hashCode_ = (size_t)code;
		  return hashCode_;
	 }

	 //print the constants (for each variable separately)
     ostream & print(ostream & out) {
		  IntArray * varConstants;
          int varCnt = getVarCount();
		  for(int var=1;var<=varCnt;var++) {
			   varConstants = (*varConstantsArr_)[var];
			   out<<"var "<<var<<" => ";
			   if(varConstants == NULL) {
					out<<"NULL"<<endl;
					continue;
			   }
			 for(int i=0;i<varConstants->size();i++) { 
				  out<<(*varConstants)[i]<<" ";
			 }
			 out<<endl;
		  }
		  return out;
	 }
	 
 private:
 
 //Each element in the array is the array of constants for a particular variable
 //(i.e. subdomain for that variable) participating in this hypercube. 
 //Element at index 0 is ignored
  Array<IntArray *> *varConstantsArr_;
  
  //IntArray *constants_;
  //Array<int> *constants_;
  //varConstantsArr_;
  
  int hashCode_;
  bool dirty_;
}; 


class HashHyperCube
{
 public:
  size_t operator()(HyperCube *hc) const  { 
	   return hc->getHashCode();
  }
};


class EqualHyperCube
{
 public:
  bool operator()(HyperCube * const & hc1, HyperCube * const & hc2) const {
	   const int *items1, *items2;
	   int constCnt1, constCnt2;
	   int varCnt1 = hc1->getVarCount();
	   int varCnt2 = hc2->getVarCount();
	  
	   if(varCnt1 != varCnt2)
			return false;

	   for(int var=1;var<=varCnt1;var++) {

		 constCnt1 = hc1->getVarConstantCount(var);
		 constCnt2 = hc2->getVarConstantCount(var);

		 if(constCnt1 != constCnt2)
			  return false;

	     items1 = hc1->getVarConstants(var)->getItems();
		 items2 = hc2->getVarConstants(var)->getItems();
		 
		 bool same = (memcmp(items1, items2,constCnt1*sizeof(int)) == 0);
		 if(!same)
          return false;
	   }

	   return true;
  }
};

typedef HashArray<HyperCube*, HashHyperCube, EqualHyperCube> HyperCubeHashArray;

typedef hash_map<HyperCube*, HyperCube*, HashHyperCube, EqualHyperCube> HyperCubeToHyperCubeMap;
typedef hash_map<HyperCube*, HyperCubeHashArray*, HashHyperCube, EqualHyperCube> HyperCubeToHyperCubeHashArrayMap;
typedef hash_map<HyperCube*, Array<HyperCube *> *, HashHyperCube, EqualHyperCube> HyperCubeToHyperCubeArrayMap;

#endif

