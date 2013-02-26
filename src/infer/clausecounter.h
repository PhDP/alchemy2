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
#ifndef _CLAUSE_COUNTER_H_DEC_2007
#define _CLAUSE_COUNTER_H_DEC_2007

#include <ext/hash_set>
#include <ostream>
#include "inferutil.h"

using namespace std;
using namespace __gnu_cxx;

/******************************************************************************/
// Clause Counter
/******************************************************************************/
//typedef hash_map<int, int, HashInt, EqualInt> IntToIntMap;

//this class is for maintaining the counts of for an indexed set of
//clauses (where each entity has a unique id)
class ClauseCounter {
	 
	 public:
	 ClauseCounter() {
	   idToIndex_ = new IntToIntMap();
	   clauseIds_ = new Array<int>;
	   clauseCntsArr_ = new Array<Array<double>*>;
	   //clauseCntsArr_ = new Array<Array<int>*>;
	   dirty_ = true;
	 }

	 ~ClauseCounter(){ 
		  delete idToIndex_;
	      delete clauseIds_;
		  for(int i=0;i<clauseCntsArr_->size();i++) {
		   delete (*clauseCntsArr_)[i];
		  }
		  delete clauseCntsArr_;
	 }

	 int getNumClauses(){ return clauseIds_->size();}
	 

	 const Array<int> * getClauseIds()  const { return clauseIds_;}
	 const Array<double>*  getClauseCounts(int index) const { return (*clauseCntsArr_)[index];}
	 //const Array<int>*  getClauseCounts(int index) const { return (*clauseCntsArr_)[index];}
	 int  getClauseId(int index) { return (*clauseIds_)[index];}

	 void incrementCount(int clauseId, int predId, int clausePredCnt, double cnt) {
		  int index;
		  int id = clauseId;
		  Array<double> *clauseCnts;
		  //Array<int> *clauseCnts;
		  IntToIntMap::iterator itr;
		  itr = idToIndex_->find(id);
		  
		  if(itr == idToIndex_->end()) {
		   index = clauseIds_->size();
		   (*idToIndex_)[id] = index;
		   clauseIds_->append(id);
		  
		   clauseCnts = new Array<double>();
		   //clauseCnts = new Array<int>();
		   clauseCnts->growToSize(clausePredCnt,0);
		   clauseCntsArr_->append(clauseCnts);
		  
		  } else {
		   index = (*idToIndex_)[id];
		  }

		  //increment the count for this clauseId/predId combination
		  clauseCnts = (*clauseCntsArr_)[index];
          //(*clauseCnts)[predId] += (int)cnt;
          (*clauseCnts)[predId] += cnt;
		  dirty_ = true;
	 }
	
	 void multiplyCounts(double scale) {	  
		  Array<double> * clauseCnts;
		  for(int index=0;index<clauseIds_->size();index++) {
			   clauseCnts = (*clauseCntsArr_)[index];
			   for(int predId=0;predId<clauseCnts->size();predId++) {
			    (*clauseCnts)[predId] = (*clauseCnts)[predId] * scale;
			   }
		  }
	 }


	 size_t getHashCode() {	  
		  if(!dirty_)
			   return hashCode_;
		  IntToIntMap::iterator itr;
		  Array<double> * clauseCnts;
		  //Array<int> * clauseCnts;
          int code = 1;
          //double code = 1;
		  for(int index=0;index<clauseIds_->size();index++) {
			   code = 31*code + (*clauseIds_)[index];
			   clauseCnts = (*clauseCntsArr_)[index];
			   for(int predId=0;predId<clauseCnts->size();predId++) {
			    code = 31*code + (int)(*clauseCnts)[predId];
			   }
		  }
		  dirty_ = false;
		  hashCode_ = (size_t)code;
		  return hashCode_;
	 }

     ostream & print(ostream & out) {
		  for(int index=0;index<clauseIds_->size();index++) {
			   out<<(*clauseIds_)[index]<<" : ";
			   printArray(*(*clauseCntsArr_)[index], out);
               out<<" ** ";
		  }
          return out;
	 }


	 private:
		IntToIntMap *idToIndex_;
		
		//NOTE: Assumes that clauseIds_ are stored in a pre-defined order for
		//all the ClauseCounter instances. This is critical for making
		//sure that hashCode/Equal functions on two instances of this
		//class work as expected
		
		Array<int> * clauseIds_;
		Array<Array<double>*> * clauseCntsArr_;
		//Array<Array<int>*> * clauseCntsArr_;
		bool dirty_;
		size_t hashCode_;
};


class HashClauseCounter
{
 public:
  size_t operator()(ClauseCounter *cc) const  { 
	   return cc->getHashCode();
  }
};


class EqualClauseCounter
{
 public:
  bool operator()(ClauseCounter * const & cc1, ClauseCounter * const & cc2) const {
	   bool same;
	   const int *items1, *items2;
	   //const int *cnts1, *cnts2;
	   //const double *cnts1, *cnts2;
	   const Array<double> *cnts1, *cnts2;
	   
	   int size1, size2;
	   
	   size1 = cc1->getNumClauses();
	   size2 = cc2->getNumClauses();

	   if(size1 != size2) return false;

	   //first check if the clause ids match
       items1 = (cc1->getClauseIds())->getItems();
       items2 = (cc2->getClauseIds())->getItems();

	   same = memcmp(items1, items2, size1*sizeof(int))==0;
	   if(!same)
			return same;

	   double epsilon = 1e-6;
	   //now check if the clause cnts match for each of the indices 
	   //(no need check the index sizes again - they must be same at this point)
       for(int index=0;index<size1;index++) {
			cnts1 = cc1->getClauseCounts(index);
			cnts2 = cc2->getClauseCounts(index);
			for(int i=0;i<cnts1->size();i++) {
				 if(((*cnts1)[i] + epsilon < (*cnts2)[i])|| 
				    ((*cnts1)[i] - epsilon > (*cnts2)[i]))
					  return false;
			}

			/*
			cnts1 = (cc1->getClauseCounts(index))->getItems();
            cnts2 = (cc2->getClauseCounts(index))->getItems();
	        int cntsSize = (cc1->getClauseCounts(index))->size();
			//same = memcmp(cnts1, cnts2, cntsSize*sizeof(int))==0;
			same = memcmp(cnts1, cnts2, cntsSize*sizeof(double))==0;
			
			if(!same)
				 return same;
		   */
			
	   }
	   return true;
  }
};

#endif
