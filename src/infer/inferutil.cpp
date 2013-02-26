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
#include "inferutil.h"

using namespace __gnu_cxx;
using namespace std;

enum Operation {
	 UNION,
	 SUBTRACT,
	 INTERSECT
};

void operateArrayKeepOrder(IntArray * const & arr1, IntArray * const & arr2, Operation oper);

/***************************************************************************/
// Mehods
/***************************************************************************/

/*
//Note: KeepOrder suffix for the method below means the following:
//It is assumed that both arrays are sorted in increasing order
//and method maintains this increasing order in the output
 * */

void intersectHashArray(IntHashArray * const & arr1, IntHashArray * const & arr2) {
	 intersectHashArrayKeepOrder(arr1, arr2);
}

void intersectHashArrayKeepOrder(IntHashArray * const & arr1, IntHashArray * const & arr2) {
	 
	 IntHashArray *arrIntersect = new IntHashArray();
	 int elem;
	 for(int i=0;i<arr1->size();i++) {
		  elem = (*arr1)[i];
		  if(arr2->find(elem) >= 0){ 
			   arrIntersect->append(elem);
		  }
	 }
	 (*arr1) = (*arrIntersect);
	 arr1->compress();
	 delete arrIntersect;
}


//perform the union of two arrays
void unionArrayKeepOrder(IntArray * const & arr1, IntArray * const & arr2) {
	 Operation oper = UNION;
	 operateArrayKeepOrder(arr1, arr2, oper);
}

//perform the intersection
void intersectArrayKeepOrder(IntArray * const & arr1, IntArray * const & arr2) {
	 Operation oper = INTERSECT;
	 operateArrayKeepOrder(arr1, arr2, oper);
}

//perform the subtraction
void subtractArrayKeepOrder(IntArray * const & arr1, IntArray * const & arr2) {
	 Operation oper = SUBTRACT;
	 operateArrayKeepOrder(arr1, arr2, oper);
}


//perform the desired operation on the arrays (union/intersection/subtraction)
//and store the result in first one
void operateArrayKeepOrder(IntArray * const & arr1, IntArray * const & arr2, Operation oper) {
	 IntArray *arrResult = new IntArray();
	 int elem1, elem2;
	 
	 int index1 = 0;
	 int index2 = 0;
	 int size1 = arr1->size();
	 int size2 = arr2->size();
	 while(index1 < size1 && index2 < size2) {
	       elem1 = (*arr1)[index1];
		   elem2 = (*arr2)[index2];
		   if(elem1 == elem2) {
				if(oper == UNION || oper == INTERSECT) {
					 arrResult->append(elem1);
				}
				index1++;
				index2++;
		   }
		   if(elem1 < elem2) {
				if(oper == UNION || oper == SUBTRACT) {
					 arrResult->append(elem1);
				}
				index1++;
		   }
		   if(elem2 < elem1) {
				if(oper == UNION) {
					 arrResult->append(elem2);
				}
				index2++;
		   }
	 }
	 
	 if(oper == UNION || oper == SUBTRACT) {
	  for(;index1<size1;index1++) {
		  arrResult->append((*arr1)[index1]);
		}
	 }
	 
	 if(oper == UNION) {
	  for(;index2<size2;index2++) {
		  arrResult->append((*arr2)[index2]);
	   }
	 }

	 (*arr1) = (*arrResult);
	 arr1->compress();
	 delete arrResult;
}

/* return the index of the variable assuming array is ordered in increasing order 
 * uses binary search*/
int orderedArrayFind(IntArray * const & arr, int elem) {
	 if(arr->size() == 0)
		  return -1;
	 int start = 0;
	 int end = arr->size()-1; 
	 int index;

	 if(elem < (*arr)[start] || elem > (*arr)[end])
		  return -1;
	 
	 while(start <= end) {
	        index = (start + end)/2;
		    if(elem == (*arr)[index])
				 return index;
			if(elem < (*arr)[index])
				 end = index-1;
			if(elem > (*arr)[index])
				 start = index+1;
	 }
	 return -1;
}

// check whether two ordered arrays are the same
bool orderedArrayAreSame(IntArray * const & arr1, IntArray * const & arr2) {
	 if(arr1->size() != arr2->size())
		  return false;

	 int size = arr1->size();
	 for(int i=0;i<size;i++) {
		  if((*arr1)[i] != (*arr2)[i])
			   return false;
	 }
	 return true;
}

//sort the index array based on the array of values given
void sortAscending(IntArray * const & indices, IntArray * const & vals) {
	 int size = vals->size();
	 indices->clear();
	 for(int i=0;i<size;i++) {
		  indices->append(i);
	 }
     for(int i=0;i<size;i++) {
      for(int j=i+1;j<size;j++) {
		   int indexi = (*indices)[i]; 
		   int indexj = (*indices)[j];
		   if((*vals)[indexi] > (*vals)[indexj]) {
				int tmp = (*indices)[i];
				(*indices)[i] = (*indices)[j];
				(*indices)[j] = tmp;
		   }
	  }
	 }
}


//Print the array
ostream& printArray(IntArray & array, ostream& out)
{
  return printArray(array, 0, out);
}

//Print the array
ostream& printArray(IntArray & array, int beginIndex, ostream& out)
{
  char delimiter = ' ';
  for (int i = beginIndex; i < array.size(); i++)
  {
    out<<array[i]<<delimiter;
  }
  return out;
}


/* Get the Maximum variable id appearing in the clause
 * Assumes there are no constants in the clause 
 * - should be in clause.h ideally */
int getMaxVarId(Clause * const & clause) {
     Predicate *pred;
	 const Term *term;
	 int maxVarId = 0;  
	 for(int predno=0;predno<clause->getNumPredicates();predno++) {
	     pred = new Predicate(*(clause->getPredicate(predno)));
         for (int termno = 0; termno < pred->getNumTerms(); termno++) {
            term = pred->getTerm(termno);
            int varId = -term->getId();
            assert(varId > 0);
			if(varId > maxVarId)
				 maxVarId = varId;
		  }
	 }
     return maxVarId;
}

//check whether any variable is repeated in the predicate's arguments
bool containsDuplicateVariable(Predicate * const & pred) {
	   IntHashArray *uniqueVarIds = new IntHashArray();
       int termCnt = pred->getNumTerms();
	   for (int termno = 0; termno < termCnt; termno++)
       {
	        const Term *term = pred->getTerm(termno);
			if(term->isConstant())
				 continue;
			int varId = term->getId();
			if(uniqueVarIds->find(varId) >= 0) {
				 delete uniqueVarIds;
				 return true;
			} else {
				 uniqueVarIds->append(varId);
			}
	   }
	   delete uniqueVarIds;
	   return false;
}


