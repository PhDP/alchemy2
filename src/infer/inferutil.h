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
#ifndef _INFER_UTIL_H_MAY_2008
#define _INFER_UTIL_H_MAY_2008

#include <stdlib.h>
#include <assert.h>
#include <iostream>

#include "array.h"
#include "hashint.h"
#include "clause.h"

#include <ext/hash_map>
#include <ext/hash_set>


using namespace __gnu_cxx;
using namespace std;

const int MAX_INT = 2147483647;

/***************************************************************************/
// Typedefs
/***************************************************************************/

typedef hash_map<int, int, HashInt, EqualInt> IntToIntMap;
typedef Array<int> IntArray;
typedef hash_map<int,IntHashArray*,HashInt,EqualInt> IntToIntHashArrayMap;

enum HyperCubeCreateType {
     Basic,
	 DT
};


/***************************************************************************/
// Mehods
/***************************************************************************/
 
/* Ideally, should be in array.h */


/*intersect the given int hash arrays and store the result in the
 * first array */ 
void intersectHashArray(IntHashArray * const & arr1, IntHashArray * const & arr2); 
void intersectHashArrayKeepOrder(IntHashArray * const & arr1, IntHashArray * const & arr2); 


//get the union of two Int Arrays into first Int Array
void unionArrayKeepOrder(IntArray * const & arr1, IntArray * const & arr2);

//get the union of two Int Arrays into first Int Array
void intersectArrayKeepOrder(IntArray * const & arr1, IntArray * const & arr2); 

// subtract the second int array from first int array and store the result in first one 
void subtractArrayKeepOrder(IntArray * const & arr1, IntArray * const & arr2); 

// return the index of the variable assuming array is ordered in increasing order
int orderedArrayFind(IntArray * const & arr, int elem);

// check whether two ordered arrays are the same
bool orderedArrayAreSame(IntArray * const & arr1, IntArray * const & arr2);

//sort the index array based on the array of values given
void sortAscending(IntArray * const & indices, IntArray * const & vals);

//Print the array
ostream& printArray(IntArray & array, ostream& out);

//Print the array
ostream& printArray(IntArray & array, int beginIndex, ostream& out);


int getMaxVarId(Clause * const & clause);

//check whether any variable is repeated in the predicate's arguments
bool containsDuplicateVariable(Predicate * const & pred);


#endif

