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
/* This file contains methods for joining together the predicate hypercubes
 * to construct the corresponding clause hypercubes */

#include <ext/hash_set>
#include <ext/hash_map>

#include "util.h"
#include "hypercube.h"
#include "hypercubereverseindex.h"
#include "partialgndhctuple.h"
#include "inferutil.h"
#include "superclause.h"


using namespace std;
using namespace __gnu_cxx;


//Get array of hypercubes (corresponding to pred) joining with the given partially ground tuple
Array<HyperCube *> * getJoinHyperCubes(PartialGndHCTuple * const & pgTuple,
									   HyperCubeReverseIndex * const & predHCReverseIndex,
									   Predicate * const & pred) 
{
	   IntToIntMap *joinVarIdToPredTermNo = new IntToIntMap();
	   IntToIntMap::iterator itr; 
	   
	   IntHashArray *joinHCIndices, *hcIndices;
       IntHashArray *seenVarIds = new IntHashArray();
        	   
       for(int termno=0; termno<pred->getNumTerms();termno++) {
         const Term * term = pred->getTerm(termno);
         int varId = -term->getId();
		 
		 //only the first occurrence of each varId should be considered
		 if(seenVarIds->find(varId) >= 0)
			  continue;

		 seenVarIds->append(varId);
         //join needs to be done on the variables which have been gnded
		 if(!pgTuple->isVarGnded(varId))
			  continue;
		 (*joinVarIdToPredTermNo)[varId] = termno+1;
	   }
	   
	   delete seenVarIds;

	   //if there is no variable on which it is to be joined, simply return
	   //all the hypercubes
	   if(joinVarIdToPredTermNo->empty()) {
			delete joinVarIdToPredTermNo;
			return predHCReverseIndex->getAllHyperCubes();
	   }
       
	   for(itr = joinVarIdToPredTermNo->begin();
		   itr != joinVarIdToPredTermNo->end();
		   itr++) {
			
			int joinVarId = itr->first;
			int predTermNo = itr->second;
			int refJoinVarId = pgTuple->getReferenceVarId(joinVarId);
			IntArray *varConstants = pgTuple->getVarConstants(refJoinVarId);
			
	
			/*
			cout<<"varId = "<<joinVarId<<", termno = "<<predTermNo<<endl;
            cout<<"Var constants =>"<<endl;
			printArray(*varConstants,cout);
			cout<<endl;
			*/

			//retrieve the hypercube indices which correspond to the given set of
			//var constants
			int refPredTermNo = predHCReverseIndex->getReferenceVarId(predTermNo);
			hcIndices = predHCReverseIndex->getHyperCubeIndices(varConstants,refPredTermNo);
            //cout<<"Size of hcIndices "<<hcIndices->size()<<endl;

			if(itr == joinVarIdToPredTermNo->begin()) {
			     joinHCIndices = hcIndices;
			} else {
			     //intersectArray(joinHCIndices,hcIndices);
			     //while(true) {
				  //cout<<"LALAL LALALALA "<<endl;
				  intersectHashArray(joinHCIndices,hcIndices);
				 //}
			     delete hcIndices;
			}
	   }

	   Array<HyperCube *> * joinHyperCubes = predHCReverseIndex->getHyperCubes(joinHCIndices);
	   delete joinHCIndices;
	   delete joinVarIdToPredTermNo;
	   
	   return joinHyperCubes;
}

//get the clause hypercubes given:
//(a) the clause 
//(b) hypercubes for the constituent predicates
Array<HyperCube *> *getClauseHyperCubes(Clause * const & clause, 
		                                SuperClause * const & superClause,
		                                Array<HyperCubeReverseIndex *> * const & predHCReverseIndexArr,
		                                Array<HyperCubeRefinement *> * const & hyperCubeRefinementArr,
										IntHashArray * const & unknownPredIds,
										IntArrayHashArray * const & neqConstraints,
										bool useCT,
										Domain * const & domain) {

	 Array<HyperCube *> * clauseHyperCubes = NULL;
	 
	 Array<HyperCube *> * joinHyperCubes, *hyperCubes;
	 PartialGndHCTuple *pgTuple, *joinedPGTuple;
	 HyperCubeReverseIndex *predHCReverseIndex;
     HyperCube *hyperCube, *joinHyperCube;
	 Predicate * pred;

	 Array<int> *varIdToCanonicalVarId;
	 Array<PartialGndHCTuple *> *pgTuples = new Array<PartialGndHCTuple *>();
	 const Array<Predicate *> *preds = clause->getPredicates();

	 pgTuple = new PartialGndHCTuple(clause, useCT);
	 pgTuples->append(pgTuple);

	 //we return the array of hypercubes if superclause argument is null
	 if(!superClause)
	  clauseHyperCubes = new Array<HyperCube *>();
	
	 int tupleCnt = 0;
	 while(pgTuples->size() > 0) {
         pgTuple = pgTuples->removeLastItem();
		 int pindex = pgTuple->getFirstUnGndPredIndex();
         assert(pindex >= 0);
		 pred = (*preds)[pindex];
		 //cout<<"Dealing with the predicate : "<<endl;
		 //pred->print(cout,domain);
		 predHCReverseIndex = (*predHCReverseIndexArr)[pindex];
		 //cout<<"Number of hypercubes in the join index = ";
		 //cout<<predHCReverseIndex->getNumHyperCubes()<<endl;
		 joinHyperCubes = getJoinHyperCubes(pgTuple,predHCReverseIndex,pred);
		 //cout<<"Size of join hypercubes = "<<joinHyperCubes->size()<<endl;
         for(int i=0;i<joinHyperCubes->size();i++) {
			  joinHyperCube = (*joinHyperCubes)[i];
			  
			  /*
			  cout<<"Join hypercube is ="<<endl;
			  joinHyperCube->print(cout);
			  cout<<endl;
			  */

			  joinedPGTuple = new PartialGndHCTuple(pgTuple);
			  joinedPGTuple->join(joinHyperCube, pindex);
			  if(joinedPGTuple->isGnded()) {
				 joinedPGTuple->createConstrainedHyperCubes(neqConstraints);
				 joinedPGTuple->updateHyperCubeRefinement(hyperCubeRefinementArr,unknownPredIds);
                 hyperCubes = joinedPGTuple->getConstrainedHyperCubes();

				 //hyperCube = new HyperCube(joinedPGTuple->getHyperCube());
                 //hyperCube = joinedPGTuple->getHyperCube();
				 
				 for(int hindex=0;hindex<hyperCubes->size();hindex++) {
					hyperCube = (*hyperCubes)[hindex];
				    if(superClause == NULL) {
				     clauseHyperCubes->append(hyperCube);
				    } else {
					 varIdToCanonicalVarId = superClause->getVarIdToCanonicalVarId();
		
					 /*
					 cout<<"Before canonicalizing.. "<<endl;
	                 hyperCube->print(cout);
					 cout<<endl;
					 */
					 
					 int cnt = hyperCube->canonicalize(varIdToCanonicalVarId);
					 /*
					 cout<<"After canonicalizing.. "<<endl;
	                 hyperCube->print(cout);
					 cout<<endl;
					 */
					 //double wt = clause->getWt();
					 //superClause->addNewHyperCubeAndIncrementCount((*hyperCubes)[i],wt);
					 superClause->addNewHyperCubeAndIncrementCount(hyperCube,cnt);
					 tupleCnt++;
				   }
				 } 
			    delete joinedPGTuple;
			  } else {
				 pgTuples->append(joinedPGTuple);
			  }
		 }

		 delete joinHyperCubes;
		 delete pgTuple;
	 }

     delete pgTuples;
	 //cout<<"Number of hypercubes using Raw Calculuations = "<<tupleCnt<<endl;
	 return clauseHyperCubes;
}

