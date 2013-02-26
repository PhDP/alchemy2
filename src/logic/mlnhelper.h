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
#ifndef MLNHELPER_H_NOV_4_2005
#define MLNHELPER_H_NOV_4_2005

//auxiliary data structures used by MLN


  //index is the index of an IndexClause object in an array.
struct IndexClause 
{
  IndexClause(const int& iindex, Clause* const & cclause)
    : index(iindex), clause(cclause) {}
  ~IndexClause() {}
  int index; 
  Clause* clause; 
};

class HashIndexClause
{
 public:
  size_t operator()(IndexClause* const& ic) const 
  { return ic->clause->hashCode(); }
};

class EqualIndexClause
{
 public:
  bool operator()(IndexClause* const & ic1, IndexClause* const & ic2) const
  { return ic1->clause->same(ic2->clause); }
};

typedef HashArray<IndexClause*, HashIndexClause, EqualIndexClause> 
  IndexClauseHashArray;

/////////////////////////////////////////////////////////////////////////////


  //Used to quickly access a clause in predIdToClausesMap_. predId is an index
  //of predIdToClausesMap_, while clauseIndex is an index of the clause array 
  //predIdToClausesMap_[predId];
struct PredIdClauseIndex 
{
  PredIdClauseIndex(const int& ppredId, int* const & cclauseIndex)
    : predId(ppredId), clauseIndex(cclauseIndex) {}
  ~PredIdClauseIndex() {}
  int predId; 
  int* clauseIndex; 
};

  //Used to quickly access a clause in formAndClauses_. formulaIndex an index of
  //formAndClauses_ and clauseIndex is an index of forAndClauses_[formulaIndex].
struct FormulaClauseIndexes 
{ 
  FormulaClauseIndexes(int* const & fformulaIndex, int* const & cclauseIndex)
    : formulaIndex(fformulaIndex), clauseIndex(cclauseIndex) {}
  ~FormulaClauseIndexes() {}
  int* formulaIndex; 
  int* clauseIndex; 
};

  //Each MLNClauseInfo corresponds to a clause. index is the index of 
  //MLNClauseInfo in clauseInfos_. It should be the same as the index of its 
  //corresponding clause in clauses_. 
  //predIdsClauseIndexes and formulaClauseIndexes help one to quickly find the 
  //occurrences of the clause in predIdToClausesMap_, and formAndClauses_.
struct MLNClauseInfo
{
  MLNClauseInfo(const int& iindex, const double& ppriorMean) 
    : index(iindex), priorMean(ppriorMean) {}
  MLNClauseInfo(const int& iindex) : index(iindex), priorMean(0) {}

  ~MLNClauseInfo() { predIdsClauseIndexes.deleteItemsAndClear();
                     formulaClauseIndexes.deleteItemsAndClear(); }
  void compress() { predIdsClauseIndexes.compress(); 
                    formulaClauseIndexes.compress(); }
  int index;
  Array<PredIdClauseIndex*> predIdsClauseIndexes;
  Array<FormulaClauseIndexes*> formulaClauseIndexes;
  double priorMean;
};

/////////////////////////////////////////////////////////////////////////////

struct FormulaAndClauses 
{ 
  FormulaAndClauses(const string& fformula, const int& iindex, 
                    const bool& hhasExist, const bool& ttiedClauses,
                    const bool& iisConjunction)
    : formula(fformula), indexClauses(new IndexClauseHashArray), index(iindex), 
      hasExist(hhasExist), tiedClauses(ttiedClauses),
      isConjunction(iisConjunction), numPreds(-1), isHard(false), priorMean(0),
      wt(0), isExistUnique(false) {}
  ~FormulaAndClauses() { indexClauses->deleteItemsAndClear();
                         delete indexClauses; }
  string formula; 
  IndexClauseHashArray* indexClauses; 
  int index; 
  bool hasExist;
  bool tiedClauses;
  bool isConjunction;
  int numPreds;
  bool isHard;
  double priorMean;
  double wt;
  bool isExistUnique; //contains existentially and uniquely quant. vars
};

class HashFormulaAndClauses
{
 public:
  size_t operator()(const FormulaAndClauses* const& f) const
  { return hash<char const *>()(f->formula.c_str()); }
};

class EqualHashFormulaAndClauses
{
 public:
  bool operator()(const FormulaAndClauses* const & f1, 
                  const FormulaAndClauses* const & f2) const
  { return (f1->formula.compare(f2->formula)==0); }
};

typedef HashArray<FormulaAndClauses*, HashFormulaAndClauses, 
  EqualHashFormulaAndClauses> FormulaAndClausesArray;


#endif
