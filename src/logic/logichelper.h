#ifndef LOGIC_HELPER_H_NOV_7_2008
#define LOGIC_HELPER_H_NOV_7_2008


#include "array.h"

// This function generate random solutionution to a propositional clause.
// The clause is specificed by two parameters: 1. Number of variables in the clause -- clauseSize; 2. It's a conjunction or disjunction clause, bConjOrDisj.
// Two possible state could be set for the clause: satisfied or not, by bSat.
// Result values are stored in solution.
void GetClauseRandomSolution(int clause_size, bool satisfied, bool conjunction_disjunction, Array<bool>* solution);

void RecordAtomMakeClauseFalse(const Array<int>& clause, const Array<bool>& atoms, bool conj_disj, Array<int>* record);
void RecordAtomMakeClauseTrue(const Array<int>& clause, const Array<bool>& atoms, bool conj_disj, Array<int>* record);




#endif
