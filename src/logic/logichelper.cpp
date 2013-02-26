#include <assert.h>
#include <math.h>

#include "logichelper.h"

/*
#include <climits>
#include <stdlib.h>
#include <ostream>
using namespace std;
*/

// The discrete clauses defined in this library are of the same clause format in alchemy.

void GetClauseRandomSolution(int clause_size, bool satisfied, bool conjunction_disjunction, Array<bool>* solution)
{
	solution->clear();
	if (satisfied && conjunction_disjunction)
	{
		solution->growToSize(clause_size, true); // A conjunctive clause could only be satisfied when every variable is true.
	}
	else if(!satisfied && !conjunction_disjunction)
	{
		solution->growToSize(clause_size, false); // A disjunctive clause could only be unsatisfied when every variable is false.
	}
	else // For cases where conjunctive clause not to be satisfied and disjunctive clause to be satisfied.
	{
		bool good_solution = false;
		while (!good_solution)
		{
			solution->clear();
			for(int i = 0; i < clause_size; i++)
			{
				bool b = random()%2;
				if (((b && satisfied) || (!b && !satisfied)))
				{
					good_solution = true;
				}
				solution->append(b);
			}				
		}
	}
}

void RecordAtomMakeClauseFalse(const Array<int>& clause, const Array<bool>& atoms, bool conj_disj, Array<int>* record) {
        assert(record);
	if (conj_disj) {  // If conjunction, any false atom will do.
		// If there already exists false atom(s), no flip is necessary, otherwise, randomly pick an atom to flip.
		for (int i = 0; i < clause.size(); ++i) {
			if (atoms[abs(clause[i])] != (clause[i] > 0))  // False atom
			  return;
		}
		record->append(abs(clause[random() % clause.size()]));
	} else {  // If disjunction, every atom has to be false.
		for (int i = 0; i < clause.size(); ++i) {
			if(atoms[abs(clause[i])] == (clause[i] > 0)) {  // It's true literal, record its index as it is to be flipped.
				record->append(abs(clause[i]));
			}
		}
	}
}

void RecordAtomMakeClauseTrue(const Array<int>& clause, const Array<bool>& atoms, bool conj_disj, Array<int>* record) {
        assert(record);
	if (!conj_disj) {  // If disjunction, any true atom will do.
		// If there already exists true atom(s), no flip is necessary, otherwise, randomly pick an atom to flip.
		for (int i = 0; i < clause.size(); ++i) {
			if (atoms[abs(clause[i])] == (clause[i] > 0))  // True atom
			   return;
		}
		record->append(abs(clause[random() % clause.size()]));
	} else {  // If conjunction, every atom has to be true.
		for (int i = 0; i < clause.size(); ++i) {
			if(atoms[abs(clause[i])] != (clause[i] > 0)) {  // It's an false literal, record its index as it is to be flipped.
				record->append(abs(clause[i]));
			}
		}
	}
}
