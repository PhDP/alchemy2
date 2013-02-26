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
#ifndef VARIABLESTATE_H_
#define VARIABLESTATE_H_

#include "mrf.h"
#include "timer.h"

const int NOVALUE = 100000000;
const int DISANDCONT = 200000000;
const int MULTIDIS = 300000000;
const double NOSOL = 1234567890;
const bool vsdebug = false;

/**
 * Represents the state of propositional variables and clauses. Some of this
 * code is based on the MaxWalkSat package of Kautz et al.
 * 
 * All inference algorithms should have a VariableState to access the
 * information needed in its predicates and clauses.
 * 
 * Each atom has its own index starting at 1. The negation of an atom with
 * index a is represented by -a (this is why the indices do not start at 0).
 * Each clause has its own index starting at 0.
 * 
 * A VariableState is either eager or lazy. Eager states build an MRF upfront
 * based on an MLN and a domain. Thus, all ground clauses and predicates are
 * in memory after building the MRF. Lazy states activate atoms and clauses as
 * they are needed from the MLN and domain. An atom is activated if it is in
 * an unsatisfied clause with the assumption of all atoms being false or if it
 * is looked at during inference (it is flipped or the cost of flipping it is
 * computed). Active clauses are those which contain active atoms.
 */
class VariableState
{
 public:

  /**
   * Constructor for a VariableState. The hard clause weight is set. In lazy
   * mode, the initial active atoms and clauses are retrieved and in eager
   * mode, the MRF is built and the atoms and clauses are retrieved from it.
   * In addition, all information arrays are filled with information.
   * 
   * @param unknownQueries Query predicates with unknown values used to build
   * MRF in eager mode.
   * @param knownQueries Query predicates with known values used to build MRF
   * in eager mode.
   * @param knownQueryValues Truth values of the known query predicates.
   * @param allPredGndingsAreQueries Array used to build MRF in eager mode.
   * @param mln mln and domain are used to build MRF in eager state and to
   * retrieve active atoms in lazy state.
   * @param domain mln and domain are used to build MRF in eager state and to
   * retrieve active atoms in lazy state.
   * @param lazy Flag stating whether lazy mode is used or not.
   */
  VariableState(GroundPredicateHashArray* const& unknownQueries,
                GroundPredicateHashArray* const& knownQueries,
                Array<TruthValue>* const & knownQueryValues,
                const Array<int>* const & allPredGndingsAreQueries,
                const bool& markHardGndClauses,
                const bool& trackParentClauseWts,
                const MLN* const & mln, const Domain* const & domain,
                const bool& lazy)
  {
    stillActivating_ = true;
    breakHardClauses_ = false;
      // By default MaxWalkSAT mode
    inferenceMode_ = MODE_MWS;
    Timer timer;
    double startTime = timer.time();

    this->mln_ = (MLN*)mln;
    this->domain_ = (Domain*)domain;
    this->lazy_ = lazy;

      // Instantiate information
    baseNumAtoms_ = 0;
    activeAtoms_ = 0;
    numFalseClauses_ = 0;
    costOfFalseClauses_ = 0.0;
    lowCost_ = LDBL_MAX;
    lowBad_ = INT_MAX;

      // Clauses and preds are stored in gndClauses_ and gndPreds_
    gndClauses_ = new GroundClauseHashArray;
    gndPreds_ = new Array<GroundPredicate*>;

      // Set the hard clause weight
    setHardClauseWeight();

      // Lazy version: Produce state with initial active atoms and clauses
    if (lazy_)
    {
        // Set number of non-evidence atoms from domain
      domain_->computeNumNonEvidAtoms();
      numNonEvAtoms_ = domain_->getNumNonEvidenceAtoms();
        // Unknown preds are treated as false
      domain_->getDB()->setPerformingInference(true);
      clauseLimit_ = INT_MAX;
      noApprox_ = false;
      haveDeactivated_ = false;

        ///// Get initial active atoms and clauses /////
        // Get initial set of active atoms (atoms in unsat. clauses)
        // Assumption is: all atoms are initially false except those in blocks

        // One atom in each block is set to true and activated
      initBlocksRandom();

      //bool ignoreActivePreds = false;
      bool ignoreActivePreds = true;
      cout << "Getting initial active atoms ... " << endl;
      getActiveClauses(newClauses_, ignoreActivePreds);
      cout << "done." << endl;
      int defaultCnt = newClauses_.size();
      long double defaultCost = 0;

      for (int i = 0; i < defaultCnt; i++)
      {
        if (newClauses_[i]->isHardClause())
          defaultCost += hardWt_;
        else
          defaultCost += abs(newClauses_[i]->getWt());
      }

        // Clear ground clauses in the ground preds
      for (int i = 0; i < gndPredHashArray_.size(); i++)
        gndPredHashArray_[i]->removeGndClauses();

        // Delete new clauses
      for (int i = 0; i < newClauses_.size(); i++)
        delete newClauses_[i];
      newClauses_.clear();

      baseNumAtoms_ = gndPredHashArray_.size();
      cout << "Number of Baseatoms = " << baseNumAtoms_ << endl;
      cout << "Default => Cost\t" << "******\t" << " Clause Cnt\t" << endl;
      cout << "           " << defaultCost << "\t" << "******\t" << defaultCnt
           << "\t" << endl << endl;

        // Set base atoms as active in DB
      for (int i = 0; i < baseNumAtoms_; i++)
      {
        domain_->getDB()->setActiveStatus(gndPredHashArray_[i], true);
        activeAtoms_++;     
      }

        // Get the initial set of active clauses
      ignoreActivePreds = false;
      cout << "Getting initial active clauses ... ";
      getActiveClauses(newClauses_, ignoreActivePreds);      
      cout << "done." << endl;
	} // End lazy version
      // Eager version: Use KBMC to produce the state
    else
    {
      unePreds_ = unknownQueries;
      knePreds_ = knownQueries;
      knePredValues_ = knownQueryValues;

        // MRF is built on known and unknown queries
      int size = 0;
      if (unknownQueries) size += unknownQueries->size();
      if (knownQueries) size += knownQueries->size();
      GroundPredicateHashArray* queries = new GroundPredicateHashArray(size);
      if (unknownQueries) queries->append(unknownQueries);
      if (knownQueries) queries->append(knownQueries);
      mrf_ = new MRF(queries, allPredGndingsAreQueries, domain_,
                     domain_->getDB(), mln_, markHardGndClauses,
                     trackParentClauseWts, -1);

        //delete to save space. Can be deleted because no more gndClauses are
        //appended to gndPreds beyond this point
      mrf_->deleteGndPredsGndClauseSets();
        //do not delete the intArrRep in gndPreds_;
      delete queries;

        // Put ground clauses in newClauses_
      newClauses_ = *(Array<GroundClause*>*)mrf_->getGndClauses();
        // Put ground preds in the hash array
      //const Array<GroundPredicate*>* gndPreds = mrf_->getGndPreds();
      const GroundPredicateHashArray* gndPreds = mrf_->getGndPreds();
      for (int i = 0; i < gndPreds->size(); i++)
        gndPredHashArray_.append((*gndPreds)[i]);
    
        // baseNumAtoms_ are all atoms in eager version
      baseNumAtoms_ = gndPredHashArray_.size();        
    } // End eager version
    
      // At this point, ground clauses are held in newClauses_
      // and ground predicates are held in gndPredHashArray_
      // for both versions
    
      // Add the clauses and preds and fill info arrays
    bool initial = true;
    addNewClauses(initial);
    
	cout << "[VS] ";
	Timer::printTime(cout,timer.time()-startTime);
	cout << endl;
	cout << ">>> DONE: Initial num. of clauses: " << getNumClauses() << endl;
  }

  /**
   * Destructor. MRF is deleted in eager version.
   */ 
  ~VariableState()
  {
    if (lazy_)
    {
      if (gndClauses_)
        for (int i = 0; i < gndClauses_->size(); i++)
          delete (*gndClauses_)[i];

      for (int i = 0; i < gndPredHashArray_.size(); i++)
      {
        gndPredHashArray_[i]->removeGndClauses();
        delete gndPredHashArray_[i];
      }
    }
    else
    {
        // MRF from eager version is deleted
      if (mrf_)
      {
        delete mrf_;
        mrf_ = NULL;
      }
      //if (unePreds_) delete unePreds_;
      //if (knePreds_) delete knePreds_;
      //if (knePredValues_) delete knePredValues_;
    }
    //delete gndClauses_;
    //delete gndPreds_;
    //delete mln_;
  }


  /**
   * New clauses are added to the state. If not the initialization, then
   * makecost and breakcost are updated for the new atoms.
   * 
   * @param initial If true, this is the first time new clauses have been
   * added.
   */
  void addNewClauses(bool initial)
  {
    if (initial) addNewClauses(ADD_CLAUSE_INITIAL, newClauses_);
    else addNewClauses(ADD_CLAUSE_REGULAR, newClauses_);
  }

  /**
   * Information about the state is reset and initialized.
   */
  void init()
  {
      // Initialize info arrays
    initMakeBreakCostWatch();
  }

  /**
   * State is re-initialized with all new clauses and atoms.
   */
/*
  void reinit()
  {
    clause_.clearAndCompress();
    clauseCost_.clearAndCompress();
    falseClause_.clearAndCompress();
    whereFalse_.clearAndCompress();
    numTrueLits_.clearAndCompress();
    watch1_.clearAndCompress();
    watch2_.clearAndCompress();
    isSatisfied_.clearAndCompress();
    deadClause_.clearAndCompress();
    threshold_.clearAndCompress();

    //newClauses_.append(gndClauses_);
    for (int i = 0; i < gndClauses_->size(); i++)
      newClauses_.append((*gndClauses_)[i]);

    gndClauses_->clearAndCompress();
    gndPreds_->clearAndCompress();
    for (int i = 0; i < occurence_.size(); i++)
      occurence_[i].clearAndCompress();
    occurence_.clearAndCompress();
    
      // Add the clauses and preds and fill info arrays
    bool initial = true;
    addNewClauses(initial);
    baseNumAtoms_ = gndPredHashArray_.size();
    init();    
  }
*/  
  /**
   * Makes a random truth assigment to all (active) atoms. Blocks are
   * taken into account: exactly one atom in the block is set to true
   * and the others are set to false.
   */
  void initRandom()
  {
      // If performing SampleSat, do a run of UP
  	if (inferenceMode_ == MODE_SAMPLESAT)
  	{
      unitPropagation();
    }

      // Set one in each block to true randomly
    initBlocksRandom();

    for (int i = 1; i <= gndPreds_->size(); i++)
    {
      if (fixedAtom_[i] != 0)
      {
        continue;
      }

        // Blocks are initialized above
      if (getBlockIndex(i - 1) >= 0)
      {
      	if (vsdebug) cout << "Atom " << i << " in block" << endl;
      	continue;
      }
        // Not fixed and not in block
      else
      {
        if (vsdebug) cout << "Atom " << i << " not in block" << endl;
			
          // only need to activate inactive true atoms
        if (!lazy_ || isActive(i))
        {
          bool isTrue = random() % 2;
          bool activate = false;
          setValueOfAtom(i, isTrue, activate, -1);
        }
        else if (inferenceMode_ == MODE_SAMPLESAT)
        {
          bool isInPrevUnsat = false;
          Array<int> oca = getPosOccurenceArray(i);
          for (int k = 0; k < oca.size(); k++)
            if (!prevSatisfiedClause_[oca[k]])
            {
              isInPrevUnsat = true;
              break;
            }
            
          if (!isInPrevUnsat)
          {
            oca = getNegOccurenceArray(i);
            for (int k = 0; k < oca.size(); k++)
              if (!prevSatisfiedClause_[oca[k]])
              {
                isInPrevUnsat = true;
                break;
              }
          }
          if (isInPrevUnsat)
          {
            bool isTrue = random() % 2;
            if (isTrue)
            {
              if (activateAtom(i, false, false))
                setValueOfAtom(i, true, false, -1);
            }
          }
        }
      }
    }

    init();	
  }

  /**
   * Sets one atom in each block to true and the rest to false.
   */
  void initBlocksRandom()
  {
    if (vsdebug)
    {
      cout << "Initializing blocks randomly" << endl;
      cout << "Num. of blocks: " << domain_->getNumPredBlocks() << endl;
    }
      // For each block: select one to set to true
    for (int i = 0; i < domain_->getNumPredBlocks(); i++)
    {
      int trueFixedAtom = getTrueFixedAtomInBlock(i);
        // True fixed atom in the block: set others to false
      if (trueFixedAtom >= 0)
      {
        if (vsdebug)
        {
          cout << "True fixed atom " << trueFixedAtom  << " in block "
               << i << endl;
        }
        setOthersInBlockToFalse(trueFixedAtom, i);
        continue;
      }

        // If evidence atom exists, then all others are false
      if (domain_->getBlockEvidence(i))
      {
        if (vsdebug) cout << "Block evidence in block " << i << endl;
          // If first argument is -1, then all are set to false
        setOthersInBlockToFalse(-1, i);
        continue;
      }

        // Pick one pred from block at random
      bool ok = false;
      while (!ok)
      {
        const Predicate* pred = domain_->getRandomPredInBlock(i);
        GroundPredicate* gndPred = new GroundPredicate((Predicate*)pred);
        int atomIdx = gndPredHashArray_.find(gndPred);
        delete pred;

        assert(lazy_ || atomIdx >= 0);        
          // Pred not yet present: add it
        if (atomIdx == -1)
        {
          atomIdx = gndPredHashArray_.append(gndPred);
          bool initial = false;
          addNewClauses(initial);
          ok = true;
        }
          // Pred already present: check that it's not fixed
        else
        {
          delete gndPred;
          if (fixedAtom_[atomIdx + 1] == 0)
          {
            if (vsdebug) cout << "Atom " << atomIdx + 1 
                              << " chosen in block " << i << endl;
            ok = true;
          }
          else
          {
            if (vsdebug) cout << "Atom " << atomIdx + 1 
                              << " is fixed to " << fixedAtom_[atomIdx + 1]
                              << endl;
              // Atom is fixed
            continue;
          }
        }
        bool activate = false;
        setValueOfAtom(atomIdx + 1, true, activate, i);
        setOthersInBlockToFalse(atomIdx, i);
      }
    }
    if (vsdebug) cout << "Done initializing blocks randomly" << endl;
  }      

  /**
   * Resets the makeCost, breakCost and watch arrays.
   */
  void resetMakeBreakCostWatch()
  {
  	  // Reset info concerning true lits, false clauses, etc.
  	for (int i = 0; i < getNumClauses(); i++) numTrueLits_[i] = 0;
    numFalseClauses_ = 0;
    costOfFalseClauses_ = 0.0;
    lowCost_ = LDBL_MAX;
    lowBad_ = INT_MAX;

    assert(makeCost_.size() == breakCost_.size());
      // Reset make and break cost to 0
    for (int i = 0; i < makeCost_.size(); i++)
    {
      makeCost_[i] = breakCost_[i] = 0.0;
    }
  }


  /**
   * Re-initializes the makeCost, breakCost and watch arrays based on the
   * current variable assignment.
   */
  void initMakeBreakCostWatch()
  {
    resetMakeBreakCostWatch();
    initMakeBreakCostWatch(0);
  }

  /**
   * Initialize makeCost and breakCost and watch arrays based on the current
   * variable assignment.
   * 
   * @param startClause All clauses with index of this or greater are
   * initialized.
   */
  void initMakeBreakCostWatch(const int& startClause)
  {
    int theTrueLit = -1;
      // Initialize breakCost and makeCost in the following:
    for (int i = startClause; i < getNumClauses(); i++)
    {
        // Don't look at dead or satisfied clauses
      if (deadClause_[i] || isSatisfied_[i]) continue;

      int trueLit1 = 0;
      int trueLit2 = 0;
      long double cost = clauseCost_[i];
      numTrueLits_[i] = 0;
      for (int j = 0; j < getClauseSize(i); j++)
      {
        if (isTrueLiteral(clause_[i][j]))
        { // ij is true lit
          numTrueLits_[i]++;
          theTrueLit = abs(clause_[i][j]);
          if (!trueLit1) trueLit1 = theTrueLit;
          else if (trueLit1 && !trueLit2) trueLit2 = theTrueLit;
        }
      }

        // Unsatisfied positive-weighted clauses or
        // Satisfied negative-weighted clauses
      if ((numTrueLits_[i] == 0 && cost >= 0) ||
          (numTrueLits_[i] > 0 && cost < 0))
      {
        whereFalse_[i] = numFalseClauses_;
        falseClause_[numFalseClauses_] = i;
        numFalseClauses_++;
        costOfFalseClauses_ += abs(cost);
        if (highestCost_ == abs(cost)) {eqHighest_ = true; numHighest_++;}

          // Unsat. pos. clause: increase makeCost_ of all atoms
        if (numTrueLits_[i] == 0)
          for (int j = 0; j < getClauseSize(i); j++)
          {
            makeCost_[abs(clause_[i][j])] += cost;
          }

          // Sat. neg. clause: increase makeCost_ if one true literal
        if (numTrueLits_[i] == 1)
        {
            // Subtract neg. cost
          makeCost_[theTrueLit] -= cost;
          watch1_[i] = theTrueLit;
        }
        else if (numTrueLits_[i] > 1)
        {
          watch1_[i] = trueLit1;
          watch2_[i] = trueLit2;
        }
      }
        // Pos. clauses satisfied by one true literal
      else if (numTrueLits_[i] == 1 && cost >= 0)
      {
        breakCost_[theTrueLit] += cost;
        watch1_[i] = theTrueLit;
      }
        // Pos. clauses satisfied by 2 or more true literals
      else if (cost >= 0)
      { /*if (numtruelit[i] == 2)*/
        watch1_[i] = trueLit1;
        watch2_[i] = trueLit2;
      }
        // Unsat. neg. clauses: increase breakCost of all atoms
      else if (numTrueLits_[i] == 0 && cost < 0)
      {
        for (int j = 0; j < getClauseSize(i); j++)
          breakCost_[abs(clause_[i][j])] -= cost;
      }
    } // for all clauses
  }

  int getNumAtoms() { return gndPreds_->size(); }
  
  int getNumClauses() { return gndClauses_->size(); }
  
  int getNumDeadClauses()
  { 
    int count = 0;
    for (int i = 0; i < deadClause_.size(); i++)
      if (deadClause_[i]) count++;
    return count;
  }

  /**
   * Returns the absolute index of a random atom. If in lazy mode, a
   * non-evidence atom which is not in memory could be picked. In this case,
   * the atom is added to the set in memory and its index is returned.
   */
  int getIndexOfRandomAtom()
  {
      // Lazy: pick a random non-evidence atom from domain
    if (lazy_)
    {
      Predicate* pred = domain_->getNonEvidenceAtom(random() % numNonEvAtoms_);
      GroundPredicate* gndPred = new GroundPredicate(pred);
      delete pred;

      int idx = gndPredHashArray_.find(gndPred);
        // Already present, then return index
      if (idx >= 0)
      {
        delete gndPred;
        return idx + 1;
      }
        // Not present, add it to the state
      else
      {
        if (vsdebug)
        {
          cout << "Adding randomly ";
          gndPred->print(cout, domain_);
          cout << " to the state" << endl; 
        }
        gndPredHashArray_.append(gndPred);
        bool initial = false;
        addNewClauses(initial);
          // index to return is the last element
        return gndPredHashArray_.size();
      }
    }
      // Eager: pick an atom from the state
    else
    {
      int numAtoms = getNumAtoms();
      if (numAtoms == 0) return NOVALUE;
      return random()%numAtoms + 1;
    }
  }

  /**
   * Returns the absolute index of a random atom in a random unsatisfied
   * pos. clause or the absolute index of a random true literal in a random
   * satisfied clause.
   */
  int getIndexOfAtomInRandomFalseClause()
  {
    if (numFalseClauses_ == 0) return NOVALUE;
    int clauseIdx = falseClause_[random()%numFalseClauses_];

      // Pos. clause: return index of random atom
    if (clauseCost_[clauseIdx] >= 0)
    {
        // Q: Can all atoms in a clause be fixed?
      while (true)
      {
        int i = random()%getClauseSize(clauseIdx);
        if (!fixedAtom_[abs(clause_[clauseIdx][i])])
          return abs(clause_[clauseIdx][i]);
      }
    }
      // Neg. clause: find random true lit
    else
      return getRandomTrueLitInClause(clauseIdx);
  }
  
  /**
   * Returns the index of a random unsatisfied pos. clause or a random
   * satisfied neg. clause.
   */
  int getRandomFalseClauseIndex()
  {
    if (numFalseClauses_ == 0) return NOVALUE;
    return falseClause_[random()%numFalseClauses_];
  }
  
  /**
   * Returns the cost of the unsatisfied positive and satisfied negative
   * clauses.
   */
  long double getCostOfFalseClauses()
  {
    return costOfFalseClauses_;
  }
  
  /**
   * Returns the number of the unsatisfied positive and
   * satisfied negative clauses.
   */
  int getNumFalseClauses()
  {
    return numFalseClauses_;
  }

  /**
   * Returns the truth value of an atom.
   * 
   * @atomIdx Index of atom whose truth value is returned.
   * @return Truth value of atom.
   */
  bool getValueOfAtom(const int& atomIdx)
  {
    return atom_[atomIdx];
  }

  /**
   * Returns the truth value of an atom in the best state.
   * 
   * @atomIdx Index of atom whose truth value is returned.
   * @return Truth value of atom.
   */
  bool getValueOfLowAtom(const int& atomIdx)
  {
    return lowAtom_[atomIdx];
  }

  /**
   * Sets the truth value of an atom. The truth value is propagated to the
   * database.
   * 
   * @param atomIdx Index of atom whose value is to be set.
   * @param value Truth value being set.
   * @param activate If true, atom is activated if not active.
   * @param blockIdx Indicates which block the atom being flipped is in. If not
   * in one, this is -1.
   */
  void setValueOfAtom(const int& atomIdx, const bool& value,
                      const bool& activate, const int& blockIdx)
  {
      // If atom already has this value, then do nothing
    if (atom_[atomIdx] == value) return;
    if (vsdebug) cout << "Setting value of atom " << atomIdx 
                      << " to " << value << endl;
      // Propagate assignment to DB
    GroundPredicate* p = gndPredHashArray_[atomIdx - 1];
    if (value)
      domain_->getDB()->setValue(p, TRUE);
    else
      domain_->getDB()->setValue(p, FALSE);

      // If not active, then activate it
    if (activate && lazy_ && !isActive(atomIdx))
    {
      bool ignoreActivePreds = false;
      bool groundOnly = false;
      activateAtom(atomIdx, ignoreActivePreds, groundOnly);
    }
    atom_[atomIdx] = value;

      // If in block and setting to true, tell the domain
    if (blockIdx > -1 && value)
    {
      Predicate* pred = p->createEquivalentPredicate(domain_);
      domain_->setTruePredInBlock(blockIdx, pred);
      if (vsdebug)
      {
        cout << "Set true pred in block " << blockIdx << " to ";
        pred->printWithStrVar(cout, domain_);
        cout << endl;
      }
    }
  }

  /**
   * Retrieves the negative occurence array of an atom.
   */
  Array<int>& getNegOccurenceArray(const int& atomIdx)
  {
    int litIdx = 2*atomIdx;
    return getOccurenceArray(litIdx);
  }

  /**
   * Retrieves the positive occurence array of an atom.
   */
  Array<int>& getPosOccurenceArray(const int& atomIdx)
  {
    int litIdx = 2*atomIdx - 1;
    return getOccurenceArray(litIdx);
  }

  /**
   * Flip the truth value of an atom and update info arrays.
   * 
   * @param toFlip Index of atom to be flipped.
   * @param blockIdx Index of block in which atom is, or -1 if not in one.
   */
  void flipAtom(const int& toFlip, const int& blockIdx)
  {
    bool toFlipValue = getValueOfAtom(toFlip);
    register int clauseIdx;
    int sign;
    int oppSign;
    int litIdx;
    if (toFlipValue)
      sign = 1;
    else
      sign = 0;
    oppSign = sign ^ 1;

    flipAtomValue(toFlip, blockIdx);
      // Update all clauses in which the atom occurs as a true literal
    litIdx = 2*toFlip - sign;
    Array<int>& posOccArray = getOccurenceArray(litIdx);
    for (int i = 0; i < posOccArray.size(); i++)
    {
      clauseIdx = posOccArray[i];
        // Don't look at dead or satisfied clauses
      if (deadClause_[clauseIdx] || isSatisfied_[clauseIdx]) continue;

	  // The true lit became a false lit
	  int numTrueLits = decrementNumTrueLits(clauseIdx);
      long double cost = getClauseCost(clauseIdx);
      int watch1 = getWatch1(clauseIdx);
      int watch2 = getWatch2(clauseIdx);
        // 1. If last true lit was flipped, then we have to update
        // the makecost / breakcost info accordingly        
      if (numTrueLits == 0)
      {
          // Pos. clause
        if (cost >= 0)
        {
            // Add this clause as false in the state
          addFalseClause(clauseIdx);
            // Decrease toFlip's breakcost (add neg. cost)
          addBreakCost(toFlip, -cost);
            // Increase makecost of all vars in clause (add pos. cost)
          addMakeCostToAtomsInClause(clauseIdx, cost);
        }
          // Neg. clause
        else
        {
          assert(cost < 0);
            // Remove this clause as false in the state
          removeFalseClause(clauseIdx);
            // Increase breakcost of all vars in clause (add pos. cost)
          addBreakCostToAtomsInClause(clauseIdx, -cost);        
            // Decrease toFlip's makecost (add neg. cost)
          addMakeCost(toFlip, cost);
        }
      }
        // 2. If there is now one true lit left, then move watch2
        // up to watch1 and increase the breakcost / makecost of watch1
      else if (numTrueLits == 1)
      {
        if (watch1 == toFlip)
        {
          assert(watch1 != watch2);
          setWatch1(clauseIdx, watch2);
          watch1 = getWatch1(clauseIdx);
        }

          // Pos. clause: Increase toFlip's breakcost (add pos. cost)
        if (cost >= 0)
        {
          addBreakCost(watch1, cost);
        }
          // Neg. clause: Increase toFlip's makecost (add pos. cost)
        else
        {
          assert(cost < 0);
          addMakeCost(watch1, -cost);
        }
      }
        // 3. If there are 2 or more true lits left, then we have to
        // find a new true lit to watch if one was flipped
      else
      { /* numtruelit[clauseIdx] >= 2 */
          // If watch1[clauseIdx] has been flipped
        if (watch1 == toFlip)
        {
            // find a different true literal to watch
          int diffTrueLit = getTrueLiteralOtherThan(clauseIdx, watch1, watch2);
          setWatch1(clauseIdx, diffTrueLit);
        }
          // If watch2[clauseIdx] has been flipped
        else if (watch2 == toFlip)
        {
            // find a different true literal to watch
          int diffTrueLit = getTrueLiteralOtherThan(clauseIdx, watch1, watch2);
          setWatch2(clauseIdx, diffTrueLit);
        }
      }
    }
        
      // Update all clauses in which the atom occurs as a false literal
    litIdx = 2*toFlip - oppSign;
    Array<int>& negOccArray = getOccurenceArray(litIdx);
    for (int i = 0; i < negOccArray.size(); i++)
    {
      clauseIdx = negOccArray[i];
        // Don't look at dead or satisfied clauses
      if (deadClause_[clauseIdx] || isSatisfied_[clauseIdx]) continue;

        // The false lit became a true lit
      int numTrueLits = incrementNumTrueLits(clauseIdx);
      long double cost = getClauseCost(clauseIdx);
      int watch1 = getWatch1(clauseIdx);

        // 1. If this is the only true lit, then we have to update
        // the makecost / breakcost info accordingly        
      if (numTrueLits == 1)
      {
          // Pos. clause
        if (cost >= 0)
        {
            // Remove this clause as false in the state
          removeFalseClause(clauseIdx);
            // Increase toFlip's breakcost (add pos. cost)
          addBreakCost(toFlip, cost);        
            // Decrease makecost of all vars in clause (add neg. cost)
          addMakeCostToAtomsInClause(clauseIdx, -cost);
        }
          // Neg. clause
        else
        {
          assert(cost < 0);
            // Add this clause as false in the state
          addFalseClause(clauseIdx);
            // Decrease breakcost of all vars in clause (add neg. cost)
          addBreakCostToAtomsInClause(clauseIdx, cost);
            // Increase toFlip's makecost (add pos. cost)
          addMakeCost(toFlip, -cost);
        }
          // Watch this atom
        setWatch1(clauseIdx, toFlip);
      }
        // 2. If there are now exactly 2 true lits, then watch second atom
        // and update breakcost
      else
      if (numTrueLits == 2)
      {
        if (cost >= 0)
        {
            // Pos. clause
            // Decrease breakcost of first atom being watched (add neg. cost)
          addBreakCost(watch1, -cost);
        }
        else
        {
            // Neg. clause
          assert(cost < 0);
            // Decrease makecost of first atom being watched (add neg. cost)
          addMakeCost(watch1, cost);
        }
        
          // Watch second atom
        setWatch2(clauseIdx, toFlip);
      }
    }
  }

  /**
   * Flips the truth value of an atom. If in lazy mode, the truth value
   * is propagated to the database.
   * 
   * @param atomIdx Index of atom to be flipped.
   * @param blockIdx Index of block in which the atom is in, or -1 if not in one.
   */  
  void flipAtomValue(const int& atomIdx, const int& blockIdx)
  {
    bool opposite = !atom_[atomIdx];
    bool activate = true;
    setValueOfAtom(atomIdx, opposite, activate, blockIdx);
  }

  /**
   * Calculates the improvement achieved by flipping an atom. This is the cost
   * of clauses which flipping the atom makes good minus the cost of clauses
   * which flipping the atom makes bad. In lazy mode, if the atom is not
   * active, then the atom is activated and the makecost and breakcost are
   * updated.
   * 
   * @param atomIdx Index of atom to flip.
   * 
   * @return Improvement achieved by flipping the atom.
   */
  long double getImprovementByFlipping(const int& atomIdx)
  {
    if (!breakHardClauses_ && breakCost_[atomIdx] >= hardWt_)
      return -LDBL_MAX;
    long double improvement = makeCost_[atomIdx] - breakCost_[atomIdx];
    return improvement;
  }
  
  /**
   * Set the active status of an atom to true in the database if in lazy mode.
   * 
   * @atomIdx Index of the atom to be set to active.
   */
  void setActive(const int& atomIdx)
  {
    if (lazy_)
    {
      Predicate* p =
        gndPredHashArray_[atomIdx - 1]->createEquivalentPredicate(domain_);
      domain_->getDB()->setActiveStatus(p, true);
      activeAtoms_++;
      delete p;
    }
  }

  /**
   * If in lazy mode, an atom is activated and all clauses activated by this
   * atom are added to the state. If in eager mode, nothing happens.
   * 
   * @param atomIdx Index of atom to be activated.
   * @param ignoreActivePreds If true, active atoms are ignored when getting
   * active clauses (just gets unsatisfied clauses).
   * @param groundOnly If true, atom itself is not activated, just its active
   * clauses are retrieved.
   */
  bool activateAtom(const int& atomIdx, const bool& ignoreActivePreds,
                    const bool& groundOnly)
  {
    /***
     * . [IF MCSAT] 
     * - check neg-wt clauses
     * . if alive, fix
     * . else mark dead
     * - if not fixed, add all pos-wt, killClauses
     * - else add dead
     * . [ELSE] add all
     ***/
    if (lazy_ && !isActive(atomIdx))
    {
      if (vsdebug)
      {
        cout << "\tActivating ";
        gndPredHashArray_[atomIdx-1]->print(cout,domain_);
        cout << " " <<
        domain_->getDB()->getActiveStatus(gndPredHashArray_[atomIdx-1]) << endl;
      }

        // Set atom to true first, see if avoid getting clauses already active
      bool needToFlipBack = false;
      if (!atom_[atomIdx])
      {
        bool activate = false;
        setValueOfAtom(atomIdx, true, activate, -1);
        updateMakeBreakCostAfterFlip(atomIdx);
        needToFlipBack = true;
      } 

        // gather active clauses
      Predicate* p =
        gndPredHashArray_[atomIdx - 1]->createEquivalentPredicate(domain_);

      bool isFixed = false;
      Array<GroundClause*> unsatClauses;
      Array<GroundClause*> deadClauses;
      Array<GroundClause*> toAddClauses;

        // Using newClauses_ causes problem since it's shared across the board
      getActiveClauses(p, unsatClauses, true, ignoreActivePreds);

        // if MC-SAT: check neg-wt or unit clauses and see if can fix		  
      if (useThreshold_)
      {
          // first append deadClauses, then pos-wt
        for (int ni = 0; ni < unsatClauses.size(); ni++)
        {
          GroundClause* c = unsatClauses[ni];
          if (c->getWt() < 0)
          {  // Neg. weighted clause
              // since newClauses_ excludes clauses in memory, c must be
              // goodInPrevious
            double threshold = RAND_MAX*(1 - exp(c->getWt()));
            if (random() <= threshold)
            {
                // fix
              isFixed = true;

                // append dead and fixed clauses first, so fixed atoms do not
                // activate them unnecessarily
              if (deadClauses.size() > 0)
                addNewClauses(ADD_CLAUSE_DEAD, deadClauses);

              Array<GroundClause*> fixedClauses;
              fixedClauses.append(c);
              addNewClauses(ADD_CLAUSE_SAT, fixedClauses);

                // --- fix the atoms
              for (int pi = 0; pi < c->getNumGroundPredicates(); pi++)
              {
                int lit = c->getGroundPredicateIndex(pi);
                fixAtom(abs(lit), (lit < 0));
              }

                // - can quit now
              break;
            }
            else
            {
                // dead
              deadClauses.append(c);
            }
          }
          else if (c->getNumGroundPredicates() == 1)
          {  // Unit clause
            double threshold = RAND_MAX*(1 - exp(-c->getWt()));
            if (random() <= threshold)
            {
                // fix
              isFixed = true;

                // append dead and fixed clauses first, so fixed atoms do not
                // activate them unnecessarily
              if (deadClauses.size() > 0)
                addNewClauses(ADD_CLAUSE_DEAD, deadClauses);

              Array<GroundClause*> fixedClauses;
              fixedClauses.append(c);
              addNewClauses(ADD_CLAUSE_SAT, fixedClauses);

                // --- fix the atoms
              int lit = c->getGroundPredicateIndex(0);
              fixAtom(abs(lit), (lit > 0)); // this is positive wt

                // - can quit now
              break;
            }
            else
            {
                // dead
              deadClauses.append(c);
            }
          }
          else toAddClauses.append(c);
        }
      }

        // Add the clauses and preds and fill info arrays
      if (!isFixed)
      {
        if (deadClauses.size() > 0)
          addNewClauses(ADD_CLAUSE_DEAD, deadClauses);

        if (useThreshold_)
          addNewClauses(ADD_CLAUSE_REGULAR, toAddClauses);
        else // lazy-sat
          addNewClauses(ADD_CLAUSE_REGULAR, unsatClauses);

          // To avoid append already-active: if not fixed, might need to
          // flip back
        if (needToFlipBack)
        {
          bool activate = false;
          setValueOfAtom(atomIdx, false, activate, -1);
          updateMakeBreakCostAfterFlip(atomIdx);
        }

        if (!groundOnly)
        {
            // Set active status in db
          domain_->getDB()->setActiveStatus(p, true);
          activeAtoms_++;
        }
      }

      delete p;
      unsatClauses.clear();
      deadClauses.clear();
      toAddClauses.clear();

      return !isFixed; // if !isFixed, then activated
    }
    return true;
  }

  /**
   * Set the inference mode currently being used
   */
  void setInferenceMode(int mode)
  {
    if (mode == inferenceMode_) return;
    inferenceMode_ = mode;
    if (inferenceMode_ == MODE_MWS)
    {		  
      for (int i = 0; i < gndClauses_->size(); i++)
      {
        if ((*gndClauses_)[i]->isHardClause())
          clauseCost_[i] = hardWt_;
        else
          clauseCost_[i] = (*gndClauses_)[i]->getWt();
      }
      initMakeBreakCostWatch();
    }
    else
    {
      if (inferenceMode_ == MODE_HARD) eliminateSoftClauses();
      else if (inferenceMode_ == MODE_SAMPLESAT) makeUnitCosts();
    }
  }

  /**
   * Returns true, if an atom is fixed to true or false.
   */
  bool isFixedAtom(int atomIdx)
  {
    return fixedAtom_[atomIdx];
  }

  bool isSatisfiedClause(const int& clauseIdx) 
  {
    return isSatisfied_[clauseIdx];
  }

  void setSatisfiedClause(const int& clauseIdx) 
  {
    isSatisfied_[clauseIdx] = true;
  }

  // set prevStat array for mcsat test-M
  void updatePrevSatisfied()
  {
    prevSatisfiedClause_.clearAndCompress();
    prevSatisfiedClause_.growToSize(clause_.size(),false);
    for (int i = 0; i < clause_.size(); i++)
	{
      long double cost = clauseCost_[i];
      bool isTrue = false;
      for (int j = 0; j < getClauseSize(i); j++)
      {
        if (isTrueLiteral(clause_[i][j]))
        { // ij is true lit
          isTrue = true;
          break;
        }
      }
      prevSatisfiedClause_[i] = ((isTrue && cost>0) || (!isTrue && cost<0));
    }
  }

  /**
   * Perform unit propagation on the current clauses.
   */
  void unitPropagation()
  {
    // ---------------------------------- //
    // check neg-wt clauses
    // ---------------------------------- //
    for (int i = 0; i < getNumClauses(); i++)
    {
        // Don't look at dead clauses
        // pre_wt_unsat -> dead		
      if (isDeadClause(i) || isSatisfiedClause(i) || getClauseCost(i) >= 0)
        continue;

        // fix all atoms
      Array<int> atoms = getAtomsInClause(i);
      for (int j = 0; j < atoms.size(); j++)
      {
        int lit = atoms[j];
        bool value = (lit < 0)? true : false;
        fixAtom(abs(lit), value);
      }
    }

    // ---------------------------------- //
    // prop among pos-wt clauses
    // ---------------------------------- //
    bool done = false;
      // We are done when no new fixed atoms
    while (!done)
    {
      if (vsdebug) cout << endl << endl;
      done = true;

        // Note that getNumClauses() could change between the runs
        // for fixAtom may activate
      for (int ci = 0; ci < getNumClauses(); ci++)
      {
          // skip if irrelevant
        if (isDeadClause(ci) || isSatisfiedClause(ci) || getClauseCost(ci) <= 0)
          continue;

          // check if satisfied, and number of non-fixed atoms
        int numNonfixedAtoms = 0;
        int nonfixedAtom = 0;
	  
        bool isSat = false;
        for (int li = 0; li < getClauseSize(ci); li++)
        {
          int lit = getAtomInClause(li,ci);
          int fixedValue = fixedAtom_[abs(lit)];

          if (fixedValue==0)
          {
            numNonfixedAtoms++;
            nonfixedAtom = lit;
            continue;
          }

          if ((fixedValue == 1 && lit > 0) || (fixedValue == -1 && lit < 0))
          {
            isSat = true;
            break;
          }
        }
        
        if (isSat) setSatisfiedClause(ci);
        else if (numNonfixedAtoms == 1)
        {
          fixAtom(abs(nonfixedAtom), (nonfixedAtom > 0) ? true : false);
          done = false;
        }
      }
    }
      // Done with UP, save this state
    saveLowState();
    if (vsdebug) cout << ">>> [vs.unitpropagation] DONE" << endl;
  }

  /**
   * Adds new clauses to the state.
   * 
   * @param addType How the clauses are being added (see ADD_CLAUSE_*)
   * @param clauses Array of ground clauses to be added to the state.
   */
  void addNewClauses(int addType, Array<GroundClause*> & clauses)
  {
    if (vsdebug)
      cout << "Adding " << clauses.size() << " new clauses.." << endl;

      // Must be in mc-sat to add dead or sat
    if (!useThreshold_ &&
        (addType == ADD_CLAUSE_DEAD || addType == ADD_CLAUSE_SAT))
    {
      cout << ">>> [ERR] add_dead/sat but useThreshold_ is false" << endl;
      exit(0);
    }

      // Store the old number of clauses and atoms
    int oldNumClauses = getNumClauses();
    int oldNumAtoms = getNumAtoms();

      //gndClauses_->append(clauses);
    for (int i = 0; i < clauses.size(); i++)
    {
      gndClauses_->append(clauses[i]);
      clauses[i]->appendToGndPreds(&gndPredHashArray_);
    }

    gndPreds_->growToSize(gndPredHashArray_.size());

    int numAtoms = getNumAtoms();
    int numClauses = getNumClauses();
      // If no new atoms or clauses have been added, then do nothing
    if (numAtoms == oldNumAtoms && numClauses == oldNumClauses) return;

    if (vsdebug) cout << "Clauses: " << numClauses << endl;
      // atomIdx starts at 1
    atom_.growToSize(numAtoms + 1, false);

    makeCost_.growToSize(numAtoms + 1, 0.0);
    breakCost_.growToSize(numAtoms + 1, 0.0);
    lowAtom_.growToSize(numAtoms + 1, false);
    fixedAtom_.growToSize(numAtoms + 1, 0);

	  // Copy ground preds to gndPreds_ and set values in atom and lowAtom
    for (int i = oldNumAtoms; i < gndPredHashArray_.size(); i++)
    {
      (*gndPreds_)[i] = gndPredHashArray_[i];

      if (vsdebug)
      {
        cout << "New pred " << i + 1 << ": ";
        (*gndPreds_)[i]->print(cout, domain_);
        cout << endl;
      }

      lowAtom_[i + 1] = atom_[i + 1] = 
        (domain_->getDB()->getValue((*gndPreds_)[i]) == TRUE) ? true : false;
    }
    clauses.clear();

    clause_.growToSize(numClauses);
    clauseCost_.growToSize(numClauses);
    falseClause_.growToSize(numClauses);
    whereFalse_.growToSize(numClauses);
    numTrueLits_.growToSize(numClauses);
    watch1_.growToSize(numClauses);
    watch2_.growToSize(numClauses);
    isSatisfied_.growToSize(numClauses, false);
    deadClause_.growToSize(numClauses, false);
    threshold_.growToSize(numClauses, false);
  
    occurence_.growToSize(2*numAtoms + 1);

    IntBoolPairItr itr;
    IntBoolPair *clauseFrequencies;

    for (int i = oldNumClauses; i < numClauses; i++)
    {
      GroundClause* gndClause = (*gndClauses_)[i];

      if (vsdebug)
      {
        cout << "New clause " << i << ": ";
        gndClause->print(cout, domain_, &gndPredHashArray_);
        clauseFrequencies = gndClause->getClauseFrequencies();
        cout << " from fo clauses" << endl;
        if (clauseFrequencies)
        {
          for (itr = clauseFrequencies->begin();
               itr != clauseFrequencies->end(); itr++)
          {
            int clauseno = itr->first;
            int frequency = itr->second.first;
            bool invertWt = itr->second.second;
            cout << clauseno << " " << frequency << " invertWt " << invertWt
                 << endl;
          }
        }
      }

        // Set thresholds for clause selection
      if (gndClause->isHardClause()) threshold_[i] = RAND_MAX;
      else
      {
        double w = gndClause->getWt();
        threshold_[i] = RAND_MAX*(1 - exp(-abs(w)));
        if (vsdebug)
        {
          cout << "Weight: " << w << endl;            
        }
      }
      if (vsdebug)
        cout << "Threshold: " << threshold_[i] << endl;            

      int numGndPreds = gndClause->getNumGroundPredicates();
      clause_[i].growToSize(numGndPreds);
      for (int j = 0; j < numGndPreds; j++) 
      {
        int lit = gndClause->getGroundPredicateIndex(j);
        clause_[i][j] = lit;
        int litIdx = 2*abs(lit) - (lit > 0);
        occurence_[litIdx].append(i);
      }

        // Hard clause weight has been previously determined
      if (inferenceMode_==MODE_SAMPLESAT)
      {
        if (gndClause->getWt()>0) clauseCost_[i] = 1;
        else clauseCost_[i] = -1;
      }
      else if (gndClause->isHardClause())
        clauseCost_[i] = hardWt_;
      else
        clauseCost_[i] = gndClause->getWt();

        // filter hard clauses
      if (inferenceMode_ == MODE_HARD && !gndClause->isHardClause())
      {
          // mark dead
        deadClause_[i] = true;
      }
    }

    if (addType == ADD_CLAUSE_DEAD)
    {
        // set to dead: no need for initMakeBreakCost, won't impact anyone
      for (int i = oldNumClauses; i < numClauses; i++)
      {
        deadClause_[i] = true;
      }
    }
    else if (addType == ADD_CLAUSE_SAT)
    {
        // set to dead: no need for initMakeBreakCost, won't impact anyone
      for (int i = oldNumClauses; i < numClauses; i++)
      {
        isSatisfied_[i]=true;
      }
    }
    else if (addType == ADD_CLAUSE_REGULAR)
    {
      if (useThreshold_)
        killClauses(oldNumClauses);
      else
        initMakeBreakCostWatch(oldNumClauses);
    }

    if (vsdebug) 
      cout << "Done adding new clauses.." << endl;
  }

  /**
   * Updates isSatisfiedClause after an atom is fixed.
   * ASSUME: fixed val already applied, and updateMakeBreakCost called
   * 
   * @param toFix Index of atom fixed
   */
  void updateSatisfiedClauses(const int& toFix)
  {
      // Already flipped
    bool toFlipValue = getValueOfAtom(toFix);
      
    register int clauseIdx;
    int sign;
    int litIdx;
    if (toFlipValue)
      sign = 1;
    else
      sign = 0;

      // Set isSatisfied: all clauses in which the atom occurs as a true literal
    litIdx = 2*toFix - sign;
    Array<int>& posOccArray = getOccurenceArray(litIdx);
    for (int i = 0; i < posOccArray.size(); i++)
    {
      clauseIdx = posOccArray[i];
        // Don't look at dead clauses
        // ignore already sat
      if (deadClause_[clauseIdx] || isSatisfied_[clauseIdx]) continue;

      if (getClauseCost(clauseIdx) < 0) 
      {
        cout << "ERR: in MC-SAT, active neg-wt clause (" << clauseIdx
             << ") is sat by fixed "<<endl;
        exit(0);
      }
      isSatisfied_[clauseIdx] = true;  
    }
  }

  /**
   * Updates cost after flip.
   * Assume: already activated and flipped.
   * Called by flipAtom and fixAtom (if fixed to the opposite val)
   */
  void updateMakeBreakCostAfterFlip(const int& toFlip)
  {	  
      // Already flipped
    bool toFlipValue = !getValueOfAtom(toFlip);

    register int clauseIdx;
    int sign;
    int oppSign;
    int litIdx;
    if (toFlipValue)
      sign = 1;
    else
      sign = 0;
    oppSign = sign ^ 1;

      // Update all clauses in which the atom occurs as a true literal
    litIdx = 2*toFlip - sign;
    Array<int>& posOccArray = getOccurenceArray(litIdx);

    for (int i = 0; i < posOccArray.size(); i++)
    {
      clauseIdx = posOccArray[i];
        // Don't look at dead clauses and sat clauses
      if (deadClause_[clauseIdx] || isSatisfied_[clauseIdx]) continue;

        // The true lit became a false lit
      int numTrueLits = decrementNumTrueLits(clauseIdx);
      long double cost = getClauseCost(clauseIdx);
      int watch1 = getWatch1(clauseIdx);
      int watch2 = getWatch2(clauseIdx);

        // 1. If last true lit was flipped, then we have to update
        // the makecost / breakcost info accordingly        
      if (numTrueLits == 0)
      {
          // Pos. clause
        if (cost >= 0)
        {
            // Add this clause as false in the state
          addFalseClause(clauseIdx);
            // Decrease toFlip's breakcost (add neg. cost)
          addBreakCost(toFlip, -cost);
            // Increase makecost of all vars in clause (add pos. cost)
          addMakeCostToAtomsInClause(clauseIdx, cost);
        }
          // Neg. clause
        else
        {
          assert(cost < 0);
            // Remove this clause as false in the state
          removeFalseClause(clauseIdx);
            // Increase breakcost of all vars in clause (add pos. cost)
          addBreakCostToAtomsInClause(clauseIdx, -cost);        
            // Decrease toFlip's makecost (add neg. cost)
          addMakeCost(toFlip, cost);
        }
      }
        // 2. If there is now one true lit left, then move watch2
        // up to watch1 and increase the breakcost / makecost of watch1
      else if (numTrueLits == 1)
      {
        if (watch1 == toFlip)
        {
          assert(watch1 != watch2);
          setWatch1(clauseIdx, watch2);
          watch1 = getWatch1(clauseIdx);
        }
          // Pos. clause: Increase toFlip's breakcost (add pos. cost)
        if (cost >= 0)
        {
          addBreakCost(watch1, cost);
        }
            // Neg. clause: Increase toFlip's makecost (add pos. cost)
        else
        {
          assert(cost < 0);
          addMakeCost(watch1, -cost);
        }
      }
        // 3. If there are 2 or more true lits left, then we have to
        // find a new true lit to watch if one was flipped
      else
      {
          /* numtruelit[clauseIdx] >= 2 */
          // If watch1[clauseIdx] has been flipped
        if (watch1 == toFlip)
        {
            // find a different true literal to watch
          int diffTrueLit = getTrueLiteralOtherThan(clauseIdx, watch1, watch2);
          setWatch1(clauseIdx, diffTrueLit);
        }
          // If watch2[clauseIdx] has been flipped
        else if (watch2 == toFlip)
        {
            // find a different true literal to watch
          int diffTrueLit = getTrueLiteralOtherThan(clauseIdx, watch1, watch2);
          setWatch2(clauseIdx, diffTrueLit);
        }
      }
    }

      // Update all clauses in which the atom occurs as a false literal
    litIdx = 2*toFlip - oppSign;
    Array<int>& negOccArray = getOccurenceArray(litIdx);
    for (int i = 0; i < negOccArray.size(); i++)
    {
      clauseIdx = negOccArray[i];
        // Don't look at dead clauses or satisfied clauses
      if (deadClause_[clauseIdx] || isSatisfied_[clauseIdx]) continue;

        // The false lit became a true lit
      int numTrueLits = incrementNumTrueLits(clauseIdx);
      long double cost = getClauseCost(clauseIdx);
      int watch1 = getWatch1(clauseIdx);

        // 1. If this is the only true lit, then we have to update
        // the makecost / breakcost info accordingly        
      if (numTrueLits == 1)
      {
          // Pos. clause
        if (cost >= 0)
        {
            // Remove this clause as false in the state
          removeFalseClause(clauseIdx);
            // Increase toFlip's breakcost (add pos. cost)
          addBreakCost(toFlip, cost);        
            // Decrease makecost of all vars in clause (add neg. cost)
          addMakeCostToAtomsInClause(clauseIdx, -cost);
        }
          // Neg. clause
        else
        {
          assert(cost < 0);
            // Add this clause as false in the state
          addFalseClause(clauseIdx);
            // Decrease breakcost of all vars in clause (add neg. cost)
          addBreakCostToAtomsInClause(clauseIdx, cost);
            // Increase toFlip's makecost (add pos. cost)
          addMakeCost(toFlip, -cost);
        }
          // Watch this atom
        setWatch1(clauseIdx, toFlip);
      }
        // 2. If there are now exactly 2 true lits, then watch second atom
        // and update breakcost
      else if (numTrueLits == 2)
      {
        if (cost >= 0)
        {
            // Pos. clause
            // Decrease breakcost of first atom being watched (add neg. cost)
          addBreakCost(watch1, -cost);
        }
        else
        {
            // Neg. clause
          assert(cost < 0);
            // Decrease makecost of first atom being watched (add neg. cost)
          addMakeCost(watch1, cost);
        }
          // Watch second atom
        setWatch2(clauseIdx, toFlip);
      }
    }
  }


  /**
   * Checks if an atom is active.
   * 
   * @param atomIdx Index of atom to be checked.
   * @return true, if atom is active, otherwise false.
   */  
  bool isActive(const int& atomIdx)
  {
    return domain_->getDB()->getActiveStatus(gndPredHashArray_[atomIdx-1]);
  }
  
  /**
   * Checks if an atom is active.
   * 
   * @param pred Predicate to be checked.
   * @return true, if atom is active, otherwise false.
   */  
  bool isActive(const Predicate* pred)
  {
    return domain_->getDB()->getActiveStatus(pred);
  }

  /**
   * Retrieves the occurence array of an atom.
   */
  Array<int>& getOccurenceArray(const int& idx)
  {
    return occurence_[idx];
  }
  
  /**
   * Increments the number of true literals in a clause.
   */
  int incrementNumTrueLits(const int& clauseIdx)
  {
    return ++numTrueLits_[clauseIdx];
  }

  /**
   * Decrements the number of true literals in a clause.
   */
  int decrementNumTrueLits(const int& clauseIdx)
  {
    return --numTrueLits_[clauseIdx];
  }

  /**
   * Retrieves the number of true literals in a clause.
   */
  int getNumTrueLits(const int& clauseIdx)
  {
    return numTrueLits_[clauseIdx];
  }

  /**
   * Retrieves the cost associated with a clause.
   */
  long double getClauseCost(const int& clauseIdx)
  {
    return clauseCost_[clauseIdx];
  }
  
  /**
   * Retrieves the atoms in a clauses as an Array of ints.
   */
  Array<int>& getAtomsInClause(const int& clauseIdx)
  {
    return clause_[clauseIdx];
  }

  /**
   * Marks a clause as false in the state.
   */
  void addFalseClause(const int& clauseIdx)
  {
    falseClause_[numFalseClauses_] = clauseIdx;
    whereFalse_[clauseIdx] = numFalseClauses_;
    numFalseClauses_++;
    costOfFalseClauses_ += abs(clauseCost_[clauseIdx]);
  }
  
  /**
   * Unmarks a clause as false in the state.
   */
  void removeFalseClause(const int& clauseIdx)
  {
    numFalseClauses_--;
    falseClause_[whereFalse_[clauseIdx]] = falseClause_[numFalseClauses_];
    whereFalse_[falseClause_[numFalseClauses_]] = whereFalse_[clauseIdx];
    costOfFalseClauses_ -= abs(clauseCost_[clauseIdx]);
  }

  /**
   * Increases the breakcost of an atom.
   */
  void addBreakCost(const int& atomIdx, const long double& cost)
  {
    breakCost_[atomIdx] += cost;
  }

  /**
   * Decreases the breakcost of an atom.
   */
  void subtractBreakCost(const int& atomIdx, const long double& cost)
  {
    breakCost_[atomIdx] -= cost;
  }

  /**
   * Increases breakCost of all atoms in a given clause.
   * 
   * @param clauseIdx Index of clause whose atoms' breakCost is altered.
   * @param cost Cost to be added to atoms' breakCost.
   */
  void addBreakCostToAtomsInClause(const int& clauseIdx,
                                   const long double& cost)
  {
    register int size = getClauseSize(clauseIdx);
    for (int i = 0; i < size; i++)
    {
      register int lit = clause_[clauseIdx][i];
      breakCost_[abs(lit)] += cost;
    }
  }

  /**
   * Decreases breakCost of all atoms in a given clause.
   * 
   * @param clauseIdx Index of clause whose atoms' breakCost is altered.
   * @param cost Cost to be subtracted from atoms' breakCost.
   */
  void subtractBreakCostFromAtomsInClause(const int& clauseIdx,
                                          const long double& cost)
  {
    register int size = getClauseSize(clauseIdx);
    for (int i = 0; i < size; i++)
    {
      register int lit = clause_[clauseIdx][i];
      breakCost_[abs(lit)] -= cost;
    }
  }

  /**
   * Increases makeCost of an atom.
   * 
   * @param atomIdx Index of atom whose makeCost is altered.
   * @param cost Cost to be added to atom's makeCost.
   */
  void addMakeCost(const int& atomIdx, const long double& cost)
  {
    makeCost_[atomIdx] += cost;
  }

  /**
   * Decreases makeCost of an atom.
   * 
   * @param atomIdx Index of atom whose makeCost is altered.
   * @param cost Cost to be subtracted from atom's makeCost.
   */
  void subtractMakeCost(const int& atomIdx, const long double& cost)
  {
    makeCost_[atomIdx] -= cost;
  }

  /**
   * Increases makeCost of all atoms in a given clause.
   * 
   * @param clauseIdx Index of clause whose atoms' makeCost is altered.
   * @param cost Cost to be added to atoms' makeCost.
   */
  void addMakeCostToAtomsInClause(const int& clauseIdx,
                                  const long double& cost)
  {
    register int size = getClauseSize(clauseIdx);
    for (int i = 0; i < size; i++)
    {
      register int lit = clause_[clauseIdx][i];
      makeCost_[abs(lit)] += cost;
    }
  }

  /**
   * Decreases makeCost of all atoms in a given clause.
   * 
   * @param clauseIdx Index of clause whose atoms' makeCost is altered.
   * @param cost Cost to be subtracted from atoms' makeCost.
   */
  void subtractMakeCostFromAtomsInClause(const int& clauseIdx,
                                         const long double& cost)
  {
    register int size = getClauseSize(clauseIdx);
    for (int i = 0; i < size; i++)
    {
      register int lit = clause_[clauseIdx][i];
      makeCost_[abs(lit)] -= cost;
    }
  }

  /**
   * Retrieves a true literal in a clause other than the two given.
   * 
   * @param clauseIdx Index of clause from which literal is retrieved.
   * @param atomIdx1 Index of first atom excluded from search.
   * @param atomIdx2 Index of second atom excluded from search.
   * 
   * @return Index of atom found.
   */
  const int getTrueLiteralOtherThan(const int& clauseIdx,
                                    const int& atomIdx1,
                                    const int& atomIdx2)
  {
    register int size = getClauseSize(clauseIdx);
    for (int i = 0; i < size; i++)
    {
      register int lit = clause_[clauseIdx][i];
      register int v = abs(lit);
      if (isTrueLiteral(lit) && v != atomIdx1 && v != atomIdx2)
        return v;
    }
      // If we're here, then no other true lit exists
    assert(false);
    return -1;
  }
  
  /**
   * Checks if a literal is true (Negated and false or non-negated an true).
   */
  const bool isTrueLiteral(const int& literal)
  {
    return ((literal > 0) == atom_[abs(literal)]);
  }

  /**
   * Retrieves the absolute index of the nth atom in a clause.
   */
  const int getAtomInClause(const int& atomIdxInClause, const int& clauseIdx)
  {
    return clause_[clauseIdx][atomIdxInClause];
  }

  /**
   * Retrieves the absolute index of a random atom in a clause.
   */
  const int getRandomAtomInClause(const int& clauseIdx)
  {
    return clause_[clauseIdx][random()%getClauseSize(clauseIdx)];
  }

  /**
   * Retrieves the index of a random true literal in a clause.
   * 
   * @param clauseIdx Index of clause from which the literal is retrieved.
   * @return Index of the atom retrieved.
   */  
  const int getRandomTrueLitInClause(const int& clauseIdx)
  {
    assert(numTrueLits_[clauseIdx] > 0);
    int trueLit = random()%numTrueLits_[clauseIdx];
    int whichTrueLit = 0;
    for (int i = 0; i < getClauseSize(clauseIdx); i++)
    {
      int lit = clause_[clauseIdx][i];
      int atm = abs(lit);
        // True literal
      if (isTrueLiteral(lit))
        if (trueLit == whichTrueLit++)
          return atm;
    }
      // If we're here, then no other true lit exists
    assert(false);
    return -1;
  }

  const double getMaxClauseWeight()
  {
    double maxWeight = 0.0;
    for (int i = 0; i < getNumClauses(); i++)
    {
      double weight = abs(clauseCost_[i]);
      if (weight > maxWeight) maxWeight = weight;
    }
    return maxWeight;
  }
  
  const long double getMakeCost(const int& atomIdx)
  {
    return makeCost_[atomIdx];
  }
   
  const long double getBreakCost(const int& atomIdx)
  {
    return breakCost_[atomIdx];
  }
   
  const int getClauseSize(const int& clauseIdx)
  {
    return clause_[clauseIdx].size();
  }

  const int getWatch1(const int& clauseIdx)
  {
    return watch1_[clauseIdx];
  }

  void setWatch1(const int& clauseIdx, const int& atomIdx)
  {
    watch1_[clauseIdx] = atomIdx;
  }
  
  const int getWatch2(const int& clauseIdx)
  {
    return watch2_[clauseIdx];
  }

  void setWatch2(const int& clauseIdx, const int& atomIdx)
  {
    watch2_[clauseIdx] = atomIdx;
  }
  
  /** 
   * Returns the index of the block which the atom with index atomIdx
   * is in. If not in any, returns -1.
   */
  const int getBlockIndex(const int& atomIdx)
  {
    const GroundPredicate* gndPred = (*gndPreds_)[atomIdx];
    return domain_->getBlock(gndPred);
  }

  /**
   * Returns the cost of bad clauses in the optimum state.
   */
  const long double getLowCost()
  {
    return lowCost_; 
  }

  /**
   * Returns the number of bad clauses in the optimum state.
   */
  const int getLowBad()
  {
    return lowBad_;
  }

  /**
   * Turns all costs into units. Positive costs are converted to 1, negative
   * costs are converted to -1
   */
  void makeUnitCosts()
  {
    for (int i = 0; i < clauseCost_.size(); i++)
    {
      if (clauseCost_[i] >= 0) clauseCost_[i] = 1.0;
      else
      {
        assert(clauseCost_[i] < 0);
        clauseCost_[i] = -1.0;
      }
    }
    if (vsdebug) cout << "Made unit costs" << endl;
    initMakeBreakCostWatch();
  }

  /**
   * Save current assignment as an optimum state.
   */
  void saveLowState()
  {
    if (vsdebug) cout << "Saving low state: " << endl;
    for (int i = 1; i <= getNumAtoms(); i++)
    {
      lowAtom_[i] = atom_[i];
      if (vsdebug) cout << lowAtom_[i] << endl;
    }
    lowCost_ = costOfFalseClauses_;
    lowBad_ = numFalseClauses_;
  }

  /**
   * Returns index of an atom if it is true and fixed in block, otherwise -1    
   */
  int getTrueFixedAtomInBlock(const int& blockIdx)
  {
    const Predicate* truePred = domain_->getTruePredInBlock(blockIdx);
    if (truePred)
    {
      if (vsdebug)
      {
        cout << "True pred in block " << blockIdx << ": ";
        truePred->printWithStrVar(cout, domain_);
        cout << endl;
      }
    
      GroundPredicate* trueGndPred = new GroundPredicate((Predicate*)truePred);

      int atomIdx = gndPredHashArray_.find(trueGndPred);
      delete trueGndPred;
      if (atomIdx > -1 && fixedAtom_[atomIdx + 1] > 0)
        return atomIdx;
    }
    return -1;
  }

  const GroundPredicateHashArray* getGndPredHashArrayPtr() const
  {
    return &gndPredHashArray_;
  }

  const GroundPredicateHashArray* getUnePreds() const
  {
    return unePreds_;
  }

  const GroundPredicateHashArray* getKnePreds() const
  {
    return knePreds_;
  }

  const Array<TruthValue>* getKnePredValues() const
  {
    return knePredValues_;
  }

  /**
   * Sets the weight of a ground clause to the sum of its parents' weights.
   */
  void setGndClausesWtsToSumOfParentWts()
  {
    for (int i = 0; i < gndClauses_->size(); i++)
    {
      GroundClause* gndClause = (*gndClauses_)[i];
      gndClause->setWtToSumOfParentWts(mln_);
      if (gndClause->isHardClause())
        clauseCost_[i] = hardWt_;
      else
        clauseCost_[i] = gndClause->getWt();
      
      if (vsdebug) cout << "Setting cost of clause " << i << " to "
                        << clauseCost_[i] << endl;

        // Set thresholds for clause selection
      if (gndClause->isHardClause()) threshold_[i] = RAND_MAX;
      else
      {
        double w = gndClause->getWt();
        threshold_[i] = RAND_MAX*(1 - exp(-abs(w)));
        if (vsdebug)
        {
          cout << "Weight: " << w << endl;            
        }
      }
      if (vsdebug)
        cout << "Threshold: " << threshold_[i] << endl;
    }
  }

  /**
   * Gets the number of (true or false) clause groundings in this state. If
   * eager, the first order clause frequencies in the mrf are used. If lazy, the
   * f.o. clause frequencies computed when activating are used.
   * 
   * @param numGndings Array being filled with values.
   * @param clauseCnt Number of first-order clauses (size of numGndings)
   * @param tv If true, number of true groundings are retrieved, otherwise
   * false groundings are retrieved.
   */
  void getNumClauseGndings(Array<double>* const & numGndings, bool tv)
  {
      // If lazy, then all groundings of pos. clauses not materialized in this
      // state are true (non-materialized groundings of neg. clauses are false),
      // so if tv is false and clause is pos. (or tv is true and clause is
      // neg.), then eager and lazy counts are equivalent. If tv is
      // true and clause is pos. (or tv is false and clause is neg.), then we
      // need to keep track of true *and* false groundings and
      // then subtract from total groundings of the f.o. clause.
      
      // This holds the false groundings when tv is true and lazy.
    Array<double> lazyFalseGndings(numGndings->size(), 0);
    Array<double> lazyTrueGndings(numGndings->size(), 0);

    IntBoolPairItr itr;
    IntBoolPair *clauseFrequencies;
    
      // numGndings should have been initialized with non-negative values
    int clauseCnt = numGndings->size();
    assert(clauseCnt == mln_->getNumClauses());
    for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
      assert ((*numGndings)[clauseno] >= 0);
    
    for (int i = 0; i < gndClauses_->size(); i++)
    {
      GroundClause *gndClause = (*gndClauses_)[i];
      int satLitcnt = 0;
      for (int j = 0; j < gndClause->getNumGroundPredicates(); j++)
      {
        int lit = gndClause->getGroundPredicateIndex(j);
        if (isTrueLiteral(lit)) satLitcnt++;
      }

      clauseFrequencies = gndClause->getClauseFrequencies();
      for (itr = clauseFrequencies->begin();
           itr != clauseFrequencies->end(); itr++)
      {
        int clauseno = itr->first;
        int frequency = itr->second.first;
        bool invertWt = itr->second.second;
          // If flipped unit clause, then we want opposite kind of grounding
        if (invertWt)
        {
            // Want true and is true => don't count it
          if (tv && satLitcnt > 0)
          {
              // Lazy: We need to keep track of false groundings also
            if (lazy_) lazyFalseGndings[clauseno] += frequency;
            continue;
          }
            // Want false and is false => don't count it
          if (!tv && satLitcnt == 0)
          {
              // Lazy: We need to keep track of false groundings also
            if (lazy_) lazyTrueGndings[clauseno] += frequency;
            continue;
          }
        }
        else
        {
            // Want true and is false => don't count it
          if (tv && satLitcnt == 0)
          {
              // Lazy: We need to keep track of false groundings also
            if (lazy_) lazyFalseGndings[clauseno] += frequency;
            continue;
          }
            // Want false and is true => don't count it
          if (!tv && satLitcnt > 0)
          {
              // Lazy: We need to keep track of false groundings also
            if (lazy_) lazyTrueGndings[clauseno] += frequency;
            continue;
          }
        }
        (*numGndings)[clauseno] += frequency;
      }
    }
    
      // Getting true counts in lazy: we have to add the remaining groundings
      // not materialized (they are by definition true groundings if clause is
      // positive, or false groundings if clause is negative)
    if (lazy_)
    {
      for (int c = 0; c < mln_->getNumClauses(); c++)
      {
        const Clause* clause = mln_->getClause(c);
          // Getting true counts and positive clause
        if (tv && clause->getWt() >= 0)
        {
          double totalGndings = domain_->getNumNonEvidGroundings(c);
          assert(totalGndings >= (*numGndings)[c] + lazyFalseGndings[c]);
          double remainingTrueGndings = totalGndings - lazyFalseGndings[c] -
                                        (*numGndings)[c];
          (*numGndings)[c] += remainingTrueGndings;
        }
          // Getting false counts and negative clause
        else if (!tv && clause->getWt() < 0)
        {
          double totalGndings = domain_->getNumNonEvidGroundings(c);
          assert(totalGndings >= (*numGndings)[c] + lazyTrueGndings[c]);
          double remainingFalseGndings = totalGndings - lazyTrueGndings[c] -
                                        (*numGndings)[c];
          (*numGndings)[c] += remainingFalseGndings;
        }
      }
    }

  }

  /**
   * Gets the number of (true or false) clause groundings in this state. If
   * eager, the first order clause frequencies in the mrf are used.
   * 
   * @param numGndings Will hold the number of groundings for each first-order
   * clause.
   * @param clauseCnt Number of first-order clauses whose groundings are being
   * counted.
   * @param tv If true, true groundings are counted, otherwise false groundings.
   * @param unknownPred If pred is marked as unknown, it is ignored in the count
   */
  void getNumClauseGndingsWithUnknown(double numGndings[], int clauseCnt,
                                      bool tv,
                                      const Array<bool>* const& unknownPred)
  {
    // TODO: lazy version
    assert(unknownPred->size() == getNumAtoms());
    IntBoolPairItr itr;
    IntBoolPair *clauseFrequencies;
    
    for (int clauseno = 0; clauseno < clauseCnt; clauseno++)
      numGndings[clauseno] = 0;
    
    for (int i = 0; i < gndClauses_->size(); i++)
    {
      GroundClause *gndClause = (*gndClauses_)[i];
      int satLitcnt = 0;
      bool unknown = false;
      for (int j = 0; j < gndClause->getNumGroundPredicates(); j++)
      {
        int lit = gndClause->getGroundPredicateIndex(j);
        if ((*unknownPred)[abs(lit) - 1])
        {
          unknown = true;
          continue;
        }
        if (isTrueLiteral(lit)) satLitcnt++;
      }

      clauseFrequencies = gndClause->getClauseFrequencies();
      for (itr = clauseFrequencies->begin();
           itr != clauseFrequencies->end(); itr++)
      {
        int clauseno = itr->first;
        int frequency = itr->second.first;
        bool invertWt = itr->second.second;
          // If flipped unit clause, then we want opposite kind of grounding
        if (invertWt)
        {
            // Want true and is true or unknown => don't count it
          if (tv && (satLitcnt > 0 || unknown))
            continue;
            // Want false and is false => don't count it
          if (!tv && satLitcnt == 0)
            continue;
        }
        else
        {
            // Want true and is false => don't count it
          if (tv && satLitcnt == 0)
            continue;
            // Want false and is true => don't count it
          if (!tv && (satLitcnt > 0 || unknown))
            continue;
        }
        numGndings[clauseno] += frequency;
      }
    }
  }

  /**
   * Sets the truth values of all atoms in a block except for the one given.
   * 
   * @param atomIdx Absolute index of atom in block exempt from being set to
   * false.
   * @param blockIdx Index of block whose atoms are set to false.
   */
  void setOthersInBlockToFalse(const int& atomIdx, const int& blockIdx)
  {
    if (vsdebug)
    {
      cout << "Set all in block " << blockIdx << " to false except "
           << atomIdx << endl;
    }
    int blockSize = domain_->getBlockSize(blockIdx);
    for (int i = 0; i < blockSize; i++)
    {
      const Predicate* pred = domain_->getPredInBlock(i, blockIdx);
      GroundPredicate* gndPred = new GroundPredicate((Predicate*)pred);      
      int idx = gndPredHashArray_.find(gndPred);

      if (vsdebug)
      {
        cout << "Gnd pred in block ";
        gndPred->print(cout, domain_);
        cout << " (idx " << idx << ")" << endl;
      }

      delete gndPred;
      delete pred;
      
        // Pred already present: check that it's not fixed
      if (idx >= 0)
      {
          // Atom not the one specified and not fixed
        if (idx != atomIdx && fixedAtom_[idx + 1] == 0)
        {
          if (vsdebug)
            cout << "Set " << idx + 1 << " to false" << endl;

          bool activate = true;
          setValueOfAtom(idx + 1, false, activate, -1);
        }
      }
    }
  }


  /**
   * Fixes an atom to a truth value. This means the atom can not change
   * its truth value again. If the atom has already been fixed to the
   * opposite value, then we have a contradiction and the program terminates.
   * 
   * @param atomIdx Index of atom to be fixed.
   * @param value Truth value to which the atom is fixed.
   */
  void fixAtom(const int& atomIdx, const bool& value)
  {
    assert(atomIdx > 0);
    if (vsdebug)
    {
      cout << "Fixing ";
      (*gndPreds_)[atomIdx - 1]->print(cout, domain_);
      cout << " to " << value << endl;	  
    }

      // Must be in mc-sat
    if (!useThreshold_)
    {
      cout << ">>> [ERR] useThreshold_ is false" << endl;
      exit(0);
    }

      // If already fixed to opp. sense, then contradiction
    if ((fixedAtom_[atomIdx] == 1 && value == false) ||
        (fixedAtom_[atomIdx] == -1 && value == true))
    {
      cout << "Contradiction: Tried to fix atom " << atomIdx <<
      " to true and false ... exiting." << endl;
      exit(0);
    }

      // already fixed
    if (fixedAtom_[atomIdx] != 0) return;

    fixedAtom_[atomIdx] = (value) ? 1 : -1;	
    if (atom_[atomIdx] != value) 
    {
      bool activate = false;
      int blockIdx =  getBlockIndex(atomIdx - 1);
      setValueOfAtom(atomIdx, value, activate, blockIdx);
      updateMakeBreakCostAfterFlip(atomIdx);
    }
    
      // update isSat
    updateSatisfiedClauses(atomIdx);

      // If in a block and fixing to true, fix all others to false
    if (value)
    {
      int blockIdx = getBlockIndex(atomIdx - 1);
      if (blockIdx > -1)
      {
        int blockSize = domain_->getBlockSize(blockIdx);
        for (int i = 0; i < blockSize; i++)
        {
          const Predicate* pred = domain_->getPredInBlock(i, blockIdx);
          GroundPredicate* gndPred = new GroundPredicate((Predicate*)pred);      
          int idx = gndPredHashArray_.find(gndPred);
          delete gndPred;
          delete pred;
      
            // Pred already present: check that it's not fixed
          if (idx >= 0)
          {
              // Atom not the one specified
            if (idx != (atomIdx - 1))
            {
                // If already fixed to true, then contradiction
              if (fixedAtom_[idx + 1] == 1)
              {
                cout << "Contradiction: Tried to fix atom " << idx + 1 <<
                " to true and false ... exiting." << endl;
                exit(0);
              }
                // already fixed
              if (fixedAtom_[idx + 1] == -1) continue;
              if (vsdebug)
              {
                cout << "Fixing ";
                (*gndPreds_)[idx]->print(cout, domain_);
                cout << " to 0" << endl;
              }
              fixedAtom_[idx + 1] = -1; 
              if (atom_[idx + 1] != false) 
              {
                bool activate = false;
                setValueOfAtom(idx + 1, value, activate, blockIdx);
                updateMakeBreakCostAfterFlip(idx + 1);
              }
                // update isSat
              updateSatisfiedClauses(idx + 1);
            }
          }
        }
      }
    }

      // activate clauses falsified
      // need to updateMakeBreakCost for new clauses only, which aren't sat
    if (lazy_ && !isActive(atomIdx) && value)
    {
        // if inactive fixed to true, need to activate falsified
        // inactive clauses gather active clauses
      Predicate* p =
        gndPredHashArray_[atomIdx - 1]->createEquivalentPredicate(domain_);
		
      bool ignoreActivePreds = false; // only get currently unsat ones
      Array<GroundClause*> unsatClauses;
      getActiveClauses(p, unsatClauses, true, ignoreActivePreds);		
		
        // filter out clauses already added (can also add directly, but
        // will elicit warning)
      addNewClauses(ADD_CLAUSE_REGULAR, unsatClauses);

        // Set active status in db
      domain_->getDB()->setActiveStatus(p, true);
      activeAtoms_++;

      delete p;
    }
  }

  /**
   * Simplifies a clause using atoms which have been fixed. If clause is
   * satisfied from the fixed atoms, this is marked in isSatisfied_ and an
   * empty array is returned. If clause is empty and not satisfied, then a
   * contradiction has occured. Otherwise, the simplified clause is returned.
   * 
   * Returned array should be deleted.
   * 
   * @param clauseIdx Index of the clause to be simplified
   * @return Simplified clause
   */
  Array<int>* simplifyClauseFromFixedAtoms(const int& clauseIdx)
  {
    Array<int>* returnArray = new Array<int>;
      // If already satisfied from fixed atoms, then return empty array
    if (isSatisfied_[clauseIdx]) return returnArray;

      // Keeps track of pos. clause being satisfied or 
      // neg. clause being unsatisfied due to fixed atoms
    bool isGood = (clauseCost_[clauseIdx] >= 0) ? false : true;
      // Keeps track of all atoms being fixed to false in a pos. clause
    bool allFalseAtoms = (clauseCost_[clauseIdx] >= 0) ? true : false;
      // Check each literal in clause
    for (int i = 0; i < getClauseSize(clauseIdx); i++)
    {
      int lit = clause_[clauseIdx][i];
      int fixedValue = fixedAtom_[abs(lit)];

      if (clauseCost_[clauseIdx] >= 0)
      { // Pos. clause: check if clause is satisfied
        if ((fixedValue == 1 && lit > 0) ||
            (fixedValue == -1 && lit < 0))
        { // True fixed lit
          isGood = true;
          allFalseAtoms = false;
          returnArray->clear();
          break;
        }
        else if (fixedValue == 0)
        { // Lit not fixed
          allFalseAtoms = false;
          returnArray->append(lit);
        }
      }
      else
      { // Neg. clause:
        assert(clauseCost_[clauseIdx] < 0);
        if ((fixedValue == 1 && lit > 0) ||
            (fixedValue == -1 && lit < 0))
        { // True fixed lit
          cout << "Contradiction: Tried to fix atom " << abs(lit) <<
          " to true in a negative clause ... exiting." << endl;
          exit(0);
        }
        else
        { // False fixed lit or non-fixed lit
          returnArray->append(lit);
            // Non-fixed lit
          if (fixedValue == 0) isGood = false;          
        }
      }
    }
    if (allFalseAtoms)
    {
      cout << "Contradiction: All atoms in clause " << clauseIdx <<
      " fixed to false ... exiting." << endl;
      exit(0);
    }
    if (isGood) isSatisfied_[clauseIdx] = true;
    return returnArray;
  }

  /**
   * Checks if a clause is dead.
   * 
   * @param clauseIdx Index of clause being checked.
   * @return True, if clause is dead, otherwise false.
   */
  const bool isDeadClause(const int& clauseIdx)
  {
    return deadClause_[clauseIdx];
  }

  /**
   * Marks soft clauses as dead.
   */
  void eliminateSoftClauses()
  {
    bool atLeastOneDead = false;
    for (int i = 0; i < getNumClauses(); i++)
    {
      if (!(*gndClauses_)[i]->isHardClause())
      {
        atLeastOneDead = true;
        deadClause_[i] = true;
      }
    } 

	if (atLeastOneDead) initMakeBreakCostWatch();
  }

  /**
   * Marks clauses as dead which were not good in the previous iteration of
   * inference or are not picked according to a weighted coin flip. 
   *
   * @param startClause All clauses with index of this or greater are
   * looked at to be killed.
   */
  void killClauses(const int& startClause)
  {
      // for hard init, no need to killClauses
    if (inferenceMode_ != MODE_HARD)
    {
      for (int i = startClause; i < getNumClauses(); i++)
      {
        GroundClause* clause = (*gndClauses_)[i];
        if ((clauseGoodInPrevious(i)) &&
            (clause->isHardClause() || random() <= threshold_[i]))
        {
          if (vsdebug)
          {
            cout << "Keeping clause "<< i << " ";
            clause->print(cout, domain_, &gndPredHashArray_);
            cout << endl;
          }
          deadClause_[i] = false;
        }
        else
        {
          deadClause_[i] = true;
        }
      }
    }
    
    initMakeBreakCostWatch(startClause);
  }

  
  /**
   * Checks if a clause was good in the previous iteration of inference, i.e.
   * if it is positive and satisfied or negative and unsatisfied.
   * 
   * @param clauseIdx Index of clause being checked.
   * @return true, if clause was good, otherwise false.
   */
  const bool clauseGoodInPrevious(const int& clauseIdx)
  {
      // satisfied: wt-sat
    return (clauseIdx >= prevSatisfiedClause_.size() ||
            prevSatisfiedClause_[clauseIdx]);
  }

  /**
   * Resets all dead clauses to be alive again.
   */
  void resetDeadClauses()
  {
    for (int i = 0; i < deadClause_.size(); i++)
      deadClause_[i] = false;
    initMakeBreakCostWatch();
  }
   
  /**
   * Resets all fixed atoms to be not fixed again.
   */
  void resetFixedAtoms()
  {
    for (int i = 0; i < fixedAtom_.size(); i++)
      fixedAtom_[i] = 0;
    for (int i = 0; i < isSatisfied_.size(); i++)
      isSatisfied_[i] = false;
  }

  void setLazy(const bool& l) { lazy_ = l; }
  const bool getLazy() { return lazy_; }

  void setUseThreshold(const bool& t) { useThreshold_ = t;}
  const bool getUseThreshold() { return useThreshold_; }
  
  long double getHardWt() { return hardWt_; }
  
  const Domain* getDomain() { return domain_; }

  const MLN* getMLN() { return mln_; }

  /**
   * Prints the best state found to a stream.
   * 
   * @param out Stream to which the state is printed.
   */
  void printLowState(ostream& out)
  {
    for (int i = 0; i < getNumAtoms(); i++)
    {
      (*gndPreds_)[i]->print(out, domain_);
      out << " " << lowAtom_[i + 1] << endl;
    }
  }

  /**
   * Prints a ground predicate to a stream.
   * 
   * @param predIndex Index of predicate to be printed.
   * @param out Stream to which predicate is printed.
   */
  void printGndPred(const int& predIndex, ostream& out)
  {
    (*gndPreds_)[predIndex]->print(out, domain_);
  }

  /**
   * Finds and returns the index of a ground predicate.
   * 
   * @param gndPred GroundPredicate to be found.
   * @return Index of the ground predicate, if present; otherwise, -1.
   */
  int getIndexOfGroundPredicate(GroundPredicate* const & gndPred)
  {
    return gndPredHashArray_.find(gndPred);
  }
  
  /**
   * Sets a GroundPredicate to be evidence and sets its truth value. If it is
   * already present as evidence with the given truth value, then nothing
   * happens. If the predicate was a query, then additional clauses may be
   * eliminated. reinit() should be called after this in order to ensure that
   * the clause and atom information is correct.
   * 
   * @param predicate GroundPredicate to be set as evidence.
   * @param trueEvidence The truth value of the predicate is set to this.
   */
  void setAsEvidence(const GroundPredicate* const & predicate,
                     const bool& trueEvidence)
  {
    if (vsdebug)
    {
      cout << "Setting to evidence " ;
      predicate->print(cout, domain_);
      cout << endl;
    }
    Database* db = domain_->getDB();
    int atomIdx = gndPredHashArray_.find((GroundPredicate*)predicate);
      // If already evidence, then check its truth value
    if (atomIdx <= 0)
    {
        // If predicate already evidence with same truth value, then do nothing
      if (db->getValue(predicate) == trueEvidence)
        return;
        
        // Changing truth value of evidence
      if (trueEvidence)
        db->setValue(predicate, TRUE);
      else
        db->setValue(predicate, FALSE);
    }
    else
    {
      Array<int> gndClauseIndexes;      
      int deleted;
      gndClauseIndexes = getNegOccurenceArray(atomIdx + 1);
      gndClauseIndexes.bubbleSort();
        // Keep track of how many clauses deleted, because the indices are
        // are adjusted when we remove an element from HashArray
      deleted = 0;
      for (int i = 0; i < gndClauseIndexes.size(); i++)
      {
          // Atom appears neg. in these clauses, so remove the atom from
          // these clauses if true evidence, or remove clause if false evidence
          // or a unit clause
        if (!trueEvidence ||
            (*gndClauses_)[gndClauseIndexes[i]]->getNumGroundPredicates() == 1)
        {          
          if (vsdebug)
          {
            cout << "Deleting ground clause " << gndClauseIndexes[i] << " ";
            cout << endl;
          }
            // Real index is old index adjusted one lower for every element
            // deleted up until now
          delete (*gndClauses_)[gndClauseIndexes[i] - deleted];
          gndClauses_->removeItem(gndClauseIndexes[i] - deleted);
          deleted++;
        }
        else
        {
          if (vsdebug)
          {
            cout << "Removing gnd pred " << -(atomIdx + 1)
                 << " from ground clause " << gndClauseIndexes[i] << endl;
          }
          (*gndClauses_)[gndClauseIndexes[i]]->removeGndPred(-(atomIdx + 1));
        }
      }

      gndClauseIndexes = getPosOccurenceArray(atomIdx + 1);
      gndClauseIndexes.bubbleSort();
        // Keep track of how many clauses deleted, because the indices are
        // are adjusted when we remove an element from HashArray
      deleted = 0;
      for (int i = 0; i < gndClauseIndexes.size(); i++)
      {
          // Atom appears pos. in these clauses, so remove the atom from
          // these clauses if false evidence, or remove clause if true evidence
          // or a unit clause
        if (trueEvidence ||
            (*gndClauses_)[gndClauseIndexes[i]]->getNumGroundPredicates() == 1)
        {
          if (vsdebug)
          {
            cout << "Deleting ground clause " << gndClauseIndexes[i] << " ";
            cout << endl;
          }
            // Real index is old index adjusted one lower for every element
            // deleted up until now
          delete (*gndClauses_)[gndClauseIndexes[i] - deleted];
          gndClauses_->removeItem(gndClauseIndexes[i] - deleted);
          deleted++;
        }
        else
        {
          if (vsdebug)
          {
            cout << "Removing gnd pred " << -(atomIdx + 1)
                 << " from ground clause " << gndClauseIndexes[i] << endl;
          }
          (*gndClauses_)[gndClauseIndexes[i]]->removeGndPred(atomIdx + 1);
        }
      }
      
      gndPredHashArray_.removeItemFastDisorder(atomIdx);
      gndPredHashArray_.compress();
      gndPreds_->removeItemFastDisorder(atomIdx);
      gndPreds_->compress();
        // By removing a pred, the pred at the end of the array gets the
        // index of the pred deleted, so we have to update to the new index
        // in all clauses
      int oldIdx = gndPredHashArray_.size();
      replaceAtomIndexInAllClauses(oldIdx, atomIdx);      
    }
  }

  /**
   * Sets a GroundPredicate to be query. If it is already present as query,
   * then nothing happens. If the predicate was evidence, then additional
   * clauses may be added. reinit() should be called after this in order to
   * ensure that the clause and atom information is correct.
   * 
   * @param predicate GroundPredicate to be set as a query.
   */
  void setAsQuery(const GroundPredicate* const & predicate)
  {
    if (vsdebug)
    {
      cout << "Setting to query " ;
      predicate->print(cout, domain_);
      cout << endl;
    }
    Database* db = domain_->getDB();
      // If already non-evidence, then do nothing
    if (gndPredHashArray_.contains((GroundPredicate*)predicate))
      return;
    else
    {
        // Evidence -> query
        // Add predicate to query set and get clauses
      gndPredHashArray_.append((GroundPredicate*)predicate);
      Predicate* p = predicate->createEquivalentPredicate(domain_);
      db->setEvidenceStatus(p, false);
      bool ignoreActivePreds = true;
      getActiveClauses(p, newClauses_, true, ignoreActivePreds);
    }
  }

  void printNetwork(ostream& out)
  {
    for (int i = 0; i < gndClauses_->size(); i++)
    {
      (*gndClauses_)[i]->print(out, domain_, &gndPredHashArray_);
      out << endl;
    }
  }
  
  
  ////////////// BEGIN: MCMC Functions //////////////
  
  /**
   * Gets a pointer to a GroundPredicate.
   * 
   * @param index Index of the GroundPredicate to be retrieved.
   * @return Pointer to the GroundPredicate 
   */
  GroundPredicate* getGndPred(const int& index)
  {
    return (*gndPreds_)[index];
  }

  /**
   * Gets a pointer to a GroundClause.
   * 
   * @param index Index of the GroundClause to be retrieved.
   * @return Pointer to the GroundClause 
   */
  GroundClause* getGndClause(const int& index)
  {
    return (*gndClauses_)[index];
  }

  /**
   * The atom assignment in the best state is saved to the ground predicates.
   */
  void saveLowStateToGndPreds()
  {
    for (int i = 0; i < getNumAtoms(); i++)
      (*gndPreds_)[i]->setTruthValue(lowAtom_[i + 1]);
  }

  /**
   * The atom assignment in the best state is saved to the database.
   */
  void saveLowStateToDB()
  {
    for (int i = 0; i < getNumAtoms(); i++)
    {
      GroundPredicate* p = gndPredHashArray_[i];
      bool value = lowAtom_[i + 1];
      if (value)
      {
        domain_->getDB()->setValue(p, TRUE);
      }
      else
      {
        domain_->getDB()->setValue(p, FALSE);
      }
    }
  }

  /**
   * Gets the index of a GroundPredicate in this state.
   * 
   * @param gndPred GroundPredicate whose index is being found.
   * @return Index of the GroundPredicate, if found; otherwise -1.
   */
  const int getGndPredIndex(GroundPredicate* const& gndPred)
  {
    return gndPreds_->find(gndPred);
  }

     
  ////////////// END: MCMC Functions //////////////


  ////////////// BEGIN: Lazy Functions //////////////
 
  /** 
   * Gets clauses and weights activated by the predicate inputPred,
   * if active is true. If false, inactive clauses (and their weights)
   * containing inputPred are retrieved. If inputPred is NULL, then all
   * active (or inactive) clauses and their weights are retrieved.
   * 
   * @param inputPred Only clauses containing this Predicate are looked at.
   * If NULL, then all active clauses are retrieved.
   * @param activeClauses New active clauses are put here.
   * @param active If true, active clauses are retrieved, otherwise inactive.
   * @param ignoreActivePreds If true, active preds are not taken into account.
   * This results in the retrieval of all unsatisfied clauses.
   */
  void getActiveClauses(Predicate *inputPred,
                        Array<GroundClause*>& activeClauses,
                        bool const & active,
                        bool const & ignoreActivePreds)
  {
    Timer timer;
    double currTime;

    Clause *fclause;
    GroundClause* newClause;
    int clauseCnt;
    GroundClauseHashArray clauseHashArray;

    Array<GroundClause*>* newClauses = new Array<GroundClause*>; 
  
    const Array<IndexClause*>* indexClauses = NULL;
      
      // inputPred is null: all active clauses should be retrieved
    if (inputPred == NULL)
    {
      clauseCnt = mln_->getNumClauses();
    }
      // Otherwise, look at all first order clauses containing the pred
    else
    {
      if (domain_->getDB()->getDeactivatedStatus(inputPred)) return;
      int predId = inputPred->getId();
      indexClauses = mln_->getClausesContainingPred(predId);
      clauseCnt = indexClauses->size();
    }

      // Look at each first-order clause and get active groundings
    int clauseno = 0;

    while (clauseno < clauseCnt)
    {
      if (inputPred)
        fclause = (Clause *) (*indexClauses)[clauseno]->clause;           
      else
        fclause = (Clause *) mln_->getClause(clauseno);

      if (vsdebug)
      {
        cout << "Getting active clauses for FO clause: ";
        fclause->print(cout, domain_);
        cout << endl;
      }

      currTime = timer.time();

      const double* parentWtPtr = NULL;
      if (!fclause->isHardClause()) parentWtPtr = fclause->getWtPtr();
      const int clauseId = mln_->findClauseIdx(fclause);
      newClauses->clear();

      if (stillActivating_)
        stillActivating_ = fclause->getActiveClauses(inputPred, domain_,
                                                     newClauses,
                                                     &gndPredHashArray_,
                                                     ignoreActivePreds);

      for (int i = 0; i < newClauses->size(); i++)
      {
        long double wt = fclause->getWt();
        newClause = (*newClauses)[i];

          // If already active, ignore; assume all gndClauses_ are active
        if (gndClauses_->find(newClause) >= 0)
        {
          delete newClause;
          continue;
        }

        bool invertWt = false;
          // We want to normalize soft unit clauses to all be positives
        if (!fclause->isHardClause() &&
            newClause->getNumGroundPredicates() == 1 &&
            !newClause->getGroundPredicateSense(0))
        {
          newClause->setGroundPredicateSense(0, true);
          newClause->setWt(-newClause->getWt());
          wt = -wt;
          invertWt = true;
          int addToIndex = gndClauses_->find(newClause);
          if (addToIndex >= 0)
          {
            (*gndClauses_)[addToIndex]->addWt(wt);
            if (parentWtPtr)
              (*gndClauses_)[addToIndex]->incrementClauseFrequency(clauseId, 1,
                                                                   invertWt);
            delete newClause;
            continue;
          }
        }
        
        int pos = clauseHashArray.find(newClause);
          // If clause already present, then just add weight
        if (pos >= 0)
        {
          if (vsdebug)
          {
            cout << "Adding weight " << wt << " to clause ";
            clauseHashArray[pos]->print(cout, domain_, &gndPredHashArray_);
            cout << endl;
          }
          clauseHashArray[pos]->addWt(wt);
          if (parentWtPtr)
            clauseHashArray[pos]->incrementClauseFrequency(clauseId, 1,
                                                           invertWt);
          delete newClause;
          continue;
        }

          // If here, then clause is not yet present        
        newClause->setWt(wt);
        if (parentWtPtr)
          newClause->incrementClauseFrequency(clauseId, 1, invertWt);

        if (vsdebug)
        {
          cout << "Appending clause ";
          newClause->print(cout, domain_, &gndPredHashArray_);
          cout << endl;
        }
        clauseHashArray.append(newClause);
      }
      clauseno++; 

    } //while (clauseno < clauseCnt)

    for (int i = 0; i < clauseHashArray.size(); i++)
    {
      newClause = clauseHashArray[i];
      activeClauses.append(newClause);
    }

    newClauses->clear();
    delete newClauses;

    clauseHashArray.clear();
  }

  /**
   * Get all the active clauses in the database.
   * 
   * @param allClauses Active clauses are retrieved into this Array.
   * @param ignoreActivePreds If true, active preds are ignored; this means
   * all unsatisfied clauses are retrieved.
   */
  void getActiveClauses(Array<GroundClause*> &allClauses,
                        bool const & ignoreActivePreds)
  {
    getActiveClauses(NULL, allClauses, true, ignoreActivePreds);
  }
  
  int getNumActiveAtoms()
  {
    return activeAtoms_; 
  }

  const void setBreakHardClauses(const bool& breakHardClauses)
  {
    breakHardClauses_ = breakHardClauses;
  }
  ////////////// END: Lazy Functions //////////////

  
 private:
 
  /**
   * Sets the hard clause weight to the sum of all soft clause weights in
   * the domain (not in the MRF) plus a constant (10). This is then the hard
   * clause weight when no soft clauses are present.
   */
  void setHardClauseWeight()
  {
      // Soft weights are summed up to determine hard weight
    long double sumSoftWts = 0.0;
      // Determine hard clause weight
    int clauseCnt = mln_->getNumClauses();    
      // Sum up the soft weights of all grounded clauses
    for (int i = 0; i < clauseCnt; i++)
    {
      Clause* fclause = (Clause *) mln_->getClause(i);
        // Skip hard clauses
      if (fclause->isHardClause()) continue;
        // Weight could be negative
      long double wt = abs(fclause->getWt());
      long double numGndings = fclause->getNumGroundings(domain_);
      sumSoftWts += wt*numGndings;
    }
    assert(sumSoftWts >= 0);
      // Add constant so weight isn't zero if no soft clauses present
    hardWt_ = sumSoftWts + 10.0;
    cout << "Set hard weight to " << hardWt_ << endl;
  }

  /**
   * Replaces the old index of an atom with a new one in all ground clauses
   * in which the old index occurs.
   * 
   * @param oldIdx Old index of atom.
   * @param newIdx New index of atom.
   */
  void replaceAtomIndexInAllClauses(const int& oldIdx, const int& newIdx)
  {
    Array<int> gndClauseIndexes;

    gndClauseIndexes = getNegOccurenceArray(oldIdx + 1);
    for (int i = 0; i < gndClauseIndexes.size(); i++)
    {
        // Atom appears neg. in these clauses, so it appears as -(oldIdx + 1)
      if ((*gndClauses_)[gndClauseIndexes[i]])
        (*gndClauses_)[gndClauseIndexes[i]]->changeGndPredIndex(-(oldIdx + 1),
                                                                -(newIdx + 1));
    }

    gndClauseIndexes = getPosOccurenceArray(oldIdx + 1);
    for (int i = 0; i < gndClauseIndexes.size(); i++)
    {
        // Atom appears pos. in these clauses, so it appears as (oldIdx + 1)
      if ((*gndClauses_)[gndClauseIndexes[i]])
        (*gndClauses_)[gndClauseIndexes[i]]->changeGndPredIndex(oldIdx + 1,
                                                                newIdx + 1);
    }
  }

 public:
   // Inference mode: decide how getActiveClauses and addNewClauses
   // filter and set cost
  const static int MODE_MWS = 0;
  const static int MODE_HARD = 1;
  const static int MODE_SAMPLESAT = 2;

  const static int ADD_CLAUSE_INITIAL = 0;
  const static int ADD_CLAUSE_REGULAR = 1;
  const static int ADD_CLAUSE_DEAD = 2;
  const static int ADD_CLAUSE_SAT = 3;    // the clauses are sat by fixed atoms

 private:
    // If true, this is a lazy variable state, else eager.
  bool lazy_;

    // Weight used for hard clauses (sum of soft weights + constant)
  long double hardWt_;
  
    // mln and domain are used to build MRF in eager state and to
    // retrieve active atoms in lazy state.
  MLN* mln_;
  Domain* domain_;

    // Eager version: Pointer to gndPreds_ and gndClauses_ in MRF
    // Lazy version: Holds active atoms and clauses
  Array<GroundPredicate*>* gndPreds_;
  
  //Array<GroundClause*>* gndClauses_;
  GroundClauseHashArray* gndClauses_;
  
    // Predicates corresponding to the groundings of the unknown non-evidence
    // predicates
  GroundPredicateHashArray* unePreds_;

    // Predicates corresponding to the groundings of the known non-evidence
    // predicates
  GroundPredicateHashArray* knePreds_;
    // Actual truth values of ground known non-evidence preds
  Array<TruthValue>* knePredValues_;  
  

    ////////// BEGIN: Information used only by lazy version //////////
    // Number of distinct atoms in the first set of unsatisfied clauses
  int baseNumAtoms_;
    // If true, atoms are not deactivated when mem. is full
  bool noApprox_;
    // Indicates whether deactivation of atoms has taken place yet
  bool haveDeactivated_;
    // Max. amount of memory to use
  int memLimit_;
    // Max. amount of clauses memory can hold
  int clauseLimit_;
    ////////// END: Information used only by lazy version //////////


    ////////// BEGIN: Information used only by eager version //////////
    // MRF is used with eager states. If lazy, this stays NULL.
  MRF* mrf_;
    ////////// END: Information used only by eager version //////////

    // Holds the new active clauses
  Array<GroundClause*> newClauses_;
    // Holds the new gnd preds
  Array<GroundPredicate*> newPreds_;
    // Holds the ground predicates in a hash array.
    // Fast access is needed for comparing preds when activating clauses.
  GroundPredicateHashArray gndPredHashArray_;

    // Clauses to be satisfied
    // Indexed as clause_[clause_num][literal_num]
  Array<Array<int> > clause_;
    // Cost of each clause (can be negative)
  Array<long double> clauseCost_;
    // Highest cost of false clause
  long double highestCost_;
    // If true, more than one clause has highest cost
  bool eqHighest_;
    // Number of clauses with highest cost
  int numHighest_;
    // Clauses which are pos. and unsatisfied or neg. and satisfied
  Array<int> falseClause_;
    // Where each clause is listed in falseClause_
  Array<int> whereFalse_;
    // Number of true literals in each clause
  Array<int> numTrueLits_;
    // watch1_[c] contains the id of the first atom which c is watching
  Array<int> watch1_;
    // watch2_[c] contains the id of the second atom which c is watching
  Array<int> watch2_;
    // Which clauses are satisfied by fixed atoms
    // wt_sat (pos-wt & sat || neg-wt & unsat)
  Array<bool> isSatisfied_;
    // Clauses which are not to be considered
  Array<bool> deadClause_;
    // Use threshold to exclude clauses from the state?
  bool useThreshold_;
    // Pre-computed thresholds for each clause
  Array<long double> threshold_;

    // Holds the index of clauses in which each literal occurs
    // Indexed as occurence_[2*abs(lit) - (lit > 0)][occurence_num]
  Array<Array<int> > occurence_;

    // Current assigment of atoms
  Array<bool> atom_;
  //Array<bool> atomEvi_;
    // Cost of clauses which would become satisfied by flipping each atom
  Array<long double> makeCost_;
    // Cost of clauses which would become unsatisfied by flipping each atom
  Array<long double> breakCost_;
    // Indicates if an atom is fixed to a value (0 = not fixed, -1 = false,
    // 1 = true)
  Array<int> fixedAtom_;

    // Assigment of atoms producing lowest cost so far
  Array<bool> lowAtom_;
    // Cost of false clauses in the currently best state
  long double lowCost_;
    // Number of false clauses in the currently best state
  int lowBad_;

    // Current no. of unsatisfied clauses
  int numFalseClauses_;
    // Cost associated with the number of false clauses
  long double costOfFalseClauses_;
  
    // Number of active atoms in state.  
  int activeAtoms_;

    // For MC-SAT: True if clause satisfied in last iteration
  Array<bool> prevSatisfiedClause_;
  
    // Number of non-evidence atoms in the domain
  int numNonEvAtoms_;

    // Inference mode: decide how getActiveClauses and addNewClauses
    // filter and set cost
  int inferenceMode_;
  
    // Indicates if atoms are still being activated. If main memory is full,
    // then activating atoms is disabled and inference is run in the partial
    // network
  bool stillActivating_;
  
    // If true, hard clauses can be broken
  bool breakHardClauses_;

};

#endif /*VARIABLESTATE_H_*/
