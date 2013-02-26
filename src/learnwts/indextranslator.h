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
#ifndef INDEXTRANSLATOR_H_NOV_29_2005
#define INDEXTRANSLATOR_H_NOV_29_2005


#include "mln.h"

//NOTE: "domain index" and "database index" are used interchangeably

struct IdxDiv { int idx; double div; };


class IndexTranslator
{
 public:
  IndexTranslator(const Array<MLN*>* const & mlns, 
                  const Array<Domain*>* const & domains) 
    : mlns_(mlns), domains_(domains), cIdxToCFIdxsPerDomain_(NULL), 
      wtsPerDomain_(NULL), gradsPerDomain_(NULL),
      clauseOrdering_(NULL), existFormulaOrdering_(NULL) 
  { createClauseIdxToClauseFormulaIdxsMap(); }


  ~IndexTranslator() 
  {
    if (cIdxToCFIdxsPerDomain_) 
    {
      for (int i = 0; i < cIdxToCFIdxsPerDomain_->size(); i++)
        (*cIdxToCFIdxsPerDomain_)[i].deleteItemsAndClear();
      delete cIdxToCFIdxsPerDomain_;
    }
    if (wtsPerDomain_) delete wtsPerDomain_;
    if (gradsPerDomain_) delete gradsPerDomain_;
    if (clauseOrdering_) 
    { clauseOrdering_->deleteItemsAndClear(); delete clauseOrdering_; }
    if (existFormulaOrdering_) delete existFormulaOrdering_;
  }    


    //Returns true if there are multiple databases, and there are formulas
    //whose CNF contains a different number of clauses for each database;
    //otherwise returns false. 
    //ASSUMPTION: the MLNs were formed from parsing the same .mln file
  static bool needIndexTranslator(const Array<MLN*>& mlns,
                                  const Array<Domain*>& domains)
  {
    if (mlns.size() <= 1) return false;

    for (int i = 1; i < mlns.size(); i++)
      if (mlns[0]->getNumClauses() != mlns[i]->getNumClauses()) return true;

    //there are more than one MLN and all MLNs have the same number of clauses
    //now check whether all tie-able formulas have one clause in their CNF

    bool allTiedFormulasHaveOneClause = true;
    for (int i = 0; i < mlns.size(); i++)
    {
      const FormulaAndClausesArray* fnca = mlns[i]->getFormulaAndClausesArray();
      for (int j = 0; j < fnca->size(); j++)
      {
        if (formulaWtTiedToClauseWts((*fnca)[j]) && 
            (*fnca)[j]->indexClauses->size() > 1)
        { 
          allTiedFormulasHaveOneClause = false;
          break;
        }
      }
      if (!allTiedFormulasHaveOneClause) break;
    }

    if (allTiedFormulasHaveOneClause) return false;

    //all MLNs have the same number of clauses, and a tie-able formula
    //has more than one clause
    
    //now check that the tie-able clauses are exactly the same across databases

    for (int i = 0; i < mlns[0]->getNumClauses(); i++)
    {
        //if this is a clause whose wt can be tied to the formula wt
      if (tieClauseWtToFormulaWt(mlns[0], i))
      {
        ostringstream oss0;
        mlns[0]->getClause(i)->printWithoutWtWithStrVar(oss0, domains[0]);
        for (int j = 1; j < mlns.size(); j++)
        {
          ostringstream ossj;
          mlns[j]->getClause(i)->printWithoutWtWithStrVar(ossj, domains[j]);
          if (oss0.str().compare(ossj.str()) != 0) return true;          
        }
      }
    }
    
    return false;
  }


  const Array<Array<Array<IdxDiv>*> >* 
  getClauseIdxToClauseFormulaIdxsPerDomain() const 
  { return cIdxToCFIdxsPerDomain_; }

  Array<Array<double> >* getWtsPerDomain() const { return wtsPerDomain_; }

  Array<Array<double> >* getGradsPerDomain() const { return gradsPerDomain_; }

  const ClauseHashArray* getClauseOrdering() const { return clauseOrdering_; }

  const StringHashArray* getExistFormulaOrdering() const 
  { return existFormulaOrdering_; }

  int getNumClausesAndExistFormulas() const 
  { return clauseOrdering_->size() + existFormulaOrdering_->size(); }

  
    //check that the sizes of cIdxToCFIdxsPerDomain_[d], wtsPerDomain_[d], and
    //gradsPerDomain_[d] is equal to the num of clauses in mln_[d] + numAdded
  bool checkCIdxWtsGradsSize(const int& numAdded)
  {
    for (int i = 0; i < mlns_->size(); i++) // for each domain
    {
      int sz = (*mlns_)[i]->getNumClauses() + numAdded;
      if ((*cIdxToCFIdxsPerDomain_)[i].size() != sz) return false;
      if ((*wtsPerDomain_)[i].size() != sz) return false;
      if ((*gradsPerDomain_)[i].size() != sz) return false;
    }
    return true;
  }


    //Size of wts must be sum of the sizes of clauseOrdering_, 
    //existFormulaOrdering_, appendedClauses, and appendedFormulas.
  void updateClauseFormulaWtsInMLNs(const Array<double>& wts,
                                    const Array<Clause*>*const& appendedClauses,
                                    const Array<string>*const& appendedFormulas)
  {
    assert(appendedClauses == NULL || appendedFormulas == NULL);
    Domain* dom0 = (*domains_)[0];
    for (int i = 0; i < mlns_->size(); i++) // for each domain
    {
      MLN* mln = (*mlns_)[i];
      Domain* dom = (*domains_)[i];

      for (int j = 0; j < clauseOrdering_->size(); j++)
      {
        Clause* c = (*clauseOrdering_)[j];
        if (i > 0 && c->containsConstants()) c->translateConstants(dom0, dom);
        Clause* cc = (Clause*) mln->findClause(c);
        if (i > 0 && c->containsConstants()) c->translateConstants(dom, dom0);
        cc->setWt(wts[j]);
      }

      int offset = clauseOrdering_->size();

      for (int j = 0; j < existFormulaOrdering_->size(); j++)
        mln->setFormulaWt((*existFormulaOrdering_)[j], wts[offset+j]);

      offset += existFormulaOrdering_->size();

      if (appendedClauses)
      {
        for (int j = 0; j < appendedClauses->size(); j++)
        {
          Clause* c = (*appendedClauses)[j];
          if (i > 0 && c->containsConstants()) c->translateConstants(dom0, dom);
          Clause* cc = (Clause*) mln->findClause(c);
          if (i > 0 && c->containsConstants()) c->translateConstants(dom, dom0);
          cc->setWt(wts[offset+j]);
        }
      }
      else
      if (appendedFormulas)
      {
        for (int j = 0; j < appendedFormulas->size(); j++)
          mln->setFormulaWt((*appendedFormulas)[j], wts[offset+j]);
      }       
    }
  }


  void getClauseFormulaWts(Array<double>& wts) const
  {
    wts.clear();
    wts.growToSize(clauseOrdering_->size() + existFormulaOrdering_->size());
    int a = 0;
    for (int j = 0; j < clauseOrdering_->size(); j++)
      wts[a++] = (*mlns_)[0]->findClause((*clauseOrdering_)[j])->getWt();

    for (int j = 0; j < existFormulaOrdering_->size(); j++)
      wts[a++] = (*mlns_)[0]->getFormulaWt((*existFormulaOrdering_)[j]);    
  }


    //There must be more than one MLN in mlns. The size of mlns_ must be the 
    //same as the size of domains_.
  void createClauseIdxToClauseFormulaIdxsMap()
  {
    assert(mlns_->size() == domains_->size());
    if (mlns_->size() <= 1) return;
    if (cIdxToCFIdxsPerDomain_) 
    {
      for (int i = 0; i < cIdxToCFIdxsPerDomain_->size(); i++)
        (*cIdxToCFIdxsPerDomain_)[i].deleteItemsAndClear();
      delete cIdxToCFIdxsPerDomain_;
    }
    if (wtsPerDomain_) delete wtsPerDomain_;
    if (gradsPerDomain_) delete gradsPerDomain_;
    if (clauseOrdering_)     
    { clauseOrdering_->deleteItemsAndClear(); delete clauseOrdering_; }
    if (existFormulaOrdering_) delete existFormulaOrdering_;

    trackClauseConstants();
    clauseOrdering_ = new ClauseHashArray;
    existFormulaOrdering_ = new StringHashArray;
    createClausesFormulasOrdering((*mlns_)[0], clauseOrdering_, 
                                  existFormulaOrdering_);

    cIdxToCFIdxsPerDomain_ = new Array<Array<Array<IdxDiv>*> >;
    wtsPerDomain_ = new Array<Array<double> >;
    gradsPerDomain_ = new Array<Array<double> >;
    cIdxToCFIdxsPerDomain_->growToSize(mlns_->size());
    wtsPerDomain_->growToSize(mlns_->size());
    gradsPerDomain_->growToSize(mlns_->size());

    Domain* dom0 = (*domains_)[0];
    for (int i = 0; i < mlns_->size(); i++) // for each domain
    {
      MLN* mln = (*mlns_)[i];
      Domain* dom = (*domains_)[i];
      const FormulaAndClausesArray* fnca = mln->getFormulaAndClausesArray();
      Array<Array<IdxDiv>*>& cIdxToCFIdxsMap = (*cIdxToCFIdxsPerDomain_)[i];
      cIdxToCFIdxsMap.growToSize(mln->getNumClauses());
      (*wtsPerDomain_)[i].growToSize(mln->getNumClauses());
      (*gradsPerDomain_)[i].growToSize(mln->getNumClauses());

      for (int j = 0; j < mln->getNumClauses(); j++) // for each clause in MLN
      {
        Array<IdxDiv>* cfIdxs = new Array<IdxDiv>;
        cIdxToCFIdxsMap[j] = cfIdxs;

          //if we don't need to tie the clause weight to formula weight as when
          //this clause is in the CNF of a non-existentially quant. and
          //non-existentially and uniquely quant. formula.
        if (noTieClauseWtToFormulaWt(mln,j))
        {
          Clause* c = (Clause*) mln->getClause(j);
          if (i > 0 && c->containsConstants()) c->translateConstants(dom, dom0);
          int idx = clauseOrdering_->find(c);
          cfIdxs->append(IdxDiv());
          cfIdxs->lastItem().idx = idx;
          cfIdxs->lastItem().div = 1;
          if (i > 0 && c->containsConstants()) c->translateConstants(dom0, dom);
          if (idx < 0) 
          { 
            cout <<"ERROR: in IndexTranslator::"
                 << "createClauseIdxToClauseFormulaIdxsMap(): clause ";
            c->printWithoutWtWithStrVar(cout, dom); cout << "not found!"<<endl;
            exit(-1);
          }
        }

          //if we must tie the clause weight to formula weight as when
          //the clause is in the CNF of an existentially quantified formula,
          //or an existentially and uniquely quantified formula,
        if (tieClauseWtToFormulaWt(mln, j))
        {
          const Array<FormulaClauseIndexes*>& fciArr 
            = mln->getMLNClauseInfo(j)->formulaClauseIndexes;

          for (int k = 0; k < fciArr.size(); k++)
          {
            int fidx = *(fciArr[k]->formulaIndex);
            if (!formulaWtTiedToClauseWts((*fnca)[fidx])) continue;
            string formula = (*fnca)[fidx]->formula;
            int numClauses = (*fnca)[fidx]->indexClauses->size();
            int idx = existFormulaOrdering_->find(formula);
            cfIdxs->append(IdxDiv());
                //exist. quant. formulas comes after 'normal' clauses
            cfIdxs->lastItem().idx = idx + clauseOrdering_->size();

            if (mln->isExistClause(j))
              cfIdxs->lastItem().div = numClauses;
            else
            if (mln->isExistUniqueClause(j))
            {
                //this clause states that there is at least one true value
              Clause* firstClause 
                = (*((*fnca)[fidx]->indexClauses))[0]->clause;
              int numPreds = firstClause->getNumPredicates();
              Clause* c = (Clause*) mln->getClause(j);
              if (c->same(firstClause))
                cfIdxs->lastItem().div = (numPreds > 1) ? numPreds+1 : 1;
              else
                cfIdxs->lastItem().div = (numPreds*numPreds-1)*0.5;
            }
            else
              assert(false);


            if (idx < 0)
            { 
              cout << "ERROR: in createClauseIdxToClauseFormulaIdxsMap(): "
                   << "formula " << formula << "not found!" << endl;
              exit(-1);
            }
          }
        }

        cfIdxs->compress();

      } // for each clause in MLN
    } // for each DB
  }


  void appendClauseIdxToClauseFormulaIdxs(const int& numCFIdxs,
                                          const int& numCIdxsPerCFIdx,
                                          const int& domainIdx)
  {
    if (numCFIdxs == 0) return;
    Array<Array<IdxDiv>*>& cidxs = (*cIdxToCFIdxsPerDomain_)[domainIdx];
    Array<double>& wts = (*wtsPerDomain_)[domainIdx];
    Array<double>& grads = (*gradsPerDomain_)[domainIdx];
    int numClausesFormulas = getNumClausesAndExistFormulas();
    
    for (int j = 0; j < numCFIdxs; j++)
    {
      int cfidx = numClausesFormulas + j;
      for (int k = 0; k < numCIdxsPerCFIdx; k++)
      {
        IdxDiv idiv;
        idiv.idx = cfidx;
        idiv.div = 1;
        Array<IdxDiv>* arr = new Array<IdxDiv>(1);
        arr->append(idiv);
        cidxs.append(arr);
        
        wts.append(0.0);
        grads.append(0.0);          
      }
    }
  }

  
  void appendClauseIdxToClauseFormulaIdxs(const int& numCFIdxs,
                                          const int& numCIdxsPerCFIdx)
  {
    if (numCFIdxs == 0) return;
    for (int i = 0; i < cIdxToCFIdxsPerDomain_->size(); i++)
      appendClauseIdxToClauseFormulaIdxs(numCFIdxs, numCIdxsPerCFIdx, i);
  }


  void removeClauseIdxToClauseFormulaIdxs(const int& numCFIdxs,
                                          const int& numCIdxsPerCFIdx,
                                          const int& domainIdx)
  {
    if (numCFIdxs == 0) return;
    Array<Array<IdxDiv>*>& cidxs = (*cIdxToCFIdxsPerDomain_)[domainIdx];
    Array<double>& wts = (*wtsPerDomain_)[domainIdx];
    Array<double>& grads = (*gradsPerDomain_)[domainIdx];
    
    for (int j = 0; j < numCFIdxs; j++)
      for (int k = 0; k < numCIdxsPerCFIdx; k++)
      {
        delete cidxs.removeLastItem();
        wts.removeLastItem();
        grads.removeLastItem();          
      }    
  }


  void removeClauseIdxToClauseFormulaIdxs(const int& numCFIdxs,
                                          const int& numCIdxsPerCFIdx)
  {
    if (numCFIdxs == 0) return;
    for (int i = 0; i < cIdxToCFIdxsPerDomain_->size(); i++)
      removeClauseIdxToClauseFormulaIdxs(numCFIdxs, numCIdxsPerCFIdx, i);
  }


  void setPriorMeans(Array<double>& priorMeans)
  {
    priorMeans.clear();
    priorMeans.growToSize(clauseOrdering_->size() + 
                          existFormulaOrdering_->size());
    int m = 0;
    MLN* mln0 = (*mlns_)[0];

    for (int i = 0; i < clauseOrdering_->size(); i++)
    {
      Clause* c = (*clauseOrdering_)[i]; 
      int cidx = mln0->findClauseIdx(c);      
      priorMeans[m++] = mln0->getMLNClauseInfo(cidx)->priorMean;
    }
    
    const FormulaAndClausesArray* fnca = mln0->getFormulaAndClausesArray();
    for (int i = 0; i < existFormulaOrdering_->size(); i++)
    {
      FormulaAndClauses tmp((*existFormulaOrdering_)[i], 0, false, false,
                            false);
      int a = fnca->find(&tmp);
      assert(a >= 0);
      priorMeans[m++] = (*fnca)[a]->priorMean;
    }
  }


    //Maps each clause index in the array to the clause/formula index, and
    //replaces the former with the latter.
  void getClauseFormulaIndexes(Array<int>& clauseIndexes, const int& domainIdx)
  {
    for (int i = 0; i < clauseIndexes.size(); i++)
    {
      int cidx = clauseIndexes[i];
      Array<IdxDiv>* idxDivs = (*cIdxToCFIdxsPerDomain_)[domainIdx][cidx];
      for (int j = 0; j < idxDivs->size(); j++)
        if ((*idxDivs)[j].idx < clauseOrdering_->size())
        {
            //only get the index of clauses (and not those of exist.
            //quant. formulas)
          assert((*idxDivs)[j].div == 1);
          clauseIndexes[i] = (*idxDivs)[j].idx;
        }
    }
  }


    //assign the weights belonging to clauses (and none of those belonging to
    //existentially/exitentially and uniquely  quantified formulas) to the MLNs
  void assignNonTiedClauseWtsToMLNs(const double* const & wts)
  {
    for (int i = 0; i < mlns_->size(); i++)
    {
      MLN* mln = (*mlns_)[i];
      for (int j = 0; j < mln->getNumClauses(); j++)
      {
        if (!noTieClauseWtToFormulaWt(mln,j)) continue;

        Clause* c = ((Clause*) mln->getClause(j));
        assert(inClauseOrdering(c, i));
        Array<IdxDiv>* idxDivs =(*cIdxToCFIdxsPerDomain_)[i][j];
        double wt = 0;
        for (int k = 0; k < idxDivs->size(); k++)
        {
            //don't assign the weights of existential formulas to clauses
          if ((*idxDivs)[k].idx < clauseOrdering_->size())
            wt += wts[ (*idxDivs)[k].idx ] / (*idxDivs)[k].div;
        }
        c->setWt(wt);
      }
    }
  }

  
  void printClauseFormulaWts(ostream& out, const bool& includeIdx)
  {
    Array<double> wts;
    getClauseFormulaWts(wts);
      //print those clauses that appear in a non-exist. quant. formula
    Array<Clause*> tmp;
    for (int i = 0; i < clauseOrdering_->size(); i++) 
    {
      (*clauseOrdering_)[i]->setWt(wts[i]);
      tmp.append((*clauseOrdering_)[i]);
    }
    Clause::sortByLen(tmp);
      
    out.setf(ios_base::left, ios_base::adjustfield);
    for (int i = 0; i < tmp.size(); i++)
    {
      if (includeIdx) { out << i << ":  "; out.width(14); }
      else            { out.width(10); }
      out << tmp[i]->getWt(); out.width(0); out << " ";
      tmp[i]->printWithoutWtWithStrVar(out, (*domains_)[0]); out << endl; 
    }
    out.width(0);

    for (int i = 0; i < existFormulaOrdering_->size(); i++)
    {
      int idx = i+clauseOrdering_->size();
      if (includeIdx) { out << idx <<":  ";out.width(14);}
      else            { out.width(10); }
      out << wts[idx]; out.width(0); 
      out << " " << (*existFormulaOrdering_)[i] << endl;      
    }
    out.width(0);
  }


    //used by VotedPerceptron
  void setRelevantClausesFormulas(Array<bool>& relevantClausesFormulas,
                                  const Array<bool>& relevantClausesInMLN0)
  {
    const Array<bool>& relevantClauses = relevantClausesInMLN0;

    if (relevantClausesFormulas.size() < relevantClauses.size())
      relevantClausesFormulas.growToSize(relevantClauses.size());
    else
    if (relevantClausesFormulas.size() > relevantClauses.size())
      relevantClausesFormulas.shrinkToSize(relevantClauses.size());
      
    for (int i = 0; i < relevantClausesFormulas.size(); i++)
      relevantClausesFormulas[i] = false;

     
    for (int i = 0; i < clauseOrdering_->size(); i++)
    {
      int cidx = (*mlns_)[0]->findClauseIdx((*clauseOrdering_)[i]);
      relevantClausesFormulas[i] = relevantClauses[cidx];
    }

    int offset = clauseOrdering_->size();      
      
    for (int i = 0; i < existFormulaOrdering_->size(); i++)
    {
      const IndexClauseHashArray* icha 
        = (*mlns_)[0]->getClausesOfFormula((*existFormulaOrdering_)[i]);
      for (int j = 0; j < icha->size(); j++)
      {
        int cidx = (*mlns_)[0]->findClauseIdx((*icha)[j]->clause);
        relevantClausesFormulas[offset+i] = relevantClauses[cidx];
        if (relevantClausesFormulas[offset+i]) break;
      }
    }
  }


  void printRelevantClausesFormulas(ostream& out,
                                    const Array<bool>& relevantClausesFormulas)
  {
    for (int i = 0; i < clauseOrdering_->size(); i++)
    {
      if (relevantClausesFormulas[i])
      {
        out << i << ": ";
        (*clauseOrdering_)[i]->printWithoutWtWithStrVar(out, (*domains_)[0]);
        out << endl;
      }
    }
    
    int offset = clauseOrdering_->size();

    for (int i = 0; i < existFormulaOrdering_->size(); i++)
    {
      if (relevantClausesFormulas[offset+i])
        out << offset+i << ": " << (*existFormulaOrdering_)[i] << endl;
    }
  }


 private:
  void trackClauseConstants()
  {
    for (int i = 0; i < mlns_->size(); i++)
      for (int j = 0; j < (*mlns_)[i]->getNumClauses(); j++)
        if (noTieClauseWtToFormulaWt((*mlns_)[i], j))
        {
          Clause* c = (Clause*) (*mlns_)[i]->getClause(j);
          if (c->getAuxClauseData() == NULL) c->newAuxClauseData();
          c->trackConstants();
        }
  }


  void createClausesFormulasOrdering(MLN* const & mln, 
                                     ClauseHashArray* const & clauses,
                                     StringHashArray* const & existFormulas)
  {
    for (int i = 0; i < mln->getNumClauses(); i++)
      if (noTieClauseWtToFormulaWt(mln,i))
      { 
        Clause* c = new Clause( *((Clause*)mln->getClause(i)) );
        clauses->append(c);
      }

    const FormulaAndClausesArray* fnca = mln->getFormulaAndClausesArray();
    for (int i = 0; i < fnca->size(); i++)
      if (formulaWtTiedToClauseWts((*fnca)[i])) 
        existFormulas->append((*fnca)[i]->formula);
  }


  bool inClauseOrdering(Clause* const & c, const int& domainIdx)
  {
    assert(c->getAuxClauseData());
    if (domainIdx > 0 && c->containsConstants()) 
      c->translateConstants((*domains_)[domainIdx], (*domains_)[0]);
    bool in = (clauseOrdering_->find(c) >= 0);
    if (domainIdx > 0 && c->containsConstants()) 
      c->translateConstants((*domains_)[0], (*domains_)[domainIdx]);
    return in;
  }
  

    //Returns true if the clause weights must be tied to the formulas (as when
    //the clause belongs to the CNF of an existentially quantified or 
    //existentially and uniquely quantified formula
    //i is the clause's index in mln
  static bool tieClauseWtToFormulaWt(const MLN* const & mln, const int& i)
  { return (mln->isExistClause(i) || mln->isExistUniqueClause(i)); }


  static bool formulaWtTiedToClauseWts(const FormulaAndClauses* const & fnc)
  { return (fnc->hasExist || fnc->isExistUnique); }


    //Returns true if the clause weights need not be tied to the formulas 
    //(as when the clause belongs to the CNF of a non-existentially quantified 
    //and non-existentially and uniquely quantified formula
    //i is the clause's index in mln
  bool noTieClauseWtToFormulaWt(const MLN* const & mln, const int& i)
  { return mln->clauseInNonExistAndNonExistUniqueFormulaCNF(i); }

 private:
  const Array<MLN*>* mlns_; // not owned by this object; do not delete;
  const Array<Domain*>* domains_; // not owned by this object; do not delete;

    //cIdxToCFIdxsPerDomain_[d] maps clause indexes (cIdx) in an MLN to 
    //clause/formula indexes (CFIdxs) for domain d.
    //Used when the clauses across multiple domains do not line up perfectly
    //as when an existential formula has a different number of CNF clauses for
    //different domains. Used in PseudoLogLikelihood.
  Array<Array<Array<IdxDiv>*> >* cIdxToCFIdxsPerDomain_;
  Array<Array<double> >* wtsPerDomain_; //Used in PseudoLogLikelihood
  Array<Array<double> >* gradsPerDomain_; //Used in PseudoLogLikelihood

  ClauseHashArray* clauseOrdering_; 
  StringHashArray* existFormulaOrdering_;
};


#endif
