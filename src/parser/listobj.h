/*
 * All of the documentation and software included in the
 * Alchemy Software is copyrighted by Stanley Kok, Parag
 * Singla, Matthew Richardson, Pedro Domingos, Marc
 * Sumner, Hoifung Poon, and Daniel Lowd.
 * 
 * Copyright [2004-07] Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, and Daniel Lowd. All rights reserved.
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
 * Poon, and Daniel Lowd in the Department of Computer Science and
 * Engineering at the University of Washington".
 * 
 * 4. Your publications acknowledge the use or
 * contribution made by the Software to your research
 * using the following citation(s): 
 * Stanley Kok, Parag Singla, Matthew Richardson and
 * Pedro Domingos (2005). "The Alchemy System for
 * Statistical Relational AI", Technical Report,
 * Department of Computer Science and Engineering,
 * University of Washington, Seattle, WA.
 * http://www.cs.washington.edu/ai/alchemy.
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
#ifndef LISTOBJ_H_JUN_21_2005
#define LISTOBJ_H_JUN_21_2005

#include <cstdio>
#include <string>
#include <list>
#include <ostream>
#include <ext/hash_map>
using namespace __gnu_cxx;
#include "array.h"
#include "hashstring.h"
#include "domain.h"
#include "arraysaccessor.h"

typedef hash_map<const string, int, HashString, EqualString> VarTypeMap;


// Modified from http://aima.cs.berkeley.edu/lisp/doc/overview-LOGIC.html 
// Represent a first-order formula in prefix form e.g. (AND (OR (P x)) (Q y)).
// Egs of expression in prefix form (! P), (^ P Q), (v P Q), (=> P Q), 
// (<=> P Q), (P x) (a predicate or function), (FORALL (x) (P x)), 
// (EXIST (x) (P x)) 

class ListObj
{
 public:
  ListObj() 
  {
    str_ = new char[1];
    str_[0] = '\0';
  }
    
  
  ListObj(const char* const & str)
  {
    if (str != NULL)
    {
      str_ = new char[strlen(str)+1];
      strcpy(str_, str);
    }
    else
    {
      str_ = new char[1];
      str_[0] = '\0';
    }
  }


  ListObj(const ListObj* const & p)
  {
    for (int i = 0; i < p->size(); i++)
      list_.append(new ListObj((*p)[i]));   
    str_ = new char[strlen(p->str_)+1];
    strcpy(str_, p->str_);
       
    assert(p->list_.size() == list_.size());

  }


  ~ListObj() { delete [] str_; list_.deleteItemsAndClear(); }


  void clear() { list_.clear(); delete [] str_; str_=new char[1];str_[0]='\0'; }

  void deleteAndClear() { list_.deleteItemsAndClear(); clear(); }

  void clearList() { list_.clear(); }

  ostream& print(ostream& out) const
  {
    if (isStr()) 
      out << str_;
    else
    {
      out << "(";
      for (int i = 0; i < list_.size(); i++)
      {
        list_[i]->print(out); 
        out << ((i<list_.size()-1)? " " : "");
      }
      out << ")";
    }
    return out;
  }


    // caller should not modify the returned Array
  const Array<ListObj*>& getList() const { return list_; }

    // caller should not modify or delete the returned string
  const char* getStr() const { return str_; }

    // if isStr(), replace str_ with the provided str
    // caller is responsible for deleting str if necessary
  void setStrIfStr(const char* str)
  {
    if (isStr())
    {
      delete [] str_;
      str_ = new char[strlen(str)+1];
      strcpy(str_,str);
    }
  }

  bool isStr() const { return (strlen(str_) > 0); }

  bool isList() const { return !isStr(); }

  void append(ListObj* const & p) { list_.append(p); }  

  void append(const char* const & s) { list_.append(new ListObj(s)); }

  ListObj* removeLast()
  {
    ListObj* lo = list_.removeLastItem();
    list_.compress();
    return lo;
  }

  int size() const { return list_.size(); }

  ListObj* operator[](const int& index) const { return list_[index]; } 


  static ListObj* toCNF(const ListObj* const & p, const ListObj* const & vars, 
                        const Domain* const & domain, 
                        const VarTypeMap * const& vtMap, bool& hasExist)
  {
    ListObj* p1 = eliminateImplication(p);

    const char* opp = op(p1)->str_;
    
    if (strcmp(opp, "!")==0)
    {
      ListObj* p2 = moveNotInwards(arg1(p1));
      if (isLiteralClause(p2)) { delete p1; return p2; }

      ListObj* ret = toCNF(p2, vars, domain, vtMap, hasExist);
      delete p2; delete p1;
      return ret;
    }
    
    if (strcmp(opp, "^")==0)
    {
      ListObj* n = NULL;
      int numConj = 0;
      ListObj* firstConj = NULL;

      for (int i = 1; i < p1->size(); i++) // for each arg
      {        
        ListObj* conj = toCNF((*p1)[i], vars, domain, vtMap, hasExist);

        bool isConj = strcmp(op(conj)->str_, "^")==0;
        int jbeg = (isConj) ? 1 : 0;
        int jmax = (isConj) ? conj->size() : 1;
        ListObj* conjj;

        for (int j = jbeg; j < jmax; j++) 
        {
          if (isConj) conjj = (*conj)[j];
          else        conjj = conj;

          numConj++;
          if (numConj==1) { firstConj = conjj; continue; }
          else 
          if (numConj==2) {n=new ListObj; n->append("^"); n->append(firstConj);}
          n->append(conjj);
        }

        if (isConj) { delete (*conj)[0]; conj->clearList(); delete conj; }

      }

      delete p1;
      if (numConj == 1) { return firstConj; }        
      return n;
    }
    
    if (strcmp(opp, "v")==0)
    {
      ListObj* n = new ListObj;
      for (int i = 1; i < p1->size(); i++) // for each arg
        n->append(toCNF((*p1)[i], vars, domain, vtMap, hasExist));
      ListObj* ret = mergeDisjuncts(n, n->size());

      delete n; delete p1;
      return ret;
    }
    
    if (strcmp(opp, "FORALL")==0)
    {
      const ListObj* arg1var = arg1(p1);

      ListObj* newVar = new ListObj;
      for (int i = 0; i < arg1var->size(); i++)
        newVar->append(newVariable((*arg1var)[i]));
      assert(newVar->size() == arg1var->size());
      
      const ListObj* aarg2 = new ListObj(arg2(p1));

      for (int i = 0; i < newVar->size(); i++)
        for (int j = 0; j < aarg2->size(); j++)
          (*aarg2)[j]->sub((*arg1var)[i], (*newVar)[i]);

      ListObj* app = new ListObj;
      for (int i = 0; i < newVar->size(); i++)  app->append((*newVar)[i]);
      for (int i = 0; i < vars->size(); i++)    app->append((*vars)[i]);
      
      newVar->clearList();
      delete newVar;
      ListObj* ret = toCNF(aarg2, app, domain, vtMap, hasExist);
      delete aarg2; delete app; delete p1;
      return ret;
    }

    if (strcmp(opp, "EXIST")==0)
    {
      hasExist = true;
      const ListObj* a2 = arg2(p1); // formula
      const ListObj* a1 = arg1(p1); // variables

        // without a domain, there are no constants to ground an existentially
        // quantified clause, so we must skolemize
      if (domain == NULL)
      {
        ListObj* sk = skolemize(a2, a1, vars);
        ListObj* c  = toCNF(sk, vars, domain, vtMap, hasExist);
        delete sk; delete p1;
        return c;
      }
        // get the variable names
      Array<string> varNames;
      for (int i = 0; i < a1->size(); i++)
        varNames.append((*a1)[i]->getStr());

      ArraysAccessor<int> arrAccessor;
      for (int i = 0; i < varNames.size(); i++)
      {        
        VarTypeMap::const_iterator mit = vtMap->find(varNames[i]);
        assert(mit != vtMap->end());
        int typeId = (*mit).second;
        assert(typeId > 0);
        const Array<int>* constArr = domain->getConstantsByType(typeId);
        if (constArr->size() == 0)
        {
          cout << "You must declare constants for type " 
               << domain->getTypeName(typeId) 
               << " in order to ground the variable " << varNames[i] << endl;
          exit(-1);
        }
        arrAccessor.appendArray(constArr);
      }

      if (varNames.size()==0){cout<<"ERROR: EXIST has no vars."<<endl;exit(-1);}

      ListObj* newFormula = new ListObj;
      newFormula->append("v");
      Array<int> constIds;
      while (arrAccessor.getNextCombination(constIds))
      {
        ListObj* replacedFormula = new ListObj(a2);

        assert(constIds.size() == varNames.size());
        for (int i = 0; i < varNames.size(); i++)
        {
          const char* constName = domain->getConstantName(constIds[i]);
          replace(replacedFormula, varNames[i].c_str(), constName, domain);
        }
        newFormula->append(replacedFormula);
      }
      //arrAccessor.deleteArraysAndClear();
      //newFormula->print(cout); cout << endl;

      delete p1;
      ListObj* c  = toCNF(newFormula, vars, domain, vtMap, hasExist);

      delete newFormula; 
      return c;
    }

    return p1;
  }


  void replaceAsterisk(const Array<bool>& bArr, int& idx)
  {
    if (!isStr())
    {
      const char* opp = op(this)->str_;

      if (strcmp(opp, "*")==0)
      {
        assert(size()==2);
          // if false, replace the * with !, 
          // else set (*this) to be the formula without the !
        if (!bArr[idx]) 
          (*this)[0]->setStrIfStr("!");
        else            
        { 
          ListObj* unNegatedFormula = new ListObj((*this)[1]); 
          assert(!unNegatedFormula->isStr());
          deleteAndClear();
          for (int i = 0; i < unNegatedFormula->size(); i++)
            append((*unNegatedFormula)[i]);
          unNegatedFormula->clear();
          delete unNegatedFormula;
        }
        idx++;
        if (idx == bArr.size()) return;
      }

      for (int i = 0; i < size(); i++)
      {
        (*this)[i]->replaceAsterisk(bArr, idx);
        if (idx == bArr.size()) return;
      }
    }
  }


    // removes $ and .10 from $x.10
  void cleanUpVars()
  {
    if (isStr())
    {
      if (str_[0] != '$') return;
      char* dot = strrchr(str_, '.'); 
      if (dot==NULL) return; // dealing with skolem constant of function
      int numChar = dot-str_-1;
      char* tmp = new char[numChar+1];
      strncpy(tmp, str_+1, numChar);
      tmp[numChar] = '\0';
      delete [] str_;
      str_ = tmp;
    }
    else
      for (int i = 0; i < list_.size(); i++)
        list_[i]->cleanUpVars(); 
  }


  void removeRedundantPredicates()
  {
    const char* opp = op(this)->str_;

    if (strcmp(opp,"^")==0)
    {
      for (int i = 1; i < size(); i++) // for each clause
        removeRedundantPredicates((*this)[i]);      
    }
    else
      removeRedundantPredicates(this);
  }


  void removeRedundantClauses()
  {
    //for each clause, remove redundant predicates
    //for each clause
    //  compare it to every other clause c
    //  if it contains c, remove it from *this

    const char* opp = op(this)->str_;

    if (strcmp(opp,"^")==0)
    {
      for (int i = 1; i < size(); i++) // for each clause
        removeRedundantPredicates((*this)[i]);
      
      Array<ListObj*> clauses, redundantClauses;
      for (int i = 1; i < size(); i++)
      {
        bool keep = true;
        for (int j = 1; j < size(); j++)
        {
          if (i==j) continue;
          if (subsetOf((*this)[j],(*this)[i])) 
            if ((*this)[i]->size() != (*this)[j]->size() || i < j) 
            { keep = false; break; }
        }
      if (keep) clauses.append((*this)[i]);
      else      redundantClauses.append((*this)[i]);
      }

      ListObj* andOper = (*this)[0];
      assert(strcmp(andOper->str_,"^")==0);
        //delete the operator and redundant clauses
      for (int i = 0; i < redundantClauses.size(); i++)
        delete redundantClauses[i];
      clear();
      append(andOper);
      for (int i = 0; i < clauses.size(); i++) append(clauses[i]);
    }
    else
      removeRedundantPredicates(this);
  }


  static void replace(const ListObj* const & p, const char* const & varName,
                      const char* const& constName, const Domain* const& domain)
  {
    const char* opp = op(p)->str_;
    
      //if there is variable of the same name that is existentially/universally
      //in an inner scope, then return;
    if (strcmp(opp, "EXIST")==0 || strcmp(opp, "FORALL")==0)
      if (containsVarName((*p)[1], varName)) return;

      // if we are currently looking at a predicate or a function
    const bool isPredFunc = (   domain->isPredicate(opp) 
                             || domain->isFunction(opp));

    for (int i = 1; i < p->size(); i++)
    {
        // if we are currently looking at a predicate/function 
        // AND the param is a possible var
      if (isPredFunc && (*p)[i]->isStr())
      {
          //if we found the variable
        if (strcmp((*p)[i]->getStr(), varName)==0)
          (*p)[i]->setStrIfStr(constName);
      }
      else
        replace((*p)[i], varName, constName, domain);
    }
  }


  void replace(const char* const & oldop, const char* const & newop)
  {
    for (int i = 0; i < size(); i++)
    {
      if (strcmp((*this)[i]->getStr(), oldop)==0)
        (*this)[i]->setStrIfStr(newop);
      else
        (*this)[i]->replace(oldop, newop);
    }
  }

  
  bool hasOrOp()
  {
    if (!isList()) return false;
    return (strcmp(op(this)->str_,"v")==0);
  }

  
  bool hasAndOp()
  {
    if (!isList()) return false;
    return (strcmp(op(this)->str_,"^")==0);
  }


 private:
  static bool same(const ListObj* const & lo1, const ListObj* const & lo2)
  {
    if (lo1->isStr() && lo2->isStr()) 
    {
      if (strcmp(lo1->str_, lo2->str_)==0) return true;
      return false;
    }

    if (lo1->isList() && lo2->isList())
    {
      if (lo1->size() == lo2->size())
      {
        for (int i = 0; i < lo1->size(); i++)
          if (!same((*lo1)[i],(*lo2)[i])) return false;
        return true;
      }
      return false;        
    }
    return false;
  }


  void sub(const ListObj* const & oldlo, const ListObj* const & newlo)
  {
    if (isStr()) 
    {
      if (same(this, oldlo))
      {
        if (newlo->isStr()) 
        {
          delete [] str_;
          str_ = new char[strlen(newlo->str_)+1];
          strcpy(str_, newlo->str_);
        }
        else
        {
          delete [] str_;
          str_ = new char[1];
          str_[0] = '\0';
          list_.deleteItemsAndClear();

          for (int i = 0; i < newlo->size(); i++)
            list_.append(new ListObj((*newlo)[i]));
        }
      }
      return;
    }
    
    for (int i = 0; i < list_.size(); i++)
    {
      if (same(list_[i],oldlo)) {delete list_[i]; list_[i]=new ListObj(newlo); }
      else                       list_[i]->sub(oldlo, newlo);
    }
  }


  static const ListObj* op(const ListObj* const & exp)
  {
    assert(exp->isList());
    if (exp->size() == 0) return exp;
    return (*exp)[0];
  }


  static const ListObj* arg1(const ListObj* const & exp)
  { 
    if (exp->isList() && exp->size() > 1) return (*exp)[1];
    else                                  return null_;
  }
  

  static const ListObj* arg2(const ListObj* const & exp)  
  { 
    if (exp->isList() && exp->size() > 2) return (*exp)[2];
    else                                  return null_;
  }
  


  static bool isConnective(const char* const & c)
  {
    if (strcmp(c, "^")==0 || strcmp(c, "v")==0 || strcmp(c, "!")==0 || 
        strcmp(c, "=>")==0 || strcmp(c, "<=>")==0 )
      return true;
    return false;
  }


  static bool isQuantifier(const char* const & q)
  {
    if (strcmp(q, "FORALL")==0 || strcmp(q, "EXIST")==0)
      return true;
    return false;
  }
  

    // An atomic clause has no connectives or quantifiers.
  static bool isAtomicClause(const ListObj* const & sentence)
  {
      //op(sentence) is a temporary object. Assigning opp as in 
      //"const char* opp = op(sentence).str_" cause opp to be mangled.
    const char* opp = op(sentence)->str_;
    if (!isConnective(opp) && !isQuantifier(opp)) return true;
    return false;      
  }


  static bool isNegatedClause(const ListObj* const & sentence)
  {
    if (strcmp(op(sentence)->str_, "!")==0) return true;
    return false;
  }


  static bool isLiteralClause(const ListObj* const & sentence)
  {
    if (isAtomicClause(sentence) || (isNegatedClause(sentence) 
                                     && isAtomicClause(arg1(sentence))))
      return true;
    return false;
  }


  static ListObj* eliminateImplication(const ListObj* const & p)
  {
    if (isLiteralClause(p)) return new ListObj(p);

    const char* opp = op(p)->str_;

    if (strcmp(opp, "=>")==0)
    {
      ListObj* n1 = new ListObj; 
      n1->append("!"); n1->append(new ListObj(arg1(p)));
      ListObj* n = new ListObj; 
      n->append("v"); n->append(new ListObj(arg2(p))); n->append(n1);
      return n;
    }
    else
    if (strcmp(opp, "<=>")==0)
    {        
      ListObj* n1a = new ListObj;
      n1a->append("!"); n1a->append(new ListObj(arg1(p)));
      ListObj* na = new ListObj; 
      na->append("v"); na->append(n1a); na->append(new ListObj(arg2(p)));

      ListObj* n1b = new ListObj; 
      n1b->append("!"); n1b->append(new ListObj(arg2(p)));
      ListObj* nb = new ListObj; 
      nb->append("v"); nb->append(new ListObj(arg1(p)));
      nb->append(n1b); 
      
      ListObj* n = new ListObj; n->append("^"); n->append(nb); n->append(na);
      return n;
    }
    else
    {
      ListObj* n = new ListObj; n->append(opp);
      for (int i = 1; i < p->size(); i++) // for each arg
        n->append(eliminateImplication((*p)[i]));
      return n;
    }
  }


  static ListObj* moveNotInwards(const ListObj* const & p)
  {
    const char* opp = op(p)->str_;

    if (strcmp(opp, "!")==0) return new ListObj((*p)[1]);

    if (strcmp(opp, "^")==0)
    {
      ListObj* n = NULL;
      int numDisj = 0;
      ListObj* firstDisj = NULL;

      for (int i = 1; i < p->size(); i++) //for each arg
      {
        numDisj++;
        if (numDisj == 1) { firstDisj = moveNotInwards((*p)[i]); continue; }
        else 
        if (numDisj == 2) {n=new ListObj; n->append("v"); n->append(firstDisj);}
        n->append(moveNotInwards((*p)[i]));
      }
      if (numDisj == 1) return firstDisj;
      return n;
    }

    if (strcmp(opp, "v")==0)
    {
      ListObj* n = NULL;
      int numConj = 0;
      ListObj* firstConj = NULL;
      
      for (int i = 1; i < p->size(); i++) //for each arg
      {
        numConj++;
        if (numConj == 1) { firstConj = moveNotInwards((*p)[i]); continue; }
        else 
        if (numConj == 2) { n=new ListObj; n->append("^");n->append(firstConj);}
        n->append(moveNotInwards((*p)[i]));
      }

      if (numConj == 1) return firstConj;     
      return n;
    }

    if (strcmp(opp, "FORALL")==0)
    {
      ListObj* n = new ListObj;
      n->append("EXIST"); 
      n->append(new ListObj((*p)[1]));
      n->append(moveNotInwards((*p)[2]));
      return n;
    }

    if (strcmp(opp, "EXIST")==0)
    {
      ListObj* n = new ListObj;
      n->append("FORALL");
      n->append(new ListObj((*p)[1]));
      n->append(moveNotInwards((*p)[2]));
      return n;
    }    

    ListObj* n = new ListObj; n->append("!");  n->append(new ListObj(p));
    return n;
  }


  static ListObj* disjunctionAndAppend(const ListObj* const & p1, 
                                       const ListObj* const & p2)
  {
    bool p1IsDisjunc = (strcmp(op(p1)->str_, "v")==0);
    int ibeg = (p1IsDisjunc) ? 1 : 0;
    int iend = (p1IsDisjunc) ? p1->size() : 1;

    bool p2IsDisjunc = (strcmp(op(p2)->str_, "v")==0);
    int jbeg = (p2IsDisjunc) ? 1 : 0;
    int jend = (p2IsDisjunc) ? p2->size() : 1;

    ListObj* n = new ListObj();
    if (p1->size()+p2->size() > 1) n->append("v");

    for (int i = ibeg; i < iend; i++)
    {
      if (p1IsDisjunc) n->append(new ListObj((*p1)[i]));
      else             n->append(new ListObj(p1));
    }

    for (int j = jbeg; j < jend; j++)  
    {
      if (p2IsDisjunc) n->append(new ListObj((*p2)[j]));
      else             n->append(new ListObj(p2));
    }
    return n;
  }


  static ListObj* mergeDisjuncts(const ListObj* const & ddisjuncts, 
                                 const int& numDisj)
  {
    assert(ddisjuncts->size() >= 1);    
    if (numDisj == 1) return new ListObj((*ddisjuncts)[ ddisjuncts->size()-1 ]);
  
    ListObj* y = mergeDisjuncts(ddisjuncts, numDisj-1);
    ListObj* x = new ListObj((*ddisjuncts)[ ddisjuncts->size()-numDisj ]);
    int numConj = 0;
    ListObj* result = NULL;
    ListObj* firstConj = NULL;
    ListObj* yj, * xi;

    bool yIsConj = strcmp(op(y)->str_, "^")==0;
    int jbeg = (yIsConj) ? 1 : 0;
    int jmax = (yIsConj) ? y->size() : 1;

    bool xIsConj = strcmp(op(x)->str_, "^")==0;
    int xbeg = (xIsConj) ? 1 : 0;
    int xmax = (xIsConj) ? x->size() : 1;

    for (int j = jbeg; j < jmax; j++)
    {
      if (yIsConj) yj = (*y)[j];
      else         yj = y;

      for (int i = xbeg; i < xmax; i++)
      {
        if (xIsConj) xi = (*x)[i];
        else         xi = x;

        numConj++;
        if (numConj == 1) 
        { 
          firstConj = disjunctionAndAppend(xi, yj); 
          continue; 
        }
        else
        if (numConj == 2) 
        { 
          result = new ListObj(); 
          result->append("^"); 
          result->append(firstConj);
        }

        result->append(disjunctionAndAppend(xi, yj));
      }
    }

    delete y; delete x;
    if (numConj == 1) return firstConj;
    return result;
  }


  static bool isVariable(const ListObj* const & x) 
  {
    if (strlen(x->str_) == 0) return false;
    return ((x->str_)[0] == '$'); 
  }


  static ListObj* newVariable(const ListObj* const & var)
  {
    assert(var->isStr());
    string s;
    if (!isVariable(var)) s.append("$");
    char buf[128];
    sprintf(buf, "%d", ++newVarCounter_);
    s.append(var->str_).append(".").append(buf);
    return new ListObj(s.c_str());
  }


  static ListObj* skolemConstant(const ListObj* const & name)
  {
    assert(name->isStr());
    char buf[1024];
    sprintf(buf, "$%s_%d", name->str_, ++newVarCounter_);
    return new ListObj(buf);
  }
  
    //only used in skolemize()
  static ListObj* cons(ListObj* const & p1, const ListObj* const & p2)
  {
    assert(p2->isList());
    ListObj* n = new ListObj;
    n->append(p1);
    for (int i = 0; i < p2->size(); i++) 
      n->append(new ListObj((*p2)[i]));
    return n;
  }


  static ListObj* skolemize(const ListObj* const & p, 
                            const ListObj* const & vars, 
                            const ListObj* const & outsideVars)
  {
    ListObj* pcopy = new ListObj(p);
    
    ListObj* skolemSym = new ListObj; //skolem constants or functions
    if (outsideVars->size() == 0)
    {
      for (int i = 0; i < vars->size(); i++) 
        skolemSym->append(skolemConstant((*vars)[i]));
    }
    else
    {
      for (int i = 0; i < vars->size(); i++) 
        skolemSym->append( cons(skolemConstant((*vars)[i]),outsideVars) );
    }
    assert(skolemSym->size() == vars->size());
    
    for (int i = 0; i < pcopy->size(); i++)
      for (int j = 0; j < skolemSym->size(); j++)
        (*pcopy)[i]->sub((*vars)[j], (*skolemSym)[j]);

    delete skolemSym;

    return pcopy;
  }


  static bool containsVarName(const ListObj* const & vars, 
                              const char* const & varName)
  {
    for (int i = 0; i < vars->size(); i++)
      if (strcmp((*vars)[i]->getStr(), varName)==0) return true;
    return false;
  }


    // remove redundant predicates from a clause with two or more predicates
  static void removeRedundantPredicates(ListObj* const & p)
  {
    const char* opp = op(p)->str_;

      // if not a clause with two or more predicates
    if (strcmp(opp,"v")!=0) return; 

    Array<ListObj*> preds, redundantPreds;   
    for (int i = 1; i < p->size(); i++)
    {
      bool unique = true;
      for (int j = i+1; j < p->size(); j++)
        if (same((*p)[i],(*p)[j])) { unique = false; break; }
      if (unique) preds.append((*p)[i]);
      else        redundantPreds.append((*p)[i]);
    }
    
    ListObj* orOper = (*p)[0];
    assert(strcmp(orOper->str_,"v")==0);
    for (int i = 0; i < redundantPreds.size(); i++) 
      delete redundantPreds[i];
    p->clear();
      // If clause has reduced to a unit clause, we can't just append
    if (preds.size() > 1)
    {
      p->append(orOper);
      for (int i = 0; i < preds.size(); i++) p->append(preds[i]);
    }
    else if (preds.size() == 1)
    {
      p->str_ = preds[0]->str_;
      p->list_ = preds[0]->list_;
    }
  }


  static bool subsetOf(const ListObj* const & p, const ListObj* const & q)
  {
    if (p->size() > q->size()) return false;
    for (int i = 0; i < p->size(); i++)
    {
      bool sameAsOneElt = false;
      for (int j = 0; j < q->size(); j++)
        if (same((*p)[i],(*q)[j])) { sameAsOneElt = true; break; }
      if (!sameAsOneElt) return false;
    }
    return true;
  }


 private:
    // If the formula is (v (P x) (Q y)), 
    // list_[0]: is a ListObj with str_ == "v"
    // list_[1]: is a ListObj lo1 s.t lo1.list_[0].str_ == "P", 
    //                                lo1.list_[1].str_ == "x"

  Array<ListObj*> list_;

  char* str_;

  static int newVarCounter_;
  static ListObj* null_;
};


inline
ostream& operator<<(ostream& out, const ListObj& p) { return p.print(out); }


#endif
