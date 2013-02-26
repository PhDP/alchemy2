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
#ifndef POLYNOMIAL_H_JAN_6_2007
#define POLYNOMIAL_H_JAN_6_2007

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include "array.h"

using namespace std;

typedef map<int, double> VarCont;

struct DoubleRange
{
  DoubleRange()
  {
    low = 0;
    high = 0;
    iType = 0;
  }
  double low;
  double high;
    // Type of the double range:
    // 0 - inward range, 1 - outward range, 2 - 1 point range, 3 - empty range.
  int iType; 
};

template <typename Type>
bool ArrayCompare(Array<Type> const & a1, Array<Type> const & a2);

class PolyNomial 
{
public:	
  PolyNomial();
  PolyNomial(const PolyNomial& pl);
  PolyNomial(PolyNomial& pl);
  ~PolyNomial();

  void Clear();

    // Explicit copy function.
  void Copy(PolyNomial& pl);
  void Copy(const PolyNomial& pl);

    // operator=
  void operator=(PolyNomial& pl);
  void operator=(const PolyNomial& pl);

    // Compute the value of the polynomial, given values for the variables.
  double ComputePlValue(const Array<double>& vars);
  double ComputePlValue();

    // Compute the polynomial's first-order derivative (a polynomial) at the given variable.
  PolyNomial GetGradient(int varInIdx);
    // Compute the polynomial's first-order derivative values for each variable.
  Array<double> GetGradient();

    // Normalize the polynomial to canonical base form, merge items with similar bases, etc.
  void Normalize();
    // Normalize the polynomial to log Gaussian format (x-a)^2 and return true, if not single-var quadratic, return false.
  bool NormalizeGaussian();

    // i.o functions.
  void ReadFrom(istream& is);  // Read polynomial from input stream.
  void PrintVars(ostream& os) const;  // Print the polynomial variables to stream.
  void PrintTo(ostream& os) const;  // Print the polynomial to stream.
  void PrintCoef(ostream& os, double coef) const;  // Print the coefficient value, dealing with the two number sign cases.

    // Polynomial computation functions.
  void MultiplyConst(double coef);  // Multiply the polynomial with a constant value.
    // Assume the polynomial is representing a single variable Gaussian
    // this function returns the paramters to the Gaussian distribution represented by the polynomial.
  void GetGaussianPara(double* miu, double* stdev);
  void AddPl(const PolyNomial& pl);

    // Return the highest order of the given variable.
  double GetHighestOrder(int varInIdx);

    // Add variable to polynomial and keep the polynomial unchanged.
  bool AddVar(const string& var)
  {
    if (varsearch_.find(var) != varsearch_.end())
    {
      return false;
    }
    varsearch_[var] = numVar_;
    numVar_++;		
    var_.append(var);
    varValue_.append(0.0);  // Default value for the new variable is set to 0.
    return true;
  }

  bool AddItem(const VarCont& item, double coef)
  {
    VarCont::const_iterator cit = item.begin();
    for (; cit != item.end(); ++cit)
    {
      int var_id = cit->first;
      if (var_id >= numVar_) return false;
    }
    numItems_ ++;
    coef_.append(coef);
    items_.append(item);
    strItems_.append(GenerateItemString(item));
    Normalize();
    return true;
  }

  void SetVarName(const Array<string>& varName)
  {
    assert(varName.size()== numVar_);
    var_.clear();
    var_.append(varName);

    varsearch_.clear();
    for (int i = 0; i < var_.size(); ++i)
    {
      varsearch_[var_[i]] = i;
    }
  }

    // Set values of variables.
  void SetVarValue(const Array<double>& varValue)
  {
    assert(varValue.size()== numVar_);
    varValue_.clear();
    varValue_.append(varValue);
  }

  void SetConstantValue(double val)
  {
    constValue_ = val;
  }

  double GetConstantValue()
  {
    return constValue_;
  }

    // Reduce the polynomial to the polynomial over one specific variable,
    // given the values of the rest.
  void ReduceToOneVar(const Array<double>& vars, int varIdx);
  
  void ReduceToOneVar(int varIdx);

    // Generate the string of specific item in the polynomial. 
  string GenerateItemString(const VarCont& vc);

    // Solve quadratic polynomial constraint.
  DoubleRange SolvePoly2(double d);  //assume it's a GREATER THAN (">=") polynomial constraints
  DoubleRange SolvePoly2(double a, double b, double c, double d);  //c != 0, constraint type a + b*x +c*x^2 >= d
	
    // Do quadratic optimization of the quadratic polynomials, and return the value of optimized polynomial.
  double QuadraticOptimization();

  const Array<double>& GetVarValue() const {return varValue_;}
  double GetVarValue(int varIdx) {return varValue_[varIdx];}
  void ClearVarValue() {varValue_.clear();}
  void AppendVarValue(double value) {varValue_.append(value);}
  const int GetVarNum() const {return numVar_;}
  
  const string& GetVarAt(int varIdx) const {return var_[varIdx];}

  int GetVarIdx(string& var) const
  {
    map<string, int>::const_iterator citer = varsearch_.find(var);
    if (citer != varsearch_.end()) 
      return citer->second;
    else return -1;
  }

  const double& GetCoef(int varIdx) const { return coef_[varIdx];}
  const Array<string> & GetVar() const {return var_;}
  const map<string, int>& GetVarSearch() const{ return varsearch_;}

    // Compute the value of the given item, given the values of variables.
  double ComputeItem(int itemIdx);
  bool IsQuadratic();  // Check if the polynomial is a log Gaussian term, i.e., if it's a quadratic polynomial.
  int Optimize2(Array<double>* vars, double* optValue);  // Optimizing the polynomial when it's a quadratic one. Values returned: assignment to variables and the optimum value of the polynomial.
    // General optimization interface for arbitrary polynomial is implemented in function minimize in lbfgsp.h

    // Parameters for quadratic polynomials.
  struct QuadraticPolyPara
  {
    double a; //constant
    double b; // coefficient of 1 order term
    double c; // coefficient of 2 order term
  };
  
  Array<double> varValue_;  // Holder for value assignment to variables.

public:	

  int numVar_;
  int numItems_;
  double constValue_;	
  QuadraticPolyPara gpara_;  // Parameters for quadratic polynomial representations.
  Array<string> strItems_;  // Holder for each item's string representation.	
  Array<double> coef_;  // From constant to higher orders, coef_[0] is always 1. ??
  Array<VarCont> items_;  // Coefficient descriptors of different terms in a polynomial, each element in this vector is a map from variable index to its order value in the term.
  Array<string> var_;  // We assume within a polynomial different variables have different names.
  map<string, int> varsearch_;  // Map from variable string to its index.
};

#endif
