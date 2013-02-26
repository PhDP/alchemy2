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
#include "Polynomial.h"

#include <cstdio>
#include <iomanip>
#include <math.h>

#include "complex.h"

#define DOUBLE_ZERO_THRESHOLD 0.000000001
#define NOOPT -1
double PBIGVALUE = 1010101;

bool pldbg = false;

bool DoubleEqual(const double& d1, const double& d2, const double& threshold) {
	if (fabs(d1 - d2) < threshold) {
		return true;
	}
	return false;
}

// Analytical solution to quadratic polynomials exists. This function could return the optimal solution to quadratic polynomials with one or two variables.
// Return value: indicator variable for solution status, NOOPT
int PolyNomial::Optimize2(Array<double>* vars, double* optValue)
{
	vars->clear();
	if (1 == numVar_)  // There is only 1 variable in the polynomial.
	{
		double a = 0, b = 0;  // Get parameters for the Gaussian representation.
		for(int i = 0; i < items_.size(); ++i) {
			map<int, double>::const_iterator citer = items_[i].begin();
			if (DoubleEqual(citer->second, 2.0, DOUBLE_ZERO_THRESHOLD)) {
				a = coef_[i];
			}

			if (DoubleEqual(citer->second, 1.0, DOUBLE_ZERO_THRESHOLD)) {
				b = coef_[i];
			}
		}

		// ERROR: The quadratic term is empty.
		if (fabs(a) < DOUBLE_ZERO_THRESHOLD) {
			cout << "The quadratic term is empty" << endl;
			exit(0);
		}

		vars->append(-b/(2*a));
		*optValue = ComputePlValue(*vars);
		return 0;
	}

	// Two-var quadratic polynomial representation: a*x1^2 + b*x2^2 + c*x1*x2 + d*x1 + e*x2.
	double a = 0, b = 0, c = 0, d = 0, e = 0;
	//normalize and assume each polynomial are 2-var poly.
	for(int i = 0; i < items_.size(); ++i) {
		map<int, double>::const_iterator citer = items_[i].begin();
		if (citer->first == 0 && DoubleEqual(citer->second, 2.0, DOUBLE_ZERO_THRESHOLD)) //x1^2
		{
			a = coef_[i];
		}
		if (citer->first == 1 && DoubleEqual(citer->second, 2.0, DOUBLE_ZERO_THRESHOLD)) //x2^2
		{
			b = coef_[i];
		}
		if (citer->first == 0 && DoubleEqual(citer->second, 1.0, DOUBLE_ZERO_THRESHOLD) && items_[i].size() == 2) //x1*x2
		{
			c = coef_[i];
		}
		if (citer->first == 0 && DoubleEqual(citer->second, 1.0, DOUBLE_ZERO_THRESHOLD) && items_[i].size() == 1) //x1
		{
			d = coef_[i];
		}
		if (citer->first == 1 && DoubleEqual(citer->second, 1.0, DOUBLE_ZERO_THRESHOLD) && items_[i].size() == 1)//x2
		{
			e = coef_[i];
		}
		if (citer->first == 1 && DoubleEqual(citer->second, 1.0, DOUBLE_ZERO_THRESHOLD) && items_[i].size() > 1) //wrong
		{
			cout << "shouldn't be here," << endl;
		}
	}

	vars->clear();
	double delta = 4 * a * b - c * c;
	//cout << "delta: " << delta << endl;
	// delta == 0, |a| > 0, |b| > 0
	if (fabs(delta) < DOUBLE_ZERO_THRESHOLD) {
		if (fabs(a)> DOUBLE_ZERO_THRESHOLD) {
			vars->append(-d / 2 * a);
			vars->append(0.0);
		} else if (fabs(a)> DOUBLE_ZERO_THRESHOLD) {
			vars->append(0.0);
			vars->append(-e / 2 * b);
		} else {
			// In this case, there is no appropriate optimum
			PrintTo(cout); cout << endl;
			cout << "delta = 0, no optimums2" << endl;
			exit(0);
		}
		*optValue = ComputePlValue(*vars);
		return NOOPT;
	}

	vars->append((e * c-2 * b * d)/delta);
	vars->append((d * c-2 * a * e)/delta);
	*optValue = ComputePlValue(*vars);
	return 0;
}

// Solve quadratic constraint which is of the form: f(x) > 0
// c != 0, constraint type a + b*x +c*x^2 >= d
DoubleRange PolyNomial::SolvePoly2(double a, double b, double c, double d)
{
	DoubleRange r;
	if (c == 0)
	{
		if (b == 0)
		{
			if(a - d >= 0) //any value will do
			{
				r.low = - PBIGVALUE;
				r.high = PBIGVALUE;
				r.iType = 0;
			}
			else
			{
				r.iType = 3;
			}
		}
		else
		{
			if (b > 0)
			{
				r.low = d/b;
				r.high = PBIGVALUE;
				r.iType = 0;
			}
			else
			{
				r.high = d/b;
				r.low = -PBIGVALUE;
				r.iType = 0;
			}
  		}
		return r;
	}
	
	double za = a/c;
	double zb = b/c;
	double zd = d/c;
	double delta = zd - za + zb * zb/4;

	if (delta < 0) {   
		r.iType = 3;
	} else if (delta == 0) {
		r.low = r.high = -zb/2;
		r.iType = 2;
	} else {
		r.low = -zb/2 - sqrt(delta);
		r.high = -zb/2 + sqrt(delta);
		if (c > 0) {
			r.iType = 1;
		} else {
			r.iType = 0;
		}
	}
	return r;
}

DoubleRange PolyNomial::SolvePoly2(double d)
{
	DoubleRange r;
	if (!NormalizeGaussian())
	{
		cout << "unsolvable" << endl;
		r.low = 1;
		r.high = 0;
		r.iType = 3;
		return r;
	}

	return SolvePoly2(gpara_.a, gpara_.b, gpara_.c, d);
}

PolyNomial::PolyNomial()
{
	numVar_ = 0;
	constValue_ = 0;
	numItems_ = 0;
}

PolyNomial::PolyNomial(const PolyNomial& pl)
{
	numVar_ = pl.numVar_;
	var_ = pl.var_;
	coef_ =  pl.coef_;
	numItems_ = pl.numItems_;
	strItems_ = pl.strItems_;
	items_ = pl.items_;
	constValue_ = pl.constValue_;
	varValue_ = pl.varValue_;

	varsearch_.clear();

	const map<string, int>& varsearch = pl.varsearch_;
	map<string, int>::const_iterator iter = varsearch.begin(); 
	for(;iter != varsearch.end(); ++iter) {
		varsearch_.insert(map<string, int>::value_type(iter->first, iter->second));
	}
}

PolyNomial::PolyNomial(PolyNomial& pl)
{
	numVar_ = pl.numVar_;
	var_ = pl.var_;
	coef_ =  pl.coef_;
	numItems_ = pl.numItems_;
	strItems_ = pl.strItems_;
	items_ = pl.items_;
	constValue_ = pl.constValue_;
	varValue_ = pl.varValue_;

	varsearch_.clear();

	const map<string, int>& varsearch = pl.varsearch_;
	map<string, int>::const_iterator iter = varsearch.begin(); 
	for(;iter != varsearch.end(); ++iter)
	{
		varsearch_.insert(map<string, int>::value_type(iter->first, iter->second));
	}
}

double PolyNomial::QuadraticOptimization()
{
	if (!NormalizeGaussian())
	{
		cout << "unsolvable" << endl;
		return PBIGVALUE;
	}

	return -gpara_.b / (2 * gpara_.c);
}

void PolyNomial::Copy(PolyNomial& pl)
{
	numVar_ = pl.numVar_;
	var_ = pl.var_;
	coef_ =  pl.coef_;
	numItems_ = pl.numItems_;
	strItems_ = pl.strItems_;
	items_ = pl.items_;
	constValue_ = pl.constValue_;
	varValue_ = pl.varValue_;

	varsearch_.clear();

	const map<string, int>& varsearch = pl.varsearch_;
	map<string, int>::const_iterator iter = varsearch.begin(); 
	for(;iter != varsearch.end(); ++iter)
	{
		varsearch_.insert(map<string, int>::value_type(iter->first, iter->second));
	}
}


void PolyNomial::Copy(const PolyNomial& pl)
{
	numVar_ = pl.numVar_;
	var_ = pl.var_;
	coef_ =  pl.coef_;
	numItems_ = pl.numItems_;
	strItems_ = pl.strItems_;
	items_ = pl.items_;
	constValue_ = pl.constValue_;
	varValue_ = pl.varValue_;

	varsearch_.clear();

	const map<string, int>& varsearch = pl.varsearch_;
	map<string, int>::const_iterator iter = varsearch.begin(); 
	for(;iter != varsearch.end(); ++iter)
	{
		varsearch_.insert(map<string, int>::value_type(iter->first, iter->second));
	}
}

void PolyNomial::operator=(PolyNomial& pl)
{
	Copy(pl);
}

void PolyNomial::operator=(const PolyNomial& pl)
{
	Copy(pl);
}

PolyNomial::~PolyNomial()
{

}

void PolyNomial::Clear()
{
	numVar_ = 0;
	constValue_ = 0;
	var_.clear();
	varsearch_.clear();
	coef_.clear();
	varValue_.clear();
	items_.clear(); 
	strItems_.clear();
	numItems_ = 0;
}

double PolyNomial::ComputeItem(int itemIdx)
{
	double v = coef_[itemIdx];
	map<int, double>::const_iterator citer;
	for(citer = items_[itemIdx].begin(); citer !=items_[itemIdx].end(); ++citer)
	{
		v = v * pow(varValue_[citer->first], double(citer->second));
	}
	return v;
}

void PolyNomial::ReduceToOneVar(int varIdx)
{
	if (pldbg)
	{
		cout << "Entering ReduceToOneVar" << endl;
	}
	
	assert(0<= varIdx && varIdx < var_.size());

	map<double, double> new_order_coef;
	double new_const_value = constValue_;	

	// Record the order and new coefficients for the given variable. Leave the old stuff unchanged.
	for(int i = 0; i < items_.size(); ++i)
	{
		VarCont& vc = items_[i];
		// This polynomial item does not contain the variable.
		VarCont::const_iterator cit = vc.find(varIdx);
		if ( cit == vc.end())
		{
			new_const_value += ComputeItem(i);
			//remove corresponding item
		} else {
			double coef = coef_[i];
			double order = cit->second;
			map<int, double>::const_iterator citer;
			for(citer = vc.begin(); citer != vc.end(); ++citer)
			{
				if (citer->first != varIdx)
				{
					coef *= pow(varValue_[citer->first], citer->second);
				}
			}
			new_order_coef[order] += coef;
		}
	}

	// Clear and reconstruct the old stuff.
	constValue_ = new_const_value;
	numVar_ = 1;
	double v = varValue_[varIdx];
	varValue_.clear();
	varValue_.append(v);
	string str = var_[varIdx];
	var_.clear();
	var_.append(str);
	varsearch_.clear();
	varsearch_[str] = 0;

	numItems_ = new_order_coef.size();
	items_.clear();
	strItems_.clear();
	coef_.clear();
	map<double, double>::const_iterator cit = new_order_coef.begin();
	for (; cit != new_order_coef.end(); ++cit) {
		coef_.append(cit->second);
		VarCont item;
		item[0] = cit->first;
		items_.append(item);
		strItems_.append(GenerateItemString(item));
	}

	if (pldbg)
	{
		cout << "Leaving ReduceToOneVar" << endl;
	}	
}


//void PolyNomial::ReduceToOneVar(int varIdx)
//{
//	cout << "Entering ReduceToOneVar" << endl;
//	assert(0<= varIdx && varIdx < var_.size());
//
//	for(int i = 0; i < items_.size(); i++)
//	{
//		cout << "Checking item " << i << ", " << strItems_[i] << endl;
//		VarCont& vc = items_[i];
//		// This polynomial item does not contain the variable.
//		if (vc.find(varIdx) == vc.end())
//		{
//			cout << "Not contain var " << var_[varIdx] << ", removing item .." << endl;
//			constValue_ += ComputeItem(i);
//			//remove corresponding item
//			items_.removeItem(i);
//			coef_.removeItem(i);
//			strItems_.removeItem(i);
//			numItems_ --;
//			i--;
//		} else {
//			cout << "Contain var " << var_[varIdx] << endl;
//			map<int, double>::iterator citer;
//			double coef = coef_[i];
//			for(citer = vc.begin(); citer != vc.end(); )
//			{
//				if (citer->first == varIdx)
//				{
//					++citer;
//				}
//				else
//				{
//					coef *= pow(varValue_[citer->first], double(citer->second));
//					vc.erase(citer++);					
//				}
//			}
//
//			double pow = vc.begin()->second;
//			vc.clear();
//			vc.insert(map<int,double>::value_type(0, pow));
//
//			coef_[i] = coef;
//			strItems_[i] = GenerateItemString(vc);
//		}
//	}
//
//	numVar_ = 1;
//	double v = varValue_[varIdx];
//	varValue_.shrinkToSize(1);
//	varValue_[0] = v;
//	string str = var_[varIdx];
//	var_.shrinkToSize(1);
//	var_[0] = str;
//	varsearch_.clear();
//	varsearch_.insert(map<string,int>::value_type(str, 0));
//	Normalize();
//	numItems_ = items_.size();	
//	cout << "Leaving ReduceToOneVar" << endl;
//}



void PolyNomial::ReduceToOneVar(const Array<double>& vars, int varIdx)
{
	SetVarValue(vars);
	ReduceToOneVar(varIdx);
}

double PolyNomial::ComputePlValue(const Array<double>& vars)
{
	assert(vars.size() == numVar_);
	
	Array<double> varValueBak_;

	varValueBak_.copyFrom(varValue_);

	varValue_.copyFrom(vars);

	double plValue = ComputePlValue();

	varValue_.copyFrom(varValueBak_);

	return plValue;
}

double PolyNomial::ComputePlValue()
{
	double plValue = 0;
	plValue = constValue_;
	for(int i = 0; i < items_.size(); i++)
	{
		plValue += ComputeItem(i);
	}
	return plValue;
}

// Get the highest order number for the variable at varInIdx.
double PolyNomial::GetHighestOrder(int varInIdx) {
	// find the highest order item
	double highest = 0;
	int size = items_.size();
	for(int i = 0; i < size; i++)
	{
		VarCont::const_iterator citer;
		for(citer = items_[i].begin(); citer != items_[i].end(); ++citer)
		{
			if(citer->first == varInIdx)
			{
				if (citer->second > highest)
				{
					highest = citer->second;
				}
			}
		}
	}
	return highest;
}

void PolyNomial::Normalize()
{
	if (pldbg)
	{
		cout << "Entering Normalize" << endl;
	}
	
	map<string, int> tmpCont;
	map<string, int>::const_iterator citer;
	for(int i = 0; i < numItems_; i++)
	{
		if ((citer = tmpCont.find(strItems_[i])) != tmpCont.end())
		{
			coef_[citer->second] += coef_[i];
			strItems_.removeItem(i);
			items_.removeItem(i);
			coef_.removeItem(i);
			numItems_ --;
			i--;
		} else {
			tmpCont.insert(map<string, int>::value_type(strItems_[i], i));
		}
	}
	if (pldbg)
	{
		cout << "Leaving Normalize" << endl;
	}	
}

bool PolyNomial::NormalizeGaussian()
{
	//do nth here
	if (!IsQuadratic())
	{
		return false;
	}
	
	Normalize();

	if (numItems_ > 2)
	{
		return false;
	}

	if (numItems_ == 1)
	{
		VarCont& vc = items_[0];
		VarCont::const_iterator citer = vc.begin();
		if (citer->second != 2)
		{
			return false;
		}

		gpara_.a = constValue_; 
		gpara_.b = 0;
		gpara_.c = coef_[0];
		assert(gpara_.c != 0);
	}
	else //numitems == 2
	{
		gpara_.a = constValue_;
		for(int i = 0; i < numItems_; i++)
		{
			VarCont& vc = items_[i];
			assert(vc.size() == 1);

			VarCont::const_iterator citer = vc.begin();

			if (citer->second == 1)
			{
				gpara_.b = coef_[i];
			}
			else
			{
				gpara_.c = coef_[i];
			}

		}
	}
	return true;
}


string PolyNomial::GenerateItemString(const VarCont& vc)
{
	string str;
	VarCont::const_iterator citer;
	for(citer = vc.begin(); citer != vc.end(); ++citer) {
	   char a[100];	   
	   sprintf(a, "(%d,%f)", citer->first, citer->second);
	   str += string(a);
	}
	return str;
}

// Read polynomial from an input stream of the form:
// CONST_VALUE/var_1;var_2;...;var_n/coef_1,(var_idx_1,var_order_1)...()/coef_2,(var_idx_2,var_order_2)...()/
void PolyNomial::ReadFrom(istream& is)
{
	Clear();
	string str;
	getline(is, str, '/');
	constValue_ = atof(str.c_str());
	getline(is, str, '/');
	stringstream ss(str);
	string strVar;
	int i = 0;
	while(getline(ss, strVar, ';'))
	{
		assert(varsearch_.find(strVar) == varsearch_.end());
		var_.append(strVar);
		varsearch_.insert(map<string, int>::value_type(strVar, i));
		i++;
	}

	i = 0;
	map<string, int> strngItems;
	while (getline(is, str, '/'))
	{
		VarCont vc;
		int varIdx;
		double varOrder;
		stringstream as(str);
		string strTmp;
		getline(as, strTmp, ',');

		double coef = atof(strTmp.c_str());

		while (getline(as, strTmp, '('))
		{
			getline(as, strTmp, ',');
			varIdx = atoi(strTmp.c_str());
			getline(as, strTmp, ')');
			varOrder  = atof(strTmp.c_str());
			vc.insert(map<int,double>::value_type(varIdx, varOrder));
		}

		string strItem = GenerateItemString(vc);
		map<string, int>::const_iterator citer = strngItems.find(strItem);
		if (citer != strngItems.end())
		{
			int itemIdx = citer->second;
			coef_[itemIdx] += coef;	//dedup & merge
			continue;
		}

		strItems_.append(strItem);
		strngItems.insert(map<string, int>::value_type(strItem, i++));

		coef_.append(coef);
		items_.append(vc);  
	}

	numVar_ = var_.size();
	varValue_.growToSize(numVar_, 0);
	numItems_ = strItems_.size();
}

void PolyNomial::PrintVars(ostream& os) const
{
	for(int i = 0; i < numVar_; i++) {
		os << var_[i] << ":" << varValue_[i] << "\t";
	}
	os << endl;
}

void PolyNomial::PrintTo(ostream& os) const
{
	os.setf(ios::fixed, ios::floatfield);
	os.setf(ios::showpoint);

	os << setprecision(20) << constValue_  << '/';
	for(int i = 0; i < var_.size(); i++) {
		os << var_[i] << ';';
	}
	os << '/';
	for(int i = 0; i < items_.size(); i++) {
		os << setprecision(20) <<coef_[i] << ',';
		map<int,double>::const_iterator citer;
		for(citer = items_[i].begin(); citer != items_[i].end(); ++citer) {
			os << '(' << citer->first << ',' << setprecision(20) <<citer->second << ')' ;
		}
		os << '/';
	}
}

// Check if it's quadratic.
bool PolyNomial::IsQuadratic()
{
	if (numVar_ != 1)
	{
		return false;
	}

	double highorder = -1;
	for(int i = 0; i < items_.size(); i++)
	{
		VarCont& vc = items_[i];
		map<int, double>::const_iterator citer = vc.begin();
		if (citer->second > highorder)
		{
			highorder = citer->second;
		}
	}

	if (highorder != 2)
	{
		return false;
	}

	return true;
}

void PolyNomial::PrintCoef(ostream& os, double coef) const
{
	if (coef < 0)
	{
		os << setprecision(20) <<coef;
	} else {
		os << "+"<< setprecision(20)<< coef;
	}
}


void PolyNomial::MultiplyConst(double coef)
{
	constValue_ *= coef;
	for(int i = 0; i < numItems_; i++)
	{
	   coef_[i] *= coef;
	}
}

//Assume two pls have the same set of variables.
void PolyNomial::AddPl(const PolyNomial& pl)
{				
	if (pldbg)
	{
		cout << "Entering AddPl" << endl;
	}
	if (numVar_ == 0 && pl.numVar_ > 0)
	{
		numVar_ = pl.numVar_;
		constValue_ += pl.constValue_;
		var_.copyFrom(pl.var_);
		items_.copyFrom(pl.items_);
		coef_.copyFrom(pl.coef_);
		varValue_.copyFrom(pl.varValue_);
		strItems_.copyFrom(pl.strItems_);
		map<string,int>::const_iterator citer = pl.varsearch_.begin();
		for(;citer!= pl.varsearch_.end(); ++citer)
		{
			varsearch_.insert(map<string,int>::value_type(citer->first, citer->second));
		}
	} else if (pl.numVar_ == 0) {
		constValue_ += pl.constValue_;
	} else if (numVar_ > 0 && pl.numVar_> 0) {
		constValue_ += pl.constValue_;
		map<string, int> itemCont;
		for(int i = 0; i < strItems_.size(); i++) {
			itemCont.insert(map<string, int>::value_type(strItems_[i],i));
		}
		map<string,int>::const_iterator citer;
		for(int i = 0; i < pl.strItems_.size(); i++)
		{
			if ((citer = itemCont.find(pl.strItems_[i])) != itemCont.end())
			{
				coef_[citer->second] += pl.coef_[i];
			} else {
				items_.append(pl.items_[i]);
				strItems_.append(pl.strItems_[i]);
				coef_.append(pl.coef_[i]);
			}
		}
		numItems_ = items_.size();
	}
	if (pldbg)
	{
		cout << "Leaving AddPl" << endl;
	}
}
void PolyNomial::GetGaussianPara(double* miu, double* stdev)
{	
	if (numVar_ == 1)  // Only 1 variable
	{
		double a = 0, b = 0;
		for(int i = 0; i < items_.size(); i++)
		{
			map<int, double>::const_iterator citer = items_[i].begin();
			if (citer->second == 2) {
				a = coef_[i];
			}
			if (citer->second == 1)	{
				b = coef_[i];
			}
		}
		assert(a < 0);
		*miu = -b/(2.0 * a);
		*stdev = sqrt(-1.0/(2.0 * a));		
	} else {
		cout << "can not handle multivariate Gaussian now" << endl;
		exit(0);
	}
}

PolyNomial PolyNomial::GetGradient(int varInIdx)
{
	PolyNomial gradient;
	gradient.Copy(*this);
	//gradient.PrintTo(cout);	cout << endl;
	gradient.ReduceToOneVar(varValue_, varInIdx);
	//gradient.PrintTo(cout); cout << endl;
	gradient.constValue_ = 0;
	
	if(0 == numItems_)
	{
		gradient.Clear();
		return gradient;		
	}

	for (int i = 0; i < gradient.numItems_; i++)
	{
		//double coef = gradient.coef_[i];
		
		map<int,double>::iterator it = gradient.items_[i].begin();
		if(it->second == 1)
		{
			double coef = gradient.coef_[i];
			gradient.items_.removeItem(i);
			gradient.coef_.removeItem(i);
			gradient.strItems_.removeItem(i);
			gradient.numItems_ --;
			i--;
			gradient.constValue_ += coef;
		} else {
			gradient.coef_[i] *= it->second;
			it->second = it->second - 1;
			gradient.strItems_[i] = gradient.GenerateItemString(gradient.items_[i]);
		}		
	}
	return gradient;
}
Array<double> PolyNomial::GetGradient()
{
	Array<double> var;

	for(int i = 0; i < numVar_; i++) {
		var.append(GetGradient(i).ComputePlValue());
	}
	return var;
}

// This routine solves order 3 polynomial equation, assumed to be of the following canonical format:
//	x^3   +   ax^2   +   bx   +   c   =   0;  
//	Let x=y-a/3, we have:
//	y^3   +   py   +   q   =   0;  
//Then we introduce the following auxilliary variables:
//	w1   =   (-1+sqrt(-3))/2;
//	w2   =   (-1-sqrt(-3))/2;  
//	z1   ��   pow((-q/2 + sqrt(q*q/4+p*p*p/27)), (1/3))
//	z2   ��   pow((-q/2-sqrt(q*q/4+p*p*p/27)), (1/3))
//	y1   =   z1 + z2;  
//	y2   =   w1*z1 + w2*z2;  
//	y3   =   w2*z1 + s2*z2;  
//	Then compute x from y1 y2 and y3.
//  Solving function for order 3 polynomials.
Array<Complex> solvePolyEq3(double a, double b, double c)
{
	Array<Complex> root;
	Complex w1,w2;

	w1.Set(-0.5,sqrt(3.0) / 2);
	w2.Set(-0.5,-sqrt(3.0) / 2);

	double p = b-(1.0/3) * a * a;
	double q = 2.0 / 27.0 * a * a * a - a * b / 3.0 + c;

	double delta = q * q / 4.0 + p * p * p / 27.0;  

	// if delta < 0, too complex, we don't handle it now
	if (delta < 0)
	{
		// if delta < 0, too complex, we don't handle it now
		Complex cDelta(delta, 0.0);
		Array<Complex> ar = cDelta.Root(2);

		// pick the same square root here
		Complex zz1 = ar[0] - q / 2;
		Complex zz2 = ar[0] * (-1.0) - q/2;

		// pick the cubic root with 120 degree phase discrepancy
		ar.clear();
		ar = zz1.Root(3);
		Complex z1 = ar[0];

		ar.clear();
		ar = zz2.Root(3);
		Complex z2 = ar[2];

		Complex y1 = z1 + z2;
		Complex y2 = w1 * z1 + w2 * z2;
		Complex y3 = w2 * z1 + w1 * z2;

		Complex x1 = y1 - a / 3;
		Complex x2 = y2 - a / 3;
		Complex x3 = y3 - a / 3;

		root.append(x1);
		root.append(x2);
		root.append(x3);
	} else {
		double sqrtDelta = sqrt(delta);
		double v1 = q * (-0.5) + sqrtDelta;
		double v2 = q * (-0.5) - sqrtDelta;

		double z1;
		if (v1 >= 0)
		{
			z1 = pow(v1, 1.0 / 3.0);
		}
		else 
		{
			z1 = - pow(-v1, 1.0 / 3.0);
		}

		double z2 = pow(q * (-0.5) - sqrtDelta, 1.0 / 3.0);
		if (v2 >= 0)
		{
			z2 = pow(v2, 1.0 / 3.0);
		}
		else 
		{
			z2 = - pow(-v2, 1.0 / 3.0);
		}

		Complex y1; y1.Set(z1 + z2, 0);
		Complex y2 = w1 * z1 + w2 * z2;
		Complex y3 = w2 * z1 + w1 * z2;			

		root.append(y1 - a / 3);
		root.append(y2 - a / 3);
		root.append(y3 - a / 3);
	}

	return root;
}
