/*
 * LogDouble.cpp
 *
 *  Created on: Nov 5, 2011
 *      Author: Vibhav Gogate
 *				The University of Texas at Dallas
 *				All rights reserved
 */

#include "logdouble.h"

// Default Constructor
LogDouble::LogDouble() :
		value(0.0), is_zero(true) {
}
// Initialize the logdouble using a long double
LogDouble::LogDouble(const long double d, bool logspace) {
	if (logspace) {
		value = d;
		is_zero = false;
	} else {
		if (d <= 0) {
			value = 0.0;
			is_zero = true;
		} else {
			value = log(d);
			is_zero = false;
		}
	}
}
//Copy constructor
LogDouble::LogDouble(const LogDouble& dub) :
		value(dub.value), is_zero(dub.is_zero) {
}

LogDouble& LogDouble::operator=(const LogDouble& other) {
	value = other.value;
	is_zero = other.is_zero;
	return *this;
}

void LogDouble::setValue(long double val) {
	value = val;
	if(val!=0)
		is_zero = false;
}

LogDouble LogDouble::operator+(const LogDouble& other) const {
	LogDouble out(other);
	if (is_zero){
		return out;
	}
	if (out.is_zero){
		out=*this;
		return out;
	}
	if (value > out.value) {
		out.value = log(1 + exp(out.value - value)) + value;
	} else {
		out.value += log(1 + exp(value - out.value));
	}
	return out;
}
LogDouble& LogDouble::operator+=(const LogDouble& other) {
	if (is_zero) {
		*this = other;
		return *this;
	}
	if (other.is_zero)
		return *this;
	if (value > other.value) {
		value += log(1 + exp(other.value - value));
	} else {
		value = log(1 + exp(value - other.value)) + other.value;
	}
	return *this;
}
LogDouble LogDouble::operator*(const LogDouble& other) const {
	if (is_zero || other.is_zero)
	{
		return LogDouble();
	}
	LogDouble out(other);
	out.value += value;
	return out;
}
LogDouble& LogDouble::operator*=(const LogDouble& other) {
	if (is_zero || other.is_zero){
		value=0.0;
		is_zero=true;
		return *this;
	}
	value += other.value;
	return *this;
}
bool LogDouble::operator <(const LogDouble& other) const {
	if (other.is_zero)
		return false;
	if (is_zero)
		return true;
	return value < other.value;
}
bool LogDouble::operator >(const LogDouble& other) const {
	if (other.is_zero)
		return true;
	if (is_zero)
		return false;
	return value > other.value;
}
LogDouble LogDouble::operator/(const LogDouble& other) const {
	if (is_zero || other.is_zero)
	{
		return LogDouble();
	}
	LogDouble out(value);
	out.value -= other.value;
	return out;
}

void LogDouble::LDPower(LogDouble in,int n,LogDouble& out)  {
	if (n==0)
	{
		return;
	}
	if(in.is_zero)
	{
		out.is_zero = true;
		return;
	}
	out=LogDouble(in.value*n);
	
	//for(unsigned int i=0;i<n;i++)
		//out = out*in;
	
}


LogDouble LogDouble::Factorial(int n)
{
	if(n==0)
		return LogDouble(1.0,false);
	LogDouble fact(n,false);
	for(unsigned int i=n-1;i>1;i--)
	{
		LogDouble t (i,false);
		fact*=t;
	}
	return fact;
}

LogDouble LogDouble::Binomial(int n,int r)
{
	
	LogDouble factorial = Factorial(n)/(Factorial(n-r)*Factorial(r));	
	return factorial;
}
