/*
 * LogDouble.h
 *
 *  Created on: Nov 3, 2011
 *      Author: Vibhav Gogate
 *				The University of Texas at Dallas
 *				All rights reserved
 */

#ifndef LOGDOUBLE_H_
#define LOGDOUBLE_H_

#include <cmath>
#include <iostream>
#ifndef _MSC_VER
   #include <limits>
#endif
using namespace std;

// To avoid underflows and overflows, we maintain real numbers in the log space
// The following class provides the necessary functionality and operators for real numbers
// in the log-space

struct LogDouble {
	// public variables
	long double value;
	bool is_zero;
	bool is_inf;

	// Default constructor, value=0.0 and is_zero=false
	LogDouble();
	// Initialize using a long double value d
	// Set the variable in_logspace to true if d is already in log-space
	// Other set it to false and the function will initialize value=log(d)
	LogDouble(const long double d, bool in_logspace=true);
	// Copy Constructor
	LogDouble(const LogDouble& dub);
	// Operators for manipulating Log Doubles:
	// "=","+","+=","*","*=","-","-=","/","/=","<","<=",">",">=","=="
	LogDouble& operator=(const LogDouble& other);
	void setValue(long double val);
	LogDouble operator+(const LogDouble& other) const;
	LogDouble& operator+=(const LogDouble& other);
	LogDouble operator*(const LogDouble& other) const;
	LogDouble& operator*=(const LogDouble& other);
	bool operator < (const LogDouble& other) const;
	bool operator > (const LogDouble& other) const;
	LogDouble operator/(const LogDouble& other) const;
	void printValue()
	{
		long double actVal = exp(value);
		if(actVal < numeric_limits<double>::max() && actVal > numeric_limits<double>::min())
			cout<<" (Actual) "<<actVal;
		else
			cout<<" (Log Space) "<<value;
	}	

	double getPrintableValue(bool& logspace)
	{
		if(is_zero)
			return 0;
		long double actVal = exp(value);
		if(actVal < numeric_limits<double>::max() && actVal > numeric_limits<double>::min())
		{
			return actVal;
		}
		else
		{
			logspace=true;
			return value;
		}
	}

	friend ostream& operator<<(ostream& out, const LogDouble& ld) // output
	{
		if (ld.is_zero)
			out<<0.0;
		else
			out << exp(ld.value);
	    return out;
	}

	friend istream& operator>>(istream& in, LogDouble& ld) // input
	{
	    long double x;
	    in>>x;
	    ld=LogDouble(x);
	    return in;
	}
	static void LDPower(LogDouble in,int n,LogDouble& out);
	static LogDouble Binomial(int n,int r);
	static LogDouble Factorial(int n);
};

#endif /* LOGDOUBLE_H_ */
