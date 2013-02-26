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
#include <iostream>
#include <cmath>
#include "random.h"
#include "array.h"
class Complex
{

public:

	Complex() : _real(0), _imag(0) {}
	explicit Complex( double r) : _real(r), _imag(0) {}
	Complex(double r, double i) : _real(r), _imag(i) {}

	Complex& operator+=(const double& d)
	{
		_real += d;
		return *this;
	}

	Complex& operator+=(const Complex& c)
	{
		_real += c._real;
		_imag += c._imag;
		return *this;
	}

	Complex& operator-=(const double &d)
	{
		_real -= d;
		return *this;
	}

	Complex& operator-=(const Complex& c)
	{
		_real -= c._real;
		_imag -= c._imag;
		return *this;
	}

	Complex& operator*=(const double& d)
	{
		_real *= d;
		_imag *= d;
		return *this;
	}

	Complex& operator*=(const Complex& c)
	{
		double re = _real;
		double im = _imag;
		_real = re * c._real - im * c._imag;
		_imag = re * c._imag + im * c._real;
		return *this;
	}

	Complex& operator/=(const double& d)
	{
		_real /= d;
		_imag /= d;
		return *this;
	}

	Complex& operator/=(const Complex& c)
	{
		double re = _real;
		double im = _imag;
		double d = c._real * c._real + c._imag * c._imag;
		_real = (re * c._real + im * c._imag) / d;
		_imag = (im * c._real - re * c._imag) / d;
		return *this;
	}

	Complex Conj() const
	{
		return Complex(_real, -_imag);
	}

	Complex TransToRadiusAngle() const
	{
		if (_real == 0 && _imag == 0)
		{
			return Complex(0,0);
		}

		double r = sqrt(_real*_real + _imag*_imag);
		double sinx = _imag / r;
		double alpha = asin(sinx); // alpha will be in -pi/2 to pi/2
		if (alpha >=0) // means in 1st or 2nd quardrant
		{
			if (_real < 0)
			{
				alpha = PI - alpha;
			}
		}
		else 
		{
			if (_real < 0)
			{
				alpha = -alpha + PI;
			}
			else
			{
				alpha = 2*PI + alpha;
			}
		}

		return Complex(r, alpha);		
	}

	Array<Complex> Sqrt() const
	{
		//
		Array<Complex> roots;

		if (_real == 0 && _imag == 0)
		{				  
			roots.append(Complex(0,0));
			return roots;
		}

		// transform to radius and angle form 

		Complex ra = TransToRadiusAngle();
		double rRoot = sqrt(ra.Real());
		

		for(int i = 0; i < 2; i++)
		{
			double angle = (ra.Imag() + 2*double(i)*PI)/2;
			double real = rRoot * cos(angle);
			double imag = rRoot * sin(angle);
			roots.append(Complex(real,imag)); 
		}

		return roots;				
	}

	Array<Complex> Root(int k) const
	{
		Array<Complex> roots;
		if (_real == 0 && _imag == 0)
		{
			roots.append(Complex(0,0));
			return roots;
		}
		Complex ra = TransToRadiusAngle();
		double rRoot = sqrt(ra.Real());
		

		for(int i = 0; i < k; i++)
		{
			double angle = (ra.Imag() + 2*double(i)*PI)/double(k);
			double real = rRoot * cos(angle);
			double imag = rRoot * sin(angle);
			roots.append(Complex(real,imag)); 
		}

		return roots;				
	}
	double Real() const { return _real; }
	double Imag() const { return _imag; }
	void Real(const double& re) { _real = re ; }
	void Imag(const double& im) { _imag = im ; }
	void Set(const double& re, const double& im){ _real = re; _imag = im ; }
	double Modsq() const { return _real*_real + _imag * _imag ; }
	double Mod() const { return sqrt(_real*_real + _imag * _imag); }

private:
	double _real;
	double _imag;

};

inline Complex operator+(const Complex& c)
{
	return Complex(c.Real(), c.Imag());
}

inline Complex operator-(const Complex& c)
{
	return Complex(-c.Real(), -c.Imag());
}

inline Complex operator+(const Complex& c, const double& d)
{
	return Complex(c.Real() + d, c.Imag());
}

inline Complex operator+(const double& d, const Complex& c)
{
	return Complex(d + c.Real(), c.Imag());
}

inline Complex operator+(const Complex& c1, const Complex& c2)
{
	return Complex(c1.Real() + c2.Real(), c1.Imag() + c2.Imag());
}

inline Complex operator-(const Complex& c, const double& d)
{
	return Complex(c.Real() - d, c.Imag());
}

inline Complex operator-(const double& d, const Complex& c)
{
	return Complex(d - c.Real(), -c.Imag());
}

inline Complex operator-(const Complex& c1, const Complex& c2)
{
	return Complex(c1.Real() - c2.Real(), c1.Imag() - c2.Imag());
}


inline Complex operator*(const Complex& c, const double& d)
{
	return Complex(c.Real() * d, c.Imag() * d);
}

inline Complex operator*(const double& d, const Complex& c)
{
	return Complex(c.Real() * d, c.Imag() * d);
}

inline Complex operator*(const Complex& c1, const Complex& c2)
{
	double real = c1.Real() * c2.Real() - c1.Imag() * c2.Imag();
	double imag = c1.Real() * c2.Imag() + c1.Imag() * c2.Real();
	return Complex(real, imag);
}

inline Complex operator/(const Complex& c, const double& d)
{
	return Complex(c.Real() / d, c.Imag() / d);
}

inline Complex operator/(const double& d, const Complex& c)
{
	double dd = c.Real() * c.Real() + c.Imag() * c.Imag();
	return Complex((d * c.Real())/dd, (-d * c.Imag())/dd);
}

inline Complex operator/(const Complex& c1, const Complex& c2)
{
	double d = c2.Real() * c2.Real() + c2.Imag() * c2.Imag();
	double real = (c1.Real() * c2.Real() + c1.Imag() * c2.Imag()) / d;
	double imag = (c1.Imag() * c2.Real() - c1.Real() * c2.Imag()) / d;
	return Complex(real, imag);
}

inline double real(const Complex &c)
{
	return c.Real();
}

inline double imag(const Complex &c)
{
	return c.Imag();
}

inline double abs(const Complex &c)
{
	return sqrt(c.Real() * c.Real() + c.Imag() * c.Imag());
}

inline double norm(const Complex &c)
{
	return c.Real() * c.Real() + c.Imag() * c.Imag();
}

inline Complex conj(const Complex &c)
{
	return Complex(c.Real(), -c.Imag());
}

inline ostream &operator<<(ostream &os, const Complex &c)
{
	os<<c.Real()<<"+"<<c.Imag()<<"i";
	return os;
} 
