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
#include "random.h"
#include "math.h"

bool ExtRandom::deviate_available = false;

double ExtRandom::second_deviate = 0;

double ExtRandom::uniformRandom()
{
  double u;
  u = (double)rand() / (double)RAND_MAX;
  return u;
}

double ExtRandom::gaussRandom(double mu, double sd) 
{
    // Ref. Numerical Recipes in C - Chapter 7.2: Gaussian Deviates.
  if (deviate_available)
  {
    deviate_available = false;
    return mu + second_deviate * sd;
  }
  else
  {
    double v1, v2, rsq;
    do
    {
      v1 = 2.0 * uniformRandom() - 1.0;
      v2 = 2.0 * uniformRandom() - 1.0;
      rsq = v1*v1 + v2*v2;
    }
    while (rsq >= 1.0 || rsq == 0.0);
    double fac = sqrt(-2.0*log(rsq)/rsq);
    second_deviate = v1 * fac;
    deviate_available = true;
    return mu + v2 * fac * sd;
  }
}

double ExtRandom::expRandom(double lambda) 
{
    // Ref. Numerical Recipes in C - Chapter 7.2: Exponential Deviates.
  double num;
  do
    num = uniformRandom();
  while (num == 0.0);
  return -log(num) / lambda;
}

double ExtRandom::ComputeGauss(double mu, double sigma, double x)
{
  double re;
  re = x - mu;
  re *= re;
  re = -re/(2*sigma*sigma);
  re = exp(re);

  re = re/(sigma * sqrt(2*PI));

  return re;
}


double ExtRandom::ComputeLnGauss(double mu, double sigma, double v1, double v2,
                                 double x)
{
  double re;
  //double factor;

  if ( x < v1 || x > v2 ) // out of range
  {
    return -BigValue;
  }

  // factor = GaussianIntegral(mu, sigma,v1,v2);  // normalization factor
  // re =  ComputeGauss(mu,sigma,x);
  // re = re / factor;
  re = x - mu;
  re *= re;
  re = re / (2 * sigma * sigma);

  re = -re;
  //re = - log ( sigma * sqrt(2 * PI)) - re;

  return re;	
}

double ExtRandom::GaussianIntegral(double mu, double sigma, double v1,
                                   double v2)
{
  double re = 0;
  double x;
  double delta;
  delta = (v2-v1)/INTEGRALSTEP;
  int i;
  for ( i = 0; i < INTEGRALSTEP; i++ )
  {
    x = v1 + (double)i * delta;
    re = re + ComputeGauss(mu, sigma,x) * delta;
  }
  return re;
}

void ExtRandom::GaussianParaLearning(double &mu, double &sigma, double &v1,
                                     double &v2, vector<double> Tdata)
{
  int NumOfTData = Tdata.size();
  int i;
  vector<double>::iterator fiter = Tdata.begin();
  v1 = Tdata[0];
  v2 = Tdata[0];

  mu = 0;
  for (i = 0; i < NumOfTData; i++)
  {
    mu += *(fiter+i);
    if (v1 > Tdata[i])
    {
      v1 = Tdata[i];
    }

    if (v2 < Tdata[i])
    {
      v2 = Tdata[i];
    }
  }
  mu /= NumOfTData;

  sigma = 0;
  for (i = 0; i < NumOfTData; i++)
  {
    sigma += pow(*(fiter+i)-mu, 2);
  }
  sigma /= (double)(NumOfTData-1);
  sigma = sqrt(sigma);
}
