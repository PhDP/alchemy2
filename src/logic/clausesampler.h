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
#ifndef CLAUSESAMPLER_H_NOV_16_2005
#define CLAUSESAMPLER_H_NOV_16_2005

#include <cmath>
#include "database.h"
#include "truefalsegroundingsstore.h"

  // expl and logl aren't available in cygwin / windows, so use exp and log
#ifndef expl
# define expl exp
# define logl log
#endif

  // test for convergence once every this many samples
const int NUM_SAMPLES_TEST_CONV = 200; 

class Clause;
class VarsGroundedType;


class ClauseSampler
{
 public:
  ClauseSampler(const double& delta, const double& epsilon,
                const int& minSamples, const int& maxSamples) 
    : delta_(delta), epsilon_(epsilon), minSamples_(minSamples),
      maxSamples_(maxSamples) { random_.init(-1); }
  
  ~ClauseSampler() {}
  
  void setEpsilon(const double& epsilon) { epsilon_ = epsilon; }
  double getEpsilon() const { return epsilon_; }

  void setDelta(const double& delta) { delta_ = delta; }
  double getDelta() const { return delta_; }

  void setMinSamples(const int& min) { minSamples_ = min; }
  int getMinSamples() const { return minSamples_; }

  void setMaxSamples(const int& max) { maxSamples_ = max; }
  int getMaxSamples() const { return maxSamples_; }


  double computeNumSamples(const int& numPreds) const
  { 
    assert(numPreds > 1);
    double n =  9/2.0 * (numPreds-1) * 1/(epsilon_*epsilon_) * logl(2.0/delta_);
    if (minSamples_ >= 0 && n < minSamples_)     n = minSamples_;
    else if (maxSamples_ >= 0 && n > maxSamples_) n = maxSamples_;
    return n;
  }

    //clause must have more than one predicate and been canonicalized
  double estimateNumTrueGroundings(Clause* const& clause, 
                                   const Predicate* const & flippedGndPred,
                                   const Domain* const& domain,
                                   double numSamples=-1);

 private:
  void getNumGroundingWithEachPredGrounded(const Clause* const & clause,
                                          Array<double>& gndingsWithPredGnded,
                                        const Array<VarsGroundedType*>& vgtArr);
  
  void getProbInfo(const Clause* const & clause, const Database* const & db, 
                   const Predicate* const & flippedGndPred,
                   const Array<VarsGroundedType*>* const & vgtArr, 
                   TrueFalseGroundingsStore* const & tfGndingsStore,
                   Array<float>& corrProbLimit, double& sumTrueGnds, 
                   Array<double>& numTrueGndsPerPos);  


  int choosePredPos(const Array<float>& corrProbLimit)
  {
    float r = random_.random();    
    for (int predPos = 0; predPos < corrProbLimit.size(); predPos++)
      if (r < corrProbLimit[predPos]) return predPos;
    return -1;
  }


  Array<int>* chooseSample(const Clause* const & clause, 
                           const Array<VarsGroundedType*>* const & vgtArr,
                           const Domain* const & domain, const int& predPos,
                           const Predicate* const & flippedGndPred);


  void groundClause(const Array<VarsGroundedType*>& vgtArr, 
                    const Array<int>& samp);

  
  int testSampleMembership(Clause* const & clause, 
                           Array<VarsGroundedType*>* const & vgtArr,
                           const Database* const & db, 
                           const Array<int>& samp,
                           const double& sumTrueGnds);


  double getNumSamplesNeeded(const double& sigma, const double& gamma,
                             const double& epsilon)
  {
    double p = (gamma+1)/2.0;
    
      // inverse of normal function
    double a[] = {-3.969683028665376e+01,  2.209460984245205e+02,
                  -2.759285104469687e+02,  1.383577518672690e+02,
                  -3.066479806614716e+01,  2.506628277459239e+00};
    double b[] = {-5.447609879822406e+01,  1.615858368580409e+02,
                  -1.556989798598866e+02,  6.680131188771972e+01,
                  -1.328068155288572e+01 };
    
    double c[] = {-7.784894002430293e-03, -3.223964580411365e-01,
                  -2.400758277161838e+00, -2.549732539343734e+00,
                  4.374664141464968e+00,  2.938163982698783e+00};
    
    double d[] = {7.784695709041462e-03, 3.224671290700398e-01,
                  2.445134137142996e+00, 3.754408661907416e+00};
    
      // Define break-points.
    double plow  = 0.02425;
    double phigh = 1 - plow;
    double invNorm;
    
    if (p < plow)  // Rational approximation for lower region:
    {
      double q = sqrt(-2*log(p));
      invNorm = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                 ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    else
    if (phigh < p) // Rational approximation for upper region:
    {
      double q  = sqrt(-2*log(1-p));
      invNorm =  -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
                  ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
    }
    else     // Rational approximation for central region:
    {
      double q = p - 0.5;
      double r = q*q;
      invNorm = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
                (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
    }
      // compute number of samples needed
    double val = invNorm * sigma / epsilon;
    return val*val;
  }


  
 private:
  double delta_;
  double epsilon_;

  int minSamples_;
  int maxSamples_;

  Random random_;
};


#endif
