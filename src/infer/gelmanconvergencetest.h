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
#ifndef GELMANCONVERGENCETEST_H_SEP_30_2005
#define GELMANCONVERGENCETEST_H_SEP_30_2005

#include "meanvariance.h"

/**
 * 
 */
class GelmanConvergenceTest
{
 public:
  GelmanConvergenceTest(const int& numChains) :
    numChains_(numChains), 
    withinChainMeanVars_(new MeanVariance[numChains]),
    numSamples_(0) {}

  ~GelmanConvergenceTest() { delete [] withinChainMeanVars_; }


    // values is array of size numChains_
  void appendNewValues(const double* const & values) 
  {
    for (int i = 0; i < numChains_; i++)
      withinChainMeanVars_[i].appendValue(values[i]);
    numSamples_++;
  }

    // values is array of size numChains_
  void appendNewValues(const bool* const & values) 
  {
    for (int i = 0; i < numChains_; i++)
      withinChainMeanVars_[i].appendValue((double)values[i]);
    numSamples_++;
  }


  double getConvergenceScore()
  {
    betweenChainsMeanVar_.reset();
    double totalW = 0;
    for (int i = 0; i < numChains_; i++) 
    {
      betweenChainsMeanVar_.appendValue( withinChainMeanVars_[i].getMean() );
      totalW += withinChainMeanVars_[i].getVariance();
    }
    int numValues = withinChainMeanVars_[0].getNumValues();

    double B = betweenChainsMeanVar_.getVariance() * numValues;
    double W = totalW / numChains_;

      // score as stated in "Probability and Statistics", DeGroot and Schervish
    double score = B/W;
    return score; 
  }
  

  int getNumSamplesAdded() { return numSamples_; }


  static bool checkConvergenceOfAll(GelmanConvergenceTest* tests[], 
                                    const int& numTests,
                                    const bool& print=false) 
  {
    // threshold as stated in "Probability and Statistics", DeGroot & Schervish 
    double threshold = 1 + 0.44 * tests[0]->getNumSamplesAdded();
    int maxItem     = -1;
    double maxScore = -1;
    int numbad      = 0;

    for (int f = 0; f < numTests; f++) 
    {
      double score = tests[f]->getConvergenceScore();

      if (!finite(score)) { numbad++; continue; }

      if (score > threshold) 
      {
        if (print)
          cout << " Item " << f << "'s score of " << score << " > threshold of "
               << threshold << endl; 
        return false;
      }

      if (score > maxScore) 
      {
        maxScore = score;
        maxItem = f;
      }
    }

    if (numbad == numTests)
    {
      if (print) cout << " All scores were inf or Nan!" << endl;
      return false;
    }
  
    // at this point all scores are less than the threshold

    if (print)
      cout << " max item is " << maxItem << " with score " << maxScore 
           << " < threshold of " << threshold << endl; 

    return true;
  }


 private:  
  int numChains_;
  MeanVariance* withinChainMeanVars_;
  int numSamples_;
  MeanVariance betweenChainsMeanVar_;

};


#endif
