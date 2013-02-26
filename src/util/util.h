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
#ifndef UTIL_H_OCT_17_2005
#define UTIL_H_OCT_17_2005

#include <string>
using namespace std;
#include <sstream>
#include <sys/time.h>
#include <sys/resource.h>

class Util
{
 public:
	
  static bool substr(const string& pred, unsigned int& cur, string& outstr,
                     const char* const & delimiter)
  {
    outstr.clear();
    while (cur < pred.length() && isspace(pred.at(cur))) cur++;
    string::size_type dlt = pred.find(delimiter, cur);
    if (dlt == string::npos) return false;
    unsigned int cur2 = dlt-1;
    while (cur2 >= 0 && isspace(pred.at(cur2))) cur2--;
    if (cur2 < cur) return false;
    outstr = pred.substr(cur, cur2+1-cur);
    cur = dlt+1;
    return true;
  }


  static string trim(const string& str)
  {
    if (str.empty()) return str;
    int lt = 0; 
    while (lt < (int)str.length() && isspace(str.at(lt))) lt++;
    int rt = str.length()-1;
    while (rt >= 0 && isspace(str.at(rt))) rt--;
    return str.substr(lt, rt+1-lt);
  }
  
  /* this function converts an integer to a string */
  static string intToString(int i) {
	 ostringstream myStream;
	 myStream<<i;
	 return myStream.str();
  }

  /* this function converts a double to a string */
  static string doubleToString(double i) {
	 ostringstream myStream;
	 myStream<<i;
	 return myStream.str();
  }

  static int factorial(int n)
  {
    int f = 1;
    for (int i = 1; i <= n; i++)
      f = f*i;
    return f;
  }

  static int permute(int n, int k)
  {
    int result = 1;
    if (n == 0 || k == 0)
      return result;
    for (int i = n; i > (n-k); i--)
      result = result*i;
    return result;
  }

  static double elapsed_seconds()
  {
    double answer;

    static struct rusage prog_rusage;
    static long prev_rusage_seconds = 0;
    static long prev_rusage_micro_seconds = 0;

    getrusage(0, &prog_rusage);
    answer = (double)(prog_rusage.ru_utime.tv_sec-prev_rusage_seconds)
           + ((double)(prog_rusage.ru_utime.tv_usec-prev_rusage_micro_seconds))
           / 1000000.0 ;
    prev_rusage_seconds = prog_rusage.ru_utime.tv_sec ;
    prev_rusage_micro_seconds = prog_rusage.ru_utime.tv_usec ;
    return answer;
  }
  
};

#endif
