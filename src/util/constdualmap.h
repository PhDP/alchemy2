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
#ifndef CONSTDUALMAP_H_JUN_21_2005
#define CONSTDUALMAP_H_JUN_21_2005


#include <limits>
#include <ext/hash_map>
//using namespace __gnu_cxx;
#include "array.h"
#include "strint.h"


  // StrInt stores the constant name & its type id
typedef hash_map<const StrInt*, int, HashStrInt, EqualStrInt> StrIntToIntMap;


  // Maps int to StrInt* and vice versa.
class ConstDualMap
{  
 public:
  ConstDualMap()  
  {
    intToStrIntArr_ = new Array<const StrInt*>;
    strIntToIntMap_ = new StrIntToIntMap;
  }


  ~ConstDualMap() 
  {
    for (int i = 0; i < intToStrIntArr_->size(); i++)
      delete (*intToStrIntArr_)[i];
    delete intToStrIntArr_;
    delete strIntToIntMap_;
  }


    // Returns const char* corresponding to i or NULL if there is no such char*.
    // The returned const char* should not be deleted.
  const char* getStr(const int& i)
  { 
    if (0 <= i && i < intToStrIntArr_->size())
      return (*intToStrIntArr_)[i]->str_;
    return NULL;
  }


    // Returns the ints_ of StrInt corresponding to i; returns NULL if there is
    // no such StrInt. Should not be deleted by caller.
  Array<int>* getInt2(const int& i)
  { 
    if (0 <= i && i < intToStrIntArr_->size())
      return (*intToStrIntArr_)[i]->ints_;
    return NULL;
  }

    // Returns the ints_ of StrInt corresponding to i; returns NULL if there is
    // no such StrInt. Caller should delete s if required, but not returned Array.
  Array<int>* getInt2(const char* const & s)
  { 
    int i = getInt(s);
    if (i < 0) return NULL;
    return getInt2(i);
  }


    // Returns int corresponding to str or -1 if there is no such str.
    // Caller should delete s if required.
    // Making this function const causes the compiler to complain.
  int getInt(const char* const & s)
  {
    StrInt str(s);
    StrIntToIntMap::iterator it;
    if ((it=strIntToIntMap_->find(&str)) == strIntToIntMap_->end())
      return -1;
    return (*it).second;
  }


    // Returns corresponding int (which increases by one each time addType() is 
    // called), or -1 is str has been added before.
    // Caller should delete s if required.
  int insert(const char* const & s, const int& ii)
  {
      // See if string already present
    Array<int>* oldInts = getInt2(s);
    if (oldInts != NULL)
    {
      if (!oldInts->contains(ii))
        oldInts->append(ii);
      StrInt str(s);
      StrIntToIntMap::iterator it;
      it = strIntToIntMap_->find(&str);
      return (*it).second;
    }

    StrInt* strInt = new StrInt(s,ii);

    StrIntToIntMap::iterator it;
    if ((it=strIntToIntMap_->find(strInt)) != strIntToIntMap_->end())
    {
      cout << "Warning: In ConstDualMap::insert(), tried to insert duplicate " 
           << strInt->str_ << ", prev id " << (*it).second << endl;
      delete strInt;
      return -1;
    }
    
    if (((int)intToStrIntArr_->size()) >= numeric_limits<int>::max()) 
    {
      cout << "Error: In ConstDualMap::insert(), reach int max limit when "
           << "inserting " << strInt->str_ << endl;
      delete strInt;
      exit(-1);
    }

    intToStrIntArr_->append(strInt);
    int i = intToStrIntArr_->size()-1;
    (*strIntToIntMap_)[strInt] = i;
    return i;
  }

  
  int getNumInt() const  { return intToStrIntArr_->size(); }

    // Caller should not delete the returned Array<StrInt*>*.
  const Array<const StrInt*>* getIntToStrIntArr() const  
  { 
    return intToStrIntArr_; 
  }

    // Caller should delete the returned Array<const char*>* but not its 
    // contents
  const Array<const char*>* getIntToStrArr() const  
  { 
    Array<const char*>* a = new Array<const char*>;
    for (int i = 0; i < intToStrIntArr_->size(); i++)
      a->append((*intToStrIntArr_)[i]->str_);
    return a; 
  }


    // Caller should delete the returned Array<int>*.
/*
  const Array<int>* getIntToInt2Arr() const  
  { 
    Array<int>* a = new Array<int>;
    for (int i = 0; i < intToStrIntArr_->size(); i++)
      a->append((*intToStrIntArr_)[i]->int_);
    return a; 
  }
*/  
  void compress() { intToStrIntArr_->compress(); }
 
 private:
  Array<const StrInt*>* intToStrIntArr_;
  StrIntToIntMap* strIntToIntMap_;
};

#endif
