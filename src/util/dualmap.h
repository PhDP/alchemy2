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
#ifndef DUALMAP_H_JUN_21_2005
#define DUALMAP_H_JUN_21_2005
  
  
#include <limits>
#include <ext/hash_map>
using namespace __gnu_cxx;
#include "array.h"
#include "equalstr.h"


typedef hash_map<const char*, int, hash<const char*>, EqualStr> StrToIntMap;


  // Maps int to const char* and vice versa.
class DualMap
{ 
 public:
  DualMap() : intToStrArr_(new Array<const char*>), 
              strToIntMap_(new StrToIntMap) {}

  
  DualMap(const DualMap& dm)
  { 
    intToStrArr_ = new Array<const char*>; 
    strToIntMap_ = new StrToIntMap;
    for (int i = 0; i < dm.intToStrArr_->size(); i++)
      insert((*(dm.intToStrArr_))[i]);
    compress();
  }
  

  ~DualMap() 
  {
    for (int i = 0; i < intToStrArr_->size(); i++)
      delete [] (*intToStrArr_)[i];
    delete intToStrArr_;
    delete strToIntMap_;
  }


    // Returns const char* corresponding to i or NULL if there is no such char*.
    // The returned const char* should not be deleted.
  const char* getStr(const int& i) const
  { 
    if (0<= i && i < intToStrArr_->size())
      return (*intToStrArr_)[i];
    return NULL;
  }


    // Returns int corresponding to str or -1 if there is no such str.
    // Caller should delete str if required.
    // Making this function const causes the compiler to complain.
  int getInt(const char* const & str) const
  {
    StrToIntMap::iterator it;
    if ((it=strToIntMap_->find(str)) == strToIntMap_->end())
      return -1;
    return (*it).second;
  }


    // Returns corresponding int (which increases by one each time addType() is 
    // called), or -1 is str has been added before.
    // Caller should delete str if required.
  int insert(const char* const & str)
  {
    StrToIntMap::iterator it;
    if ((it=strToIntMap_->find(str)) != strToIntMap_->end())
    {
      cout << "Warning: In DualMap::insert(), tried to insert duplicate " 
           << str << ", prev id " << (*it).second << endl;
      return -1;
    }
    
    if (((int)intToStrArr_->size()) >= numeric_limits<int>::max()) 
    {
      cout << "Error: In DualMap::insert(), reach int max limit when inserting "
           << str << endl;
      exit(-1);
    }

    char* s = new char[strlen(str)+1];
    strcpy(s,str);
    intToStrArr_->append(s);
    int i = intToStrArr_->size()-1;
    (*strToIntMap_)[s] = i;
    return i;
  }

  
  int getNumInt() const  { return intToStrArr_->size(); }

    // Caller should not delete the returned Array<const char*>*.
  const Array<const char*>* getIntToStrArr() const  { return intToStrArr_; }
  
  void compress() { intToStrArr_->compress(); }


 private:
  Array<const char*>* intToStrArr_;
  StrToIntMap* strToIntMap_;
};

#endif
