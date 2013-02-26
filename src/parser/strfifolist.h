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
#ifndef STRFIFOLIST_H_JUN_26_2005
#define STRFIFOLIST_H_JUN_26_2005

#include <list>
#include <ostream>
using namespace std;


  // Stores the tokens returned by follex. The first token is the latest one.
class StrFifoList
{
 public: 
  StrFifoList() { fifoList_ = new list<const char*>; }

  ~StrFifoList() 
  {
    list<const char*>::iterator it = fifoList_->begin();
    for (; it != fifoList_->end(); it++)  delete [] (*it);
    delete fifoList_; 
  }

  
  void deleteContentsAndClear()
  {
    list<const char*>::iterator it = fifoList_->begin();
    for (; it != fifoList_->end(); it++)  delete [] (*it);
    fifoList_->clear();
  }


    // Caller is responsible for deleting str.
  void add(const char* const & str)
  {
    //commented out: don't impose a max
    //if (fifoList_->size() == STRFIFOLIST_MAXSIZE) fifoList_->pop_back();
    char* s = new char[strlen(str)+1];
    strcpy(s,str);
    fifoList_->push_front(s);
  }


  void addLast(const char* const & str)
  {
    //commented out: don't impose a max
    //if (fifoList_->size() == STRFIFOLIST_MAXSIZE) fifoList_->pop_back();
    char* s = new char[strlen(str)+1];
    strcpy(s,str);
    fifoList_->push_back(s);
  }


    // Caller should not delete returned const char*
  const char* getItem(const int& idx)
  {
    if ((unsigned int)idx > fifoList_->size()-1) return NULL;
    list<const char*>::iterator it = fifoList_->begin();
    int i = 0;
    for (; it != fifoList_->end(); it++)
      if (i++ == idx) return (*it);
    return NULL;
  }
  

  unsigned int getSize() const { return fifoList_->size(); }

    // the caller should not modify or delete the returned pointer
  const char* back() const { return fifoList_->back(); }

    // caller is responsible to deleting returned pointer
  const char* removeLast() 
  {
    const char* last = fifoList_->back();
    fifoList_->pop_back();
    return last;
  }


  ostream& print(ostream& out) const 
  {
    list<const char*>::iterator it = fifoList_->begin();
    for (; it != fifoList_->end(); it++)
      out << (*it) << " ";
    return out;
  } 
  
 private:
  list<const char*>* fifoList_;
  
};


inline
ostream& operator<<(ostream& out, const StrFifoList& l) { return l.print(out); }


#endif
