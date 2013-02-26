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
#ifndef MULTIDARRAY_H_SEP_12_2005
#define MULTIDARRAY_H_SEP_12_2005

#include "array.h"


template <typename Type> 
class MultDArray
{
 public:
    // Caller should delete dim if required
  MultDArray(const Array<int>* const & dim)
  {
    arr_ = new Array<Type>;
    int n = 1;
    for (int i = 0; i < dim->size(); i++) n *= (*dim)[i];
    arr_->growToSize(n);

    multiplier_ = new Array<int>;
    for (int i = 0; i < dim->size(); i++)
    {
      n /= (*dim)[i];
      multiplier_->append(n);
    }
  }


  ~MultDArray() 
  { 
    if (multiplier_) delete multiplier_; 
    if (arr_) delete arr_; 
  }    


  const Array<Type>* get1DArray() { return arr_; }
  

  Type getItem(const Array<int>* const & indexes) const
  { return (*arr_)[getIndex(indexes)]; }


  void setItem(const Array<int>* const & indexes, const Type& item)
  { (*arr_)[getIndex(indexes)] = item; }


  void addItem(const Array<int>* const & indexes, const Type& item)
  { (*arr_)[getIndex(indexes)] += item; }


 private:
  int getIndex(const Array<int>* const & indexes) const
  {
    assert(indexes->size() == multiplier_->size());
    int idx = 0;
    for (int i = 0; i < indexes->size(); i++)
      idx += (*indexes)[i] * (*multiplier_)[i];
    return idx;
  }


 private:
  Array<int>* multiplier_;
  Array<Type>* arr_;

};

#endif
