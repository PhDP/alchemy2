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
#ifndef HASHLIST_H_JUL_20_2005
#define HASHLIST_H_JUL_20_2005

#include <list>
#include <ext/hash_map>
using namespace __gnu_cxx;


  // A list that is backed up by a hash_map
template <typename Type, class HashFn, class EqualFn> 
class HashList
{
 public:
  HashList() 
    : map_(new hash_map<Type, typename list<Type>::iterator, HashFn, EqualFn>),
      list_(new list<Type>) {}
  ~HashList() { delete map_; delete list_; }


  Type& front() { return list_->front(); }
  Type& back( ) { return list_->back();  }

  typename list<Type>::iterator begin() { return list_->begin(); }
  typename list<Type>::iterator end()   { return list_->end(); }
  typename list<Type>::reverse_iterator rbegin() { return list_->rbegin(); }
  typename list<Type>::reverse_iterator rend()   { return list_->rend(); }

    // the list and its contents should not be modified
  const list<Type>* getList() const { return list_; }

  bool contains(const Type& val) const  
  { return (map_->find(val) != map_->end()); }


  Type* find(const Type& val) const
  {
    typename 
      hash_map<Type, typename list<Type>::iterator, HashFn, EqualFn>::iterator  
      mit = map_->find(val);
    
    if (mit == map_->end()) return NULL;
    return (Type*)(&((*mit).first));
  }


  typename list<Type>::iterator findListIterator(const Type& val) const
  {
    typename 
      hash_map<Type, typename list<Type>::iterator, HashFn, EqualFn>::iterator  
      mit = map_->find(val);
    
    if (mit == map_->end()) return NULL;
    return (Type*)(&((*mit).first));
  }


  void clear() { map_->clear(); list_->clear(); }


  void deleteContentsAndClear()
  {    
    typename list<Type>::iterator it = list_->begin();
    for (; it != list_->end(); it++) delete *it;
    map_->clear();
    list_->clear();
    assert(map_->size() == list_->size());
  }


  bool empty() const { return list_->empty(); }


    // returns true if val is in HashList and is erased; otherwise returns false
  bool erase(const Type& val, const bool& deleteVal=false)
  {
    typename 
      hash_map<Type, typename list<Type>::iterator, HashFn, EqualFn>::iterator  
      mit = map_->find(val);

    if (mit == map_->end()) return false;
    typename list<Type>::iterator lit = (*mit).second;
    list_->erase(lit);
    map_->erase(mit);
    assert(map_->size() == list_->size());
    if (deleteVal) delete *lit;
    return true;
  }


  bool eraseAndDelete(const Type& val)  { return erase(val, true); }


    // returns true if the val is inserted, returns false otherwise
  bool insert(typename list<Type>::iterator loc, const Type& val)
  {
    if (contains(val)) return false;
    list_->insert(loc, val);
    (*map_)[val] = --loc;
    assert(map_->size() == list_->size());
    return true;
  }
  

  int maxSize() const 
  {
    if (map_->max_size() < list_->max_size()) return map_->max_size();
    return list_->max_size();
  }


  void popBack()
  {
    Type& val = list_->back();
    list_->pop_back();
    typename 
      hash_map<Type,typename list<Type>::iterator, HashFn, EqualFn>::iterator 
      mit = map_->find(val);
    assert(mit != map_->end());
    map_->erase(mit);
    assert(map_->size() == list_->size());
  }


  void popFront()
  {
    Type& val = list_->front();
    list_->pop_front();
    typename 
      hash_map<Type,typename list<Type>::iterator, HashFn, EqualFn>::iterator 
      mit = map_->find(val);
    assert(mit != map_->end());
    map_->erase(mit);
    assert(map_->size() == list_->size());
  }


  bool pushBack(const Type& val)
  {
    if (contains(val)) return false;
    list_->push_back(val);
    typename list<Type>::iterator lit = list_->end();
    (*map_)[val] = --lit;
    assert(map_->size() == list_->size());
    return true;
  }


  bool pushFront(const Type& val)
  {
    if (contains(val)) return false;
    list_->push_front(val);
    (*map_)[val] = list_->begin();
    assert(map_->size() == list_->size());
    return true;
  }
  
  
  void reverse() { list_->reverse(); }

  int size() { assert(map_->size() == list_->size()); return list_->size(); }


 private:
    // note that there are two copies of Type, one in map_ and one in list_
  hash_map<Type, typename list<Type>::iterator, HashFn, EqualFn>* map_;
  list<Type>* list_;
};


#endif
