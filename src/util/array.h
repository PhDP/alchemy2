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
#ifndef ARRAY_H_JUN_21_2005
#define ARRAY_H_JUN_21_2005

#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include "random.h"
using namespace std;


template <typename Type> 
class Array
{
 public:
  Array(const int& initSize=0)
  {
    items_    = NULL;
    maxItems_ = 0;
    numItems_ = 0;
    if (initSize > 0) allocateMemory(initSize);
  }

  Array(const int& initSize, const Type& fill)
  {
    items_    = NULL;
    maxItems_ = 0;
    numItems_ = 0;
    growToSize(initSize, fill);
  }
  
  
    // creates an array which points to this data
    // Array owns the data
  Array(Type* data, int numdata) 
  {
    items_    = data;
    maxItems_ = numdata;
    numItems_ = numdata;
  }

  
  Array(const Array<Type>& x)
  { 
    if (x.numItems_ == 0) 
    {
      items_    = NULL;
      maxItems_ = 0;
      numItems_ = 0;
      allocateMemory(0);
      return;
    }
    items_ = new Type[x.numItems_];
    //commented out: Type may not be a basic type and thus need deep copying
    //memcpy(items_, x.items_, x.numItems_*sizeof(Type));
    for (int i = 0; i < x.numItems_; i++) items_[i] = x.items_[i];
    numItems_ = x.numItems_;
    maxItems_ = numItems_;
  }


  Array<Type>& operator=(const Array<Type>& x)
  { 
    if (items_ == x.items_) {
      return *this;
    }

    if (items_) {
      delete [] items_;
    }

    // The below code is mostly copied from the copy constructor.
    if (x.numItems_ == 0) 
    {
      items_    = NULL;
      maxItems_ = 0;
      numItems_ = 0;
      allocateMemory(0);
      return *this;
    }
    items_ = new Type[x.numItems_];
    //commented out: Type may not be a basic type and thus need deep copying
    //memcpy(items_, x.items_, x.numItems_*sizeof(Type));
    for (int i = 0; i < x.numItems_; i++) items_[i] = x.items_[i];
    numItems_ = x.numItems_;
    maxItems_ = numItems_;

    return *this;
  }


    // does not delete each individual item
  ~Array() { if (items_) delete [] items_; }


    // returns index of appended item
  int append(Type item)
  {
    if (numItems_ == maxItems_)
      allocateMemory(maxItems_*2+1);
    items_[numItems_++] = item;
    return numItems_-1;
  }


  void append(const Array<Type>& newItems)
  {
    int origSize = numItems_;
    int newItemsSize = newItems.numItems_;
    growToSize(origSize+newItemsSize);
    //commented out: Type may not be a basic type and thus need deep copying
    //memcpy(items_+origSize, newItems.items_, newItemsSize*sizeof(Type));
    for (int i = 0; i<newItemsSize;i++) items_[origSize+i] = newItems.items_[i];

  }


  void append(const Array<Type>* const & items)  { append(*items); }

 
  void copyFrom(const Array<Type>& items) { clear(); append(items); }


  bool growToSize(const int& newSize)
  {
    if (newSize <= numItems_) return false;
    
    if (newSize <= maxItems_)
      numItems_ = newSize;
    else
    {
      allocateMemory(newSize);
      numItems_ = newSize;
    }
    return true;
  }


  bool growToSize(const int& newSize, Type fill)
  {
    if (newSize <= numItems_) return false;
    int origNumItems = numItems_;
    growToSize(newSize);
    for (int i = origNumItems; i < numItems_; i++)
      items_[i] = fill;
    return true;
  }

  
  bool shrinkToSize(const int& newSize)
  {
    if (newSize >= numItems_) return false;
    allocateMemory(newSize);
    numItems_ = newSize;
    return true;
  }
  

  void clear()  { numItems_ = 0; }
  void clearAndCompress()  { numItems_ = 0; compress(); }


    // delete the contents of items_ array and set numItems = 0;
    // items_ is not deleted,and maxItems_ stay the same
  void deleteItemsAndClear()
  {
    for (int i = 0; i < numItems_; i++)
      if (items_[i]) delete items_[i];
    clear();
  }


    // delete the contents of items_ array and set numItems = 0;
    // items_ is not deleted,and maxItems_ stay the same
  void deleteItemsAndClearCompress() 
  { 
    for (int i = 0; i < numItems_; i++)
      if (items_[i]) delete items_[i];
    clear();
    compress();
  }


  int size() const { return numItems_; }

  int maxItems() const { return maxItems_; }

  bool empty() const { return size()==0; } 


  Type& item(const int& index) const 
  {
    assert(index < numItems_);
    return items_[index];
  }


  Type& operator[](const int& index) const { return item(index); } 


  Type& lastItem() const { return items_[numItems_-1]; }


    // caller should not delete the returned pointer nor modify the object
  const Type* getItems() const { return items_; }

  
  //comment out: randomOneOf no longer a static function
  //Type& randomItem() const 
  //{
  //  assert(numItems_>0); 
  //  return items_[Random::randomOneOf(numItems_)];
  //}


    // linear search
  int find(Type const & item) const
  {
    for (int i = 0; i < numItems_; i++)
      if (items_[i] == item) return i;
    return -1;
  }
  

    // linear search
  bool contains(Type const & item) const { return find(item) >= 0;}


    // Appends item only if it is not already in the array 
    // (by checking item1 == item2). Involves slow linear search.
    // Returns true of items is unique and inserted.
  bool appendUnique(Type item)
  {
    if (!contains(item))
    {
      append(item);
      return true;
    }
    return false;
  }


  void appendUnique(Array<Type>& items)
  {
    for (int i = 0; i < items.size(); i++)
      appendUnique(items.item(i));
  }


  void appendUnique(Array<Type>* items) { appendUnique(*items); }


  void removeAllNull()
  {
      // remove all entries with value NULL, and reindex so they're contiguous  
      // (e.g. [3 4 7 NULL 3 6 NULL] becomes [3 4 7 3 6]
    int currentIndex = 0;
    for (int i = 0; i < numItems_; i++)
    {
      if (item(i) != NULL)
        item(currentIndex++) = item(i);
    }
    numItems_ = currentIndex;
  }
  

  Type removeItem(const int& index)
  {
    Type removedItem = item(index);
    for (int i = index+1; i < numItems_; i++) items_[i-1] = items_[i];
    numItems_--;
    return removedItem;

    //commented out: Type may not be a basic type and thus need deep copying
      // move everything past the item back
    //int numItemsToMove = numItems_ - index - 1;
    //if (numItemsToMove > 0)
    //  memmove(&items_[index], &items_[index+1], numItemsToMove*sizeof(Type));
    //numItems_--;
  }


  Type removeLastItem() { return removeItem(numItems_-1); }


    // removes the item, but does not leave the array in the same order as 
    // before
  Type removeItemFastDisorder(const int& index) 
  {
    Type removedItem = item(index);
    item(index) = item(numItems_-1);
    numItems_--;
    return removedItem;
  }


    // resizes the array to be exactly the number of elements in it
  void compress()  { if (maxItems_ > numItems_)  allocateMemory(numItems_); }


    // reorder randomly
  void shuffle() 
  {
    for (int i = 0; i < numItems_; i++) 
    {
        // swap item i with a random item
      int swapwith = Random::randomOneOf(numItems_);
      Type tmp;
      tmp = items_[i];
      items_[i] = items_[swapwith];
      items_[swapwith] = tmp;
    }
  }


    // returns the index of largest item
  int getMaxIndex() 
  {
    assert(numItems_ > 0);
    int max = 0;
    for (int i = 1; i < numItems_; i++) 
    {
      if (items_[i] > items_[max])
        max = i;
    }
    return max;
  }


    // returns the value of largest item
  Type getMaxValue() 
  {
    assert(numItems_ > 0);
    Type max = items_[0];
    for (int i = 1; i < numItems_; i++) 
    {
      if (items_[i] > max)
        max = items_[i];
    }
    return max;
  }    


    // sort in increasing order
  void quicksort()  { quicksort(0,numItems_-1); }

    // sort in decreasing order
  void rquicksort() { rquicksort(0,numItems_-1); }
  
    // sort in increasing order
  void bubbleSort()
  {
    Type tempItem;
    for (int i = 0; i < numItems_; i++)
      for (int j = i+1; j < numItems_; j++)
        if (items_[i] > items_[j])
        {
          tempItem = items_[i];
          items_[i] = items_[j];
          items_[j] = tempItem;
        }
  }


    // sort in increasing order
  void rbubbleSort()
  {
    Type tempItem;
    for (int i = 0; i < numItems_; i++)
      for (int j = i+1; j < numItems_; j++)
        if (items_[i] < items_[j])
        {
          tempItem = items_[i];
          items_[i] = items_[j];
          items_[j] = tempItem;
        }
  }



 private:
  void truncate(const int& newNumItems) 
  {
    assert(newNumItems <= numItems_);
    numItems_ = newNumItems;
  }


    // if newSize < numItems, numItems will be shrunk to newSize
  void allocateMemory(const int& newSize)
  {
    if (newSize > 0)
    {
      if (newSize >= numItems_)
      {
        Type* tempItems = new Type[newSize];
        if (numItems_ > 0) 
        {
          //commented out: Type may not be a basic type & thus need deep copying
          //memcpy(tempItems, items_, numItems_*sizeof(Type));
          for (int i = 0; i < numItems_; i++) tempItems[i] = items_[i];
        }
        if (items_) delete [] items_;
        items_ = tempItems;
        maxItems_ = newSize;
      }
      else  // newSize < numItems_
      {
        Type* tempItems = new Type[newSize];
        if (numItems_ > 0) 
        {
          //commented out: Type may not be a basic type & thus need deep copying
          //memcpy(tempItems, items_, newSize*sizeof(Type));
          for (int i = 0; i < newSize; i++) tempItems[i] = items_[i];
        }
        if (items_) delete [] items_;
        items_ = tempItems;
        maxItems_ = newSize;
        numItems_ = newSize;
      }
    }
    else // newSize==0
    {
      if (items_) delete [] items_;
      items_ = NULL;
      maxItems_ = 0;
      numItems_ = 0;
    }
  }


  void quicksort(const int& l, const int& r)
  {
    if (l >= r) return;
    Type tmp = items_[l];
    items_[l] = items_[(l+r)/2];
    items_[(l+r)/2] = tmp;

    int last = l;
    for (int i = l+1; i <= r; i++)
      if (items_[i] < items_[l])
      {
        ++last;
        tmp = items_[last];
        items_[last] = items_[i];
        items_[i] = tmp;
      }
    
    tmp = items_[l];
    items_[l] = items_[last];
    items_[last] = tmp;
    quicksort(l, last-1);
    quicksort(last+1, r);  
  }


  void rquicksort(const int& l, const int& r)
  {
    if (l >= r) return;
    Type tmp = items_[l];
    items_[l] = items_[(l+r)/2];
    items_[(l+r)/2] = tmp;

    int last = l;
    for (int i = l+1; i <= r; i++)
      if (items_[i] > items_[l])
      {
        ++last;
        Type tmp = items_[last];
        items_[last] = items_[i];
        items_[i] = tmp;
      }
    
    tmp = items_[l];
    items_[l] = items_[last];
    items_[last] = tmp;
    rquicksort(l, last-1);
    rquicksort(last+1, r);  
  }

     
 private:
  Type* items_;
  int maxItems_;
  int numItems_;
};


// Useful functions

  // between is the filler string to be used between items
template<typename T> 
void writeArray(const Array<T*>& array, ostream& out, char* between, 
                const char& delimiter)
{
  out << array.size() << delimiter;
  for (int i = 0; i < array.size(); i++)
  {
    array.item(i)->save(out);
    if (between) out << between;
  }
}


template<typename T> 
void writeArray(const Array<T*>& array, ostream& out, char* between=NULL)
{ writeArray(array, out, between, '\n'); }


template<typename T> 
void readArray(Array<T*>& array, istream& in)
{
  int numItems;
  in >> numItems;
  for (int i = 0; i < numItems && in.good(); i++)
  {
    T* x = new T;
    x->load(in);
    array.append(x);
  }
}


template<typename T> 
void writeArray(const Array<T>& array, ostream& out, char* between, 
                const char& delimiter)
{
  out << array.Size() << delimiter;
  for (int i = 0; i < array.Size(); i++)  
  {
    array.item(i).save(out);
    if (between) out << between;
  }  
}


template<typename T> 
void writeArray(const Array<T>& array, ostream& out, char* between=NULL)
{ writeArray(array, out, between, '\n'); }


template<typename T> 
void readArray(Array<T>& array, istream& in) 
{
  int numItems;
  in >> numItems;
  for (int i = 0; i < numItems && in.good(); i++) 
  {
    T x;
    x.load(in);
    array.append(x);
  }
}


template<typename T> 
void writeBasicArray(const Array<T>& array, ostream& out, char* between,
                     const char& delimiter)
{
  out << array.Size() << delimiter;
  for (int i = 0; i < array.Size(); i++)
  {
    out << array.item(i) << " ";
    if (between) out << between;
  }
}


template<typename T> 
void writeBasicArray(const Array<T>& array, ostream& out, char* between=NULL)
{ writeBasicArray(array, out, between, '\n'); }


template<typename T> 
void readBasicArray(Array<T>& array, istream& in)
{
  int numItems=0;
  in >> numItems;
  for (int i = 0; i < numItems && in.good(); i++)
  {
    T x;
    in >> x;
    array.append(x);
  }
  if (!in.good())
    cout << "ERROR: readBasicArray(). Input not good." << endl;
}

  //Print the array
template <typename type>
ostream& printArray(const Array<type>& array, ostream& out)
{
  return printArray(array, 0, out);
}

  //Print the array
template <typename type>
ostream& printArray(const Array<type>& array, int beginIndex, ostream& out)
{
  char delimiter = ' ';
  for (int i = beginIndex; i < array.size(); i++)
  {
    out << array[i] << delimiter;
  }
  return out;
}

  //get the total number of elements in the array of arrays
template <typename type> 
int getNumArrayArrayElements(const Array<Array<type>*> & elementsArr)
{
  int numElements = 0;
  for (int i = 0; i < elementsArr.size(); i++)
  {
    numElements += elementsArr[i]->size();
  }
  return numElements;
}


#endif
