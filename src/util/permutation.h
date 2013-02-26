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
#ifndef PERMUTATION_H
#define PERMUTATION_H
#include "array.h"
#include "arraysaccessor.h"

/* The Permutation class iterates through all permutations of a set of
 * objects.  In order to avoid equivalent permutations, equality must be
 * defined among the objects.
 */
template<typename Type>
class Permutation
{
public:
    Permutation()
        : hasNext_(false)
    { /* NOP */ }

    Permutation(const Array<Type>& vals) 
    { setValues(vals); }


    void setValues(const Array<Type>& vals) 
    {
        values_ = vals;

        perm_.clear();
        perm_.growToSize(values_.size());
        filled_.clear();
        filled_.growToSize(values_.size());

        // Keep track of duplicates.  We don't want to double-count
        // equivalent permutations.
        prevDuplicate_.clear();
        numLaterDuplicates_.clear();
        for (int i = 0; i < values_.size(); i++) {
            int prevDuplicate = -1;
            for (int j = 0; j < i; j++) {
                if (values_[j] == values_[i]) {
                    prevDuplicate = j;
                    numLaterDuplicates_[j]++;
                }
            }
            prevDuplicate_.append(prevDuplicate);
            numLaterDuplicates_.append(0);
        }

        reset();
    }

    void reset()
    {
        skips_.clear();
        skips_.growToSize(values_.size());

        // Return to first permutation
        for (int i = 0; i < skips_.size(); i++) { 
            skips_[i] = 0;
        }

        for (int i = 0; i < values_.size(); i++) {
            perm_[i] = values_[i];
        }

        // NOTE: this could be wrong if we have 0 values.  Boundary case.
        hasNext_ = true;
    }

    ~Permutation() { /* NOP */ }

    int size() {
        return values_.size();
    }

    bool hasNext() {
        return hasNext_;
    }

    // Advance to the next permutation.
    void next() {

        // Advance to next permutation
        hasNext_ = false;
        for (int i = size() - 2; i >= 0; i--) {
            skips_[i]++;
            int numSkips = skips_[i];
            // Add up all previous skips
            int prevIndex = i;
            while (prevDuplicate_[prevIndex] >= 0) {
                numSkips += skips_[prevDuplicate_[prevIndex]];
                prevIndex = prevDuplicate_[prevIndex];
            }
            if (numSkips < size() - i - numLaterDuplicates_[i]) {
                hasNext_ = true;
                break;
            } else {
                skips_[i] = 0;
            }
        }

        // Mark all slots as empty
        for (int i = 0; i < size(); i++) {
            filled_[i] = false;
        }

        for (int i = 0; i < size(); i++) {
            int numSkips = skips_[i];
            if (prevDuplicate_[i] >= 0) {
                numSkips += skips_[prevDuplicate_[i]];
            }
            int currIndex = 0;
            while (numSkips > 0 || filled_[currIndex]) {
                if (!filled_[currIndex]) {
                    numSkips--;
                }
                currIndex++;
                assert (currIndex < size()); 
            }
            perm_[currIndex] = values_[i];
            filled_[currIndex] = true;
        }
    }

    Type operator[](int idx) {
        return perm_[idx];
    }


private:
    bool hasNext_;
    Array<int> skips_;
    Array<Type> values_;
    Array<bool> filled_;
    Array<Type> perm_;
    Array<int> prevDuplicate_;
    Array<int> numLaterDuplicates_;
};
#endif // ndef PERMUTATION_H
