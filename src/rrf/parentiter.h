#ifndef PARENTITER_H
#define PARENTITER_H

#include "array.h"
#include "permutation.h"
#include "arraysaccessor.h"

class ParentIter
{
public:
    ParentIter()
    {
        /* NOP */
    }

    ~ParentIter()
    {
        // Default destructor should work...
    }

    ParentIter(const ParentIter& other)
        : size_(other.size_), fillerValues_(other.fillerValues_),
          reqValues_(other.reqValues_), reqGround_(other.reqGround_),
          fillerGround_(), firstTime_(other.firstTime_)
    {
        setupFillerIter();
        while (fillerGround_.getCombinationIndex() <
                other.fillerGround_.getCombinationIndex()) {
            fillerGround_.next();
        }
    }

    // TODO: expand this to allow for objects of different types
    ParentIter(const Array<int>& required, const Array<int>* filler, int size)
        : size_(size), fillerValues_(*filler), reqValues_(required),
          reqGround_(), firstTime_(true)
    {
        Array<int> requiredWithSpaces(required);
        for (int i = required.size(); i < size; i++) {
            requiredWithSpaces.append(-2);
        }
        reqGround_.setValues(requiredWithSpaces);
        setupFillerIter();
    }

    void reset()
    {
        reqGround_.reset();
        setupFillerIter();
        firstTime_ = true;
    }

    bool hasNextGrounding() {
        return fillerGround_.hasNextCombination() || reqGround_.hasNext();
    }


    bool getNextGrounding(Array<int>& grounding)
    {
        grounding.clear();
        grounding.growToSize(size_);

        if (firstTime_ || !fillerGround_.hasNextCombination()) {
            reqGround_.next();
            setupFillerIter();
        }
        firstTime_ = false;

        // Advance to next grounding
        Array<int> fillerGrounding;
        fillerGround_.getNextCombination(fillerGrounding);

        // Copy to the passed-in array
        int fillerIndex = 0;
        for (int i = 0; i < size_; i++) {
            if (reqGround_[i] < 0) {
                grounding[i] = fillerGrounding[fillerIndex++];
            } else {
                grounding[i] = reqGround_[i];
            }
        }

        return hasNextGrounding();
    }


protected:
    void setupFillerIter()
    {
        // There's some trickiness in here to avoid duplicate groundings.
        // Consider a parent grounding that must have variable A,
        // but also has one free variable that can be anaything.  If we
        // considerthe permutations (A,?) and (?,A), we must be careful
        // not to double-count (A,A) when iterating through all values of
        // the wildcard, '?'.  Therefore, this code ensures that all
        // wildcards that take on a specific value must *precede* all
        // non-wildcards of that value.  This makes groundings unique.
        Array<int> seenSoFar;
        fillerGround_.deleteArraysAndClear();
        for (int i = 0; i < size_; i++) {
            if (reqGround_[i] >= 0) {
                seenSoFar.append(reqGround_[i]);
            } else {
                Array<int>* valList = new Array<int>;
                for (int fi = 0; fi < fillerValues_.size(); fi++) {
                    bool seen = false;
                    for (int si = 0; si < seenSoFar.size(); si++) {
                        if (seenSoFar[si] == fillerValues_[fi]) {
                            seen = true;
                            break;
                        }
                    }
                    if (!seen) {
                        valList->append(fillerValues_[fi]);
                    }
                }

                fillerGround_.appendArray(valList);
            }
        }
    }


private:
    int size_;
    Array<int> fillerValues_;
    Array<int> reqValues_;
    Permutation<int> reqGround_;
    ArraysAccessor<int> fillerGround_;
    bool firstTime_;
};

#endif // ndef PARENTITER_H
