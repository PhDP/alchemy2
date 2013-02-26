#ifndef PARENTITER2_H
#define PARENTITER2_H

#include "parentiter.h"

/* ParentIter2 generalizes ParentIter, but allows multiple types.
 * The purpose is to iterate over all groundings of a parfeature
 * that include a certain set of constants in the arguments.
 * The user specifies the type of each argument, the constants required
 * to appear (possibly none), and the list of valid constants in the
 * domain for each type.
 */
class ParentIter2
{
public:
    /* Constructor
     *
     * Arguments:
     *   required -- for each type, a lists of constants that must appear
     *               as arguments
     *   filler -- for each type, a list of constants that could appear
     *   types  -- the type for each argument.  The length of this array
     *             is the number of constants that will appear in the
     *             final configuration
     */
    ParentIter2(const Array<Array<int> >& required, 
            const Array<Array<int>*>* filler, const Array<int>& types)
        : firstTime_(true), size_(types.size())
    {
        // Count number of arguments of each type
        for (int i = 0; i < types.size(); i++) {
            while (numPerType_.size() <= types[i]) {
                argMap_.append(Array<int>() );
                numPerType_.append(0);
            }
            numPerType_[types[i]]++;
            argMap_[types[i]].append(i);
        }

        // Create one iterator for each type
        for (int i = 0; i < numPerType_.size(); i++) {
            piters_.append(ParentIter(required[i], 
                        (*filler)[i], numPerType_[i]));
        }
    }


    bool hasNextGrounding()
    {
        // See if there are any more configurations
        bool hasNext = false;
        for (int i = 0; i < piters_.size(); i++) {
            if (piters_[i].hasNextGrounding()) {
                hasNext = true;
                break;
            }
        }
        return hasNext;
    }


    bool getNextGrounding(Array<int>& grounding) 
    {
        grounding.clear();
        grounding.growToSize(size_);

        if (firstTime_) {

            pgrounds_.clear();
            pgrounds_.growToSize(piters_.size());

            // Get the first grounding
            for (int i = 0; i < piters_.size(); i++) {
                piters_[i].getNextGrounding(pgrounds_[i]);
            }
            firstTime_ = false;
        } else {
            // Advance to the next configuration
            for (int i = 0; i < piters_.size(); i++) {
                if (numPerType_[i] == 0) {
                    continue;
                }
                if (piters_[i].hasNextGrounding()) {
                    // Increment this configuration, and we're done
                    piters_[i].getNextGrounding(pgrounds_[i]);
                    break;
                } else {
                    // Go back to first configuration
                    piters_[i].reset();
                    piters_[i].getNextGrounding(pgrounds_[i]);
                }
            }
        }

        // Copy to the passed-in array
        for (int i = 0; i < numPerType_.size(); i++) {
            for (int j = 0; j < numPerType_[i]; j++) {
                grounding[argMap_[i][j]] = pgrounds_[i][j];
            }
        }

        return hasNextGrounding();
    }

private:
    Array<ParentIter> piters_;
    Array<Array<int> > argMap_;
    Array<Array<int> > pgrounds_;
    Array<int> numPerType_;
    bool firstTime_;
    int size_;
};

#endif // ndef PARENTITER2_H
