#ifndef FEATURE_H_NOV_08
#define FEATURE_H_NOV_08

#define KVAL 100
#define NORM 0
#define REC_HACK 0

#define NO_NORMALIZE 1
#define SIGMOID 0
#define SMOOTH_MAX 0
#define SMOOTH_MAX2 0

#include <stdlib.h>
inline double frand() {
    return rand()/(double)RAND_MAX;
}

#include "predicate.h"
#include "predicatetemplate.h"
#include "database.h"
#include <math.h>
#include "parentiter2.h"
using namespace __gnu_cxx;

// Forward declaration
class RRF;
class GroundRRF;
class GroundFeature;
class PredicateGroundFeature;
class ClausalGroundFeature;
class ConstantGroundFeature;
class RecursiveGroundFeature;

inline double sigmoid(double x)
{
    if (x < -100.0) {
        return 0;
    } else if (x > 100.0) {
        return 1;
    } else {
       return 1.0/(1.0 + exp(-x));
    }
}


/* Feature -- an abstract class to represent an RRF feature.
 * This could reference other, child features, or be a ground predicate,
 * depending on the derived class.
 */
class Feature
{
public:
    Feature(const char* name = NULL) : name_(NULL)
    { 
        if (name != NULL) {
            setName(name);
        }
    }

    virtual ~Feature() {
        delete [] name_;
    }

    const char* getName() const {
        return name_;
    }

    void setName(const char* name) {
        delete [] name_;
        name_ = new char[strlen(name)+1]; 
        strcpy(name_, name); 
    }

    void setId(int id) {   id_ = id; }
    int  getId()       { return id_; }

    /* Compute the value of one grounding of this feature, given
     * a particular world.
     */
    virtual double getValue(const Array<int>& grounding, Database* db) {
#if 1
        return getCachedValue(grounding, db);
#else
        return computeValue(grounding, db);
#endif
    }

    virtual double getLogValue(const Array<int>& grounding, Database* db) {
#if 1
        return getCachedLogValue(grounding, db);
#else
        return computeLogValue(grounding, db);
#endif
    }

    double getCachedValue(const Array<int>& grounding, Database* db) {
        int index = getGroundingIndex(grounding, db);

        // Expand the array as necessary
        while (cacheValid_.size() <= index) {
            cacheValid_.append(false);
            cacheValue_.append(1.0);
        }

        if (!cacheValid_[index]) {
            cacheValue_[index] = computeValue(grounding, db);
            cacheValid_[index] = true;
        }

        return cacheValue_[index];
    }


    double getCachedLogValue(const Array<int>& grounding, Database* db) {
        int index = getGroundingIndex(grounding, db);

        // Expand the array as necessary
        while (cacheLogValid_.size() <= index) {
            cacheLogValid_.append(false);
            cacheLogValue_.append(1.0);
        }

        if (!cacheLogValid_[index]) {
            cacheLogValue_[index] = computeLogValue(grounding, db);
            cacheLogValid_[index] = true;
        }

        return cacheLogValue_[index];
    }

    // For caching and updating partial derivatives of each weight
    // Used by getPsuedoCounts in gfeature.h/cpp.
    double getCount(int w) {
        if (w < counts_.size()) {
            return counts_[w];
        } else {
            return 0.0;
        }
    }

    void setCount(int w, double val) {
        while (w >= counts_.size()) {
            counts_.append(0.0);
        }
        counts_[w] = val;
    }


    void invalidateAll() {
        for (int i = 0; i < cacheValid_.size(); i++) {
            cacheValid_[i] = false;
        }
        for (int i = 0; i < cacheLogValid_.size(); i++) {
            cacheLogValid_[i] = false;
        }
    }

    virtual void invalidate(const Array<int>& fgrounding, Database* db) {

        // Invalidate this grounding of this feature, and all
        // dependent groundings of its parents
        int index = getGroundingIndex(fgrounding, db);

        bool wasValid = (cacheValid_.size() > index && cacheValid_[index])
            || (cacheLogValid_.size() > index && cacheLogValid_[index]);

        if (cacheValid_.size() > index) { 
            cacheValid_[index] = false; 
        }
        if (cacheLogValid_.size() > index) { 
            cacheLogValid_[index] = false; 
        }

        if (wasValid) {
            for (int i = 0; i < parents_.size(); i++) {
                parents_[i]->invalidateChild(id_, fgrounding, db);
            }
        }
    }

    virtual void invalidateChild(int feature, const Array<int>& grounding,
            Database* db) {
        // This should never be called.  It should only be called for
        // recursive features, which override this definition.
        assert(false);
        return;
    }
    
    virtual double computeValue(const Array<int>& grounding, Database* db) = 0;
    virtual double computeLogValue(const Array<int>& grounding, Database* db) = 0;

    /* Compute the partial derivative of this feature with respect
     * to a particular weight coefficient in a particular feature.
     *
     * As in getValue(), there's no caching.
     */
    virtual double getPartialDeriv(int featureIndex, int weightIndex,
            const Array<int>& grounding, Database* db) {

        return computePartialDeriv(featureIndex, weightIndex, grounding, db);
    }

    // Access term types
    int getNumTerms() const         { return termTypes_.size(); }
    void addTermType(int type)      { termTypes_.append(type); }
    void setTermType(int idx, int type) { termTypes_[idx] = type; }
    int getTermType(int idx) const  { return termTypes_[idx]; }
    Array<int> getTermTypes() const { return termTypes_; }

    // Get and set weights (only applies to RecursiveFeatures)
    virtual int getNumWeights() const  { return 0; }
    virtual double getWeight(int idx) { assert(false); return 0.0; }
    virtual void setWeight(int idx, double weight) { assert(false); }


    /* The following two methods are used by GroundRRF for constructing
     * a tree (or graph) of ground features.
     */

    // Utility method
    int getGroundingIndex(const Array<int>& grounding, const Database* db) 
        const {

        const Domain* domain = db->getDomain();

        int index = 0;
        for (int i = 0; i < getNumTerms(); i++) {
            index *= domain->getNumConstantsByType(getTermType(i));
            index += grounding[i];
        }

        return index;
    }

    // For building a ground feature tree
    virtual GroundFeature* constructGroundFeature(GroundRRF* rrf,
            const Array<int>& grounding, Database* db) 
    { return NULL; }

    virtual void print(ostream& out) const
    { 
        // List feature name
        if (name_ == NULL) {
            out << "f" << id_ << "(";
        } else {
            out << name_ <<  "(";
        }

        // List all terms
        for (int i = 0; i < getNumTerms(); i++) {
            if (i > 0) {
                out << ",";
            }
            out << (char)('a' + i);
            out << ':' << getTermType(i);
        }
        out << ")";
    }

    void addParent(Feature* parent) {
        parents_.append(parent);
    }


protected:
    virtual double computePartialDeriv(int featureIndex, int weightIndex,
            const Array<int>& grounding, Database* db) 
    { return 0; }


protected:
    // Type for each of our terms 
    Array<int> termTypes_;

    // Unique index associated with this feature. 
    int id_;

    // Name of the feature
    char* name_;

    Array<Feature*> parents_;
    Array<bool>   cacheValid_;
    Array<double> cacheValue_;
    Array<bool>   cacheLogValid_;
    Array<double> cacheLogValue_;

    // For caching and updating partial derivatives of each weight
    // Used by getPsuedoCounts in gfeature.h/cpp.
    Array<double> counts_;
};



/* A feature whose value is simply that of a predicate.
 * Reference no other features.
 */
class PredicateFeature : public Feature
{
public:
    PredicateFeature(const PredicateTemplate* predTemplate) 
        : pred_(predTemplate)
    { 
        // Set up placeholder terms and term types
        for (int i = 0; i < predTemplate->getNumTerms(); i++) {
            pred_.appendTerm(new Term(-1));
            addTermType(predTemplate->getTermTypeAsInt(i));
        }

        setName(pred_.getName());
    }

    virtual double computeValue(const Array<int>& grounding, Database* db) 
    { 
        assert(grounding.size() == pred_.getNumTerms());

        // Substitute in constants 
        // (all terms in the predicate are independent)
        for (int i = 0; i < grounding.size(); i++) {
            pred_.setTermToConstant(i, grounding[i]);
        }
        // Query database
        return db->getValue(&pred_); 
    }

    virtual double computeLogValue(const Array<int>& grounding, Database* db)
    {
        // This shouldn't be called!
        assert(false);
        return 0.0;
    }

    virtual double getPartialDeriv(int fi, int wi, 
            const Array<int>& grounding, Database* db) 
        { return 0; }

    virtual void print(ostream& out) const {
#if 0
        Feature::print(out);
        out << " = " << pred_.getName() << "(";

        for (int i = 0; i < getNumTerms(); i++) {
            if (i > 0) {
                out << ",";
            }
            out << (char)('a' + i);
        }
        out << ")";
#endif
    }


    // EFFICIENCY HACK!
    Predicate* getPredicate() const { return &pred_; }

    inline virtual GroundFeature* constructGroundFeature(GroundRRF* rrf,
            const Array<int>& grounding, Database* db);

private:
    // mutable -- for efficiency hack above
    mutable Predicate pred_;
};



/* A feature with a fixed, numerical value.
 */
class ConstantFeature : public Feature
{
public:
    ConstantFeature(double value=1) 
        : value_(value)
    { setName(""); }

    // Return fixed value
    virtual double getValue(const Array<int>& grounding, Database* db) 
    { return value_; }

    virtual double computeValue(const Array<int>& grounding, Database* db) 
    { return value_; }

    virtual double computeLogValue(const Array<int>& grounding, Database* db) 
    { return log(value_); }

    // Derivative of a constant is always zero
    virtual double getPartialDeriv(int fi, int wi, 
            const Array<int>& grounding, Database* db) 
    { return 0; }

    virtual inline GroundFeature* constructGroundFeature(GroundRRF* rrf,
            const Array<int>& grounding, Database* db);

    virtual void print(ostream& out) const {
        /*
        Feature::print(out);
        out << " = " << value_;
        out << endl;
        */
    }

private:
    double value_;
};



/* An RRF feature whose value is a log linear function of other features.
 */
class RecursiveFeature : public Feature
{
public:
    RecursiveFeature(const char* name, bool logDerivs = false, 
            bool normalize = true) 
        : doDerivsOfLog_(logDerivs), 
#if NO_NORMALIZE
        normalize_(false), 
#else
        normalize_(normalize), 
#endif
          cachedZinvalid_(true), cachedValuesInvalid_(true)
    {
        setName(name);
    }

    // TODO: copy constructur, operator=, and destructor...

    // Get and set weights
    virtual int getNumWeights() const { return getNumChildren(); }
    virtual double getWeight(int idx) { 
        assert(idx >= 0 && idx < getNumWeights());
        return weights_[idx]; 
    }

    virtual void setWeight(int idx, double weight) { 
        assert(idx >= 0 && idx < getNumWeights());
        cachedZinvalid_ = true;
        cachedValuesInvalid_ = true;
        weights_[idx] = weight; 
    }

    void addChild(Feature* feature, double weight, const Array<int>& map) {
        children_.append(feature);
        weights_.append(weight);
        termMap_.append(map);
        feature->addParent(this);
    }

    Feature* getChild(int i) {
        return children_[i];
    }

    int getNumChildren() const { return children_.size(); }


    virtual void invalidateChild(int feature, const Array<int>& fgrounding, 
            Database* db) {

        // Hack to make the top-level work
        if (getNumTerms() == 0) {
            Array<int> nullGrounding;
            invalidate(nullGrounding, db);
            return;
        }

        const Domain* domain = db->getDomain();

        // Iterate over all groundings of this feature, to see
        // which ones are invalid.
        Array<Array<int> > uniqVals;
        for (int i = 0; i < fgrounding.size(); i++) {
            int type = getTermType(i);
            while (uniqVals.size() <= type) {
                uniqVals.append(Array<int>());
            }
            if (!uniqVals[type].contains(fgrounding[i])) {
                uniqVals[type].append(fgrounding[i]);
            }
        }
        //ParentIter iter(uniqVals[type], *(domain->getConstantsByType(1)), 
        //        getNumTerms());
        ParentIter2 iter(uniqVals, domain->getConstantsByType(), 
                getTermTypes());

        Array<int> grounding;

        // Mark each relevant grounding and all its parents as invalid
        while (iter.hasNextGrounding()) {
            iter.getNextGrounding(grounding);
            invalidate(grounding, db);
        }
    }


    virtual double computePartialDeriv(int fi, int wi, 
            const Array<int>& grounding, Database* db)
    {
        // TODO: optimize by caching...
        double totalPartial = 0.0;

        if (numGroundings_.size() == 0) {
            cacheNumGroundings(db);
        }

        // First, compute the partial derivative of the weight sum
        if (fi == id_) {

            // We own the weight; partial is simply the sum of the 
            // corresponding child values.  (i.e., d w*x/ dw = x)
            totalPartial = childGroundSum(wi, grounding, db);

            if (normalize_) {
                totalPartial -= getNorm(wi);
            }
            // getNorm(wi) results from taking the derivative of the 
            // 1/Z component.
        } else {

            // Sum partial derivatives for each grounding of each child,
            // according to the chain rule of partial derivatives.
            for (int i = 0; i < getNumChildren(); i++) 
            {

                ArraysAccessor<int>* groundingIter = 
                    getChildGroundingIter(i, grounding, db);

                Array<int> childGrounding;
                int j = 0;
                do {
                    groundingIter->getNextCombination(childGrounding);
                    double childValue = children_[i]->getPartialDeriv(fi, wi, 
                                childGrounding, db);
#if 0
                    // DEBUG
                    cout << "RRF Child " << i << "," << j <<  "," << fi 
                        << "," << wi << ": " << childValue << endl;
#endif
                    totalPartial += weights_[i] * childValue;
                    j++;
                } while (groundingIter->hasNextCombination());

#if NORM
                // Normalize
                totalPartial /= j;
#endif

                releaseChildGroundingIter(i, groundingIter);
            }
        }

        // Multiply by value to get true partial:
        // d e^f(x)/dx = e^f(x) * d f(x)/dx
        if (doDerivsOfLog_) {
            return totalPartial;
        } else {
#if SIGMOID
            double val = getValue(grounding, db);
            return val * (1.0 - val) * totalPartial;
#else
            return getValue(grounding, db) * totalPartial;
#endif
        }
    }

    virtual double computeValue(const Array<int>& grounding, Database* db) 
    {
#if SIGMOID
        double totalValue = 0.0;

        if (numGroundings_.size() == 0) {
            cacheNumGroundings(db);
        }

        for (int i = 0; i < getNumChildren(); i++) {
            totalValue += weights_[i] * childGroundSum(i, grounding, db);
        }

        return sigmoid(totalValue);
#else
        return exp(computeLogValue(grounding, db));
#endif
    }

    virtual double computeLogValue(const Array<int>& grounding, Database* db) 
    {
        double totalValue = 0.0;

        if (numGroundings_.size() == 0) {
            cacheNumGroundings(db);
        }

        for (int i = 0; i < getNumChildren(); i++) {
            totalValue += weights_[i] * childGroundSum(i, grounding, db);
        }

        return totalValue - getLogZ();
    }

    virtual inline GroundFeature* constructGroundFeature(GroundRRF* rrf,
            const Array<int>& grounding, Database* db);

    virtual void print(ostream& out) const {

        Feature::print(out);
        out << " = exp(";

        // TODO: we may want to allow special groundings like f12(A,A),
        // even when A is a free parameter.  Currently, that's broken.
        char currFreeVar = 'a' + getNumTerms();

        // Print out weighted sum of child features
        for (int i = 0; i < getNumChildren(); i++) {
            if (i > 0) {
                out << " + ";
            }
            
            assert(children_[i]->getName() != NULL);
            out << weights_[i] << " " << children_[i]->getName();
            
            if (strlen(children_[i]->getName()) > 0) {

                out << "(";
                for (int j = 0; j < termMap_[i].size(); j++) {
                    if (j > 0) {
                        out << ",";
                    }
                    if (termMap_[i][j] < 0) {
                        out << currFreeVar;
                        //currFreeVar++;
                    } else {
                        out << (char)('a' + termMap_[i][j]);
                    }
                }
                out << ")";
            }
        }
        out << ")\n";
    }

//protected:
    // Utility method: compute the value of the given child over
    // all groundings of the child consistent with the parent grounding
    // and the term mapping.
    double childGroundSum(int childIndex, const Array<int>& parentGrounding,
            Database* db)
    {
        double totalValue = 0.0;
        ArraysAccessor<int>* groundingIter = 
            getChildGroundingIter(childIndex, parentGrounding, db);

#if NORM
        // Normalize over the number of groundings
        int numGroundings = 0;
#endif
        Array<int> childGrounding;
        do {
            groundingIter->getNextCombination(childGrounding);
            totalValue += children_[childIndex]->getValue(childGrounding, db);
#if NORM
            numGroundings++;
#endif
        } while (groundingIter->hasNextCombination());

        releaseChildGroundingIter(childIndex, groundingIter);
#if NORM
        return totalValue/numGroundings;
#else
        return totalValue;
#endif
    }

    // Returns an iterator through all groundings of the specified
    // child that are consistent with the provided grounding of the 
    // parent (this).  Uses db to obtain a list of constants for each
    // type.  Must be freed by releaseChildGroundingIter().
    ArraysAccessor<int>* getChildGroundingIter(int childId, 
            const Array<int>& grounding, Database* db) 
    {
        const Domain* domain = db->getDomain();
        ArraysAccessor<int>* groundingIter = new ArraysAccessor<int>;

        // Setup object to iterate over all groundings of this child.
        // Respect mappings from parent terms to child terms;
        // all others are free to be any constant of the proper type.
        for (int term = 0; term < children_[childId]->getNumTerms(); 
                term++) {

            // TODO: Make this work with multiple types!
            if (termMap_[childId][term] < 0) {
                // Iterate over all constants of the appropriate type
                int type = children_[childId]->getTermType(term);
                groundingIter->appendArray(domain->getConstantsByType(type));
            } else {
                // Use the passed in value
                Array<int>* singleton = new Array<int>;
                singleton->append(grounding[termMap_[childId][term]]);
                groundingIter->appendArray(singleton);
            }
        }

        return groundingIter;
    }

    // Frees a grounding iterator created by the above function.
    void releaseChildGroundingIter(int childId,
            ArraysAccessor<int>* groundingIter) 
    {
        // Delete singleton arrays we created for the groundingIter
        for (int j = 0; j < children_[childId]->getNumTerms(); j++) {
            if (termMap_[childId][j] >= 0) {
                delete groundingIter->getArray(j);
            }
        }

        delete groundingIter;
    }

    double getLogZ() const {
#if SIGMOID
        return 0.0;
#else

        if (!normalize_) { 
            return 0.0; 
        }

        // This array must already be set up.  HACK-ish
        assert(numGroundings_.size() > 0);

        if (cachedZinvalid_) {
            cachedLogZ_ = 0.0;
            for (int i = 0; i < getNumChildren(); i++) {
#if SMOOTH_MAX
                // Similar to dividing by the largest possible value,
                // but softened so the derivative is continuous.
                double sigma_w = sigmoid(weights_[i]);
                // Compute sigma(weights_[i]) robustly, to avoid NaN
                cachedLogZ_ += weights_[i] * sigma_w * numGroundings_[i];
#elif SMOOTH_MAX2
                if (weights_[i] < -10.0) {
                    cachedLogZ_ += 0.0;
                } else if (weights_[i] > 10.0) {
                    cachedLogZ_ += weights_[i] * numGroundings_[i];
                } else {
                    cachedLogZ_ += numGroundings_[i]
                        * log(1.0 + exp(KVAL * weights_[i])) / KVAL;
                }
#else
                if (weights_[i] > 0.0) {
                    cachedLogZ_ += weights_[i] * numGroundings_[i];
                }
#endif
            }

            cachedZinvalid_ = false;
        }
        return cachedLogZ_;
#endif
    }

    virtual double getZ() const {
        return exp(getLogZ());
    }

    virtual double getNorm(int idx) {
#if SIGMOID
        return 0.0;
#else

        if (!normalize_) {
            return 0.0;
        }

        // This array must already be set up.  HACK-ish
        assert(numGroundings_.size() > 0);

        if (cachedValuesInvalid_) {
            cachedNormalizers_.growToSize(getNumWeights());
            for (int i = 0; i < getNumWeights(); i++) {

#if SMOOTH_MAX
                // Continuous hack to avoid the discontinuity with using
                // the maximum value
                double w_i = getWeight(i);
                double sigma_w = sigmoid(w_i);
                
                //norm = sigma_w * (1.0 - sigma_w) * numGroundings_[i];
                double norm = sigma_w * (1.0 + w_i * (1.0 - sigma_w)) 
                    * numGroundings_[i];
                cachedNormalizers_[i] = norm;
#elif SMOOTH_MAX2
                cachedNormalizers_[i] = sigmoid(getWeight(i)) 
                    * numGroundings_[i];
#else
                if (getWeight(i) < 0.0) {
                    cachedNormalizers_[i] = 0.0;
                } else {
                    cachedNormalizers_[i] = 1.0 * numGroundings_[i];
                }
#endif
            }
            cachedValuesInvalid_ = false;
        }
        return cachedNormalizers_[idx];
#endif
    }

    void cacheNumGroundings(Database* db)
    {
        const Domain* domain = db->getDomain();

        numGroundings_.clear();
        for (int i = 0; i < getNumChildren(); i++) {
            int numGroundings = 1;
            for (int term = 0; term < children_[i]->getNumTerms(); term++) {
                if (termMap_[i][term] < 0) {
                    int type = children_[i]->getTermType(term);
                    numGroundings *= domain->getNumConstantsByType(type);
                } 
            }
            numGroundings_.append(numGroundings);
        }
    }


protected:
    Array<Feature*>    children_;
    Array<double>      weights_;
    Array<ArraysAccessor<int>* > groundingIters_;

    /* Map of parent feature's terms to child feature's terms.
     *
     * The outer array indexes over all children.
     * The inner array indexes over the terms of a child.
     * A value i >= 0 maps the child's terms to the parent's ith
     * term.  A value i < 0 maps the child's term to a new 
     * term.  This has the effect of summing over all values.
     */
    Array<Array<int> > termMap_;

    // Keep track of the number of groundings for each feature.
    // TODO: this may vary if we switch databases halfway through...
    // ...how can we avoid this?
    Array<int> numGroundings_;

    // When taking derivatives, use the log of this feature value, not
    // the feature value itself.  This should be true for the top feature
    // (and only the top feature).
    bool doDerivsOfLog_; 

    bool normalize_;

    mutable bool   cachedZinvalid_;
    mutable double cachedLogZ_;
    mutable bool   cachedValuesInvalid_;
    mutable Array<double> cachedNormalizers_;
};


/* A ClausalFeature is like a RecursiveFeature
 * except it only has PredicateFeatures as its children.
 *
 * WARNING: inheritance used for code sharing, but if you
 * add non-predicate features as children of a ClausalFeature...
 * everything will break.
 */
class ClausalFeature : public RecursiveFeature
{
public:
    // BUG: segfaults when name is NULL... why?
    ClausalFeature(const char* name = NULL) 
        : RecursiveFeature(name), cachedZinvalid_(true), 
          cachedValuesInvalid_(true)
    { /* NOP */ }
        
    virtual double computePartialDeriv(int fi, int wi,
            const Array<int>& grounding, Database* db) 
    {
        if (fi != id_) {
            return 0.0;
        }
        
        // Find the ground value of the specified term.
        // There had better only be a *single* grounding, 
        // or we're in trouble... 

        double totalValue = 
            getChildValue(wi, grounding, db) - getNorm(wi);
        // getNorm(wi) results from taking the derivative of the 1/Z component.

#if SIGMOID
        double val = getValue(grounding, db);
        return val * (1.0 - val) * totalValue;
#else
        return getValue(grounding, db) * totalValue;
#endif
    }

    virtual double computeLogValue(const Array<int>& grounding, Database* db)
    {
        return computeLogValueRaw(grounding, db) - getLogZ();
    }

    virtual double computeValue(const Array<int>& grounding, Database* db) 
    {
#if SIGMOID
        double total = 0.0;
        for (int childId = 0; childId < getNumChildren(); childId++) {
            total += weights_[childId] 
                * getChildValue(childId, grounding, db);
        }

        return sigmoid(total);
#else
        double logValue = computeLogValue(grounding, db);
        // HACK to avoid underflow
        if (logValue < -100.0) {
            return 0.0;
        } else {
            return exp(logValue); 
        }
#endif
    }

    virtual void setWeight(int idx, double weight) { 
        RecursiveFeature::setWeight(idx, weight);
        cachedZinvalid_ = true;
        cachedValuesInvalid_ = true;
    }

    virtual double getNorm(int idx) {
#if SIGMOID || NO_NORMALIZE
        return 0.0;
#else
        if (cachedValuesInvalid_) {
            cachedNormalizers_.growToSize(getNumWeights());
            for (int i = 0; i < getNumWeights(); i++) {
#if 0
                // Normalize based on all possible values
                double norm = 1.0/(1.0 + exp(-getWeight(i)));
#elif SMOOTH_MAX
                // Continuous hack to avoid the discontinuity with using
                // the maximum value
                double w_i = getWeight(i);
                double norm;
                double sigma_w = sigmoid(w_i);

                //norm = sigma_w * (1.0 - sigma_w);
                norm = sigma_w * (1.0 + w_i * (1.0 - sigma_w));
#elif SMOOTH_MAX2
                // Continuous hack to avoid the discontinuity with using
                // the maximum value
                double w_i = getWeight(i);
                double norm;
                double sigma_w = sigmoid(w_i);

                norm = sigma_w;
#else
                // Normalize based on maximum possible value
                double norm;
                if (getWeight(i) < 0.0) {
                    norm = 0.0;
                } else {
                    //norm = 1.0/getZ();
                    norm = 1.0;
                }
#endif
                cachedNormalizers_[i] = norm;
            }
            cachedValuesInvalid_ = false;
        }
        return cachedNormalizers_[idx];
#endif
    }

    void setWeightsFromExample(const Array<int>& grounding, Database* db, 
            const Array<int>& queryPreds, int numWeights)
    {
        if (numWeights < 0) {
            numWeights = getNumChildren();
        }

        // Keep going until one of the query predicates has a decent weight
        bool nonZeroQueryPredicate = false;
        while (!nonZeroQueryPredicate) {

            int numNonZero = 0;
            for (int i = 0; i < getNumChildren(); i++) {

                if (frand() < ((double)numWeights - numNonZero)
                        /((double)getNumChildren() - i)) {

                    if (getChildValue(i, grounding, db)) {
                        setWeight(i, 1.0 + frand() );
                    } else {
                        setWeight(i, -1.0 - frand() );
                    }
                    
                    int featureId = children_[i]->getId();
                    if (queryPreds.find(featureId) != -1) {
                        nonZeroQueryPredicate = true;
                    }
                    numNonZero++;
                } else {
                    //setWeight(i, (frand() - 0.5)/100.0 );
                    setWeight(i, 0.0);
                }
            }
        }
    }

    virtual GroundFeature* constructGroundFeature(GroundRRF* rrf,
            const Array<int>& grounding, Database* db);

    double getLogZ() const {
#if SIGMOID || NO_NORMALIZE
        return 0.0;
#else
        if (cachedZinvalid_) {
            cachedLogZ_ = 0.0;
            for (int i = 0; i < getNumChildren(); i++) {
#if 0
                // Normalize based on all possible values
                if (weights_[i] > 100.0) {
                    // Very large weights dominate the 1.0 term
                    cachedLogZ_ += weights_[i];
                } else if (weights_[i] > -100.0) {
                    // Very tiny weights contribute nothing
                    cachedLogZ_ += log(1.0 + exp(weights_[i]));
                }
#elif SMOOTH_MAX
                // Similar to dividing by the largest possible value,
                // but softened so the derivative is continuous.
                double sigma_w = sigmoid(weights_[i]);
                // Compute sigma(weights_[i]) robustly, to avoid NaN
                cachedLogZ_ += weights_[i] * sigma_w;
#elif SMOOTH_MAX2
                if (weights_[i] < -10.0) {
                    /* NOP */
                } else if (weights_[i] > 10.0) {
                    cachedLogZ_ += weights_[i];
                } else {
                    cachedLogZ_ += log(1.0 + exp(KVAL*weights_[i]))/KVAL;
                }
#else
                // Normalize by dividing by the largest possible value
                if (weights_[i] > 0.0) {
                    cachedLogZ_ += weights_[i];
                }
#endif
            }
            cachedZinvalid_ = false;
        }
        return cachedLogZ_;
#endif
    }

    virtual double getZ() const {
        return exp(getLogZ());
    }

protected:

    double computeLogValueRaw(const Array<int>& grounding, 
            Database* db) const {

        double total = 0.0;
        for (int childId = 0; childId < getNumChildren(); childId++) {
            total += weights_[childId] 
                * getChildValue(childId, grounding, db);
        }
        return total;
    }


    double getChildValue(int childId, const Array<int>& grounding, 
            Database* db) const
    {
#if 0
        Predicate* pred = ((PredicateFeature*)children_[childId])
                              ->getPredicate();
        for (int i = 0; i < children_[childId]->getNumTerms(); i++) {
            pred->setTermToConstant(i, grounding[termMap_[childId][i]]);
        }

        // NOTE: This had better be 1 or 0, not 2 (UNKNOWN)
        return (int)db->getValue(pred);
#else
        Array<int> cgrounding;
        for (int i = 0; i < children_[childId]->getNumTerms(); i++) {
            cgrounding.append(grounding[termMap_[childId][i]]);
        }
        return children_[childId]->getValue(cgrounding, db);
#endif
    }

    mutable bool   cachedZinvalid_;
    mutable double cachedLogZ_;
    mutable bool   cachedValuesInvalid_;
    mutable Array<double> cachedNormalizers_;
};

#include "gfeature.h"

#endif
