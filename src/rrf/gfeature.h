#ifndef GFEATURE_H_APR_06
#define GFEATURE_H_APR_06

#include "feature.h"

class GroundFeature
{
public:
    GroundFeature()
        : dirtyDeriv_(true), dirtyValue_(true), dirtyLogValue_(true)
    { /* NOP */ }

    virtual ~GroundFeature() { /* NOP */ }

    double getValue() {

        if (dirtyValue_) {
            cachedValue_ = computeValue();
            dirtyValue_ = false;
        }

        return cachedValue_;
    }

    double getLogValue() {

        if (dirtyLogValue_) {
            cachedLogValue_ = computeLogValue();
            dirtyLogValue_ = false;
        }
        
        return cachedLogValue_;
    }

    double getCounts(int featureIndex, int weightIndex) {
        return computePartialDeriv(featureIndex, weightIndex);
    }

    virtual double computeValue() { return 0.0; }
    virtual double computeLogValue() { return 0.0; }
    virtual double computePartialDeriv(int featureIndex, int weightIndex) 
        { return 0.0; }

    inline double getDeriv();

    // This is for computing the partial derivative of a feature with
    // respect to a child feature.  It's only relevant for recursive features.
    virtual double getChildDeriv(int i) {
        assert(false);
        return 0.0;
    }

    // Clear cached value/partial derivative
    inline void setDirty();
    
    inline void updateParents(double oldValue, double newValue);
    inline void updateParentCounts(double oldValue, double newValue);

    void addParent(RecursiveGroundFeature* parent, int parentIndex) {
        parents_.append(parent);
        parentIndices_.append(parentIndex);
    }
    // TODO: remove parent?  clear parent?

protected:
    bool dirtyDeriv_;
    bool dirtyValue_;
    bool dirtyLogValue_;
    double cachedDeriv_;
    double cachedValue_;
    double cachedLogValue_;
    Array<RecursiveGroundFeature*> parents_;
    Array<int> parentIndices_;
};


class RecursiveGroundFeature : public GroundFeature
{
public:
    RecursiveGroundFeature(RecursiveFeature* featureTemplate, 
            bool doDerivsOfLog=false)
        : featureTemplate_(featureTemplate), doDerivsOfLog_(doDerivsOfLog)
    { /* NOP */ }

    virtual ~RecursiveGroundFeature() { /* NOP */ }

    virtual double computeValue()
    {
#if SIGMOID
        computeLogValue();
        return sigmoid(cachedSum_);
#else
        return exp(computeLogValue());
#endif
    }

    virtual double computeLogValue()
    {
        cachedSum_ = 0.0;
        subSums_.growToSize(children_.size());

        for (int i = 0; i < children_.size(); i++) {
            double currFeatureTotal = 0.0;
            for (int j = 0; j < children_[i].size(); j++) {
                //cout << children_[i][j]->getLogValue() << endl;
                //cout << children_[i][j]->getValue() << endl;
                currFeatureTotal += children_[i][j]->getValue();
            }
#if NORM
            cachedSum_ += featureTemplate_->getWeight(i) * currFeatureTotal 
                / children_[i].size(); // Normalize by the number of groundings
#else
            cachedSum_ += featureTemplate_->getWeight(i) * currFeatureTotal;
#endif
            subSums_[i] = currFeatureTotal;
        }

        RecursiveFeature* recFeature = (RecursiveFeature*)featureTemplate_;
        // DEBUG
        //cout << featureTemplate_->getId() << endl;
        //cout << "cachedSum_ = " << cachedSum_ << endl;
        return cachedSum_ - recFeature->getLogZ();
    }

    virtual void update(int idx, double oldValue, double newValue)
    {
        dirtyDeriv_ = true;
        // TODO -- what's the right way to do this?  Possible bug.
        if (!dirtyValue_ || !dirtyLogValue_) {
            double oldCachedValue = cachedValue_;
#if NORM
            cachedSum_ += featureTemplate_->getWeight(idx) 
                * (newValue - oldValue)
                /children_[idx].size(); // Normalize by number of children
#else
            cachedSum_ += featureTemplate_->getWeight(idx) 
                * (newValue - oldValue);
#endif
            subSums_[idx] += newValue - oldValue;

            RecursiveFeature* recFeature = (RecursiveFeature*)featureTemplate_;
            cachedLogValue_ = cachedSum_ - recFeature->getLogZ();
#if SIGMOID
            cachedValue_ = sigmoid(cachedSum_);
#else
            cachedValue_ = exp(cachedLogValue_);
#endif
            updateParents(oldCachedValue, cachedValue_);
        }
    }

    virtual void updateCounts(int idx, double oldValue, double newValue)
    {
        if ((dirtyValue_ && dirtyLogValue_) || dirtyDeriv_) {
           cout << "ERROR: cache must be up-to-date when calling "
               << "updateCounts!\n";
           return;
        }

        double counts = featureTemplate_->getCount(idx);
        //int id = featureTemplate_->getId();

        // Figure out what the partial derivative with respect
        // to the idx'd weight used to be
        double totalPartial = 0.0;
        {
#if 0
            // HACK DEBUG -- verify that our nubmers are consistent here
            double totalValue = 0.0;
            for (int w = 0; w < children_.size(); w++) {
                double sum = 0.0;
                for (int j = 0; j < children_[w].size(); j++) {
                    sum += children_[w][j]->getValue();
                }
                totalValue += sum * featureTemplate_->getWeight(w);
            }
            cout << "Old:      " << cachedSum_ << endl;
            cout << "New:      " << totalValue << endl;
            totalValue += (oldValue - newValue) 
                * featureTemplate_->getWeight(idx);
            cout << "Computed: " << totalValue << endl;
            // END CONSISTENCY CHECK
#endif

            RecursiveFeature* recTemplate 
                = (RecursiveFeature*)featureTemplate_;
            totalPartial = subSums_[idx] - recTemplate->getNorm(idx);
            
            if (!doDerivsOfLog_) {
#if SIGMOID
                totalPartial *= getValue() * (1.0 - getValue());
#else
                totalPartial *= getValue();
#endif
            }
            //cout << "totalPartial(sum2) = " << totalPartial << endl;
        }

        // DEBUG
        //cout << "deriv = " << cachedDeriv_ << endl;
        //cout << "oldcounts = " << counts << endl;

        // Use it to subtract old counts 
        // (using old derivative and old partial derivative here)
        counts -= cachedDeriv_ * totalPartial;

        // DEBUG
        //cout << "newcounts = " << counts << endl;

        // Update value
        double oldCachedValue = cachedValue_;
        cachedSum_ += featureTemplate_->getWeight(idx) 
            * (newValue - oldValue);
        subSums_[idx] += newValue - oldValue;

        RecursiveFeature* recFeature = (RecursiveFeature*)featureTemplate_;
        cachedLogValue_ = cachedSum_ - recFeature->getLogZ();
#if SIGMOID
        cachedValue_ = sigmoid(cachedSum_);
#else
        cachedValue_ = exp(cachedLogValue_);
#endif
        
        // Update parent values (and derivatives)
        updateParentCounts(oldCachedValue, cachedValue_);

        // Update our derivative (depends on parent derivs)
        cachedDeriv_ = 0.0;
        for (int i = 0; i < parents_.size(); i++) {
            cachedDeriv_ += parents_[i]->getDeriv()
                * parents_[i]->getChildDeriv(parentIndices_[i]);
        }

        if (parents_.size() == 0.0) {
            cachedDeriv_ = 1.0;
        }

        // Add in updated counts
        {
            RecursiveFeature* recTemplate 
                = (RecursiveFeature*)featureTemplate_;
            totalPartial = subSums_[idx] - recTemplate->getNorm(idx);
            
            if (!doDerivsOfLog_) {
#if SIGMOID
                totalPartial *= getValue() * (1.0 - getValue());
#else
                totalPartial *= getValue();
#endif
            }
        }

        counts += cachedDeriv_ * totalPartial;
        featureTemplate_->setCount(idx,counts);
    }


    virtual double getChildDeriv(int i)
    {
        double totalPartial = featureTemplate_->getWeight(i);
#if NORM
        totalPartial /= children_[i].size();
#endif

        if (doDerivsOfLog_) {
            return totalPartial;
        } else {
            return getValue() * totalPartial;
        }
    }

    virtual void addChild(int childIndex, GroundFeature* child)
    {
        while (children_.size() <= childIndex) {
            children_.append(Array<GroundFeature*>());
        }

        children_[childIndex].append(child);
        child->addParent(this, childIndex);
    }

    virtual double computePartialDeriv(int fi, int wi) 
    {
        double totalPartial = 0.0;

        if (fi == featureTemplate_->getId()) {

            if (dirtyLogValue_ && dirtyValue_) {
                computeLogValue();
            }

            RecursiveFeature* recTemplate 
                = (RecursiveFeature*)featureTemplate_;
            totalPartial = subSums_[wi] - recTemplate->getNorm(wi);
#if NORM
            totalPartial /= children_[wi].size(); // Normalize by number of groundings
#endif
        } else {

            // We should no longer be using this code!
            assert(false);

            // Sum partial derivatives for each grounding of each child,
            // according to the chain rule of partial derivatives.
            for (int i = 0; i < children_.size(); i++) 
            {
                double featureTotal = 0.0;
                for (int j = 0; j < children_[i].size(); j++) {
                    featureTotal += children_[i][j]->getCounts(fi, wi);
                }
#if NORM
                totalPartial += featureTotal * featureTemplate_->getWeight(i)
                    / children_[i].size(); // Normalize by number of groundings
#else
                totalPartial += featureTotal * featureTemplate_->getWeight(i);
#endif
            }
        }

        if (doDerivsOfLog_) {
            return totalPartial;
        } else {
#if SIGMOID
            return getValue() * (1.0 - getValue()) * totalPartial;
#else
            return getValue() * totalPartial;
#endif
        }
    }

protected:
    Array<Array<GroundFeature*> > children_;
    RecursiveFeature* featureTemplate_;
    bool doDerivsOfLog_;
    double cachedSum_;
    Array<double> subSums_;
};


class ClausalGroundFeature : public RecursiveGroundFeature
{
public:
    ClausalGroundFeature(ClausalFeature* featureTemplate)
        : RecursiveGroundFeature(featureTemplate)
    { /* NOP */ }

    virtual ~ClausalGroundFeature() { /* NOP */ }

#if 0
    virtual double computeLogValue() {
        // TODO -- fix this!
        assert(false);
    }
#endif

    // TODO -- implement computeLogValue()

#if 0 // TODO: add subSums_ to this, or something.
    virtual double computeValue() {

        cachedSum_ = 0.0;
        for (int i = 0; i < children_.size(); i++) {
            cachedSum_ += featureTemplate_->getWeight(i) 
                * children_[i][0]->getValue();
        }
#if SIGMOID
        return sigmoid(cachedSum_);
#else
        ClausalFeature* clausalTemplate = (ClausalFeature*)featureTemplate_;
        return exp(cachedSum_ - clausalTemplate->getLogZ());
#endif
    }
#endif
    
#if 0
    virtual double getChildDeriv(int i)
    {
        // NOTE: this should never actually get called, since the only 
        // children are predicate features that have no parameters.
        // But if there were children with paramters, this is how the
        // method should be implemented.
        assert(false);
        double weight_i = featureTemplate_->getWeight(i);
        double ret = getValue() * (weight_i - 
                - 1.0/(1.0 + exp(-children_[i][0]->getValue())));
        return ret;
    }

    virtual double computePartialDeriv(int fi, int wi) 
    {
        ClausalFeature* clausalTemplate = (ClausalFeature*)featureTemplate_;

        if (fi != clausalTemplate->getId()) {
            return 0.0;
        }

        assert(children_[wi].size() == 1);
#if SIGMOID
        return getValue() * (1.0 - getValue()) 
            * children_[wi][0]->getValue();
#else
        return getValue() * (children_[wi][0]->getValue() - 
                clausalTemplate->getNorm(wi));
#endif
    }

    virtual void update(int idx, double oldValue, double newValue)
    {
        dirtyDeriv_ = true;
        if (!dirtyValue_) {
            double oldCachedValue = cachedValue_;
#if 0
            double oldCachedSum = cachedSum_;
#endif
            cachedSum_ += featureTemplate_->getWeight(idx) 
                * (newValue - oldValue);
#if SIGMOID
            cachedValue_ = sigmoid(cachedSum_);
#else
            ClausalFeature* clausalTemplate = 
                (ClausalFeature*) featureTemplate_;
            cachedValue_ = exp(cachedSum_ - clausalTemplate->getLogZ());
#endif
#if 0
            double ourCachedSum_ = cachedSum_;
            computeValue();
            if (fabs(ourCachedSum_ - cachedSum_) > 0.00001) {
                cout << "Sums are different!\n";
                cout << ourCachedSum_ << endl;
                cout << cachedSum_ << endl;
                cout << "Old: " << oldCachedSum << endl;
                cout << oldValue << " -> " << newValue << endl;
                cout << featureTemplate_->getWeight(idx) << endl;
            }
#endif
            updateParents(oldCachedValue, cachedValue_);
        } 
    }
#endif
};


class PredicateGroundFeature : public GroundFeature
{
public:
    PredicateGroundFeature(const Predicate& pred)
        : groundPred_(pred), predValue_(false)
    { /* NOP */ }

    virtual ~PredicateGroundFeature() { /* NOP */ }

    Predicate* getPredicate() { return &groundPred_; }
    virtual double computeValue() { return predValue_ ? 1.0 : 0.0; }
    void setValue(bool val)   
    { 
        if (predValue_ != val) {
            predValue_ = val; 
            setDirty();  
        }
    }

    void setValueAndUpdate(bool val)
    {
        if (predValue_ != val) {
            predValue_ = val;

            if (val) {
                cachedValue_ = 1.0;
                updateParents(0.0, 1.0);
            } else {
                cachedValue_ = 0.0;
                updateParents(1.0, 0.0);
            }
        }
    }

    void setValueAndUpdateCounts(bool val)
    {
        if (predValue_ != val) {
            predValue_ = val;

            if (val) {
                cachedValue_ = 1.0;
                updateParentCounts(0.0, 1.0);
            } else {
                cachedValue_ = 0.0;
                updateParentCounts(1.0, 0.0);
            }
        }
    }

    virtual double computePartialDeriv(int featureIndex, int weightIndex) 
        { return 0; }

private:
    Predicate groundPred_;
    bool predValue_;
};


inline GroundFeature* PredicateFeature::constructGroundFeature(
        GroundRRF* rrf, const Array<int>& grounding, Database* db)
{
    PredicateGroundFeature* ret = new PredicateGroundFeature(pred_);
    ret->setValue(computeValue(grounding, db));
    return ret;
}

class ConstantGroundFeature : public GroundFeature
{
public:
    ConstantGroundFeature(double value)
        : val_(value)
    { /* NOP */ }

    virtual ~ConstantGroundFeature() { /* NOP */ }

    virtual double computeValue() { return val_; }
    virtual double computePartialDeriv(int featureIndex, int weightIndex) 
        { return 0; }

protected:
    double val_;
};

inline GroundFeature* ConstantFeature::constructGroundFeature(
        GroundRRF* rrf, const Array<int>& grounding, Database* db) 
{
    return new ConstantGroundFeature(value_);
}

class GroundRRF
{
public:
    GroundRRF(RRF* rrf, Database* db);

    double getValue() { return root_->getValue(); }

    double getLogValue() { return root_->getLogValue(); }

    // DEBUG
    //double getLogPseudoLikelihood(const Array<int>& queryPreds, const Domain* domain);
    double getLogPseudoLikelihood(const Array<int>& queryPreds);
    double getLogPseudoLikelihood(const Array<Predicate*>& queryPreds);

    void getCounts(Array<double>& counts/*, Database* db*/);

    void getPseudoCounts(Array<double>& counts, const Array<int>& queryPreds,
            double samplingFrac);
    void getPseudoCountsFast(Array<double>& counts, 
            const Array<int>& queryPreds, double samplingFrac);

    int getNumPredicateGroundings(int predIdx) {
        // NOTE: assuming that predicate 0 is equality and therefore skipped,
        // and that predicates appear first.  This could be wrong.
        return allFeatures_[predIdx-1].size();
    }

    bool getPredicateValue(int predIdx, int groundIdx) {
        // NOTE: assuming that predicate 0 is equality and therefore skipped,
        // and that predicates appear first.  This could be wrong.
        if (allFeatures_[predIdx-1][groundIdx] != NULL) {
            return (((PredicateGroundFeature*)allFeatures_[predIdx-1][groundIdx])
                ->getValue() != 0.0);
        } else {
            return false;
        }
    }

    void setPredicateValue(int predIdx, int groundIdx, bool value) {
        // NOTE: assuming that predicate 0 is equality and therefore skipped,
        // and that predicates appear first.  This could be wrong.

        if (allFeatures_[predIdx-1][groundIdx] != NULL) {
            ((PredicateGroundFeature*)allFeatures_[predIdx-1][groundIdx])
                ->setValue(value);
        }
    }

    void setPredicateAndUpdate(int predIdx, int groundIdx, bool value) {
        // NOTE: assuming that predicate 0 is equality and therefore skipped,
        // and that predicates appear first.  This could be wrong.
        if (allFeatures_[predIdx-1][groundIdx] != NULL) {
            ((PredicateGroundFeature*)allFeatures_[predIdx-1][groundIdx])
                ->setValueAndUpdate(value);
                //->setValue(value);
        }
    }

    void setPredicateAndUpdateCounts(int predIdx, int groundIdx, bool value) {
        // NOTE: assuming that predicate 0 is equality and therefore skipped,
        // and that predicates appear first.  This could be wrong.
        if (allFeatures_[predIdx-1][groundIdx] != NULL) {
            ((PredicateGroundFeature*)allFeatures_[predIdx-1][groundIdx])
                ->setValueAndUpdateCounts(value);
        } else {
            // HACK DEBUG
            cout << "Feature was NULL! Possible bug.\n";
        }
    }

    GroundFeature* getGroundFeature(Feature* feature, Array<int>& grounding)
    {
        int featureIdx = feature->getId();
        int groundIdx  = feature->getGroundingIndex(grounding, db_);

        while (allFeatures_[featureIdx].size() <= groundIdx) {
            allFeatures_[featureIdx].append((GroundFeature*)NULL);
        }

        if (allFeatures_[featureIdx][groundIdx] == NULL) {
            allFeatures_[featureIdx][groundIdx] 
                = feature->constructGroundFeature(this, grounding, db_);
        }

        return allFeatures_[featureIdx][groundIdx];
    }

    int getNumGroundings(int featureId) const {
        return allFeatures_[featureId].size();
    }

    void dirtyAll() {
        for (int i = 0; i < allFeatures_.size(); i++) {
            for (int j = 0; j < allFeatures_[i].size(); j++) {
                if (allFeatures_[i][j] != NULL) {
                    allFeatures_[i][j]->setDirty();
                }
            }
        }
    }

private:
    GroundFeature* root_;
    RRF* rrf_;
    Database* db_;
    int numCounts_;

    Array<Array<GroundFeature*> > allFeatures_;
};


inline GroundFeature* RecursiveFeature::constructGroundFeature(
        GroundRRF* rrf, const Array<int>& grounding, Database* db) 
{
    RecursiveGroundFeature* ret = new RecursiveGroundFeature(
            this, doDerivsOfLog_);

    if (numGroundings_.size() == 0) {
        cacheNumGroundings(db);
    }

    // Add all children (recurses as necessary)
    for (int i = 0; i < getNumChildren(); i++) {

        Array<int> childGrounding;
        ArraysAccessor<int>* groundingIter = 
            getChildGroundingIter(i, grounding, db);

        do {
            groundingIter->getNextCombination(childGrounding);
            ret->addChild(i, rrf->getGroundFeature(children_[i], 
                        childGrounding));
        } while (groundingIter->hasNextCombination());

        releaseChildGroundingIter(i, groundingIter);
    }

    return ret;
}

inline GroundFeature*  ClausalFeature::constructGroundFeature(
        GroundRRF* rrf, const Array<int>& grounding, Database* db) 
{
    ClausalGroundFeature* ret = new ClausalGroundFeature(this);

    // Add all children (recurses as necessary)
    for (int i = 0; i < getNumChildren(); i++) {

        Array<int> childGrounding;
        ArraysAccessor<int>* groundingIter = 
            getChildGroundingIter(i, grounding, db);

        do {
            groundingIter->getNextCombination(childGrounding);
            ret->addChild(i, rrf->getGroundFeature(children_[i], 
                        childGrounding));
        } while (groundingIter->hasNextCombination());

        releaseChildGroundingIter(i, groundingIter);
    }

    return ret;
}


double GroundFeature::getDeriv() 
{
    if (dirtyDeriv_) {
        if (parents_.size() == 0) {
            cachedDeriv_ = 1.0;
        } else {
            cachedDeriv_ = 0.0;
            for (int i = 0; i < parents_.size(); i++) {
                cachedDeriv_ += parents_[i]->getDeriv()
                    * parents_[i]->getChildDeriv(parentIndices_[i]);
            }
        }
        dirtyDeriv_ = false;
    }

    return cachedDeriv_;
}

void GroundFeature::setDirty()
{ 
    if (!dirtyDeriv_ || !dirtyValue_ || !dirtyLogValue_) {
        dirtyDeriv_ = true; 
        dirtyValue_ = true;
        dirtyLogValue_ = true;
        for (int i = 0; i < parents_.size(); i++) {
            parents_[i]->setDirty();
        }
    }
}

void GroundFeature::updateParents(double oldValue, double newValue)
{
    for (int i = 0; i < parents_.size(); i++) {
        parents_[i]->update(parentIndices_[i], oldValue, newValue);
    }
}

void GroundFeature::updateParentCounts(double oldValue, double newValue)
{
    for (int i = 0; i < parents_.size(); i++) {
        parents_[i]->updateCounts(parentIndices_[i], oldValue, newValue);
    }
}

#endif // ndef GFEATURE_H_APR_06
