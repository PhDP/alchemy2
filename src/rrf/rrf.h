#ifndef RRF_H
#define RRF_H
#define SIGMASQ 10000
/*
 * rrf.h -- Defines the RRF class.
 *
 * RRF stands for Recursive Random Field, a probabilistic model that
 * generalizes Markov Logic Networks (MLNs).  In an RRF,
 *
 *     P(World) = 1/Z exp( sum_i w_i f_i )
 *
 * That is, an RRF is a log linear model (just like an MLN).  The
 * difference is in how the features, f_i are defined.  In an MLN,
 * each is typically a clause, e.g.:  (X_1 v X_2 v X_4)
 * In an RRF, each feature is itself an RRF:
 *
 *     f_i(x) = 1/Z_i exp( sum_j w_j f_j )
 *
 * ...or a ground predicate...
 *   
 *     f_i(x) = Pred(x)
 *
 * We can therefore think of the entire probability distribution as
 * one top-level feature:
 *
 *     P(World) = f_0
 *
 * Each feature may have one or more terms, represented here as x.
 * Each term in a child feature may be bound to a term of the
 * parent feature, or may be left unspecified, in which case we sum over
 * all groundings.  (This is similar to summing over all groundings in an
 * MLN.)
 */

// TODO: const correctness

#include "feature.h"
#include "predicate.h"
#include "predicatetemplate.h"
#include <math.h>
#include <iostream>
using namespace __gnu_cxx;

class RRF
{
private:
    // Private, unimplemented copy-constructor so it won't get called
    // accidentally.
    RRF(const RRF& other);

public:

    RRF() { maxFeatureId_ = -1; topFeature_ = NULL; }

    ~RRF() {
        featureArray_.deleteItemsAndClear();
    }

    // To be called *after* all features have been allocated.
    // TODO -- mmake this more flexible and intelligent
    void load(istream& in, Domain* domain);

    // Create a two-level RRF with specified number of features.
    // 
    // of each arity.  The first dimension is 1.  So, for the
    // array {2, 1, 4}, this will create 2 features with arity 1
    // (\forall X), 1 of arity 2 (\forall X,Y), and 4 of arity 3
    // (\forall X,Y,Z).
    //
    // HACK: For now, we assume only a single type.
    // Caller takes ownership over returned RRF
    static RRF* generateTwoLevel(Domain* domain);

    // Add a feature with randomly selected term mappings
    void addRandomFeature(int numChildren, Array<int> queryPreds,
        Array<int> typeArity, Domain* domain);

    // Add a feature with every possible term mapping of every predicate
    void addCompleteFeature(Array<int> typeArity, Domain* domain, 
            const Array<int>& queryPreds, int numChildren);

    double getValue(Database* db) { 
        Array<int> emptyGrounding;
        // DEBUG
        //invalidateAll();
        return topFeature_->getValue(emptyGrounding, db); 
    }

    double getLogValue(Database* db) { 
        Array<int> emptyGrounding;
        // DEBUG
        //invalidateAll();
        return topFeature_->getLogValue(emptyGrounding, db); 
    }

    // TODO: should all of this be log likelihood...?
    double getExactLikelihood(Database* db) {
        return getValue(db) / getExactZ(db->getDomain());
    }

    double getExactConditionalLikelihood(Database* db, 
            const Array<int>& queryPreds) {
        return getValue(db) / getExactZ(db->getDomain(), queryPreds, db);
    }

    // Get log pseudo likelihood
    double getLogPseudoLikelihood(Database* db, const Array<int>& queryPreds);

    double getLogPseudoLikelihood(Database* db, const Array<Predicate*>& 
            queryPreds);

#if 0
    // Get log likelihood via Gibbs sampling 
    // (estimated as product of marginals)
    double getGibbsLogLikelihood(Database* db, const Array<int>& queryPreds);
#endif

    double getWeightLogLikelihood(double sigmaSq) const {

        double ll = 0.0;
        for (int i = 0; i < getNumWeights(); i++) {
            ll -= getWeight(i) * getWeight(i) / sigmaSq;
        }
        return ll;
    }

    void changedPredicate(const Predicate* pred, Database* db) {
        Array<int> grounding;
        for (int i = 0; i < pred->getNumTerms(); i++) {
            grounding.append(pred->getTerm(i)->getId());
        }
        featureArray_[pred->getId()-1]->invalidate(grounding, db);
    }

    void invalidateAll() {

        for (int i = 0; i < getNumFeatures(); i++) {
            featureArray_[i]->invalidateAll();
        }
    }

#if 0
    double getPredicateLogLikelihood(const Predicate* pred, Database* db)
    {

        const Domain* domain = db->getDomain();
        TruthValue originalValue = db->getValue(pred);

        RecursiveFeature* root = (RecursiveFeature*)topFeature_;

        double posWeightSum = 0.0;
        double negWeightSum = 0.0;
        Array<int> nullGrounding;
        for (int wi = 0; wi < topFeature_->getNumWeights(); wi++)
        {
            ArraysAccessor<int>* iter 
                = root->getChildGroundingIter(wi, nullGrounding, db);

#if 0
            Array<int> grounding;
#else
            // Prepare array once.  Efficiency hack.
            Array<int> grounding(root->getChild(wi)->getNumTerms());
            for (int i = 0; i < root->getChild(wi)->getNumTerms(); i++) {
                grounding.append(-1);
            }
#endif
            while (iter->hasNextCombination()) {
#if 0
                iter->getNextCombination(grounding);
#else
                for (int i = 0; i < grounding.size(); i++) {
                    grounding[i] = iter->getItem(i);
                }
                iter->next();
#endif
                // HACK DEBUG
                //posWeightSum += root->getWeight(wi) * 1.0;
                //continue;

                bool allTermsPresent = true;
                for (int i = 0; i < pred->getNumTerms(); i++) {
                    int id = pred->getTerm(i)->getId();
                    bool termPresent = false;
                    for (int j = 0; j < grounding.size(); j++) {
                        if (grounding[j] == id) {
                            termPresent = true;
                            break;
                        }
                    }

                    if (!termPresent) {
                        allTermsPresent = false;
                        break;
                    }
                }

                if (allTermsPresent) {
                    db->setValue(pred, TRUE);
                    posWeightSum += root->getWeight(wi) * 
                        featureArray_[domain->getNumPredicates()+wi]
                        ->computeValue(grounding, db);
                    db->setValue(pred, FALSE);
                    negWeightSum += root->getWeight(wi) * 
                        featureArray_[domain->getNumPredicates()+wi]
                        ->computeValue(grounding, db);
                }
            }
        }

        db->setValue(pred, originalValue);

        if (originalValue == TRUE) {
            return -log(1.0 + exp(negWeightSum - posWeightSum));
        } else {
            return -log(1.0 + exp(posWeightSum - negWeightSum));
        }
    }
#endif


    // Note: SLOW!  This computes Z, the partition function, by looping
    // over all 2^n possible worlds, where n is the number of ground 
    // predicates.
    double getExactZ(const Domain* domain); 

    double getExactZ(const Domain* domain, const Array<int>& queryPreds,
        Database* origDb);

    
    void getCounts(Array<double>& counts, Database* db) 
    {
        Array<int> emptyGrounding;
        counts.clear();

        // For each weight of each feature, gather counts
        for (int i = 0; i < featureArray_.size(); i++) {
            for (int j = 0; j < featureArray_[i]->getNumWeights(); j++) {
                counts.append(topFeature_->getPartialDeriv(
                        i, j, emptyGrounding, db));
            }
        }
    }

    void getPseudoCounts(Array<double>& counts, Database* db, 
            Array<Predicate*> queryPreds)
    {
        // TODO...
    }

    void setWeight(int idx, double wt) {
        // TODO: optimize?
        for (int i = 0; i < featureArray_.size(); i++) {
            int numWeights = featureArray_[i]->getNumWeights();
            if (idx < numWeights) {
                featureArray_[i]->setWeight(idx, wt);
                return;
            }
            idx -= numWeights;
        }

        // We should not reach here unless given invalid input.
        // TODO: better error reporting?
        assert(false);
    }

    double getWeight(int idx) const {
        // TODO: optimize?
        for (int i = 0; i < featureArray_.size(); i++) {
            int numWeights = featureArray_[i]->getNumWeights();
            if (idx < numWeights) {
                return featureArray_[i]->getWeight(idx);
            }
            idx -= numWeights;
        }

        // We should not reach here unless given invalid input.
        // TODO: better error reporting?
        assert(false);
        return 0;
    }

    Feature* getRoot() {
        return topFeature_;
    }

    Feature* getFeature(int idx) {
        assert(idx < featureArray_.size());
        return featureArray_[idx];
    }

    const Feature* getFeature(int idx) const {
        assert(idx < featureArray_.size());
        return featureArray_[idx];
    }

    Feature* getFeature(const char* name) {
        for (int i = 0; i < featureArray_.size(); i++) {
            if (!strcmp(featureArray_[i]->getName(), name)) {
                return featureArray_[i];
            }
        }
        return NULL;
    }

    const Feature* getFeature(const char* name) const {
        for (int i = 0; i < featureArray_.size(); i++) {
            if (!strcmp(featureArray_[i]->getName(), name)) {
                return featureArray_[i];
            }
        }
        return NULL;
    }

    int getNumFeatures() const {
        return featureArray_.size();
    }

    int getNumWeights() const {

        int numWeights = 0;
        for (int i = 0; i < featureArray_.size(); i++) {
            numWeights += featureArray_[i]->getNumWeights();
        }
        return numWeights;
    }

   /*
    void learnWeights();
    void inferMissing(Database* db);
    ...
   */

    Array<Feature*> getFeatureArray() { return featureArray_; }

protected:
    Feature* topFeature_;
    int maxFeatureId_;
    Array<Feature*> featureArray_;
};


inline ostream& operator<<(ostream& out, const RRF* rrf)
{
    out << "RRF {\n";
    for (int i = 0;
            i < rrf->getNumFeatures(); i++) {
        rrf->getFeature(i)->print(out);
    }
    out << "}\n";
    return out;
}

#endif
