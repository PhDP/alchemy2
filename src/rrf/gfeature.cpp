#include "gfeature.h"
#include "rrf.h"

GroundRRF::GroundRRF(RRF* rrf, Database* db)
    : rrf_(rrf), db_(db) 
{
    numCounts_ = 0;
    int numFeatures = rrf_->getNumFeatures();
    for (int i = 0; i < numFeatures; i++) {
        allFeatures_.append(Array<GroundFeature*>());
        numCounts_ += rrf_->getFeature(i)->getNumWeights();
    }

    Array<int> nullGrounding;
    root_ = rrf->getRoot()->constructGroundFeature(this, nullGrounding, db);
    allFeatures_[rrf->getRoot()->getId()].append(root_);
}


void GroundRRF::getCounts(Array<double>& counts)
{
    counts.growToSize(numCounts_);

    // For each weight of each feature, gather counts
    int countIdx = 0;
    for (int i = 0; i < allFeatures_.size(); i++) {
        for (int j = 0; j < rrf_->getFeature(i)->getNumWeights(); j++) {

            double totalCounts = 0.0;
            for (int k = 0; k < allFeatures_[i].size(); k++) {
                if (allFeatures_[i][k] == NULL) {
                    continue;
                }
                totalCounts += allFeatures_[i][k]->getDeriv()
                    * allFeatures_[i][k]->computePartialDeriv(i,j);

#if 0
                cout << "Deriv " << i << " " << k << ": ";
                cout << allFeatures_[i][k]->getDeriv() << endl;
                cout << allFeatures_[i][k]->computePartialDeriv(i,j) << endl;
#endif
            }
#if 0
            double oldCounts = root_->getCounts(i,j);
            if (2.0*fabs((oldCounts - totalCounts)/(oldCounts + totalCounts))
                    > 0.01) {
                cout << "Counts differ " << i << " " << j << ": \n";
                cout << "Groundings: " << allFeatures_[i].size() << endl;
                cout << "Deriv: " << allFeatures_[i][0]->getDeriv() << endl;
                cout << "Partial: " << allFeatures_[i][0]->computePartialDeriv(i,j) << endl;
                cout << oldCounts << endl;
                cout << totalCounts << endl;
            }
                    
            counts.append(root_->getCounts(i, j));
#endif
            counts[countIdx++] = totalCounts;
        }
    }
}


void GroundRRF::getPseudoCounts(Array<double>& counts, 
        const Array<int>& queryPreds, double samplingFrac)
{
    // Set-up counts array
    counts.growToSize(numCounts_);
    for (int i = 0; i < numCounts_; i++) {
        counts[i] = 0.0;
    }

    // Get statistics for the true values; we'll need them for reference
    double trueSum = getLogValue();

    if (isnan(trueSum) || isinf(trueSum)) {
        for (int i = 0; i < numCounts_; i++) {
            counts[i] = trueSum;
        }
        return;
    }

    Array<double> trueCounts;
    getCounts(trueCounts);

    double falseSum;
    Array<double> falseCounts;

    // Sum the effect from each ground predicate
    for (int q = 0; q < queryPreds.size(); q++) {

        int i = queryPreds[q];
        for (int j = 0; j < getNumPredicateGroundings(i); j++) {

            if (frand() < samplingFrac) {

            // Corrupt the value of a single ground predicate
            bool origValue = getPredicateValue(i,j);
            setPredicateAndUpdate(i,j,!origValue);

            // Compute its contribution to the overall gradient
            falseSum = getLogValue();
            getCounts(falseCounts);

            double prob_false = sigmoid(falseSum - trueSum);
            for (int c = 0; c < numCounts_; c++) {
                counts[c] += prob_false * (trueCounts[c] - falseCounts[c]);
            }

            // Reset the predicate's value
            setPredicateAndUpdate(i,j,origValue);
            }
        }
    }

#if 0
    for (int c = 0; c < numCounts_; c++) {
        cout << "lpl counts[" << c << "] = " << counts[c] << endl;
    }
#endif

    // Renormalize, so that we're independent of the sampling fraction
    for (int i = 0; i < numCounts_; i++) {
        counts[i] /= samplingFrac;
    }

    // Renormalize, because pll is actually *average* pll
    double numQueryPreds = 0.0;
    for (int q = 0; q < queryPreds.size(); q++) {
        numQueryPreds += getNumPredicateGroundings(queryPreds[q]);
    }
    for (int i = 0; i < numCounts_; i++) {
        counts[i] /= numQueryPreds;
    }
}


void GroundRRF::getPseudoCountsFast(Array<double>& counts, 
        const Array<int>& queryPreds, double samplingFrac)
{
#if 0 // THIS IS BROKEN -- don't use it!
    // Set-up counts array
    counts.growToSize(numCounts_);
    for (int i = 0; i < numCounts_; i++) {
        counts[i] = 0.0;
    }

    // Get statistics for the true values; we'll need them for reference
    double trueSum = getLogValue();
    Array<double> trueCounts;
    getCounts(trueCounts);

    double falseSum;
    Array<double> falseCounts;
    Array<double> falseCounts2;
#if 0
    Array<double> counts1;
    Array<double> counts2;
    counts1.growToSize(numCounts_);
    counts2.growToSize(numCounts_);
    for (int i = 0; i < numCounts_; i++) {
        counts1[i] = 0.0;
        counts2[i] = 0.0;
    }
#endif

    // Sum the effect from each ground predicate
    for (int q = 0; q < queryPreds.size(); q++) {

        int i = queryPreds[q];
        for (int j = 0; j < getNumPredicateGroundings(i); j++) {

            if (frand() < samplingFrac) {

            int c = 0;
#if 0
            for (int f = 0; f < allFeatures_.size(); f++) {
                Feature* feat = rrf_->getFeature(f);
                for (int w = 0; w < feat->getNumWeights(); w++) {
                    cout << "pre[" << c <<  "] = " 
                        << feat->getCount(w) << endl;
                    feat->setCount(w,0.0);
                    c++;
                }
            }
#endif

            // Corrupt the value of a single ground predicate
            bool origValue = getPredicateValue(i,j);
            setPredicateAndUpdateCounts(i,j,!origValue);

            // Compute its contribution to the overall gradient
            falseSum = getLogValue();
            double prob_false = sigmoid(falseSum - trueSum);

            // c indexes into the count array.
            c = 0;
            for (int f = 0; f < allFeatures_.size(); f++) {
                Feature* feat = rrf_->getFeature(f);
                for (int w = 0; w < feat->getNumWeights(); w++) {
                    // DEBUG
                    if (c == 2) {
                        cout << "pfalse = " << prob_false << endl;
                        cout << "-featCount = " << -feat->getCount(w) << endl;
                    }
                    counts[c] += prob_false * -feat->getCount(w);
#if 0
                    //cout << (trueCounts[c] - feat->getCount(w)) 
                    //    << " " << falseCounts[c] << " " << falseCounts2[c] << endl;

                    cout << prob_false*(trueCounts[c] - falseCounts[c]) << " " 
                        << prob_false*(trueCounts[c] - falseCounts2[c]) << endl;
#endif
                    c++;
                }
            }

#if 0
            // DEBUG
            c = 0;
            for (int f = 0; f < allFeatures_.size(); f++) {
                Feature* feat = rrf_->getFeature(f);
                for (int w = 0; w < feat->getNumWeights(); w++) {
                    cout << "mid[" << c <<  "] = " 
                        << feat->getCount(w) << endl;
                    c++;
                }
            }
#endif

            // Reset the predicate's value
            //setPredicateAndUpdateCounts(i,j,origValue);
            setPredicateAndUpdateCounts(i,j,origValue);
#if 0
            c = 0;
            for (int f = 0; f < allFeatures_.size(); f++) {
                Feature* feat = rrf_->getFeature(f);
                for (int w = 0; w < feat->getNumWeights(); w++) {
                    cout << "post[" << c <<  "] = " 
                        << feat->getCount(w) << endl;
                    //feat->setCount(w,0.0);
                    c++;
                }
            }
#endif
            }
        }
    }

    for (int c = 0; c < numCounts_; c++) {
        cout << "lpl1 counts[" << c << "] = " << counts[c] << endl;
    }

#if 0
    for (int c = 0; c < numCounts_; c++) {
        cout << c << ": " << counts1[c] << " ; " << counts2[c] << endl;
    }
#endif

    // Renormalize, so that we're independent of the sampling fraction
    for (int i = 0; i < queryPreds.size(); i++) {
        counts[i] /= samplingFrac;
    }
#endif
}


double GroundRRF::getLogPseudoLikelihood(const Array<Predicate*>& queryPreds)
{
    // LPL = log pseudo-likelihood, 
    //       the conditional probability of each query predicate 
    //       conditioned on all data.
    double lpl = 0.0;
    double logTruthProb = getLogValue();

    for (int q = 0; q < queryPreds.size(); q++) {

        int i = queryPreds[q]->getId();
        Array<int> grounding;
        for (int ti = 0; ti < queryPreds[q]->getNumTerms(); ti++) {
            grounding.append(queryPreds[q]->getTerm(ti)->getId());
        }
        int j = rrf_->getFeature(i)->getGroundingIndex(grounding, db_);

        if (allFeatures_[i][j] == NULL) {
            lpl += log(0.5);
        } else {
            double currValue = allFeatures_[i][j]->getValue();
            double newValue = (currValue == 0.0) ? 1.0 : 0.0;
            //setPredicateValue(i, j, newValue);
            setPredicateAndUpdate(i, j, newValue);
            double logUntruthProb = getLogValue();
            //setPredicateValue(i, j, currValue);
            setPredicateAndUpdate(i, j, currValue);
            
            if (fabs(logUntruthProb - logTruthProb) > 100) {
                lpl += (logUntruthProb - logTruthProb);
            } else {
                lpl += -log(1.0 + exp(logUntruthProb - logTruthProb));
            }
            //lpl += logTruthProb - log(exp(logTruthProb) + exp(logUntruthProb));
        }
    }

    return lpl/queryPreds.size();
}


double GroundRRF::getLogPseudoLikelihood(const Array<int>& queryPreds)
{
    // LPL = log pseudo-likelihood, 
    //       the conditional probability of each query predicate 
    //       conditioned on all data.
    double lpl = 0.0;
    dirtyAll();
    double logTruthProb = getLogValue();

    // DEBUG
    //cout << "logTruthProb = " << logTruthProb << endl;

    int totalPreds = 0;
    for (int q = 0; q < queryPreds.size(); q++) {

        int i = queryPreds[q]; 
        for (int j = 0; j < allFeatures_[i-1].size(); j++) {
            totalPreds++;
            if (allFeatures_[i-1][j] == NULL) {
                lpl += log(0.5);
            } else {
                double currValue = allFeatures_[i-1][j]->getValue();
                double newValue = (currValue == 0.0) ? 1.0 : 0.0;
                //setPredicateValue(i, j, newValue);
                setPredicateAndUpdate(i, j, newValue);
                double logUntruthProb = getLogValue();
#if 0
                // DEBUG
                if (isnan(logUntruthProb)) {
                    cout << "\nuntruthProb = " << logUntruthProb;
                } else {
                    cout << ".";
                }
#endif
                //setPredicateValue(i, j, currValue);
                setPredicateAndUpdate(i, j, currValue);


                double curr_lpl;
                if (logUntruthProb - logTruthProb > 100) {
                    curr_lpl = logTruthProb - logUntruthProb;
                } else if (logTruthProb - logUntruthProb > 100) {
                    curr_lpl = 0.0;
                } else {
                    curr_lpl = -log(1.0 + exp(logUntruthProb - logTruthProb));
                }
                lpl += curr_lpl;
#if 0
                // DEBUG
                //((PredicateGroundFeature*)allFeatures_[i-1][j])
                //    ->getPredicate()->printWithStrVar(cout,domain);
                cout << " = " << curr_lpl << endl;
                cout << logUntruthProb << " ; " << logTruthProb << endl;
                lpl += logTruthProb - log(exp(logTruthProb) + exp(logUntruthProb));
#endif
            }
        }
    }

    // DEBUG
    //cout << endl;

    return lpl/totalPreds;
}
