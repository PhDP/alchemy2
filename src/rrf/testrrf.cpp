#include "mln.h"
#include "fol.h"
#include "rrf.h"
#include "random.h"
#include <math.h>
#include "arguments.h"
#include "infer.h"  // For extractPredNames()
#include "timer.h"
#include "gfeature.h"
#include "lbfgsr.h"

#define BURN_IN 100
//#define BURN_IN 10

// DEBUG
int llcount = 0;

//#define MAX_TRAIN_ITERS 10000
#define CMP_RRF 0

enum inferenceMethod { INF_ICM, INF_GIBBS, INF_EXACT, INF_PSEUDO, INF_MAX };

char* afilename = NULL;
char* nonEvidPredsStr = "all";
char* aInferenceMethodStr = "p";
int   aGibbsIters = 100;
int   nf = 20;
int   infMethod = INF_PSEUDO;
double aAlpha = 0.1;
bool aRec = false;
bool aGround = false;
int  aNumFeatures = 10;
int  aPredsPerFeature = 5;
double aSigmaSq = 100;
char* atestConstant = NULL;
char* amodelFilename = NULL;
char* aoutputFilename = NULL;
int aMaxTrainIters = 10000;
bool aVerbose = false;
int  aSeed = 1;
bool aTrainTopLevel = false;
bool aTrainBottomLevel = false;
bool aInfer = false;
bool aLPL = false;
double aSamplingFrac = 1.0;
bool aPseudoFast = false;


ARGS ARGS::Args[] =
{
    ARGS("i", ARGS::Req, afilename, "input .mln file"),
    ARGS("weights", ARGS::Opt, amodelFilename,
            "Filename containing initial weights\n"),
    ARGS("r", ARGS::Opt, aoutputFilename,
            "File for saving learned model\n"),
    ARGS("ne", ARGS::Opt, nonEvidPredsStr, 
       //"(for discriminative learning only) "
       "first-order non-evidence predicates (comma-separated with no space)"),
    ARGS("inf", ARGS::Opt, aInferenceMethodStr, "Inference method to use\n"
            "i = ICM; g = Gibbs sampling; p = pseudo-likelihood; e = Exact\n"),
    ARGS("frac", ARGS::Opt, aSamplingFrac, 
           "Fraction of ground predicates to consider in psuedo-likelihood\n"),
    ARGS("iters", ARGS::Opt, aGibbsIters, 
            "Number of iterations to use for Gibbs sampling\n"),
    ARGS("nf", ARGS::Opt, aNumFeatures, "Number of features to use\n"),
    ARGS("alpha", ARGS::Opt, aAlpha, "Learning rate\n"),
    ARGS("sigmasq", ARGS::Opt, aSigmaSq, "Large-weight penalty\n"),
    ARGS("rec", ARGS::Opt, aRec, "Compute counts recursively\n"),
    ARGS("ppf", ARGS::Opt, aPredsPerFeature, 
            "Number of ground predicates per feature\n"),
    ARGS("ground", ARGS::Opt, aGround, 
            "Compute counts using fully ground tree\n"),
    ARGS("test", ARGS::Opt, atestConstant,
            "Test constant, whose query preds are unknown\n"),
    ARGS("trainiters", ARGS::Opt, aMaxTrainIters,
            "Maximum number of training iterations\n"),
    ARGS("v", ARGS::Tog, aVerbose,
            "Verbose output\n"),
    ARGS("seed", ARGS::Opt, aSeed, "Random seed\n"),
    ARGS("top", ARGS::Tog, aTrainTopLevel,
            "Only train top-level weights (feature coefficients)\n"),
    ARGS("bottom", ARGS::Tog, aTrainBottomLevel,
            "Only train lower-level (within-feature) weights\n"),
    ARGS("infer", ARGS::Tog, aInfer, 
            "Run inference and produce per-constant probabilities\n"),
    ARGS("lpl", ARGS::Tog, aLPL, 
            "Produces pseudo-log-likelihood for query predicates\n"),
    ARGS("f", ARGS::Tog, aPseudoFast,
      "Use optimized (and broken) pseudo-likelihood gradient computation.\n"),
    ARGS()
};

// HACKISH
Array<TruthValue> truthValues;
void saveState(Database* db, Array<Predicate*>& queryPreds)
{
    truthValues.clear();
    for (int i = 0; i < queryPreds.size(); i++) {
        truthValues.append(db->getValue(queryPreds[i]));
    }
}

void restoreState(Database* db, Array<Predicate*>& queryPreds)
{
    for (int i = 0; i < queryPreds.size(); i++) {
        db->setValue(queryPreds[i], truthValues[i]);
    }
}


void runICMInference(RRF* rrf, Database* db, Array<Predicate*>& queryPreds)
{
    // Initialize randomly
    for (int i = 0; i < queryPreds.size(); i++) {
        TruthValue tv = (frand() < 0.5) ? TRUE : FALSE;
        db->setValue(queryPreds[i], tv);
    }
    rrf->invalidateAll();

    // Toggle predicates one by one, until convergence
    bool changed = true;
    while (changed) {
        changed = false;

        double currLikelihood = rrf->getLogValue(db);
        for (int i = 0; i < queryPreds.size(); i++) {

            // Try toggling the value.
            //cout << *queryPreds[i] << endl;
            TruthValue oldValue = db->getValue(queryPreds[i]);
            TruthValue newValue = (oldValue == TRUE) ? FALSE : TRUE;

            db->setValue(queryPreds[i], newValue);
            rrf->changedPredicate(queryPreds[i], db);
            double newLikelihood = rrf->getLogValue(db);
            if (newLikelihood > currLikelihood) {
                // If it yields higher likelihood, keep it.
                currLikelihood = newLikelihood;
                changed = true;
            } else {
                // Otherwise, revert the change.
                db->setValue(queryPreds[i], oldValue);
                rrf->changedPredicate(queryPreds[i], db);
            }
        }
    }
}


void runICMInference(GroundRRF* grrf, const Array<int>& queryPreds,
        Database* db, RRF* rrf)
{
    // Initialize randomly
    for (int q = 0; q < queryPreds.size(); q++) {
#if CMP_RRF
        Array<Predicate*> queryGroundings;
        Predicate::createAllGroundings(queryPreds[q], db->getDomain(), 
                queryGroundings);
#endif
        for (int j = 0; j < grrf->getNumGroundings(queryPreds[q]-1); j++) {
            bool tv = (frand() < 0.5);
            grrf->setPredicateValue(queryPreds[q], j, tv);
#if CMP_RRF
            db->setValue(queryGroundings[j], tv ? TRUE : FALSE);
#endif
        }
    }

    // Toggle predicates one by one, until convergence
    bool changed = true;
    while (changed) {
        changed = false;

        double currLikelihood = grrf->getLogValue();
#if CMP_RRF
        double rrfLikelihood = rrf->getLogValue(db);
        if (fabs(rrfLikelihood - currLikelihood) > 0.00001) {
            cout << "Likelihoods differ!\n";
        }
#endif
        for (int q = 0; q < queryPreds.size(); q++) {

            int i = queryPreds[q];
#if CMP_RRF 
            Array<Predicate*> queryGroundings;
            Predicate::createAllGroundings(queryPreds[q], db->getDomain(), 
                    queryGroundings);
#endif
            for (int j = 0; j < grrf->getNumGroundings(i-1); j++) {

                // Try toggling the value.
                bool oldValue = grrf->getPredicateValue(i, j);
                grrf->setPredicateAndUpdate(i, j, !oldValue);
                //grrf->setPredicateValue(i, j, !oldValue);

                double newLikelihood = grrf->getLogValue();
#if CMP_RRF 
                db->setValue(queryGroundings[j], oldValue ? FALSE : TRUE);
                double rrfLikelihood = rrf->getLogValue(db);
                if (fabs(rrfLikelihood - newLikelihood) > 0.00001) {
                    cout << "Toggled likelihoods differ!\n";
                }
#endif
                if (newLikelihood > currLikelihood) {
                    // If it yields higher likelihood, keep it.
                    currLikelihood = newLikelihood;
                    changed = true;
                } else {
                    // Otherwise, revert the change.
                    grrf->setPredicateAndUpdate(i, j, oldValue);
                    //grrf->setPredicateValue(i, j, oldValue);
#if CMP_RRF
                    db->setValue(queryGroundings[j], oldValue ? TRUE : FALSE);
#endif
                }
            }
        }
    }
}


void getExactCounts(RRF* rrf, Database* db, Array<Predicate*>& queryPreds,
        Array<double>& counts)
{
    // Start with no counts
    rrf->getCounts(counts, db);
    for (int i = 0; i < counts.size(); i++) {
        counts[i] = 0.0;
    }

    // The below was cut and pasted from RRF:getExactZ()
    ArraysAccessor<TruthValue> predValues;
    Array<TruthValue> truthValues;
    truthValues.append(TRUE);
    truthValues.append(FALSE);

    for (int i = 0; i < queryPreds.size(); i++) {
        predValues.appendArray(&truthValues);
    }

    double Z = 0.0;

    // Sum value for all worlds (keeping the evidence pred values the same)
    do {
        Array<TruthValue> truthValues;
        predValues.getNextCombination(truthValues);
        for (int i = 0; i < truthValues.size(); i++) {
            db->setValue(queryPreds[i], truthValues[i]);
        }

        double likelihood = rrf->getValue(db);
        Z += likelihood;

        // Collect counts
        Array<double> newCounts;
        rrf->getCounts(newCounts, db);
        for (int j = 0; j < newCounts.size(); j++) {
            counts[j] += likelihood * newCounts[j];
        }

        //cout << "truthValues.size() = " << truthValues.size() << endl;
        //cout << "Likelihood " << index << " = " << likelihood << endl;
    } while (predValues.hasNextCombination());

    // Normalize
    for (int j = 0; j < counts.size(); j++) {
        counts[j] /= Z;
    }
}

void getGibbsCounts(RRF* rrf, Database* db, Array<Predicate*>& queryPreds,
        Array<double>& counts)
{
    // Get initial likelihood and counts
    double oldLikelihood = rrf->getValue(db);
    rrf->getCounts(counts, db);

    // TODO -- optimize...
    for (int iter = 0; iter < aGibbsIters; iter++) {
        for (int i = 0; i < queryPreds.size(); i++) {

            // Sample i'th predicate, conditioned on the rest
            TruthValue oldValue = db->getValue(queryPreds[i]);
            TruthValue newValue = (oldValue == TRUE) ? FALSE : TRUE;

            db->setValue(queryPreds[i], newValue);
            rrf->changedPredicate(queryPreds[i], db);
            double newLikelihood = rrf->getValue(db);
            double prob = newLikelihood/(oldLikelihood + newLikelihood);
            if (frand() < prob) {
                oldLikelihood = newLikelihood;
            } else {
                db->setValue(queryPreds[i], oldValue);
                rrf->changedPredicate(queryPreds[i], db);
            }
        }

        // Collect counts
        Array<double> newCounts;
        rrf->getCounts(newCounts, db);
        for (int j = 0; j < newCounts.size(); j++) {
            counts[j] += newCounts[j];
        }
    }

    // Normalize
    for (int j = 0; j < counts.size(); j++) {
        counts[j] /= ((double)aGibbsIters+1.0);
    }
}

void runGibbs(RRF* rrf, Database* db, Array<Predicate*>& queryPreds,
        Array<double>& counts)
{
    // Get initial likelihood and counts
    double oldLikelihood = rrf->getValue(db);
    counts.clear();

    // TODO -- optimize...
    for (int iter = 0; iter < aGibbsIters; iter++) {
        for (int i = 0; i < queryPreds.size(); i++) {

            // Sample i'th predicate, conditioned on the rest
            TruthValue oldValue = db->getValue(queryPreds[i]);
            TruthValue newValue = (oldValue == TRUE) ? FALSE : TRUE;
            TruthValue finalValue;

            db->setValue(queryPreds[i], newValue);
            rrf->changedPredicate(queryPreds[i], db);
            double newLikelihood = rrf->getValue(db);
            double prob = newLikelihood/(oldLikelihood + newLikelihood);
            if (frand() < prob) {
                oldLikelihood = newLikelihood;
                finalValue = newValue;
            } else {
                db->setValue(queryPreds[i], oldValue);
                rrf->changedPredicate(queryPreds[i], db);
                finalValue = oldValue;
            }

            // Keep track of counts
            if (counts.size() <= i) {
                counts.append(0.5);
            }
            if (finalValue == TRUE) {
                counts[i]++;
            }
        }
    }

    // Normalize
    for (int j = 0; j < counts.size(); j++) {
        counts[j] /= ((double)aGibbsIters+1.0);
    }
}

double getGibbsLogLikelihood(RRF* rrf, Database* db, 
        Array<Predicate*>& queryPreds)
{
    Array<double> counts;

    Array<TruthValue> trueValues;

    for (int i = 0; i < queryPreds.size(); i++) {
        trueValues.append(db->getValue(queryPreds[i]));
        //db->setValue(queryPreds[i], (frand() > 0.5) ? TRUE : FALSE);
        counts.append(0.1);
    }
    double total = 0.2;

    // Initialize using ICM
    runICMInference(rrf, db, queryPreds);


    rrf->invalidateAll();
    double oldLikelihood = rrf->getLogValue(db);

    for (int iter = 0; iter < aGibbsIters + BURN_IN; iter++) {
        for (int i = 0; i < queryPreds.size(); i++) {

            // Sample i'th predicate, conditioned on the rest
            TruthValue oldValue = db->getValue(queryPreds[i]);
            TruthValue newValue = (oldValue == TRUE) ? FALSE : TRUE;
            TruthValue finalValue;
            
            db->setValue(queryPreds[i], newValue);
            rrf->changedPredicate(queryPreds[i], db);
            double newLikelihood = rrf->getLogValue(db);
            double prob = 1.0/(1.0 + exp(oldLikelihood - newLikelihood));
            if (oldLikelihood - newLikelihood > 100) {
                prob = 0.0;
            } else if (oldLikelihood - newLikelihood < -100) {
                prob = 1.0;
            }
            if (frand() < prob) {
                oldLikelihood = newLikelihood;
                finalValue = newValue;
            } else {
                db->setValue(queryPreds[i], oldValue);
                rrf->changedPredicate(queryPreds[i], db);
                finalValue = oldValue;
            }

            if (iter >= BURN_IN && finalValue) {
                counts[i]++;
            }
        }

        if (iter >= BURN_IN) {
            total++;
        }
    }

    double ll = 0.0;

    // Normalize
    for (int j = 0; j < counts.size(); j++) {

        if (trueValues[j]) {
            ll += log(counts[j]/total);
        } else {
            ll += log(1.0 - counts[j]/total);
        }

        cout << *(queryPreds[j]) << " = " << trueValues[j] 
            << " (" << (counts[j]/total) << ")\n";
    }

    //cout << right/(right + wrong) << endl;
    //cout << naive/(right + wrong) << endl;

    return ll;
}


double getGibbsLogLikelihood2(RRF* rrf, Database* db, 
        Array<Predicate*>& queryPreds)
{
    Array<int> queryPredId;
    Array<int> queryGroundId;

    Array<double> counts;
    Array<TruthValue> trueValues;

    double ll = 0.0;

    for (int i = 0; i < queryPreds.size(); i++) {
        trueValues.append(db->getValue(queryPreds[i]));

        // Store identity of query predicates
        queryPredId.append(queryPreds[i]->getId());
        Array<int> grounding;
        for (int j = 0; j < queryPreds[i]->getNumTerms(); j++) {
            grounding.append(queryPreds[i]->getTerm(j)->getId());
        }
        queryGroundId.append(rrf->getFeature(queryPredId[i]-1)->
                getGroundingIndex(grounding, db));

        counts.append(0.1);
    }
    double total = 0.2;


    for (int trial= 0; trial < 10; trial++) {

        Array<int> qpreds;
        qpreds.append(queryPreds[0]->getId());
        GroundRRF* grrf = new GroundRRF(rrf, db);
        grrf->dirtyAll();

        double best_ll = -1.0e100;
        for (int init = 0; init < 10; init++) {
            runICMInference(grrf, qpreds, db, rrf);
            if (grrf->getLogValue() > best_ll) {
                best_ll = grrf->getLogValue();
                saveState(db, queryPreds);
            }
        }
        restoreState(db, queryPreds);

        double oldLikelihood = grrf->getLogValue();
        for (int iter = 0; iter < aGibbsIters + BURN_IN; iter++) {
            for (int i = 0; i < queryPreds.size(); i++) {

                // Sample i'th predicate, conditioned on the rest
                bool oldValue = grrf->getPredicateValue(
                        queryPredId[i], queryGroundId[i]);
                bool finalValue;

                grrf->setPredicateAndUpdate(queryPredId[i], queryGroundId[i], 
                        !oldValue);
                double newLikelihood = grrf->getLogValue();
                //double prob = newLikelihood/(oldLikelihood + newLikelihood);
                double prob = 1.0/(1.0 + exp(oldLikelihood - newLikelihood));
                if (oldLikelihood - newLikelihood > 100) {
                    prob = 0.0;
                } else if (oldLikelihood - newLikelihood < -100) {
                    prob = 1.0;
                }
                //cout << prob << endl;
                if (frand() < prob) {
                    oldLikelihood = newLikelihood;
                    finalValue = !oldValue;
                } else {
                    grrf->setPredicateAndUpdate(queryPredId[i], 
                            queryGroundId[i], oldValue);
                    finalValue = oldValue;
                }

                if (iter >= BURN_IN && finalValue) {
                    counts[i]++;
                }
            }

            if (iter >= BURN_IN) {
                total++;
            }
        }
    }

    // Normalize
    for (int j = 0; j < counts.size(); j++) {

        double curr_ll;
        if (trueValues[j]) {
            curr_ll = log(counts[j]/total);
        } else {
            curr_ll = log(1.0 - counts[j]/total);
        }
        ll += curr_ll;

        queryPreds[j]->printWithStrVar(cout, db->getDomain());
        cout << " = " << trueValues[j] 
            << " (" << (counts[j]/total) << ")  " << curr_ll << endl;

    }

    return ll;
}


#define USE_LOGS 1

void getGibbsCounts(GroundRRF* grrf, const Array<int>& queryPreds,
        Array<double>& counts, Database* db, RRF* rrf)
{
    // Get initial likelihood and counts
#if USE_LOGS
    double oldLikelihood = grrf->getLogValue();
#else
    double oldLikelihood = grrf->getValue();
#endif
    grrf->getCounts(counts);

#if CMP_RRF
    Array<double> rrfCounts;

    for (int q = 0; q < queryPreds.size(); q++) {
        int i = queryPreds[q];
        Array<Predicate*> queryGroundings;
        Predicate::createAllGroundings(i, db->getDomain(), queryGroundings);
        for (int j = 0; j < queryGroundings.size(); j++) {
            bool value = grrf->getPredicateValue(i,j);
            db->setValue(queryGroundings[j], value ? TRUE : FALSE);
        }
    }
#endif

    Array<double> newCounts;
    for (int iter = 0; iter < aGibbsIters; iter++) {
        for (int q = 0; q < queryPreds.size(); q++) {
            int i = queryPreds[q];
            // NOTE: assuming that predicate 0 is equality and therefore 
            // skipped.  This is the source of the (i-1).
            for (int j = 0; j < grrf->getNumGroundings(i-1); j++) {

                // Sample predicate, conditioned on the rest
                bool oldValue = grrf->getPredicateValue(i,j);
                grrf->setPredicateAndUpdate(i,j,!oldValue);
                //grrf->setPredicateValue(i,j,!oldValue);
#if CMP_RRF
                db->setValue(queryGroundings[j], oldValue ? FALSE : TRUE);
#endif

#if USE_LOGS
                double newLikelihood = grrf->getLogValue();
#else
                double newLikelihood = grrf->getValue();
#endif

#if CMP_RRF  // TODO -- fix this...
                double rrfLikelihood = rrf->getValue(db);
                if (fabs(rrfLikelihood - newLikelihood) > 0.00001) {
                    cout << "Different likelihoods: \n";
                    cout << rrfLikelihood << endl;
                    cout << newLikelihood << endl;
                }
#endif

#if USE_LOGS
                double prob = 1.0/(1.0 + exp(oldLikelihood - newLikelihood));
                if (oldLikelihood - newLikelihood > 100) {
                    prob = 0.0;
                } else if (oldLikelihood - newLikelihood < -100) {
                    prob = 1.0;
                }
#else
                double prob = newLikelihood/(oldLikelihood + newLikelihood);
#endif

                if (frand() < prob) {
                    oldLikelihood = newLikelihood;
                } else {
                    grrf->setPredicateAndUpdate(i,j,oldValue);
                    //grrf->setPredicateValue(i,j,oldValue);
#if CMP_RRF
                    db->setValue(queryGroundings[j], oldValue ? TRUE : FALSE);
#endif
                }

#if CMP_RRF
                db->setValue(queryGroundings[j], 
                        grrf->getPredicateValue(i,j) ? TRUE : FALSE);
#endif
            }
        }

        if (grrf == NULL) {
            cout << "grrf is somehow NULL!\n";
        }

        // Collect counts
        grrf->dirtyAll();
        grrf->getCounts(newCounts);
#if CMP_RRF
        rrf->getCounts(rrfCounts, db);
#endif
        for (int c = 0; c < newCounts.size(); c++) {
#if CMP_RRF
            if (fabs(newCounts[c] - rrfCounts[c]) > 0.000001) {
                cout << "Different counts for " << c << ":\n";
                cout << "Ground: " << newCounts[c] << endl;
                cout << "True:   " << rrfCounts[c] << endl;
            } 
#endif
            counts[c] += newCounts[c];
        }
    }

    // Normalize
    for (int c = 0; c < counts.size(); c++) {
        counts[c] /= ((double)aGibbsIters+1.0);
    }
}


void runGibbs(GroundRRF* grrf, const Array<int>& queryPreds,
        Array<double>& counts, Database* db)
{
    // Get initial likelihood and counts
    double oldLikelihood = grrf->getLogValue();

    counts.clear();
    for (int iter = 0; iter < aGibbsIters; iter++) {
        int countIndex = 0;
        for (int q = 0; q < queryPreds.size(); q++) {
            int i = queryPreds[q];
            // NOTE: assuming that predicate 0 is equality and therefore 
            // skipped.  This is the source of the (i-1).
            for (int j = 0; j < grrf->getNumGroundings(i-1); j++) {

                // Sample predicate, conditioned on the rest
                bool oldValue = grrf->getPredicateValue(i,j);
                grrf->setPredicateAndUpdate(i,j,!oldValue);
                //grrf->setPredicateValue(i,j,!oldValue);
                bool finalValue;

                double newLikelihood = grrf->getLogValue();
                double prob = 1.0/(1.0 + exp(oldLikelihood - newLikelihood));
                if (oldLikelihood - newLikelihood > 100) {
                    prob = 0.0;
                } else if (oldLikelihood - newLikelihood < -100) {
                    prob = 1.0;
                }

                if (frand() < prob) {
                    oldLikelihood = newLikelihood;
                    finalValue = !oldValue;
                } else {
                    grrf->setPredicateAndUpdate(i,j,oldValue);
                    finalValue = oldValue;
                }

                if (counts.size() <= countIndex) {
                    counts.append(0.5);
                }
                if (finalValue) {
                    counts[countIndex]++;
                }
                countIndex++;
            }
        }

        if (grrf == NULL) {
            cout << "grrf is somehow NULL!\n";
        }
    }

    // Normalize
    for (int c = 0; c < counts.size(); c++) {
        counts[c] /= ((double)aGibbsIters+1.0);
    }
}

void scrambleWeights(RRF* rrf) 
{
    // Set initial (random) weights
    Array<Feature*> features = rrf->getFeatureArray();
    for (int i = 0; i < features.size(); i++) {
        for (int j = 0; j < features[i]->getNumWeights(); j++) {
            // Initialize to a random weight in (-0.5, +0.5)
            double randWt = frand() - 0.5;
            features[i]->setWeight(j, randWt);
        }
    }
}

void trainRRF(RRF* rrf, Database* db, const Array<int>& queryPreds)
{
    // Set initial (random) weights
    // scrambleWeights(rrf);

    // Initial preparations for recursive counts computation
    Array<Predicate*> allGroundings;
    Array<TruthValue> truePredValues;

    if (aRec || CMP_RRF) {
        const Domain* domain = db->getDomain();
        for (int i = 0; i < queryPreds.size(); i++) {

            // Create all groundings for this predicate
            Array<Predicate*> predGroundings;
            Predicate::createAllGroundings(queryPreds[i], domain, predGroundings);
            allGroundings.append(predGroundings);

            // Save all original truth values (assuming complete data for now)
            for (int j = 0; j < predGroundings.size(); j++)  {
                truePredValues.append(db->getValue(predGroundings[j]));
            }
        }
    }

    // Initial preparations for full ground tree counts computation
    GroundRRF* grrf = NULL;
#if 1
    GroundRRF* trueGrrf = NULL;
#else
    Array<Array<bool> > trueValues;
#endif

    if (aGround) {
        grrf = new GroundRRF(rrf, db);
#if 1
        trueGrrf = new GroundRRF(rrf, db);
#else
        // Save true values of all query predicates
        for (int q = 0; q < queryPreds.size(); q++) {
            int i = queryPreds[q];
            trueValues.append(Array<bool>());
            for (int j = 0; j < grrf->getNumGroundings(i-1); j++) {
                trueValues[q].append(grrf->getPredicateValue(i,j));
            }
        }
#endif
    }


    if (infMethod == INF_PSEUDO) {

        // Print out initial lpl
        cout << "ground lpl = " 
            << trueGrrf->getLogPseudoLikelihood(queryPreds) << endl;

        // Get weights
        int numWeights = rrf->getNumWeights();
        double wts[numWeights+1];
        for (int i = 0; i < numWeights; i++) {
            wts[i+1] = rrf->getWeight(i);
        }


        int minWt = 0;
        int maxWt = numWeights;

        if (aTrainBottomLevel) {
            minWt = rrf->getRoot()->getNumWeights();
        }
        if (aTrainTopLevel) {
            maxWt = rrf->getRoot()->getNumWeights();
        }

        // Solve
        int iter;
        bool error;
        LBFGSR solver(aMaxTrainIters, 1e-5, aSamplingFrac, rrf, trueGrrf, 
                queryPreds, maxWt - minWt, minWt, aSigmaSq, aPseudoFast);
        solver.minimize(wts+minWt, iter, error);
        cout << "Iters required: " << iter << endl;

        // Save weights
        for (int i = 0; i < numWeights; i++) {
            rrf->setWeight(i, wts[i+1]);
        }

        trueGrrf->dirtyAll();

        // Print out final lpl
        cout << "ground lpl = " 
            << trueGrrf->getLogPseudoLikelihood(queryPreds) << endl;
    } else {


    // Loop until convergence
    Timer timer;
    double lastWriteSec = timer.time();
    for (int iter = 0; iter < aMaxTrainIters; iter++) {

        double begSec = timer.time();

#if 0
        for (int i = 0; i < rrf->getNumFeatures(); i++) {
            cout << "F" << i << ":";
            Feature* feature = rrf->getFeature(i);
            for (int j = 0; j < feature->getNumWeights(); j++) {
                cout << " " << feature->getWeight(j);
            }
            cout << "\n";
        }
#endif

        // Print out log likelihood
        if (aRec) {
            cout << "rec lpl = " << rrf->getLogPseudoLikelihood(db, queryPreds) 
                << endl;
        }
        if (aGround) {
#if 0
            cout << "ground lpl = " << grrf->getLogPseudoLikelihood(queryPreds) 
                << endl;
#else
            cout << "ground lpl = " 
                << trueGrrf->getLogPseudoLikelihood(queryPreds) << endl;
            cout << "ground ll = "
                << trueGrrf->getLogValue() << endl;
#if 0
            // DEBUG
            if (llcount++ % 100 == 0) {
                cout << "ll = " << log(rrf->getExactConditionalLikelihood(db, queryPreds)) << endl;
            }
#endif
#endif
#if CMP_RRF
            cout << "rec lpl = " << rrf->getLogPseudoLikelihood(db, queryPreds) 
                << endl;
#endif
        }
        cout << "wll = " << rrf->getWeightLogLikelihood(aSigmaSq) << endl;

        // Run inference
        if (aRec) {
            if (infMethod != INF_EXACT) {
                runICMInference(rrf, db, allGroundings);
            }
        }
        if (aGround) {
            if (infMethod == INF_ICM) {
                runICMInference(grrf, queryPreds, db, rrf);
            }
        }
        
        // Get inferred counts
        Array<double> groundInfCounts;
        Array<double> recInfCounts;

        if (infMethod == INF_GIBBS) {
            if (aRec)    { getGibbsCounts(rrf, db, allGroundings, recInfCounts); }
            if (aGround) { getGibbsCounts(grrf, queryPreds, groundInfCounts, 
                    db, rrf); }
        } else if (infMethod == INF_ICM) {
            if (aRec)    { rrf->getCounts(recInfCounts, db); }
            if (aGround) { grrf->getCounts(groundInfCounts); }
        } else if (infMethod == INF_EXACT) {
            if (aRec)    { getExactCounts(rrf, db, allGroundings, recInfCounts); }
            // TODO -- implement this
            if (aGround) { cout << "ERROR: exact inference not supported for ground network.\n";  }
        } else if (infMethod == INF_PSEUDO) {
            // No need to do any expectation!
        } else {
            cout << "ERROR: unknown inference method " << infMethod << endl;
        }

        // Reset database
#if 0
        if (aGround) {
            for (int q = 0; q < queryPreds.size(); q++) {
                for (int j = 0; j < trueValues[q].size(); j++) {
                    grrf->setPredicateValue(queryPreds[q], j, trueValues[q][j]);
                }
            }
        }
#endif

        if (aRec || CMP_RRF) 
        {
            // Reset database
            for (int i = 0; i < allGroundings.size(); i++) {
                db->setValue(allGroundings[i], truePredValues[i]);
            }
            rrf->invalidateAll();
        }

        // Get true counts
        Array<double> groundTrueCounts;
        Array<double> recTrueCounts;

        if (aGround) {
#if 0
            grrf->getCounts(groundTrueCounts);
#else
            if (infMethod == INF_PSEUDO) {
                if (aPseudoFast) {
                    trueGrrf->getPseudoCountsFast(groundTrueCounts, 
                            queryPreds, aSamplingFrac);
                } else {
                    trueGrrf->getPseudoCounts(groundTrueCounts, queryPreds,
                            aSamplingFrac);
                }
            } else {
                trueGrrf->getCounts(groundTrueCounts);
            }
#endif
        }
        if (aRec) {
            if (infMethod == INF_PSEUDO) {
                // TODO... implement this properly
                cout << "ERROR: pseudo-likelihood requires ground RRF.\n";
            } else {
                rrf->getCounts(recTrueCounts, db);
            }
        }

        // Compute and follow gradient for one step
        int iMin = 0;
        int iMax = aRec ? recTrueCounts.size() : groundTrueCounts.size();

        if (aTrainBottomLevel) {
            iMin = rrf->getRoot()->getNumWeights();
        }
        if (aTrainTopLevel) {
            iMax = rrf->getRoot()->getNumWeights();
        }


        Array<double> gradient;
        for (int i = 0; i < iMin; i++) {
            gradient.append(0.0);
        }

        // DEBUG
        double sqlength = 0.0;
        double weightsum = 0.0;
        int maxIndex = 0;
        double maxTotal = 0.0;
        for (int i = iMin; i < iMax; i++) 
        {
#if CMP_RRF
            // NOTE: no support for psuedo-log-likelihood
            if (aRec && aGround) {
                if (2.0*fabs((recTrueCounts[i] - groundTrueCounts[i])/
                        (recTrueCounts[i] + groundTrueCounts[i])) > 0.01) {
                    cout << "True counts of " << i << " differ.\n";
                }
                if (2.0*fabs(recInfCounts[i] - groundInfCounts[i]) > 0.0001) {
                    cout << "Inferred counts of " << i << " differ:\n";
                    cout << "Rec:    " << recInfCounts[i] << endl;
                    cout << "Ground: " << groundInfCounts[i] << endl;
                }

                double recDiff = recInfCounts[i] - recTrueCounts[i];
                double groundDiff = groundInfCounts[i] - groundTrueCounts[i];
                if (2.0*fabs(recDiff - groundDiff) > 0.0001) {
                    cout << "Gradient differs on " << i << endl;
                    cout << "Rec:    " << recDiff << endl;
                    cout << "Ground: " << groundDiff << endl;
                }
            }
#endif
            double diff;
            if (infMethod == INF_PSEUDO) {
                diff = groundTrueCounts[i];
            } else if (aRec) {
                diff = recTrueCounts[i] - recInfCounts[i];
            } else {
                diff = groundTrueCounts[i] - groundInfCounts[i];
            }
            double weight = rrf->getWeight(i);
            double total = diff - weight/aSigmaSq;
            gradient.append(total);
            sqlength += total * total;

            // DEBUG
            if (fabs(total) > fabs(maxTotal)) {
                maxTotal = total;
                maxIndex = i;
            }
            weightsum += fabs(weight);
        }
        //double norm = sqrt(sqlength);

#if 0
        // DEBUG
        cout << "largest change = " << maxTotal*aAlpha << " (" << maxIndex << ")\n";
        cout << "sqlength = " << sqlength*aAlpha*aAlpha << endl;
        cout << "length = " << sqrt(sqlength*aAlpha*aAlpha) << endl;
        cout << "sum of weights = " << weightsum << endl;
#endif

        for (int i = iMin; i < iMax; i++) 
        {
            // Take current weight, add normalized(?) gradient
            double weight = rrf->getWeight(i);
            rrf->setWeight(i, weight + gradient[i]*aAlpha
                   // *weightsum
                   // /norm
                    );
        }

        // Weights changed; cache invalid!
        if (aGround) {
            grrf->dirtyAll();
#if 1
            trueGrrf->dirtyAll();
#endif
        }

        rrf->invalidateAll();

        // Print out time for this iteration
        cout << "Time: ";
        Timer::printTime(cout, timer.time() - begSec); cout << endl;

        if (aVerbose) {
            cout << rrf;
        }

        if (timer.time() - lastWriteSec > 5) {
            ofstream fout(aoutputFilename);
            fout << rrf;
            fout.close();
            lastWriteSec = timer.time();
        }
    }
    }

    // One final save.
    ofstream fout(aoutputFilename);
    fout << rrf;

    if (aVerbose) {
        cout << "Final model:\n";
        cout << rrf;
    }
    fout.close();
}


void inferRRF(RRF* rrf, Database* db, const Array<int>& queryPreds)
{
    // Initial preparations for recursive counts computation
    Array<Predicate*> allGroundings;
    Array<TruthValue> truePredValues;
    const Domain* domain = db->getDomain();

    //if (aRec || CMP_RRF) 
    if (1) {
        for (int i = 0; i < queryPreds.size(); i++) {

            // Create all groundings for this predicate
            Array<Predicate*> predGroundings;
            Predicate::createAllGroundings(queryPreds[i], domain, predGroundings);
            allGroundings.append(predGroundings);

            // Save all original truth values (assuming complete data for now)
            for (int j = 0; j < predGroundings.size(); j++)  {
                truePredValues.append(db->getValue(predGroundings[j]));
            }
        }
    }

    // Initial preparations for full ground tree counts computation
    GroundRRF* grrf = NULL;
#if 1
    GroundRRF* trueGrrf = NULL;
#else
    Array<Array<bool> > trueValues;
#endif

    if (aGround) {
        grrf = new GroundRRF(rrf, db);
#if 1
        trueGrrf = new GroundRRF(rrf, db);
#else
        // Save true values of all query predicates
        for (int q = 0; q < queryPreds.size(); q++) {
            int i = queryPreds[q];
            trueValues.append(Array<bool>());
            for (int j = 0; j < grrf->getNumGroundings(i-1); j++) {
                trueValues[q].append(grrf->getPredicateValue(i,j));
            }
        }
#endif
    }

    // Run inference
    if (aRec) {
        runICMInference(rrf, db, allGroundings);
    }

    if (aGround) {
        runICMInference(grrf, queryPreds, db, rrf);
    }

    Array<double> recInfCounts;
    Array<double> groundInfCounts;
    
    getGibbsLogLikelihood2(rrf, db, allGroundings);

#if 0
    if (infMethod == INF_GIBBS) {
        if (aRec)    { runGibbs(rrf, db, allGroundings, recInfCounts); }
        if (aGround) { runGibbs(grrf, queryPreds, groundInfCounts, db); }
    } else if (infMethod == INF_ICM) {
        //if (aRec)    { rrf->getCounts(recInfCounts, db); }
        //if (aGround) { grrf->getCounts(groundInfCounts); }
        // TODO!
    } 

    // Print out all probabilities
    for (int i = 0; i < allGroundings.size(); i++) {
        allGroundings[i]->printWithStrVar(cout, domain);
        if (aGround) {
            cout << " " << groundInfCounts[i];
        }
        if (aRec) {
            cout << " " << recInfCounts[i];
        }
        cout << endl;
    }
#endif
}



int main(int argc, char* argv[]) 
{
    ARGS::parse(argc, argv, &cout);

    switch(aInferenceMethodStr[0]) {
        case 'i': infMethod = INF_ICM;
                  break;
        case 'g': infMethod = INF_GIBBS;
                  break;
        case 'e': infMethod = INF_EXACT;
                  break;
        case 'p': infMethod = INF_PSEUDO;
                  break;
        default:  cout << "ERROR: Unknown inference method "
                  << aInferenceMethodStr << endl;
                  exit(-1);
    }

    if (!aRec && !aGround) {
        aGround = true;
    }

    srand(aSeed);

    //// Load domain ////
    MLN* mln = new MLN; // required, but unused...
    Domain* domain = new Domain;
    StringHashArray* openWorldPredNames = new StringHashArray();
    //StringHashArray* closedWorldPredNames = new StringHashArray();
    StringHashArray* queryPredNames = new StringHashArray();
    runYYParser(mln, domain, afilename, true, 
            openWorldPredNames, queryPredNames, false, 
            false, 0, false, 
            NULL, true, false);
    //delete mln; mln = NULL;

    // Get a list of query predicates
    Array<int> queryPreds;
    if (!strcmp(nonEvidPredsStr, "all")) {

        // Skip the equality predicate
        for (int i = 1; i < domain->getNumPredicates(); i++) {
            queryPreds.append(i);
        }
    } else {
        StringHashArray nonEvidPredNames;
        if(!extractPredNames(nonEvidPredsStr, NULL, nonEvidPredNames)) {
            cout << "ERROR: failed to extract non-evidence predicate names.\n";
            return -1;
        }

        for (int i = 0; i < nonEvidPredNames.size(); i++) {
            int predId = domain->getPredicateId( nonEvidPredNames[i].c_str() );
            if (predId < 0) {
                cout << "Ignoring unrecognized query predicate \""
                    << nonEvidPredNames[i] << "\".\n";
            } else {
                queryPreds.append(predId);
            }
        }
    }


    RRF* rrf;

    if (amodelFilename) {
        // DEBUG
        cout << "Loading weights from file.\n";
        rrf = new RRF();
        ifstream modelIn(amodelFilename);
        rrf->load(modelIn, domain);
        modelIn.close();

    } else {
#if 0
        //// Create RRF ////
        rrf = RRF::generateTwoLevel(domain);
        Array<int> typeArity;
        typeArity.append(0);
        typeArity.append(1);
        for (int i = 0; i < aNumFeatures; i++) {
            rrf->addCompleteFeature(typeArity, domain, queryPreds, 
                    aPredsPerFeature);
        }
#else
        cout << "ERROR: no longer supported!\n";
        exit(-1);
#endif
    }

#if 0
    if (atestConstant) {
        // HACK -- do leave-one-out experiments
        int testId = domain->getConstantId(atestConstant);
        if (testId < 0) {
            cout << "Error: unknown constant \"" << atestConstant << "\"\n";
        }

        // Make list of all query preds involving the held-out constant
        Array<Predicate*> groundQueryPreds;
        Array<int> req;
        req.append(testId);

        // For each query predicate...
        for (int i = 0; i < queryPreds.size(); i++) {

            // Iterate through all groundings that include the held-out
            // constant.

            cout << "req = ";
            for (int j = 0; j < req.size(); j++) {
                cout << req[j] << " ";
            }
            cout << endl;
            // Broken...
           // ParentIter2 iter(req, domain->getConstantsByType(), 
            //    domain->getPredicateTemplate(queryPreds[i])->getNumTerms());
            Array<int> grounding;
            while (iter.hasNextGrounding()) {

                // Copy the groundings to an instance of a predicate
                iter.getNextGrounding(grounding);
                Predicate* pred = new Predicate(
                        domain->getPredicateTemplate(queryPreds[i]));
                for (int j = 0; j < grounding.size(); j++) {
                    pred->setTermToConstant(j, grounding[j]);
                }

                // Add the predicate to the list of all test preds
                groundQueryPreds.append(pred);
            }
        }

        double lpl = rrf->getLogPseudoLikelihood(domain->getDB(),
                groundQueryPreds);

        // Compute log likelihood of test preds
        double ll;
        if (aGround) {
            ll = getGibbsLogLikelihood2(rrf, domain->getDB(), 
                groundQueryPreds);
        } else {
            ll = getGibbsLogLikelihood(rrf, domain->getDB(), 
                groundQueryPreds);
        }
        cout << "lpl = " << lpl << endl;
        cout << "ll = " << ll << endl;

    } else 
#endif
    if (aInfer) {
        inferRRF(rrf, domain->getDB(), queryPreds);
    } else if (aLPL) {
        if (aGround) {
            GroundRRF* grrf = new GroundRRF(rrf, domain->getDB());
            cout << "lpl = " << grrf->getLogPseudoLikelihood(queryPreds) 
                << endl;
        } else {
            cout << "ERROR: pll for non-ground RRF is not yet supported!\n";
        }
    } else {
        if (aoutputFilename == NULL) {
            cout << "Warning: you should specify an output model filename "
                << "using the -r option.\n";
        }
        trainRRF(rrf, domain->getDB(), queryPreds);
    }

    return 0;
}
