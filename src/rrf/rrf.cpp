#include "rrf.h"
#include "domain.h"

RRF* RRF::generateTwoLevel(Domain* domain) 
{
    RRF* rrf = new RRF();

    // Create predicate (bottom level) features
    // HACK -- ignore equality predicate, 0
    for (int pred = 1; pred < domain->getNumPredicates(); pred++) {
        PredicateFeature* predFeat = new PredicateFeature(
                domain->getPredicateTemplate(pred));
        rrf->maxFeatureId_++;
        predFeat->setId(rrf->maxFeatureId_);
        rrf->featureArray_.append(predFeat);
    }

    // Create top-level feature
    rrf->topFeature_ = new RecursiveFeature("root", true, false);
    rrf->maxFeatureId_++;
    rrf->topFeature_->setId(rrf->maxFeatureId_);
    rrf->featureArray_.append(rrf->topFeature_);

    return rrf;
}

void RRF::addCompleteFeature(Array<int> typeArity, Domain* domain, 
        const Array<int>& queryPreds, int numChildren)
{
    char clausalName[100];
    sprintf(clausalName, "f%d", ++maxFeatureId_);
    ClausalFeature* feature = new ClausalFeature(clausalName);
    //RecursiveFeature* feature = new RecursiveFeature(clausalName);
    feature->setId(maxFeatureId_);
    Array<Array<int> > termIndicesByType;

    // For each type...
    for (int i = 0; i < typeArity.size(); i++) {

        // Keep a list of all indices for terms of the current type
        Array<int> currTypeTermIndices;
        for (int j = 0; j < typeArity[i]; j++) {

            // Add a term of the appropriate type to the current
            // feature, and record its index in the array.
            feature->addTermType(i);
            currTypeTermIndices.append(feature->getNumTerms() - 1);
        }

        // Append this array to the list that includes all types
        termIndicesByType.append(currTypeTermIndices);
    }

    // Add all groundings of each predicate
    // HACK: start at 1 to avoid equality predicate
    for (int pred = 1; pred < domain->getNumPredicates(); pred++) {

        const PredicateTemplate* predTemp = 
                domain->getPredicateTemplate(pred);

        // Construct iterator over all term mappings of 
        // this predicate.
        // DEBUG
        //cout << predTemp->getNumTerms() << endl;
        ArraysAccessor<int> mapIter;
        for (int i = 0; i < predTemp->getNumTerms(); i++) {
            int termType = predTemp->getTermTypeAsInt(i);
            mapIter.appendArray(&(termIndicesByType[termType]));
            //Array<int>* copiedArray 
            //    = new Array<int>(termIndicesByType[termType]);
            //mapIter.appendArray(copiedArray);
        }

        // Add the predicate feature with each possible mapping.
        // Assumption: feature 0 is the top level feature;
        //   next features are the base predicate features.
        Feature* predFeature = featureArray_[pred-1];
        Array<int> mapping;
        do {
            mapIter.getNextCombination(mapping);
            feature->addChild(predFeature, 0, mapping);
        } while (mapIter.hasNextCombination());
    }

    Array<int> mapping;
    for (int i = 0; i < feature->getNumTerms(); i++) {
        mapping.append(-i - 1);
    }
    ((RecursiveFeature*)topFeature_)->addChild(feature, 0.0, mapping);


    // Use random example to set weights
    Array<int> grounding;
    for (int i = 0; i < feature->getNumTerms(); i++) {
        int type = feature->getTermType(i);
        
        const Array<int>* constants = domain->getConstantsByType(type);
        grounding.append((*constants)[rand() % constants->size()]);
    }
    feature->setWeightsFromExample(grounding, domain->getDB(), queryPreds, numChildren);
    topFeature_->setWeight(topFeature_->getNumWeights() - 1, frand());

    featureArray_.append(feature);
}

void RRF::addRandomFeature(int numChildren, Array<int> queryPreds,
        Array<int> typeArity, Domain* domain)
{
    ClausalFeature* feature = NULL;

    Array<Array<int> > termIndicesByType;

    bool queryPredChosen = false;
    while (!queryPredChosen) {

        feature = new ClausalFeature;

        // For each type...
        for (int i = 0; i < typeArity.size(); i++) {

            // Keep a list of all indices for terms of the current type
            Array<int> currTypeTermIndices;
            for (int j = 0; j < typeArity[i]; j++) {

                // Add a term of the appropriate type to the current
                // feature, and record its index in the array.
                feature->addTermType(i);
                currTypeTermIndices.append(feature->getNumTerms() - 1);
            }

            // Append this array to the list that includes all types
            termIndicesByType.append(currTypeTermIndices);
        }

        // Add the specified number of random predicates
        for (int i = 0; i < numChildren; i++) {

            int pred = (int)(rand() % (domain->getNumPredicates()-1)) + 1; 
            const PredicateTemplate* predTemp = 
                    domain->getPredicateTemplate(pred);

            // Add the feature with a random mapping
            Array<int> mapping;
            for (int i = 0; i < predTemp->getNumTerms(); i++) {
                int termType = predTemp->getTermTypeAsInt(i);
                int numTerms = termIndicesByType[termType].size();
                mapping.append(termIndicesByType[termType][rand() % numTerms]);
            }

            feature->addChild(featureArray_[pred], 0.0, mapping);

            // Check to see if we have a query predicate
            if (queryPreds.find(pred) != -1) {
                queryPredChosen = true;
            }
        }
    }

    feature->setId(++maxFeatureId_);

    Array<int> mapping;
    for (int i = 0; i < feature->getNumTerms(); i++) {
        mapping.append(-i - 1);
    }
    ((RecursiveFeature*)topFeature_)->addChild(feature, 0.0, mapping);

    // Use random example to set weights
    Array<int> grounding;
    for (int i = 0; i < feature->getNumTerms(); i++) {
        int type = feature->getTermType(i);
        
        const Array<int>* constants = domain->getConstantsByType(type);
        grounding.append((*constants)[rand() % constants->size()]);
    }
    feature->setWeightsFromExample(grounding, domain->getDB(), queryPreds, 
            numChildren);
    topFeature_->setWeight(topFeature_->getNumWeights() - 1, 1.0 + frand());

    featureArray_.append(feature);
}

double RRF::getLogPseudoLikelihood(Database* db, const Array<int>& queryPreds)
{
    Array<Predicate*> allPreds;
    for (int i = 0; i < queryPreds.size(); i++) {
        Array<Predicate*> predArray;
        Predicate::createAllGroundings(queryPreds[i], db->getDomain(), 
                predArray);
        allPreds.append(predArray);
    }

    return getLogPseudoLikelihood(db, allPreds);
}


double RRF::getLogPseudoLikelihood(Database* db, 
        const Array<Predicate*>& allPreds)
{

    // LPL = log pseudo-likelihood, 
    //       the conditional probability of each query predicate 
    //       conditioned on all data.
    double lpl = 0.0;
#if 1
    invalidateAll();
    double logTruthProb = getLogValue(db);
#endif

    for (int i = 0; i < allPreds.size(); i++) {
#if 1
        TruthValue currValue = db->getValue(allPreds[i]);
        assert(currValue != UNKNOWN);
        TruthValue newValue = (currValue == TRUE) ? FALSE : TRUE;
        db->setValue(allPreds[i], newValue);
        changedPredicate(allPreds[i], db);
        double logUntruthProb = getLogValue(db);
        db->setValue(allPreds[i], currValue);
        changedPredicate(allPreds[i], db);

        //cout << log(truthProb / (truthProb + untruthProb));

        if (fabs(logUntruthProb - logTruthProb) > 100) {
            lpl += logUntruthProb - logTruthProb;
        } else {
            lpl += -log(1.0 + exp(logUntruthProb - logTruthProb));
        }
#else
        lpl += getPredicateLogLikelihood(allPreds[i], db);
#endif
    }

    return lpl/allPreds.size();
}


double RRF::getExactZ(const Domain* domain, const Array<int>& queryPreds,
        Database* origDb)
{
    ArraysAccessor<TruthValue> predValues;
    Array<TruthValue> truthValues;
    truthValues.append(TRUE);
    truthValues.append(FALSE);

    Array<Predicate*> allPreds;
    for (int i = 0; i < queryPreds.size(); i++) {
        Array<Predicate*> predArray;
        Predicate::createAllGroundings(queryPreds[i], domain, predArray);
        allPreds.append(predArray);
    }

    for (int i = 0; i < allPreds.size(); i++) {
        predValues.appendArray(&truthValues);
    }

    double Z = 0.0;

    // Create database, using closed-world assumption
    Array<bool> closedWorld;
    for (int i = 0; i < domain->getNumPredicates(); i++) {
        closedWorld.append(true);
    }
    Database* db = new Database(domain, closedWorld, true);

    // Fill in values from the original db
    if (origDb != NULL) {
        for (int i = 1; i < domain->getNumPredicates(); i++) {
            Array<Predicate*> predArray;
            Predicate::createAllGroundings(i, domain, predArray);
            for (int j = 0; j < predArray.size(); j++) {
                db->setValue(predArray[j], origDb->getValue(predArray[j]));
            }
            predArray.deleteItemsAndClear();
        }
    }

    // Sum value for all worlds (keeping the evidence pred values the same)
    do {
        Array<TruthValue> truthValues;
        predValues.getNextCombination(truthValues);
        for (int i = 0; i < truthValues.size(); i++) {
            db->setValue(allPreds[i], truthValues[i]);
        }

        Array<int> emptyGrounding;
        Z += topFeature_->getValue(emptyGrounding, db);
    } while (predValues.hasNextCombination());

    allPreds.deleteItemsAndClear();

    delete db;

    return Z;
}

double RRF::getExactZ(const Domain* domain) 
{
    Array<int> allPreds;
    for (int i = 1; i < domain->getNumPredicates(); i++) {
        allPreds.append(i);
    }
    return getExactZ(domain, allPreds, NULL);
}

inline bool isnamebegin(char c) {
    return (isalpha(c) || c == '_');
}

inline bool isname(char c) {
    return (isalnum(c) || c == '_');
}

char* nextChild(char* buf, const Array<char*>& vars, 
        double& weight, string& name, Array<int>& grounding)
{
    char* pbuf;
    bool isConstant = false;

    // Read weight
    while (isspace(*buf)) { buf++; }
    pbuf = buf;
    while (!isspace(*pbuf) && *pbuf != ')' && *pbuf) { pbuf++; }
    if (*pbuf == ')') {
        isConstant = true;
    } 
    *pbuf = '\0';
    weight = atof(buf);
    buf = pbuf + 1;

#if PARSE_DEBUG
    cout << "Read weight: " << weight << endl;
#endif

    while (!isConstant && isspace(*buf)) { buf++; }
    if (*buf == ')') {
        isConstant = true;
    }
    if (!isnamebegin(*buf)) {
        isConstant = true;
    }

    if (!isConstant) {
        // Read name
        char* pbuf = buf;
        while (isname(*pbuf)) { pbuf++; }
        if (*pbuf != '(') {
            cout << "ERROR: expected '('; found '" << *pbuf << "'.\n";
            return NULL;
        }
        // Mark the end of the name with a NULL and save it to a string
        *pbuf = '\0';
        name = buf;
        buf = pbuf + 1;

#if PARSE_DEBUG
        cout << "Read name: " << name << endl;
#endif

        // Read arguments.  No spaces allowed in between!
        while (isspace(*buf)) { buf++; }
        if (*buf != ')') {
            pbuf = buf;
            bool argsDone = false;
            while (!argsDone) {
                while (isname(*pbuf)) { pbuf++; }
                if (*pbuf != ',' && *pbuf != ')') {
                    cout << "ERROR: expected ','; found '" << *pbuf << "'.\n";
                    return NULL;
                }

                if (*pbuf == ')') {
                    argsDone = true;
                }

                // Mark the end of the argument with a NULL
                *pbuf = '\0';
                int varIndex = -1;
                for (int i = 0; i < vars.size(); i++) {
                    if (!strcmp(buf, vars[i])) {
                        varIndex = i;
                        break;
                    }
                }
                grounding.append(varIndex);

                // Advance past the NULL, to the next argument
                pbuf++;
                buf = pbuf;
            }
        } else {
            // Advance past ')'
            buf++;
        }
    } else {
        name = "";
    }

    // Advance to next child, so we're ready to parse it.
    // Return NULL if there are no more children.
    while (isspace(*buf)) { buf++; }
    if (*buf == '+') {
        return ++buf;
    } else if (*buf == ')' || *buf == '\n' || *buf == '\r' || *buf == '\0') {
        return NULL;
    } else {
        // TODO: distinguish between NULL and error?
        cout << "ERROR: read \'" << *buf << "\'; expected '+' or ')'.\n";
        return NULL;
    }
}


// Utility function
void readFeature(Array<string>& types, Array<double>& weights, 
        Array<string>& children, Array<Array<int> >& mapping, 
        string& name, char* buf)
{
    types.clear();
    weights.clear();
    children.clear();
    mapping.clear();

    // ASSUMPTION: no spaces in between arguments
    char* pbuf;

    // Get feature name
    while (isspace(*buf)) { buf++; }
    pbuf = buf;
    while (isname(*pbuf)) { pbuf++; }
    if (*pbuf != '(') { 
        cout << "ERROR: expected ')'; read '" << *pbuf << "'.\n";
        return;
    }
    *pbuf = '\0';
    name = buf;

#if PARSE_DEBUG
    cout << "Read feature name: " << name << endl;
#endif

    // Read in feature arguments
    pbuf++;
    buf = pbuf;
    Array<char*> vars;
    bool done = false;
    while (isspace(*pbuf)) { ++pbuf; }
    if (*pbuf == ')') {
        done = true;
    }

    // Read name/type pairs
    while (*buf && !done) {

        char* p2buf;

        //
        // READ NAME
        //
        while (isname(*pbuf))  { ++pbuf; }

        // Read spaces, ':', and more spaces
        p2buf = pbuf;
        while (isspace(*p2buf)) { ++p2buf; }
        if (*p2buf != ':') {
            cout << "ERROR: expected ':'; read '" << *p2buf << "'\n";
            return;
        }
        p2buf++;
        while (isspace(*p2buf)) { ++p2buf; }

        // Save name in array and update buffer pointers
        *pbuf = '\0';
        vars.append(buf);
        pbuf = buf = p2buf;

        //
        // READ TYPE
        //
        while(isname(*pbuf))  { ++pbuf; }

        // Read spaces, ',' or ')', and more spaces
        p2buf = pbuf;
        while(isspace(*p2buf)) { ++p2buf; }
        if ((*p2buf != ',') && (*p2buf != ')')) {
            cout << "ERROR: expected ',' or ')'; read '" << *p2buf << "'.\n";
            return;
        }
        if (*p2buf == ')') {
            done = true;
        }
        p2buf++;
        while (isspace(*p2buf)) { ++p2buf; }

        // Save type in array and update buffer pointers
        *pbuf = '\0';
        types.append(buf);
        pbuf = buf = p2buf;
    }

    // Skip the " = exp(" part
    while (*buf && *buf != '(') { buf++; }
    buf++;

    while (buf) {

        double weight;
        Array<int> grounding;
        string childName;

        // Read the next child
        buf = nextChild(buf, vars, weight, childName, grounding);

        weights.append(weight);
        children.append(childName);
        mapping.append(grounding);
    }
}



void RRF::load(istream& in, Domain* domain)
{
#if PARSE_DEBUG
    cout << "Loading RRF from file...\n";
#endif

    // HACK -- this is fragile, and presumes a very rigid input format.
    // TODO -- fix and generalize.

    // Create predicate (bottom level) features
    // HACK -- ignore equality predicate, 0
    for (int pred = 1; pred < domain->getNumPredicates(); pred++) {
        PredicateFeature* predFeat = new PredicateFeature(
                domain->getPredicateTemplate(pred));
        maxFeatureId_++;
        predFeat->setId(maxFeatureId_);
        featureArray_.append(predFeat);
    }

    // Create a single constant feature (only one is needed)
    {
        ConstantFeature* conFeature = new ConstantFeature(1.0);
        maxFeatureId_++;
        conFeature->setId(maxFeatureId_);
        featureArray_.append(conFeature);
    }

#if PARSE_DEBUG
    cout << "Created predicate features.\n";
#endif
    int firstNonPredFeature = maxFeatureId_+1;


    char buffer[10240];
    char* buf = buffer;

    // Skip first line, which should be "RRF {\n"
    in.getline(buf, 10240);

    // Keep track of children, so feature can refer to later features
    // as necessary, rather than being defined in a specific order.
    Array<Array<string> > allChildren;
    Array<Array<double> > allWeights;
    Array<Array<Array<int> > > allGroundings;

    // Each of the following lines should contain a feature.
    // Read in and set all feature weights.
    while (1) {

        Array<string> terms;
        Array<double> weights;
        Array<string> children;
        Array<Array<int> > groundings;
        string name;

        in.getline(buf, 10240);

        // Stop if we hit a closing '}'
        while (isspace(*buf))     { buf++; }
        if (buf[0] == '}' || !in) { break; }

        readFeature(terms, weights, children, groundings, name, buf);
        RecursiveFeature* feat;
        if (!name.compare("root") || !name.compare("f0")) {
            feat = new RecursiveFeature(name.c_str(), true, false);
            topFeature_ = feat;
        } else {
            feat = new RecursiveFeature(name.c_str(), false, true);
            //feat = new ClausalFeature(name.c_str());
        }
        feat->setId(++maxFeatureId_);
        featureArray_.append(feat);
        for (int j = 0; j < terms.size(); j++) {
            int type;
            type = atoi(terms[j].c_str());
            if (type < 1) {
                type = domain->getTypeId(terms[j].c_str());
            }
            if (type < 1) {
                cout << "ERROR: Unknown term type '" << terms[j] 
                    << "' for feature '" << name << "'.\n";
                type = 1;
            }
            feat->addTermType(type);
        }

#if PARSE_DEBUG
        cout << "Read feature: " << name << endl;
        cout << "Terms: ";
        for (int j = 0; j < terms.size(); j++) {
            cout << terms[j] << endl;
        }
#endif

        // We have to postpone this, because the referenced features may not
        // have been defined yet.
        allChildren.append(children);
        allWeights.append(weights);
        allGroundings.append(groundings);
    }

    // After reading in all features, construct parent/child relationships
    // TODO: check for cycles!
    for (int i = 0; i < getNumFeatures() - firstNonPredFeature; i++) {
        for (int j = 0; j < allChildren[i].size(); j++) {
            Feature* child = getFeature(allChildren[i][j].c_str());
            if (child == NULL) {
                cout << "ERROR: unknown feature " << allChildren[i][j] <<
                    " referenced by feature " << 
                    getFeature(i + firstNonPredFeature)->getName() << ".\n";
            } 
            RecursiveFeature* parent = (RecursiveFeature*)getFeature(
                    i + firstNonPredFeature);
            if (parent == NULL) {
                cout << "ERROR: could not find index " 
                    << (i + firstNonPredFeature) << " feature!\n";
            }
            parent->addChild(child, allWeights[i][j], allGroundings[i][j]);
        }
    }

    // DEBUG
    cout << this << endl;
    return;
}
