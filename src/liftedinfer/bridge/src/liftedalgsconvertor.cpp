#include "liftedalgsconvertor.h"
#include "queryupdater.h"

//constructor for weight learning where information needs to be converted back and forth
//between the alchemy format and the lifted inference format (I.E. WEIGHT LEARNING)
LiftedAlgsConvertor::LiftedAlgsConvertor(VariableState* state):state_(state),mln_(state_->getMLN()),
domain_(state_->getDomain())
{
	lvrMln = new LvrMLN();
	//convertInputMLNToLifted(state_->getMLN(),state_->getDomain());
	bool mlncontainsGrndClause=false;
	convertInputMLNToLifted(mln_,domain_,mlncontainsGrndClause);
	//mln_= state->getMLN();
	//domain_ = domain;
	//intilialize the query updater to contain all atoms in the state
	const GroundPredicateHashArray* grndPredHashArray = state_->getGndPredHashArrayPtr();
	vector<vector<int> > queryIntRep(grndPredHashArray->size());
	//cout<<grndPredHashArray->size()<<endl;
	for(int i=0;i<grndPredHashArray->size();i++)
	{
		//const Array<unsigned int*> intRep = (*grndPredHashArray)[i]->getIntArrRep();
		int sz = (*grndPredHashArray)[i]->getIntArrRep()->size();
		vector<int> tmp(sz);
		for(int j=0;j<sz;j++)
		{
			tmp[j]=(*((*grndPredHashArray)[i]->getIntArrRep()))[j];
		}
		queryIntRep[i] = tmp;
		//(*grndPredHashArray)[i]->print(cout,domain_);
		//cout<<endl;
		//LFileDump::dumpQueriesToFile(queryIntRep);
	}
	//no need for query strings or results file since never written to file
	vector<string> qstrings;
 	LvrQueryUpdater::createInstance(queryIntRep,qstrings,"");
 	liftedAlgsHandler = NULL;
}

//constructor for marginal estimations where information only needs to be converted
//from alchemy format to lifted inference format
LiftedAlgsConvertor::LiftedAlgsConvertor(MLN* mln, Domain* domain,GroundPredicateHashArray queries):mln_(mln),
		domain_(domain)
{
	state_ = NULL;
	lvrMln = NULL;
	queries_ = queries;
	liftedAlgsHandler = NULL;
}

LiftedAlgsConvertor::~LiftedAlgsConvertor()
{
	delete lvrMln;
	if(liftedAlgsHandler)
		delete liftedAlgsHandler;
}

void LiftedAlgsConvertor::doConversion(vector<vector<int> >& queriesIntRep, vector<vector<int> > & evidenceIntRep, 
	vector<string>& queryStrings,bool removeSelfJoins)
{
	lvrMln = new LvrMLN();
	bool mlncontainsGrndClause=false;
	convertInputMLNToLifted(mln_,domain_,mlncontainsGrndClause);
	//process the evidence
	addEvidence(evidenceIntRep);
	//LFileDump::dumpEvidenceToFile(evidenceIntRep);
	if(evidenceIntRep.size() > 0)
	{
		cout<<"Normalizing evidence..."<<endl;
		lvrMln->preprocessEvidence(evidenceIntRep);
		//for(int i=0;i<lvrMln->clauses.size();i++)
			//lvrMln->clauses[i]->print();
		//cout<<endl;
	}
	else if(mlncontainsGrndClause)
	{
		//the input mln has a "+" term, needs to be normalized
		LNormalizer ln(*lvrMln);
		ln.normalizeClauses(lvrMln->clauses);
		for(unsigned int i=0;i<lvrMln->clauses.size();i++)
		{
			for(unsigned int j=0;j<lvrMln->clauses[i]->atoms.size();j++)
			{
				lvrMln->clauses[i]->atoms[j]->symbol->parentId =
						lvrMln->clauses[i]->atoms[j]->symbol->id;
			}
		}
	}
	//cout<<"Converted MLN"<<endl;
	//for(int i=0;i<lvrMln->clauses.size();i++)
		//lvrMln->clauses[i]->print();
	//cout<<endl;
	queriesIntRep.resize(queries_.size());
	for(int i=0;i<queries_.size();i++)
	{
		for(int j=0;j<queries_[i]->getIntArrRep()->size();j++)
		{
			queriesIntRep[i].push_back((*queries_[i]->getIntArrRep())[j]);
		}
		string s=queries_[i]->getPredString(domain_);
		queryStrings.push_back(s);
	}
	//LFileDump::dumpQueriesToFile(queriesIntRep);
	if(removeSelfJoins)
	{
		bool isSelfJoined = false;
		for(int i=0;i<lvrMln->clauses.size();i++)
		{
			if(lvrMln->clauses[i]->isSelfJoined())
			{
				isSelfJoined = true;
				break;
			}
		}
		if(isSelfJoined)
		{
			cout<<"Preprocessing normalized MLN..."<<endl;
			processSelfJoins();
		}
	}


	liftedAlgsHandler = new LiftedAlgsHandler(*lvrMln);
	//LFileDump::dumpMLNToFile(lvrMln);
	cout<<"done writing dump files ("<<MLNDUMPFILENAME<<","<<SYMBOLSDUMPFILENAME<<","<<
	QUERIESDUMPFILENAME<<","<<EVIDENCEDUMPFILENAME<<")"<<endl;
	//for(int i=0;i<lvrMln->clauses.size();i++)
		//lvrMln->clauses[i]->print();
}

void LiftedAlgsConvertor::updateMLNWeights(vector<double> weights)
{
	for(unsigned int i=0;i<lvrMln->clauses.size();i++)
	{
		//get the new weight based on the original clause index
		lvrMln->clauses[i]->weight =
				LogDouble(weights[lvrMln->clauses[i]->originalClauseIndex],false);
	}
}


void LiftedAlgsConvertor::doLBGMar(LvrParams* params)
{
	vector<vector<int> > evidenceIntRep;
	vector<vector<int> > queriesIntRep;
	vector<string> queryStrings;
	doConversion(queriesIntRep,evidenceIntRep,queryStrings,true);
	if(params->operationMode != ECLUSTERING)
	{
		//create the singleton instance for queries
		LvrQueryUpdater::createInstance(queriesIntRep,queryStrings,
				params->resultFile,&evidenceIntRep);
	}
	//obtain the string representations of all queries to display
	
	liftedAlgsHandler->LBGApproxMarginals(params);
}

void LiftedAlgsConvertor::doLISApproxZ(LvrParams* params)
{
	vector<vector<int> > evidenceIntRep;
	vector<vector<int> > queriesIntRep;
	vector<string> queryStrings;
	doConversion(queriesIntRep,evidenceIntRep,queryStrings,false);
	liftedAlgsHandler->LISApproxPartition(params);
}

void LiftedAlgsConvertor::doLISApproxMar(LvrParams* params)
{
	vector<vector<int> > evidenceIntRep;
	vector<vector<int> > queriesIntRep;
	vector<string> queryStrings;
	doConversion(queriesIntRep,evidenceIntRep,queryStrings,false);
	//create the singleton instance for queries
	LvrQueryUpdater::createInstance(queriesIntRep,queryStrings,
			params->resultFile,&evidenceIntRep);
	liftedAlgsHandler->LISApproxMarginals(params);

}

void LiftedAlgsConvertor::doWMPTPApproxZ(LvrParams* params)
{
	vector<vector<int> > evidenceIntRep;
	lvrMln = new LvrMLN();
	bool mlncontainsGrndClause=false;
	convertInputMLNToLifted(mln_,domain_,mlncontainsGrndClause);
	//process the evidence
	addEvidence(evidenceIntRep);
	//cout<<"Converted MLN"<<endl;
	//for(int i=0;i<lvrMln->clauses.size();i++)
		//lvrMln->clauses[i]->print();
	//cout<<endl;
	liftedAlgsHandler = new LiftedAlgsHandler(*lvrMln);
	liftedAlgsHandler->WMPTPZApprox(params);
}

void LiftedAlgsConvertor::doWMPTPExactZ(LvrParams* params)
{
	vector<vector<int> > evidenceIntRep;
	lvrMln = new LvrMLN();
	bool mlncontainsGrndClause=false;
	convertInputMLNToLifted(mln_,domain_,mlncontainsGrndClause);
	//process the evidence
	addEvidence(evidenceIntRep);
	//cout<<"Converted MLN"<<endl;
	//for(int i=0;i<lvrMln->clauses.size();i++)
		//lvrMln->clauses[i]->print();
	//cout<<endl;
	vector<vector<int> > queriesIntRep(queries_.size());
	vector<string> queryStrings;
	for(int i=0;i<queries_.size();i++)
	{
		for(int j=0;j<queries_[i]->getIntArrRep()->size();j++)
		{
			queriesIntRep[i].push_back((*queries_[i]->getIntArrRep())[j]);
		}
		string s=queries_[i]->getPredString(domain_);
		queryStrings.push_back(s);
	}
	//LFileDump::dumpQueriesToFile(queriesIntRep);
	if(queriesIntRep.size() > 0)
	{
		//initialize the query updater
		LvrQueryUpdater::createInstance(queriesIntRep,queryStrings,
				params->resultFile,&evidenceIntRep);
	}
	liftedAlgsHandler = new LiftedAlgsHandler(*lvrMln);
	liftedAlgsHandler->WMPTPZExact(params);
}


PredicateSymbol* LiftedAlgsConvertor::createPredicateSymbol(Predicate* pred)
{
	int predicateId = pred->getId();
	vector<int> var_types(pred->getNumTerms());
	string symbolName(pred->getName());
	PredicateSymbol* pSymbol = new PredicateSymbol(predicateId,symbolName,var_types,LogDouble(1,false),LogDouble(1,false),predicateId);
	return pSymbol;
}

void LiftedAlgsConvertor::addEvidence(vector<vector<int> >& evidenceIntRep)
{
	set<int> completedIds;
	for(int i=0;i<mln_->getNumClauses();i++)
	{
		for(int j=0;j<mln_->getClause(i)->getNumPredicates();j++)
		{
			if(completedIds.count(mln_->getClause(i)->getPredicate(j)->getId()) == 0)
			{
				completedIds.insert(mln_->getClause(i)->getPredicate(j)->getId());
				Array<Predicate*>* indexedGndings = new Array<Predicate*>();
				domain_->getDB()->getTrueGndings(mln_->getClause(i)->getPredicate(j)->getId()
						,indexedGndings);
				for(int k=0;k<indexedGndings->size();k++)
				{
					GroundPredicate* gpTmp = new GroundPredicate((*indexedGndings)[k]);
					vector<int> tmpArr(gpTmp->getIntArrRep()->size());
					for(unsigned int jj=0;jj<tmpArr.size();jj++)
						tmpArr[jj] = (*gpTmp->getIntArrRep())[jj];
					evidenceIntRep.push_back(tmpArr);
					PredicateSymbol* ps = createPredicateSymbol((*indexedGndings)[k]);
					vector<LvrTerm*> atomTerms((*indexedGndings)[k]->getNumTerms());
					for (int termno = 0; termno < (*indexedGndings)[k]->getNumTerms(); termno++)
					{
							vector<int> domValues(1);
							const Term* term = (*indexedGndings)[k]->getTerm(termno);
							domValues[0] = term->getId();
							LvrTerm* lvrTerm = new LvrTerm(0,domValues);
							atomTerms[termno]=lvrTerm;
					}
					Atom* atom = new Atom(ps,atomTerms);
					WClause* nClause = new WClause();
					nClause->atoms.push_back(atom);
					nClause->sign.push_back(false);
					lvrMln->clauses.push_back(nClause);
					//evidenceAtoms.push_back(atom);
				}
			}
		}
	}
}


WClause* LiftedAlgsConvertor::convertClauseToLifted(Clause* clause, const Domain* domain, vector<PredicateSymbol*>& symbolsToAdd,bool& containsGroundedClause)
{
   const Array<Predicate*>* preds = clause->getPredicates();
   map<int,LvrTerm*> mapOfTerms;
   WClause* nClause = new WClause();
   nClause->atoms = vector<Atom*>(preds->size());
   nClause->sign = vector<bool>(preds->size());
   nClause->weight = LogDouble(clause->getWt(),false);
   set<int> completedIds;
	for (int i = 0; i < preds->size(); i++)
	{
		nClause->sign[i] = !((*preds)[i]->getSense());
		int predicateId = (*preds)[i]->getId();
		vector<int> var_types((*preds)[i]->getNumTerms());
		string symbolName((*preds)[i]->getName());
		PredicateSymbol* pSymbol = new PredicateSymbol(predicateId,symbolName,var_types,LogDouble(1,false),LogDouble(1,false),predicateId);
		const PredicateTemplate *ptemplate = (*preds)[i]->getTemplate();
		vector<LvrTerm*> atomTerms((*preds)[i]->getNumTerms());
		for (int termno = 0; termno < (*preds)[i]->getNumTerms(); termno++)
		{
			 const Term* term = (*preds)[i]->getTerm(termno);
			 if(term->getType()==Term::CONSTANT)
			 {
				 int id = term->getId();
				 vector<int> domValues(1);
				 domValues[0] = id;
				 LvrTerm* lvrTerm = new LvrTerm(0,domValues);
				 atomTerms[termno]=lvrTerm;
				 //contains a "+" in the term
				 containsGroundedClause=true;
				 continue;
			 }
			 int varId = term->getId();
			 map<int,LvrTerm*>::iterator mapIt = mapOfTerms.find(varId);
			 if(mapIt==mapOfTerms.end())
			 {
				 //new term, add it
				 int varTypeId = ptemplate->getTermTypeAsInt(termno);
				 const Array<int>* constants = domain->getConstantsByType(varTypeId);
				 vector<int> domValues;
				 for(unsigned int j=0;j<constants->size();j++)
				 {
					 domValues.push_back((*constants)[j]);
				 }
				 LvrTerm* lvrTerm = new LvrTerm(0,domValues);
				 mapOfTerms.insert(pair<int,LvrTerm*> (varId,lvrTerm));
				 atomTerms[termno]=lvrTerm;
			 }
			 else
			 {
				 //use existing term
				 atomTerms[termno] = mapIt->second;
			 }
		}
		//sort the domains of all terms
		for(unsigned int jj=0;jj<atomTerms.size();jj++)
			sort(atomTerms[jj]->domain.begin(),atomTerms[jj]->domain.end());
		Atom* atom = new Atom(pSymbol,atomTerms);
		nClause->atoms[i] = atom;
		if(completedIds.count(pSymbol->id) == 0)
		{
			//ADD TO LvrMLN
			PredicateSymbol* nSymbol = new PredicateSymbol(predicateId,symbolName,var_types,LogDouble(1,false),LogDouble(1,false),predicateId);
			symbolsToAdd.push_back(nSymbol);
		}
	}
	return nClause;
}

void LiftedAlgsConvertor::convertInputMLNToLifted(const MLN* mln, const Domain* domain,
		bool& containsGroundedClause)
{
  cout<<"Converting MLN..."<<endl;
  const FormulaAndClausesArray* fca = mln->getFormulaAndClausesArray();
  int clauseIndex = 0;
  for (int i = 0; i < fca->size(); i++)
  {
	IndexClauseHashArray* indexClauses = (*fca)[i]->indexClauses;
	for (int j = 0; j < indexClauses->size(); j++)
	{
	  int idx = (*indexClauses)[j]->index;
	  Clause* c = (*indexClauses)[j]->clause;
	  vector<PredicateSymbol*> symbolsToAdd;
	  WClause* lCls = convertClauseToLifted(c,domain,symbolsToAdd,containsGroundedClause);
	  lCls->originalClauseIndex = clauseIndex++;
	  lvrMln->clauses.push_back(lCls);
	  for(unsigned int jj=0;jj<symbolsToAdd.size();jj++)
		  lvrMln->symbols.push_back(symbolsToAdd[jj]);
	}
  }
}

void LiftedAlgsConvertor::processSelfJoins()
{
	LNormalizer ln(*lvrMln);
	ln.removeNPSelfJoins(lvrMln->clauses);
}

void LiftedAlgsConvertor::updateState()
{
	const Domain* domain = state_->getDomain();
	const GroundPredicateHashArray* grndPredHashArray = state_->getGndPredHashArrayPtr();
	for(int i=0;i<grndPredHashArray->size();i++)
	{
		//const Array<unsigned int*> intRep = (*grndPredHashArray)[i]->getIntArrRep();
		int idx = state_->getIndexOfGroundPredicate((*grndPredHashArray)[i]);
		if(idx==-1)
			continue;
		int blockIdx = domain->getBlock((*grndPredHashArray)[i]);
		int sz = (*grndPredHashArray)[i]->getIntArrRep()->size();
		vector<int> tmp(sz);
		for(int j=0;j<sz;j++)
		{
			tmp[j]=(*((*grndPredHashArray)[i]->getIntArrRep()))[j];
		}
		int val = LvrQueryUpdater::Instance()->getCurrentSampleValue(tmp);
		if(val==-1)
			continue;
		bool bVal=false;
		if(val > 0)
			bVal=true;
		if(bVal)
			cout<<1;
		else
			cout<<0;
	    state_->setValueOfAtom(idx+1,bVal,false,blockIdx);
	}
	cout<<endl;
}
void LiftedAlgsConvertor::processWeightLearningInput(LvrParams* params,bool initialization)
{
	params->isWeightLearning = true;
	/*
	if(params->algorithm == LISM)
	{
		if(!liftedAlgsHandler)
			liftedAlgsHandler = new LiftedAlgsHandler(*lvrMln);
		if(initialization && params->samplingMode == EINFORMED)
		{
			liftedAlgsHandler->buildProposal();
		}
		liftedAlgsHandler->LISApproxMarginals(params);
	}
	 */
	if(params->algorithm == LVG)
	{
		//preprocess self joins
		if(initialization)
		{
			bool isSelfJoined = false;
			for(int i=0;i<lvrMln->clauses.size();i++)
			{
				if(lvrMln->clauses[i]->isSelfJoined())
				{
					isSelfJoined = true;
					break;
				}
			}
			if(isSelfJoined)
			{
				cout<<"Preprocessing converted MLN..."<<endl;
				processSelfJoins();
			}
		}
		params->operationMode = ESAMPLER;
		if(!liftedAlgsHandler)
			liftedAlgsHandler = new LiftedAlgsHandler(*lvrMln);
		liftedAlgsHandler->LBGApproxMarginals(params);
	}
}

