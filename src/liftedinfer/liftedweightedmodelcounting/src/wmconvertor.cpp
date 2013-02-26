#include "wmconvertor.h"
#include <sstream>
using namespace std;

	
AtomLocation* LWMConvertor::getMatchedLocation(vector<AtomLocation*> atomLocs, vector<AtomLocation*> termLocs,
	WClause* clause)
{
	for(int n=0;n<termLocs.size();n++)
	{
		int iter =0;
		for(int k=0;k<clause->atoms.size();k++)
		{
			for(int m=0;m<clause->atoms[k]->terms.size();m++)
			{
				if(termLocs.at(n)->isEqual(atomLocs.at(iter++)))
				{
					//clause index is not really needed, can return -1 in its place
					return (new AtomLocation(-1,k,m));
				}
			}
		}
	}
	return NULL;
}
void LWMConvertor::appendNewPredicate(vector<LvrTerm*> newTerms, vector<vector<AtomLocation*> > termLocs, 
	vector<vector<AtomLocation*> > atomLocs,vector<WClause*>& clauses,PredicateSymbol* ps)
{	
	for(int i=0;i<clauses.size();i++)
	{
		vector<LvrTerm*> newSymbolTerms(newTerms.size());
		//for each new term, match its location in clause i's atom locations
		for(int j=0;j<termLocs.size();j++)
		{
			/*AtomLocation* match = NULL;
			for(int n=0;n<termLocs[j].size();n++)
			{
				int iter =0;
				for(int k=0;k<clauses[i]->atoms.size();k++)
				{
					for(int m=0;m<clauses[i]->atoms[k]->terms.size();m++)
					{
						if(termLocs[j].at(n)->isEqual(atomLocs[i].at(iter++)))
						{
							match = new AtomLocation(i,k,m);
							break;
						}
						if(match)
							break;
					}
					if(match)
						break;
				}
				if(match)
					break;
			}
			*/
			AtomLocation* match = getMatchedLocation(atomLocs[i],termLocs[j],clauses[i]);
			if(match == NULL)
			{
				//use a new copy of the term j, no reuse required for new predicate symbol
				newSymbolTerms[j] = LvrMLN::create_new_term(newTerms[j]);
			}
			else
			{
				newSymbolTerms[j] = clauses[i]->atoms[match->atomIndex]->terms[match->termIndex];				
			}
		}
		
		//Atom* newSymAtom = new Atom(ps,newSymbolTerms);
		Atom* newSymAtom = new Atom(LvrMLN::create_new_symbol(ps),newSymbolTerms);
		clauses[i]->atoms.push_back(newSymAtom);
		clauses[i]->sign.push_back(false);
	}
	
}

void LWMConvertor::appendNewPredicate1(vector<LvrTerm*> newTerms, vector<vector<AtomLocation*> > termLocs,vector<WClause*>& clauses,PredicateSymbol* ps)
{
	for(int i=0;i<clauses.size();i++)
	{
		vector<LvrTerm*> newSymbolTerms(newTerms.size());
		//for each new term, match its location in clause i's atom locations
		for(int j=0;j<termLocs.size();j++)
		{
			//check if term occurs in clause i
			AtomLocation* foundLocation = NULL;
			for(int k=0;k<termLocs[j].size();k++)
			{				
				//check which element matches current location
				for(int m=0;m<clauses[i]->atoms.size();m++)
				{
					for(int n=0;n<clauses[i]->atoms[m]->terms.size();n++)
					{
						AtomLocation newLocation(i,m,n);
						if(termLocs[j].at(k)->isEqual(&newLocation))
						{
							//found the correct term location
							foundLocation = termLocs[j].at(k);
							break;
						}
					}
					if(foundLocation!=NULL)
						break;
				}
			}
			if(foundLocation == NULL)
			{
				//use a new copy of the term j, no reuse required for new predicate symbol
				newSymbolTerms[j] = LvrMLN::create_new_term(newTerms[j]);
			}
			else
			{
				//new predicate symbol needs to reuse term from found atom location
				newSymbolTerms[j] = clauses[i]->atoms[foundLocation->atomIndex]->terms[foundLocation->termIndex];
			}
		}
		Atom* newSymAtom = new Atom(LvrMLN::create_new_symbol(ps),newSymbolTerms);
		clauses[i]->atoms.push_back(newSymAtom);
		clauses[i]->sign.push_back(true);
	}
}



void LWMConvertor::getNewAtomTerms(int formulaIndex, vector<LvrTerm*>& terms,vector<vector<AtomLocation*> >& locationIndex)
{
	int start = mln.formulas[formulaIndex]->MLNClauseStartIndex;
	int end = mln.formulas[formulaIndex]->MLNClauseEndIndex;
	for(int i=start;i<end;i++)
	{	
		int relLocation = i - start;
		vector<LvrTerm*> completedTerms;
		int currIndex = locationIndex.size();
		for(int j=0;j<mln.clauses[i]->atoms.size();j++)
		{
			for(int k=0;k<mln.clauses[i]->atoms[j]->terms.size();k++)
			{
				LvrTerm* term = mln.clauses[i]->atoms[j]->terms[k];
				bool found =false;
				for(int m=0;m<completedTerms.size();m++)
				{
					if(completedTerms[m]==term)
					{
						found=true;
						//locationIndex[currIndex+m].push_back(new AtomLocation(i,j));
						locationIndex[currIndex+m].push_back(new AtomLocation(relLocation,j,k));
						break;
					}
				}
				if(!found)
				{
					completedTerms.push_back(term);
					//AtomLocation* atLoc = new AtomLocation(i,j);
					AtomLocation* atLoc = new AtomLocation(relLocation,j,k);
					vector<AtomLocation*> tempLoc;
					tempLoc.push_back(atLoc);
					locationIndex.push_back(tempLoc);
				}
			}
		}
		for(int k=0;k<completedTerms.size();k++)
			terms.push_back(completedTerms[k]);
	}
}

//recursive function to convert !(C1 ^ C2 ^...CN) to CNF
vector<WClause*> LWMConvertor::generateAtomCombinations(int index,int startIndex,vector<vector<AtomLocation*> >& locations)
{
	if(index == startIndex)
	{
		vector<WClause*> tmpClauses;
		for(int i=0;i<mln.clauses[index]->atoms.size();i++)
		{
			Atom* newAtom = LvrMLN::create_new_atom(mln.clauses[index]->atoms[i]);
			vector<AtomLocation*> aLocVec;
			for(int k=0;k<newAtom->terms.size();k++)
			{
				AtomLocation* aLoc = new AtomLocation(index-startIndex,i,k);
				aLocVec.push_back(aLoc);
			}
			locations.push_back(aLocVec);
			WClause* newClause = new WClause();
			newClause->atoms.push_back(newAtom);
			newClause->sign.push_back(!(mln.clauses[index]->sign[i]));
			newClause->satisfied = mln.clauses[index]->satisfied;
			tmpClauses.push_back(newClause);
		}
		return tmpClauses;
	}

	vector<WClause*> allCombClauses = generateAtomCombinations(index-1,startIndex,locations);
	vector<WClause*> augmentedClauses;
	
	for(int i=0;i<mln.clauses[index]->atoms.size();i++)
	{
		for(int j=0;j<allCombClauses.size();j++)
		{
			WClause* newClause = LvrMLN::create_new_clause(allCombClauses[j]);
			//make a new copy of location vector associated with j
			vector<AtomLocation*> newLocVec;
			for(int k=0;k<locations[j].size();k++)
				newLocVec.push_back(locations[j].at(k));
			//augment the combination
			Atom* nAtom = LvrMLN::create_new_atom(mln.clauses[index]->atoms[i]);
			bool sign = !(mln.clauses[index]->sign[i]);
			newClause->atoms.push_back(nAtom);
			newClause->sign.push_back(sign);
			augmentedClauses.push_back(newClause);

			//augment the location vector
			for(int k=0;k<nAtom->terms.size();k++)
			{
			AtomLocation* aLoc = new AtomLocation(index-startIndex,i,k);
			newLocVec.push_back(aLoc);
			}
			locations.push_back(newLocVec);
		}
	}
	locations.erase(locations.begin(),locations.begin()+allCombClauses.size());
	return augmentedClauses;
}
int LWMConvertor::toSymbolIndex(int predId)
{
	for(int i=0;i<mln.symbols.size();i++)
	{
		if(mln.symbols[i]->id == predId)
			return i;
	}
	return 0;
}

void LWMConvertor::updateMLN(int formulaIndex,int symId,vector<WClause*>& newClausestoAdd)
{
	vector<LvrTerm*> newAtomTerms;
	vector<vector<AtomLocation*> > newAtomTermLocs;
	//gather the terms for new predicate and their locations in the original LvrMLN
	getNewAtomTerms(formulaIndex,newAtomTerms,newAtomTermLocs);
	int newAtomSize = 1;
	for(unsigned int i=0;i<newAtomTerms.size();i++)
		newAtomSize *= newAtomTerms[i]->domain.size();
	double v = exp(mln.clauses[formulaIndex]->weight.value);
	double v1 = newAtomSize*v;
	conversionFactor *= LogDouble(v1,true);
	vector<vector<AtomLocation*> > atomLocs;
	int formulaEndIndex = mln.formulas[formulaIndex]->MLNClauseEndIndex-1;
	int formulaStartIndex = mln.formulas[formulaIndex]->MLNClauseStartIndex;
	//create new clauses which represents (C1 ^ C2...)->NSM [ !(C1 ^ C2...)v NSM ]
	vector<WClause*> tempClauses = generateAtomCombinations(formulaEndIndex,formulaStartIndex,atomLocs);
	int symbolIndex = toSymbolIndex(symId);
	appendNewPredicate(newAtomTerms, newAtomTermLocs, atomLocs,tempClauses,mln.symbols[symbolIndex]);
	
	//create n new clauses: NSM->(C1 ^ C2...) [!NSM V C1, !NSM V C2,...]
	vector<WClause*> tempClauses1(formulaEndIndex-formulaStartIndex+1);
	int iter = 0;
	for(int i=formulaStartIndex;i<=formulaEndIndex;i++)
	{
		tempClauses1[iter++]=LvrMLN::create_new_clause(mln.clauses[i]);
	}
	appendNewPredicate1(newAtomTerms, newAtomTermLocs,tempClauses1,mln.symbols[symbolIndex]);
	
	for(int i=0;i<tempClauses.size();i++)
	{
		newClausestoAdd.push_back(tempClauses[i]);
	}

	for(int i=0;i<tempClauses1.size();i++)
	{
		newClausestoAdd.push_back(tempClauses1[i]);
	}

}

void LWMConvertor::convertMLN()
{
	mln.preprocessEvidenceWMPTP();
	int originalMLNSize = mln.clauses.size();
	vector<WClause*> newClausestoAdd;
	mln.formulas.clear();
	for(int i=0;i<mln.clauses.size();i++)
	{
		Formula* f = new Formula(i,i+1,mln.clauses[i]->weight);
		mln.formulas.push_back(f);
	}
	for(int i=0;i<mln.clauses.size();i++)
	{
		set<LvrTerm*> terms;
		//get the weight for this formula;
		//LogDouble wt = mln.clauses[i]->weight;
		double v = exp(mln.clauses[i]->weight.value);
		double v1 = exp(-1*v);
		LogDouble wt = LogDouble(v1,false);
		for(int k=0;k<mln.clauses[i]->atoms.size();k++)
		{
			for(int m=0;m<mln.clauses[i]->atoms[k]->terms.size();m++)
			{
				terms.insert(mln.clauses[i]->atoms[k]->terms[m]);
			}
		}
		//create a new predicate symbol for each formula
		//int curr_num_symbols = mln.symbols.size();
		int curr_num_symbols = mln.getMaxPredicateId();
		vector<int> var_type(terms.size());
		for(int m=0;m<var_type.size();m++)
			var_type[m]=0;
		//new symbol has weight(NSM)=1 and weight(!NSM)=weight of formula
		string newString(NEWPREDICATEPREFIX);
		stringstream wtStr;
		wtStr<<curr_num_symbols;
		newString.append(wtStr.str());
		PredicateSymbol* p = new PredicateSymbol(curr_num_symbols, newString, var_type,LogDouble(1,false),wt,mln.getMaxNormId());
		//mln.setMaxPredicateId(mln.getMaxPredicateId()+1);
		mln.symbols.push_back(p);
		mln.setMaxPredicateId(mln.getMaxPredicateId()+1);
		updateMLN(i,curr_num_symbols,newClausestoAdd);
	}
	//erase all LvrMLN clauses and formulas
	mln.clauses.clear();
	mln.formulas.clear();
	//update the mln with augmented clauses
	mln.clauses = newClausestoAdd;
	for(int jj=0;jj<mln.clauses.size();jj++)
	{
		if(mln.clauses[jj]->satisfied)
		{
			for(int kk=0;kk<mln.clauses[jj]->atoms.size();kk++)
			{
				if(mln.clauses[jj]->atoms[kk]->symbol->symbol.find(NEWPREDICATEPREFIX)!=string::npos)
				{
					mln.clauses[jj]->removeAtom(kk);
					break;
				}
			}
		}
	}
#ifdef __DEBUG_PRINT__
	cout<<"*************WMCONVERTOR Output****************"<<endl;
	for(int i=0;i<mln.clauses.size();i++)
		mln.clauses[i]->print();
	cout<<"********************************************"<<endl;
#endif
	mln.conversionFactor = conversionFactor;
}

