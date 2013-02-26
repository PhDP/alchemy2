#include "lvrmln.h"
#include "cleanuputils.h"
#include "normalizer.h"
#include "hashalgorithm.h"


#include <sstream>
using namespace std;



void Atom::print(bool verbose)
{
	if(verbose)
	{
		cout<<symbol->symbol.c_str()<<"["<< symbol->id<<","<<symbol->parentId<<"{"<<symbol->pweight<<","<<symbol->nweight<<"}]"<<" (";
		for(unsigned int j=0;j<terms.size();j++){
			if((int)terms[j]->domain.size()==1){
				cout<<terms[j]->domain[0];
			}
			else{
				cout<<terms[j]<<"[ #"<<terms[j]->domain.size()<<" ]";
			}
			if(j!=terms.size()-1)
				cout<<",";
		}
		cout<<")";
	}
	else
	{
		cout<<symbol->symbol.c_str()<<symbol->id<<"(";
		for(unsigned int j=0;j<terms.size();j++){
			if((int)terms[j]->domain.size()==1){
				cout<<terms[j]->domain[0];
			}
			else{
				cout<<terms[j]<<"#"<<terms[j]->domain.size();
			}
			/*if(terms[j]->origDomainSize!=-1)
			{
				cout<<"D#"<<terms[j]->origDomainSize;
			}
			*/
			if(j!=terms.size()-1)
				cout<<",";
		}
		cout<<")";


	}
}


bool Atom::isConstant()
{
	bool constantDomain = true;
	for(unsigned int i=0;i<terms.size();i++)
	{
		if(terms[i]->domain.size() > 1)
		{
			constantDomain = false;
			break;
		}
	}
	return constantDomain;
}

bool Atom::isSingletonAtom(int& variableIndex)
{
	/*
	if(terms.size() == 1)
	{
		variableIndex = 0;
		return true;
	}
	*/

	variableIndex = -1;
	int numVariables = 0;
	for(unsigned int i=0;i<terms.size();i++)
	{
		if(terms[i]->domain.size() > 1)
		{
			//variable
			numVariables++;
			if(numVariables > 1)
				break;
			//store index of the first variable
			if(variableIndex == -1)
				variableIndex = i;
		}
	}
	//singleton if no variables or 1 variable
	if(numVariables == 1)
		return true;
	else
		return false;
}

int Atom::getNumberOfGroundings()
{
	int numberOfGroundings=1;
	for(unsigned int i=0;i<terms.size();i++)
	{
		numberOfGroundings *= terms[i]->domain.size();
	}
	return numberOfGroundings;
}


void WClause::findSelfJoins(map<int,vector<int> >& selfJoinedAtoms)
{
	for(unsigned int i=0;i<atoms.size();i++)
	{
		int id = atoms[i]->symbol->id;
		if(selfJoinedAtoms.find(id) == selfJoinedAtoms.end())
		{
			vector<int> tempPos;
			tempPos.push_back(i);
			for(unsigned int j=i+1;j<atoms.size();j++)
			{
				if(atoms[i]->symbol->id == atoms[j]->symbol->id)
				{
					//self joined
					tempPos.push_back(j);
				}
			}
			if(tempPos.size()>1)
			{
				//self joined
				selfJoinedAtoms.insert(pair<int,vector<int> >(id,tempPos));
			}
		}
	}
}

void WClause::find_decomposer(vector<LvrTerm*>& terms, vector<vector<int> >& positions)
{
	map<int,vector<int> > selfJoinedAtoms;
	findSelfJoins(selfJoinedAtoms);
	set<LvrTerm*> clause_terms;

	for(unsigned int i=0;i<atoms.size();i++){
		//if any constant in clause, can never have a decomposer
		if(atoms[i]->isConstant())
			return;
		for(unsigned int j=0;j<atoms[i]->terms.size();j++){
			clause_terms.insert(atoms[i]->terms[j]);
		}
	}

	vector<vector<int> > clause_positions(clause_terms.size());
	for(unsigned int i=0;i<atoms.size();i++){
		if(atoms[i]->isConstant())
			continue;
		for(unsigned int j=0;j<atoms[i]->terms.size();j++){
			int k=distance(clause_terms.begin(),clause_terms.find(atoms[i]->terms[j]));
			clause_positions[k].push_back(j);
		}
	}
	int k=0;
	for(set<LvrTerm*>::iterator i=clause_terms.begin();i!=clause_terms.end();i++){
		//check if term is a variable and not a constant
		if((*i)->domain.size() > 1)
		{
			// Check if the term participates in all atoms and unifies with self joined atoms. If it does, it is decomposer
			if (clause_positions[k].size()==atoms.size()){
				//check for self joins
				bool unified = true;
				for(map<int,vector<int> >::iterator it = selfJoinedAtoms.begin();it!=selfJoinedAtoms.end();it++)
				{
					//check for unification with all self joined atoms
					
					int pos = clause_positions[k].at(it->second[0]);
					for(unsigned int m=1;m<it->second.size();m++)
					{
						if(clause_positions[k].at(m)!=pos)
						{
							unified = false;
							break;
						}
					}
					if(!unified)
						break;
				}
				if(unified)
				{
					terms.push_back(*i);
					positions.push_back(clause_positions[k]);
				}
			}
		}
		++k;
	}

}


WClause::~WClause()
{
	//each terms reference count
	set<LvrTerm*> deletionRef;
	for(unsigned int i=0;i<atoms.size();i++)
	{
		for(unsigned int j=0;j<atoms[i]->terms.size();j++)
		{
			if(deletionRef.count(atoms[i]->terms[j]) > 0)
			{
				atoms[i]->terms.erase(atoms[i]->terms.begin()+j);
				j--;
			}
			else
			{
				deletionRef.insert(atoms[i]->terms[j]);
				removeItem(atoms[i]->terms,j);
				j--;
			}
		}
	}

	cleanup(atoms);
	
	sign.clear();
}

int WClause::getNumberOfGroundedClauses()
{
	set<LvrTerm*> completedTerms;
	int totalSize = 1;
	for(unsigned int j=0;j<atoms.size();j++)
	{
		for(unsigned int k=0;k<atoms[j]->terms.size();k++)
		{
			completedTerms.insert(atoms[j]->terms[k]);
		}
	}
	for(set<LvrTerm*>::iterator it=completedTerms.begin();it!=completedTerms.end();it++)
	{
		totalSize *= (*it)->domain.size();
	}
	return totalSize;
}

void WClause::removeAtom(int index)
{
	//terms are shared in the same clause, check if the atom has shared terms
	for(unsigned int i=0;i<atoms[index]->terms.size();i++)
	{
		bool shared = false;
		for(unsigned int j=0;j<atoms.size();j++)
		{
			if(j==index)
				continue;
			for(unsigned int k=0;k<atoms[j]->terms.size();k++)
			{
				if(atoms[index]->terms[i] == atoms[j]->terms[k])
				{
					shared = true;
					break;
				}
			}
			if(shared)
				break;
		}
		if(!shared)
		{
			removeItem(atoms[index]->terms,i);
			i--;
		}
	}
	removeItem(atoms,index);
	sign.erase(sign.begin()+index);
}

void WClause::removeAtomWithWeightUpdation(int index)
{
	//terms are shared in the same clause, check if the atom has shared terms
	for(unsigned int i=0;i<atoms[index]->terms.size();i++)
	{
		bool shared = false;
		for(unsigned int j=0;j<atoms.size();j++)
		{
			if(j==index)
				continue;
			for(unsigned int k=0;k<atoms[j]->terms.size();k++)
			{
				if(atoms[index]->terms[i] == atoms[j]->terms[k])
				{
					shared = true;
					break;
				}
			}
			if(shared)
				break;
		}
		if(!shared)
		{
			//update clause weight
			weight = weight*LogDouble(atoms[index]->terms[i]->domain.size(),false);
			removeItem(atoms[index]->terms,i);
			i--;
		}
	}
	removeItem(atoms,index);
	sign.erase(sign.begin()+index);
}

bool WClause::isSelfJoinedOnAtom(Atom* atom)
{
	int count=0;
	for(unsigned int i=0;i<atoms.size();i++)
	{
		if(atoms[i]->symbol->id==atom->symbol->id)
			count++;
	}
	if(count > 1)
		return true;
	else
		return false;

}


Atom::~Atom()
{
	//delete the symbol
	delete symbol;
}

LvrTerm* LvrMLN::create_new_term(LvrTerm* term)
{
	int type = term->type;
	vector<int> domain;
	for(unsigned int k=0;k<term->domain.size();k++)
		domain.push_back(term->domain[k]);
	return new LvrTerm(type,domain);
}

Atom* LvrMLN::create_new_atom(Atom* atom)
{
	vector<LvrTerm*> terms;
	for(unsigned int i=0;i<atom->terms.size();i++)
	{
		LvrTerm* newTerm = create_new_term(atom->terms[i]);
		newTerm->origDomainSize = atom->terms[i]->origDomainSize;
		for(unsigned int j=0;j<atom->terms[i]->origDomain.size();j++)
			newTerm->origDomain.push_back(atom->terms[i]->origDomain[j]);
		terms.push_back(newTerm);
	}
	//Atom* newAtom = new Atom(atom->symbol,terms);
	Atom* newAtom = new Atom(create_new_symbol(atom->symbol),terms);
	return newAtom;

}

PredicateSymbol* LvrMLN::create_new_symbol(PredicateSymbol* symbol)
{
	vector<int> var_types;
	for(unsigned int i=0;i<symbol->variable_types.size();i++)
		var_types.push_back(symbol->variable_types[i]);
	PredicateSymbol* newSymbol = new PredicateSymbol(symbol->id,symbol->symbol,var_types,
		symbol->pweight,symbol->nweight,symbol->normParentId);
	//USED IN CACHING
	newSymbol->parentId = symbol->parentId;
	return newSymbol;
}

void LvrMLN::copyAllClauses(vector<WClause*> origClauses, vector<WClause*>& newClauses)
{
	for(unsigned int i=0;i<origClauses.size();i++)
	{
		WClause* newClause = create_new_clause(origClauses[i]);
		newClauses.push_back(newClause);
	}
}

bool WClause::isPropositional()
{
	for(unsigned int i=0;i<atoms.size();i++)
	{
		if(!atoms[i]->isConstant())
			return false;
	}
	return true;
}


LvrMLN::~LvrMLN()
{
	cleanup(symbols);
	cleanup(clauses);
	//clauses.clear();
	cleanup(formulas);
}

void LvrMLN::preprocessEvidenceWMPTP()
{
	LNormalizer ln(*this);
	//TODO-integrate WM ptp
	ln.normalizeClauses(clauses,true,true);
	vector<WClause*> evidence;
	vector<int> positions;
	for(unsigned int i=0;i<clauses.size();i++)
	{
		//mln.clauses[i]->print();
		//cout<<mln.clauses[i]->weight<<endl;
		if(clauses[i]->weight.is_zero && clauses[i]->atoms.size()==1)
		{
			evidence.push_back(clauses[i]);
			positions.push_back(i);
		}
	}
	
	for(unsigned int i=0;i<clauses.size();i++)
	{
		if(find(positions.begin(),positions.end(),i)!=positions.end())
			continue;
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			for(unsigned int k=0;k<evidence.size();k++)
			{
				if(evidence[k]->atoms[0]->symbol->id == clauses[i]->atoms[j]->symbol->id)
				{
					//check sign
					if(evidence[k]->sign[0] == clauses[i]->sign[j])
					{
						clauses[i]->satisfied = true;
					}
					clauses[i]->removeAtom(j);
					j--;
					break;
				}

			}
		}
	}
	
	for(unsigned int i=0;i<positions.size();i++)
	{
		removeItem(clauses,positions[i]-i);
	}

	//vector<bool> completedIds(getMaxPredicateId());	
	cleanup(symbols);
	//create a new set of symbols
	set<int> completedSymbIds;
	int maxId = 0;
	for(unsigned int i=0;i<clauses.size();i++)
	{
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			//if(!completedIds[clauses[i]->atoms[j]->symbol->id])
			if(completedSymbIds.count(clauses[i]->atoms[j]->symbol->id)==0)
			{
				//completedIds[clauses[i]->atoms[j]->symbol->id] = true;
				completedSymbIds.insert(clauses[i]->atoms[j]->symbol->id);
				symbols.push_back(LvrMLN::create_new_symbol(clauses[i]->atoms[j]->symbol));
				if(clauses[i]->atoms[j]->symbol->id > maxId)
					maxId = clauses[i]->atoms[j]->symbol->id;

			}
		}
	}
	max_predicate_id = maxId+1;
}
void LvrMLN::putWeightsOnClauses()
{
	for(unsigned int i=0;i<formulas.size();i++)
	{
		//cout<<mln.formulas[i]->weight<<endl;
		//FOR NOW EACH FORMULA HAS EXACTLY ONE CLAUSE
		clauses[i]->weight = formulas[i]->weight;
	}
}
WClause* LvrMLN::create_new_clause(WClause* clause) {
	WClause* new_clause = new WClause();
	new_clause->sign = clause->sign;
	new_clause->satisfied = clause->satisfied;
	//NEW
	new_clause->weight = clause->weight;
	//added for weight learning efficiency
	new_clause->originalClauseIndex = clause->originalClauseIndex;

	//if atoms have common terms their relationship must be maintained when new clause is created
	vector<LvrTerm*> newTerms;
	vector<LvrTerm*> oldTerms;
	for (size_t i = 0; i < clause->atoms.size(); i++) 
	{
		for (size_t j = 0; j < clause->atoms[i]->terms.size(); j++)
		{
			int termPosition=-1;
			for(unsigned int m=0;m<oldTerms.size();m++)
			{
				if(oldTerms[m]==clause->atoms[i]->terms[j])
				{
					termPosition = m;
				}
			}
			if(termPosition==-1)
			{
				LvrTerm* term = new LvrTerm();
				term->type = clause->atoms[i]->terms[j]->type;
				term->origDomainSize = clause->atoms[i]->terms[j]->origDomainSize;
				for(unsigned int k=0;k<clause->atoms[i]->terms[j]->domain.size();k++)
					term->domain.push_back(clause->atoms[i]->terms[j]->domain[k]);
				for(unsigned int k=0;k<clause->atoms[i]->terms[j]->origDomain.size();k++)
					term->origDomain.push_back(clause->atoms[i]->terms[j]->origDomain[k]);
				newTerms.push_back(term);
				oldTerms.push_back(clause->atoms[i]->terms[j]);
			}
			else
			{
				newTerms.push_back(newTerms[termPosition]);
				oldTerms.push_back(clause->atoms[i]->terms[j]);
			}
		}
	}
	int ind=0;
	new_clause->atoms = vector<Atom*>(clause->atoms.size());
	for (size_t i = 0; i < new_clause->atoms.size(); i++) {
		new_clause->atoms[i] = new Atom();
		new_clause->atoms[i]->symbol = create_new_symbol(clause->atoms[i]->symbol);
		new_clause->atoms[i]->terms = vector<LvrTerm*>(clause->atoms[i]->terms.size());
		for (size_t j = 0; j < new_clause->atoms[i]->terms.size(); j++) {
			new_clause->atoms[i]->terms[j] = newTerms[ind];
			ind++;
		}
	}
	return new_clause;
}

void LvrMLN::preprocessEvidence(vector<vector<int> > evidenceIntRep)
{
	map<int,int> evidenceHashMap;
	for(unsigned int i=0;i<evidenceIntRep.size();i++)
	{
		evidenceHashMap.insert(pair<int,int>(LvrHashAlgorithm::DJBHash(evidenceIntRep[i]),1));
	}
	doPreprocessEvidence(evidenceHashMap);
}

void LvrMLN::preprocessEvidence(vector<Atom*> evidenceAtoms)
{
	map<int,int> evidenceHashMap;
	for(unsigned int i=0;i<evidenceAtoms.size();i++)
	{
		int key = LvrHashAlgorithm::convertToHash(evidenceAtoms[i]);
		evidenceHashMap.insert(pair<int,int>(key,1));
	}
	doPreprocessEvidence(evidenceHashMap);
}


void LvrMLN::doPreprocessEvidence(map<int,int> evidenceHashMap)
{
	LNormalizer ln(*this);
	ln.normalizeClauses(clauses);
	for(unsigned int i=0;i<clauses.size();i++)
	{
		vector<int> clauseHashCodes;
		int removed=0;
		for(int j=0;j<clauses[i]->atoms.size();j++)
		{
			if(!clauses[i]->atoms[j]->isConstant())
				continue;
			int key = LvrHashAlgorithm::convertToHash(clauses[i]->atoms[j]);
			if(evidenceHashMap.find(key)!=evidenceHashMap.end())
			{
				//evidence atom is found
				if(!clauses[i]->sign[j])
					clauses[i]->satisfied = true;
				clauses[i]->removeAtom(j);
				j--;
			}
		}
		if(clauses[i]->satisfied)
			continue;

		//remove all self joins
		vector<int> removeInds;
		for(int j=0;j<clauses[i]->atoms.size();j++)
		{
			if(!clauses[i]->atoms[j]->isConstant())
				continue;
			for(int k=0;k<clauses[i]->atoms.size();k++)
			{
				if(k==j || !clauses[i]->atoms[k]->isConstant())
					continue;
				if(clauses[i]->atoms[j]->symbol->id==clauses[i]->atoms[k]->symbol->id 
					&& clauses[i]->sign[j]!=clauses[i]->sign[k])
				{
					clauses[i]->satisfied;
					break;
				}
			}
			if(clauses[i]->satisfied)
				break;
		}
	}
	/*
	//For each atom, append id to its predicate symbol string
	for(unsigned int i=0;i<clauses.size();i++)
	{
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
		{
			//clauses[i]->atoms[j]->symbol->parentId = clauses[i]->atoms[j]->symbol->id;
			stringstream st;
			st<<clauses[i]->atoms[j]->symbol->id;
			clauses[i]->atoms[j]->symbol->symbol.append(st.str());
		}
		if(clauses[i]->atoms.size()==0 && clauses[i]->satisfied && clauses[i]->weight.is_zero)
		{
			removeItem(clauses,i);
			i--;
		}
	}
	*/
	//remove all clauses with 0 weight and empty
	for(unsigned int i=0;i<clauses.size();i++)
	{
		//make each atom's parent as itself
		for(unsigned int j=0;j<clauses[i]->atoms.size();j++)
			clauses[i]->atoms[j]->symbol->parentId = clauses[i]->atoms[j]->symbol->id;
		if(clauses[i]->atoms.size()==0 && clauses[i]->satisfied && clauses[i]->weight.is_zero)
		{
			removeItem(clauses,i);
			i--;
		}
	}
}



