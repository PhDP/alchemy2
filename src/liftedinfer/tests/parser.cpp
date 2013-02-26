#include "parser.h"
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
using namespace std;
#include "cleanuputils.h"
#include "stringconversionutils.h"

LMParser::~LMParser()
{
	cleanup(domainList);
	predicateDomainMap.clear();
}

bool LMParser::isTermConstant(string term)
{
	//if starts with a capital letter or number, it is taken as a constant
	if((term[0] >= 65 && term[0] <= 90)|| (term[0] >= 48 && term[0] <= 57))
		return true;
	else
		return false;
}

LvrTerm* LMParser::create_new_term(int domainSize)
{
	//map domain to an integer domain for ease of manipulation
	vector<int> iDomain(domainSize);
	for(int k=0;k<domainSize;k++)
		iDomain[k] = k;
	//create a new term
	LvrTerm* term = new LvrTerm(0,iDomain);
	return term;
}

//TO DO when displaying convert internal integer representation of domain value to string
string LMParser::convert_atom_to_string(Atom* atom)
{
	//get the index into the domain list
	//vector<int> domainIndexList = predicateDomainMap[atom->symbol->id];
	
	//for(int i=0;i<atom->terms.size();i++)
	//{
		//int domainIndex = domainIndexList[i];
	//}
	return string(" ");
}

WClause* LMParser::create_new_clause(vector<int> predicateSymbolIndex,vector<bool> sign,
						   vector<vector<LvrTerm*> > iTermsList)
{
	int numAtoms = predicateSymbolIndex.size();
	WClause* clause = new WClause();
	clause->atoms = vector<Atom*>(numAtoms);
	clause->satisfied =false;
	clause->sign = sign;
	for(int i=0;i<numAtoms;i++)
	{
		PredicateSymbol* symbol = LvrMLN::create_new_symbol(mln.symbols[predicateSymbolIndex[i]]);
		vector<LvrTerm*> terms = iTermsList[i];
		Atom* atom = new Atom(symbol,terms);
		clause->atoms[i]= atom;
	}
	return clause;
}

//void LMParser::parseCNFString(string formula,vector<PredicateSymbol*> predicateList)
void LMParser::parseCNFString(string formula,fstream& filestr)
{
	//apend multiple lines
	char* buf = new char[1024];
	size_t pos;
	while(1)
	{
		pos = formula.find(WEIGHTSEPARATOR);
		if(pos!=string::npos)
			break;
		else
		{
			filestr.getline(buf,1024);
			string s(buf);
			s.erase( remove(s.begin(), s.end(), ' '), s.end() );	
			formula.append(s);
		}
	}
	delete[] buf;
	//extract the weight
	string weight = formula.substr(pos+2);
	stringstream convert(weight);
	double wt;
	convert >> wt;
	formula = formula.substr(0,pos);
	//cout<<wt<<endl;
	vector<string> clauseStrings;
	LStringConversionUtils::tokenize(formula,clauseStrings,ANDOPERATOR);
	vector<WClause*> CNF;
	for(int i=0;i<clauseStrings.size();i++)
	{
		//If clause starts with a paranthesis
		if(clauseStrings[i].find(LEFTPRNTH)==0)
		{
			//eliminate the first and last paranthesis
			clauseStrings[i] = clauseStrings[i].substr(1,clauseStrings[i].size()-2);
		}
		
		vector<string> atomStrings;
		LStringConversionUtils::tokenize(clauseStrings[i],atomStrings,OROPERATOR);
		vector<vector<string> > sTermsList(atomStrings.size());
		//sign of an atom
		vector<bool> sign(atomStrings.size());
		//index into the predicate symbols
		vector<int> predicateSymbolIndex(atomStrings.size());
		
		for(int j=0;j<atomStrings.size();j++)
		{
			//find opening and closing braces
			int startpos=atomStrings[j].find(LEFTPRNTH);
			string predicateName = atomStrings[j].substr(0,startpos);
			if(predicateName.find(NOTOPERATOR)==0)
			{
				//negation
				predicateName = predicateName.substr(1,predicateName.size());
				sign[j]=true;
			}
			for(int k=0;k<mln.symbols.size();k++)
			{
				//found the predicate
				if(mln.symbols[k]->symbol.compare(predicateName)==0)
				{
					predicateSymbolIndex[j] = k;
					break;
				}
			}
			
			int endpos=atomStrings[j].find(RIGHTPRNTH);
			string termsString = atomStrings[j].substr(startpos+1,endpos-startpos-1);
			vector<string> terms;
			LStringConversionUtils::tokenize(termsString,terms,COMMASEPARATOR);
			sTermsList[j]=terms;
			//check if the number of terms is equal to the declared predicate
			if(terms.size()!=mln.symbols[predicateSymbolIndex[j]]->variable_types.size())
			{
				cout<<"Error! Number/domain of terms in the predicate delcaration does not match in formula::"<<predicateName.c_str()<<endl;
				exit(-1);
			}
		}
		//create required terms
		vector<vector<LvrTerm*> > iTermsList(atomStrings.size());
		for(int j=0;j<atomStrings.size();j++)
		{
			//for each term of atom i, check if it has already appeared in previous atoms of clause
			vector<LvrTerm*> iTerms(sTermsList[j].size());
			for(int k=0;k<sTermsList[j].size();k++)
			{
				int domainIndex = predicateDomainMap[predicateSymbolIndex[j]].at(k);
				//if term is a constant must be a unique term
				if(isTermConstant(sTermsList[j].at(k)))
				{
					//find the id of the term
					int id=-1;
					for(int m=0;m<domainList[domainIndex]->values.size();m++)
					{
						if(domainList[domainIndex]->values[m].compare(sTermsList[j].at(k))==0)
						{
							id=m;
							break;
						}
					}
					if(id==-1)
					{
						cout<<"Constant does not match predicate's domain::"<<domainList[domainIndex]->name<<endl;
						exit(-1);
					}
					iTerms[k] = new LvrTerm(0,id);
				}
				else
				{
					int domainSize = domainList[domainIndex]->values.size();
					bool isExistingTerm = false;
					int atomIndex=-1;
					int termIndex=-1;
					//check in term lists for atoms 0 to j;
					for(int m=0;m<j;m++)
					{
						for(int n=0;n<sTermsList[m].size();n++)
						{
							if((sTermsList[m])[n].compare((sTermsList[j])[k])==0)
							{
								//check if the domains of the matched variables are the same
								int atomSymbolIndex1 = predicateSymbolIndex[m];
								int atomId1 = mln.symbols[atomSymbolIndex1]->id;
								int domainListIndex1 = predicateDomainMap[atomId1].at(n);

								int atomSymbolIndex2 = predicateSymbolIndex[j];
								int atomId2 = mln.symbols[atomSymbolIndex2]->id;
								int domainListIndex2 = predicateDomainMap[atomId2].at(k);
								if(domainList[domainListIndex1]->name.compare(domainList[domainListIndex2]->name)!=0)
								{
									cout<<"Error! variables do not match type ::"<<atomStrings[j]<<"("<<domainList[domainListIndex1]->name<<", "<<domainList[domainListIndex2]->name<<")"<<endl;
									exit(-1);
								}
								//variable is repeated, use the term created for atom m, term n
								isExistingTerm = true;
								atomIndex = m;
								termIndex = n;								
								break;
							}
						}
						if(isExistingTerm)
							break;
					}
					if(isExistingTerm)
					{
						//use the terms created for previous atoms
						iTerms[k] = (iTermsList[atomIndex])[termIndex];
					}
					else
					{
						//create a new LvrTerm
						iTerms[k] = create_new_term(domainSize);
					}
				}
			}
			iTermsList[j] = iTerms;
		}//j atoms
		WClause* newClause = create_new_clause(predicateSymbolIndex,sign,iTermsList);
		newClause->weight = LogDouble(wt,false);
		CNF.push_back(newClause);
		
	}//i clauses
	int formulaStartIndex = mln.clauses.size();

	for(int i=0;i<CNF.size();i++)
	{
		mln.clauses.push_back(CNF[i]);
	}
	int formulaEndIndex = mln.clauses.size();
	mln.formulas.push_back(new Formula(formulaStartIndex,formulaEndIndex,LogDouble(wt,false)));
}

void LMParser::parseDomainString(string line,fstream& filestr)
{
	size_t pos= line.find(EQUALSTO);
	if(pos!=string::npos)
	{
		//new domain name found, create a domain string
		string domainString;
		string domainName = line.substr(0,pos);
		size_t startpos = line.find(LEFTFLOWER);
		size_t endpos = line.find(RIGHTFLOWER);
		if(startpos!=string::npos && endpos!=string::npos)
		{
			domainString.append(line.substr(startpos+1,endpos-startpos-1));
		}
		else if(startpos!=string::npos)
		{
			domainString.append(line.substr(startpos+1));
			bool done =false;
			char* buf = new char[1024];
			while(!done)
			{
				filestr.getline(buf,1024);
				string nextLine(buf);
				nextLine.erase( remove(nextLine.begin(), nextLine.end(), ' '), nextLine.end() );
				size_t endpos = nextLine.find("}");
				if(endpos!=string::npos)
				{
					domainString.append(nextLine.substr(0,endpos));
					done =true;
				}
				else
				{
					domainString.append(nextLine);
				}
			}
			delete [] buf;
		}
		vector<string> domainValues;
		LStringConversionUtils::tokenize(domainString,domainValues,COMMASEPARATOR);
		PDomain* domain=new PDomain(domainName,domainValues);
		domainList.push_back(domain);
	}
}

void LMParser::parsePredicateString(string line)
{
	size_t pos;
	int startpos = line.find(LEFTPRNTH);
	int endpos = line.find(RIGHTPRNTH);
	string symbolName = line.substr(0,startpos);
	symbolName = symbolName.substr(0, symbolName.find(LEFTPRNTH));
	string symbolTerms = line.substr(startpos+1,endpos-startpos-1);
	vector<string> termNames;
	

	LStringConversionUtils::tokenize(symbolTerms,termNames,COMMASEPARATOR);
	//create a new predicate symbol
	vector<int> var_types(termNames.size());
	for(int m=0;m<var_types.size();m++)
		var_types[m] = 0;
	PredicateSymbol* p = new PredicateSymbol(predicateId,symbolName,var_types,LogDouble(1,false),LogDouble(1,false),predicateId);
	//predicateList.push_back(p);
	mln.symbols.push_back(p);

	//Build the map for this predicate;
	//For predicateid, generate a vector of domainIds that index Domains
	vector<int> domainIndex;
	for(int i=0;i<termNames.size();i++)
	{
		int matchingIndex = -1;
		for(int j=0;j<domainList.size();j++)
		{
			if(termNames[i].compare(domainList[j]->name)==0)
			{
				matchingIndex = j;
				break;
			}
		}
		if(matchingIndex == -1)
		{
			cout<<"Error! PDomain does not exist for predicate::"<<symbolName.c_str()<<endl;
			exit(-1);
		}
		domainIndex.push_back(matchingIndex);
	}
	predicateDomainMap.insert(pair<int,vector<int> >(predicateId,domainIndex));
	//increment predicateid
	predicateId++;
}

void LMParser::parseInputMLNFile(string filename)
{
	fstream filestr(filename.c_str());
	if(filestr == NULL)
	{
		cout<<"Error reading LvrMLN File"<<endl;
		exit(-1);
	}
	char* buf = new char[1024];
	int cnt=0;
	ParsingState state;
	size_t pos;
	while(filestr)
	{
		filestr.getline(buf,1024);
		string line(buf);
		line.erase( remove(line.begin(), line.end(), ' '), line.end() );
		if(line.size()==0)
			continue;
		PDomain* domain;
		if(line.find(DOMAINSTART)!=string::npos)
		{
			//set the sate to parsing domains
			state = Domains;
			continue;
		}
		else if(line.find(PREDICATESTART)!=string::npos)
		{
			//set the sate to parsing domains
			state = Predicates;
			continue;
		}
		else if(line.find(FORMULASTART)!=string::npos)
		{
			state = Formulas;
			continue;
		}

		switch(state)
		{
			case Domains:
				{
					parseDomainString(line,filestr);
				}
				break;
			case Predicates:
				{
					parsePredicateString(line);
				}
				break;
			case Formulas:
				{
					//if all forumlas are seperated as clauses (CNF)
					parseCNFString(line,filestr);
				}
				break;
			default:
				break;

		}

	}
	delete[] buf;
	filestr.close();
#ifdef __DEBUG_PRINT__
	cout<<"*************Parser Output****************"<<endl;
	for(int i=0;i<mln.clauses.size();i++)
		mln.clauses[i]->print();
	for(int i=0;i<mln.formulas.size();i++)
		cout<<"["<<mln.formulas[i]->MLNClauseStartIndex<<","<<mln.formulas[i]->MLNClauseEndIndex<<"] "<<mln.formulas[i]->weight<<endl;
	cout<<"********************************************"<<endl;
#endif
}

bool LMParser::checkTermsValidity(int predicateId,vector<string> terms, vector<int>& matchedIndexList)
{
	if(terms.size()!=predicateDomainMap[predicateId].size())
		return false;
	for(int i=0;i<predicateDomainMap[predicateId].size();i++)
	{
		if(!isTermConstant(terms[i]))
		{
			cout<<"Error! should enter constants in DB file::"<<terms[i].c_str()<<endl;
			exit(-1);
		}
		int domainIndex = predicateDomainMap[predicateId].at(i);
		bool found=false;
		for(int j=0;j<domainList[domainIndex]->values.size();j++)
		{
			if(domainList[domainIndex]->values[j].compare(terms[i])==0)
			{
				found =true;
				matchedIndexList.push_back(j);
			}
		}
		if(!found)
			return false;
	}
	return true;
}


void LMParser::parseDB(string filename)
{
	ifstream filestr(filename.c_str());
	if(filestr == NULL)
	{
		cout<<"Error reading DB File"<<endl;
		return;
	}
	char* buf = new char[1024];
	int cnt=0;
	size_t pos;
	while(filestr)
	{
		filestr.getline(buf,1024);
		string line(buf);
		if(line.size()==0)
			continue;
		pos = line.find(LEFTPRNTH);
		bool sign =false;
		int start = 0;
		size_t posS = line.find(NOTOPERATOR);
		if(posS!=string::npos)
		{
			sign=true;
			start = posS+1;
		}
		string predicateName = line.substr(start,pos-start);
		//find the predicate id
		int id = -1;
		int index = -1;
		for(int i=0;i<mln.symbols.size();i++)
		{
			if(mln.symbols[i]->symbol.compare(predicateName)==0)
			{
				id = mln.symbols[i]->id;
				index = i;
				break;
			}
		}
		if(id==-1)
		{
			cout<<"Error! Predicate in DB file not found::"<<predicateName<<endl;
			exit(-1);
		}
		int endpos = line.find(RIGHTPRNTH);
		vector<string> terms;
		LStringConversionUtils::tokenize(line.substr(pos+1,endpos-pos-1),terms,COMMASEPARATOR);
		if(terms.size() > mln.symbols[index]->variable_types.size())
		{
			cout<<"Error! Wrong terms in Predicate in DB file::"<<predicateName<<endl;
			exit(-1);
		}
		vector<int> matchedIndexList;
		//check type of each term
		if(!checkTermsValidity(id,terms,matchedIndexList))
		{
			cout<<"Error! Wrong value of term in Predicate in DB file::"<<predicateName<<endl;
			exit(-1);
		}
		vector<LvrTerm*> termList(matchedIndexList.size());
		for(int i=0;i<matchedIndexList.size();i++)
		{
			LvrTerm* term = new LvrTerm(0,matchedIndexList[i]);
			termList[i]=term;
		}
		Atom* atom = new Atom(LvrMLN::create_new_symbol(mln.symbols[id]),termList);
		WClause* newClause = new WClause();
		newClause->atoms.push_back(atom);
		newClause->satisfied = false;
		newClause->sign.push_back(sign);
		newClause->weight =  LogDouble(0,false);
		//append the formula
		int formulaStartIndex = mln.clauses.size();
		mln.clauses.push_back(newClause);
		int formulaEndIndex = mln.clauses.size();
		//hard formula i.e. evidence
		mln.formulas.push_back(new Formula(formulaStartIndex,formulaEndIndex,LogDouble(0,false)));
	}
	delete[] buf;
	filestr.close();
}
