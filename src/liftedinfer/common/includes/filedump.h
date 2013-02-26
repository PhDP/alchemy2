#ifndef __LFILEDUMP
#define __LFILEDUMP
#include "lvrmln.h"
#include<sstream>
#include <fstream>
using namespace std;
#include "stringconversionutils.h"

#define MLNDUMPFILENAME "mlndump.dat"
#define	SYMBOLSDUMPFILENAME "symbolsdump.dat"
#define	QUERIESDUMPFILENAME "queriesdump.dat"
#define	EVIDENCEDUMPFILENAME "evidencedump.dat"
#define PROPOSALDUMPFILENAME "proposaldump.dat"
#define PROPOSALISOLATEDTERMSDUMPFILE "proposalITdump.dat"
#define CLAUSESDUMPFILE "clausesdump.dat"

struct LSymbolDump
{
	int sign;
	int nTerms;
	int id;
	int pid;
	int nid;
	string name;
	double pval;
	double nval;
};

struct LFileDump
{
	static void printClause(ofstream* out,WClause* clause)
	{
		if(clause->satisfied==true)
			(*out)<<1<<":";
		else
			(*out)<<0<<":";
		(*out)<<clause->atoms.size()<<":";
		for(unsigned int i=0;i<clause->atoms.size();i++)
		{
			if(clause->sign[i])
				(*out)<<0<<"#";
			else
				(*out)<<1<<"#";
			(*out)<<clause->atoms[i]->terms.size()<<"#"<<clause->atoms[i]->symbol->id<<"#"<<clause->atoms[i]->symbol->parentId<<"#"<<
				clause->atoms[i]->symbol->normParentId<<"#"<<clause->atoms[i]->symbol->symbol.c_str()<<"#"<<clause->atoms[i]->symbol->pweight<<"#"<<
				clause->atoms[i]->symbol->nweight;
			if(i<clause->atoms.size()-1)
				(*out)<<"|";
		}
		(*out)<<"@";
		vector<LvrTerm*> terms;
		for(unsigned int i=0;i<clause->atoms.size();i++)
		{
			for(unsigned int j=0;j<clause->atoms[i]->terms.size();j++)
			{
				if(find(terms.begin(),terms.end(),clause->atoms[i]->terms[j])==terms.end())
				{
					terms.push_back(clause->atoms[i]->terms[j]);
				}
			}
		}

		for(unsigned int i=0;i<terms.size();i++)
		{
			//print the term domain
			for(unsigned int jj=0;jj<terms[i]->domain.size();jj++)
			{
				(*out)<<terms[i]->domain[jj];
				if(jj<terms[i]->domain.size()-1)
					(*out)<<",";
			}
			(*out)<<"&";
			string atomTermsStr;
			for(unsigned int j=0;j<clause->atoms.size();j++)
			{
				string s;
				for(unsigned int k=0;k<clause->atoms[j]->terms.size();k++)
				{
					if(terms[i]==clause->atoms[j]->terms[k])
					{
						string s1;
						stringstream st1(s1);
						st1<<j;
						string s2;
						stringstream st2(s2);
						st2<<k;
						s.append(st1.str());
						s.append(":");
						s.append(st2.str());
						s.append(",");
						//(*out)<<j<<":"<<k;
						//(*out)<<",";
					}
				}
				if(s.length()>0)
					atomTermsStr.append(s);
			}
			(*out)<<atomTermsStr.substr(0,atomTermsStr.length()-1);
			if(i<terms.size()-1)
				(*out)<<"%";
		}
		(*out)<<"@";
		(*out)<<clause->weight;
	}
	static void dumpSymbols(LvrMLN* mln,string symbolsdumpfile)
	{
		ofstream* out = new ofstream(symbolsdumpfile.c_str());
		(*out)<<mln->getMaxPredicateId()<<endl;
		(*out)<<mln->getMaxDegree()<<endl;
		for(unsigned int i=0;i<mln->symbols.size();i++)
		{
			(*out)<<-1<<"#";
			(*out)<<mln->symbols[i]->variable_types.size()<<"#"<<mln->symbols[i]->id<<"#"<<mln->symbols[i]->parentId<<"#"<<
				mln->symbols[i]->normParentId<<"#"<<mln->symbols[i]->symbol.c_str()<<"#"<<mln->symbols[i]->pweight<<"#"<<
				mln->symbols[i]->nweight;
			(*out)<<endl;
		}
		out->close();
		delete out;
	}

	static void dumpMLNToFile(LvrMLN* mln,string mlndumpfile=MLNDUMPFILENAME,string symbolsdumpfile = SYMBOLSDUMPFILENAME)
	{
		ofstream* out = new ofstream(mlndumpfile.c_str());
		
		//write all clauses in a specific format
		for(unsigned int i=0;i<mln->clauses.size();i++)
		{
			printClause(out,mln->clauses[i]);
			(*out)<<endl;
		}
		out->close();
		dumpSymbols(mln,symbolsdumpfile);
		delete out;
	}

	static void dumpClausesToFile(vector<WClause*> clauses,string clausesdumpfile = CLAUSESDUMPFILE)
	{
		ofstream* out = new ofstream(clausesdumpfile.c_str());
		
		//write all clauses in a specific format
		for(unsigned int i=0;i<clauses.size();i++)
		{
			printClause(out,clauses[i]);
			(*out)<<endl;
		}
		out->close();
		delete out;
	}

	static void readDumpedclausesFromFile(vector<WClause*>& clauses,string clausesdumpfile = CLAUSESDUMPFILE)
	{
		fstream filestr(clausesdumpfile.c_str());
		if(filestr == NULL)
		{
			cout<<"Error reading LvrMLN File"<<endl;
			exit(-1);
		}
		//char* buf = new char[1024];
		while(filestr)
		{
			//filestr.getline(buf,1024);
			string line;
			getline(filestr,line);
			if(line.size()==0)
				continue;
			bool sat = false;
			//check for empty clauses
			{
				if(line[2] == '0')
				{
					if(line[0]=='1')
					{
						sat=true;
					}
					string str = line.substr(6);
					stringstream sstr(str);
					double wt;
					sstr >> wt;
					WClause* empClause = new WClause();
					empClause->satisfied = sat;
					empClause->weight = LogDouble(wt,false);
					clauses.push_back(empClause);
					continue;
				}
			}

			vector<string> tokens;
			LStringConversionUtils::tokenize(line,tokens,"@");
			int pos = tokens[0].find(":");
			int sVal;
			stringstream str(tokens[0].substr(0,pos));
			str>>sVal;
			if(sVal==1)
				sat=true;
			//string remaining = tokens[0].substr(0,pos+1);
			//string nAtomsS = tokens[0].substr(0,pos);
			int p1 = tokens[0].rfind(":");

			string atomsS = tokens[0].substr(p1+1);
			vector<string> atomNames;
			LStringConversionUtils::tokenize(atomsS,atomNames,"|");
			vector<LSymbolDump*> sList;
			vector<PredicateSymbol*> symbols;
			for(unsigned int i=0;i<atomNames.size();i++)
			{
				LSymbolDump* ns = new LSymbolDump();
				convertToDump(atomNames[i],ns,symbols);
				sList.push_back(ns);
			}

			vector<string> toksTerms;
			LStringConversionUtils::tokenize(tokens[1],toksTerms,"%");
			vector<vector<LvrTerm*> > atomTerms(sList.size());
			for(unsigned int i=0;i<atomTerms.size();i++)
			{
				atomTerms[i].resize(sList[i]->nTerms);
			}
			for(unsigned int i=0;i<toksTerms.size();i++)
			{
				int pos = toksTerms[i].find("&");
				string s1 = toksTerms[i].substr(0,pos);
				string s2 = toksTerms[i].substr(pos+1);
				vector<string> dValsS;
				vector<int> dVals;
				LStringConversionUtils::tokenize(s1,dValsS,",");
				LStringConversionUtils::toIntArr(dValsS,dVals);
				LvrTerm* term = new LvrTerm(0,dVals);

				vector<string> postoks;
				LStringConversionUtils::tokenize(s2,postoks,",");
				for(unsigned int j=0;j<postoks.size();j++)
				{
					int pos = postoks[j].find(":");
					string s1 = postoks[j].substr(0,pos);
					string s2 = postoks[j].substr(pos+1);
					int apos = LStringConversionUtils::toInt(s1);
					int tpos = LStringConversionUtils::toInt(s2);
					atomTerms[apos].at(tpos) = term;
				}
			}
			WClause* clause = new WClause();
			//create the clause
			for(unsigned int i=0;i<atomTerms.size();i++)
			{
				Atom* atom = new Atom(symbols[i],atomTerms[i]);
				clause->atoms.push_back(atom);
				if(sList[i]->sign == 0)
					clause->sign.push_back(true);
				else
					clause->sign.push_back(false);
			}
			clause->weight = LogDouble(LStringConversionUtils::toDouble(tokens[2]),false);
			clause->satisfied=sat;
			clauses.push_back(clause);
		}
		filestr.close();
	}

	static void dumpQueriesToFile(vector<vector<int> > intRep,string queriesdumpfile = QUERIESDUMPFILENAME)
	{
		ofstream* out = new ofstream(queriesdumpfile.c_str());
		for(unsigned int i=0;i<intRep.size();i++)
		{
			for(unsigned int j=0;j<intRep[i].size();j++)
			{
				(*out)<<intRep[i].at(j);
				if(j<intRep[i].size()-1)
					(*out)<<" ";
			}
			(*out)<<endl;
		}
		out->close();
		delete out;
	}

	static void readQueriesFromFile(vector<vector<int> >& intRep,string queriesdumpfile = QUERIESDUMPFILENAME )
	{
		fstream filestr(queriesdumpfile.c_str());
		if(filestr == NULL)
		{
			cout<<"Error reading LvrMLN File"<<endl;
			exit(-1);
		}
		char* buf = new char[1024];
		while(filestr)
		{
			filestr.getline(buf,1024);
			string line(buf);
			if(line.size()==0)
				continue;
			vector<string> tokens;
			LStringConversionUtils::tokenize(line,tokens," ");
			vector<int> tmpArr;
			LStringConversionUtils::toIntArr(tokens,tmpArr);
			intRep.push_back(tmpArr);
		}
	}
	static void dumpEvidenceToFile(vector<vector<int> > intRep,string evidencedumpfile = EVIDENCEDUMPFILENAME)
	{
		ofstream* out = new ofstream(evidencedumpfile.c_str());
		for(unsigned int i=0;i<intRep.size();i++)
		{
			for(unsigned int j=0;j<intRep[i].size();j++)
			{
				(*out)<<intRep[i].at(j);
				if(j<intRep[i].size()-1)
					(*out)<<" ";
			}
			(*out)<<endl;
		}
		out->close();
		delete out;
	}

	static void readEvidenceFromFile(vector<vector<int> >& intRep,string evidencedumpfile = EVIDENCEDUMPFILENAME)
	{
		fstream filestr(evidencedumpfile.c_str());
		if(filestr == NULL)
		{
			cout<<"Error reading LvrMLN File"<<endl;
			exit(-1);
		}
		char* buf = new char[1024];
		while(filestr)
		{
			filestr.getline(buf,1024);
			string line(buf);
			if(line.size()==0)
				continue;
			vector<string> tokens;
			LStringConversionUtils::tokenize(line,tokens," ");
			vector<int> tmpArr;
			LStringConversionUtils::toIntArr(tokens,tmpArr);
			intRep.push_back(tmpArr);
		}
	}

	static void convertToDump(string s, LSymbolDump* sd,vector<PredicateSymbol*>& symbols)
	{
		vector<string> params;
		LStringConversionUtils::tokenize(s,params,"#");
		sd->sign = LStringConversionUtils::toInt(params[0]);
		sd->nTerms = LStringConversionUtils::toInt(params[1]);
		sd->id = LStringConversionUtils::toInt(params[2]);
		sd->pid = LStringConversionUtils::toInt(params[3]);
		sd->nid = LStringConversionUtils::toInt(params[4]);
		sd->name = params[5];
		sd->pval = LStringConversionUtils::toDouble(params[6]);
		sd->nval = LStringConversionUtils::toDouble(params[7]);
		vector<int> var_types(sd->nTerms);
		PredicateSymbol *ps = new PredicateSymbol(sd->id,sd->name,var_types,LogDouble(sd->pval,false),LogDouble(sd->nval,false),sd->nid);
		ps->parentId = sd->pid;
		symbols.push_back(ps);
	}
	static void readDumpedMLNFromFile(LvrMLN& mln,string mlndumpfile = MLNDUMPFILENAME)
	{
		fstream filestr(mlndumpfile.c_str());
		if(filestr == NULL)
		{
			cout<<"Error reading LvrMLN File"<<endl;
			exit(-1);
		}
		//char* buf = new char[1024];
		while(filestr)
		{
			//filestr.getline(buf,1024);
			string line;
			getline(filestr,line);
			if(line.size()==0)
				continue;
			bool sat = false;
			//check for empty clauses
			{
				if(line[2] == '0')
				{
					if(line[0]=='1')
					{
						sat=true;
					}
					string str = line.substr(6);
					stringstream sstr(str);
					double wt;
					sstr >> wt;
					WClause* empClause = new WClause();
					empClause->satisfied = sat;
					empClause->weight = LogDouble(wt,false);
					mln.clauses.push_back(empClause);
					continue;
				}
			}

			vector<string> tokens;
			LStringConversionUtils::tokenize(line,tokens,"@");
			int pos = tokens[0].find(":");
			int sVal;
			stringstream str(tokens[0].substr(0,pos));
			str>>sVal;
			if(sVal==1)
				sat=true;
			//string remaining = tokens[0].substr(0,pos+1);
			//string nAtomsS = tokens[0].substr(0,pos);
			int p1 = tokens[0].rfind(":");

			string atomsS = tokens[0].substr(p1+1);
			vector<string> atomNames;
			LStringConversionUtils::tokenize(atomsS,atomNames,"|");
			vector<LSymbolDump*> sList;
			vector<PredicateSymbol*> symbols;
			for(unsigned int i=0;i<atomNames.size();i++)
			{
				LSymbolDump* ns = new LSymbolDump();
				convertToDump(atomNames[i],ns,symbols);
				sList.push_back(ns);
			}

			vector<string> toksTerms;
			LStringConversionUtils::tokenize(tokens[1],toksTerms,"%");
			vector<vector<LvrTerm*> > atomTerms(sList.size());
			for(unsigned int i=0;i<atomTerms.size();i++)
			{
				atomTerms[i].resize(sList[i]->nTerms);
			}
			for(unsigned int i=0;i<toksTerms.size();i++)
			{
				int pos = toksTerms[i].find("&");
				string s1 = toksTerms[i].substr(0,pos);
				string s2 = toksTerms[i].substr(pos+1);
				vector<string> dValsS;
				vector<int> dVals;
				LStringConversionUtils::tokenize(s1,dValsS,",");
				LStringConversionUtils::toIntArr(dValsS,dVals);
				LvrTerm* term = new LvrTerm(0,dVals);

				vector<string> postoks;
				LStringConversionUtils::tokenize(s2,postoks,",");
				for(unsigned int j=0;j<postoks.size();j++)
				{
					int pos = postoks[j].find(":");
					string s1 = postoks[j].substr(0,pos);
					string s2 = postoks[j].substr(pos+1);
					int apos = LStringConversionUtils::toInt(s1);
					int tpos = LStringConversionUtils::toInt(s2);
					atomTerms[apos].at(tpos) = term;
				}
			}
			WClause* clause = new WClause();
			//create the clause
			for(unsigned int i=0;i<atomTerms.size();i++)
			{
				Atom* atom = new Atom(symbols[i],atomTerms[i]);
				clause->atoms.push_back(atom);
				if(sList[i]->sign == 0)
					clause->sign.push_back(true);
				else
					clause->sign.push_back(false);
			}
			clause->weight = LogDouble(LStringConversionUtils::toDouble(tokens[2]),false);
			clause->satisfied=sat;
			mln.clauses.push_back(clause);
		}
		filestr.close();
		readDumpedSymbolsFromFile(mln);
	}

	static void readDumpedSymbolsFromFile(LvrMLN& mln,string symbolsdumpfile=SYMBOLSDUMPFILENAME)
	{
		fstream filestr(symbolsdumpfile.c_str());
		if(filestr == NULL)
		{
			cout<<"Error reading LvrMLN File"<<endl;
			exit(-1);
		}
		char* buf = new char[1024];
		filestr.getline(buf,1024);
		string line0(buf);
		mln.max_predicate_id = LStringConversionUtils::toInt(line0);
		filestr.getline(buf,1024);
		string line1(buf);
		mln.maxDegree = LStringConversionUtils::toInt(line1);
		while(filestr)
		{
			filestr.getline(buf,1024);
			string line(buf);
			if(line.size()==0)
				continue;
		
			vector<LSymbolDump*> sList;
			vector<PredicateSymbol*> symbols;
			LSymbolDump* ns = new LSymbolDump();
			convertToDump(line,ns,mln.symbols);
			sList.push_back(ns);
		}
		filestr.close();
	}

	static void dumpProposalAtom(Atom* atom, ofstream* out)
	{
		(*out)<<atom->symbol->symbol.c_str()<<"#"<<atom->symbol->id<<"#"<<atom->symbol->parentId<<"#";
		for(unsigned int i=0;i<atom->terms.size();i++)
		{
			for(unsigned int j=0;j<atom->terms[i]->domain.size();j++)
			{
				(*out)<<atom->terms[i]->domain[j];
				if(j<atom->terms[i]->domain.size()-1)
					(*out)<<",";
			}
			if(i < atom->terms.size()-1)
				(*out)<<"!";
		}
		(*out)<<"$";
	}

	static void dumpProposalToFile(vector<Atom*> atoms,vector<vector<Atom*> > parents, vector<int> ids, string proposaldumpfile=PROPOSALDUMPFILENAME)
	{
		ofstream* out = new ofstream(proposaldumpfile.c_str());
		for(unsigned int a=0;a<atoms.size();a++)
		{
			(*out)<<ids[a]<<"$";
			Atom* atom = atoms[a];
			dumpProposalAtom(atom,out);
			for(unsigned int i=0;i<parents[a].size();i++)
			{
				dumpProposalAtom(parents[a].at(i),out);
			}
			(*out)<<endl;
		}
		(*out).close();
		delete out;

	}

	static void readDumpedProposalAtom(vector<Atom*>& atoms,string atomString)
	{
		vector<string> tokens;
		LStringConversionUtils::tokenize(atomString,tokens,"#");
		string symbol = tokens[0];
		int id = LStringConversionUtils::toInt(tokens[1]);
		int parentid = LStringConversionUtils::toInt(tokens[2]);
		vector<string> termsS;
		LStringConversionUtils::tokenize(tokens[3],termsS,"!");
		vector<LvrTerm*> termDoms(termsS.size());
		for(unsigned int i=0;i<termsS.size();i++)
		{
			vector<int> vals;
			vector<string> domS;
			LStringConversionUtils::tokenize(termsS[i],domS,",");
			LStringConversionUtils::toIntArr(domS,vals);
			LvrTerm* term = new LvrTerm(0,vals);
			termDoms[i] = term;

		}
		vector<int> var_types(termDoms.size());
		PredicateSymbol *ps = new PredicateSymbol(id,symbol,var_types,LogDouble(1,false),LogDouble(1,false));
		ps->parentId = parentid;
		Atom* atom = new Atom(ps,termDoms);
		atoms.push_back(atom);
	}

	static void readDumpedProposalFile(vector<Atom*>& atoms,vector<vector<Atom*> >& parents, vector<int>& ids,string proposaldumpfile=PROPOSALDUMPFILENAME)
	{
		fstream filestr(proposaldumpfile.c_str());
		if(filestr == NULL)
		{
			cout<<"Error reading LvrMLN File"<<endl;
			exit(-1);
		}
		char* buf = new char[1024];
		while(filestr)
		{
			filestr.getline(buf,1024);
			string line(buf);
			if(line.size()==0)
				continue;
			vector<string> tokens;
			LStringConversionUtils::tokenize(line,tokens,"$");
			int id = LStringConversionUtils::toInt(tokens[0]);
			readDumpedProposalAtom(atoms,tokens[1]);
			vector<Atom*> tmpParents;
			for(unsigned int i=2;i<tokens.size();i++)
			{
				if(tokens[i].size() > 0)
				{
					readDumpedProposalAtom(tmpParents,tokens[i]);
				}
			}
			ids.push_back(id);
			parents.push_back(tmpParents);
		}
		filestr.close();
	}

	static void dumpProposalITTermsFile(vector<vector<bool> > isolatedTerms,string proposalITTermsDum = PROPOSALISOLATEDTERMSDUMPFILE)
	{
		ofstream out(proposalITTermsDum.c_str());
		for(unsigned int i=0;i<isolatedTerms.size();i++)
		{
			for(unsigned int j=0;j<isolatedTerms[i].size();j++)
			{
				if(isolatedTerms[i].at(j))
					out<<1;
				else
					out<<0;
				if(j<isolatedTerms[i].size()-1)
					out<<",";
			}
			out<<endl;
		}
		out.close();
	}

	static void readDumpedProposalIsolatedTermsFile(vector<vector<bool> >& isolatedTerms,string proposaldumpfile=PROPOSALISOLATEDTERMSDUMPFILE)
	{
		fstream filestr(proposaldumpfile.c_str());
		if(filestr == NULL)
		{
			cout<<"Error reading LvrMLN File"<<endl;
			exit(-1);
		}
		char* buf = new char[1024];
		while(filestr)
		{
			filestr.getline(buf,1024);
			string line(buf);
			if(line.size()==0)
				continue;
			vector<bool> tmpisolatedTerms;
			vector<string> tokens;
			LStringConversionUtils::tokenize(line,tokens,",");
			for(unsigned int i=0;i<tokens.size();i++)
			{
				int val;
				stringstream st(tokens[i]);
				st>>val;
				if(val > 0)
				{
					tmpisolatedTerms.push_back(true);
				}
				else
					tmpisolatedTerms.push_back(false);
			}
			isolatedTerms.push_back(tmpisolatedTerms);
		}
		filestr.close();
	}

};
#endif
