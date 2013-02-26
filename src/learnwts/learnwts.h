/*
 * All of the documentation and software included in the
 * Alchemy Software is copyrighted by Stanley Kok, Parag
 * Singla, Matthew Richardson, Pedro Domingos, Marc
 * Sumner, Hoifung Poon, Daniel Lowd, and Jue Wang.
 * 
 * Copyright [2004-09] Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang. All rights reserved.
 * 
 * Contact: Pedro Domingos, University of Washington
 * (pedrod@cs.washington.edu).
 * 
 * Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that
 * the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the
 * following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the
 * above copyright notice, this list of conditions and the
 * following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 * 
 * 3. All advertising materials mentioning features or use
 * of this software must display the following
 * acknowledgment: "This product includes software
 * developed by Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang in the Department of
 * Computer Science and Engineering at the University of
 * Washington".
 * 
 * 4. Your publications acknowledge the use or
 * contribution made by the Software to your research
 * using the following citation(s): 
 * Stanley Kok, Parag Singla, Matthew Richardson and
 * Pedro Domingos (2005). "The Alchemy System for
 * Statistical Relational AI", Technical Report,
 * Department of Computer Science and Engineering,
 * University of Washington, Seattle, WA.
 * http://alchemy.cs.washington.edu.
 * 
 * 5. Neither the name of the University of Washington nor
 * the names of its contributors may be used to endorse or
 * promote products derived from this software without
 * specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF WASHINGTON
 * AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OF WASHINGTON OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
#ifndef LEARNWTS_H_NOV_23_2005
#define LEARNWTS_H_NOV_23_2005

#include <sys/time.h>
#include "util.h"
#include "timer.h"
#include "fol.h"
#include "mln.h"
#include "indextranslator.h"


extern const char* ZZ_TMP_FILE_POSTFIX; //defined in folhelper.h
const bool DOMAINS_SHARE_DATA_STRUCT = true;

/**
 * Extracts file names from a list of files separated by a comma with
 * no space.
 * 
 * @param namesStr List of files from which the names are extracted
 * @param namesArray Array which contains the file names after extraction
 */
void extractFileNames(const char* const & namesStr, Array<string>& namesArray)
{
  if (namesStr == NULL) return;
  string s(namesStr);
  s = Util::trim(s);
  if (s.length() == 0) return;
  s.append(",");
  string::size_type cur = 0;
  string::size_type comma;
  string name;
  while (true)
  {
    comma = s.find(",", cur);
    if (comma == string::npos) return;
    name = s.substr(cur, comma-cur);
    namesArray.append(name);
    cur = comma+1;
  }
}

/**
 * Extracts a string of arguments separated by white space.
 * 
 * @param argsStr List of arguments being extracted
 * @param argc Number of arguments extracted
 * @param argv Array which contains the arguments after extraction
 */
void extractArgs(const char* const & argsStr, int& argc, char** argv)
{
  argc = 0;
  if (argsStr == NULL) return;
  string s(argsStr);
  s = Util::trim(s);
  if (s.length() == 0) return;
  s.append(" ");
  string::size_type cur = 0;
  string::size_type blank;
  string arg;

  while (true)
  {
    blank = s.find(" ", cur);
    if (blank == string::npos) return;
    arg = s.substr(cur, blank-cur);
    arg = Util::trim(arg);
    memset(argv[argc], '\0', 500);
    arg.copy(argv[argc], arg.length());
    argc++;
    cur = blank + 1;
  }
}

void createDomainAndMLN(Array<Domain*>& domains, Array<MLN*>& mlns,
                        const string& inMLNFile, ostringstream& constFiles,
                        ostringstream& dbFiles,
                        const StringHashArray* const & nonEvidPredNames,
                        const bool& addUnitClauses, const double& priorMean,
                        const bool& checkPlusTypes, const bool& mwsLazy,
                        const bool& allPredsExceptQueriesAreCW,
                        const StringHashArray* const & owPredNames,
                        const StringHashArray* const & cwPredNames)
{
  string::size_type bslash = inMLNFile.rfind("/");
  string tmp = (bslash == string::npos) ? 
               inMLNFile:inMLNFile.substr(bslash+1,inMLNFile.length()-bslash-1);
  char tmpInMLN[100];  
  sprintf(tmpInMLN, "%s%d%s",  tmp.c_str(), getpid(), ZZ_TMP_FILE_POSTFIX);

  ofstream out(tmpInMLN);
  ifstream in(inMLNFile.c_str());
  if (!out.good()) { cout<<"ERROR: failed to open "<<tmpInMLN <<endl; exit(-1);}
  if (!in.good())  { cout<<"ERROR: failed to open "<<inMLNFile<<endl; exit(-1);}

  string buffer;
  while(getline(in, buffer)) out << buffer << endl;
  in.close();

  out << constFiles.str() << endl 
      << dbFiles.str() << endl;
  out.close();
  
  // read the formulas from the input MLN
  Domain* domain = new Domain;
  MLN* mln = new MLN();
  
    // Unknown evidence atoms are filled in by EM
  bool warnAboutDupGndPreds = true;
  bool mustHaveWtOrFullStop = false;
  bool flipWtsOfFlippedClause = false;
  Domain* domain0 = (checkPlusTypes) ? domains[0] : NULL;

  bool ok = runYYParser(mln, domain, tmpInMLN, allPredsExceptQueriesAreCW, 
                        owPredNames, cwPredNames, nonEvidPredNames,
                        addUnitClauses,  warnAboutDupGndPreds, priorMean,
                        mustHaveWtOrFullStop, domain0, mwsLazy,
                        flipWtsOfFlippedClause);

  unlink(tmpInMLN);
  if (!ok) exit(-1);
  domains.append(domain);
  mlns.append(mln);
}


void createDomainsAndMLNs(Array<Domain*>& domains, Array<MLN*>& mlns, 
                          const bool& multipleDatabases,
                          const string& inMLNFile,
                          const Array<string>& constFilesArr,
                          const Array<string>& dbFilesArr,
                          const StringHashArray* const & nonEvidPredNames,
                          const bool& addUnitClauses, const double& priorMean,
                          const bool& mwsLazy,
                          const bool& allPredsExceptQueriesAreCW,
                          const StringHashArray* const & owPredNames,
                          const StringHashArray* const & cwPredNames)
{
  if (!multipleDatabases)
  {
    ostringstream constFilesStream, dbFilesStream;
    for (int i = 0; i < constFilesArr.size(); i++) 
      constFilesStream << "#include \"" << constFilesArr[i] << "\"" << endl;
    for (int i = 0; i < dbFilesArr.size(); i++)    
      dbFilesStream << "#include \"" << dbFilesArr[i] << "\"" << endl;
    createDomainAndMLN(domains, mlns, inMLNFile, constFilesStream, 
                       dbFilesStream, nonEvidPredNames,
                       addUnitClauses, priorMean, false, mwsLazy,
                       allPredsExceptQueriesAreCW, owPredNames, cwPredNames);
  }
  else
  {   //if multiple databases
    for (int i = 0; i < dbFilesArr.size(); i++) // for each domain
    {
      cout << "parsing MLN and creating domain " << i << "..." << endl;
      ostringstream constFilesStream, dbFilesStream;
      if (constFilesArr.size() > 0)
        constFilesStream << "#include \"" << constFilesArr[i] << "\"" << endl;
      dbFilesStream    << "#include \"" << dbFilesArr[i]    << "\"" << endl;
      
      bool checkPlusTypes = (i > 0);

      createDomainAndMLN(domains, mlns, inMLNFile, constFilesStream,
                         dbFilesStream, nonEvidPredNames,
                         addUnitClauses, priorMean, checkPlusTypes, mwsLazy,
                         allPredsExceptQueriesAreCW, owPredNames, cwPredNames);

        // let the domains share data structures
      if (DOMAINS_SHARE_DATA_STRUCT && i > 0)
      {
		  // Add new constants to base domain as external
          // Assumption: all domains have same number and ordering of types
        int numTypes = domains[i]->getNumTypes();
        for (int j = 0; j < numTypes; j++)
        {
          const char* typeName = domains[i]->getTypeName(j);
          const Array<int>* constantsByType = domains[i]->getConstantsByType(j);
          for (int k = 0; k < constantsByType->size(); k++)
          {
            const char* constantName =
              domains[i]->getConstantName((*constantsByType)[k]);
              // Add external constant to base domain
            if (!domains[0]->isConstant(constantName))
            {
              domains[0]->addExternalConstant(constantName, typeName);
            }
          }
        }
 
          // Share all predicate templates
        const ClauseHashArray* carr = mlns[i]->getClauses();
        for (int j = 0; j < carr->size(); j++)
        {
          Clause* c = (*carr)[j];

            // Share all predicate templates
          for (int k = 0; k < c->getNumPredicates(); k++)
          {
            Predicate* p = c->getPredicate(k);
            const PredicateTemplate* t 
              = domains[0]->getPredicateTemplate(p->getName());
            assert(t);
            p->setTemplate((PredicateTemplate*)t);
          }
        } // for all clauses

          // All mlns can use the following data structures from mln0
        ((Domain*)domains[i])->replaceTypeDualMap((
                                         DualMap*)domains[0]->getTypeDualMap());
/*
        ((Domain*)domains[i])->replaceStrToPredTemplateMapAndPredDualMap(
                  (StrToPredTemplateMap*) domains[0]->getStrToPredTemplateMap(),
                  (DualMap*) domains[0]->getPredDualMap());
*/
        ((Domain*)domains[i])->replaceStrToFuncTemplateMapAndFuncDualMap(
                  (StrToFuncTemplateMap*) domains[0]->getStrToFuncTemplateMap(),
                  (DualMap*) domains[0]->getFuncDualMap());
        ((Domain*)domains[i])->replaceEqualPredTemplate(
                        (PredicateTemplate*)domains[0]->getEqualPredTemplate());
        ((Domain*)domains[i])->replaceFuncSet(
                                        (FunctionSet*)domains[0]->getFuncSet());

      }
    } // for each domain

	// Reorder constants in all domains
	domains[0]->reorderConstants(mlns[0]);
    for (int i = 1; i < domains.size(); i++)
    {
      ((Domain*)domains[i])->reorderConstants(
                                   (ConstDualMap*)domains[0]->getConstDualMap(),
                          (Array<Array<int>*>*)domains[0]->getConstantsByType(),
                  (Array<Array<int>*>*)domains[0]->getExternalConstantsByType(),
                                              mlns[i]);
    }

      // Share clauses across all mlns as external clauses.
      // This happens when per-constant is used and dbs have diff. constants
    Array<Array<bool>*>* externalClausesPerMLN = new Array<Array<bool>*>;
    externalClausesPerMLN->growToSize(mlns.size());
    (*externalClausesPerMLN)[0] = new Array<bool>;
    (*externalClausesPerMLN)[0]->growToSize(mlns[0]->getNumClauses(), false);
    for (int i = 1; i < mlns.size(); i++)
    {
      const ClauseHashArray* carr = mlns[i]->getClauses();
      (*externalClausesPerMLN)[i] = new Array<bool>;
      for (int j = 0; j < carr->size(); j++)
      {
        Clause* c = (*carr)[j];
        if (!mlns[0]->containsClause(c))
        {
          string formulaString = mlns[i]->getParentFormula(j, 0);
          bool hasExist = mlns[i]->isExistClause(j);
          if (mlns[0]->appendExternalClause(formulaString, hasExist,
                                            new Clause(*c), domains[0], false,
                                            false, false))
          {
            (*externalClausesPerMLN)[0]->append(false);
          }
        }
      } // for all clauses
    } // for all mlns except base

    const ClauseHashArray* carr = mlns[0]->getClauses();
    for (int j = 0; j < carr->size(); j++)
    {
      Clause* c = (*carr)[j];
        // Do not add existentially qualified clauses
      //if (mlns[i]->isExistClause(j)) continue;
      for (int k = 1; k < mlns.size(); k++)
      {
        if (mlns[k]->containsClause(c))
          (*externalClausesPerMLN)[k]->append(false);
        else
          (*externalClausesPerMLN)[k]->append(true);
      }
    } // for all clauses

    for (int i = 1; i < mlns.size(); i++)
    {
      mlns[i]->replaceClauses(new
        ClauseHashArray(*(ClauseHashArray*)mlns[0]->getClauses()));
      mlns[i]->replaceMLNClauseInfos(new
        Array<MLNClauseInfo*>(*(Array<MLNClauseInfo*>*)mlns[0]->
          getMLNClauseInfos()));
      mlns[i]->replacePredIdToClausesMap(new
        Array<Array<IndexClause*>*>(*(Array<Array<IndexClause*>*>*)mlns[0]->
          getPredIdToClausesMap()));
      mlns[i]->replaceFormulaAndClausesArray(new
        FormulaAndClausesArray(*(FormulaAndClausesArray*)mlns[0]->
          getFormulaAndClausesArray()));
      mlns[i]->replaceExternalClause((*externalClausesPerMLN)[i]);
    }

    delete externalClausesPerMLN;

    for (int i = 1; i < mlns.size(); i++)
      mlns[i]->rehashClauses();
 } // if multiple databases
  

  //commented out: not true when there are domains with different constants
  //int numClauses = mlns[0]->getNumClauses();
  //for (int i = 1; i < mlns.size(); i++) 
  //  assert(mlns[i]->getNumClauses() == numClauses);
  //numClauses = 0; //avoid compilation warning
}

void assignWtsAndOutputMLN(ostream& out, Array<MLN*>& mlns, 
                           Array<Domain*>& domains, const Array<double>& wts, 
                           IndexTranslator* const& indexTrans)
{
    //assign the optimal weights belonging to clauses (and none of those 
    //belonging to existentially quantified formulas) to the MLNs
  double* wwts = (double*) wts.getItems();
  indexTrans->assignNonTiedClauseWtsToMLNs(++wwts);


    // output the predicate declaration
  out << "//predicate declarations" << endl;
  domains[0]->printPredicateTemplates(out);
  out << endl;

    // output the function declarations
  if (domains[0]->getNumFunctions() > 0) 
  {
    out << "//function declarations" << endl;
    domains[0]->printFunctionTemplates(out);
    out << endl;
  }

  mlns[0]->printMLNNonExistFormulas(out, domains[0]);

  const ClauseHashArray* clauseOrdering = indexTrans->getClauseOrdering();
  const StringHashArray* exFormOrdering = indexTrans->getExistFormulaOrdering();
  for (int i = 0; i < exFormOrdering->size(); i++)
  {
      // output the original formula and its weight
    out.width(0); out << "// "; out.width(6); 
    out << wts[1+clauseOrdering->size()+i] << "  " <<(*exFormOrdering)[i]<<endl;
    out << wts[1+clauseOrdering->size()+i] << "  " <<(*exFormOrdering)[i]<<endl;
    out << endl;
  }
}


void assignWtsAndOutputMLN(ostream& out, Array<MLN*>& mlns, 
                           Array<Domain*>& domains, const Array<double>& wts)
{
    // assign the optimal weights to the clauses in all MLNs
  for (int i = 0; i < mlns.size(); i++)
  {
    MLN* mln = mlns[i];
    const ClauseHashArray* clauses = mln->getClauses();
    for (int j = 0; j < clauses->size(); j++)
      (*clauses)[j]->setWt(wts[j+1]);
  }

  for (int i = 0; i < mlns.size(); i++)
  {
    MLN* mln = mlns[i];
    Domain* domain = domains[i];
    const Array<Clause*>* hclauses = mln->getHybridClauses();
    for (int j = 0; j < hclauses->size(); j++)
    {
      double variance = mln->getHybridClauseVariance(j, domain);
      double w = 1.0/(2.0*variance);
      (*hclauses)[j]->setWt(w);
    }
  }

    // output the predicate declaration
  out << "//predicate declarations" << endl;
  domains[0]->printPredicateTemplates(out);
  out << endl;

    // output the function declarations
  if (domains[0]->getNumFunctions() > 0) 
  {
    // output the function declarations
    out << "//function declarations" << endl;
    domains[0]->printFunctionTemplates(out);
    out << endl;
  }
  mlns[0]->printMLN(out, domains[0]);
}


void deleteDomains(Array<Domain*>& domains)
{
  for (int i = 0; i < domains.size(); i++) 
  {
    if (DOMAINS_SHARE_DATA_STRUCT && i > 0)
    {
      ((Domain*)domains[i])->setTypeDualMap(NULL);
      ((Domain*)domains[i])->setStrToPredTemplateMapAndPredDualMap(NULL, NULL);
      ((Domain*)domains[i])->setStrToFuncTemplateMapAndFuncDualMap(NULL, NULL);
      ((Domain*)domains[i])->setEqualPredTemplate(NULL);
      ((Domain*)domains[i])->setFuncSet(NULL);

        // Need this since it is pointing to domain0's copy
	  ((Domain*)domains[i])->setConstDualMap(NULL);
	  ((Domain*)domains[i])->setConstantsByType(NULL);
    }
    delete domains[i];
  }
}


#endif
