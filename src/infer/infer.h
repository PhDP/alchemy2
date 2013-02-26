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
#ifndef _INFER_H_OCT_30_2005
#define _INFER_H_OCT_30_2005

/**
 * Collection of functions useful for performing inference.
 */

#include "util.h"
#include "fol.h"
#include "mrf.h"
#include "learnwts.h"
#include "inferenceargs.h"
#include "maxwalksat.h"
#include "mcsat.h"
#include "gibbssampler.h"
#include "simulatedtempering.h"
#include "bp.h"
#include "variablestate.h"
#include "hvariablestate.h"
#include "hmcsat.h"
#include "lbfgsp.h"

// Variables for holding inference command line args are in inferenceargs.h

char* aevidenceFiles  = NULL;
char* aresultsFile    = NULL;
char* aqueryPredsStr  = NULL;
char* aqueryFile      = NULL;

int anumerator = 50;
int adenominator = 100;

char* aMaxSeconds = NULL;
bool aStartPt = true;

int saInterval = 100;

bool aPrintSamplePerIteration = false;

string queryPredsStr, queryFile;
GroundPredicateHashArray queries;
GroundPredicateHashArray knownQueries;

/**
 * Creates query predicates from a string of predicates
 * separated by a comma without a space. This method also
 * serves to print out ground preds and their probs / truth
 * values to a stream if printToFile is set to true.
 * 
 * @param queryPredsStr Comma-separated string of predicate names which is
 * parsed.
 * @param domain Domain from which predicates are generated.
 * @param db Database containing truth values of the preds.
 * @param queries Array being filled with groundings of preds specified in
 * queryPredsStr (those with unknown values in db).
 * @param knownQueries Array being filled with groundings of preds specified in
 * queryPredsStr (those with known values in db).
 * @param allPredGndingsAreQueries Boolean array indicating which predicates
 * have all their gndings as queries (indexed by the predicate's id in domain)
 * @param printToFile If true, atom, and their truth values, are printed to
 * file, not collected in HashArray.
 * @param out stream to which the groundings are printed when printToFile is set
 * to true.
 * @param amapPos If true, only positive atoms are printed out. Otherwise, all
 * atoms with 0/1 are printed. Used in conjunction with printToFile.
 * @param trueQueries If specified, probabilities from the ground preds here are
 * written to out. Used in conjunction with printToFile.
 * 
 * @return false, if a pred name or constant is not defined, or if a pred has
 * the wrong number of terms; otherwise, true. 
 */
bool createComLineQueryPreds(const string& queryPredsStr,
                             const Domain* const & domain,
                             Database* const & db,
                             GroundPredicateHashArray* const & queries,
                             GroundPredicateHashArray* const & knownQueries,
                             Array<int>* const & allPredGndingsAreQueries,
                             bool printToFile, ostream& out, bool amapPos,
                            const GroundPredicateHashArray* const & trueQueries,
                             const Array<double>* const & trueProbs,
                             Array<Array<Predicate* >* >* queryConjs)
{
  if (queryPredsStr.length() == 0) return true;
  string predConjs = Util::trim(queryPredsStr);

    //replace the comma or semi-colon between query predicates with '\n'
  int balparen = 0;
  for (unsigned int i = 0; i < predConjs.length(); i++)
  {
    if (predConjs.at(i)=='(')                     balparen++;
    else if (predConjs.at(i)==')')                balparen--;
    else if ((predConjs.at(i)==';' || predConjs.at(i)==',') &&
             balparen==0) predConjs.at(i) = '\n';
  }
  
  bool ret = true;
  string predConj;
  istringstream iss(predConjs);
  char delimit[2]; delimit[1] = '\0';

    // for each query formula
  while (getline(iss, predConj))
  {
      // replace the ^ with '\n'
    int balparen = 0;
    for (unsigned int i = 0; i < predConj.length(); i++)
    {
      if (predConj.at(i)=='(')                     balparen++;
      else if (predConj.at(i)==')')                balparen--;
      else if (predConj.at(i)=='^' && balparen==0) predConj.at(i) = '\n';
    }
    bool onlyPredName;
    unsigned int cur;
    int termId, varIdCnt = 0;
    hash_map<string, int, HashString, EqualString> varToId;
    hash_map<string, int, HashString, EqualString>::iterator it;
    Array<VarsTypeId*>* vtiArr;
    string pred, predName, term;
    istringstream iss2(predConj);
    const PredicateTemplate* ptemplate;
    int predicate = 0;

    Array<Predicate* >* predArray = new Array<Predicate*>;
    if (queryConjs) queryConjs->append(predArray);
      // for each query pred on command line
    while (getline(iss2, pred))
    {
      pred = Util::trim(pred);
      onlyPredName = false;
      varToId.clear();
      varIdCnt = 0;
      cur = 0;
      bool negated = false;

        // find if pred is negated
      if (pred.at(0) == '!')
      {
        negated = true;
        pred.at(0) = ' ';
        pred = Util::trim(pred);
      }
        // get predicate name
      if (!Util::substr(pred, cur, predName, "("))
      {
        predName = pred;
        onlyPredName = true;
      }
    
        // Predicate must be in the domain
      ptemplate = domain->getPredicateTemplate(predName.c_str());
      if (ptemplate == NULL)
      {
        cout << "ERROR: Cannot find command line query predicate" << predName 
             << " in domain." << endl;
        ret = false;
        continue;
      }
      Predicate ppred(ptemplate);

        // if the terms of the query predicate are also specified
      if (!onlyPredName)
      {
          // get term name
        for (int i = 0; i < 2; i++)
        {
          if (i == 0) delimit[0] = ',';
          else        delimit[0] = ')';
          while (Util::substr(pred, cur, term, delimit))
          {
              // this is a constant
            if (isupper(term.at(0)) || term.at(0) == '"' || isdigit(term.at(0)))
            {
              termId = domain->getConstantId(term.c_str());
              if (termId < 0) 
              {
                cout <<"ERROR: Cannot find constant "<<term<<" in database"<<endl;
                ret = false;
              }
            }
            else
            {   // it is a variable        
              if ((it=varToId.find(term)) == varToId.end()) 
              {
                termId = --varIdCnt;
                varToId[term] = varIdCnt; 
              }
              else
                termId = (*it).second;
            }
            ppred.appendTerm(new Term(termId, (void*)&ppred, true));
          }
        }
      }
      else
      {   // if only the predicate name is specified
          // HACK DEBUG
        //(*allPredGndingsAreQueries)[ptemplate->getId()] = true;
        for (int i = 0; i < ptemplate->getNumTerms(); i++)
          ppred.appendTerm(new Term(--varIdCnt, (void*)&ppred, true));
      }  

        // Check if number of terms is correct
      if (ppred.getNumTerms() != ptemplate->getNumTerms())
      {
        cout << "ERROR: " << predName << " requires " << ptemplate->getNumTerms()
             << " terms but given " << ppred.getNumTerms() << endl;
        ret = false;
      }
      if (!ret) continue;
    
      ///////////////////// create all groundings of predicate ///////////////
      vtiArr = NULL;
      ppred.createVarsTypeIdArr(vtiArr);

        // If a ground predicate was specified on command line
      if (vtiArr->size() <= 1)
      {
        assert(ppred.isGrounded());
        assert(!db->isClosedWorld(ppred.getId()));
        TruthValue tv = db->getValue(&ppred);

        if (negated) ppred.setSense(false);
        predArray->append(new Predicate(ppred));
        
        GroundPredicate* gndPred = new GroundPredicate(&ppred);
          // If just printing to file, then all values must be known
        if (printToFile) assert(tv != UNKNOWN);
        if (tv == UNKNOWN)
        {
          if (queries->append(gndPred) < 0) delete gndPred;
        }
        else
        {
            // If just printing to file
          if (printToFile)
          {
              // If trueQueries is given as argument, then get prob. from there
            if (trueQueries)
            {
              double prob = 0.0;
              if (domain->getDB()->getEvidenceStatus(&ppred))
              {
                  // Don't print out evidence atoms
                continue;
                //prob = (tv == TRUE) ? 1.0 : 0.0;
              }
              else
              {
                int found = trueQueries->find(gndPred);
                if (found >= 0) prob = (*trueProbs)[found];
                else
                    // Uniform smoothing
                  prob = (prob*10000+1/2.0)/(10000+1.0);
              
              }
              gndPred->print(out, domain); out << " " << prob << endl;
            }
            else
            {
              if (amapPos) //if show postive ground query predicates only
              {
      		    if (tv == TRUE)
                {
      	  	      ppred.printWithStrVar(out, domain);
      	  	      out << endl;
                }
              }
              else //print all ground query predicates
              {
                ppred.printWithStrVar(out, domain);
                out << " " << tv << endl;
              }
            }
            delete gndPred;
          }
          else // Building queries for HashArray
          {
            //if (tv == TRUE) gndPred->setProbTrue(1);
            //else            gndPred->setProbTrue(0);
            if (knownQueries->append(gndPred) < 0) delete gndPred;  
          }
        }
      }
      else // Variables need to be grounded
      {
        ArraysAccessor<int> acc;
        for (int i = 1; i < vtiArr->size(); i++)
        {
          const Array<int>* cons =
            domain->getConstantsByType((*vtiArr)[i]->typeId);
          acc.appendArray(cons);
        } 

          // form all groundings of the predicate
        Array<int> constIds;
        while (acc.getNextCombination(constIds))
        {
          assert(constIds.size() == vtiArr->size()-1);
          for (int j = 0; j < constIds.size(); j++)
          {
            Array<Term*>& terms = (*vtiArr)[j+1]->vars;
            for (int k = 0; k < terms.size(); k++)
              terms[k]->setId(constIds[j]);
          }

          // at this point the predicate is grounded
          assert(!db->isClosedWorld(ppred.getId()));
 
          TruthValue tv = db->getValue(&ppred);        
          GroundPredicate* gndPred = new GroundPredicate(&ppred);

            // If just printing to file, then all values must be known
          if (printToFile) assert(tv != UNKNOWN);
          if (tv == UNKNOWN)
          {
            if (queries->append(gndPred) < 0) delete gndPred;
          }
          else
          {
              // If just printing to file
            if (printToFile)
            {
                // If trueQueries is given as argument, then get prob. from there
              if (trueQueries)
              {
                double prob = 0.0;
                if (domain->getDB()->getEvidenceStatus(&ppred))
                {
                    // Don't print out evidence atoms
                  continue;
                  //prob = (tv == TRUE) ? 1.0 : 0.0;
                }
                else
                {
                  int found = trueQueries->find(gndPred);
                  if (found >= 0) prob = (*trueProbs)[found];
                  else
                      // Uniform smoothing
                    prob = (prob*10000+1/2.0)/(10000+1.0);
                }
                  // Uniform smoothing
                //prob = (prob*10000+1/2.0)/(10000+1.0);
                gndPred->print(out, domain); out << " " << prob << endl;
              }
              else
              {
                if (amapPos) //if show postive ground query predicates only
                {
                  if (tv == TRUE)
                  {
                    ppred.printWithStrVar(out, domain);
                    out << endl;
                  }
                }
                else //print all ground query predicates
                {
                  ppred.printWithStrVar(out, domain);
                  out << " " << tv << endl;
                }
              }
              delete gndPred;          
            }
            else // Building queries
            {
                //if (tv == TRUE) gndPred->setProbTrue(1);
          	    //else            gndPred->setProbTrue(0);
              if (knownQueries->append(gndPred) < 0) delete gndPred;  
            }
          }
        }
      }
      
      ppred.deleteVarsTypeIdArr(vtiArr);
      predicate++;
    } // for each query pred on command line
    if (predArray->size() == 0)
    {
      if (queryConjs) queryConjs->removeLastItem();
      delete predArray;
    }    
  } // for each query formula

  if (!printToFile)
  {
  	queries->compress();
  	knownQueries->compress();
  }
  
  return ret;
}

/**
 * Calls createComLineQueryPreds (see above) with trueQueries set to NULL
 * for backwards compatibility.
 */
bool createComLineQueryPreds(const string& queryPredsStr,
                             const Domain* const & domain,
                             Database* const & db,
                             GroundPredicateHashArray* const & queries,
                             GroundPredicateHashArray* const & knownQueries,
                             Array<int>* const & allPredGndingsAreQueries,
                             Array<Array<Predicate* >* >* queryConjs)
{
  return createComLineQueryPreds(queryPredsStr, domain, db,
                             	 queries, knownQueries,
                                 allPredGndingsAreQueries,
                             	 false, cout, false, NULL, NULL, queryConjs);
}

/**
 * Extracts predicate names specified in a string and in a file
 * and stores them in a HashArray. Note that preds is a copy.
 * 
 * @param preds string of comma-separated predicate names to be parsed
 * @param queryFile name of file containing predicate names
 * @param predNames Stores the parsed predicate names
 * 
 * @return false, if queryFile can not be opened; otherwise, true.
 */
bool extractPredNames(string preds, const string* queryFile, 
                      StringHashArray& predNames)
{ 
  predNames.clear();

    // first extract the query pred names specified on command line
  string::size_type cur = 0, ws, ltparen;
  string qpred, predName;
  
  if (preds.length() > 0)
  {
    preds.append(" "); //terminate preds with a whitespace
    
      //replace the comma or semi-colon between query predicates with ' '
    int balparen = 0;
    for (unsigned int i = 0; i < preds.length(); i++)
    {
      if (preds.at(i) == '(')      balparen++;
      else if (preds.at(i) == ')') balparen--;
      else if ((preds.at(i) == ',' || preds.at(i) == ';') &&
               balparen == 0) preds.at(i) = ' ';
    }
    
    while (preds.at(cur) == ' ') cur++;
    while (true)
    {
      ws = preds.find(" ", cur);
      if (ws == string::npos) break;
      qpred = preds.substr(cur,ws-cur+1);
      cur = ws+1;
      while (cur < preds.length() &&
             (preds.at(cur) == ' ' || preds.at(cur) == '^' ||
              preds.at(cur) == '!')) cur++;
      ltparen = qpred.find("(",0);
      
      if (ltparen == string::npos) 
      { 
        ws = qpred.find(" ");
        if (ws != string::npos) qpred = qpred.substr(0,ws);
        predName = qpred; 
      }
      else
        predName = qpred.substr(0,ltparen);
      
      predNames.append(predName);
    }
  }

  if (queryFile == NULL || queryFile->length() == 0) return true;

    // next extract query predicates specified in query file  
  ifstream in((*queryFile).c_str());
  if (!in.good())
  {
    cout << "ERROR: unable to open " << *queryFile << endl;
    return false;
  }
  string buffer;
  while (getline(in, buffer))
  {
    cur = 0;
    while (cur < buffer.length() &&
           (buffer.at(cur) == ' ' || buffer.at(cur) == '^' ||
            buffer.at(cur) == '!')) cur++;
    ltparen = buffer.find("(", cur);
    if (ltparen == string::npos) continue;
    predName = buffer.substr(cur, ltparen-cur);
    predNames.append(predName);
  }

  in.close();
  return true;
}

/**
 * Get the ids of the predicates represented by the string
 */
IntHashArray *getPredicateIds(string & predsStr, Domain * const & domain)
{
  IntHashArray *predIds = new IntHashArray();
  StringHashArray predNames;
  extractPredNames(predsStr,NULL,predNames);
  for(int i=0;i<predNames.size();i++)
  {
    string predName = predNames[i];
    int predId = domain->getPredicateId(predName.c_str());
    if(predId < 0)
      continue;
    predIds->append(predId);
  }
  return predIds;
}


/**
 * Returns the first character of the truth value
 * TRUE, FALSE or UKNOWN
 */
char getTruthValueFirstChar(const TruthValue& tv)
{
  if (tv == TRUE)    return 'T';
  if (tv == FALSE)   return 'F';
  if (tv == UNKNOWN) return 'U';
  assert(false);
  exit(-1);
  return '#'; //avoid compilation warning
}

/**
 * Sets the seed for the srand() function according to time of day.
 */
void setsrand()
{
  struct timeval tv;
  struct timezone tzp;
  gettimeofday(&tv,&tzp);
  unsigned int seed = (( tv.tv_sec & 0177 ) * 1000000) + tv.tv_usec;
  srand(seed);
} 


  //copy srcFile to dstFile, & append '#include "dbFiles"' to latter
void copyFileAndAppendDbFile(const string& srcFile, string& dstFile, 
                             const Array<string>& dbFilesArr,
                             const Array<string>& constFilesArr)
{
  ofstream out(dstFile.c_str());
  ifstream in(srcFile.c_str());
  if (!out.good()) { cout<<"ERROR: failed to open "<<dstFile<<endl;exit(-1);}
  if (!in.good()) { cout<<"ERROR: failed to open "<<srcFile<<endl;exit(-1);}
  
  string buffer;
  while(getline(in, buffer)) out << buffer << endl;
  in.close();

  out << endl;
  for (int i = 0; i < constFilesArr.size(); i++) 
    out << "#include \"" << constFilesArr[i] << "\"" << endl;
  out << endl;
  for (int i = 0; i < dbFilesArr.size(); i++)    
    out << "#include \"" << dbFilesArr[i] << "\"" << endl;
  out.close();
}


bool checkQueryPredsNotInClosedWorldPreds(const StringHashArray& qpredNames,
                                          const StringHashArray& cwPredNames)
{
  bool ok = true;
  for (int i = 0; i < qpredNames.size(); i++)
    if (cwPredNames.contains(qpredNames[i]))
    {
      cout << "ERROR: query predicate " << qpredNames[i] 
           << " cannot be specified as closed world" << endl; 
      ok = false;
    }
  return ok;
}

/**
 * Creates query predicates from a file of predicates
 * separated by newline character. This method also
 * serves to print out ground preds and their probs / truth
 * values to a stream if printToFile is set to true.
 * 
 * @param queryFile Location of file containing predicate
 * names which is parsed.
 * @param domain Domain from which predicates are generated.
 * @param db Database containing truth values of the preds.
 * @param queries Array being filled with groundings of preds
 * specified in queryPredsStr (those with unknown values in db).
 * @param knownQueries Array being filled with groundings of preds
 * specified in queryPredsStr (those with known values in db).
 * @param printToFile If true, atom, and their truth values,
 * are printed to file, not collected in HashArray.
 * @param out stream to which the groundings are printed when
 * printToFile is set to true.
 * @param amapPos If true, only positive atoms are printed out.
 * Otherwise, all atoms with 0/1 are printed. Used in
 * conjunction printToFile.
 * @param trueQueries If specified, probabilities from the ground
 * preds here are written to out. Used in conjunction with
 * printToFile.
 * 
 * @return false, if a pred name or constant is not defined, or if
 * a pred has the wrong number of terms; otherwise, true. 
 */
bool createQueryFilePreds(const string& queryFile,
                          const Domain* const & domain,
                          Database* const & db,
                          GroundPredicateHashArray* const &queries,
                          GroundPredicateHashArray* const &knownQueries,
                          bool printToFile, ostream& out, bool amapPos,
                          const GroundPredicateHashArray* const &trueQueries,
                          const Array<double>* const & trueProbs,
                          Array<Array<Predicate* >* >* queryConjs)
{
  if (queryFile.length() == 0) return true;

  bool ret = true;
  string predConj;
  ifstream in(queryFile.c_str());
  char delimit[2]; delimit[1] = '\0';

    // for each query formula
  while (getline(in, predConj))
  {
      // replace the ^ with '\n'
    int balparen = 0;
    for (unsigned int i = 0; i < predConj.length(); i++)
    {
      if (predConj.at(i)=='(')                     balparen++;
      else if (predConj.at(i)==')')                balparen--;
      else if (predConj.at(i)=='^' && balparen==0) predConj.at(i) = '\n';
    }

    bool onlyPredName;
    unsigned int cur;
    int termId, varIdCnt = 0;
    hash_map<string, int, HashString, EqualString> varToId;
    hash_map<string, int, HashString, EqualString>::iterator it;
    Array<VarsTypeId*>* vtiArr;
    string pred, predName, term;
    istringstream iss2(predConj);
    const PredicateTemplate* ptemplate;
    int predicate = 0;

    Array<Predicate* >* predArray = new Array<Predicate*>;
    queryConjs->append(predArray);
      // for each query pred on command line
    while (getline(iss2, pred))
    {
      pred = Util::trim(pred);
      onlyPredName = false;
      varToId.clear();
      varIdCnt = 0;
      cur = 0;
      bool negated = false;

        // find if pred is negated
      if (pred.at(0) == '!')
      {
        negated = true;
        pred.at(0) = ' ';
        pred = Util::trim(pred);
      }
        // get predicate name
      if (!Util::substr(pred, cur, predName, "("))
      {
        predName = pred;
        onlyPredName = true;
      }
    
        // Predicate must be in the domain
      ptemplate = domain->getPredicateTemplate(predName.c_str());
      if (ptemplate == NULL)
      {
        cout << "ERROR: Cannot find command line query predicate" << predName 
             << " in domain." << endl;
        ret = false;
        continue;
      }
      Predicate ppred(ptemplate);

        // if the terms of the query predicate are also specified
      if (!onlyPredName)
      {
          // get term name
        for (int i = 0; i < 2; i++)
        {
          if (i == 0) delimit[0] = ',';
          else        delimit[0] = ')';
          while (Util::substr(pred, cur, term, delimit))
          {
              // this is a constant
            if (isupper(term.at(0)) || term.at(0) == '"' || isdigit(term.at(0)))
            {
              termId = domain->getConstantId(term.c_str());
              if (termId < 0) 
              {
                cout <<"ERROR: Cannot find constant "<<term<<" in database"<<endl;
                ret = false;
              }
            }
            else
            {   // it is a variable        
              if ((it=varToId.find(term)) == varToId.end()) 
              {
                termId = --varIdCnt;
                varToId[term] = varIdCnt; 
              }
              else
                termId = (*it).second;
            }
            ppred.appendTerm(new Term(termId, (void*)&ppred, true));
          }
        }
      }
      else
      {   // if only the predicate name is specified
          // HACK DEBUG
        //(*allPredGndingsAreQueries)[ptemplate->getId()] = true;
        for (int i = 0; i < ptemplate->getNumTerms(); i++)
          ppred.appendTerm(new Term(--varIdCnt, (void*)&ppred, true));
      }  

        // Check if number of terms is correct
      if (ppred.getNumTerms() != ptemplate->getNumTerms())
      {
        cout << "ERROR: " << predName << " requires " << ptemplate->getNumTerms()
             << " terms but given " << ppred.getNumTerms() << endl;
        ret = false;
      }
      if (!ret) continue;
    
      ///////////////////// create all groundings of predicate ///////////////
      vtiArr = NULL;
      ppred.createVarsTypeIdArr(vtiArr);

        // If a ground predicate was specified on command line
      if (vtiArr->size() <= 1)
      {
        assert(ppred.isGrounded());
        assert(!db->isClosedWorld(ppred.getId()));
        TruthValue tv = db->getValue(&ppred);

        if (negated) ppred.setSense(false);
        predArray->append(new Predicate(ppred));
        
        GroundPredicate* gndPred = new GroundPredicate(&ppred);
          // If just printing to file, then all values must be known
        if (printToFile) assert(tv != UNKNOWN);
        if (tv == UNKNOWN)
        {
          if (queries->append(gndPred) < 0) delete gndPred;
        }
        else
        {
            // If just printing to file
          if (printToFile)
          {
              // If trueQueries is given as argument, then get prob. from there
            if (trueQueries)
            {
              double prob = 0.0;
              if (domain->getDB()->getEvidenceStatus(&ppred))
              {
                  // Don't print out evidence atoms
                continue;
                //prob = (tv == TRUE) ? 1.0 : 0.0;
              }
              else
              {
                int found = trueQueries->find(gndPred);
                if (found >= 0) prob = (*trueProbs)[found];
                else
                    // Uniform smoothing
                  prob = (prob*10000+1/2.0)/(10000+1.0);
              
              }
              gndPred->print(out, domain); out << " " << prob << endl;
            }
            else
            {
              if (amapPos) //if show postive ground query predicates only
              {
                if (tv == TRUE)
                {
                  ppred.printWithStrVar(out, domain);
                  out << endl;
                }
              }
              else //print all ground query predicates
              {
                ppred.printWithStrVar(out, domain);
                out << " " << tv << endl;
              }
            }
            delete gndPred;
          }
          else // Building queries for HashArray
          {
            //if (tv == TRUE) gndPred->setProbTrue(1);
            //else            gndPred->setProbTrue(0);
            if (knownQueries->append(gndPred) < 0) delete gndPred;  
          }
        }      
      }
      else // Variables need to be grounded
      {
        ArraysAccessor<int> acc;
        for (int i = 1; i < vtiArr->size(); i++)
        {
          const Array<int>* cons =
            domain->getConstantsByType((*vtiArr)[i]->typeId);
          acc.appendArray(cons);
        } 

          // form all groundings of the predicate
        Array<int> constIds;
        while (acc.getNextCombination(constIds))
        {
          assert(constIds.size() == vtiArr->size()-1);
          for (int j = 0; j < constIds.size(); j++)
          {
            Array<Term*>& terms = (*vtiArr)[j+1]->vars;
            for (int k = 0; k < terms.size(); k++)
              terms[k]->setId(constIds[j]);
          }

          // at this point the predicate is grounded
          assert(!db->isClosedWorld(ppred.getId()));
 
          TruthValue tv = db->getValue(&ppred);        
          GroundPredicate* gndPred = new GroundPredicate(&ppred);

            // If just printing to file, then all values must be known
          if (printToFile) assert(tv != UNKNOWN);
          if (tv == UNKNOWN)
          {
            if (queries->append(gndPred) < 0) delete gndPred;
          }
          else
          {
              // If just printing to file
            if (printToFile)
            {
                // If trueQueries is given as argument, then get prob. from there
              if (trueQueries)
              {
                double prob = 0.0;
                if (domain->getDB()->getEvidenceStatus(&ppred))
                {
                    // Don't print out evidence atoms
                  continue;
                  //prob = (tv == TRUE) ? 1.0 : 0.0;
                }
                else
                {
                  int found = trueQueries->find(gndPred);
                  if (found >= 0) prob = (*trueProbs)[found];
                  else
                      // Uniform smoothing
                    prob = (prob*10000+1/2.0)/(10000+1.0);
                }
                  // Uniform smoothing
                //prob = (prob*10000+1/2.0)/(10000+1.0);
                gndPred->print(out, domain); out << " " << prob << endl;
              }
              else
              {
                if (amapPos) //if show postive ground query predicates only
                {
                  if (tv == TRUE)
                  {
                    ppred.printWithStrVar(out, domain);
                    out << endl;
                  }
                }
                else //print all ground query predicates
                {
                  ppred.printWithStrVar(out, domain);
                  out << " " << tv << endl;
                }
              }
              delete gndPred;          
            }
            else // Building queries
            {
                //if (tv == TRUE) gndPred->setProbTrue(1);
                //else            gndPred->setProbTrue(0);
              if (knownQueries->append(gndPred) < 0) delete gndPred;  
            }
          }
        }
      }
      
      ppred.deleteVarsTypeIdArr(vtiArr);
      predicate++;
    } // for each query pred on command line
    if (predArray->size() == 0)
    {
      if (queryConjs) queryConjs->removeLastItem();
      delete predArray;
    }    
  } // fore each query formula

  if (!printToFile)
  {
    queries->compress();
    knownQueries->compress();
  }
  
  in.close();
  return ret;
}

/**
 * Calls createQueryFilePreds (see above) with trueQueries set to NULL
 * for backwards compatibility.
 */
bool createQueryFilePreds(const string& queryFile, const Domain* const & domain,
                          Database* const & db,
                          GroundPredicateHashArray* const &queries,
                          GroundPredicateHashArray* const &knownQueries,
                          Array<Array<Predicate* >* >* queryConjs)
{
  return createQueryFilePreds(queryFile, domain, db, queries, knownQueries,
                              false, cout, false, NULL, NULL, queryConjs);
}


/**
 * Parses the .mln and .db files and builds an inference procedure.
 * 
 * @param inference Inference built from the files and parameters.
 * @return -1, if an error occurred; otherwise 1
 */
int buildInference(Inference*& inference, Domain*& domain,
                   bool const &aisQueryEvidence, Array<Predicate *> &queryPreds,
                   Array<TruthValue> &queryPredValues)
{
  string inMLNFile, wkMLNFile, evidenceFile;

  StringHashArray queryPredNames;
  StringHashArray owPredNames;
  StringHashArray cwPredNames;
  MLN* mln = NULL;
  Array<string> constFilesArr;
  Array<string> evidenceFilesArr;

  
  //the second .mln file to the last one in ainMLNFiles _may_ be used 
  //to hold constants, so they are held in constFilesArr. They will be
  //included into the first .mln file.

    //extract .mln, .db file names
  extractFileNames(ainMLNFiles, constFilesArr);
  assert(constFilesArr.size() >= 1);
  inMLNFile.append(constFilesArr[0]);
  constFilesArr.removeItem(0);
  extractFileNames(aevidenceFiles, evidenceFilesArr);
  
  if (aqueryPredsStr) queryPredsStr.append(aqueryPredsStr);
  if (aqueryFile) queryFile.append(aqueryFile);

  if (queryPredsStr.length() == 0 && queryFile.length() == 0)
  { cout << "No query predicates specified" << endl; return -1; }

  if (agibbsInfer && agibbsTestConvergence && amcmcNumChains < 2) 
  {
    cout << "ERROR: If testing for convergence, there must be at least 2 "
         << "MCMC chains in Gibbs sampling" 
         << endl; return -1;
  }

  if (adecisionInfer && !abpInfer && !aefbpInfer)
  {
    aefbpInfer = true;
  }
  else if (!asimtpInfer && !amapPos && !amapAll && !agibbsInfer &&
           !amcsatInfer && !aHybrid && !aSA && !abpInfer && !aefbpInfer &&
           !aoutputNetwork)
  {
      // If nothing specified, use MC-SAT
    amcsatInfer = true;
  }

    //extract names of all query predicates
  if (queryPredsStr.length() > 0 || queryFile.length() > 0)
  {
    if (!extractPredNames(queryPredsStr, &queryFile, queryPredNames)) return -1;
  }

  if (amwsMaxSteps <= 0)
  { cout << "ERROR: mwsMaxSteps must be positive" << endl; return -1; }

  if (amwsTries <= 0)
  { cout << "ERROR: mwsTries must be positive" << endl; return -1; }

    //extract names of open-world evidence predicates
  if (aOpenWorldPredsStr)
  {
    if (!extractPredNames(string(aOpenWorldPredsStr), NULL, owPredNames)) 
      return -1;
    assert(owPredNames.size() > 0);
  }

    //extract names of closed-world non-evidence predicates
  if (aClosedWorldPredsStr)
  {
    if (!extractPredNames(string(aClosedWorldPredsStr), NULL, cwPredNames)) 
      return -1;
    assert(cwPredNames.size() > 0);
    if (!checkQueryPredsNotInClosedWorldPreds(queryPredNames, cwPredNames))
      return -1;
  }

  // TODO: Check if query atom in -o -> error

  // TODO: Check if atom in -c and -o -> error


  // TODO: Check if ev. atom in -c or
  // non-evidence in -o -> warning (this is default)


    // Set SampleSat parameters
  SampleSatParams* ssparams = new SampleSatParams;
  ssparams->lateSa = assLateSa;
  ssparams->saRatio = assSaRatio;
  ssparams->saTemp = assSaTemp;

    // Set MaxWalksat parameters
  MaxWalksatParams* mwsparams = new MaxWalksatParams;
  mwsparams->ssParams = ssparams;
  mwsparams->maxSteps = amwsMaxSteps;
  mwsparams->maxTries = amwsTries;
  mwsparams->targetCost = amwsTargetWt;
  mwsparams->hard = amwsHard;
    // numSolutions only applies when used in SampleSat.
    // When just MWS, this is set to 1
  mwsparams->numSolutions = amwsNumSolutions;
  mwsparams->heuristic = amwsHeuristic;
  mwsparams->tabuLength = amwsTabuLength;
  mwsparams->lazyLowState = amwsLazyLowState;

    // Set MC-SAT parameters
  MCSatParams* msparams = new MCSatParams;
  msparams->mwsParams = mwsparams;
    // MC-SAT needs only one chain
  msparams->numChains          = 1;
  msparams->burnMinSteps       = amcmcBurnMinSteps;
  msparams->burnMaxSteps       = amcmcBurnMaxSteps;
  msparams->minSteps           = amcmcMinSteps;
  msparams->maxSteps           = amcmcMaxSteps;
  msparams->maxSeconds         = amcmcMaxSeconds;
  //msparams->numStepsEveryMCSat = amcsatNumStepsEveryMCSat;

    // Set Gibbs parameters
  GibbsParams* gibbsparams = new GibbsParams;
  gibbsparams->mwsParams    = mwsparams;
  gibbsparams->numChains    = amcmcNumChains;
  gibbsparams->burnMinSteps = amcmcBurnMinSteps;
  gibbsparams->burnMaxSteps = amcmcBurnMaxSteps;
  gibbsparams->minSteps     = amcmcMinSteps;
  gibbsparams->maxSteps     = amcmcMaxSteps;
  gibbsparams->maxSeconds   = amcmcMaxSeconds;

  gibbsparams->gamma           = 1 - agibbsDelta;
  gibbsparams->epsilonError    = agibbsEpsilonError;
  gibbsparams->fracConverged   = agibbsFracConverged;
  gibbsparams->walksatType     = agibbsWalksatType;
  gibbsparams->testConvergence = agibbsTestConvergence;
  gibbsparams->samplesPerTest  = agibbsSamplesPerTest;
  
    // Set Sim. Tempering parameters
  SimulatedTemperingParams* stparams = new SimulatedTemperingParams;
  stparams->mwsParams    = mwsparams;
  stparams->numChains    = amcmcNumChains;
  stparams->burnMinSteps = amcmcBurnMinSteps;
  stparams->burnMaxSteps = amcmcBurnMaxSteps;
  stparams->minSteps     = amcmcMinSteps;
  stparams->maxSteps     = amcmcMaxSteps;
  stparams->maxSeconds   = amcmcMaxSeconds;

  stparams->subInterval = asimtpSubInterval;
  stparams->numST       = asimtpNumST;
  stparams->numSwap     = asimtpNumSwap;

    // Set BP parameters
  BPParams* bpparams = new BPParams;
  bpparams->maxSteps               = amcmcMaxSteps;
  bpparams->maxSeconds             = amcmcMaxSeconds;
  bpparams->lifted                 = aliftedInfer;
  bpparams->useHC                  = auseHC;
  bpparams->useCT                  = auseCT;
  bpparams->convergenceThresh      = abpConvergenceThresh;
  bpparams->convergeRequiredItrCnt = abpConvergeRequiredItrCnt;
  bpparams->implicitRep            = !aexplicitRep;
  bpparams->hcCreateType           = (HyperCubeCreateType)ahcCreateType;
  bpparams->hcCreateNoise          = ahcCreateNoise;
  bpparams->lncIter                = alncIter;
  bpparams->outputNetwork          = aoutputNetwork;  
  if (anoHCPredsStr)
  {
  	(bpparams->noHCPredsStr).append(anoHCPredsStr);
  }

  //////////////////// read in clauses & evidence predicates //////////////////

  cout << "Reading formulas and evidence predicates..." << endl;

    // Copy inMLNFile to workingMLNFile & app '#include "evid.db"'
  string::size_type bslash = inMLNFile.rfind("/");
  string tmp = (bslash == string::npos) ? 
               inMLNFile:inMLNFile.substr(bslash+1,inMLNFile.length()-bslash-1);
  char buf[100];
  sprintf(buf, "%s%d%s", tmp.c_str(), getpid(), ZZ_TMP_FILE_POSTFIX);
  wkMLNFile = buf;
  copyFileAndAppendDbFile(inMLNFile, wkMLNFile,
                          evidenceFilesArr, constFilesArr);

    // Parse wkMLNFile, and create the domain, MLN, database
  domain = new Domain;
  mln = new MLN();
  bool addUnitClauses = false;
  bool mustHaveWtOrFullStop = true;
  bool warnAboutDupGndPreds = true;
  bool flipWtsOfFlippedClause = true;
  bool allPredsExceptQueriesAreCW = false;
  //bool allPredsExceptQueriesAreCW = owPredNames.empty();
  Domain* forCheckingPlusTypes = NULL;

    // Parse as if lazy inference is set to true to set evidence atoms in DB
    // If lazy is not used, this is removed from DB
  cout<<wkMLNFile.c_str()<<endl;
  if (!runYYParser(mln, domain, wkMLNFile.c_str(), allPredsExceptQueriesAreCW, 
                   &owPredNames, &cwPredNames, &queryPredNames, addUnitClauses, 
                   warnAboutDupGndPreds, 0, mustHaveWtOrFullStop, 
                   forCheckingPlusTypes, true, flipWtsOfFlippedClause))
  {
    unlink(wkMLNFile.c_str());
    return -1;
  }

  unlink(wkMLNFile.c_str());
  const FormulaAndClausesArray* fca = mln->getFormulaAndClausesArray();
  for (int i = 0; i < fca->size(); i++)
  {
    IndexClauseHashArray* indexClauses = (*fca)[i]->indexClauses;
    for (int j = 0; j < indexClauses->size(); j++)
    {
      int idx = (*indexClauses)[j]->index;
      Clause* c = (*indexClauses)[j]->clause;
      cout << "idx " << idx << ": ";
      c->printWithWtAndStrVar(cout, domain);
      cout << endl;
    }
  }
    //////////////////////////// run inference /////////////////////////////////

    ///////////////////////// read & create query predicates ///////////////////
  Array<int> allPredGndingsAreQueries;
  Array<Array<Predicate* >* >* queryFormulas =  new Array<Array<Predicate*> *>;

    // Eager inference: Build the queries for the mrf
    // Lazy version evaluates the query string / file when printing out
  if (!aLazy)
  {
    if (queryFile.length() > 0)
    {
      cout << "Reading query predicates that are specified in query file..."
           << endl;
      bool ok = createQueryFilePreds(queryFile, domain, domain->getDB(),
                                     &queries, &knownQueries, queryFormulas);
      if (!ok) { cout<<"Failed to create query predicates."<<endl; exit(-1); }
    }

    allPredGndingsAreQueries.growToSize(domain->getNumPredicates(), false);
    if (queryPredsStr.length() > 0)
    {
      // unePreds = unknown non-evidence predicates
      // nePreds  = known non-evidence predicates
      GroundPredicateHashArray unePreds;
      GroundPredicateHashArray knePreds;
      bool ok = createComLineQueryPreds(queryPredsStr, domain, 
                                  domain->getDB(), &unePreds, &knePreds,
                                  &allPredGndingsAreQueries, queryFormulas);
      if (!ok) { cout<<"Failed to create query predicates."<<endl; exit(-1); }
      Array<Predicate*>* indexedGndings = new Array<Predicate*>();
      domain->getDB()->getIndexedGndings(indexedGndings,mln->getClause(0)->getPredicate(0),true,false);
      for(int i=0;i<indexedGndings->size();i++)
      {
    	  (*indexedGndings)[i]->print(cout,domain);
    	  cout<<endl;
      }
      for(int i=0;i<queryFormulas->size();i++)
      {
    	  for(int j=0;j<(*queryFormulas)[i]->size();j++)
    	  {
    		  (*(*queryFormulas)[i])[j]->print(cout,domain);
    		  cout<<"  ";
    	  }
    	  cout<<endl;
      }
      if (aisQueryEvidence)
      {
        // If the isQueryEvidence flag is set, then all query predicates
        // that are specified in the database are assumed to be the
        // set of queries we're actually interested in, while all
        // unspecified query predicates are assumed to be false evidence.
        // This is useful for doing inference with canopies,
        // e.g., (Singla & Domingos, 2005).

        // All unknown queries are actually false evidence 
        queryPredValues.clear();
        for (int predno = 0; predno < unePreds.size(); predno++) 
          domain->getDB()->setValue(unePreds[predno], FALSE);
        knownQueries = unePreds;

        // All known queries are actually unknown 
        for (int predno = 0; predno < knePreds.size(); predno++) 
        {
          TruthValue origValue 
              = domain->getDB()->setValue(knePreds[predno], UNKNOWN);
          queryPredValues.append(origValue);
        }
        queries = knePreds;
      }
      else
      {
        queries = unePreds;
        knownQueries = knePreds;
      }
    }
  }
  for(int i=0;i<queries.size();i++)
  {
	  queries[i]->print(cout,domain);
	  cout<<endl;
  }
  bool trackClauseTrueCnts = false;
  VariableState* state = NULL;
  HVariableState* hstate = NULL;
  FactorGraph* factorGraph = NULL;
  
  if (abpInfer || aefbpInfer || aoutputNetwork)
  {
    factorGraph = new FactorGraph(bpparams->lifted, bpparams->useHC,
                                  bpparams->useCT, bpparams->implicitRep,
                                  bpparams->hcCreateType,
                                  bpparams->hcCreateNoise, bpparams->lncIter,
                                  bpparams->noHCPredsStr,mln, domain,
                                  queryFormulas);
    inference = new BP(factorGraph, bpparams, queryFormulas);
  }
  else
  {
      // Create inference algorithm and state based on queries and mln / domain
    bool markHardGndClauses = true;
    bool trackParentClauseWts = false;
    if (aHybrid)
    {
        // Create inference algorithm and state based on queries and mln / domain
      bool markHardGndClauses = true;
      bool trackParentClauseWts = false;
	  hstate = new HVariableState(&queries, NULL, NULL,
		                          &allPredGndingsAreQueries,
		                          markHardGndClauses,
	                              trackParentClauseWts,
		                          mln, domain, aLazy);

	  //hstate->LoadContGroundedMLN(aContPartFile);
        // Ground out cont. formulas
      hstate->LoadContMLN();

	  //hstate->WriteGndPredIdxMap(aGndPredIdxMapFile);
	  hstate->bMaxOnly_ = (amapPos || amapAll);
    }
    else
    {
        // Lazy version: queries and allPredGndingsAreQueries are empty,
        // markHardGndClause and trackParentClauseWts are not used
      state = new VariableState(&queries, NULL, NULL,
                                &allPredGndingsAreQueries,
                                markHardGndClauses,
                                trackParentClauseWts,
                                mln, domain, aLazy);
    }

    if (aStartPt && aHybrid)
    {
      //hstate->LoadDisEviValuesFromRst(atestdis);
      hstate->loadDisEviValues();
      for (int i = 0; i < evidenceFilesArr.size(); i++)
        hstate->LoadContEviValues(evidenceFilesArr[i]);
      //hstate->loadContEviValues();
      hstate->setInitFromEvi(true);
    }

      // MAP inference, MC-SAT, Gibbs or Sim. Temp.
    if (amapPos || amapAll || amcsatInfer || agibbsInfer || asimtpInfer ||
        aHybrid || aSA)
    {
      if ((amapPos || amapAll) && !aHybrid)
      { // MaxWalkSat
          // When standalone MWS, numSolutions is always 1
          // (maybe there is a better way to do this?)
        mwsparams->numSolutions = 1;
        inference = new MaxWalkSat(state, aSeed, trackClauseTrueCnts, mwsparams);
      }
      else if ((amapPos || amapAll) && aHybrid && !amcsatInfer && !aSA)
	  {
        hstate->setProposalStdev(aProposalStdev);
          // maximizing all the numeric terms individually by l-bfgs, 
          // and cache the optimal solution for each term
        hstate->optimizeIndividualNumTerm();
        mwsparams->numSolutions = 1;
        mwsparams->heuristic = HMWS;
        inference = new HMaxWalkSat(hstate, aSeed, trackClauseTrueCnts,
                                    mwsparams);		
        HMaxWalkSat* p = (HMaxWalkSat*)inference;
        p->SetNoisePra(anumerator, adenominator);
        if (aMaxSeconds)
        {
          p->SetMaxSeconds(atof(aMaxSeconds));
        }

        if (aStartPt)
        {
		  hstate->setInitFromEvi(true);
		}
	  }
	  else if ((amapPos || amapAll) && aHybrid && !amcsatInfer && aSA)
	  {
		hstate->setProposalStdev(aProposalStdev);
		mwsparams->numSolutions = 1;
		mwsparams->heuristic = HSA;
		inference = new HMaxWalkSat(hstate, aSeed, trackClauseTrueCnts,
                                    mwsparams);

		HMaxWalkSat* p = (HMaxWalkSat*)inference;
		//p->SetNoisePra(anumerator, adenominator);
		p->setHeuristic(HSA);
		p->setSATempDownRatio(aSATempDownRatio);
		p->SetMaxSeconds(atof(aMaxSeconds));
		p->SetSAInterval(saInterval);
		if (aStartPt)
		{
          hstate->setInitFromEvi(true);
		}
	  }
	  else if (aHybrid && amcsatInfer)
	  {
		cout << "Creating HMCSAT instance." << endl;
		hstate->setProposalStdev(aProposalStdev);
		if (aContSamples == NULL)
		{
			//cout << "Numeric sample file is error." <<endl;
		}
		inference = new HMCSAT(hstate, aSeed, trackClauseTrueCnts, msparams);
		HMCSAT* p = (HMCSAT*) inference;
		p->SetContSampleFile(aContSamples);
        p->SetPrintVarsPerSample(aPrintSamplePerIteration);
	  }
      else if (amcsatInfer && !aHybrid)
      { // MC-SAT
        inference = new MCSAT(state, aSeed, trackClauseTrueCnts, msparams,
                              queryFormulas);
      }
      else if (asimtpInfer && !aHybrid)
      { // Simulated Tempering
          // When MWS is used in Sim. Temp., numSolutions is always 1
          // (maybe there is a better way to do this?)
        mwsparams->numSolutions = 1;
        inference = new SimulatedTempering(state, aSeed, trackClauseTrueCnts,
                                           stparams);
      }
      else if (agibbsInfer && !aHybrid)
      { // Gibbs sampling
          // When MWS is used in Gibbs, numSolutions is always 1
          // (maybe there is a better way to do this?)
        mwsparams->numSolutions = 1;
        inference = new GibbsSampler(state, aSeed, trackClauseTrueCnts,
                                     gibbsparams);
      }
    }
  }
  return 1;
}
 
  // Typedefs
typedef hash_map<string, const Array<const char*>*, HashString, EqualString> 
StringToStrArrayMap;

#endif
