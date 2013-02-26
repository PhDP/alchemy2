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
%{
#define YYSTYPE int
#define YYDEBUG 1
%}


/*************************** Declarations ***************************/
// the token 'error' is declared by Bison
%token ZZ_NUM
%token ZZ_DOTDOTDOT
%token ZZ_STRING
%token ZZ_INCLUDE

%token ZZ_PREDICATE
%token ZZ_FUNCTION
%token ZZ_CONSTANT
%token ZZ_VARIABLE
%token ZZ_TYPE

%left '-' '+'

%token ZZ_FORALL
%token ZZ_EXIST
%left ZZ_EQUIV
%left ZZ_IMPLY
%left 'v'
%left '^'
%left '!'
%left '*' '/' '%'

%glr-parser
%expect 16
%error-verbose

// This must be in a %code {} section since bison 2.4 so it appears before tokens
// The comments following { and } must remain intact as they
// are used to insert the proper directive in replacefoly.pl
%code { //bisonopencode
#include "fol.h"
#include "follex.cpp"

  // 0: no output; 1,2: increasing order of verbosity
int folDbg = 0;
//int folDbg = 1;
//int folDbg = 2;
} //bisonclosecode

%% 

/*************************** Grammar **************************/
// printouts after "input" in the following rules cause shift/reduce conflicts
input:
  { if (folDbg >= 2) printf("input: empty\n"); }
  // empty
| input error '\n' 
  { 
    yyerrok; // tell Bison not to suppress any errors
    const char* tok; 
    while (true)
    {
      tok = zztokenList.removeLast();
      if (strcmp(tok,"\n")==0) { delete tok; break; }
      delete tok;
    }
  }
| input newline
| input include // newline is in include rule
| input type_declaration newline
| input numeric_type_declaration newline
| input predicate_declaration newline
| input function_declaration newline 
| input at //{ zzconsumeToken(zztokenList,"@"); }
  pd_not_qs predicate_definition
| input at //{ zzconsumeToken(zztokenList,"@"); }
  function_definition
| input at newline
| 
    // input weight sentence
  input  
  { 
    if (folDbg >= 2) printf("input: weight\n"); 
    zzreset();
  }
  weight
  {
    if (folDbg >= 2) printf("input: utility\n"); 
  }
  utility
  {  
    if (folDbg >= 2) printf("input: sentence\n");
      // the states should be reset because a parse error may have occurred
  }
  sentence
  continuous_part
  fullstop
  {
    ListObj* formula;
    zzassert(zzoldNewVarList.size()==0,"expected zzoldNewVarList.size()==0");
    zzassert(zzformulaListObjs.size()==1,"expected zzformulaListObjs.size()=1");
   	formula = zzformulaListObjs.top(); zzformulaListObjs.pop();

   	zzdetermineEqPredTypes(formula);//set the any unknown type of '=' predicates
    zzeqPredList.clear();
   	zzdetermineIntPredTypes(formula);//set the any unknown type of internal predicates
    zzintPredList.clear();
   	zzdetermineIntFuncTypes(formula);//set the any unknown type of internal functions
    zzintFuncList.clear();
    zzsetPlusVarTypeId();// set typeIds of variables with pluses before them
    zzcheckVarNameToIdMap();
    
    if (zzhasFullStop && zzwt != NULL)
    {
      zzerr("A weight must not be specified for a formula that is "
            "terminated with a period (i.e. the formula is hard).");
      delete zzwt; zzwt = NULL;
    }
	
      // at this point we are sure we are dealing with a formula

      //if there are errors, we can keep this formula for further processing
    if (zznumErrors == 0)
    {
        //defer converting the formula to CNF until we have read all the
        //.db files and are sure that we know all the constants to ground
        //formula's clauses
      ZZFormulaInfo* epfi 
        = new ZZFormulaInfo(formula, zzformulaStr, zzfdnumPreds, zzwt, 
                            zzdefaultWt, zzdomain, zzmln, zzvarNameToIdMap, 
                            zzplusVarMap, zznumAsterisk,
                            zzhasFullStop, zzreadHardClauseWts, 
                            zzmustHaveWtOrFullStop, zzinIndivisible,
                            zzisHybrid, zzcontPred, zzmean,
                            zzhasWeightFullStop, zzutil);
      zzformulaInfos.append(epfi); 
    }

    if (zzwt) { delete zzwt; zzwt = NULL; }
    if (zzutil) { delete zzutil; zzutil = NULL; }
  }
  newline
;


continuous_part:
  // empty 
| 
  '*'
  {
    // Indicator
    if (folDbg >= 2) printf("sentence: indicator function\n");
    zzconsumeToken(zztokenList,"*"); 
    if (folDbg >= 1) printf("* ");
    zzisHybrid = true;
    zzformulaStr.append(" * ");
  }
  numeric_term
  {
  	// Real-valued
    if (folDbg >= 2) printf("sentence: numeric term\n");
  }


at: 
  '@'    { zzconsumeToken(zztokenList,"@"); }
| at '@' { zzconsumeToken(zztokenList,"@"); }


fullstop:
  // empty 
| 
  '.' 
{ 
  if (folDbg >= 1) printf(".\n"); zzconsumeToken(zztokenList,"."); 
  zzassert(!zzhasFullStop, "expecting no full stop");
  zzhasFullStop = true;
  zzformulaStr.append(".");
}


newline:
  '\n' { if (folDbg >= 1) printf("\\n\n"); zzconsumeToken(zztokenList,"\n"); }
| 
  '\r' { if (folDbg >= 1) printf("\\r\n"); zzconsumeToken(zztokenList,"\r"); }


optnewline:
  //empty
|
  '\n' { if (folDbg >= 1) printf("\\n\n"); zzconsumeToken(zztokenList,"\n"); }
| 
  '\r' { if (folDbg >= 1) printf("\\r\n"); zzconsumeToken(zztokenList,"\r"); }



/*************************** include statement ****************************/

include: ZZ_INCLUDE ZZ_STRING nnewline
{
  const char* inc = zztokenList.removeLast();
  const char* str = zztokenList.removeLast();
  const char* nl = zztokenList.removeLast();
  zzassert(strcmp(inc,"#include")==0,"expecting #include keyword");

  if (folDbg >= 1) printf("#include %s ", str);

  string s(str);
  if (s.length() == 0)
  {
    zzerr("empty string used with #include");
    delete [] inc;
    delete [] str;
    delete [] nl;
    break;
  }

  zzassert(s.at(0) == '"' && s.at(s.length()-1) == '"', "no enclosing \"");
  s = s.substr(1, s.length()-2);
  int len = s.length();

  // if it is a .cpp file, then we are dealing with linked-in functions and predicates
  if (s.at(len-4)=='.' && s.at(len-3)=='c' && s.at(len-2)=='p' &&
	  s.at(len-1)=='p')
  {
    zzcompileFunctions(str);
    zzusingLinkedPredicates = true;    
    zzusingLinkedFunctions = true;    
    break;
  }

  FILE* newin = fopen(s.c_str(), "r" );
  if (newin)
  {
    zzinStack.push(ZZFileState(yyin, string(zzinFileName), zznumCharRead, 
                               zzline, zzcolumn)); 
    ungetc('\n', newin); // pretend that file begins with a newline
    zzline = 1;
    zzcolumn = -1;
    zzline--;
    zzinFileName = str;
    yyrestart(newin); 
    zznumCharRead = 0;

      // if it is a .db file containing ground predicates
    if ((s.at(len-3)=='.' && s.at(len-2)=='d' && s.at(len-1)=='b'))
      zzparseGroundPred = true;
  } 
  else
    zzerr("Failed to open file %s.", str);
  
  delete [] inc;
  delete [] str;
  delete [] nl;
}
;

nnewline:
  '\n' { if (folDbg >= 1) printf("\\n\n"); }
| 
  '\r' { if (folDbg >= 1) printf("\\r\n"); }


/********************** type and constant declaration ********************/

// The first time a type is declared it is not seen before, so flex returns it
// as a ZZ_VARIABLE instead of ZZ_TYPE
// is_variable_type '=' '{' constant_declarations  '}'
type_declaration: 
is_variable_type
'=' '{'  
{ 
  zzconsumeToken(zztokenList,"=");
  zzconsumeToken(zztokenList,"{");
  if (folDbg >= 1) printf("= { "); 
  if (folDbg >= 2) printf("type_declarations: constant_declarations\n");
}
constant_declarations 
'}'
{
  if (folDbg >= 1) printf("} ");           
  zzconsumeToken(zztokenList, "}");
  delete [] zztypeName;
  zztypeName = NULL;
}
;


is_variable_type:
ZZ_VARIABLE 
{
  const char* idf = zztokenList.removeLast();
  if (folDbg >= 1) printf("type t_%s ", idf);
  zzassert(!zzdomain->isType(idf), "expecting type to be undefined");
  int id = zzaddTypeToDomain(zzdomain, idf);
  zzassert(id >= 0,"expecting id >= 0");
  zzassert(zztypeName==NULL,"expecting zztypeName==NULL");
  zztypeName = new char[strlen(idf)+1];
  strcpy(zztypeName, idf);
  delete [] idf;
}
|
ZZ_TYPE 
{
  const char* idf = zztokenList.removeLast();
  if (folDbg >= 1) printf("type t_%s ", idf);
  zzassert(zzdomain->isType(idf),"expecting type to be defined");
  //zzwarn("Type %s has been declared before.",idf);
  zzassert(zztypeName==NULL,"expecting zztypeName==NULL");
  zztypeName = new char[strlen(idf)+1];
  strcpy(zztypeName, idf);
  delete [] idf;
}


// ZZ_VARIABLE is an undeclared constant.
// empty | constant_declarations constant_sep ZZ_VARIABLE
constant_declarations:   
  // empty 
|
  //{if (folDbg>=2) printf("constant_declarations: constant_declarations\n");}
  constant_declarations constant_sep 
  { if (folDbg >= 2) printf("constant_declarations: ZZ_VARIABLE\n"); }
  variable_or_string_or_constant 
  {
    const char* vsc = zztokenList.removeLast();
    if (folDbg >= 1) printf("cd_%s ", vsc);
    zzaddConstantToDomain(vsc, zztypeName);
    delete [] vsc;
  }
;


variable_or_string_or_constant: ZZ_VARIABLE | ZZ_STRING | ZZ_CONSTANT
;



constant_sep:   
  //empty
| 
  ','  { zzconsumeToken(zztokenList, ","); }
| 
  '\n' ',' '\n' 
  { 
    zzconsumeToken(zztokenList,"\n");
    zzconsumeToken(zztokenList,",");
    zzconsumeToken(zztokenList,"\n");
  }
| 
  '\n' ',' 
  { 
    zzconsumeToken(zztokenList,"\n");
    zzconsumeToken(zztokenList,",");
  }
| 
  ',' '\n'      
  { 
    zzconsumeToken(zztokenList,",");
    zzconsumeToken(zztokenList,"\n");
  }
;


numeric_type_declaration: 
is_variable_type
'=' '{'  
{
  if (folDbg >= 2) printf("numeric_type_declarations: \n");
  zzconsumeToken(zztokenList,"=");
  zzconsumeToken(zztokenList,"{");  
}
ZZ_NUM
{
  const char* numStr = zztokenList.removeLast();
  if (folDbg >= 1) printf(" %s ", numStr);
  
  double d = atof(numStr);
  zzinitialNumDeclaration = int(d);
  if (d != zzinitialNumDeclaration) zzerr("constant %s must be an integer", numStr);

  char constStr[100];
  zzcreateIntConstant(constStr, zztypeName, zzinitialNumDeclaration);
  zzaddConstantToDomain(constStr, zztypeName);

  delete [] numStr;
}
','
{
  zzconsumeToken(zztokenList, ",");
}
numeric_types
'}'
{
  zzconsumeToken(zztokenList, "}");
}
;

numeric_types:
ZZ_DOTDOTDOT ',' ZZ_NUM
{
  //const char* numStr1 = zztokenList.removeLast();
  zzconsumeToken(zztokenList, "...");
  zzconsumeToken(zztokenList, ",");
  const char* numStr2 = zztokenList.removeLast();
  if (folDbg >= 1) printf("= ... , %s } ", numStr2);
  
  double d2 = atof(numStr2);
  int i2 = int(d2);
  if (d2 != i2) zzerr("constant %s must be an integer", numStr2);
  if (zzinitialNumDeclaration > i2)
    zzerr("first constant cannot be larger than last constant %s", numStr2);

  for (int i = zzinitialNumDeclaration + 1; i <= i2; i++)
  {
    char constStr[100];
    zzcreateIntConstant(constStr, zztypeName, i);
    zzaddConstantToDomain(constStr, zztypeName);
  }

  delete [] numStr2; delete [] zztypeName;
  zztypeName = NULL;
}
|
single_numeric_types
{
  delete [] zztypeName;
  zztypeName = NULL;
}
;

// ZZ_NUM ',' single_numeric_types | ZZ_NUM
single_numeric_types:
  single_numeric_types
  ',' 
  {  
    zzconsumeToken(zztokenList, ",");
    if (folDbg >= 1) printf(", "); 
    if (folDbg >= 2) printf("single_numeric_types: single_numeric_type\n"); 
  }
  single_numeric_type
| 
  single_numeric_type
;

single_numeric_type: ZZ_NUM
{
  const char* numStr = zztokenList.removeLast();
  if (folDbg >= 1) printf(" %s ", numStr);
  
  double d = atof(numStr);
  int i = int(d);
  if (d != i) zzerr("constant %s must be an integer", numStr);

  char constStr[100];
  zzcreateIntConstant(constStr, zztypeName, i);
  zzaddConstantToDomain(constStr, zztypeName);

  delete [] numStr;
}
;

/******************* predicate and function declaration **********************/

// The first time a predicate is declared, flex returns it as a variable because
// it has not been seen before
// ZZ_VARIABLE '(' types ')'
predicate_declaration:
//  if (folDbg >= 2) printf("predicate_declaration: ZZ_VARIABLE\n");
ZZ_VARIABLE 
{
  const char* predName = zztokenList.removeLast();

  if (folDbg >= 1) printf("ZZ_PREDICATE pc_%s ", predName);  
    //predicate has not been declared a function
  zzassert(zzdomain->getFunctionId(predName) < 0, 
           "not expecting pred name to be declared as a function name");
  zzassert(zzpredTemplate==NULL,"expecting zzpredTemplate==NULL");
  zzpredTemplate = new PredicateTemplate();
  zzpredTemplate->setName(predName);
  delete [] predName;
}
'('    { if (folDbg >= 1) printf("( "); zzconsumeToken(zztokenList,"("); }
types  { if (folDbg >= 2) printf("predicate_declaration: types\n"); }
')' 
{  
  if (folDbg >= 1) printf(") "); 

  zzconsumeToken(zztokenList,")");

  zzassert(zzpredTemplate, "not expecting zzpredTemplate==NULL");
  int id = zzdomain->addPredicateTemplate(zzpredTemplate);
  zzassert(id >= 0, "expecting pred template id >= 0");
  zzpredTemplate->setId(id);
  zzpredTemplate = NULL;
}
|
ZZ_VARIABLE 
{
  const char* predName = zztokenList.removeLast();

  if (folDbg >= 1) printf("ZZ_PREDICATE pc_%s ", predName);  
    //predicate has not been declared a function
  zzassert(zzdomain->getFunctionId(predName) < 0, 
           "not expecting pred name to be declared as a function name");
  zzassert(zzpredTemplate==NULL,"expecting zzpredTemplate==NULL");
  zzpredTemplate = new PredicateTemplate();
  zzpredTemplate->setName(predName);
  delete [] predName;

	// Declare this predicate to have type "AlchemyPropositionalType"
  const char* varName = Domain::PROPOSITIONAL_TYPE;
  if (folDbg >= 1) printf("t_%s ", varName);
  if (!zzdomain->isType(varName))
  {
    int id = zzaddTypeToDomain(zzdomain, varName);
    zzassert(id >= 0, "expecting var id >= 0");
  }
  zzaddType(varName, zzpredTemplate, NULL, false, zzdomain);

  zzassert(zzpredTemplate, "not expecting zzpredTemplate==NULL");
  int templateId = zzdomain->addPredicateTemplate(zzpredTemplate);
  zzassert(templateId >= 0, "expecting pred template id >= 0");
  zzpredTemplate->setId(templateId);
  zzpredTemplate = NULL;
}
;


// The first time a function is declared, flex returns it as a variable because
// it has not been seen before
// function_return_type ZZ_VARIABLE '(' types ')'
function_declaration: 
function_return_type 
{ 
  const char* retTypeName = zztokenList.removeLast();
  if (folDbg >= 1) printf("ZZ_FUNCTION t_%s ", retTypeName); 
  if (folDbg >= 2) printf("function_declaration: ZZ_VARIABLE\n"); 

  if (!zzdomain->isType(retTypeName))
  {
    int id = zzaddTypeToDomain(zzdomain, retTypeName);
    zzassert(id >= 0, "expecting retTypeName's id >= 0");
  }

  zzassert(zzfuncTemplate==NULL, "expecting zzfuncTemplate==NULL");
  zzfuncTemplate = new FunctionTemplate();
  zzfuncTemplate->setRetTypeName(retTypeName,zzdomain);

  // We are creating a new predicate as well
  zzassert(zzpredTemplate==NULL,"expecting zzpredTemplate==NULL");
  zzpredTemplate = new PredicateTemplate();
  zzaddType(retTypeName, zzpredTemplate, NULL, false, zzdomain);

  delete [] retTypeName;
}
ZZ_VARIABLE 
{
  const char* funcName = zztokenList.removeLast();

  if (folDbg >= 1) printf("fc_%s ", funcName); 
  zzassert(zzdomain->getPredicateId(funcName) < 0, 
           "not expecting func name to be declared as pred name");
  zzfuncTemplate->setName(funcName);

  // Predicate name is PredicateTemplate::ZZ_RETURN_PREFIX + function name
  char* predName;
  predName = (char *)malloc((strlen(PredicateTemplate::ZZ_RETURN_PREFIX) +
  							strlen(funcName) + 1)*sizeof(char));
  strcpy(predName, PredicateTemplate::ZZ_RETURN_PREFIX);
  strcat(predName, funcName);
    	
  // Check that predicate has not been declared a function
  zzassert(zzdomain->getFunctionId(predName) < 0, 
           "not expecting pred name to be declared as a function name");
  zzassert(zzpredTemplate,"expecting zzpredTemplate!=NULL");
  zzpredTemplate->setName(predName);

  free(predName);
  delete [] funcName;
}
'(' 
{
  zzconsumeToken(zztokenList,"(");
  if (folDbg >= 1) printf("( "); 
  if (folDbg >= 2) printf("function_declaration: types\n"); 
}
types ')' 
{  
  zzconsumeToken(zztokenList,")");
  if (folDbg >= 1) printf(") "); 
  zzassert(zzfuncTemplate, "expecting zzfuncTemplate != NULL");
  int id = zzdomain->addFunctionTemplate(zzfuncTemplate);
  zzassert(id >= 0, "expecting function template's id >= 0");
  zzfuncTemplate->setId(id);
  zzfuncTemplate = NULL;

  zzassert(zzpredTemplate, "not expecting zzpredTemplate==NULL");
  // Predicate could already be there
  if (!zzdomain->isPredicate(zzpredTemplate->getName()))
  {
    int predId = zzdomain->addPredicateTemplate(zzpredTemplate);
    zzassert(predId >= 0, "expecting pred template id >= 0");
    zzpredTemplate->setId(predId);
  }
  zzpredTemplate = NULL;

}
;

  // tokens consumed in function_declaration
function_return_type: ZZ_TYPE | ZZ_VARIABLE
;


types:  //{ if (folDbg >= 2) printf("types: types\n"); }
  types ',' 
  { if (folDbg >= 1) printf(", "); zzconsumeToken(zztokenList,","); }
  type_code
| 
  types ','  
  { if (folDbg >= 1) printf(", "); zzconsumeToken(zztokenList,","); }
  variable_code
| type_code
| variable_code
;


type_code: ZZ_TYPE
{  
  const char* ttype = zztokenList.removeLast();
  if (folDbg >= 1) printf("t_%s ", ttype);
  zzaddType(ttype, zzpredTemplate, NULL, false, zzdomain);
  // If we are in a function definition, then add the type to the function too
  if (zzfuncTemplate)
	zzaddType(ttype, NULL, zzfuncTemplate, false, zzdomain);
  delete [] ttype;
}
exist_unique

variable_code: ZZ_VARIABLE
{
  const char* varName = zztokenList.removeLast();
  if (folDbg >= 1) printf("t_%s ", varName);
  zzassert(!zzdomain->isType(varName), "expecting varName to be not a type");
  int id = zzaddTypeToDomain(zzdomain, varName);
  zzassert(id >= 0, "expecting var id >= 0");
  zzaddType(varName, zzpredTemplate, NULL, false, zzdomain);
  // If we are in a function definition, then add the type to the function too
  if (zzfuncTemplate)
	zzaddType(varName, NULL, zzfuncTemplate, false, zzdomain);
  delete [] varName;
}
exist_unique

exist_unique:
  // No '!'
  {
  }
| '!'
  {
    zzconsumeToken(zztokenList, "!");
    if (folDbg >= 1) printf("! "); 
    if (zzfuncTemplate)
      zzerr("'!' cannot be used in a function declaration.");
    zzassert(zzpredTemplate, "'!' used outside a predicate declaration.");

    int termIdx = zzpredTemplate->getNumTerms()-1;
    zzpredTemplate->addUniqueVarIndex(termIdx);
  }


/***************** predicate and function groundings ***********************/

pd_not_qs:  
  // empty
  { 
    zztrueFalseUnknown = 't';
  }
| '!'  
  { 
    zzconsumeToken(zztokenList,"!"); 
    if (folDbg >= 1) printf("! "); 
    zztrueFalseUnknown = 'f';
  }
| '?'  
  { 
    zzconsumeToken(zztokenList,"?"); 
    if (folDbg >= 1) printf("? "); 
    zztrueFalseUnknown = 'u';
  }
;


// ZZ_PREDICATE '(' constants_in_groundings ')'
predicate_definition:  
ZZ_PREDICATE 
{
  const char* predName = zztokenList.removeLast();
  if (folDbg >= 1) printf("pd_%s ", predName); 
  
  zzcreatePred(zzpred, predName);
  if (zztrueFalseUnknown == 't')      zzpred->setTruthValue(TRUE);
  else if (zztrueFalseUnknown == 'f') zzpred->setTruthValue(FALSE);
  else if (zztrueFalseUnknown == 'u') zzpred->setTruthValue(UNKNOWN);
  else { zzassert(false, "expecting t,f,u"); }

  delete [] predName;
}
nothing_or_constants
{
  zzcheckPredNumTerm(zzpred);
  int predId = zzpred->getId();
  hash_map<int,PredicateHashArray*>::iterator it;
  if ((it=zzpredIdToGndPredMap.find(predId)) == zzpredIdToGndPredMap.end())
    zzpredIdToGndPredMap[predId] = new PredicateHashArray;
    
  if (zzrealValue != NULL)
    zzpred->setRealValue(*zzrealValue);
    
  PredicateHashArray* pha = zzpredIdToGndPredMap[predId];
  if (pha->append(zzpred) < 0)
  {
    int a = pha->find(zzpred);
    zzassert(a >= 0, "expecting ground predicate to be found");
    string origTvStr = (*pha)[a]->getTruthValueAsStr();
    (*pha)[a]->setTruthValue(zzpred->getTruthValue());
    string newTvStr = (*pha)[a]->getTruthValueAsStr();

    if (zzwarnDuplicates)
    {
      ostringstream oss;
      oss << "Duplicate ground predicate "; zzpred->print(oss, zzdomain); 
      oss << " found. ";
      if (origTvStr.compare(newTvStr) != 0)
        oss << "Changed its truthValue from " << origTvStr << " to " <<newTvStr 
            << endl;
      zzwarn(oss.str().c_str());
    }
    delete zzpred;
  }
  zzpred = NULL;
}
;

nothing_or_constants:
  // Nothing (propositional case)
{
  // Add the one constant for propositional case
  const char* constName = Domain::PROPOSITIONAL_CONSTANT;
  if (folDbg >= 1) printf("cg_%s ", constName);
  zzaddConstantToPredFunc(constName);
}
|
'(' 
{ 
  zzconsumeToken(zztokenList,"("); 
  if (folDbg >= 1) printf("( "); 
  if (folDbg >= 2) printf("predicate_definition: constants_in_groundings\n");
}
constants_in_groundings ')'
{ 
  zzconsumeToken(zztokenList,")"); 
  if (folDbg >= 1) printf(")\n"); 
}
real_value
;

real_value:
  // whitespace between pred and number was eaten up
ZZ_NUM
{
  const char* value = zztokenList.removeLast();
  if (folDbg >= 1) printf("rv_%f ", atof(value));
  if (zzrealValue) delete zzrealValue;
  zzrealValue = new double(atof(value));
  delete [] value;
}
|
  //empty
;

// function_return_constant '=' ZZ_FUNCTION '(' constants_in_groundings ')'
function_definition:
function_return_constant
'='
{ 
  zzconsumeToken(zztokenList,"="); 
  if (folDbg >= 1) printf("= "); 
  if (folDbg >= 2) printf("function_definition: constants_in_groundings\n");
}
ZZ_FUNCTION
{
  // Predicate name is PredicateTemplate::ZZ_RETURN_PREFIX + function name
  const char* funcName = zztokenList.removeLast();
  char* predName;
  predName = (char *)malloc((strlen(PredicateTemplate::ZZ_RETURN_PREFIX) +
  							strlen(funcName) + 1)*sizeof(char));
  strcpy(predName, PredicateTemplate::ZZ_RETURN_PREFIX);
  strcat(predName, funcName);

  if (folDbg >= 1) printf("pd_%s ", predName); 
  
  zzcreatePred(zzpred, predName);
  zzpred->setTruthValue(TRUE);

  char constName[100];
  char* constString;
  if (zztmpReturnNum)
  {
  	zzcreateAndCheckNumConstant(zztmpReturnConstant, zzfunc, zzpred, zzdomain, constName);
    if (constName == NULL)
    {
      constString = (char *)malloc((strlen(zztmpReturnConstant) + 1)*sizeof(char));
      strcpy(constString, zztmpReturnConstant);
    } else
    {
	  constString = (char *)malloc((strlen(constName) + 1)*sizeof(char));
      strcpy(constString, constName);
    }
  }
  else
  {
    constString = (char *)malloc((strlen(zztmpReturnConstant) + 1)*sizeof(char));  
  	strcpy(constString, zztmpReturnConstant);
  }

  //Add return constant parsed earlier
  zzaddConstantToPredFunc(constString);

  zztmpReturnNum = false;
  free(zztmpReturnConstant);
  delete [] funcName;
  free(predName);
  free(constString);
}
'(' 
{ 
  zzconsumeToken(zztokenList,"("); 
  if (folDbg >= 1) printf("( "); 
  if (folDbg >= 2) printf("function_definition: constants_in_groundings\n");
}
constants_in_groundings ')'
{ 
  zzconsumeToken(zztokenList,")"); 
  if (folDbg >= 1) printf(")\n"); 
  
  zzcheckPredNumTerm(zzpred);
  int predId = zzpred->getId();
  hash_map<int,PredicateHashArray*>::iterator it;
  if ((it=zzpredIdToGndPredMap.find(predId)) == zzpredIdToGndPredMap.end())
    zzpredIdToGndPredMap[predId] = new PredicateHashArray;
  
  PredicateHashArray* pha = zzpredIdToGndPredMap[predId];
  if (pha->append(zzpred) < 0)
  {
    int a = pha->find(zzpred);
    zzassert(a >= 0, "expecting ground predicate to be found");
    string origTvStr = (*pha)[a]->getTruthValueAsStr();
    (*pha)[a]->setTruthValue(zzpred->getTruthValue());
    string newTvStr = (*pha)[a]->getTruthValueAsStr();

    if (zzwarnDuplicates)
    {
      ostringstream oss;
      oss << "Duplicate ground predicate "; zzpred->print(oss, zzdomain); 
      oss << " found. ";
      if (origTvStr.compare(newTvStr) != 0)
        oss << "Changed its truthValue from " << origTvStr << " to " <<newTvStr 
            << endl;
      zzwarn(oss.str().c_str());
    }
    //delete zzpred;
  }

  // Insert FALSE for all other return values

  zzpred = NULL;
}

// is_constant | constants_in_groundings constant_sep is_constant
constants_in_groundings:  
  is_constant_string_num_variable
  {     
    const char* constName = zztokenList.removeLast();
    if (folDbg >= 1) printf("cg_%s ", constName);
    zzaddConstantToPredFunc(constName);
    delete [] constName;
  }
| 
  constants_in_groundings constant_sep is_constant_string_num_variable  
  {     
    const char* constName = zztokenList.removeLast();
    if (folDbg >= 1) printf("cg_%s ", constName);
    zzaddConstantToPredFunc(constName);
    delete [] constName;
  }


is_constant_string_num_variable:
  ZZ_CONSTANT
| ZZ_STRING  
  {
    if (zzconstantMustBeDeclared)
      zzerr("Constant %s must be declared before it is used.",
            zztokenList.back());
  }
| ZZ_NUM
  {
    const char* intStr = zztokenList.removeLast();
    char constName[100];
    zzcreateAndCheckNumConstant(intStr, zzfunc, zzpred, zzdomain, constName);
    if (constName == NULL) zztokenList.addLast(intStr);
    else                   zztokenList.addLast(constName);
    delete [] intStr;
  }
| ZZ_VARIABLE
  {
    if (zzconstantMustBeDeclared)
      zzerr("Constant %s must be declared before it is used",
            zztokenList.back());
  }

function_return_constant:
  ZZ_CONSTANT
  {
    const char* tmp = zztokenList.removeLast();
	zztmpReturnConstant = (char *)malloc((strlen(tmp) + 1)*sizeof(char));
  	strcpy(zztmpReturnConstant, tmp);
  	zztmpReturnNum = false;
  	delete []tmp;
  	if (folDbg >= 1) printf("ic_%s ", zztmpReturnConstant); 
  }
| ZZ_STRING  
  {
    if (zzconstantMustBeDeclared)
      zzerr("Constant %s must be declared before it is used.",
            zztokenList.back());
    const char* tmp = zztokenList.removeLast();
	zztmpReturnConstant = (char *)malloc((strlen(tmp) + 1)*sizeof(char));
  	strcpy(zztmpReturnConstant, tmp);
  	zztmpReturnNum = false;
  	delete []tmp;
  	if (folDbg >= 1) printf("ic_%s ", zztmpReturnConstant);
  }
| ZZ_NUM
  {
  	const char* tmp = zztokenList.removeLast();
  	zztmpReturnConstant = (char *)malloc((strlen(tmp) + 1)*sizeof(char));
  	strcpy(zztmpReturnConstant, tmp);
  	zztmpReturnNum = true;
  	delete []tmp;
  	if (folDbg >= 1) printf("icnum_%s ", zztmpReturnConstant); 
  }
| ZZ_VARIABLE
  {
    if (zzconstantMustBeDeclared)
      zzerr("Constant %s must be declared before it is used",
            zztokenList.back());
    const char* tmp = zztokenList.removeLast();
	zztmpReturnConstant = (char *)malloc((strlen(tmp) + 1)*sizeof(char));
  	strcpy(zztmpReturnConstant, tmp);
  	zztmpReturnNum = false;
  	delete []tmp;
  	if (folDbg >= 1) printf("ic_%s ", zztmpReturnConstant);
  }

/**************************** sentence *************************************/

weight:   
  // empty 
| 
  ZZ_NUM
  {
    const char* wt = zztokenList.removeLast();
    if (folDbg >= 1) printf("n_%f ", atof(wt));
    if (zzwt) delete zzwt;
    zzwt = new double(atof(wt));
    delete [] wt;
  }
  weight_fullstop
;

weight_fullstop:
  // empty 
| 
  '.' 
{ 
  if (folDbg >= 1) printf("."); zzconsumeToken(zztokenList,"."); 
  zzassert(!zzhasWeightFullStop, "expecting no full stop");
  zzhasWeightFullStop = true;
  zzformulaStr.append(".");
}

utility:
  // empty
|
  ':'
{
  if (folDbg >= 1) printf(":"); zzconsumeToken(zztokenList,":"); 
//  zzassert(!zzhasUtility, "expecting no utility");
//  zzhasUtility = true;
  zzformulaStr.append(":");
}
  ZZ_NUM
  {
    const char* util = zztokenList.removeLast();
    if (folDbg >= 1) printf("n_%f ", atof(util));
    if (zzutil) delete zzutil;
    zzutil = new double(atof(util));
    delete [] util;
  }

// You must explicitly specify the connectives between the "sentence"s in order
// for Bison to effect the left-associativity of the connectives
sentence:
  '['  
  { 
      // Square Brackets
    zzconsumeToken(zztokenList,"["); 
    if (folDbg >= 1) printf("[ "); 
    if (folDbg >= 2) printf("sentence: '[' sentence\n");
    zzformulaStr.append("[");
    zzinIndivisible = true;
  }
  sentence 
  ']'  
  { 
    zzconsumeToken(zztokenList,"]"); 
    if (folDbg >= 1) printf("] "); 
    zzformulaStr.append("]");
    //zzinIndivisible = false;
  }
| 
  '('  
  { 
    zzconsumeToken(zztokenList,"("); 
    if (folDbg >= 1) printf("( "); 
    if (folDbg >= 2) printf("sentence: '(' sentence\n");
    zzformulaStr.append("(");
  }
  sentence 
  ')'  
  { 
    zzconsumeToken(zztokenList,")"); 
    if (folDbg >= 1) printf(") "); 
    zzformulaStr.append(")");
  }
| 
  { if (folDbg >= 2) printf("sentence: atomic_sentence\n"); }
  atomic_sentence
|
  sentence ZZ_IMPLY
  {
    const char* imp = zztokenList.removeLast(); 
    if (folDbg >= 1) printf("=> ");
    zzformulaStr.append(" => ");
    delete [] imp;
  }
  optnewline
  sentence
  { zzcreateListObjFromTopTwo(zzformulaListObjs, "=>"); }
  %dprec 1
|  
  sentence ZZ_EQUIV
  { 
    const char* eq = zztokenList.removeLast(); 
    if (folDbg >= 1) printf("<=> ");
    zzformulaStr.append(" <=> ");
    delete [] eq;  
  }
  optnewline
  sentence
  { zzcreateListObjFromTopTwo(zzformulaListObjs, "<=>"); }
  %dprec 1
|
  sentence 'v' 
  { 
    zzconsumeToken(zztokenList,"v"); 
    if (folDbg >= 1) printf("v "); 
    zzformulaStr.append(" v ");
  }
  optnewline
  sentence
  { zzcreateListObjFromTopTwo(zzformulaListObjs, "v"); }
  %dprec 1
|  
  sentence '^' 
  { 
    zzconsumeToken(zztokenList,"^"); 
    if (folDbg >= 1) printf("^ "); 
    zzformulaStr.append(" ^ ");
  }
  optnewline
  sentence
  { zzcreateListObjFromTopTwo(zzformulaListObjs, "^"); }
  %dprec 1
| 
  // quantifier variables sentence
  { if (folDbg >= 2) printf("sentence: quantifier\n"); }
  quantifier  
  { 
    if (folDbg >= 2) printf("sentence: variables\n"); 
    zzformulaListObjs.push(new ListObj);
    pair<StringToStringMap*,int> pr(new StringToStringMap,zzscopeCounter++);
    zzoldNewVarList.push_front(pr);
  }
  variables   { if (folDbg >= 2) printf("sentence: sentence\n"); }
  sentence    
  { 
    zzcreateListObjFromTopThree(zzformulaListObjs);
    pair<StringToStringMap*, int> pr = zzoldNewVarList.front();
    zzoldNewVarList.pop_front();
    delete pr.first;
  }
  %dprec 2
| 
  '!'  // '!' sentence
  { 
    zzassert(!zzisNegated,"expecting !zzisNegated");
    zzisNegated = true;

    zzconsumeToken(zztokenList,"!");
    if (folDbg >= 1) printf("! "); 
    if (folDbg >= 2) printf("sentence: sentence\n");
    zzformulaStr.append("!");
  }
  sentence
  { zzcreateListObjFromTop(zzformulaListObjs, "!"); }
;


quantifier:
  ZZ_FORALL
  { 
    const char* fa = zztokenList.removeLast();
    if (folDbg >= 1) printf("FORALL ");
    zzformulaListObjs.push(new ListObj("FORALL"));
    zzformulaStr.append("FORALL ");
    delete [] fa;
  }
| ZZ_EXIST
  { 
    const char* ex = zztokenList.removeLast();
    if (folDbg >= 1) printf("EXIST "); 
    zzformulaListObjs.push(new ListObj("EXIST"));
    zzformulaStr.append("EXIST ");
    delete [] ex;
  }
;


//variables ',' quant_variable | quant_variable
variables:  
  { if (folDbg >= 2) printf("variables: variables\n"); }
  quant_variable ',' 
  {  
    zzconsumeToken(zztokenList,",");
    if (folDbg >= 1) printf(", "); 
    zzformulaStr.append(",");
  }
  variables
| 
  quant_variable  { zzformulaStr.append(" "); }
;


quant_variable: ZZ_VARIABLE
{
  const char* varName = zztokenList.removeLast();  
  if (folDbg >= 1) printf("v_%s ", varName); 

  //if (isupper(varName[0])) 
  if (zzisConstant(varName)) 
  {
    zzerr("Variable %s must be begin with a lowercase character.", varName);
    ((char*)varName)[0] = tolower(varName[0]);
  }

  int scopeNum = zzoldNewVarList.front().second;
  string newVarName = zzappend(varName,scopeNum);

  StringToStringMap* oldNewVarMap = zzoldNewVarList.front().first;
  (*oldNewVarMap)[varName]=newVarName;
  zzvarNameToIdMap[newVarName] = ZZVarIdType(--zzvarCounter, zzanyTypeId);
  zzformulaListObjs.top()->append(newVarName.c_str());
  zzformulaStr.append(varName); //use the old var name in the orig string
  delete [] varName;
}


// asterisk ZZ_PREDICATE '(' terms ')' | term internal_predicate_sign term
atomic_sentence:  
  asterisk ZZ_PREDICATE
  {
    const char* predName = zztokenList.removeLast();
    if (folDbg >= 1) printf("p_%s ", predName); 

    zzfdnumPreds++;
    ListObj* predlo = new ListObj;

	if (PredicateTemplate::isInternalPredicateTemplateName(predName))
	{
	  //zzinfixPredName is misused here to store internal pred. name
	  zzinfixPredName = (char *)malloc((strlen(predName)
    								  	+ 1)*sizeof(char));
	  strcpy(zzinfixPredName, predName);
	  const PredicateTemplate* t = zzdomain->getEmptyPredicateTemplate();
      zzassert(zzpred == NULL,"expecting zzpred==NULL");
      zzpred = new Predicate(t);
      predlo->append(PredicateTemplate::EMPTY_NAME);
	}
	else
	{
      zzcreatePred(zzpred, predName);
      predlo->append(predName);
	}
	
    zzformulaStr.append(predName);
    if(zzisNegated)  { zzpred->setSense(false); zzisNegated = false; }
    zzpredFuncListObjs.push(predlo);
    delete [] predName;
  }
  nothing_or_terms
| 
  // term internal_predicate_sign term
  {
    ++zzfdnumPreds;
    //zzfdisEqualPred = true;
    //const PredicateTemplate* t = zzdomain->getEqualPredicateTemplate();
    const PredicateTemplate* t = zzdomain->getEmptyPredicateTemplate();

    zzassert(zzpred == NULL, "expecting zzpred==NULL");
    zzpred = new Predicate(t);

    ListObj* predlo = new ListObj;
    //predlo->append(PredicateTemplate::EQUAL_NAME);
    predlo->append(PredicateTemplate::EMPTY_NAME);
    zzpredFuncListObjs.push(predlo);
    if (zzisNegated)  { zzpred->setSense(false); zzisNegated = false; }

    if (folDbg >= 2) printf("atomic_sentence (left): term\n"); 
  }
  term         
  internal_predicate_sign
  {
  	ListObj* predlo = zzpredFuncListObjs.top();
    predlo->replace(PredicateTemplate::EMPTY_NAME, zzinfixPredName);
    
  	  // If type known from LHS, then set the pred types accordingly
    int lTypeId = zzgetTypeId(zzpred->getTerm(0), (*predlo)[1]->getStr());
    if (lTypeId > 0)
    {
      if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
      {
        zzsetEqPredTypeName(lTypeId);
      }
      else
 	  {
 	    zzsetInternalPredTypeName(zzinfixPredName, lTypeId);
 	  }
    }
  }
  term
  {  
    ListObj* predlo = zzpredFuncListObjs.top();
    //predlo->replace(PredicateTemplate::EMPTY_NAME, zzinfixPredName);

	// types are possibly unknown
	// If '=' predicate then types are possibly unknown
	//if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0) {
      int lTypeId = zzgetTypeId(zzpred->getTerm(0), (*predlo)[1]->getStr());
      int rTypeId = zzgetTypeId(zzpred->getTerm(1), (*predlo)[2]->getStr());

      if (lTypeId > 0 && rTypeId > 0) //if both types are known
      {
        if (lTypeId != rTypeId)
          zzerr("The types on the left and right of '=' must match.");
        if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
        {
 	      zzsetEqPredTypeName(lTypeId);
 	    }
 	    else
 	    {
 	      zzsetInternalPredTypeName(zzinfixPredName, lTypeId);
 	    }
      }
      else  // if only one type is known
      if ( (lTypeId<=0 && rTypeId>0) || (lTypeId>0 && rTypeId<=0) )
      {
        int knownTypeId = (lTypeId>0) ? lTypeId : rTypeId;
        if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
        {
          zzsetEqPredTypeName(knownTypeId);
        }
        else
 	    {
 	      zzsetInternalPredTypeName(zzinfixPredName, knownTypeId);
 	    }
        const char* lvarName = (*predlo)[1]->getStr();
        const char* rvarName = (*predlo)[2]->getStr();
        const char* unknownVarName = (lTypeId>0) ?  rvarName : lvarName;
        zzvarNameToIdMap[unknownVarName].typeId_ = knownTypeId;
      }
      else // if both types are unknown
      {
          //both sides must be variables
        const char* lvarName = (*predlo)[1]->getStr();
        const char* rvarName = (*predlo)[2]->getStr();
        lTypeId = zzgetVarTypeId(lvarName);
        rTypeId = zzgetVarTypeId(rvarName);
	
        if (lTypeId > 0 && rTypeId > 0) //if both types are known
        {
          if (lTypeId != rTypeId)
            zzerr("The types of %s and %s on the left and right of "
                   "'=' must match.", lvarName, rvarName);
          if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
          {
          	zzsetEqPredTypeName(lTypeId);
          }
          else
 	      {
 	      	zzsetInternalPredTypeName(zzinfixPredName, lTypeId);
 	      }
        }
        else  // if only one type is known
        if ( (lTypeId <= 0 && rTypeId > 0) || (lTypeId > 0 && rTypeId <= 0) )
        {
          int knownTypeId = (lTypeId > 0) ? lTypeId : rTypeId;
          if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME) == 0)
          {
          	zzsetEqPredTypeName(knownTypeId);
          }
          else
 	      {
 	      	zzsetInternalPredTypeName(zzinfixPredName, knownTypeId);
 	      }
          const char* unknowVarName = (lTypeId>0) ?  rvarName : lvarName;
          zzvarNameToIdMap[unknowVarName].typeId_ = knownTypeId;
        }
        else
        {      
		  if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
          {
          	string unknownEqName = zzappend(PredicateTemplate::EQUAL_NAME, 
            	                            zzeqTypeCounter++);
          	zzeqPredList.push_back(ZZUnknownEqPredInfo(unknownEqName,lvarName,
                                                       rvarName));
          	predlo->replace(PredicateTemplate::EQUAL_NAME, unknownEqName.c_str());
          }
          else
          {
          	string unknownIntPredName =
          	  zzappendWithUnderscore(zzinfixPredName, zzintPredTypeCounter++);
          	zzintPredList.push_back(ZZUnknownIntPredInfo(unknownIntPredName, lvarName,
                                                       	 rvarName));
          	predlo->replace(zzinfixPredName, unknownIntPredName.c_str());
          }
        }
      }
	//}
	//else // Infix predicate other than '='
	//{
	  //Only left term could be unknown
	  //const char* leftTerm = (*predlo)[1]->getStr();
	//}

    zzassert(zzpredFuncListObjs.size()==1,
             "expecting zzpredFuncListObjs.size()==1");
    ListObj* topPredlo = zzpredFuncListObjs.top();
    zzpredFuncListObjs.pop();
    zzformulaListObjs.push(topPredlo);
      // This allows for the '!=' operator
    if (zzisNegated)
    {
      zzcreateListObjFromTop(zzformulaListObjs, "!");
      zzisNegated = false;
    }
    
	delete zzpred;
	zzpred = NULL;
	
	// If we have replaced a function inside the predicate
	// then we have to add the conjunction
	while (!zzfuncConjStack.empty())
	{
	  // create the second part of the conjunction
	  //zzformulaStr.append(" ^ ");
		
      ListObj* topPredlo = zzfuncConjStack.top();
	  zzfuncConjStack.pop();
      //zzformulaListObjs.push(topPredlo);
	  //zzcreateListObjFromTopTwo(zzformulaListObjs, "^");
	  zzformulaListObjs.push(topPredlo);
	  zzcreateListObjFromTopTwo(zzformulaListObjs, "v");
	} //while (!zzfuncConjStack.empty())

	if (!zzfuncConjStr.empty())
	{
	  zzformulaStr.append(zzfuncConjStr);
	  zzfuncConjStr.clear();
	}
	zzfunc = NULL;
	free(zzinfixPredName);
    zzinfixPredName = NULL;
  }
;

nothing_or_terms:
  // Nothing (propositional case)
  {
	// Add the one constant for propositional case
    const char* constName = Domain::PROPOSITIONAL_CONSTANT;
    if (folDbg >= 1) printf("c2_%s ", constName); 
    zztermIsConstant(constName, constName, false);

    zzcheckPredNumTerm(zzpred);
    delete zzpred;
    zzpred = NULL;
    if (!zzisHybrid)
    {
      zzassert(zzpredFuncListObjs.size()==1,
               "expecting zzpredFuncListObjs.size()==1");
      ListObj* predlo = zzpredFuncListObjs.top();
      zzpredFuncListObjs.pop();
      zzformulaListObjs.push(predlo);
    }
  }
|
  '(' 
  {  
    zzconsumeToken(zztokenList, "(");
    if (folDbg >= 1) printf("( "); 
    if (folDbg >= 2) printf("atomic_sentence: terms\n"); 
    zzformulaStr.append("(");
    //if (zzisHybrid) zzcontinuousStr.append("(");
  }
  terms 
  ')'
  {  
    zzconsumeToken(zztokenList, ")");
    if (folDbg >= 1) printf(") "); 

	  //If an internal pred., then need to determine type
	  //zzinfixPredName is misused here to store internal pred. name
	if (zzinfixPredName)
	{
	  ListObj* predlo = zzpredFuncListObjs.top();
      predlo->replace(PredicateTemplate::EMPTY_NAME, zzinfixPredName);
		// types are possibly unknown
		// If '=' predicate then types are possibly unknown
		//if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0) {
      int lTypeId = zzgetTypeId(zzpred->getTerm(0), (*predlo)[1]->getStr());
      int rTypeId = zzgetTypeId(zzpred->getTerm(1), (*predlo)[2]->getStr());

      if (lTypeId > 0 && rTypeId > 0) //if both types are known
      {
        if (lTypeId != rTypeId)
          zzerr("The types on the left and right of '=' must match.");
        if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
        {
 	      zzsetEqPredTypeName(lTypeId);
 	    }
 	    else
 	    {
 	      zzsetInternalPredTypeName(zzinfixPredName, lTypeId);
 	    }
      }
      else  // if only one type is known
      if ( (lTypeId<=0 && rTypeId>0) || (lTypeId>0 && rTypeId<=0) )
      {
        int knownTypeId = (lTypeId>0) ? lTypeId : rTypeId;
        if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
        {
          zzsetEqPredTypeName(knownTypeId);
        }
        else
 	    {
 	      zzsetInternalPredTypeName(zzinfixPredName, knownTypeId);
 	    }
        const char* lvarName = (*predlo)[1]->getStr();
        const char* rvarName = (*predlo)[2]->getStr();
        const char* unknownVarName = (lTypeId>0) ?  rvarName : lvarName;
        zzvarNameToIdMap[unknownVarName].typeId_ = knownTypeId;
      }
      else // if both types are unknown
      {
          //both sides must be variables
        const char* lvarName = (*predlo)[1]->getStr();
        const char* rvarName = (*predlo)[2]->getStr();
        lTypeId = zzgetVarTypeId(lvarName);
        rTypeId = zzgetVarTypeId(rvarName);
	
        if (lTypeId > 0 && rTypeId > 0) //if both types are known
        {
          if (lTypeId != rTypeId)
            zzerr("The types of %s and %s on the left and right of "
                   "'=' must match.", lvarName, rvarName);
          if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
          {
          	zzsetEqPredTypeName(lTypeId);
          }
          else
 	      {
 	      	zzsetInternalPredTypeName(zzinfixPredName, lTypeId);
 	      }
        }
        else  // if only one type is known
        if ( (lTypeId<=0 && rTypeId>0) || (lTypeId>0 && rTypeId<=0) )
        {
          int knownTypeId = (lTypeId>0) ? lTypeId : rTypeId;
          if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
          {
          	zzsetEqPredTypeName(knownTypeId);
          }
          else
 	      {
 	      	zzsetInternalPredTypeName(zzinfixPredName, knownTypeId);
 	      }
          const char* unknowVarName = (lTypeId>0) ?  rvarName : lvarName;
          zzvarNameToIdMap[unknowVarName].typeId_ = knownTypeId;
        }
        else
        {      
		  if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
          {
          	string unknownEqName = zzappend(PredicateTemplate::EQUAL_NAME, 
            	                            zzeqTypeCounter++);
          	zzeqPredList.push_back(ZZUnknownEqPredInfo(unknownEqName,lvarName,
                                                       rvarName));
          	predlo->replace(PredicateTemplate::EQUAL_NAME, unknownEqName.c_str());
          }
          else
          {
          	string unknownIntPredName =
          	  zzappendWithUnderscore(zzinfixPredName, zzintPredTypeCounter++);
          	zzintPredList.push_back(ZZUnknownIntPredInfo(unknownIntPredName, lvarName,
                                                       	 rvarName));
          	predlo->replace(zzinfixPredName, unknownIntPredName.c_str());
          }
        }
      }
	  free(zzinfixPredName);
      zzinfixPredName = NULL;
	}

    zzcheckPredNumTerm(zzpred);
    delete zzpred;
    zzpred = NULL;
    
    if (!zzisHybrid)
    {
      zzassert(zzpredFuncListObjs.size()==1,
               "expecting zzpredFuncListObjs.size()==1");
      ListObj* predlo = zzpredFuncListObjs.top();
      zzpredFuncListObjs.pop();

      if (zzisAsterisk)
      {
        zzisAsterisk = false;
        ListObj* lo = new ListObj;
        lo->append("*"); lo->append(predlo);
        zzformulaListObjs.push(lo);
      }
      else
        zzformulaListObjs.push(predlo);
	}
	else
	{
	  
	}
    zzformulaStr.append(")");

	// If we have replaced a function inside the predicate
	// then we have to add the conjunction
	while (!zzfuncConjStack.empty())
	{
		// create the second part of the conjunction
		//zzformulaStr.append(" ^ ");
      ListObj* topPredlo = zzfuncConjStack.top();
      zzfuncConjStack.pop();
      //zzformulaListObjs.push(topPredlo);
      //zzcreateListObjFromTopTwo(zzformulaListObjs, "^");
      zzformulaListObjs.push(topPredlo);
      zzcreateListObjFromTopTwo(zzformulaListObjs, "v");
	} //while (!zzfuncConjStack.empty())

	if (!zzfuncConjStr.empty())
	{
      zzformulaStr.append(zzfuncConjStr);
      zzfuncConjStr.clear();
	}
	zzfunc = NULL;
  }
;

internal_predicate_sign:
  '>'
  {
  	zzconsumeToken(zztokenList, ">");
    if (folDbg >= 1) printf("> "); 
    if (folDbg >= 2) printf("atomic_sentence (right): term\n"); 
    zzformulaStr.append(" > ");
    zzinfixPredName = (char *)malloc((strlen(PredicateTemplate::GT_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixPredName, PredicateTemplate::GT_NAME);
  	//zzcreateInternalPredTemplate(zzinfixPredName);
    //const PredicateTemplate* t = zzdomain->getPredicateTemplate(zzinfixPredName);
	//zzpred->setTemplate((PredicateTemplate*)t);
  }
|
  '<'
  {
   	zzconsumeToken(zztokenList, "<");
    if (folDbg >= 1) printf("< "); 
    if (folDbg >= 2) printf("atomic_sentence (right): term\n"); 
    zzformulaStr.append(" < "); 
    zzinfixPredName = (char *)malloc((strlen(PredicateTemplate::LT_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixPredName, PredicateTemplate::LT_NAME);
  	//zzcreateInternalPredTemplate(zzinfixPredName);
    //const PredicateTemplate* t = zzdomain->getPredicateTemplate(zzinfixPredName);
	//zzpred->setTemplate((PredicateTemplate*)t);
  }
|
  '>''='
  {
  	zzconsumeToken(zztokenList, ">");
    if (folDbg >= 1) printf(">");
    zzformulaStr.append(" >");
    zzconsumeToken(zztokenList, "=");
    if (folDbg >= 1) printf("= "); 
    if (folDbg >= 2) printf("atomic_sentence (right): term\n"); 
    zzformulaStr.append("= ");
    zzinfixPredName = (char *)malloc((strlen(PredicateTemplate::GTEQ_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixPredName, PredicateTemplate::GTEQ_NAME);
  	//zzcreateInternalPredTemplate(zzinfixPredName);
	//const PredicateTemplate* t = zzdomain->getPredicateTemplate(zzinfixPredName);
	//zzpred->setTemplate((PredicateTemplate*)t);
  }
|
  '<''='
  {
	zzconsumeToken(zztokenList, "<");
    if (folDbg >= 1) printf("<");
    zzformulaStr.append(" <");
    zzconsumeToken(zztokenList, "=");
    if (folDbg >= 1) printf("= "); 
    if (folDbg >= 2) printf("atomic_sentence (right): term\n"); 
    zzformulaStr.append("= ");
    zzinfixPredName = (char *)malloc((strlen(PredicateTemplate::LTEQ_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixPredName, PredicateTemplate::LTEQ_NAME);
  	//zzcreateInternalPredTemplate(zzinfixPredName);
    //const PredicateTemplate* t = zzdomain->getPredicateTemplate(zzinfixPredName);
	//zzpred->setTemplate((PredicateTemplate*)t);
  }
|
  '='
  {
  	zzconsumeToken(zztokenList, "=");
    if (folDbg >= 1) printf("= "); 
    if (folDbg >= 2) printf("atomic_sentence (right): term\n"); 
    zzformulaStr.append(" = ");
    zzinfixPredName = (char *)malloc((strlen(PredicateTemplate::EQUAL_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixPredName, PredicateTemplate::EQUAL_NAME);
    const PredicateTemplate* t = zzdomain->getEqualPredicateTemplate();
	zzpred->setTemplate((PredicateTemplate*)t);
  }
|
  '!''='
  {
  	zzconsumeToken(zztokenList, "!");
  	zzconsumeToken(zztokenList, "=");
    if (folDbg >= 1) printf("!= "); 
    if (folDbg >= 2) printf("atomic_sentence (right): term\n"); 
    zzformulaStr.append(" != ");
    zzinfixPredName = (char *)malloc((strlen(PredicateTemplate::EQUAL_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixPredName, PredicateTemplate::EQUAL_NAME);
    const PredicateTemplate* t = zzdomain->getEqualPredicateTemplate();
	zzpred->setTemplate((PredicateTemplate*)t);
	  // Theoretically, we could have "!(a1 != a2)" which would be "a1 = a2"
    if (zzpred->getSense())
    {
      zzpred->setSense(false);
      zzisNegated = true;
    }
    else { zzpred->setSense(true); }
  }
;


asterisk: 
  // empty
| '*'
  { 
    zzconsumeToken(zztokenList, "*");
    if (folDbg >= 1) printf("* "); 
    if (zzisNegated) zzerr("'!' and '*' cannot be used at the same time");
    zznumAsterisk++;
    zzassert(!zzisAsterisk,"expecting !zzisAsterisk");
    zzisAsterisk = true;
    zzformulaStr.append("*");
  }
;


// term ',' terms | term
terms:  
  terms  
  ',' 
  {  
    zzconsumeToken(zztokenList, ",");
    if (folDbg >= 1) printf(", "); 
    if (folDbg >= 2) printf("terms: term\n"); 
    // While parsing function, we do not want to append anything to the formula
    if (zzfunc == NULL) zzformulaStr.append(",");
  }
  term
| 
  term
  {
  	// After the first term in an internal pred., check if we can determine type
  	if (zzpred && zzpred->isEmptyPred() &&
  		zzpred->getNumTerms() == 1 && zzinfixPredName)
  	{
      ListObj* predlo = zzpredFuncListObjs.top();
      predlo->replace(PredicateTemplate::EMPTY_NAME, zzinfixPredName);
    
  	  	// If type known from term, then set the pred types accordingly
	  int lTypeId = zzgetTypeId(zzpred->getTerm(0), (*predlo)[1]->getStr());
      if (lTypeId>0)
      {
      	if (strcmp(zzinfixPredName, PredicateTemplate::EQUAL_NAME)==0)
      	{
          zzsetEqPredTypeName(lTypeId);
      	}
      	else
 	  	{
 	      zzsetInternalPredTypeName(zzinfixPredName, lTypeId);
 	  	}
      }
  	}
  }
;


term:
  ZZ_CONSTANT 
  {
    const char* constName = zztokenList.removeLast();
    if (folDbg >= 1) printf("c2_%s ", constName); 
    zztermIsConstant(constName, constName, true);
    if (zzfunc) zzfdfuncConstants.append(string(constName));
    else        zzfdconstName = constName;
    delete [] constName;
  }
|
  ZZ_STRING // the string is a constant
  {
    const char* constName = zztokenList.removeLast();
    if (folDbg >= 1) printf("c2_%s ", constName); 
      if (zzconstantMustBeDeclared)
        zzerr("Constant %s must be declared before it is used", constName);
    zztermIsConstant(constName, constName, true);
    if (zzfunc) zzfdfuncConstants.append(string(constName));
    else        zzfdconstName = constName;
    delete [] constName;
  }
|
  ZZ_NUM // the integer is a constant
  {
    const char* intStr = zztokenList.removeLast();
    if (folDbg >= 1) printf("c3_%s ", intStr);

    char constName[100];
    zzcreateAndCheckNumConstant(intStr, zzfunc, zzpred, zzdomain, constName);
    if (constName == NULL) { break; delete [] intStr; }

    zztermIsConstant(constName, intStr, true);
    if (zzfunc) zzfdfuncConstants.append(string(constName));
    else        zzfdconstName = constName;
    delete [] intStr;
  }
| 
  ZZ_VARIABLE
  {
    zztermIsVariable(folDbg);
    if (zzisPlus) zzisPlus = false;
  }
| 
  '+' ZZ_VARIABLE
  {
    zzconsumeToken(zztokenList, "+");
    if (folDbg >= 1) printf("+ "); 
    zzassert(!zzisPlus,"expecting !zzisPlus");
    zzisPlus = true;
    zzformulaStr.append("+");
    zztermIsVariable(folDbg);
    if (zzisPlus) zzisPlus = false;
  }  
|
  function_term
  {
    zzassert(zzfunc != NULL,"expecting zzfunc != NULL");
    zzcheckFuncNumTerm(zzfunc);
    zzassert(zzpred != NULL, "expecting zzpred != NULL");

      // replace function
	 
	  // Append Term funcVar<zzfuncVarCounter> to predicate where
      // function appears and set flag to add conjunction

      // 1. Make new variable
    char* varName;
    string newVarName;
    
      // Check if function mapping already occurs in formula
    ZZFuncToFuncVarMap::iterator mit;
    if ((mit = zzfuncToFuncVarMap.find(zzfunc)) != zzfuncToFuncVarMap.end())
    {
      newVarName = (*mit).second;
      varName = (char *)malloc((strlen(newVarName.c_str()) + 1) * sizeof(char));
      strcpy(varName, newVarName.c_str());
    }
    else
    {
      char funcVarCounterString[10];
      int funcVarCounterLength =
        sprintf(funcVarCounterString, "%d", zzfuncVarCounter);
      varName = (char *)malloc((strlen(ZZ_FUNCVAR_PREFIX) +
  								funcVarCounterLength + 1)*sizeof(char));
      strcpy(varName, ZZ_FUNCVAR_PREFIX);
      strcat(varName, funcVarCounterString);
      ++zzfdnumVars;
      ++zzfuncVarCounter;
      newVarName = zzgetNewVarName(varName);
      
      Function* func = new Function(*zzfunc);
      zzfuncToFuncVarMap[func] = newVarName;
    }

   	bool rightNumTerms = true;
   	bool rightType = true;
   	int exp, unexp;

      // 2. Create new predicate
      // Predicate name is PredicateTemplate::ZZ_RETURN_PREFIX + function name
	char* predName;
	predName = (char *)malloc((strlen(PredicateTemplate::ZZ_RETURN_PREFIX) +
  	    					   strlen(zzfunc->getName()) + 1)*sizeof(char));
	strcpy(predName, PredicateTemplate::ZZ_RETURN_PREFIX);
 	strcat(predName, zzfunc->getName());
    	
	// Only insert predicate declaration, if not yet declared
	if (zzdomain->getPredicateId(predName) < 0)
	{
      zzassert(zzpredTemplate==NULL,"expecting zzpredTemplate==NULL");
      zzpredTemplate = new PredicateTemplate();
      zzpredTemplate->setName(predName);
			
      // Register the types
      // First parameter is the return type
	  const char* ttype = zzfunc->getRetTypeName();
	  zzaddType(ttype, zzpredTemplate, NULL, false, zzdomain);

	  // Register the parameter types
	  for (int i = 0; i < zzfunc->getNumTerms(); i++)
	  {
		const char* ttype = zzfunc->getTermTypeAsStr(i);
		zzaddType(ttype, zzpredTemplate, NULL, false, zzdomain);
	  }
	
	  zzassert(zzpredTemplate, "not expecting zzpredTemplate==NULL");
	  int id = zzdomain->addPredicateTemplate(zzpredTemplate);
	  zzassert(id >= 0, "expecting pred template id >= 0");
	  zzpredTemplate->setId(id);
	  zzpredTemplate = NULL;
	}

	Predicate* prevPred = zzpred;
	zzpred = NULL;
    zzfdnumPreds++;
	zzcreatePred(zzpred, predName);
    
    ListObj* predlo = new ListObj;
    predlo->append(predName);
		
	// Put predicate list object in stack to be used in conjunction later
    zzfuncConjStack.push(predlo);
    //zzfuncConjStr.append(" ^ ");
    zzfuncConjStr.append(" v !");
    zzfuncConjStr.append(predName);
	zzfuncConjStr.append("(");
    zzputVariableInPred(varName, folDbg);
    zzfuncConjStack.top()->append(varName);
    zzfuncConjStr.append(varName);
      	
	// Information about the terms is in zzfunc
	for (int i = 0; i < zzfunc->getNumTerms(); i++)
	{
	  zzfuncConjStr.append(", ");
	  Term* term = new Term(*zzfunc->getTerm(i));
	  zzpred->appendTerm(term);

	  const char* name;
	  if (term->getType() == Term::VARIABLE)
	  {
		name = (*zzpredFuncListObjs.top())[i+1]->getStr();
	  }
	  else if (term->getType() == Term::CONSTANT)
	  {
		name = zzdomain->getConstantName(term->getId());
	  }
	  zzpredFuncListObjs.top()->append(name);
	  zzfuncConjStack.top()->append(name);
	  zzfuncConjStr.append(name);
	}
	zzfuncConjStr.append(")");
    zzcreateListObjFromTop(zzfuncConjStack, "!");
    
    // 4. Append new variable to function in stack or to predicate
    if (!zzfuncStack.empty())
    {
   	  Function* prevFunc = zzfuncStack.top(); 
	  zzfuncStack.pop();
    
      // check that we have not exceeded the number of terms
      if ((unexp = prevFunc->getNumTerms()) ==
       	  (exp = prevFunc->getTemplate()->getNumTerms()))
      {
       	rightNumTerms = false;
       	zzerr("Wrong number of terms for function %s. "
           	  "Expected %d but given %d", prevFunc->getName(), exp, unexp+1);
      }
      
      int varId = -1;      
      if (rightNumTerms)
      {
       	// check that the variable is of the right type
     	int typeId = prevFunc->getTermTypeAsInt(prevFunc->getNumTerms());
      	rightType = zzcheckRightTypeAndGetVarId(typeId, newVarName.c_str(),
	                                            varId);
      }

      if (rightNumTerms && rightType)
      {
     	prevFunc->appendTerm(new Term(varId, (void*)prevFunc, false));
     	zzpredFuncListObjs.pop();
    	zzpredFuncListObjs.top()->append(newVarName.c_str());
      }
    	      	  
      zzfunc = prevFunc;
    }
	else // function stack is empty, so append to predicate
    {
      // check that we have not exceeded the number of terms
      if ((unexp = prevPred->getNumTerms()) ==
      	  (exp = prevPred->getTemplate()->getNumTerms()))
      {
        rightNumTerms = false;
        zzerr("Wrong number of terms for predicate %s. "
              "Expected %d but given %d", prevPred->getName(), exp, unexp+1);
      }
        
      int varId = -1;
      if (rightNumTerms)
      {
          // check that the variable is of the right type
        int typeId = prevPred->getTermTypeAsInt(prevPred->getNumTerms());
        rightType = zzcheckRightTypeAndGetVarId(typeId, newVarName.c_str(),
        									    varId);
      }
      
	  if (rightNumTerms && rightType)
	  {
		prevPred->appendTerm(new Term(varId, (void*)prevPred, true));

   		// Pop the function from the stack
   		zzoldFuncLo = zzpredFuncListObjs.top();
   		zzpredFuncListObjs.pop();
   		zzpredFuncListObjs.top()->append(newVarName.c_str());
   		zzformulaStr.append(varName);
      }
	  zzfunc = NULL;
        	
    }
    free(varName);
	free(predName);
	//zzformulaStr.append(")");

    delete zzpred;
	zzpred = prevPred;
  }
;

  // ZZ_FUNCTION '(' terms ')' | term internal_function_sign term
function_term: 
  ZZ_FUNCTION
  {
    const char* funcName = zztokenList.removeLast();
    if (folDbg >= 1) printf("f_%s ", funcName);
    if (zzfunc != NULL) { zzfuncStack.push(zzfunc); zzfunc = NULL; }
    ++zzfdnumFuncs;
    zzfdfuncName = funcName;

    ListObj* funclo = new ListObj;

	  // Check if internal function
 	if (FunctionTemplate::isInternalFunctionTemplateName(funcName))
	{
	  zzinfixFuncName = (char *)malloc((strlen(funcName)
    								  	+ 1)*sizeof(char));
	  strcpy(zzinfixFuncName, funcName);
	  const FunctionTemplate* t;
	  if (FunctionTemplate::isInternalFunctionUnaryTemplateName(funcName))
	  	t = zzdomain->getEmptyFunctionUnaryTemplate();
	  else
	  	t = zzdomain->getEmptyFunctionBinaryTemplate();
      zzassert(zzfunc == NULL,"expecting zzfunc==NULL");
      zzfunc = new Function(t);
      funclo->append(FunctionTemplate::EMPTY_FTEMPLATE_NAME);
	}
	else
	{
	  zzcreateFunc(zzfunc, funcName); 
	  funclo->append(funcName);
	}
	
    zzpredFuncListObjs.push(funclo);
    //zzformulaStr.append(funcName);

    delete [] funcName;
  }
  '('
  {
    zzconsumeToken(zztokenList, "(");
    if (folDbg >= 1) printf("( "); 
    if (folDbg >= 2) printf("term: terms\n");
  }
  terms ')'
  {
  	  //If an internal func., then need to determine type
	if (zzinfixFuncName)
	{
	  ListObj* funclo = zzpredFuncListObjs.top();
      funclo->replace(FunctionTemplate::EMPTY_FTEMPLATE_NAME, zzinfixFuncName);
	  const FunctionTemplate* t =
	    zzgetGenericInternalFunctionTemplate(zzinfixFuncName);
		// If in a pred. (not LHS of infix pred.), then we know the return type
	  if (zzpred)
	  {
	  	((FunctionTemplate*)t)->setRetTypeId(
	  		zzpred->getTermTypeAsInt(zzpred->getNumTerms()), zzdomain);
      }
	  zzfunc->setTemplate((FunctionTemplate*)t);

		// types are possibly unknown
	  bool unknownTypes = false;
		
		// First element is return type
	  Array<int> typeIds(zzfunc->getNumTerms() + 1);
	  typeIds.append(zzfunc->getRetTypeId());
	  if (typeIds[0] <= 0) unknownTypes = true;
		
		// Then term types
	  for (int i = 1; i <= zzfunc->getNumTerms(); i++)
	  {
	  	typeIds.append(zzgetTypeId(zzfunc->getTerm(i-1), (*funclo)[i]->getStr()));
	  	if (typeIds[i] <= 0) unknownTypes = true;
	  }

		// If all types are known
	  if (!unknownTypes)
	  {
 	  	zzsetInternalFuncTypeName(zzinfixFuncName, typeIds);
	  }
      else // Not all types are known, delay processing
      {
      	Array<string> varNames(typeIds.size());
      	
      	char* varName;
    	char funcVarCounterString[10];
    	int funcVarCounterLength =
      		sprintf(funcVarCounterString, "%d", zzfuncVarCounter);
    	varName = (char *)malloc((strlen(ZZ_FUNCVAR_PREFIX) +
  								  funcVarCounterLength + 1)*sizeof(char));
  		strcpy(varName, ZZ_FUNCVAR_PREFIX);
 		strcat(varName, funcVarCounterString);
	    string newVarName = zzgetNewVarName(varName);
      	
      	varNames.append(newVarName);
      	
      	for (int i = 1; i < typeIds.size(); i++)
      	  varNames[i] = (*funclo)[i]->getStr();

		// Predicate name is PredicateTemplate::ZZ_RETURN_PREFIX + function name
		char* predName;
		predName = (char *)malloc((strlen(PredicateTemplate::ZZ_RETURN_PREFIX) +
  	    					   strlen(zzfunc->getName()) + 1)*sizeof(char));
		strcpy(predName, PredicateTemplate::ZZ_RETURN_PREFIX);
 		strcat(predName, zzinfixFuncName);
 	
        string unknownIntFuncName =
          	  zzappendWithUnderscore(predName, zzintFuncTypeCounter++);
        zzintFuncList.push_back(ZZUnknownIntFuncInfo(unknownIntFuncName, varNames));
        //funclo->replace(zzinfixFuncName, unknownIntFuncName.c_str());
		free(predName);
		free(varName);
      }
	  free(zzinfixFuncName);
      zzinfixFuncName = NULL;
	}
  	
    zzconsumeToken(zztokenList, ")");
    if (folDbg >= 1) printf(") ");
  }
|
  term internal_function_sign
  {
      // Create a function corresponding to the infix sign
    if (zzfunc != NULL) { zzfuncStack.push(zzfunc); zzfunc = NULL; }
    ++zzfdnumFuncs;
    zzfdfuncName = zzinfixFuncName;
	const FunctionTemplate* t =
      zzgetGenericInternalFunctionTemplate(zzinfixFuncName);
  	zzfunc = new Function(t);
	  // If in a pred. (not LHS of infix pred.), then we know the return type
      // It's the type of the previous variable (getNumTerms() - 1)
	if (zzpred)
      ((FunctionTemplate*)t)->setRetTypeId(
        zzpred->getTermTypeAsInt(zzpred->getNumTerms() - 1), zzdomain);

    ListObj* funclo = new ListObj;
    funclo->append(zzfdfuncName.c_str());
    zzpredFuncListObjs.push(funclo);

      // The first term was added to the previous pred / func
      // So this needs to be switched to the infix function
    moveTermToInfixFunction();
  }
  term
  {
    ListObj* funclo = zzpredFuncListObjs.top();

	  // types are possibly unknown
    bool unknownTypes = false;
	
      // First element is return type
    Array<int> typeIds(zzfunc->getNumTerms() + 1);
    typeIds.append(zzfunc->getRetTypeId());
    if (typeIds[0] <= 0) unknownTypes = true;
		
      // Then term types
	for (int i = 1; i <= zzfunc->getNumTerms(); i++)
	{
      typeIds.append(zzgetTypeId(zzfunc->getTerm(i-1), (*funclo)[i]->getStr()));
      if (typeIds[i] <= 0) unknownTypes = true;
	}

	  // If all types are known
	if (!unknownTypes)
	{
 	  zzsetInternalFuncTypeName(zzinfixFuncName, typeIds);
	}
    else // Not all types are known, delay processing
    {
      Array<string> varNames(typeIds.size());
      	
      char* varName;
      char funcVarCounterString[10];
      int funcVarCounterLength =
      		sprintf(funcVarCounterString, "%d", zzfuncVarCounter);
      varName = (char *)malloc((strlen(ZZ_FUNCVAR_PREFIX) +
  								funcVarCounterLength + 1)*sizeof(char));
  	  strcpy(varName, ZZ_FUNCVAR_PREFIX);
 	  strcat(varName, funcVarCounterString);
	  string newVarName = zzgetNewVarName(varName);
      	
      varNames.append(newVarName);
      	
      for (int i = 1; i < typeIds.size(); i++)
        varNames[i] = (*funclo)[i]->getStr();
		// Predicate name is PredicateTemplate::ZZ_RETURN_PREFIX + function name
      char* predName;
      predName = (char *)malloc((strlen(PredicateTemplate::ZZ_RETURN_PREFIX) +
  	       					     strlen(zzfunc->getName()) + 1)*sizeof(char));
      strcpy(predName, PredicateTemplate::ZZ_RETURN_PREFIX);
      strcat(predName, zzinfixFuncName);
 	
      string unknownIntFuncName =
       	zzappendWithUnderscore(predName, zzintFuncTypeCounter++);
      zzintFuncList.push_back(ZZUnknownIntFuncInfo(unknownIntFuncName, varNames));
        //funclo->replace(zzinfixFuncName, unknownIntFuncName.c_str());
      free(predName);
      free(varName);
    }
    free(zzinfixFuncName);
    zzinfixFuncName = NULL;
  }
;

internal_function_sign:
  '+'
  {
    zzconsumeToken(zztokenList, "+");
    if (folDbg >= 1) printf("+ "); 
    zzinfixFuncName = (char *)malloc((strlen(PredicateTemplate::PLUS_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixFuncName, PredicateTemplate::PLUS_NAME);
  }
|
  '-'
  {
   	zzconsumeToken(zztokenList, "-");
    if (folDbg >= 1) printf("- "); 
    zzinfixFuncName = (char *)malloc((strlen(PredicateTemplate::MINUS_NAME)
                                      + 1)*sizeof(char));
    strcpy(zzinfixFuncName, PredicateTemplate::MINUS_NAME);
  }
|
  '*'
  {
  	zzconsumeToken(zztokenList, "*");
    if (folDbg >= 1) printf("* ");
    zzinfixFuncName = (char *)malloc((strlen(PredicateTemplate::TIMES_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixFuncName, PredicateTemplate::TIMES_NAME);
  }
|
  '/'
  {
	zzconsumeToken(zztokenList, "/");
    if (folDbg >= 1) printf("/ ");
    zzinfixFuncName = (char *)malloc((strlen(PredicateTemplate::DIVIDEDBY_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixFuncName, PredicateTemplate::DIVIDEDBY_NAME);
  }
|
  '%'
  {
  	zzconsumeToken(zztokenList, "%");
    if (folDbg >= 1) printf("%% ");
    zzinfixFuncName = (char *)malloc((strlen(PredicateTemplate::MOD_NAME)
    								  + 1)*sizeof(char));
  	strcpy(zzinfixFuncName, PredicateTemplate::MOD_NAME);
  }
;

numeric_term:
  '['
  { 
      // Square Brackets
    zzconsumeToken(zztokenList,"["); 
    if (folDbg >= 1) printf("[ "); 
    if (folDbg >= 2) printf("numeric_term: '[' numeric_term\n");
    zzformulaStr.append("[");
  }
  numeric_term
  ']'  
  { 
    zzconsumeToken(zztokenList,"]"); 
    if (folDbg >= 1) printf("] "); 
    zzformulaStr.append("]");
  }
| 
  '('  
  { 
    zzconsumeToken(zztokenList,"("); 
    if (folDbg >= 1) printf("( "); 
    if (folDbg >= 2) printf("numeric_term: '(' numeric_term\n");
    zzformulaStr.append("(");
  }
  numeric_term
  ')'
  { 
    zzconsumeToken(zztokenList,")"); 
    if (folDbg >= 1) printf(") "); 
    zzformulaStr.append(")");
  }
|
  gaussian_shorthand
|
  polynomial
;

gaussian_shorthand:
  ZZ_PREDICATE
  {
    const char* predName = zztokenList.removeLast();
    if (folDbg >= 1) printf("p_%s ", predName); 

    ListObj* predlo = new ListObj;

	if (PredicateTemplate::isInternalPredicateTemplateName(predName))
	{
	  //zzinfixPredName is misused here to store internal pred. name
	  zzinfixPredName = (char *)malloc((strlen(predName)
    								  	+ 1)*sizeof(char));
	  strcpy(zzinfixPredName, predName);
	  const PredicateTemplate* t = zzdomain->getEmptyPredicateTemplate();
      zzassert(zzpred == NULL,"expecting zzpred==NULL");
      zzpred = new Predicate(t);
      predlo->append(PredicateTemplate::EMPTY_NAME);
	}
	else
	{
      zzcreatePred(zzpred, predName);
      predlo->append(predName);
	}
	
    zzformulaStr.append(predName);
    if(zzisNegated)  { zzpred->setSense(false); zzisNegated = false; }
    zzcontPred = predlo;
    delete [] predName;
  }
  nothing_or_terms
  '='  
   {
    zzconsumeToken(zztokenList,"="); 
    if (folDbg >= 1) printf("= "); 
    zzformulaStr.append(" = ");
    //if (zzisHybrid) zzcontinuousStr.append("=");
  }
  ZZ_NUM
  {
  	const char* tmp = zztokenList.removeLast();
  	//zztmpReturnConstant = (char *)malloc((strlen(tmp) + 1)*sizeof(char));
  	//strcpy(zztmpReturnConstant, tmp);
  	//zztmpReturnNum = true;
    zzformulaStr.append(tmp);
    if (zzisHybrid) zzmean = atof(tmp);
  	if (folDbg >= 1) printf("icnum_%s ", tmp);
  	delete []tmp;
  }
|
  '-'  
  {
    zzconsumeToken(zztokenList,"-");
    if (folDbg >= 1) printf("-"); 
    zzformulaStr.append("-");
  }
  '('  
  {
    zzconsumeToken(zztokenList,"("); 
    if (folDbg >= 1) printf("("); 
    zzformulaStr.append("(");
  }
  ZZ_PREDICATE
  {
    const char* predName = zztokenList.removeLast();
    if (folDbg >= 1) printf("p_%s ", predName); 

    ListObj* predlo = new ListObj;

	if (PredicateTemplate::isInternalPredicateTemplateName(predName))
	{
	  //zzinfixPredName is misused here to store internal pred. name
	  zzinfixPredName = (char *)malloc((strlen(predName)
    								  	+ 1)*sizeof(char));
	  strcpy(zzinfixPredName, predName);
	  const PredicateTemplate* t = zzdomain->getEmptyPredicateTemplate();
      zzassert(zzpred == NULL,"expecting zzpred==NULL");
      zzpred = new Predicate(t);
      predlo->append(PredicateTemplate::EMPTY_NAME);
	}
	else
	{
      zzcreatePred(zzpred, predName);
      predlo->append(predName);
	}
	
    zzformulaStr.append(predName);
    if(zzisNegated)  { zzpred->setSense(false); zzisNegated = false; }
    zzcontPred = predlo;
    delete [] predName;
  }
  nothing_or_terms
  '-'
  { 
    zzconsumeToken(zztokenList,"-"); 
    if (folDbg >= 1) printf("- "); 
    zzformulaStr.append(" - ");
  }
  ZZ_NUM
  {
  	const char* tmp = zztokenList.removeLast();
  	//zztmpReturnConstant = (char *)malloc((strlen(tmp) + 1)*sizeof(char));
  	//strcpy(zztmpReturnConstant, tmp);
  	//zztmpReturnNum = true;
    zzformulaStr.append(tmp);
    if (zzisHybrid) zzmean = atof(tmp);
  	if (folDbg >= 1) printf("icnum_%s ", tmp);
  	delete []tmp;
  }
  ')'  
  {
    zzconsumeToken(zztokenList,")");
    if (folDbg >= 1) printf(")"); 
    zzformulaStr.append(")");
  }
  '^'  
  {
    zzconsumeToken(zztokenList,"^"); 
    if (folDbg >= 1) printf("^"); 
    zzformulaStr.append("^");
  }
  ZZ_NUM // must be the number 2
  {
    zzconsumeToken(zztokenList,"2"); 
    if (folDbg >= 1) printf("2"); 
    zzformulaStr.append("2");
  }
;

polynomial:
  'x^2 + x + 3'
  {
  }
;

%%

/******************* function definitions ****************************/

  //warnDuplicates set to true by default
bool runYYParser(MLN* const & mln, Domain* const & dom,
                 const char* const & fileName, 
                 const bool& allPredsExceptQueriesAreClosedWorld,
                 const StringHashArray* const & openWorldPredNames,
                 const StringHashArray* const & closedWorldPredNames,
                 const StringHashArray* const & queryPredNames,
                 const bool& addUnitClauses, const bool& warnDuplicates,
                 const double& defaultWt, const bool& mustHaveWtOrFullStop,
                 const Domain* const & domain0, const bool& lazyInference,
                 const bool& flipWtsOfFlippedClause)
{
  zzinit();
  if (fileName) { yyin = fopen(fileName, "r" ); zzinFileName = fileName; }
  else            yyin = stdin;
  if (yyin == NULL) zzexit("Failed to open file %s.", fileName);

  signal(SIGFPE,zzsigHandler);
  signal(SIGABRT,zzsigHandler);
  signal(SIGILL,zzsigHandler);
  signal(SIGSEGV,zzsigHandler);

  zzwarnDuplicates = warnDuplicates;
  zzdefaultWt = defaultWt;
  zzmustHaveWtOrFullStop = mustHaveWtOrFullStop;
  zzflipWtsOfFlippedClause = flipWtsOfFlippedClause;

  ungetc('\n', yyin); // pretend that file begin with a newline
  zzmln = mln;
  zzdomain = dom;

  zzanyTypeId = zzdomain->getTypeId(PredicateTemplate::ANY_TYPE_NAME);
  zzassert(zzanyTypeId >= 0, "expecting zzanyTypeId >= 0");

  // Declare internally implemented predicates and functions
  //declareInternalPredicatesAndFunctions();

  zzisParsing = true;
  yyparse();
  zzisParsing = false;

  zzcheckAllTypesHaveConstants(zzdomain);
  zzcreateBlocks(zzdomain);
  
  // Insert groundings generated from internally implemented functions and predicates
  zzgenerateGroundingsFromInternalPredicatesAndFunctions();

  // Insert groundings generated from linked-in predicates
  if (zzusingLinkedPredicates) zzgenerateGroundingsFromLinkedPredicates(zzdomain);

  // Insert groundings generated from linked-in functions
  if (zzusingLinkedFunctions) zzgenerateGroundingsFromLinkedFunctions(zzdomain);

  if (zzok)
  {
      //append the formulas to MLN and set the weights of the hard clauses
    cout << "Adding clauses to MLN..." << endl;
    zzappendFormulasToMLN(zzformulaInfos, mln, domain0);

    if (addUnitClauses) zzappendUnitClausesToMLN(zzdomain, zzmln, zzdefaultWt);

    zzmln->setClauseInfoPriorMeansToClauseWts();

      // Change the constant ids so that constants of same type are ordered 
      // consecutively. Also ensure that the constants in the mln and map is 
      // consistent with the new constant ids
    zzdomain->reorderConstants(zzmln, zzpredIdToGndPredMap);

    zzmln->compress();

    Array<bool> isClosedWorldArr; 
    zzfillClosedWorldArray(isClosedWorldArr,allPredsExceptQueriesAreClosedWorld,
                           openWorldPredNames, closedWorldPredNames,
                           queryPredNames);

    Database* db = new Database(zzdomain, isClosedWorldArr, zzstoreGroundPreds);
    if (lazyInference) db->setLazyFlag();
    zzaddGndPredsToDb(db);
    zzdomain->setDB(db);
    zzdomain->compress();
  }
  
  if (zznumErrors > 0) cout << "Num of errors detected = " << zznumErrors<<endl;

  zzcleanUp();

  signal(SIGFPE,SIG_DFL);
  signal(SIGABRT,SIG_DFL);
  signal(SIGILL,SIG_DFL);
  signal(SIGSEGV,SIG_DFL);

  return zzok;
}

