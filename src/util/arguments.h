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
// *************************************************************************
// * A general mechanism to parse command line and file arguments for 
// * C++ programs.
// * <bilmes@media.mit.edu> Aug 1992
// *
// *************************************************************************
//

// TO DO:
//    1. Usage message should print required arguments first.  
//    2. Parse from file should have capability to read multi word string args.
//    3. Make error reports better.
//    4. Fix bug with empty line in parse args from file (libg++ 2.2).
//    5. negative action, action called if argument is not present.
//    6. Add ability to specify controls in documentation string.
//       i.e. to we can do:  ... "Parm X, default = %d",default);

#ifndef ARGS_h
#define ARGS_h

//#include "common.h"
  // ... OR ...
  // The following bool definition will work fine.
#include <iostream>
using namespace std;
// enum bool {true = 1, false = 0};
// inline ostream& operator << (ostream&os,bool b) 
// {
//   return (os << ((b == true) ? "T" : "F"));
// }

class ARGS;
class argsAction;
extern ostream& operator << (ostream& os, ARGS*);
extern ostream& operator << (ostream& os, ARGS&);


class unionClass 
{
  friend class ARGS;
  friend class argsAction;
  friend ostream& operator << (ostream&, ARGS*);
  friend ostream& operator << (ostream&, ARGS&);
 
 public:
  enum argType 
  { bool_t, char_t, str_t, int_t, uint_t, float_t, double_t, action_t };

 private:
  union argU 
  {
    bool         boolean;
    char         ch;
    char*        string;
    int          integer;
    unsigned int uinteger;
    float        single_prec;
    double       double_prec;
    argsAction*  action;  // dummy field.
  };

  argU* ptr;
  argType type;
  static const char* typeStr(unionClass::argType);
  friend ostream& operator << (ostream&, unionClass::argType&);
 
 public:
  unionClass(bool& b)          { ptr = (argU*)&b; type = bool_t; }
  unionClass(char& c)          { ptr = (argU*)&c; type = char_t; }
  unionClass(char*& s)         { ptr = (argU*)&s; type = str_t; }
  unionClass(int& i)           { ptr = (argU*)&i; type = int_t; }
  unionClass(unsigned int& i)  { ptr = (argU*)&i; type = uint_t; }
  unionClass(float& f)         { ptr = (argU*)&f; type = float_t; }
  unionClass(double& d)        { ptr = (argU*)&d; type = double_t; }
  unionClass(argsAction& a)    { ptr = (argU*)&a; type = action_t; }
};


extern ostream& operator << (ostream&, unionClass::argType&);

class argsAction 
{
  enum ActionRetCode { OK, ERR };
 public:
  virtual ~argsAction() {}
 private:
  ActionRetCode setArg(char*, unionClass::argType, unionClass::argU);

 protected:
  ActionRetCode setArg(char*, bool);
  ActionRetCode setArg(char*, char);
  ActionRetCode setArg(char*, char*);
  ActionRetCode setArg(char*, int);
  ActionRetCode setArg(char*, float);
  ActionRetCode setArg(char*, double);
  ActionRetCode satisfy(char*);

 public:
  virtual void act()           {}
  virtual void act(bool& b)    { b = b;}
  virtual void act(char& c)    { c = c; }
  virtual void act(char*& str) { str = str; }
  virtual void act(int& i)     { i = i; }
  virtual void act(unsigned int& ui) { ui = ui; }
  virtual void act(float& f)   { f = f; }
  virtual void act(double& d)  { d = d; }
};



class ARGS 
{
  friend class unionClass;
  friend class argsAction;

 public:

  enum argKind     { Opt, Req, Tog }; // optional, required, or toggle
  enum ArgsRetCode { ARG_MISSING, ARG_OK, ARG_ERR };
  static ARGS Args[];
  static ARGS END;

 private:

  static char* const NOFLAG;
  static char* const NOFL_FOUND;
  static const char CommentChar;  // for argument files.
  static bool usage_called;
  static int numArgs;
  static bool* found;
  static char* progName;
  static bool ignoreUnknownSwitch;

    // instance data members.
  const char* flag;        // name to match on command line, NULL when end.
  argKind arg_kind;  // optional, required, toggle
  unionClass uc;
  const char* description;
  argsAction* action;

  
  static argsAction no_action;
  static bool __args__dummy__;

  static ostream* pout; //used to print parameters

 public:
  ARGS(const char*,argKind, unionClass, const char*d=NULL,
       argsAction&a=no_action);
  ARGS(argKind, unionClass, const char*d=NULL, argsAction&a=no_action);
  ARGS(unionClass, const char*d=NULL, argsAction&a=no_action);
  ARGS(const ARGS&);
  ARGS();

  static void cleanUp();
  static ArgsRetCode parseFromCommandLine(int,char**);
  static ArgsRetCode parseFromFile(char*f="argsFile");
  static void parse(int, char**, char*&, ostream* prout=NULL);
  static void parse(int i, char**c, ostream* prout=NULL) 
  { char*d=NULL; parse(i,c,d, prout); }
  static void printMissing()        { checkMissing(true); }
  static void usage(bool force=false);
  static void ignoreUnknownFlag()   { ignoreUnknownSwitch = true; }
  static void regardUnknownFlag()   { ignoreUnknownSwitch = false; }

 private:
  void init(const char*, argKind, const char*, argsAction&);
  static void Error(const char*, const char*s2=NULL, const char*s3=NULL,
                    const char*s4=NULL, const char*s5=NULL);

  static bool noFlagP(const char *);
  static ArgsRetCode argsSwitch(ARGS*, char*, int&, bool&, const char*);
  static ARGS* findMatch(ARGS*, char*);
  static void makeFound();
  static bool checkMissing(bool printMessage=false);
  static bool boolable(char* string, bool&value);
  friend ostream& operator << (ostream&,ARGS*);
  friend ostream& operator << (ostream&,ARGS&);

};

#endif
