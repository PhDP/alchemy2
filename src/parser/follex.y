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
/* scanner for first-order language */

%{

#include "folhelper.h"

bool follexDbg = false; 
//bool follexDbg = true; 

%}


/* not scanning another file after current one */
%option noyywrap 

ZZ_SINGLE_LINE_COMMENT [/][/][^\r\n]*
ZZ_MULTI_LINE_COMMENT [/][*][^*]*[*]([*]|[^*/][^*]*[*])*[/]


/* double ' and " to avoid string syntax highlighting in xemacs */
 /*ZZ_STRING  [""]([^""\\n\r] | ([\][ntbr\f''""]))*[""]*/
ZZ_STRING  [""]([^""\\\n\r]|([\\][ntbr\f''""]))*[""]
ZZ_NOT [!]
ZZ_OR  [ ]+[v][ \r\n]+
ZZ_AND [ ]+[\^][ \r\n]+
ZZ_IMPLY [=][>]
ZZ_EQUIV [<][=][>]
ZZ_EXIST  [Ee][Xx][Ii][Ss][Tt]
ZZ_FORALL [Ff][Oo][Rr][Aa][Ll][Ll]
ZZ_ASTERISK [*]
ZZ_PLUS [+]
ZZ_MINUS [-]
ZZ_MINUS_OR_PLUS [+-]
ZZ_QS [?]
ZZ_EQ  [=]
ZZ_DOTDOTDOT [.][.][.]

ZZ_INCLUDE [#][i][n][c][l][u][d][e]

ZZ_DIGIT [0-9] 
/* double ' to avoid string syntax highlighting in xemacs */
/* ZZ_ID [a-zA-z_\-][a-zA-Z0-9_\-'']* */
ZZ_ID [a-zA-Z_\-][a-zA-Z0-9_\-'']*

%%

"{" {
  if (follexDbg) printf("LBRACE: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


"}" {
  if (follexDbg) printf("RBRACE: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


{ZZ_SINGLE_LINE_COMMENT} { 
  if (follexDbg) printf("SINGLE_LINE_COMMENT: %s\n", yytext);
  if (yytext[strlen(yytext)-1] == '\n' || yytext[strlen(yytext)-1] == '\r') 
    {
      ++zzline;
      zzcolumn = -1;
    } 
    else
      zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
}

{ZZ_MULTI_LINE_COMMENT} { 
  if (follexDbg) printf("MULTI_LINE_COMMENT: %s\n", yytext);
  int i=-1;
  while (yytext[++i] != '\0')
  {
    if (yytext[i] == '\n' || yytext[i] == '\r') 
    {
      ++zzline;
      zzcolumn = -1;
    } 
    else
      zzcolumn++;
  }
  zznumCharRead += strlen(yytext); 
}


{ZZ_INCLUDE} { 
  if (follexDbg) printf("INCLUDE: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return ZZ_INCLUDE;
}


{ZZ_NOT} {
  if (follexDbg) printf("NOT: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


{ZZ_OR} {
  if (follexDbg) printf("OR: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add("v");
  return 'v';
}


{ZZ_AND} {
  if (follexDbg) printf("AND: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add("^");
  return '^';
}


{ZZ_IMPLY} {
  if (follexDbg) printf("IMPLY: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return ZZ_IMPLY;
}


{ZZ_EQUIV} {
  if (follexDbg) printf("EQUIV: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return ZZ_EQUIV;
}


{ZZ_EXIST} { 
  if (follexDbg) printf("EXIST: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return ZZ_EXIST;
}


{ZZ_FORALL} { 
  if (follexDbg) printf("FORALL: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return ZZ_FORALL;
}


{ZZ_ASTERISK} {
  if (follexDbg) printf("ASTERISK: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


{ZZ_PLUS} {
  if (follexDbg) printf("PLUS: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


{ZZ_MINUS} {
  if (follexDbg) printf("MINUS: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}

{ZZ_QS} {
  if (follexDbg) printf("QS: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


{ZZ_MINUS}?{ZZ_DIGIT}+ {
  if (follexDbg) printf("INTEGER: %s (%d)\n", yytext, atoi(yytext));
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return ZZ_NUM;  
}


{ZZ_MINUS}?{ZZ_DIGIT}+"."{ZZ_DIGIT}+ {
  if (follexDbg) printf("FLOAT: %s (%g)\n", yytext, atof(yytext));
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return ZZ_NUM;
}


{ZZ_MINUS}?{ZZ_DIGIT}+"."?{ZZ_DIGIT}*"e"{ZZ_MINUS_OR_PLUS}?{ZZ_DIGIT}+ {
  if (follexDbg) printf("EXP_FLOAT: %s (%e)\n", yytext, atof(yytext));
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return ZZ_NUM;
}


{ZZ_STRING} {
  if (follexDbg) printf("STRING: %s \n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  if (zzdomain->isConstant(yytext)) return ZZ_CONSTANT;
  return ZZ_STRING;
}


{ZZ_ID} {
  //if (follexDbg) printf("IDENTIFIER: %s\n", yytext );
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  if (zzdomain->isType(yytext))
  { 
    if (follexDbg) printf("ZZ_TYPE: %s\n", yytext ); 
    return ZZ_TYPE;
  }
  if (zzdomain->isPredicate(yytext)) 
  { 
    if (follexDbg) printf("ZZ_PREDICATE: %s\n", yytext ); 
    return ZZ_PREDICATE; 
  }
  if (zzdomain->isFunction(yytext))  
  {  
    if (follexDbg) printf("ZZ_FUNCTION: %s\n", yytext ); 
    return ZZ_FUNCTION;
  }
  if (zzdomain->isConstant(yytext))
  {
    if (follexDbg) printf("ZZ_CONSTANT: %s\n", yytext ); 
    return ZZ_CONSTANT;
  }
  if (follexDbg) printf("ZZ_VARIABLE: %s\n", yytext ); 
  return ZZ_VARIABLE;
}


"\n"|"\r" {
  if (follexDbg) 
  {
    if (zzparseGroundPred) printf("AT: %c\n", '@');
    else                   printf("NEWLINE: %s", yytext);
  }
  ++zzline;
  zznumCharRead += 1;
  zzcolumn = -1;
  if (zzparseGroundPred) { 
  	zzafterRtParen = false;
  	zztokenList.add("@");
  	return '@';
  }
  zztokenList.add(yytext);
  return yytext[0];
}

"[" {
  if (follexDbg) printf("LEFT BRACKET: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  //if (zzparseGroundPred) zzafterRtParen = false;
  return yytext[0];
}


"]" {
  if (follexDbg) printf("RIGHT BRACKET: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  //if (zzparseGroundPred) zzafterRtParen = true;
  return yytext[0];
}


"(" {
  if (follexDbg) printf("LEFT PAREN: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  if (zzparseGroundPred) zzafterRtParen = false;
  return yytext[0];
}


")" {
  if (follexDbg) printf("RIGHT PAREN: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  if (zzparseGroundPred) zzafterRtParen = true;
  return yytext[0];
}


"," {
  if (follexDbg) printf("COMMA: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


{ZZ_EQ} {
  if (follexDbg) printf("EQUAL: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add("=");
  return '=';
}


{ZZ_DOTDOTDOT} {
  if (follexDbg) printf("DOTDOTDOT: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add("...");
  return ZZ_DOTDOTDOT;
}


"@" {
  if (follexDbg) printf("AT: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


"." {
  if (follexDbg) printf("FULLSTOP: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}

":" {
  if (follexDbg) printf("COLON: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];
}


[ \t]+ { /* eat up whitespace */
  if (follexDbg) 
  {
    if (zzparseGroundPred) printf("WS:\n");
  }
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);

  //if (zzparseGroundPred && zzafterRtParen) { zztokenList.add("@"); return '@'; }
}


. {
    // commented out so as to allow fol.y to deal with the problem
  //printf("\nERROR: Unrecognized character. %s: ln %d, col %d: %s \n", 
  //       zzinFileName.c_str(), zzline, zzcolumn+1, yytext);
  //exit(-1);
  if (follexDbg) printf("Unrecognized character: %s\n", yytext);
  zzcolumn += strlen(yytext);
  zznumCharRead += strlen(yytext);
  zztokenList.add(yytext);
  return yytext[0];

}

<<EOF>> {

  if (zzparseGroundPred) 
  {
    if (zznumEofSeen == 0) 
    {
      zznumEofSeen++;
      if (follexDbg) printf("EOF returned as @:\n");
      zztokenList.add("@"); 
      return '@'; 
    }
    else
    {
      zzparseGroundPred = false;
      zzafterRtParen = false;
      zznumEofSeen = 0;
    }
  }

  if (zznumEofSeen==0)
  {
    ungetc(EOF, yyin);
    zznumEofSeen = 1;
    if (follexDbg) printf("EOF returned as \\n:\n");    
    zztokenList.add("\n");  // pretend that file ends with newline
    return '\n';
  }
  
  zznumEofSeen = 0;
  
  if (follexDbg) printf("EOF %s:\n", yytext);
  if (false) yyunput(1,NULL); //avoid compilation warning

  fclose(yyin);
  if (zzinStack.empty()) return 0;

  ZZFileState fileState = zzinStack.top();
  zzinFileName = fileState.inFileName_;
  zznumCharRead = fileState.numCharRead_;
  zzline = fileState.line_;
  zzcolumn = fileState.column_;
  zzinStack.pop();

  FILE* previn = fileState.file_;
  if (fseek(previn, zznumCharRead-1, SEEK_SET)!=0)
  {
    printf("\nERROR: In follex.y. Failed to seek to previous position in file "
           "%s, ln %d, col %d\n", zzinFileName.c_str(), zzline, zzcolumn+1);
    exit(-1);
  }
  yyrestart(previn);
  
}
%%
