#!/usr/bin/perl

#Replace some of the C code in bison's output so that it is C++ compliant

$file = $ARGV[0];;
if (!$file) 
{
  print "Please specify bison's output file\n";
  exit(0);
}

open (IN, $file) || die("ERROR: Unable to open $file");
{
  undef $/;
  $target = <IN>;
}
close (IN) || die("ERROR: Unable to close $file");

#with bison 2.0
$target =~ s/fol.tab.c/fol.cpp/g;
$target =~ s/parser stack overflow/Parser stack overflow. Try disambiguating with parentheses/g;
$target =~ s/yymsg = YYMALLOC \(yysize\)/yymsg = (char*) (YYMALLOC (yysize))/g;
$target =~ s/return yytname\[yytoken\];/if (yytoken < 0) { char* buf = new char[1]; buf[0]='\\0'; return buf; } return yytname[yytoken];/g;
$target =~ s/syntax error/parse error/g;
$target =~ s/yynewItem->yystate.yyloc = \*yylocp;/yylocp->yydummy = 'a'; yynewItem->yystate.yyloc = *yylocp;/g;

#with bison 2.3
$target =~ s/YYLTYPE yyerrloc;/YYLTYPE yyerrloc;yyerrloc.yydummy = 'a';/g;

#Using bison 2.0, we do not need the following 3 lines.
#$target =~ s/YYMALLOC \(16 \*/(yyGLRState**) YYMALLOC (16 */g;
#$target =~ s/YYMALLOC \(yysize \*/(yyGLRStackItem*) YYMALLOC (yysize */g;
#$target =~ s/YYREALLOC \(yystack/(yyGLRState**) YYREALLOC (yystack/g;

open(MODIFIED, ">$file") || die("ERROR: Unable to write $file");
print MODIFIED $target;
close(MODIFIED) || die("ERROR: Unable to close $file");
exit(0);
