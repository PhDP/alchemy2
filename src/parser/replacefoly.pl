#!/usr/bin/env perl

#Replace some of the bison directives so they are backwards-compatible

$file = $ARGV[0];
$bisoncommand = $ARGV[1];
if (!$file) 
{
  print "Please specify bison's input file\n";
  exit(0);
}

open (IN, $file) || die("ERROR: Unable to open $file");
{
  undef $/;
  $target = <IN>;
}
close (IN) || die("ERROR: Unable to close $file");

my $bisonversion = `$bisoncommand --version | grep "bison"`;
my @sp = split(/\s/,$bisonversion);
# last element of sp is version number
my @version = split(/\./,$sp[$#sp]);
# 2.4.1 means $version[0] is 2, $version[1] is 4, $version[2] is 1

if ($version[0] < 2 || ($version[0] >= 2 && $version[1] < 4))
{
  #bison 2.3 or earlier
  $target =~ s/.*\/\/bisonopencode/%{ \/\/bisonopencode/g;
  $target =~ s/.*\/\/bisonclosecode/%} \/\/bisonclosecode/g;
}
else
{
  #bison 2.4 or later
  $target =~ s/.*\/\/bisonopencode/%code { \/\/bisonopencode/g;
  $target =~ s/.*\/\/bisonclosecode/} \/\/bisonclosecode/g;
}

open(MODIFIED, ">$file") || die("ERROR: Unable to write $file");
print MODIFIED $target;
close(MODIFIED) || die("ERROR: Unable to close $file");
exit(0);

