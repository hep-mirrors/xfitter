#!/usr/bin/perl
#
# Input minuit.out.txt file. Finds the result of minimization
#

$file = shift;

open F,"<$file";

print "set title \n"; 
print "new YOUR PDF \n"; 
print "parameters \n";


$write = 0;
$fix = "FIX ";
while (<F>){
  if (/TERMINATED/ || /CONVERGED/) {
    $write = 1;
  }	
  if (/DERIVATIVE/ && ($write==1)) {
    $write = 2;
  }
  if ( ($write == 2) && (/\s+(\d+)\s+(\w+)\s+([-\w\.]+)\s+([-\w\.]+)/) ) {
    $i = $1;
    $n = $2;
    $v = $3;
    $e = $4;
    $e /= 10.;
    $_ = $4;
    if ( /fixed/) { 
       $fix = $fix . "$i ";
       $e = 0.01;
    }
    printf "%4d   %-16s  %16s %14.6e\n",$i,"\'$n\'",$v,$e;
  }
if (/EXTERNAL/) {$write = 0;}
}

#print "\n $fix \n";




