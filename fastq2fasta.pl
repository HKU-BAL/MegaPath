#!/usr/bin/perl

use strict;
use warnings;
#use Cwd;

##### ##### ##### ##### #####

use Getopt::Std;
use vars qw( $opt_a $opt_c $opt_v);

# Usage
my $usage = "

fastq2fasta.pl - converts fastq files to fasta format.

       		      by
       		Brian J. Knaus
       		 December 2009

Copyright (c) 2009 Brian J. Knaus.
License is hereby granted for personal, academic, and non-profit use.
Commercial users should contact the author (http://brianknaus.com).

Great effort has been taken to make this software perform its said
task however, this software comes with ABSOLUTELY NO WARRANTY,
not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Usage: perl script.pl options
 required:
  -a   	fastq input file.
 optional:
  -c   	inset in nucleotides, used to remove barcodes [default = 0].
  -v   	verbose mode [T/F, default is F].

";

##### ##### ##### ##### #####
# Main.

my @temp;

while (<>){
  chomp($temp[0] = $_);		# First line is an id.
  chomp($temp[1] = <>);	# Second line is a sequence.
  chomp($temp[2] = <>);	# Third line is an id.
  chomp($temp[3] = <>);	# Fourth line is quality.

  # Prune first char.
  $temp[0] = substr($temp[0], 1);

  # Substring to inset value.
  $temp[1] = substr($temp[1], 0);

  # Print to fasta file.
  print STDOUT ">$temp[0]\n";
  print STDOUT "$temp[1]\n";
}

##### ##### ##### ##### #####
# EOF.