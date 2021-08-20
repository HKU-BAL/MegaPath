#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	extractViralAndUnmap.pl [options] <in.lsam.labels>
	
=head1 Options
	-t score threshold [40]
	-v extract Viral reads too
	-i append an "IGNORE" comment to those already mapped
	-g do not output a read with an "IGNORE" tag
	
=cut

use warnings;
use strict;
use Getopt::Long;

my $scoreT = 40;
my $outputViral;
my $appendIgnore;
my $skipIgnoreTag;

GetOptions(
	"t=f" => \$scoreT,
	"v" => \$outputViral,
	"g" => \$skipIgnoreTag,
	"i" => \$appendIgnore
);

unless (@ARGV == 1) {
	die `pod2text $0`;
}

my $lsam = $ARGV[0];

my $in;
if ($lsam eq "-") {
	$in = *STDIN;
} elsif ($lsam =~ /^(\S+)\.gz/) {
	open($in, "gzip -cd $lsam |") or die "cannot open $lsam";
} else {
	open($in, "<", "$lsam") or die "cannot open $lsam";
}

while (my $line1 = <$in>) {
	chomp $line1;
	my $line2 = <$in>;
	chomp $line2;

	my ($name1, $flag1, $score1, $seq1, $qual1, $tag1) = split "\t", $line1;
	my ($name2, $flag2, $score2, $seq2, $qual2, $tag2) = split "\t", $line2;

	if ($score1 < $scoreT or $score2 < $scoreT or ($outputViral and (index($tag1, "Viruses") != -1 or index($tag2, "Viruses") != -1))) {
		my $comment1 = "";
		my $comment2 = "";
		if ($appendIgnore and $score1 >= $scoreT) {
			$comment1 = " IGNORE";
		}
		if ($appendIgnore and $score2 >= $scoreT) {
			$comment2 = " IGNORE";
		}

		if (!$skipIgnoreTag or $tag1 ne "IGNORE") {
			print join("\n", "@".$name1."/1".$comment1, $seq1, "+", $qual1)."\n";
		}

		if (!$skipIgnoreTag or $tag2 ne "IGNORE") {
			print join("\n", "@".$name2."/2".$comment2, $seq2, "+", $qual2)."\n";
		}
	}
}

close($in) unless $lsam eq "-";