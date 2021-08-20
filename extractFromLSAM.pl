#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	extractViralAndUnmap.pl [options] <in.lsam.id>
	
=head1 Options
	-n output a list of read names instead of fastq
	-t score threshold [40]
	-v extract Viral reads too
	-i append an "IGNORE" comment to those already mapped
	-g do not output a read with an "IGNORE" tag
	-s SE mode: do not output a mapped read even if its mate is unmapped
	-c append tags to comments
	
=cut

use warnings;
use strict;
use Getopt::Long;

my $scoreT = 40;
my $outputViral;
my $appendIgnore;
my $skipIgnoreTag;
my $seMode;
my $outputTag;
my $outputName;

GetOptions(
	"t=f" => \$scoreT,
	"v" => \$outputViral,
	"g" => \$skipIgnoreTag,
	"i" => \$appendIgnore,
	"s" => \$seMode,
	"c" => \$outputTag,
	"n" => \$outputName
);

unless (@ARGV == 1) {
	die `pod2text $0`;
}

die "-c and -i cannot be set at the same time" if ($outputTag and $appendIgnore);

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

	my ($name1, $flag1, $score1, $seq1, $qual1, $tag1, @opt1) = split "\t", $line1;
	my ($name2, $flag2, $score2, $seq2, $qual2, $tag2, @opt2) = split "\t", $line2;

	my $scoreCut = $scoreT < 1 ? $scoreT * (length($seq1) + length($seq2)) : $scoreT;

	if ($score1 < $scoreCut or $score2 < $scoreCut or ($outputViral and (grep(/^Viruses$/, @opt1) or grep(/^Viruses$/, @opt2)))) {
		my $comment1 = "";
		my $comment2 = "";
		
		if ($appendIgnore and $score1 >= $scoreCut) {
			$comment1 = " IGNORE";
		}

		if ($appendIgnore and $score2 >= $scoreCut) {
			$comment2 = " IGNORE";
		}

		if ($outputTag) {
			$comment1 = $tag1;
			$comment2 = $tag2;
		}

		my $alreadyOutput = 0;

		if ((!$skipIgnoreTag or !grep(/^IGNORE$/, @opt1)) and (!$seMode or $score1 < $scoreCut or grep(/^Viruses$/, @opt1))) {
			if ($outputName) {
				print $name1."\n";
				$alreadyOutput = 1;
			} else {
				print join("\n", "@".$name1."/1".$comment1, $seq1, "+", $qual1)."\n";
			}
		}

		if ((!$skipIgnoreTag or !grep(/^IGNORE$/, @opt2)) and (!$seMode or $score2 < $scoreCut or grep(/^Viruses$/, @opt2))) {
			if ($outputName) {
				print $name2."\n" if !$alreadyOutput;
			} else {
				print join("\n", "@".$name2."/2".$comment2, $seq2, "+", $qual2)."\n";
			}
		}
	}
}

close($in) unless $lsam eq "-";