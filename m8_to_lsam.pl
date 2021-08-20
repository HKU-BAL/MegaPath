#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
	This script convert the m8 output of (AC-)DIAMOND to lsam

=head1 Usage
	m8_to_lsam.pl <in.m8 >out.lsam
=cut

use warnings;
use strict;
use Getopt::Long;

my $last_q = undef;
my @targets;
my $score = 0;

while (<>) {
	chomp;
	my @rec = split "\t", $_;
	if ($last_q ne $rec[0]) {
		if (defined $last_q) {
			my $tgs = scalar(@targets) > 0 ? join(";", @targets) : "*";
			print join("\t", $last_q, 0, $score, "*", "*", $tgs), "\n";
		}
		$last_q = $rec[0];
		@targets = ();
		$score = -1;
	}

	if ($rec[11] > $score) {
		$score = $rec[11];
	}

	my @tids = split /0x1/, $rec[1];
	foreach my $tid (@tids) {
		push(@targets, join(",", $rec[11], $tid));
	}
}

if (defined $last_q) {
	my $tgs = scalar(@targets) > 0 ? join(";", @targets) : "*";
	print join("\t", $last_q, 0, $score, "*", "*", $tgs), "\n";
}