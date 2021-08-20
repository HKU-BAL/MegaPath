#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	reassign.pl [options] in.lsam.id > out.lsam.id

=head1 Options
	-u <int>   [20]
	-c <float> [0.05]
	-t <int>   [40]
=cut

use warnings;
use Getopt::Long;
use strict;

my $ut = 20;
my $ct = 0.05;
my $st = 40;

GetOptions( 
	"u=i" => \$ut,
	"c=f" => \$ct,
	"t=i" => \$st
);

die `pod2text $0`if (@ARGV == 0);

my $in;
if ($ARGV[0] eq "-") {
	$in = *STDIN;
} elsif ($ARGV[0] =~ /^(\S+)\.gz/) {
	open($in, "gzip -cd $ARGV[0] |") or die "cannot open $ARGV[0]";
} else {
	open($in, "<", "$ARGV[0]") or die "cannot open $ARGV[0]";
}

my %count;
my %uniq_count;
my %intersect;

# 1-pass
while (<$in>) {
	chomp;
	my ($q, undef, $score, undef, undef, $t, @opt) = split /\s/;
	next if ($t eq "*" or $t eq "");
	next if ($score < $st);
	my @tagsWithScore = split ';', $t;
	my @tags;
	my $n = scalar(@tagsWithScore);

	foreach my $s (@tagsWithScore) {
		my @a = split ',', $s;
		push(@tags, $a[1]);
	}

	foreach my $s (@tags) {
		$count{$s}++;
	}

	$n = scalar(@tags);

	if ($n == 1) {
		$uniq_count{$tags[0]}++;
	} else {
		for (my $i = 0; $i < $n; $i++) {
			for (my $j = $i + 1; $j < scalar(@tags); $j++) {
				$intersect{$tags[$i].";".$tags[$j]}++;
			}
		}
	}
}

for my $s (keys %count) {
	unless (exists $uniq_count{$s}) {
		$uniq_count{$s} = 0;
	}
}

close $in unless $ARGV[0] eq "-";

my %weak_explain;
my %weak_explained;
my %explain;

sub weak_exp {
	my ($s1, $s2, $n_inter) = @_;
	# print STDERR $s1, " ", $s2, "\n";
	return 0 if $uniq_count{$s1} <= $ut * $uniq_count{$s2};
	return 0 if $count{$s1} - $n_inter <= $ct * $count{$s1};
	return 1;
};

# calculate weakly explain
foreach my $pair (keys %intersect) {
	my ($s1, $s2) = split ";", $pair;
	if (weak_exp($s1, $s2, $intersect{$pair})) {
		$weak_explain{$pair} = 1;
		$weak_explained{$s2} = 1;
	} elsif (weak_exp($s2, $s1, $intersect{$pair})) {
		$weak_explain{$s2.";".$s1} = 1;
		$weak_explained{$s1} = 1;
	}
}

# calculate explain
foreach my $pair (keys %weak_explain) {
	my ($s1, $s2) = split ";", $pair;
	if (!(exists $weak_explained{$s1})) {
		$explain{$pair} = 1;
		print STDERR $s1, " explains ", $s2, "\n";
	}
}

# 2-pass
if ($ARGV[0] eq "-") {
	$in = *STDIN;
} elsif ($ARGV[0] =~ /^(\S+)\.gz/) {
	open($in, "gzip -cd $ARGV[0] |") or die "cannot open $ARGV[0]";
} else {
	open($in, "<", "$ARGV[0]") or die "cannot open $ARGV[0]";
}

while (<$in>) {
	chomp;
	my ($q, $flag, $score, $seq, $qual, $t, @opt) = split /\s/;
	my @remain_tags;

	if ($t eq "*") {
		print join("\t", $q, $flag, $score, $seq, $qual, $t, @opt)."\n";
		next;
	}

	my @tagsWithScore = split ';', $t;
	my @tags;
	my %scores;

	my $n = scalar(@tagsWithScore);

	foreach my $s (@tagsWithScore) {
		my @a = split ',', $s;
		push(@tags, $a[1]);
		$scores{$a[1]} = $a[0];
	}

	foreach my $s (@tags) {
		my $remain = 1;
		foreach my $t (@tags) {
			if (exists $explain{$t.";".$s}) {
				$remain = 0;
				last;
			}
		}

		push(@remain_tags, $scores{$s}.",".$s) if $remain;
	}

	print join("\t", $q, $flag, $score, $seq, $qual, join(";", @remain_tags), @opt)."\n";
}

close $in unless $ARGV[0] eq "-";