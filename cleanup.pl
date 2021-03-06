#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	cleanup.pl [options] in.lsam.id > out.lsam.id

=head1 Options
	-c <int[,int...]> contaminant taxonomy ids [9606,32630]
	-t <float>        ratio to determine a read is homologous to contaminant [0.95]
	-p <float>        mark a species as contaminant-like if more than *float* reads are homologous to contaminant [0.8]
=cut

use warnings;
use Getopt::Long;
use strict;

my $str_tids = "9606,32630"
my $top = 0.95;
my $perc = 0.8;

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
	my @tags = split ';', $t;

	foreach my $s (@tags) {
		$count{$s}++;
	}

	if (scalar(@tags) == 1) {
		$uniq_count{$tags[0]}++;
	} else {
		for (my $i = 0; $i < scalar(@tags); $i++) {
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
	my @tags = split ';', $t;
	my @remain_tags;

	foreach my $s (@tags) {
		$count{$s}++;
	}

	foreach my $s (@tags) {
		my $remain = 1;
		foreach my $t (@tags) {
			if (exists $explain{$t.";".$s}) {
				$remain = 0;
				last;
			}
		}

		push(@remain_tags, $s) if $remain;
	}

	print join("\t", $q, $flag, $score, $seq, $qual, join(";", @remain_tags), @opt)."\n";
}

close $in unless $ARGV[0] eq "-";