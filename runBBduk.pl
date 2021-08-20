#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	runBBduk.pl [options] <read1.fq> [read2.fq]
	
=head1 Options
	-a            Trim adapters using BBDuk
	-q <int>      Quality fitlering threshold [10]
	-e <float>    Entropy for low-complexity filtering [0.75]
	-p <str>      output prefix
=cut

use warnings;
use strict;
use Getopt::Long;

my $bbduk = "/nas3/dhli_1/software/bbmap/bbduk2.sh";
my $adapters = "/nas3/dhli_1/software/bbmap/resources/adapters.fa";

my $qual_t = 10;
my $trim_ad;
my $entropy = 0.75;
my $out_prefix;

GetOptions(
	"a" => \$trim_ad,
	"q=i" => \$qual_t,
	"e=f" => \$entropy,
	"p=s" => \$out_prefix
);

unless (@ARGV >= 1 and defined $out_prefix) {
	die `pod2text $0`;
}

my $is_pe = (@ARGV >= 2)? 1 : 0;
my $read1 = $ARGV[0];
my $read2 = $ARGV[1] if ($is_pe);

print STDERR trim()." | ".low_compl_filter()."\n";
system(trim()." | ".low_compl_filter()) == 0 or die "Failed: $?";

sub trim {
	my @trim_cmd = ($bbduk, "qtrim=rl", "trimq=".$qual_t, "qin=33"); # WARNING hard code
	push(@trim_cmd, "threads=12");
	push(@trim_cmd, "minlength=50");
	push(@trim_cmd, "in=".$read1, "out=stdout.fq");
	push(@trim_cmd, "ref=".$adapters, "hdist=1") if ($trim_ad);
	push(@trim_cmd, "in2=".$read2) if ($is_pe);
	return join(" ", @trim_cmd);
}

sub low_compl_filter {
	my @filter_cmd = ($bbduk, "entropy=".$entropy, "in=stdin.fq", "outm=".$out_prefix.".low_compl.fq.gz");
	push(@filter_cmd, "threads=12", "qin=33"); # WARNING hard code
	push(@filter_cmd, "out=".$out_prefix.".bbduk_1.fq.gz");
	push(@filter_cmd, "out2=".$out_prefix.".bbduk_2.fq.gz") if $is_pe;
	return join(" ", @filter_cmd);
}