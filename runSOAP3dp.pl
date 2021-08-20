#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	runSOAP3dp.pl [options] <soap3-dp-index> <read1.fq> [read2.fq]
	
=head1 Options
	-p <str>      output prefix
	-b            use BWA mem for alignment
=cut

use warnings;
use strict;
use Getopt::Long;

my $soap3dp = "/nas3/dhli_1/software/SOAP3-dp/soap3-dp";
my $bwa = "bwa";

my $out_prefix;
my $use_bwa;

GetOptions( 
	"p=s" => \$out_prefix,
	"b" => \$use_bwa
);

unless (@ARGV >= 2 and defined $out_prefix) {
	die `pod2text $0`;
}

my $is_pe = (@ARGV >= 3)? 1 : 0;
my $index = $ARGV[0];
my $read1 = $ARGV[1];
my $read2 = $ARGV[2] if ($is_pe);

if ($use_bwa) {
	bwa() == 0 or die "error on running BWA";
} else {
	aln() == 0 or die "error on running SOAP3-dp";
}

sub aln {
	my @aln_cmd = ($soap3dp, ($is_pe ? "pair" : "single"));
	push(@aln_cmd, $index);
	push(@aln_cmd, $read1);
	push(@aln_cmd, $read2) if ($is_pe);
	push(@aln_cmd, "-o", $out_prefix);
	return system(@aln_cmd);
}

sub bwa {
	my @aln_cmd = ($bwa, "mem");
	push(@aln_cmd, "-t24", "-M");
	push(@aln_cmd, "-p") unless ($is_pe);
	push(@aln_cmd, $index);
	push(@aln_cmd, $read1);
	push(@aln_cmd, $read2) if ($is_pe);
	push(@aln_cmd, ">", $out_prefix.".sam");
	my $real_cmd = join(' ', @aln_cmd);
	return system($real_cmd);
}