#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	sam2cfq.pl [options] <in.sam|bam> >out.{fq|lsam}
	
=head1 Options
	-A <int>    match score [1]
	-B <int>    mismatch penalty [2]
	-O <int>    gap open penalty [3]
	-E <int>    gap extension penalty [1]
	-d <float>  dropout ratio [0.95]
	-l          output lsam
	
=cut

use warnings;
use strict;
use Getopt::Long;

my ($A, $B, $O, $E) = (1, 2, 3, 1);
my $dropout = 0.95;
my $lsam;

GetOptions(
	"A=i" => \$A,											
	"B=i" => \$B,
	"O=i" => \$O,
	"E=i" => \$E,
	"d=f" => \$dropout,
	"l" => \$lsam
);

unless (@ARGV == 1) {
	die `pod2text $0`;
}

my $sam_fn = $ARGV[0];

open(SAM, "samtools view -F0x100 $sam_fn 2>/dev/null |") or die "Fail to open $sam_fn with samtools";

while (my $line = <SAM>) {
	my ($query, $flag, $max_score, $seq, $qual, @aln) = getAln($line);
	$line = <SAM>;
	my ($query2, $flag2, $max_score2, $seq2, $qual2, @aln2) = getAln($line);

	my %all_ref = collectRef(@aln);
	my %all_ref2 = collectRef(@aln2);

	my (@pe_aln, @pe_aln2);
	$max_score = 0;
	$max_score2 = 0;

	while (@aln) {
		my $s = shift @aln;
		my $g = shift @aln;

		if (exists $all_ref2{$g}) {
			$s += $all_ref2{$g};
		}

		$max_score = $s if $s > $max_score;
		push @pe_aln, $s, fetchTaxid($g);
	}

	while (@aln2) {
		my $s = shift @aln2;
		my $g = shift @aln2;

		if (exists $all_ref{$g}) {
			$s += $all_ref{$g};
		}

		$max_score2 = $s if $s > $max_score2;
		push @pe_aln2, $s, fetchTaxid($g);
	}

	my (%final_aln, %final_aln2);

	while (@pe_aln) {
		my $s = shift @pe_aln;
		my $g = shift @pe_aln;
		$final_aln{$g} = $s if $s >= $max_score * $dropout;
	}

	while (@pe_aln2) {
		my $s = shift @pe_aln2;
		my $g = shift @pe_aln2;
		$final_aln2{$g} = $s if  $s >= $max_score2 * $dropout;
	}

	print join("\t", $query, $flag, $max_score, $seq, $qual, join(';', map{ "$final_aln{$_},$_" } (sort keys %final_aln)));
	print "*" unless keys %final_aln;
	print "\n";
	print join("\t", $query2, $flag2, $max_score2, $seq2, $qual2, join(';', map{ "$final_aln2{$_},$_" } (sort keys %final_aln2)));
	print "*" unless keys %final_aln2;
	print "\n";
}

sub collectRef {
	my @aln = @_;
	my %ret;
	while (@aln) {
		my $s = shift @aln;
		my $g = shift @aln;

		$ret{$g} = $s if (!exists $ret{$g} or $ret{$g} < $s);
	}
	return %ret;
}

sub getAln {
	my ($line) = @_;
	chomp $line;
	my ($query, $flag, $ref, undef, undef, $cigar, undef, undef, undef, $seq, $qual, @tags) = split "\t", $line;

	if ($flag & 0x10) {
		$seq = rev_comp($seq);
		$qual = reverse($qual);
	}

	my @aln = ();
	my $max_score = -99999999;

	my $ed = -1;
	my $mm = 0;
	foreach my $tag (@tags) {
		if ($tag =~ /^NM:i:(\d+)$/) {
			$ed = $1;
		} elsif ($tag =~ /^XA:Z:(.+)$/) {
			my @xa = split ";", $1;
			foreach my $a (@xa) {
				my @rec = split ",", $a;
				my @indel_rec = calc_indel($rec[2]);
				my $score = $indel_rec[2] * $A - ($rec[3] - $indel_rec[0]) * ($A + $B) - $indel_rec[1];

				push(@aln, $score, $rec[0]);
				$max_score = $score if $score > $max_score;
			}
		} elsif ($tag =~ /^XC:Z:(.+)$/) {
			my @xc = split ";", $1;
			foreach my $c (@xc) {
				my ($g, $s) = split ",", $c;
				push(@aln, $s, $g);
				$max_score = $s if $s > $max_score;
			}
		}
	}

	if (($flag & 0x4) == 0) {
		my ($indel, $indel_p, $m_len) = calc_indel($cigar);
		my $score = $m_len * $A - ($ed - $indel) * ($A + $B) - $indel_p;
		if ($ed < 0) {
			$score = $m_len * $A - $mm * ($A + $B) - $indel_p;
		}
		push(@aln, $score, $ref);
		$max_score = $score if $score > $max_score;
	}

	return ($query, $flag, $max_score, $seq, $qual, @aln);
}

sub fetchTaxid {
	my ($ref) = @_;
	die "cannot fetch taxid" unless $ref =~ /kraken:taxid\|(\d+)$/;
	return $1;
}

sub calc_indel {
	my ($cigar) = @_;
	my @cig_fields = split(/([A-Z=])/, $cigar);
	# print $cigar, "\n";
	my ($indel, $indel_p, $m_len) = (0, 0, 0);

	while (@cig_fields) {
		my $l = shift @cig_fields;
		my $op = shift @cig_fields;
		if ($op eq 'M' or $op eq 'X' or $op eq '=') {
			$m_len += $l;
		} elsif ($op eq 'I') {
			$indel += $l;
			$indel_p += $O + $E * $l;
		} elsif ($op eq 'D') {
			$indel += $l;
			$indel_p += $O + $E * $l;
		} 
	}

	return ($indel, $indel_p, $m_len);
}

sub rev_comp {
	my ($seq) = @_;
	$seq =~ tr /atcgATCG/tagcTAGC/;
	$seq = reverse($seq);
	return $seq;
}