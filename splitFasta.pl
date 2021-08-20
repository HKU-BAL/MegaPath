#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	splitFasta.pl [options] <in.fa>
	
=head1 Options
	-c <float>    Chunk size [3.9e9]
	-p <str>      Output prefix [in.fa]
=cut

use warnings;
use strict;
use Getopt::Long;

my $chunk_size = 3.9e9;
my $prefix = undef;

GetOptions(
	"c=f" => \$chunk_size,
	"p=s" => \$prefix
);

unless (@ARGV >= 1) {
	die `pod2text $0`;
}

my $fa_fn = $ARGV[0];

unless (defined($prefix)) {
	if ($fa_fn =~ /^(\S+)\.gz/) {
		$prefix = $1;
	} else {
		$prefix = $fa_fn;
	}
}

my $file_idx = 0;
my $out_file;
my $acc_size = 0;
my $seq_size = 0;
open($out_file, ">", $prefix.".".$file_idx);

split_fa();
close($out_file);

sub split_fa {
	my $in;
	if ($fa_fn eq "-") {
		$in = *STDIN;
	} elsif ($fa_fn =~ /^(\S+)\.gz/) {
		open($in, "gzip -cd $fa_fn |") or die "cannot open $fa_fn";
	} else {
		open($in, "<", "$fa_fn") or die "cannot open $fa_fn";
	}

	my $name = undef;
	my @seqs = undef;

	while (my $line = <$in>) {
		chomp $line;
		if ($line =~ /^>/) {
			write_seq($name, @seqs) if defined($name);
			$name = $line;
			@seqs = ();
			$seq_size = 0;
		} else {
			$acc_size += length($line);
			push(@seqs, $line);
		}
	}

	write_seq($name, @seqs) if defined($name);
	close($in) unless $fa_fn eq "-";
}

sub write_seq {
	my ($name, @seqs) = @_;
	$acc_size += $seq_size;
	if ($acc_size > $chunk_size) {
		$acc_size = $seq_size;
		$file_idx++;
		close($out_file);
		open($out_file, ">", $prefix.".".$file_idx);
	}
	print $out_file $name."\n";
	foreach my $seq (@seqs) {
		print $out_file $seq."\n";
	}
}