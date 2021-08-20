#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	maskLowerWithN.pl [options] <in.fa>
	
=head1 Options
	-l min len [15]
=cut

use warnings;
use strict;
use Getopt::Long;

my $min_len = 15;

GetOptions(
	"l=i" => \$min_len
);

unless (@ARGV >= 1) {
	die `pod2text $0`;
}


my $in;
if ($ARGV[0] eq "-") {
	$in = *STDIN;
} elsif ($ARGV[0] =~ /^(\S+)\.gz/) {
	open($in, "gzip -cd $ARGV[0] |") or die "cannot open $ARGV[0]";
} else {
	open($in, "<", "$ARGV[0]") or die "cannot open $ARGV[0]";
}

close($out_file);

sub split_fa {

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