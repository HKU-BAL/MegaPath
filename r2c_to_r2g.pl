#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	r2c_to_r2g.pl [options] <read2contigs.lsam> <contigs.lsam>
	
=head1 Options
	-t score threshold [40]
	
=cut

use warnings;
use strict;
use Getopt::Long;

my %contig2gid;
my $scoreT = 40;

GetOptions(
	"t=f" => \$scoreT
);

unless (@ARGV == 2) {
	die `pod2text $0`;
}

my $lsam = $ARGV[0];
my $ctglsam = $ARGV[1];
my $in_ctg;

if ($ctglsam eq "-") {
	$in_ctg = *STDIN;
} elsif ($ctglsam =~ /^(\S+)\.gz/) {
	open($in_ctg, "gzip -cd $ctglsam |") or die "cannot open $ctglsam";
} else {
	open($in_ctg, "<", "$ctglsam") or die "cannot open $ctglsam";
}

while (<$in_ctg>) {
	chomp;
	my ($name, undef, undef, undef, undef, $labels, @opts) = split "\t";
	next unless ($name =~ /^contig_(\S+)/);
	$contig2gid{$1} = $labels;
}
close($in_ctg) unless $ctglsam eq "-";

my $in;
if ($lsam eq "-") {
	$in = *STDIN;
} elsif ($lsam =~ /^(\S+)\.gz/) {
	open($in, "gzip -cd $lsam |") or die "cannot open $lsam";
} else {
	open($in, "<", "$lsam") or die "cannot open $lsam";
}


while (<$in>) {
	my ($name, $flag, $score, $seq, $qual, $ctgs, @opts) = split /\s/;
	next if grep(/^IGNORE$/, @opts);
	my @contigs = split ";", $ctgs;
	my @gids;

	foreach my $contig (@contigs) {
		next if ($contig eq "*");
		my ($score, $ctg) = split ",", $contig;
		push(@gids, $contig2gid{$ctg}) if $score > $scoreT and exists $contig2gid{$ctg};
	}

	my $gidTag = "*";
	$gidTag = join(";", @gids) if scalar(@gids) > 0;

	print join("\t", $name, $flag, $score, "*", "*", $gidTag, @opts)."\n";
}

close($in) unless $lsam eq "-";
