#!/usr/bin/env perl

=head1 Author
	Dinghua Li <dhli@cs.hku.hk>
	
=head1 Usage
	calcAccuracy.pl <fastq_sim.fq> <xxx.lsam.labels|xxx.mpa>
=cut

my %id2sp = (
	NC001586 => "Alphapapillomavirus_1",
	NC001531 => "Betapapillomavirus_1",
	YP001816783 => "Chapare_mammarenavirus",
	YP717901 => "Chipapillomavirus_1",
	NP041306 => "Deltapapillomavirus_1",
	NP899211 => "Guanarito_mammarenavirus",
	NC006577 => "Human_coronavirus_HKU1",
	NC001802 => "Human_immunodeficiency_virus_1",
	NC001722 => "Human_immunodeficiency_virus_2",
	NC002016 => "Influenza_A_virus",
	NP056657 => "Influenza_B_virus",
	YP002302328 => "Influenza_C_virus",
	NP899217 => "Junin_mammarenavirus",
	NC004355 => "Kyasanur_forest_disease_virus",
	NC003690 => "Langat_virus",
	NC005078 => "Machupo_mammarenavirus",
	NC001608 => "Marburg_marburgvirus",
	NC005062 => "Omsk_hemorrhagic_fever_virus",
	NC006313 => "Sabia_mammarenavirus",
	NC004718 => "Severe_acute_respiratory_syndrome-related_coronavirus",
	NC001549 => "Simian_immunodeficiency_virus",
	NC006432 => "Sudan_ebolavirus"
);

if (@ARGV < 2) {
	die `pod2text $0`;
}

my $ans_fq = $ARGV[0];
my $query = $ARGV[1];
my %ans;
my %correc_sp, %false_sp, %size_sp;

# parse answer
open(ANS_FQ, '<', $ans_fq) or die $!;

while (<ANS_FQ>) {
	my ($name, $ans, undef) = split /\s/;
	$name = substr $name, 1;
	my ($target, $id) = split /A/, $ans;
	$ans{$name} = $id2sp{$target};
	$size_sp{$id2sp{$target}}++;
	$correct_sp{$id2sp{$target}} = 0;
	$false_sp{$id2sp{$target}} = 0;

	<ANS_FQ>;
	<ANS_FQ>;
	<ANS_FQ>;
}

close(ANS_FQ);

open(IN, '<', $query) or die $!;

my $correct = 0;
my $false = 0;

if ($query =~ /.mpa$/) {
	# kraken MPA output
	while (<IN>) {
		my ($name, $a) = split /\s/;
		if ($a =~ /s__([\w|\-]+)$/) {

			if ($1 eq $ans{$name}) {
				$correct++;
				$correct_sp{$ans{$name}}++;
			} else {
				$false++;
				$false_sp{$ans{$name}}++;
			}
		}
	}
} else {
	# MegaPath's lsam
	while (<IN>) {
		my ($name, undef, undef, undef, undef, $a) = split /\s/;
		if ($name =~ /^@/) {
			$name = substr $name, 1;
		}
		if (!($a =~ /;/)) {
			if ($a =~ /species\|([\w|\-]+),/) {
				if ($1 eq $ans{$name}) {
					$correct++;
					$correct_sp{$ans{$name}}++;
				} else {
					$false++;
					$false_sp{$ans{$name}}++;
				}
			}
		}
	}
}

print $correct, "\t", $false, "\n";

foreach my $sp (keys %size_sp) {
	print join("\t", $sp, $size_sp{$sp}, $correct_sp{$sp}, $false_sp{$sp}), "\n";
}

print "SENSITIVITY\t", $correct * 100.0 / (scalar keys %ans), "\n";
print "FDR\t", $false * 100.0 / ($correct + $false), "\n";

close(IN);