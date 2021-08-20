#!/bin/bash
#bbduk in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified December 1, 2015

Description:  Performs quality-trimming, artifact removal, linker-trimming, adapter trimming, and spike-in removal using BBDukF.
Performs human/cat/dog/mouse/microbe removal using BBMap.
NOTE!  This program uses hard-coded paths and will only run on Genepool.
It requires 39500m RAM for mousecatdoghuman, but only 1GB or so without them.

Usage:        rqcfilter.sh in=<input file> path=<output directory> library=<frag, clip, lfpe, or clrs> rna=<t or f>


Optional parameters (and their defaults)

Primary I/O parameters:
in=<file>           Input reads.
in2=<file>          Use this if 2nd read of pairs are in a different file.
path=null           Set to the directory to use for all output files.

Reference file paths:
ref=<file,file>     Comma-delimited list of additional reference files for filtering.
artifactdb=<file>   Override default Illumina artifacts file.
rnadb=<file>        Override default rna spikein file.
dnadb=<file>        Override default dna spikein file.
ribodb=<file>       Override defailt ribosomal file.
fragadapter=<file>  Override default fragment library adapter file.
rnaadapter=<file>   Override default rna library adapter file.
lfpelinker=<file>   Override default lfpe linker file.
clrslinker=<file>   Override default clrs linker file.
cliplinker=<file>   Override default clip linker file.
phixref=<file>      Override default phiX reference file.

Output parameters:
scafstats=scaffoldStats.txt  Scaffold stats file name (how many reads matched which reference scaffold) .
kmerstats=kmerStats.txt      Kmer stats file name (duk-like output).
log=status.log               Progress log file name.
filelist=file-list.txt       Progress log file name.
stats=filterStats.txt        Overall stats file name.
ihist=ihist_merge.txt        Insert size histogram name.  Set to null to skip merging.
outribo=ribo.fq.gz           Output for ribosomal reads, if removeribo=t.
reproduceName=reproduce.sh   Name of shellscript to reproduce these results.
usetmpdir=t                  Write temp files to TMPDIR.
tmpdir=                      Override TMPDIR.

Adapter trimming parameters:
trimhdist=1         Hamming distance used for trimming.
trimhdist2=         Hamming distance used for trimming with short kmers.  If unset, trimhdist will be used.
trimk=23            Kmer length for trimming stage.
mink=11             Minimum kmer length for short kmers when trimming.
trimfragadapter=f   Trim all known Illumina adapter sequences, including TruSeq and Nextera.
trimrnaadapter=f    Trim Illumina TruSeq-RNA adapters.

Quality trimming parameters:
qtrim=f             Trim read ends to remove bases with quality below minq.  Performed AFTER looking for kmers.
                    Values: t (trim both ends), f (neither end), r (right end only), l (left end only).
trimq=10            Trim quality threshold.  Must also set qtrim for direction.
minlength=45        (ml) Reads shorter than this after trimming will be discarded.  Pairs will be discarded only if both are shorter.
mlf=0.333           (minlengthfraction) Reads shorter than this fraction of original length after trimming will be discarded.
minavgquality=5     (maq) Reads with average quality (before trimming) below this will be discarded.
maxns=0             Reads with more Ns than this will be discarded.
forcetrimmod=5      (ftm) If positive, right-trim length to be equal to zero, modulo this number.

Mapping parameters (for vertebrate contaminants):
mapk=14             Kmer length for mapping stage (9-15; longer is faster).
removehuman=f       (human) Remove human reads via mapping.
removedog=f         (dog) Remove dog reads via mapping.
removecat=f         (cat) Remove cat reads via mapping.
removemouse=f       (mouse) Remove mouse reads via mapping.
aggressive=f        Aggressively remove human reads (and cat/dog/mouse) using unmasked references.
*NOTE: If cat, dog, human, and mouse removal are all true, they will be done together.*
humanpath=          Use this to override the index for human mapping; default is a masked HG19 index.
catpath=            Use this to override the index for cat mapping.
dogpath=            Use this to override the index for dog mapping.
mousepath=          Use this to override the index for mouse mapping.
mapref=             Remove contaminants by mapping to this fasta file (or comma-delimited list).

Microbial contaminant removal parameters:
removemicrobes=f    (microbes) Remove common contaminant microbial reads via mapping, and place them in a separate file.
taxlist=            (tax) Remove these taxa from the database before filtering.  Typically, this would be the organism name or NCBI ID, or a comma-delimited list.  Organism names should have underscores instead of spaces, such as Escherichia_coli.
taxlevel=order      (level) Level to remove.  For example, 'phylum' would remove everything in the same phylum as entries in the taxlist.
taxtree=            (tree) Override location of the TaxTree file.
gitable=            Ovverride location of the gitable file.
loadgitable=f       Controls whether gi numbers may be used for taxonomy.
microberef=         Path to fasta file of microbes.

Filtering parameters (for artificial and microbial contaminants):
dna=t               Remove reads containing DNA-specific artifacts.
rna=f               Remove reads containing  RNA-specific artifacts.
phix=t              Remove  reads containing phiX kmers.
pjet=t              Remove  reads containing PJET kmers.
maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to increase sensitivity in the presence of errors.
maxbadkmers=0       (mbk) Reads with more than this many contaminant kmers will be discarded.
filterhdist=1       Hamming distance used for filtering.
filterqhdist=1      Query hamming distance used for filtering.

Ribosomal filtering parameters:
ribohdist=1         Hamming distance used for rRNA removal.
riboedist=0         Edit distance used for rRNA removal.
removeribo=f        (ribo) Remove ribosomal reads via kmer-matching, and place them in a separate file.

Other processing parameters:
threads=auto        (t) Set number of threads to use; default is number of logical processors.
library=frag        Set to 'frag', 'clip', 'lfpe', or 'clrs'.
filterk=27          Kmer length for filtering stage.
rcomp=t             Look for reverse-complements of kmers in addition to forward kmers.
nexteralmp=f        Split into different files based on Nextera LMP junction sequence.  Only for Nextera LMP, not normal Nextera.
extend=f            Extend reads during merging to allow insert size estimation of non-overlapping reads.
monitor=f           Kill this process if it crashes.  monitor=600,0.01 would kill after 600 seconds under 1% usage.
barcodefilter=crash Crash when improper barcodes are discovered.  Set to 'f' to disable or 't' to just remove improper barcodes.
barcodes=           A comma-delimited list of barcodes or files of barcodes.
pigz=t              Use pigz for compression.
unpigz=t            Use pigz for decompression.
khist=f             Set to true to generate a kmer-frequency histogram of the output data.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

*****   All additional parameters supported by BBDuk may also be used, and will be passed directly to BBDuk   *****

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"
NATIVELIBDIR="$DIR""jni/"

z="-Xmx1g"
z2="-Xms1g"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
    
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 3200m 84

#	if [ $RAM -ge 31000 ]; then
#		RAM=31000
#	elif [ 400 -ge $RAM ]; then
#		RAM=400
#	fi
	
#	local HOSTNAME=`hostname`
#	if [[ $NSLOTS == 8 ]] && [[ $HOSTNAME == sgi* ]]; then

	if [[ $NSLOTS == 8 ]]; then
		RAM=39200
	fi
	
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}

calcXmx "$@"


rqcfilter() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.7_64bit
		module load pigz
	export TZ="America/Los_Angeles" 
	fi
	local CMD="java -Djava.library.path=$NATIVELIBDIR $EA $z -cp $CP jgi.RQCFilter usejni $@"
    echo $CMD

	echo $CMD >&2
	eval $CMD
}

rqcfilter "$@"
