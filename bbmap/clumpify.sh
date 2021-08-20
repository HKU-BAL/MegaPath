#!/bin/bash
#clumpify in=<infile> out=<outfile> groups=<number>

usage(){
echo "
Written by Brian Bushnell
Last modified July 5, 2016

Description:  Sorts sequences to put similar reads near each other.
This is a wrapper for KmerSplit and KmerSort.
Works best with single-ended data.

Usage:   clumpify.sh in=<file> out=<file> groups=<number>

Input may be fasta or fastq, compressed or uncompressed.

Optional parameters (and their defaults)

in=<file>           Input file.
out=<file>          Output file.  May not be standard out.
groups=16           Use this many intermediate files (to save memory).
rcomp=t             Give read clumps the same orientation.
                    Should be disabled for paired reads.
rename=t            Add kmer information to the name.
consensus=t         Generate consensus reads from clumps.
divisor=80m         (div) Use a prime number at least this big as the divisor.
k=31                Use kmers of this length (1-31).
mincount=0          Ignore pivot kmers with count less than this.
prefilter=t         Use a prefilter if counting kmers.
overwrite=f         (ow) Set to false to force the program to abort rather 
                    than overwrite an existing file.
ziplevel=2          (zl) Set to 1 (lowest) through 9 (max) to change 
                    compression level; lower compression is faster.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding
                    the program's automatic memory detection.  -Xmx20g will 
                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                    The max is typically 85% of physical memory.

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

z="-Xmx2g"
z2="-Xms2g"
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
	freeRam 2000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

clumpify() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.7_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP clump.Clumpify $@"
	echo $CMD >&2
	eval $CMD
}

clumpify "$@"
