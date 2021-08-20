#!/bin/bash
#calctruequality in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified December 11, 2015

Description:  Calculates the observed quality scores from a sam file.
Generates matrices for use in recalibrating quality scores.

Usage:        calctruequality.sh in=<file,file,...file> path=<directory>


Parameters (and their defaults)

Input parameters:
in=<file,file>      Sam file or comma-delimited list of files.  Must use = and X cigar symbols.
reads=-1            Stop after processing this many reads (if positive).

Output parameters:
overwrite=t         (ow) Set to true to allow overwriting of existing files.
path=.              Directory to write quality matrices (within /ref subdir).
write=t             Write matrices.
showstats=t         Print a summary.

Other parameters:
t=auto
pigz=f              Use pigz to compress.  If argument is a number, that will set the number of pigz threads.
unpigz=t            Use pigz to decompress.
passes=2            Generate matrices for 1-pass recalibration only.  Max is 2.
recalqmax=42        Adjust max quality scores tracked.
loadq102=           For each recalibration matrix, enable or disable that matrix with t/f.
                    You can specify pass1 or pass2 like this: loadq102_p1=f loadq102_p2=t.

Java Parameters:
-Xmx                This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

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
}
calcXmx "$@"

calctruequality() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.7_64bit
		module load pigz
		module load samtools
	fi
	local CMD="java $EA $z -cp $CP jgi.CalcTrueQuality $@"
	echo $CMD >&2
	eval $CMD
}

calctruequality "$@"
