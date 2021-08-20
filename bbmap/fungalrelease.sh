#!/bin/bash
#fungalrelease in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified December 10, 2015

Description:  Reformats a fungal assembly for release.
Also creates contig and agp files.

Usage:  fungalrelease.sh in=<input file> out=<output file>


File parameters:
in=<file>           Input scaffolds.
out=<file>          Output scaffolds.
outc=<file>         Output contigs.
qfin=<file>         Optional quality scores input.
qfout=<file>        Optional quality scores output.
qfoutc=<file>       Optional contig quality scores output.
agp=<file>          Output AGP file.
legend=<file>       Output name legend file.
overwrite=f         (ow) Set to false to force the program to abort rather than
                    overwrite an existing file.

Processing parameters:
fastawrap=60        Wrap length for fasta lines.
tuc=t               Convert sequence to upper case.
baniupac=t          Crash on encountering a non-ACGTN base call.
mingap=10           Expand all gaps (Ns) to be at least this long.
mingapin=1          Only expand gaps that are at least this long.
sortcscaffolds=t    Sort scaffolds descending by length.
sortcontigs=f       Sort contigs descending by length.
renamescaffolds=t   Rename scaffolds to 'scaffold_#'.
scafnum=1           Number of first scaffold.
renamecontigs=f     Rename contigs to 'contig_#' instead of 'scafname_c#'.
contignum=1         Number of first contig; only used if renamecontigs=t.
minscaf=1           Only retain scaffolds at least this long.
mincontig=1         Only retain contigs at least this long.


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

z="-Xmx4g"
z2="-Xms4g"
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

fungalrelease() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.7_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP jgi.FungalRelease $@"
	echo $CMD >&2
	eval $CMD
}

fungalrelease "$@"
