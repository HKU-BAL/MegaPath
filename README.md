[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# MegaPath
Sensitive and rapid pathogen detection using metagenomic NGS data; MegaPath-Amplicon: filtering module for metagenomic amplicon data.

## Prerequisites

Requirements:

MegaPath

Memory: 100G
Storage: 250G

MegaPath-Amplicon

Memory: 400G
Storage: 450G

## Option 1: Bioconda
```
# prioritize channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mp -c bioconda megapath
conda activate mp
```

## Option 2: Conda Virtual Environment Setup
```
# prioritize channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n mp python=3.6.10
conda activate mp
# MegaPath
conda install samtools==1.10 bedtools==2.27.1 megahit==1.1.3
# MegaPath-Amplicon
conda install gatk4 pandas pysam pysam=0.16.0.1 bwa=0.7.12 pypy3.6 parallel=20191122 seqtk

# git clone MegaPath
git clone --depth 1 https://github.com/edwwlui/MegaPath

# MegaPath
make -C MegaPath/cc/
make -C MegaPath/soap4/2bwt-lib/
make -C MegaPath/soap4/

# MegaPath-Amplicon
cd MegaPath/scripts/realignment/realign/
g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp
g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp 
gcc -Wall -O3 -pipe -fPIC -shared -rdynamic -o libssw.so ssw.c ssw.h
```


## Download pre-built databases of MegaPath and/or MegaPath-Amplicon
```
# Option 1, Bioconda: cd ${CONDA_PREFIX}/MegaPath
# conda info --env can show the ${CONDA_PREFIX} in the current environment.
# Option 2, Conda virtual env: cd ./MegaPath (the git clone)
cd ${MEGAPATH_DIR}
# MegaPath db
wget -c http://www.bio8.cs.hku.hk/dataset/MegaPath/MegaPath_db.v1.0.tar.gz -O - | tar -xvz
# MegaPath-Amplicon db
wget -c http://www.bio8.cs.hku.hk/dataset/MegaPath/MegaPath-Amplicon_db.v1.0.tar.gz -O - | tar -xvz

```

## Basic usage
```
Usage: ./runMegaPath.sh -1 <read1.fq> -2 <read2.fq> [options]
    -p  output prefix [megapath]
    -t  number of threads [24]
    -c  NT alignment score cutoff [40]
    -s  SPIKE filter number of stdev [60]
    -o  SPIKE overlap [0.5]
    -L  max read length [150]
    -d  database directory
    -S  skip ribosome filtering
    -H  skip human filtering
    -A  Perform assembly & protein alignment

Usage: ./runMegaPath-Amplicon.sh -1 <read1.fq> -2 <read2.fq> [options]
    -p  output prefix [megapath-amplicon]
    -t  number of threads [24]
    -L  max read length [250]
    -d  database directory [/autofs/bal13/wwlui/bioconda/test_mp/test/MegaPath/db]
```

## Advanced usage for MegaPath-Amplicon
```
Keep BWA index in memory to significantly speed up batch-run
# load indices
bwa shm ${CONDA_PREFIX}/MegaPath/db/amplicon/GCF_000001405.39_GRCh38.p13_genomic.fna.gz &
bwa shm ${CONDA_PREFIX}/MegaPath/db/amplicon/Mycobacterium_tuberculosis_H37Rv_genome_v3.fasta &
bwa shm ${CONDA_PREFIX}/MegaPath/db/amplicon/ORAL_MICROBIOME_genomic.no_Myc.fna &
bwa shm ${CONDA_PREFIX}/MegaPath/db/amplicon/refseq.fna.no_blastid85_tb.nomyc.onlytb.gz &
# unload all indices
bwa shm -d
```
