[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# MegaPath
Sensitive and rapid pathogen detection using metagenomic NGS data

## Prerequisites

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
conda install samtools==0.1.18 bedtools==2.27.1

# git clone MegaPath
git clone --depth 1 https://github.com/edwwlui/MegaPath

make -C MegaPath/megahit/
make -C MegaPath/soap4/2bwt-lib/
make -C MegaPath/soap4/
```


## Pre-built Database Download
```
# Option 1, Bioconda: cd ${CONDA_PREFIX}/MegaPath
# conda info --env can show the ${CONDA_PREFIX} in the current environment.
# Option 2, Conda virtual env: cd ./MegaPath (the git clone)
cd ${MEGAPATH_DIR}
wget http://www.bio8.cs.hku.hk/dataset/MegaPath/MegaPath_db.v1.0.tar.gz
tar -xvzf MegaPath_db.v1.0.tar.gz
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
```
