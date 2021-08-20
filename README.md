[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

# MegaPath

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

## Option 2: Docker
```
sudo docker build -f ./Dockerfile -t mp_image . 
sudo docker run -it mp_image /bin/bash
```



## Pre-built Database Download
```
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
