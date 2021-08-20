FROM ubuntu

RUN apt-get update && \
 apt-get upgrade -y && \
 apt-get install -y wget build-essential git curl make vim python libz-dev libncurses5-dev libncursesw5-dev

RUN curl -L http://xrl.us/installperlnix | bash

COPY . /megapath
WORKDIR /megapath


RUN cd cc/ && make
RUN cd bedtools2/ && make && cd /megapath/
RUN cd megahit/ && make && cd /megapath/
RUN cd soap4/2bwt-lib/ && make && cd ../ && make && cd /megapath/
RUN cd samtools-0.1.18/ && make && cd /megapath/

#Download MegaPath databases
#wget http://www.bio8.cs.hku.hk/dataset/MegaPath/MegaPath_db.v1.0.tar.gz
#tar -xvzf MegaPath_db.v1.0.tar.gz
