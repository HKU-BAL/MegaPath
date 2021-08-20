#
#
#    Makefile
#    Soap3(gpu)
#
#    Copyright (C) 2011, HKU
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 2
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
#
####################################################################
#
#   Modification History
#
#   Date   : 6th January 2012
#   Author : Edward MK Wu
#   Change : Standardising SOAP4 package
#
####################################################################

CXX = g++
CC = g++
LIBFLAG = -lz -lm
CXXFLAGS = -w -O3 -funroll-loops -march=native -maccumulate-outgoing-args -Wno-unused-result -static-libgcc -mavx -fopenmp -std=c++0x -fomit-frame-pointer 
#CFLAGS = -w -O3 -funroll-loops -march=native -maccumulate-outgoing-args -Wno-unused-result -static-libgcc -mavx -fopenmp -std=c++0x -fomit-frame-pointer 

BWTLIB = 2bwt-lib
BWTOBJLIBS = $(BWTLIB)/BWT.o $(BWTLIB)/dictionary.o $(BWTLIB)/DNACount.o $(BWTLIB)/HSP.o $(BWTLIB)/HSPstatistic.o $(BWTLIB)/iniparser.o $(BWTLIB)/inistrlib.o $(BWTLIB)/karlin.o $(BWTLIB)/MemManager.o $(BWTLIB)/MiscUtilities.o $(BWTLIB)/QSufSort.o $(BWTLIB)/r250.o $(BWTLIB)/TextConverter.o $(BWTLIB)/Timing.o $(BWTLIB)/Socket.o

SAMLIB = samtools-0.1.18
SAMOBJLIBS = $(SAMLIB)/sam.o $(SAMLIB)/bam.o $(SAMLIB)/bgzf.o $(SAMLIB)/kstring.o $(SAMLIB)/bam_import.o $(SAMLIB)/faidx.o $(SAMLIB)/bam_pileup.o $(SAMLIB)/bam_aux.o $(SAMLIB)/sam_header.o $(SAMLIB)/razf.o

CPULIB = 2bwt-flex
CPUOBJLIB = $(CPULIB)/LT.o


all:	SOAP4
# Application
#
VERSION.h: force_look
	git rev-parse HEAD | awk '{print "#ifndef __VERSION_H__\n#define __VERSION_H__\n#define COMMIT_HASH \"" $$0 "\"\n\n#endif" }' > VERSION.h

force_look: 
	true

SOAP4:	SOAP4.o BGS-IO.o PE.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o CPU_DPfunctions.o CPU_DP.o DV-DPForBothUnalign.o DV-DPForSingleReads.o aio_thread.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o SeedPool.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB) VERSION.h
	$(CC) $(CXXFLAGS) SOAP4.o BGS-IO.o PE.o CPUfunctions.o alignment.o PEAlgnmt.o SAM.o DV-SemiDP.o OutputDPResult.o AlgnResult.o DV-DPfunctions.o CPU_DPfunctions.o CPU_DP.o DV-DPForBothUnalign.o DV-DPForSingleReads.o aio_thread.o IniParam.o UsageInterface.o FileUtilities.o IndexHandler.o QueryParser.o SeedPool.o $(SAMOBJLIBS) $(BWTOBJLIBS) $(CPUOBJLIB) $(STATICLIBZ) -o soap4 $(LIBFLAG)


clean:
	-rm -f .*.cpp *.o $(BWTOBJLIBS) $(CPUOBJLIB) || true
	-rm soap4 soap4-builder || true
