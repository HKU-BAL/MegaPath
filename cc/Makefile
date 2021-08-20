CXX = g++
CXXFLAGS = -std=c++11 -g -O2 -Wall -std=c++0x -fopenmp -lz

all:
	$(CXX) calculate_m8_cov.cpp $(CXXFLAGS) -o calculate_m8_cov
	$(CXX) cleanup.cpp $(CXXFLAGS) -o cleanup
	$(CXX) createDB.cpp $(CXXFLAGS) -o createDB
	$(CXX) deinterleave.cpp $(CXXFLAGS) -o deinterleave
	$(CXX) fastq2lsam.cpp $(CXXFLAGS) -o fastq2lsam
	$(CXX) filterCrossFamilyReads.cpp $(CXXFLAGS) -o filterCrossFamilyReads
	$(CXX) filterDB.cpp $(CXXFLAGS) -o filterDB
	$(CXX) filterExistingSeq.cpp $(CXXFLAGS) -o filterExistingSeq
	$(CXX) genCountTable.cpp $(CXXFLAGS) -o genCountTable
	$(CXX) genKrakenReport.cpp $(CXXFLAGS) -o genKrakenReport
	$(CXX) genomeCovFilter.cpp $(CXXFLAGS) -o genomeCovFilter
	$(CXX) kraken-censtruct.cpp $(CXXFLAGS) -o kraken-censtruct
	$(CXX) lsamReadFilter.cpp $(CXXFLAGS) -o lsamReadFilter
	$(CXX) m8_to_mapLen_hist.cpp $(CXXFLAGS) -o m8_to_mapLen_hist
	$(CXX) maskLowerWithN.cpp $(CXXFLAGS) -o maskLowerWithN
	$(CXX) mask_with_N.cpp $(CXXFLAGS) -o mask_with_N
	$(CXX) masonAccuracy.cpp $(CXXFLAGS) -o masonAccuracy
	$(CXX) prepareKrakenFasta.cpp $(CXXFLAGS) -o prepareKrakenFasta
	$(CXX) reassign.cpp $(CXXFLAGS) -o reassign
	$(CXX) revCompFastx.cpp $(CXXFLAGS) -o revCompFastx
	$(CXX) ribosomeDedup.cpp $(CXXFLAGS) -o ribosomeDedup
	$(CXX) sam2cfq.cpp $(CXXFLAGS) -o sam2cfq
	$(CXX) selectSameSpecieGenome.cc $(CXXFLAGS) -o selectSameSpecieGenome
	$(CXX) smart_pairing.cpp $(CXXFLAGS) -o smart_pairing
	$(CXX) split_ref_to_reads.cpp $(CXXFLAGS) -o split_ref_to_reads
	$(CXX) surpiAnn2id.cpp $(CXXFLAGS) -o surpiAnn2id
	$(CXX) taxLookupAcc.cpp $(CXXFLAGS) -o taxLookupAcc

#%: %.cpp %.h $(DEPS)
#	$(CXX) $< -o $@ $(CXXFLAGS)

