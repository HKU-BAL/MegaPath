#include <iostream>
#include <string.h>
#include <cstring>
#include <string>
#include <typeinfo>
#include <pthread.h>
#include <unistd.h>
#include <queue>
#include <vector>
#include <fstream>
#include <sstream>
#include <time.h>
#include <math.h>
#include <inttypes.h>
#include <zlib.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <stdio.h>
#include <set>
using namespace std;

void show_help_and_die(char* my_name) 
{
	fprintf(stderr, "Usage: %s <input genome> <input genomecov> <max depth stdev e.g. 100> \n", my_name);
	exit(1);
}

int main(int argc, char** argv)
{
	if(argc < 4) show_help_and_die(argv[0]);
	int maxDepthStdev = stoi(argv[3]);

	//read genome
	map<string, int> sequence2Length;
	map<string, double> sequence2MovingAverage;
	map<string, double> sequence2MovingVariance; //population
	map<string, int> sequence2Count;
	map<string, double> sequence2MaxDepth;
	map<string, double> sequence2DiffPower;
	FILE* GENOME = fopen(argv[1], "r");
	char sequenceId[256];
	int sequenceLength;
	while(fscanf(GENOME, "%s %d", sequenceId, &sequenceLength) == 2)
	{
		sequence2Length.insert(pair<string, int>(string(sequenceId), sequenceLength));
		sequence2MovingAverage.insert(pair<string, double>(string(sequenceId), 0.0));
		sequence2MovingVariance.insert(pair<string, double>(string(sequenceId), 0.0));
		sequence2DiffPower.insert(pair<string, double>(string(sequenceId), 0.0));
		sequence2Count.insert(pair<string, int>(string(sequenceId), 0));
	}
	fclose(GENOME);
	
	//read genomeCoverage to calculate variance
	FILE* GENOMECOV = fopen(argv[2], "r");
	int sequenceFrom, sequenceTo, depth;
	while(fscanf(GENOMECOV, "%s %d %d %d", sequenceId, &sequenceFrom, &sequenceTo, &depth) == 4)
	{
		if(sequence2Length.find(string(sequenceId)) != sequence2Length.end()) //A original B newSet
		{
			double averageDiff = (double)depth - sequence2MovingAverage[string(sequenceId)];
			double newMovingAverage = (double) sequence2MovingAverage[string(sequenceId)] + averageDiff * (sequenceTo - sequenceFrom) / (sequence2Count[string(sequenceId)] + (sequenceTo - sequenceFrom));
			double newDiffPower = sequence2DiffPower[string(sequenceId)] + pow(averageDiff, 2.0) * (sequenceTo - sequenceFrom) * sequence2Count[string(sequenceId)] / (sequence2Count[string(sequenceId)] + (sequenceTo - sequenceFrom));
			double newMovingVariance = newDiffPower / (sequence2Count[string(sequenceId)] + (sequenceTo - sequenceFrom));
			sequence2DiffPower[string(sequenceId)] = newDiffPower;
			sequence2Count[string(sequenceId)] = sequence2Count[string(sequenceId)] + (sequenceTo - sequenceFrom);
			sequence2MovingAverage[string(sequenceId)] = newMovingAverage;
			sequence2MovingVariance[string(sequenceId)] = newMovingVariance;
		}
	}
	fclose(GENOMECOV);

	double maxDepth;
	for(map<string, double>::iterator it = sequence2MovingAverage.begin(); it != sequence2MovingAverage.end(); it++)
	{
		maxDepth = it->second + maxDepthStdev * sqrt(sequence2MovingVariance[it->first]);
		sequence2MaxDepth.insert(pair<string, double>(it->first, maxDepth));
	}

	GENOMECOV = fopen(argv[2], "r");
	while(fscanf(GENOMECOV, "%s %d %d %d", sequenceId, &sequenceFrom, &sequenceTo, &depth) == 4)
	{
		if(sequence2MaxDepth.find(string(sequenceId)) != sequence2MaxDepth.end())
		{
			if(depth > sequence2MaxDepth[string(sequenceId)])
				fprintf(stdout, "%s\t%d\t%d\n", sequenceId, sequenceFrom, sequenceTo);
		}
	}
	fclose(GENOMECOV);
	return 0;
}
