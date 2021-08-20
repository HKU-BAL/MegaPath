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
	fprintf(stderr, "Usage: %s <input readFilter> <input lsam.id>\n", my_name);
	exit(1);
}

int main(int argc, char** argv)
{
	if(argc < 3) show_help_and_die(argv[0]);

	//read Filtered readId
	set<string> readIdSet; //should be filtered
	FILE* READFILTER = fopen(argv[1], "r");
	char readId[256];
	while(fscanf(READFILTER, "%s", readId) == 1)
	{
		readIdSet.insert(string(readId));
		if(string(readId)[string(readId).length()-2] == '/')
			readIdSet.insert(string(readId).substr(0, string(readId).length()-2));
	}
	fclose(READFILTER);
	
	//read lsam.id file and do filtering
	FILE* LSAM = fopen(argv[2], "r");
	char *ptr = NULL;
	size_t len;
	istringstream iss;
	string line, readIdS;
	while(getline(&ptr, &len, LSAM) != -1)
	{
		iss.str(ptr);
		iss >> readIdS;
		if(readIdSet.find(readIdS) == readIdSet.end())
		{
			fprintf(stdout, "%s", ptr);
		}
	}
	fclose(LSAM);
	return 0;
}
