#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <assert.h>
#include <unistd.h>
#include <math.h>
#include <unordered_map>
#include "kxseq.h"
#include "taxonomy.h"
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <unordered_set>
#include <map>

using namespace std;
using namespace kx;

TaxDB db;
char delimiters[5]="\t";
int taxIdIndex=4, salignedIndex=8;
unordered_map<int, int> taxid2clade;
unordered_map<int, int> taxid2aligned;
unordered_map<int, unordered_set<int>> taxid2child;
const static vector<string> Ranks = {"domain", "kingdom", "phylum", "class", "order", "family", "genus", "species"};
long long int totAligned=0;
void printIndent(int d)
{
	while(d--)
		printf("  ");
}

char getLevelCode(int tid)
{
	string r=db.taxInfo[tid].rank.c_str();
	if (r=="superkingdom") 
		return 'D';
	else if (find(Ranks.begin(),Ranks.end(), r) != Ranks.end()) 
		return r[0]+'A'-'a';
	else 
		return '-';
}

bool cmp(int a, int b)
{
	return taxid2clade[a]>taxid2clade[b];
}

void printOutcome(int p, int depth)
{
	//prec
	double prec=(double)(taxid2clade[p]*100)/totAligned;
	printf("%.2lf\t", prec);
	//n-clade
	printf("%d\t", taxid2clade[p]);
	//n-stay
	if(taxid2aligned.find(p)!=taxid2aligned.end())
		printf("%d\t", taxid2aligned[p]);
	else
		printf("0\t");
	//level
	printf("%c\t", getLevelCode(p));
	//taxonid
	printf("%d\t", p);
	//depth
	printf("%d\t", depth-1);
	//indent and name
	printIndent(depth);
	if(p==0)
		printf("unclassified\n");
	else
		printf("%s\n", db.taxInfo[p].name.c_str());

	
	if(taxid2child.find(p) != taxid2child.end())
	{
		vector<int> child;
		for(unordered_set<int>::iterator it=taxid2child[p].begin();it!=taxid2child[p].end();it++)
		{
			child.push_back(*it);
		}
		sort(child.begin(), child.end(), cmp);
		for(int i=0;i<child.size();i++)
			printOutcome(child[i], depth+1);
	}
}

void read(int tid, int aligned)
{	
	if(taxid2clade.find(tid)==taxid2clade.end())
		taxid2clade[tid]=aligned;
	else
		taxid2clade[tid]+=aligned;
	
	while(tid!=0 && tid!=1)
	{
		int p=db.taxInfo[tid].parent;
		taxid2child[p].insert(tid);
		tid=p;
		if(taxid2clade.find(tid)==taxid2clade.end())
		{
			taxid2clade[tid]=aligned;
		}
		else
		{
			taxid2clade[tid]+=aligned;
		}
		
	}
}

int main(int argc, char **argv)
{
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <nodes.dmp> <names.dmp> <input.japsa> taxidIndex=4 alignedIndex=8\n", argv[0]);
		exit(1);
	}

	db.readNode(argv[1]);
	db.readName(argv[2]);
	ifstream input(argv[3]);
	if(argc >= 5)
		taxIdIndex=stoi(argv[4]);
	if(argc>=6)
		salignedIndex=stoi(argv[5]);
	char str[200];
	char *ptr;
	input.getline(str, 200);
	totAligned=0;
	while(input.getline(str, 200))
	{
		ptr=strtok(str, delimiters);
		int cindex=0;
		int ctid, caligned;
		while(ptr){
			if(cindex==taxIdIndex)
				ctid=stoi(ptr);
			if(cindex==salignedIndex)
				caligned=stoi(ptr);
			ptr=strtok(NULL, delimiters);
			cindex++;
		}
		taxid2aligned[ctid]=caligned;
		totAligned+=caligned;
		read(ctid, caligned);
	}
	printf("prec\tn-clade\tn-stay\tlevel\ttaxonid\tdepth\tname\n");
	printOutcome(0, 1);
	printOutcome(1, 1);
	return 0;
}

