/*

   MemManager.c		Memory Manager

   This module provides memory management functions.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _WIN32
#include <mm_malloc.h>
#endif
#include "MiscUtilities.h"
#include "MemManager.h"

MMMaster mmMaster;

void *MMMalloc(const long long memSize) {

	void *address;

	address = MEMALIGN(memSize, MAX_ALIGN);
	if (address == NULL) {
		fprintf(stderr, "MMMalloc() : cannot allocate memory!\n");
		exit(1);
	}
	return address;

}

void MMFree(void *address) {

	FREEALIGN(address);

}

void MMMasterInitialize(const unsigned long long maxNumberOfPools, const unsigned long long maxNumberOfBulks, 
						const long long traceUnitByteAllocation, FILE* unitByteTraceFile) {

	unsigned long long i;

	mmMaster.maxTotalByteAllocated = 0;
	mmMaster.maxTotalByteDispatched = 0;
	mmMaster.currentUnitByteAllocated = 0;
	mmMaster.maxUnitByteAllocated = 0;

	mmMaster.maxNumberOfBulks = maxNumberOfBulks;
	mmMaster.maxNumberOfPools = maxNumberOfPools;
	if (maxNumberOfBulks > 0) {
		mmMaster.mmBulk = (MMBulk**) MEMALIGN(sizeof(MMBulk*) * maxNumberOfBulks, MAX_ALIGN);
		for (i=0; i<maxNumberOfBulks; i++) {
			mmMaster.mmBulk[i] = NULL;
		}
	} else {
		mmMaster.mmBulk = NULL;
	}
	if (maxNumberOfPools > 0) {
		mmMaster.mmPool = (MMPool**) MEMALIGN(sizeof(MMPool*) * maxNumberOfPools, MAX_ALIGN);
		for (i=0; i<maxNumberOfPools; i++) {
			mmMaster.mmPool[i] = NULL;
		}
	} else {
		mmMaster.mmPool = NULL;
	}

	mmMaster.traceUnitByteAllocation = traceUnitByteAllocation;
	mmMaster.unitByteTraceFile = unitByteTraceFile;

}

void MMMasterFreeAll() {

	unsigned long long i;

	for (i=0; i < mmMaster.maxNumberOfBulks; i++) {
		if (mmMaster.mmBulk[i] != NULL) {
			if (MMBulkIsActive(mmMaster.mmBulk[i])) {
				MMBulkFree(mmMaster.mmBulk[i]);
			}
			if (MMBulkFindPoolUsed(mmMaster.mmBulk[i]) == NULL) {
				MMUnitFree(mmMaster.mmBulk[i], sizeof(MMBulk));
			}
			mmMaster.mmBulk[i] = NULL;
		}
	}
	FREEALIGN(mmMaster.mmBulk);

	for (i=0; i < mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] != NULL) {
			if (MMPoolIsActive(mmMaster.mmPool[i])) {
				MMPoolFree(mmMaster.mmPool[i]);
			}
			FREEALIGN(mmMaster.mmPool[i]);
			mmMaster.mmPool[i] = NULL;
		}
	}
	FREEALIGN(mmMaster.mmPool);

}

long long MMMasterCurrentTotalByteAllocated() {

	unsigned long long i;
	long long currentTotalByteAllocated;

	// unit memory allocated
	currentTotalByteAllocated = mmMaster.currentUnitByteAllocated;

	// pool and temp memory allocated
	for (i=0; i < mmMaster.maxNumberOfPools; i++) {
        if (mmMaster.mmPool[i] != NULL && MMPoolIsActive(mmMaster.mmPool[i])) {
			currentTotalByteAllocated += MMPoolCurrentTotalByteAllocated(mmMaster.mmPool[i]);
		}
	}

	// bulk memory allocated
	for (i=0; i < mmMaster.maxNumberOfBulks; i++) {
        if (mmMaster.mmBulk[i] != NULL && MMBulkIsActive(mmMaster.mmBulk[i])) {
			currentTotalByteAllocated += MMBulkByteAllocated(mmMaster.mmBulk[i]);
		}
	}

	return currentTotalByteAllocated;

}

long long MMMasterCurrentTotalByteDispatched() {

	unsigned long long i;
	long long currentTotalByteDispatched;

	// unit memory dispatched
	currentTotalByteDispatched = mmMaster.currentUnitByteAllocated;

	// pool and temp memory dispatched
	for (i=0; i < mmMaster.maxNumberOfPools; i++) {
        if (mmMaster.mmPool[i] != NULL && MMPoolIsActive(mmMaster.mmPool[i])) {
			currentTotalByteDispatched += MMPoolCurrentTotalByteDispatched(mmMaster.mmPool[i]);
		}
	}

	// bulk memory dispatched
	for (i=0; i < mmMaster.maxNumberOfBulks; i++) {
        if (mmMaster.mmBulk[i] != NULL && MMBulkIsActive(mmMaster.mmBulk[i])) {
			currentTotalByteDispatched += MMBulkByteDispatched(mmMaster.mmBulk[i]);
		}
	}

	return currentTotalByteDispatched;

}

long long MMMasterMaxTotalByteAllocated() {

	long long currentTotalByteAllocated;

	currentTotalByteAllocated = MMMasterCurrentTotalByteAllocated();

	if (currentTotalByteAllocated > mmMaster.maxTotalByteAllocated) {
		return currentTotalByteAllocated;
	} else {
		return mmMaster.maxTotalByteAllocated;
	}

}

long long MMMasterMaxTotalByteDispatched() {

	long long currentTotalByteDispatched ;

	currentTotalByteDispatched = MMMasterCurrentTotalByteDispatched();

	if (currentTotalByteDispatched > mmMaster.maxTotalByteDispatched) {
		return currentTotalByteDispatched;
	} else {
		return mmMaster.maxTotalByteDispatched;
	}

}

void MMMasterSetMaxTotalByteAllocated() {

	long long currentTotalByteAllocated;
	
	currentTotalByteAllocated = MMMasterCurrentTotalByteAllocated();

	if (currentTotalByteAllocated > mmMaster.maxTotalByteAllocated) {
		mmMaster.maxTotalByteAllocated = currentTotalByteAllocated;
	}

}

void MMMasterSetMaxTotalByteDispatched() {

	long long currentTotalByteDispatched;
	
	currentTotalByteDispatched = MMMasterCurrentTotalByteDispatched();

	if (currentTotalByteDispatched > mmMaster.maxTotalByteDispatched) {
		mmMaster.maxTotalByteDispatched = currentTotalByteDispatched;
	}

}

void MMMasterPrintReport(FILE *output, const unsigned long long withUnitDetails, const unsigned long long withPoolDetails, const unsigned long long withBulkDetails) {

	unsigned long long i;

	fprintf(output, "Maximum amount of memory allocated:  %u\n", (unsigned int) MMMasterMaxTotalByteAllocated());
	fprintf(output, "Maximum amount of memory dispatched: %u\n", (unsigned int) MMMasterMaxTotalByteDispatched());

	if (withUnitDetails) {
		fprintf(output, "\n");
		MMUnitPrintReport(output);
	}

	if (withPoolDetails) {
		for (i=0; i<mmMaster.maxNumberOfPools; i++) {
			if (mmMaster.mmPool[i] != NULL) {
				fprintf(output, "\nPool number %llu\n", i);
				MMPoolPrintReport(mmMaster.mmPool[i], output);
			}
		}
	}

	if (withBulkDetails) {
		for (i=0; i<mmMaster.maxNumberOfBulks; i++) {
			if (mmMaster.mmBulk[i] != NULL) {
				fprintf(output, "\nBulk number %llu\n", i);
				MMBulkPrintReport(mmMaster.mmBulk[i], output);
			}
		}
	}

}

void *MMUnitAllocate(const long long memSize) {

	void *temp;

	#ifdef DEBUG
	if (memSize == 0) {
		fprintf(stderr, "MMUnitAllocate() : memSize = 0!\n");
		exit(1);
	}
	#endif

	temp = MEMALIGN(memSize, MAX_ALIGN);
	if (temp == NULL) {
		fprintf(stderr, "MMUnitAllocate() : cannot allocate memory!\n");
		exit(1);
	}

	mmMaster.currentUnitByteAllocated += memSize;
	if (mmMaster.traceUnitByteAllocation) {
		fprintf(mmMaster.unitByteTraceFile, "MMUnitAllocate        : %llu\n", memSize);
	}

	return temp;

}

void *MMUnitReallocate(void *address, const long long newMemSize, const long long oldMemSize) {

	void *temp;

	#ifdef DEBUG
	if (newMemSize == 0) {
		fprintf(stderr, "MMUnitReallocate() : newMemSize = 0!\n");
		exit(1);
	}
	if (oldMemSize == 0) {
		fprintf(stderr, "MMUnitReallocate() : oldMemSize = 0!\n");
		exit(1);
	}
	#endif

	if (mmMaster.traceUnitByteAllocation) {
		fprintf(mmMaster.unitByteTraceFile, "MMUnitReallocate\n");
	}

	temp = MMUnitAllocate(newMemSize);
	if (temp == NULL) {
		fprintf(stderr, "MMUnitReallocate() : cannot allocate memory!\n");
		exit(1);
	}
	memcpy(temp, address, min(newMemSize, oldMemSize));

	MMUnitFree(address, oldMemSize);

	return temp;

}

void MMUnitFree(void *address, const long long memSize) {

	#ifdef DEBUG
	if (address == NULL) {
		fprintf(stderr, "MMUnitFree() : address = NULL!\n");
		exit(1);
	}
	if (mmMaster.currentUnitByteAllocated < memSize) {
		fprintf(stderr, "MMUnitFree() : currentUnitByteAllocated < memSize!\n");
		exit(1);
	}
	#endif

	FREEALIGN(address);

	#ifdef RECORD_GRAND_TOTAL
	MMMasterSetMaxTotalByteAllocated();
	MMMasterSetMaxTotalByteDispatched();
	#endif

	if (mmMaster.currentUnitByteAllocated > mmMaster.maxUnitByteAllocated) {
		mmMaster.maxUnitByteAllocated = mmMaster.currentUnitByteAllocated;
	}
	mmMaster.currentUnitByteAllocated -= memSize;
	if (mmMaster.traceUnitByteAllocation) {
		fprintf(mmMaster.unitByteTraceFile, "MMUnitFree            : %llu\n", memSize);
	}

}

long long MMUnitCurrentByteAllocated() {

	return mmMaster.currentUnitByteAllocated;

}

long long MMUnitMaxByteAllocated() {

	if (mmMaster.currentUnitByteAllocated > mmMaster.maxUnitByteAllocated) {
		return mmMaster.currentUnitByteAllocated;
	} else {
		return mmMaster.maxUnitByteAllocated;
	}

}

void MMUnitPrintReport(FILE *output) {

	fprintf(output, "Maximum amount of unit memory allocated: %llu\n", MMUnitMaxByteAllocated());
	fprintf(output, "Amount of memory unit memory currently allocated: %llu\n", MMUnitCurrentByteAllocated());

}

MMPool *MMPoolCreate(const unsigned long long poolSize) {

	MMPool *mmPool;
	unsigned long long i;

	#ifdef DEBUG
	if (poolSize < sizeof(MMPool)) {
		fprintf(stderr, "MMPoolCreate() : poolSize < MMPool!\n");
		exit(1);
	}
	#endif

	if (poolSize / MAX_ALIGN * MAX_ALIGN != poolSize) {
		fprintf(stderr, "MMPoolCreate() : poolSize must be multiple of MAX_ALIGN (%d)!\n", MAX_ALIGN);	// Otherwise temp memory is not properly aligned
		exit(1);
	}

	mmPool = (MMPool*) MEMALIGN(poolSize, MAX_ALIGN);
	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolCreate() : cannot allocate memory!\n");
		exit(1);
	}

	mmPool->poolSize = poolSize;
	mmPool->poolByteDispatched = sizeof(MMPool);
	mmPool->poolByteSpillover = 0;
	mmPool->firstSpillOverAddress = NULL;
	mmPool->currentTempByteDispatched = 0;
	mmPool->currentTempByteSpillover = 0;
	mmPool->maxTotalByteDispatched = 0;

	for (i=0; i<mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] == NULL) {
			mmMaster.mmPool[i] = mmPool;
			return mmPool;
		}
	}

	fprintf(stderr, "MMPoolCreate() : number of pools > maxNumberOfPools!\n");
	exit(1);

}

unsigned int MMPoolIsActive(const MMPool *mmPool) {

	return ((mmPool->firstSpillOverAddress) != (void*)mmPool);

}
void MMPoolSetInactive(MMPool *mmPool) {

	if (mmPool->firstSpillOverAddress != NULL) {
		fprintf(stderr, "MMPoolSetInactive() : spillover memory not freed yet!\n");
		exit(1);
	}

	mmPool->firstSpillOverAddress = (void*)mmPool;
}

unsigned long long MMPoolCurrentTotalByteAllocated(const MMPool *mmPool) {

	return mmPool->poolSize + mmPool->poolByteSpillover + mmPool->currentTempByteSpillover;

}

unsigned long long MMPoolCurrentTotalByteDispatched(const MMPool *mmPool) {

	return mmPool->poolByteDispatched + mmPool->currentTempByteDispatched;

}

unsigned long long MMPoolMaxTotalByteDispatched(const MMPool *mmPool) {

	unsigned long long currentTotalByteDispatched;

	currentTotalByteDispatched = MMPoolCurrentTotalByteDispatched(mmPool);
	
	if (currentTotalByteDispatched > mmPool->maxTotalByteDispatched) {
		return currentTotalByteDispatched;
	} else {
		return mmPool->maxTotalByteDispatched;
	}

}

unsigned long long MMPoolByteAvailable(const MMPool *mmPool) {

	if (mmPool->poolSize > mmPool->poolByteDispatched + MAX_ALIGN) {
		return (mmPool->poolSize - mmPool->poolByteDispatched + MAX_ALIGN - 1) / MAX_ALIGN * MAX_ALIGN;
	} else {
		return 0;
	}

}

MMPool *MMPoolFree(MMPool *mmPool) {

	MMPool *dummyMMPool;
	unsigned long long i;
	void *temp1, *temp2;

	#ifdef DEBUG
	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolFree(): mmPool = NULL!\n");
		exit(1);
	}
	#endif

	#ifdef RECORD_GRAND_TOTAL
	MMMasterSetMaxTotalByteAllocated();
	MMMasterSetMaxTotalByteDispatched();
	#endif

	dummyMMPool = (MMPool*) MEMALIGN(sizeof(MMPool), MAX_ALIGN);
	if (dummyMMPool == NULL) {
		fprintf(stderr, "MMPoolFree() : cannot allocate memory!\n");
		exit(1);
	}

	// Free spillover memory
	temp1 = mmPool->firstSpillOverAddress;
	while (temp1 != NULL) {
		temp2 = *((void**)temp1);
		FREEALIGN(temp1);
		temp1 = temp2;
	}
	mmPool->firstSpillOverAddress = NULL;

	dummyMMPool->poolByteDispatched = mmPool->poolByteDispatched;
	dummyMMPool->poolByteSpillover = mmPool->poolByteSpillover;
	dummyMMPool->currentTempByteDispatched = mmPool->currentTempByteDispatched;
	dummyMMPool->currentTempByteSpillover = mmPool->currentTempByteSpillover;
	dummyMMPool->firstSpillOverAddress = mmPool->firstSpillOverAddress;
	dummyMMPool->maxTotalByteDispatched = mmPool->maxTotalByteDispatched;
	dummyMMPool->poolSize = mmPool->poolSize;

	MMPoolSetInactive(dummyMMPool);

	// Update master directory
	for (i=0; i<mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] == mmPool) {
			mmMaster.mmPool[i] = dummyMMPool;
			FREEALIGN(mmPool);
			return dummyMMPool;
		}
	}

	fprintf(stderr, "MMPoolFree() : cannot locate pool in master!\n");
	exit(1);

}

void MMPoolReset(MMPool *mmPool) {

	void *temp1, *temp2;

	#ifdef DEBUG
	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolReset(): mmPool = NULL!\n");
		exit(1);
	}
	#endif

	#ifdef RECORD_GRAND_TOTAL
	MMMasterSetMaxTotalByteAllocated();
	MMMasterSetMaxTotalByteDispatched();
	#endif

	// Free spillover memory
	temp1 = mmPool->firstSpillOverAddress;
	while (temp1 != NULL) {
		temp2 = *((void**)temp1);
		FREEALIGN(temp1);
		temp1 = temp2;
	}

	mmPool->poolByteDispatched = sizeof(MMPool);
	mmPool->poolByteSpillover = 0;
	mmPool->currentTempByteDispatched = 0;
	mmPool->currentTempByteSpillover = 0;
	mmPool->firstSpillOverAddress = NULL;
	mmPool->maxTotalByteDispatched = 0;

}

void MMPoolDestory(MMPool *mmPool) {

	unsigned long long i;
	MMPool *temp;

	#ifdef DEBUG
	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolDestory(): mmPool = NULL!\n");
		exit(1);
	}
	#endif

	if (MMPoolIsActive(mmPool)) {
		temp = MMPoolFree(mmPool);
	} else {
		temp = mmPool;
	}

	// Update master directory
	for (i=0; i<mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] == temp) {
			mmMaster.mmPool[i] = NULL;
			FREEALIGN(temp);
			temp = NULL;
		}
	}

	if (temp != NULL) {
		fprintf(stderr, "MMPoolDestory() : cannot locate pool in master!\n");
		exit(1);
	}

}

void *MMPoolDispatch(MMPool *mmPool, const unsigned long long memSize) {

	void **temp;
	unsigned long long totalPoolMemoryUsed, nextPoolMemoryOffset;
	unsigned long long align, skipForAlign;

	if (mmPool == NULL) {
		return MMUnitAllocate(memSize);
	}
	if (memSize == 0) {
		fprintf(stderr, "MMPoolDispatch(): memSize = 0!\n");
		exit(1);
	}

	totalPoolMemoryUsed = mmPool->poolByteDispatched - mmPool->poolByteSpillover +
						  mmPool->currentTempByteDispatched - mmPool->currentTempByteSpillover;
	nextPoolMemoryOffset = mmPool->poolByteDispatched - mmPool->poolByteSpillover;

	// Calculate the number of byte to skip in order to align the memory dispatched
	align = 1 << (BITS_IN_WORD - leadingZero(memSize - 1));
	if (align > MAX_ALIGN) {
		align = MAX_ALIGN;
	}
	if (align < MIN_ALIGN) {
		align = MIN_ALIGN;
	}
	skipForAlign = nextAlignedBoundary(nextPoolMemoryOffset, align) - nextPoolMemoryOffset;

	if (totalPoolMemoryUsed + memSize + skipForAlign <= mmPool->poolSize) {
		temp = (void**)(((char*)mmPool) + nextPoolMemoryOffset + skipForAlign);
		mmPool->poolByteDispatched += memSize + skipForAlign;
		return temp;
	} else {
		// Spillover
		// Allocate for linked list pointer as well
		temp = (void **) MEMALIGN(memSize + MAX_ALIGN, MAX_ALIGN);	// spillover memory is always aligned to MAX_ALIGN
		if (temp == NULL) {
			fprintf(stderr, "MMPoolDispatch(): cannot allocate memory!\n");
			exit(1);
		}
		// Add spillover memory to linked list
		*temp = mmPool->firstSpillOverAddress;
		mmPool->firstSpillOverAddress = temp;
		mmPool->poolByteSpillover += memSize + MAX_ALIGN;
		mmPool->poolByteDispatched += memSize + MAX_ALIGN;
		return (char*)temp + MAX_ALIGN;
	}
		
}

unsigned long long MMPoolDispatchOffset(MMPool *mmPool, const unsigned long long memSize) {

	unsigned long long totalPoolMemoryUsed, nextPoolMemoryOffset;
	unsigned long long align, skipForAlign;

	if (mmPool == NULL) {
		fprintf(stderr, "MMPoolDispatchOffset(): mmPool == NULL!\n");
		exit(1);
	}
	if (memSize == 0) {
		fprintf(stderr, "MMPoolDispatchOffset(): memSize = 0!\n");
		exit(1);
	}

	totalPoolMemoryUsed = mmPool->poolByteDispatched - mmPool->poolByteSpillover +
						  mmPool->currentTempByteDispatched - mmPool->currentTempByteSpillover;
	nextPoolMemoryOffset = mmPool->poolByteDispatched - mmPool->poolByteSpillover;

	// Calculate the number of byte to skip in order to align the memory dispatched
	align = 1 << (BITS_IN_WORD - leadingZero(memSize - 1));
	if (align > MAX_ALIGN) {
		align = MAX_ALIGN;
	}
	if (align < MIN_ALIGN) {
		align = MIN_ALIGN;
	}
	skipForAlign = nextAlignedBoundary(nextPoolMemoryOffset, align) - nextPoolMemoryOffset;

	if (totalPoolMemoryUsed + memSize + skipForAlign > mmPool->poolSize) {
		fprintf(stderr, "MMPoolDispatchOffset(): Not enough memory in memory pool!\n");
		exit(1);
	}

	mmPool->poolByteDispatched += memSize + skipForAlign;

	return nextPoolMemoryOffset + skipForAlign;

}

void MMPoolReturn(MMPool *mmPool, void *address, const unsigned long long memSize) {
	
	if (mmPool == NULL) {
		MMUnitFree(address, memSize);
	}

}

void MMPoolPrintReport(MMPool *mmPool, FILE *output) {

	fprintf(output, "Pool Size     : %llu\n", mmPool->poolSize);
	fprintf(output, "   Dispatched : %llu\n", mmPool->poolByteDispatched);
	fprintf(output, "     - Spillover             : %llu\n", mmPool->poolByteSpillover);
	fprintf(output, "Maximum amount of memory dispatched including temp memory : %llu\n", 
			MMPoolMaxTotalByteDispatched(mmPool));

}

void *MMTempDispatch(MMPool *mmPool, const unsigned long long memSize) {

	void **temp;
	unsigned long long totalPoolMemoryUsed, nextTempMemoryOffset;
	unsigned long long alignedMemSize;
	void **pointerToLastSpilloverAddress;

	if (mmPool == NULL) {
		return MMUnitAllocate(memSize);
	}
	if (memSize == 0) {
		fprintf(stderr, "MMTempDispatch(): memSize = 0!\n");
		exit(1);
	}

	alignedMemSize = nextAlignedBoundary(memSize, MAX_ALIGN);	// temp memory is always aligned to MAX_ALIGN

	totalPoolMemoryUsed = mmPool->poolByteDispatched - mmPool->poolByteSpillover +
						  mmPool->currentTempByteDispatched - mmPool->currentTempByteSpillover;
	nextTempMemoryOffset = mmPool->currentTempByteDispatched - mmPool->currentTempByteSpillover;

	if (totalPoolMemoryUsed + alignedMemSize <= mmPool->poolSize) {
		temp = (void**)(((char*)mmPool) + mmPool->poolSize - nextTempMemoryOffset - alignedMemSize);
		mmPool->currentTempByteDispatched += alignedMemSize;
		return temp;
	} else {
		// Spillover
		// Locate the last spillover memory
		pointerToLastSpilloverAddress = &(mmPool->firstSpillOverAddress);
		temp = (void**)(*pointerToLastSpilloverAddress);
		while (temp != NULL) {
			pointerToLastSpilloverAddress = temp;
			temp = (void**)*pointerToLastSpilloverAddress;
		}
		// Allocate for linked list pointer as well
		temp = (void **) MEMALIGN(memSize + MAX_ALIGN, MAX_ALIGN);
		if (temp == NULL) {
			fprintf(stderr, "MMTempDispatch(): cannot allocate memory!\n");
			exit(1);
		}
		*pointerToLastSpilloverAddress = temp;
		*temp = NULL;
		mmPool->currentTempByteDispatched += memSize + MAX_ALIGN;
		mmPool->currentTempByteSpillover += memSize + MAX_ALIGN;
		return (char*)temp + MAX_ALIGN;
	}
		
}

void MMTempReturn(MMPool *mmPool, void *address, const unsigned long long memSize) {

	void **temp;
	unsigned long long alignedMemSize;
	void **pointerToLastButOneSpillover;
	void *spilloverPointerAddress;

	if (mmPool == NULL) {
		MMUnitFree(address, memSize);
	} else {

		alignedMemSize = nextAlignedBoundary(memSize, MAX_ALIGN);

		if (address >= (void*)mmPool && address <= (void*)((char*)mmPool + mmPool->poolSize)) {
			// No need to record the global level max memory dispatched/allocated
			// because memory pool is allocated as a whole and fluctuation across pools should not be counted
			if (mmPool->poolByteDispatched + mmPool->currentTempByteDispatched > mmPool->maxTotalByteDispatched) {
				mmPool->maxTotalByteDispatched = mmPool->poolByteDispatched + mmPool->currentTempByteDispatched;
			}
			mmPool->currentTempByteDispatched -= alignedMemSize;
		} else {
			#ifdef RECORD_GRAND_TOTAL
			MMMasterSetMaxTotalByteAllocated();
			MMMasterSetMaxTotalByteDispatched();
			#endif
			// Spillover
			spilloverPointerAddress = (void*)((char*)address - MAX_ALIGN);	// MAX_ALIGN no. of bytes preceding temp address
			// Locate the last spillover memory
			pointerToLastButOneSpillover = &(mmPool->firstSpillOverAddress);
			temp = (void**)(*pointerToLastButOneSpillover);
			while (*temp != NULL) {
				pointerToLastButOneSpillover = temp;
				temp = (void**)*pointerToLastButOneSpillover;
			}
			if (*pointerToLastButOneSpillover != spilloverPointerAddress) {
				fprintf(stderr, "MMTempReturn(): address != lastSpilloverAddress! Last allocated temp memory must be freed first\n");
				exit(1);
			}
			FREEALIGN(spilloverPointerAddress);
			*pointerToLastButOneSpillover = NULL;

			if (mmPool->poolByteDispatched + mmPool->currentTempByteDispatched > mmPool->maxTotalByteDispatched) {
				mmPool->maxTotalByteDispatched = mmPool->poolByteDispatched + mmPool->currentTempByteDispatched;
			}
			mmPool->currentTempByteDispatched -= memSize + MAX_ALIGN;
			mmPool->currentTempByteSpillover -= memSize + MAX_ALIGN;
		}

	}

}

void MMTempPrintReport(MMPool *mmPool, FILE *output) {

	MMPoolPrintReport(mmPool, output);

}

MMBulk *MMBulkCreate(MMPool *mmPool, const unsigned long long itemSize, const unsigned long long itemPerAllocationInPowerOf2, 
					 const unsigned long long boundaryCushionSize, const unsigned long long directorySize) {

	unsigned long long i;
	MMBulk *mmBulk;

	#ifdef DEBUG
	if (itemSize == 0) {
		fprintf(stderr, "MMBulkCreate() : itemSize = 0!\n");
		exit(1);
	}
	if (itemPerAllocationInPowerOf2 >= BITS_IN_WORD) {
		fprintf(stderr, "MMBulkCreate() : itemPerAllocationInPowerOf2 >= BITS_IN_WORD!\n");
		exit(1);
	}
	#endif

	if (mmPool == NULL) {
		mmBulk = (MMBulk*) MMUnitAllocate(sizeof(MMBulk));
	} else {
		mmBulk = (MMBulk*) MMPoolDispatch(mmPool, sizeof(MMBulk));
	}

	mmBulk->itemSize = itemSize;
	mmBulk->itemPerAllocationInPowerOf2 = itemPerAllocationInPowerOf2;
	mmBulk->boundaryCushionSize = boundaryCushionSize;
	mmBulk->indexMask = truncateLeft(ALL_ONE_MASK,  BITS_IN_WORD - itemPerAllocationInPowerOf2);
	mmBulk->currentDirectoryEntry = 0;
	mmBulk->nextUnusedItem = 0;
	mmBulk->directorySize = directorySize;

	if (mmPool == NULL) {
		mmBulk->directory = (unsigned char **) MMUnitAllocate(sizeof(unsigned char*) * directorySize);
	} else {
		mmBulk->directory = (unsigned char **) MMPoolDispatch(mmPool, sizeof(unsigned char*) * directorySize);
	}

	//Allocate memory for the first directory entry
	mmBulk->directory[0] = (unsigned char *) MEMALIGN(boundaryCushionSize * 2 + (itemSize << itemPerAllocationInPowerOf2), MAX_ALIGN);
	if (mmBulk->directory[0] == NULL) {
		fprintf(stderr, "MMBulkCreate() : cannot allocate memory!\n");
		exit(1);
	}

	//Advance the address by boundaryCushionSize
	mmBulk->directory[0] += boundaryCushionSize;

	for (i=0; i<mmMaster.maxNumberOfBulks; i++) {
		if (mmMaster.mmBulk[i] == NULL) {
			mmMaster.mmBulk[i] = mmBulk;
			return mmBulk;
		}
	}

	fprintf(stderr, "MMBulkCreate() : number of bulks > maxNumberOfBulk!\n");
	exit(1);

}

unsigned int MMBulkIsActive(const MMBulk *mmBulk) {

	return (mmBulk->directory != (void*)mmBulk);

}

void MMBulkSetInactive(MMBulk *mmBulk) {

	if (mmBulk->directory != NULL) {
	}
	mmBulk->directory = (unsigned char **) mmBulk;

}

unsigned long long MMBulkByteAllocated(const MMBulk *mmBulk) {

	return (mmBulk->currentDirectoryEntry + 1) *
			(mmBulk->boundaryCushionSize * 2 + (mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2));

}

unsigned long long MMBulkByteDispatched(const MMBulk *mmBulk) {

	return (mmBulk->currentDirectoryEntry) *
			(mmBulk->boundaryCushionSize * 2 + (mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2)) +
			mmBulk->boundaryCushionSize * 2 +
			mmBulk->itemSize * mmBulk->nextUnusedItem;

}

unsigned long long MMBulkUnitDispatched(const MMBulk *mmBulk) {

	return mmBulk->currentDirectoryEntry * (1 << mmBulk->itemPerAllocationInPowerOf2) + mmBulk->nextUnusedItem;

}

void MMBulkFree(MMBulk *mmBulk) {

	unsigned long long i;

	#ifdef RECORD_GRAND_TOTAL
	MMMasterSetMaxTotalByteAllocated();
	MMMasterSetMaxTotalByteDispatched();
	#endif

	for (i=0; i<=mmBulk->currentDirectoryEntry; i++) {
		FREEALIGN(mmBulk->directory[i] - mmBulk->boundaryCushionSize);
	}

	if (MMBulkFindPoolUsed(mmBulk) == NULL) {
        MMUnitFree(mmBulk->directory, sizeof(unsigned char*) * mmBulk->directorySize);
	}

	mmBulk->directory = NULL;

	MMBulkSetInactive(mmBulk);

}

void MMBulkDestory(MMBulk *mmBulk) {

	unsigned long long i;
	MMBulk *temp;

	#ifdef DEBUG
	if (mmBulk == NULL) {
		fprintf(stderr, "MMBulkDestory(): mmBulk = NULL!\n");
		exit(1);
	}
	#endif

	if (MMBulkIsActive(mmBulk)) {
		MMBulkFree(mmBulk);
	}

	temp = mmBulk;

	// Update master directory
	for (i=0; i<mmMaster.maxNumberOfBulks; i++) {
		if (mmMaster.mmBulk[i] == temp) {
			mmMaster.mmBulk[i] = NULL;
			if (MMBulkFindPoolUsed(temp) == NULL) {
				MMUnitFree(temp, sizeof(MMBulk));
			}
			temp = NULL;
		}
	}

	if (temp != NULL) {
		fprintf(stderr, "MMBulkDestory() : cannot locate bulk in master!\n");
		exit(1);
	}

}
unsigned long long MMBulkDispatch(MMBulk *mmBulk) {

	if (mmBulk->nextUnusedItem >> mmBulk->itemPerAllocationInPowerOf2) {
		mmBulk->currentDirectoryEntry++;
		if (mmBulk->currentDirectoryEntry >= mmBulk->directorySize) {
			fprintf(stderr, "MMBulkDispatch() : memory directory size overflow!\n");
			exit(1);
		}
		//Allocate memory for the next directory entry
		mmBulk->directory[mmBulk->currentDirectoryEntry] = (unsigned char *) MEMALIGN(mmBulk->boundaryCushionSize * 2 + (mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2), MAX_ALIGN);
		if (mmBulk->directory[mmBulk->currentDirectoryEntry] == NULL) {
			fprintf(stderr, "MMBulkDispatch() : cannot allocate memory!\n");
			exit(1);
		}
		//Advance the address by boundaryCushionSize
		mmBulk->directory[mmBulk->currentDirectoryEntry] += mmBulk->boundaryCushionSize;
		mmBulk->nextUnusedItem = 0;
	}
	return ((mmBulk->currentDirectoryEntry << mmBulk->itemPerAllocationInPowerOf2) | mmBulk->nextUnusedItem++);

}

void *MMBulkAddress(const MMBulk *mmBulk, const unsigned long long index) {

	#ifdef DEBUG
	if (index >= (((mmBulk->currentDirectoryEntry+1) << mmBulk->itemPerAllocationInPowerOf2) | mmBulk->nextUnusedItem)) {
		fprintf(stderr, "MMBulkAddress() : index out of range!\n");
		exit(1);
	}
	#endif

	return &(mmBulk->directory[index >> mmBulk->itemPerAllocationInPowerOf2][(index & mmBulk->indexMask) * mmBulk->itemSize]);
}

MMPool *MMBulkFindPoolUsed(const MMBulk *mmBulk) {

	unsigned long long i;
	void *temp;

	for (i=0; i<mmMaster.maxNumberOfPools; i++) {
		if (mmMaster.mmPool[i] != NULL) {
			if ((void*)mmBulk >= (void*)mmMaster.mmPool[i] &&
				(void*)mmBulk <= (void*)((char*)mmMaster.mmPool[i] + mmMaster.mmPool[i]->poolSize)) {
				return mmMaster.mmPool[i];
			}
			temp = mmMaster.mmPool[i]->firstSpillOverAddress;
			while (temp != NULL) {
				if ((void*)((char*)temp + sizeof(void*)) == (void*)mmBulk) {
					return mmMaster.mmPool[i];
				}
				temp = *((void**)temp);
			}
		}
	}

	return NULL;

}

void MMBulkPrintReport(MMBulk *mmBulk, FILE *output){

	fprintf(output, "Memory allocated  : %llu\n", MMBulkByteAllocated(mmBulk));
	fprintf(output, "Memory dispatched : %llu\n", MMBulkByteDispatched(mmBulk));

}

void MMBulkSave(MMBulk *mmBulk, FILE *output) {

	unsigned long long i;

	fwrite(&mmBulk->itemSize, sizeof(unsigned int), 1, output);
	fwrite(&mmBulk->itemPerAllocationInPowerOf2, sizeof(unsigned int), 1, output);
	fwrite(&mmBulk->boundaryCushionSize, sizeof(unsigned int), 1, output);
	fwrite(&mmBulk->currentDirectoryEntry, sizeof(unsigned int), 1, output);
	fwrite(&mmBulk->nextUnusedItem, sizeof(unsigned int), 1, output);
	fwrite(&mmBulk->directorySize, sizeof(unsigned int), 1, output);

	for (i=0; i<mmBulk->currentDirectoryEntry; i++) {
		fwrite(mmBulk->directory[i], mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2, 1, output);
	}

	if (mmBulk->nextUnusedItem > 0) {
		fwrite(mmBulk->directory[i], mmBulk->itemSize * mmBulk->nextUnusedItem, 1, output);
	}

}

MMBulk *MMBulkLoad(MMPool *mmPool, FILE *input) {

	unsigned long long i;
	MMBulk *mmBulk;

	mmBulk = (MMBulk*) MMPoolDispatch(mmPool, sizeof(MMBulk));

	fread(&mmBulk->itemSize, sizeof(unsigned int), 1, input);
	fread(&mmBulk->itemPerAllocationInPowerOf2, sizeof(unsigned int), 1, input);
	fread(&mmBulk->boundaryCushionSize, sizeof(unsigned int), 1, input);
	fread(&mmBulk->currentDirectoryEntry, sizeof(unsigned int), 1, input);
	fread(&mmBulk->nextUnusedItem, sizeof(unsigned int), 1, input);
	fread(&mmBulk->directorySize, sizeof(unsigned int), 1, input);

	mmBulk->indexMask = truncateLeft(ALL_ONE_MASK,  BITS_IN_WORD - mmBulk->itemPerAllocationInPowerOf2);

	mmBulk->directory = (unsigned char**) MMPoolDispatch(mmPool, sizeof(unsigned char*) * mmBulk->directorySize);

	for (i=0; i<mmBulk->currentDirectoryEntry; i++) {
		mmBulk->directory[i] = (unsigned char *) MEMALIGN(mmBulk->boundaryCushionSize * 2 + 
									(mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2), MAX_ALIGN);
		if (mmBulk->directory[i] == NULL) {
			fprintf(stderr, "MMBulkLoad() : cannot allocate memory!\n");
			exit(1);
		}

		//Advance the address by boundaryCushionSize
		mmBulk->directory[i] += mmBulk->boundaryCushionSize;
		
		fread(mmBulk->directory[i], mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2, 1, input);
	}

	mmBulk->directory[i] = (unsigned char *) MEMALIGN(mmBulk->boundaryCushionSize * 2 + 
								(mmBulk->itemSize << mmBulk->itemPerAllocationInPowerOf2), MAX_ALIGN);
	if (mmBulk->directory[i] == NULL) {
		fprintf(stderr, "MMBulkLoad() : cannot allocate memory!\n");
		exit(1);
	}

	//Advance the address by boundaryCushionSize
	mmBulk->directory[i] += mmBulk->boundaryCushionSize;

	if (mmBulk->nextUnusedItem > 0) {
		fread(mmBulk->directory[i], mmBulk->itemSize * mmBulk->nextUnusedItem, 1, input);
	}


	for (i=0; i<mmMaster.maxNumberOfBulks; i++) {
		if (mmMaster.mmBulk[i] == NULL) {
			mmMaster.mmBulk[i] = mmBulk;
			return mmBulk;
		}
	}

	fprintf(stderr, "MMBulkLoad() : number of bulks > maxNumberOfBulk!\n");
	exit(1);

}
