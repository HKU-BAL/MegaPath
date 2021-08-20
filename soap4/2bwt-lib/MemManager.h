/*

   MemManager.h		Memory Manager

   This module provides memory management functions.

   Copyright (C) 2004, Wong Chi Kwong.

   This program is FREEALIGN software; you can redistribute it and/or
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

#ifndef __MEM_MANAGER_H__
#define __MEM_MANAGER_H__

#include "TypeNLimit.h"
#include <stdint.h>

#define MAX_ALIGN	64	// All memory except pool memory are aligned to MAX_ALIGN; pool memory is aligned to finer boundary for small memory size
#define MIN_ALIGN	1

#define RECORD_GRAND_TOTAL

//	Memory type:
//
//		unit memory:	allocation managed by malloc() individually;
//						to be used for large and less frequently accessed items
//						allocation can be freed individually at any time
//		pool memory:	pre-allocated memory pool for items with varying sizes
//						allocation cannot be freed by individually
//						to be used for small and frequently accessed items
//		temp memory:	temporary use granted from pool memory
//						allocation is allocated and freed like the items in a stack
//						pool memory allocation is disabled while temporary memory is in use
//		bulk memory:	pre-allocated memory pool for items with the same size
//						to be used for massively numbered items
//						memory address of dispatched items can be calculated by dispatch index


#ifdef DEBUG
#define Mem(mmBulk, index)  MMBulkAddress(mmBulk, index)
#else
#define Mem(mmBulk, index) 	(void*)&(mmBulk->directory[index >> mmBulk->itemPerAllocationInPowerOf2][(index & mmBulk->indexMask) * mmBulk->itemSize])
#endif

typedef struct MMPool {
	unsigned long long poolSize;						// Size of memory pool; the beginning of the pool holds the MMPool structure
	unsigned long long poolByteDispatched;			// Includes any spillover and memory skipped for align
	unsigned long long poolByteSpillover;				// Exclude spillover pointers
	unsigned long long currentTempByteDispatched;		// Includes any spillover
	unsigned long long currentTempByteSpillover;		// Exclude spillover pointers
	unsigned long long maxTotalByteDispatched;		// The max of pool memory + temp memory dispatched
	void *firstSpillOverAddress;				// if pool is freed, = address of mmPool
} MMPool;


typedef struct MMBulk {
	unsigned long long itemSize;
	unsigned long long itemPerAllocationInPowerOf2;
	unsigned long long boundaryCushionSize;			// boundary cushion is a piece of memory allocated so that the memory around items can be safely referenced
	unsigned long long indexMask;
	unsigned long long currentDirectoryEntry;
	unsigned long long nextUnusedItem;
	unsigned long long directorySize;
	unsigned char **directory;			// if bulk is freed, = NULL
} MMBulk;

typedef struct MMMaster {
	long long currentUnitByteAllocated;
	long long maxUnitByteAllocated;
	unsigned long long maxNumberOfPools;
	MMPool **mmPool;
	unsigned long long maxNumberOfBulks;
	MMBulk **mmBulk;
	long long maxTotalByteAllocated;
	long long maxTotalByteDispatched;
	int traceUnitByteAllocation;
	FILE *unitByteTraceFile;
} MMMaster;

void *MMMalloc(const long long memSize);
void MMFree(void *address);
void MMMasterInitialize(const unsigned long long maxNumberOfPools, const unsigned long long maxNumberOfBulks,
						const long long traceUnitByteAllocation, FILE *unitByteTraceFile);
void MMMasterFreeAll();
long long MMMasterCurrentTotalByteAllocated();
long long MMMasterCurrentTotalByteDispatched();
long long MMMasterMaxTotalByteAllocated();
long long MMMasterMaxTotalByteDispatched();
void MMMasterSetMaxTotalByteAllocated();
void MMMasterSetMaxTotalByteDispatched();
void MMMasterPrintReport(FILE *output, const unsigned long long withUnitDetails, const unsigned long long withPoolDetails, const unsigned long long withBulkDetails);

void *MMUnitAllocate(const long long memSize);
void *MMUnitReallocate(void *address, const long long newMemSize, const long long oldMemSize);
void MMUnitFree(void *address, const long long memSize);
long long MMUnitCurrentByteAllocated();
long long MMUnitMaxByteAllocated();
void MMUnitPrintReport(FILE *output);

MMPool *MMPoolCreate(const unsigned long long poolSize);
unsigned int MMPoolIsActive(const MMPool *mmPool);
void MMPoolSetInactive(MMPool *mmPool);
unsigned long long MMPoolCurrentTotalByteAllocated(const MMPool *mmPool);
unsigned long long MMPoolCurrentTotalByteDispatched(const MMPool *mmPool);
unsigned long long MMPoolMaxTotalByteDispatched(const MMPool *mmPool);
unsigned long long MMPoolByteAvailable(const MMPool *mmPool);
MMPool *MMPoolFree(MMPool *mmPool);
void MMPoolReset(MMPool *mmPool);
void MMPoolDestory(MMPool *mmPool);
void *MMPoolDispatch(MMPool *mmPool, const unsigned long long memSize);
unsigned long long MMPoolDispatchOffset(MMPool *mmPool, const unsigned long long memSize);
void MMPoolReturn(MMPool *mmPool, void *address, const unsigned long long memSize);		// Dummy function
void MMPoolPrintReport(MMPool *mmPool, FILE *output);

void *MMTempDispatch(MMPool *mmPool, const unsigned long long memsize);
void MMTempReturn(MMPool *mmPool, void *address, const unsigned long long memSize);
void MMTempPrintReport(MMPool *mmPool, FILE *output);

MMBulk *MMBulkCreate(MMPool *mmPool, const unsigned long long itemSize, const unsigned long long itemPerAllocationInPowerOf2, 
					 unsigned long long const boundaryCushionSize, unsigned long long const directorySize);
unsigned int MMBulkIsActive(const MMBulk *mmBulk);
void MMBulkSetInactive(MMBulk *mmBulk);
unsigned long long MMBulkByteAllocated(const MMBulk *mmBulk);
unsigned long long MMBulkByteDispatched(const MMBulk *mmBulk);
unsigned long long MMBulkUnitDispatched(const MMBulk *mmBulk);
void MMBulkFree(MMBulk *mmBulk);
void MMBulkDestory(MMBulk *mmBulk);
unsigned long long MMBulkDispatch(MMBulk *mmBulk);
void *MMBulkAddress(const MMBulk *mmBulk, const unsigned long long index);
MMPool *MMBulkFindPoolUsed(const MMBulk *mmBulk);
void MMBulkPrintReport(MMBulk *mmBulk, FILE *output);

void MMBulkSave(MMBulk *mmBulk, FILE *output);
MMBulk *MMBulkLoad(MMPool *mmPool, FILE *input);


#endif
