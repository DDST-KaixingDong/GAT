#pragma once
// MODE DEFINE

#ifdef WIN32

#else
#define USE_MULTIGPU
#endif

#ifdef WIN32

#else
#define FDMULCPU
#endif

// #define CHECK_CORRECT
// #define IfNoMAT 1

#define DATAINDEX 9 // 5 datasize = 7.5 GB
//#define IfSimQuery  0 // sim-query range-query

#define CELL_LEN 0.005 // here for fig.11 attention::: 0.01 n = 9 ;0.02 n = 8


#define MAXPOINTINNODE 20000 // here for fig.10 a little bit hard !! this is too inconvenient but no need to re-make?no so this is not recommended for tuning-parameter. 


#define ALPHA 0.2 // here for fig.13

#define KSIMILARITY 40

#define MAXDYNAMICM 200 

#ifdef WIN32
#define BIG_MEM 1024 // 1 GB
#define SMALL_MEM 512 // 256
#else
#define BIG_MEM 8000 // 
#define SMALL_MEM 2000 //
#endif


#define MAXLENGTH 1024

#define MAXTHREAD 256


#define MAX_TRAJ_SIZE 500000

#define MAXGAP 3600

//#define VITURAL_CELL_PARAM 8 

#define EPSILON 0.00025

#define N_BATCH_QUERY 2048
#define TRUE 1
#define FALSE 0


// #define NOT_COLUMN_ORIENTED

#include <stdio.h>
#include <string>
#include <math.h>
#include "QueryResult.h"
#include <cstring>
#include <thread>
#include <algorithm>
#include "MBB.h"

//决定计时函数
#ifdef WIN32
	#include "WinTimer.h"
#else
	#include <sys/time.h>
#endif

#ifdef WIN32
#else

class MyTimer
{
public:
	MyTimer() {
	};
	double iStart;
	double iEnd;

	double cpuSecond() {
		struct timeval tp;
		gettimeofday(&tp, NULL);
		return ((double)tp.tv_sec + (double)tp.tv_usec*1.e-6);
	}

	inline void start()
	{
		iStart = cpuSecond() * 1000;
	}
	inline void stop()
	{
		iEnd = cpuSecond() * 1000;
	}
	inline float elapse()
	{
		return iEnd - iStart;
	}
};
#endif



#define _CELL_BASED_STORAGE
//test:Similarity query based on naive grid，以定大小的grid来索引
//#define _SIMILARITY


//4+4+4+4=16 bytes
typedef struct Point {
	float x;
	float y;
	int time;
	int tID;
}Point;

//4+4+4=12 bytes
typedef struct SPoint {
	float x;
	float y;
	int tID; // 轨迹ID
}SPoint;

//2+2+4=8 bytes
typedef struct DPoint {
	short x;
	short y;
	int tID;
}DPoint;


typedef struct cellBasedTraj {
	int *cellNo = NULL;									
	int *startIdx = NULL;
	short *numOfPointInCell = NULL;	

	short length;
	int trajLength;

}cellBasedTraj;


// 4+4*7=32byte
typedef struct RangeQueryStateTable {//  leafnode

	void* ptr;
	int candidatePointNum;

	float xmin;
	float ymin;
	float xmax;
	float ymax;

	// cpu
	int queryID;
	int startIdxInAllPoints;

}RangeQueryStateTable;



typedef struct OffsetTable {
	int objectId;
	void *addr;
}OffsetTable;

typedef struct TaskInfoTableForSimilarity {
	int qID;
	int candTrajID;
}TaskInfoTableForSimilarity;

typedef struct intPair{
	int int_1;
	int int_2;
}intPair;



typedef unsigned char uint8_t;
typedef unsigned short uint16_t;
typedef unsigned int uint32_t;


int getIdxFromXY(int x, int y);
bool nonmem_cmp2(MBB t1, MBB t2); 

