#pragma once
#include "Cell.h"
#include "ConstDefine.h"
#include "MBB.h"
#include "QueryResult.h"
#include "Trajectory.h"
#include <iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include "cudaKernel.h"
#include <map>
#include <bitset>
#include"FVTable.h"

#define MAX_LEVEL (1048576-1)/3
extern Trajectory* tradb;

typedef struct MortonNode {		
	int level;
	int nid;
}MortonNode;

class MortonGrid
{
public:

	std::ofstream fout;
	MortonGrid(const MBB& mbb, float val_cell_size, int VITURAL_CELL_PARAM);
	MortonGrid();
	~MortonGrid();

	MBB range;
	float cell_size;
	int cellNum_axis;
	int cellnum;
	int totalPointNum;
	int trajNum;
	int VITURAL_CELL_PARAM;

	// ************* INDEX *************
	Cell* cellPtr;
	std::vector<cellBasedTraj> cellBasedTrajectory;
	std::bitset<size_t(MAX_LEVEL)> isLeaf;	


	SPoint* allPoints;
	Point*	allPointsPtrGPU;
	DPoint* allPointsDeltaEncoding;


	// **************** 1 ************ GPU
	void *baseAddrRange[2];	
	void *stateTableGPU[2];	


	// **************** 2 ************ CPU
	RangeQueryStateTable* stateTableRange[2];
	int stateTableLength[2];

	// **************** 3 ************ 
	// Memory Allocated Table (MAT)
	std::map<int, void*> nodeAddrTable[2];

	int nodeAddrTableLength[2];

	int testCnt;

	// ****************** Similarity Query  ******************
	FVTable freqVectors;

	// index related
	int addTrajectoryIntoCell(Trajectory &t);
	int WhichCellPointIn(SamplePoint p);
	int addDatasetToGrid(Trajectory* db, int traNum);
	int writeCellsToFile(int* cellNo, int cellNum, std::string file);
	int buildQuadTree(int level, int id);
	MBB generateMBBfromNode(int level, int id);

	//rangeQuery
	int rangeQueryBatch(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize);
	int rangeQueryBatchMultiThread(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize);
	int findMatchNodeInQuadTree(MortonNode node, MBB& bound, std::vector<MortonNode> *cells);

	int rangeQueryBatchGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize, RangeQueryStateTable* stateTableAllocate, int device_idx);
	int rangeQueryBatchGPUNoMAT(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize, RangeQueryStateTable* stateTableAllocate, int device_idx);
	
	int rangeQueryBatchMultiGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize);
	int findMatchNodeInQuadTreeGPU(MortonNode node, MBB& bound, std::vector<MortonNode> *cells, cudaStream_t stream, int queryID, int device_idx);
	int findMatchNodeInQuadTreeGPUNoMAT(MortonNode node, MBB& bound, std::vector<MortonNode> *cells, cudaStream_t stream, int queryID, int device_idx);
	
	//SimilarityQuery
	int SimilarityQueryBatch(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue);
	int SimilarityQueryBatchCPUParallel(Trajectory *qTra, int queryTrajNum, int *EDRdistance, int kValue);
	int SimilarityMultiThreadHandler(std::priority_queue<FDwithID, std::vector<FDwithID>, cmp>* queryQueue, Trajectory* qTra, int queryTrajNum, std::priority_queue<FDwithID, std::vector<FDwithID>, cmpBig>* EDRCalculated, int kValue, int startQueryIdx);
	int FDCalculateParallelHandeler(std::priority_queue<FDwithID, std::vector<FDwithID>, cmp> *queue, std::map<int, int>* freqVectorQ);
	int SimilarityExecuter(SPoint* queryTra, SPoint** candidateTra, int queryLength, int* candidateLength, int candSize, int *resultArray);
	int SimilarityQueryBatchOnGPU(Trajectory * qTra, int queryTrajNum, int * topKSimilarityTraj, int kValue);
	int SimilarityQueryBatchOnMultiGPU(Trajectory * qTra, int queryTrajNum, int * topKSimilarityTraj, int kValue);

};

