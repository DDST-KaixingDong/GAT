#pragma once
#include "ConstDefine.h"
#include<map>
#include<vector>
#include <queue>


typedef struct FDwithID {
	int traID;
	int FD;
}FDwithID;

struct cmp {
	bool operator()(FDwithID a, FDwithID b) {
		return(a.FD > b.FD);
	}
};

struct cmpBig {
	bool operator()(FDwithID a, FDwithID b) {
		return(a.FD < b.FD);
	}
};

class FVTable
{
public:
	std::vector<std::map<int, int>> FreqVector;

	int trajNum;
	int cellNum;
	void *FVTableGPU, *FVinfoGPU,*queryFVGPU,*FVTableOffset,*FDresultsGPU,*SubbedArrayGPU,*SubbedArrayOffsetGPU;
	int SubbedArrayJump = 0;
	size_t pitch;
	int nonZeroFVNum = 0;

	int initFVTable(int trajNum, int cellNum);
	int addPointToFVTable(int trajID, int pointNum, int cellID);
	int getCandidate(int lowerBound, int k, std::map<int, int>* freqVectorQ, int *candidateTrajID, int *candidateNum);
	double calculateFreqDist(int *freqVectorQ, int trajID);
	int findNeighbor(int cellID, int* neighborID);
	int formPriorityQueue(std::priority_queue<FDwithID, std::vector<FDwithID>, cmp> *queue, std::map<int, int>* freqVectorQ);
	int transferFVtoGPU();
	int formPriorityQueueGPU(std::priority_queue<FDwithID, std::vector<FDwithID>, cmp> *queue, std::map<int, int>* freqVectorQ);

	FVTable();
	~FVTable();
};

