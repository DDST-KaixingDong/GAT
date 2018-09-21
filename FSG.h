#pragma once
#include "ConstDefine.h"
#include "Cell.h"
#include <map>

#include "Trajectory.h"
#include "SamplePoint.h"

class FSG
{
public:
	//Grid�������������귶Χ
	MBB range;
	float cell_size; //length of a cell
	int cell_num_x, cell_num_y; //�������ж��ٸ�cell
	int cellnum; //upper(area(grid)/area(cell))����֤�ܷ�������cell
	Cell* cellPtr; //�洢cell�����
	std::ofstream fout;//�ļ�����ӿ�
	int totalPointNum; //grid�ڵ����
	int trajNum;


	SPoint* allPoints;//�洢���е������
	Point* allPointsPtrGPU;
	DPoint *allPointsDeltaEncoding;//Delta Encoding��ĵ�

	//Range Query on GPU ��
	void *baseAddrRange[2];
	void *stateTableGPU[2];
	RangeQueryStateTable* stateTableRange[2];
	std::map<int, void*> nodeAddrTable[2];
	int stateTableLength[2];
	//int nodeAddrTableLength[2];

	FSG(const MBB& mbb, float val_cell_size);
	int addTrajectoryIntoCell(Trajectory &t);
	int WhichCellPointIn(SamplePoint p);
	int addDatasetToGrid(Trajectory* db, int traNum);

	//Query����
	int rangeQueryBatchGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize, RangeQueryStateTable* stateTableAllocate, int device_idx);
	int rangeQueryBatchMultiGPU(MBB *bounds, int rangeNum, CPURangeQueryResult *ResultTable, int *resultSetSize);

	FSG();
	~FSG();
};

