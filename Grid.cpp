#include "Grid.h"

// for test

extern Trajectory* tradb;
extern void* baseAddrGPU;
using namespace std;

#ifdef WIN32
MyTimer timer;
#else

#endif


int getIdxFromXY(int x, int y)
{
	if (x < 0 || y < 0) {
		return -1; //太小超界
	}
	int lenx, leny;
	if (x == 0)
		lenx = 1;
	else
	{
		lenx = int(log2(x)) + 1;
	}
	if (y == 0)
		leny = 1;
	else
		leny = int(log2(y)) + 1;

	int result = 0;
	int xbit = 1, ybit = 1;
	for (int i = 1; i <= 2 * max(lenx, leny); i++) // 标识向左移动多少位
	{
		if ((i & 1) == 1) // 二进制末尾为1 该取列号
		{
			result += (x >> (xbit - 1) & 1) * (1 << (i - 1)); // x >> (xbit - 1) & 1 巧妙 运用位运算 移位运算
			xbit = xbit + 1;
		}
		else // 二进制末尾为0  该取行号
		{
			result += (y >> (ybit - 1) & 1) * (1 << (i - 1));
			ybit = ybit + 1;
		}
	}
	return result;
}

int Grid::GetIDfromXYTest(int x, int y) {
	return getIdxFromXY(x, y);
}

//确认无误
int Grid::WhichCellPointIn(SamplePoint p)
{
	//注意cell编号是从(xmin,ymax)开始的，而不是(xmin,ymin)
	int row = (int)((range.ymax - p.lat) / cell_size); //从0开始
	int col = (int)((p.lon - range.xmin) / cell_size); //从0开始
	return getIdxFromXY(col, row);
}

Grid::Grid()
{
	range = MBB(0, 0, 0, 0);
	cellnum = 0;
	cell_size = 0;
	cellNum_axis = 0;
	cellPtr = NULL;
	allPoints = NULL;
	allPointsPtrGPU = NULL;
}

Grid::Grid(const MBB& mbb, float val_cell_size, int VITURAL_CELL_PARAM)
{
	this->nodeNum = 0;
	this->VITURAL_CELL_PARAM = VITURAL_CELL_PARAM;
	range = mbb;
	cell_size = val_cell_size;
	int divideNumOnX = (int)((mbb.xmax - mbb.xmin) / val_cell_size) + 1; 
	int divideNumOnY = (int)((mbb.ymax - mbb.ymin) / val_cell_size) + 1;
	int maxValue = max(divideNumOnX, divideNumOnY);
	cellNum_axis = maxValue >> (int(log2(maxValue))) << (int(log2(maxValue)) + 1);

	cout << cellNum_axis << endl;

	cellnum = cellNum_axis * cellNum_axis;
	cellPtr = new Cell[cellnum];
	this->testCnt = 0;
	for (int i = 0; i <= cellNum_axis - 1; i++)
	{
		for (int j = 0; j <= cellNum_axis - 1; j++)
		{
			int cell_idx = getIdxFromXY(j, i);
			cellPtr[cell_idx].initial(i, j, MBB(range.xmin + cell_size * j, range.ymax - cell_size * (i + 1), range.xmin + cell_size * (j + 1), range.ymax - cell_size * (i)));
		}
	}
}

Grid::~Grid()
{
}






int Grid::buildQuadTree(int level, int id, QuadtreeNode* pNode, QuadtreeNode* parent)
{
	int totalLevel = int(log2(this->cellnum) / log2(4));
	int totalPoints = 0;
	for (int i = id * int(pow(4, (totalLevel - level))); i <= (id + 1) * int(pow(4, (totalLevel - level))) - 1; i++)
	{
		totalPoints += this->cellPtr[i].totalPointNum;
	}
	pNode->mbb = MBB(this->cellPtr[id * int(pow(4, (totalLevel - level)))].mbb.xmin, this->cellPtr[(id + 1) * int(pow(4, (totalLevel - level))) - 1].mbb.ymin, this->cellPtr[(id + 1) * int(pow(4, (totalLevel - level))) - 1].mbb.xmax, this->cellPtr[id * int(pow(4, (totalLevel - level)))].mbb.ymax);
	pNode->numPoints = totalPoints;
	pNode->NodeID = id;
	pNode->parent = parent;
	pNode->level = level;
	if ((totalPoints < MAXPOINTINNODE) || (level == totalLevel))
	{
		pNode->isLeaf = true;
		pNode->DL = NULL;
		pNode->DR = NULL;
		pNode->UL = NULL;
		pNode->UR = NULL;
		return 0;
	}
	else
	{
		this->nodeNum += 4;
		pNode->isLeaf = false;
		pNode->UL = (QuadtreeNode*)malloc(sizeof(QuadtreeNode));
		this->buildQuadTree(level + 1, id << 2, pNode->UL, pNode);
		pNode->UR = (QuadtreeNode*)malloc(sizeof(QuadtreeNode));
		this->buildQuadTree(level + 1, (id << 2) + 1, pNode->UR, pNode);
		pNode->DL = (QuadtreeNode*)malloc(sizeof(QuadtreeNode));
		this->buildQuadTree(level + 1, (id << 2) + 2, pNode->DL, pNode);
		pNode->DR = (QuadtreeNode*)malloc(sizeof(QuadtreeNode));
		this->buildQuadTree(level + 1, (id << 2) + 3, pNode->DR, pNode);
		return 0;
	}
}

//把轨迹t打碎成子轨迹，添加到cell里面
//这一步仅仅是把子轨迹放进了cell里面，组成一个个item
int Grid::addTrajectoryIntoCell(Trajectory& t)
{
	if (t.length == 0)
		return 1;//空轨迹
	SamplePoint p = t.points[0];
	int lastCellNo = WhichCellPointIn(p);
	int lastCellStartIdx = 0;
	int nowCellNo;
	//cell based traj生成，记得转换后free！
	vector<int>* tempCellBasedTraj = new vector<int>;
	tempCellBasedTraj->reserve(1048577);
	int tempCellNum = 0;
	for (int i = 0; i <= t.length - 1; i++)
	{
		p = t.points[i];
		nowCellNo = WhichCellPointIn(p);

		if (i == t.length - 1)
		{
			//到最后一条，发现这个cell也是上个cell就是最后一个cell了，添加之
			if (lastCellNo == nowCellNo)
			{
				tempCellNum++;
				tempCellBasedTraj->push_back(nowCellNo);
				cellPtr[nowCellNo].addSubTra(t.tid, lastCellStartIdx, i, i - lastCellStartIdx + 1);
				int vituralCellNo = nowCellNo >> VITURAL_CELL_PARAM; //虚格子
				this->freqVectors.addPointToFVTable(t.tid, i - lastCellStartIdx + 1, vituralCellNo);
			}
			//否则，上一个和这个cell都要添加
			else
			{
				tempCellNum += 2;
				tempCellBasedTraj->push_back(lastCellNo);
				tempCellBasedTraj->push_back(nowCellNo);
				cellPtr[lastCellNo].addSubTra(t.tid, lastCellStartIdx, i - 1, i - 1 - lastCellStartIdx + 1);
				int vituralCellNo = lastCellNo >> VITURAL_CELL_PARAM; //虚格子
				this->freqVectors.addPointToFVTable(t.tid, i - 1 - lastCellStartIdx + 1, vituralCellNo);
				cellPtr[nowCellNo].addSubTra(t.tid, i, i, 1);
				vituralCellNo = nowCellNo >> VITURAL_CELL_PARAM;
				this->freqVectors.addPointToFVTable(t.tid, 1, vituralCellNo);
			}
		}
		else
		{
			if (lastCellNo == nowCellNo)
				continue;
			else
			{
				// 终结一条子轨迹，开始下一条子轨迹
				//cellTra里面加一条
				tempCellNum++;
				tempCellBasedTraj->push_back(lastCellNo);
				//SubTra添加
				//printf("cell:%d\n", lastCellNo);
				cellPtr[lastCellNo].addSubTra(t.tid, lastCellStartIdx, i - 1, i - 1 - lastCellStartIdx + 1);
				int vituralCellNo = lastCellNo >> VITURAL_CELL_PARAM; //虚格子
				this->freqVectors.addPointToFVTable(t.tid, i - 1 - lastCellStartIdx + 1, vituralCellNo);
				lastCellNo = nowCellNo;
				lastCellStartIdx = i;
			}
		}
	}
	this->cellBasedTrajectory[t.tid].length = tempCellNum;
	this->cellBasedTrajectory[t.tid].cellNo = (int*)malloc(sizeof(int) * tempCellNum);
	if (this->cellBasedTrajectory[t.tid].cellNo == NULL) throw("alloc error");
	for (int i = 0; i <= tempCellNum - 1; i++)
	{
		this->cellBasedTrajectory[t.tid].cellNo[i] = tempCellBasedTraj->at(i);
	}
	this->cellBasedTrajectory[t.tid].trajLength = t.length;
	delete tempCellBasedTraj;
	return 0;
}

int Grid::addDatasetToGrid(Trajectory* db, int traNum)
{
	this->trajNum = traNum;
	//生成frequency vector
	this->freqVectors.initFVTable(traNum, (this->cellnum)>> this->VITURAL_CELL_PARAM); //注意这里是虚格子的个数
	//注意，轨迹编号从1开始
	this->cellBasedTrajectory.resize(traNum + 1); //扩大cellbasedtraj的规模，加轨迹的时候可以直接用
	int pointCount = 0;
	for (int i = 1; i <= traNum; i++)
	{
		addTrajectoryIntoCell(db[i]);
	}
	for (int i = 0; i <= cellnum - 1; i++)
	{
		cellPtr[i].buildSubTraTable();
		pointCount += cellPtr[i].totalPointNum;
	}
	this->totalPointNum = pointCount;
	this->root = (QuadtreeNode*)malloc(sizeof(QuadtreeNode));
	this->buildQuadTree(0, 0, this->root, NULL);
	this->allPoints = (SPoint*)malloc(sizeof(SPoint) * (this->totalPointNum));
	pointCount = 0;


	for (int i = 0; i <= cellnum - 1; i++)
	{
		cellPtr[i].pointRangeStart = pointCount;
		for (int j = 0; j <= cellPtr[i].subTraNum - 1; j++)
		{
			//for each subTra, add Points to AllPoints
			cellPtr[i].subTraTable[j].idxInAllPointsArray = pointCount;
			for (int k = cellPtr[i].subTraTable[j].startpID; k <= cellPtr[i].subTraTable[j].endpID; k++)
			{
				allPoints[pointCount].tID = cellPtr[i].subTraTable[j].traID;
				allPoints[pointCount].x = tradb[allPoints[pointCount].tID].points[k].lon;
				allPoints[pointCount].y = tradb[allPoints[pointCount].tID].points[k].lat;
				//allPoints[pointCount].time = tradb[allPoints[pointCount].tID].points[k].time;
				pointCount++;
			}
		}
		cellPtr[i].pointRangeEnd = pointCount - 1;
		if (cellPtr[i].pointRangeEnd - cellPtr[i].pointRangeStart + 1 != cellPtr[i].totalPointNum)
			cerr << "Grid.cpp: something wrong in total point statistic" << endl;
	}
	for (int i = 1; i <= this->trajNum; i++)
	{
		this->cellBasedTrajectory[i].startIdx = (int*)malloc(sizeof(int) * this->cellBasedTrajectory[i].length);
		this->cellBasedTrajectory[i].numOfPointInCell = (short*)malloc(sizeof(short) * this->cellBasedTrajectory[i].length);
		int* tempCntForTraj = (int*)malloc(sizeof(int) * this->cellnum);
		memset(tempCntForTraj, 0, sizeof(int) * this->cellnum);
		for (int cellidx = 0; cellidx <= this->cellBasedTrajectory[i].length - 1; cellidx++)
		{
			int nowCellID = this->cellBasedTrajectory[i].cellNo[cellidx];
			int j, cnt;
			for (j = 0 , cnt = 0; cnt <= tempCntForTraj[nowCellID]; j++)
			{
				if (this->cellPtr[nowCellID].subTraTable[j].traID == i)
				{
					cnt++;
				}
			}
			j--;
			//choose j; //选择第j项subTra，但是这里面还是没有起始点的信息，考虑回来把point的信息加进subtraTable。。。反正最后要把subTraTable给free掉
			this->cellBasedTrajectory[i].startIdx[cellidx] = this->cellPtr[nowCellID].subTraTable[j].idxInAllPointsArray;
			this->cellBasedTrajectory[i].numOfPointInCell[cellidx] = this->cellPtr[nowCellID].subTraTable[j].numOfPoint;
			tempCntForTraj[nowCellID]++;
		}
		free(tempCntForTraj);
	}
	// Transfer FV to GPU
	//this->freqVectors.transferFVtoGPU();


	////Delta Encoding的cell连续存储
	//this->allPointsDeltaEncoding = (DPoint*)malloc(sizeof(DPoint)*(this->totalPointNum));
	//pointCount = 0;
	//for (int i = 0; i <= cellnum - 1; i++) {
	//	cellPtr[i].pointRangeStart = pointCount;
	//	for (int j = 0; j <= cellPtr[i].subTraNum - 1; j++) {
	//		for (int k = cellPtr[i].subTraTable[j].startpID; k <= cellPtr[i].subTraTable[j].endpID; k++) {
	//			allPointsDeltaEncoding[pointCount].tID = cellPtr[i].subTraTable[j].traID;
	//			allPointsDeltaEncoding[pointCount].x = short(int((tradb[allPointsDeltaEncoding[pointCount].tID].points[k].lon)*1000000)-cellPtr[i].anchorPointX);
	//			allPointsDeltaEncoding[pointCount].y = short(int((tradb[allPointsDeltaEncoding[pointCount].tID].points[k].lat)*1000000)-cellPtr[i].anchorPointY);
	//			pointCount++;
	//		}
	//	}
	//	cellPtr[i].pointRangeEnd = pointCount - 1;
	//	if (cellPtr[i].pointRangeEnd - cellPtr[i].pointRangeStart + 1 != cellPtr[i].totalPointNum)
	//		cerr << "Grid.cpp: something wrong in total point statistic" << endl;
	//}

	//把生成好的allpoints放到GPU内
	//putCellDataSetIntoGPU(this->allPoints, this->allPointsPtrGPU, this->totalPointNum);


	return 0;
}

int Grid::writeCellsToFile(int* cellNo, int cellNum, string file)
// under editing....
{
	fout.open(file, ios_base::out);
	for (int i = 0; i <= cellNum - 1; i++)
	{
		int outCellIdx = cellNo[i];
		cout << outCellIdx << ": " << "[" << cellPtr[outCellIdx].mbb.xmin << "," << cellPtr[outCellIdx].mbb.xmax << "," << cellPtr[outCellIdx].mbb.ymin << "," << cellPtr[outCellIdx].mbb.ymax << "]" << endl;
		for (int j = 0; j <= cellPtr[outCellIdx].subTraNum - 1; j++)
		{
			int tid = cellPtr[outCellIdx].subTraTable[j].traID;
			int startpid = cellPtr[outCellIdx].subTraTable[j].startpID;
			int endpid = cellPtr[outCellIdx].subTraTable[j].endpID;
			for (int k = startpid; k <= endpid; k++)
			{
				cout << tradb[tid].points[k].lat << "," << tradb[tid].points[k].lon << ";";
			}
			cout << endl;
		}
	}
	return 0;
}






int Grid::rangeQueryBatch(MBB* bounds, int rangeNum, CPURangeQueryResult* ResultTable, int* resultSetSize)
{
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		ResultTable[i].resize(this->trajNum + 1);
	}
#ifdef CHECK_CORRECT

	for (int i = 0; i <= rangeNum - 1;i++)
	{
		for (int j = 0; j <= this->trajNum + 1;j++)
		{
			ResultTable[i][j] = 0;
		}
	}
#endif

	int totalLevel = int(log2(this->cellnum) / log2(4));
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		//int candidateNodeNum = 0;
		vector<QuadtreeNode*> cells;
		findMatchNodeInQuadTree(this->root, bounds[i], &cells);
		//printf("%d", cells.size());
		for (vector<QuadtreeNode*>::iterator iterV = cells.begin(); iterV != cells.end(); iterV++)
		{
			int nodeID = (*iterV)->NodeID;
			int nodeLevel = (*iterV)->level;
			int firstCellID = nodeID * int(pow(4, (totalLevel - nodeLevel)));
			int lastCellID = (nodeID + 1) * int(pow(4, (totalLevel - nodeLevel))) - 1;
			for (int cellID = firstCellID; cellID <= lastCellID; cellID++)
			{
				int anchorX = this->cellPtr[cellID].anchorPointX;
				int anchorY = this->cellPtr[cellID].anchorPointY;
				for (int idx = this->cellPtr[cellID].pointRangeStart; idx <= this->cellPtr[cellID].pointRangeEnd; idx++)
				{
					//compress
					//float realX = float(allPointsDeltaEncoding[idx].x + anchorX) / 1000000;
					//float realY = float(allPointsDeltaEncoding[idx].y + anchorY) / 1000000;
					// no compress
					float realX = allPoints[idx].x;
					float realY = allPoints[idx].y;
					int tID = allPoints[idx].tID;
					if (bounds[i].pInBox(realX, realY))
					{
						ResultTable[i][tID] = TRUE;
					}
				}
			}
		}
	}

	//for (int jobID = 0; jobID <= rangeNum - 1; jobID++)
	//{
	//	for (int traID = 1; traID <= (this->trajNum); traID++)
	//	{
	//		if (resultsReturned[jobID * (this->trajNum + 1) + traID] == 1)
	//		{
	//			out << "job " << jobID << "find" << traID << endl;
	//		}
	//	}
	//}

	return 0;
}

int Grid::rangeQueryBatchMultiThread(MBB* bounds, int rangeNum, CPURangeQueryResult* ResultTable, int* resultSetSize)
{
	/*
	// way2:
	vector<thread> threads_RQ;
	const int device_num = 16; // 这里用了16！！！ 说明MCPU 知道的 那20估计是弄错了 含着CPU 会满载么
	int nowidx = 0;
	// 整个batch平分为16份
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		int boundNum = rangeNum / device_num + 1;
		if (nowidx + boundNum >= rangeNum)
			boundNum = rangeNum - nowidx;
		threads_RQ.push_back(thread(std::mem_fn(&Grid::rangeQueryBatch), this, &bounds[nowidx], boundNum, &ResultTable[nowidx], resultSetSize));
		nowidx += boundNum;
	}
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));
	
	*/

	// way1:
	vector<thread> threads_RQ;
	for (int qID = 0; qID <= rangeNum - 1; qID++)
	{
		threads_RQ.push_back(thread(std::mem_fn(&Grid::rangeQueryBatch), this, &bounds[qID], 1, &ResultTable[qID], resultSetSize));// resultSetSize 没用到
		//threads_RQ.push_back(thread(std::mem_fn(&Grid::SimilarityMultiThreadHandler), this, queryQueue, qTra, 1, EDRCalculated, kValue, qID));
		//threads_EDR.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}

	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join)); // threads_FD已经指明

	return 0;
}

int Grid::findMatchNodeInQuadTree(QuadtreeNode* node, MBB& bound, vector<QuadtreeNode*>* cells)
{
	if (node->isLeaf)
	{
		cells->push_back(node);
	}
	else
	{
		if (bound.intersect(node->UL->mbb))
			findMatchNodeInQuadTree(node->UL, bound, cells);
		if (bound.intersect(node->UR->mbb))
			findMatchNodeInQuadTree(node->UR, bound, cells);
		if (bound.intersect(node->DL->mbb))
			findMatchNodeInQuadTree(node->DL, bound, cells);
		if (bound.intersect(node->DR->mbb))
			findMatchNodeInQuadTree(node->DR, bound, cells);
	}
	return 0;
}

int Grid::rangeQueryBatchGPU(MBB* bounds, int rangeNum, CPURangeQueryResult* ResultTable, int* resultSetSize, RangeQueryStateTable* stateTableAllocate, int device_idx)
{

#ifdef CHECK_CORRECT
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		ResultTable[i].resize(this->trajNum + 1);
	}
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		for (int j = 0; j <= this->trajNum; j++)
		{
			ResultTable[i][j] = 0;
		}
	}
#endif

	MyTimer timer;
	// timer.start();
	CUDA_CALL(cudaSetDevice(device_idx));

	this->stateTableRange[device_idx] = stateTableAllocate;
	this->stateTableLength[device_idx] = 0;

	this->nodeAddrTableLength[device_idx] = 0;

	// for each query, generate the nodes:
	cudaStream_t stream;
	cudaStreamCreate(&stream);

	// *************************** step1: filtering ***************************
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		findMatchNodeInQuadTreeGPU(root, bounds[i], NULL, stream, i, device_idx);
	}

	int maxPointNum = 0;
	for (int i = 0; i <= this->stateTableLength[device_idx] - 1; i++)
	{
		if (stateTableAllocate[i].candidatePointNum > maxPointNum)
			maxPointNum = stateTableAllocate[i].candidatePointNum;
	}
	// timer.stop();
	// cout << "Time 1:" << timer.elapse() << "ms" << endl;

	// *************************** step2: verification ***************************
	// timer.start();
	CUDA_CALL(cudaMemcpyAsync(this->stateTableGPU[device_idx], stateTableAllocate, sizeof(RangeQueryStateTable)*this->stateTableLength[device_idx],
		cudaMemcpyHostToDevice, stream));
	// timer.stop();
	// cout << "Time 2:" << timer.elapse() << "ms" << endl;

	//传递完成，开始调用kernel查询
	uint8_t* resultsReturned = (uint8_t*)malloc(sizeof(uint8_t) * (this->trajNum + 1) * rangeNum);
	// timer.start();
	cudaRangeQueryTestHandler((RangeQueryStateTable*)this->stateTableGPU[device_idx], this->stateTableLength[device_idx], resultsReturned, this->trajNum + 1, rangeNum, stream);
	//ofstream fp("queryResult(GTS).txt", ios_base::out);
#ifdef CHECK_CORRECT

	for (int jobID = 0; jobID <= rangeNum - 1; jobID++)
	{
		for (int traID = 0; traID <= this->trajNum; traID++)
		{
			if (resultsReturned[jobID * (this->trajNum + 1) + traID] == 1)
			{
				ResultTable[jobID][traID] = TRUE;
			}
		}
	}
#endif
	//for (vector<uint8_t>::iterator iter = resultsReturned.begin(); iter != resultsReturned.end(); iter++) {
	//	//cout << (*iter) << endl;
	//	//printf("%d\n", *iter);
	//}
	// timer.stop();
	// cout << "Time 3:" << timer.elapse() << "ms" << endl;

	//FILE *fp = fopen("resultQuery.txt", "w+");
	//for (int i = 0; i <= stateTableLength - 1; i++) {
	//	for (int j = 0; j <= stateTableAllocate[i].candidatePointNum - 1; j++) {

	//		if ((resultsReturned[i*maxPointNum + j]) == (uint8_t)(1)) {
	//			fprintf(fp,"%d\n", stateTableAllocate[i].startIdxInAllPoints + j);
	//			fprintf(fp,"%f,%f\n", allPoints[stateTableAllocate[i].startIdxInAllPoints + j].x, allPoints[stateTableAllocate[i].startIdxInAllPoints + j].y);
	//		}

	//	}
	//}
	//查询结束，善后，清空stateTable，清空gpu等
	cudaStreamDestroy(stream);
	return 0;
}

int Grid::rangeQueryBatchGPUWithoutMorton(MBB * bounds, int rangeNum, CPURangeQueryResult * ResultTable, int * resultSetSize, RangeQueryStateTable * stateTableAllocate, int device_idx)
{
#ifdef CHECK_CORRECT
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		ResultTable[i].resize(this->trajNum + 1);
	}
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		for (int j = 0; j <= this->trajNum; j++)
		{
			ResultTable[i][j] = 0;
		}
	}
#endif
	MyTimer timer;
	// timer.start();
	CUDA_CALL(cudaSetDevice(device_idx));
	this->stateTableRange[device_idx] = stateTableAllocate;
	this->stateTableLength[device_idx] = 0;
	this->nodeAddrTableLength[device_idx] = 0;
	vector<int> blockOffsetInData, blockOffsetNum, blockOffsetOfOffset, blockLength;
	// for each query, generate the nodes:
	cudaStream_t stream;
	cudaStreamCreate(&stream);
	for (int i = 0; i <= rangeNum - 1; i++)
	{
		findMatchNodeInQuadTreeGPUWithoutMorton(root, bounds[i], NULL, stream, i, device_idx, 
			blockOffsetInData, blockOffsetNum, blockOffsetOfOffset, blockLength);
	}
	//printf("StateTableLength:%d",this->stateTableLength);
	//stateTable中点的数目的最大值
	// FILE* testFile = fopen("freq.txt", "w+");
	int maxPointNum = 0;
	for (int i = 0; i <= this->stateTableLength[device_idx] - 1; i++)
	{
		if (stateTableAllocate[i].candidatePointNum > maxPointNum)
			maxPointNum = stateTableAllocate[i].candidatePointNum;
		//fprintf(testFile, "%d ", stateTableAllocate[i].candidatePointNum);
	}
	//交给GPU进行并行查询
	//先传递stateTable
	// timer.stop();
	// cout << "Time 1:" << timer.elapse() << "ms" << endl;

	// timer.start();
	// transfer offset of uncontinuous part
	size_t allSpaceForOffsetContinuous = blockOffsetInData.size() + blockOffsetNum.size() +
		blockOffsetOfOffset.size() + blockLength.size();
	// printf("%zd", allSpaceForOffsetContinuous);
	int* startGPUTemp = NULL;
	CUDA_CALL(cudaMalloc((void**)&startGPUTemp, sizeof(int)*allSpaceForOffsetContinuous));
	void* GPUTempPtr = (void*)startGPUTemp;
	int* blockOffsetInDataGPU=startGPUTemp, *blockOffsetNumGPU, *blockOffsetOfOffsetGPU, *blockLengthGPU;
	CUDA_CALL(cudaMemcpyAsync(startGPUTemp, &blockOffsetInData[0], sizeof(int)*blockOffsetInData.size(),cudaMemcpyHostToDevice,stream));
	startGPUTemp += blockOffsetInData.size();
	blockOffsetNumGPU = startGPUTemp;
	CUDA_CALL(cudaMemcpyAsync(startGPUTemp, &blockOffsetNum[0], sizeof(int)*blockOffsetNum.size(), cudaMemcpyHostToDevice, stream));
	startGPUTemp += blockOffsetNum.size();
	blockOffsetOfOffsetGPU = startGPUTemp;
	CUDA_CALL(cudaMemcpyAsync(startGPUTemp, &blockOffsetOfOffset[0], sizeof(int)*blockOffsetOfOffset.size(), cudaMemcpyHostToDevice, stream));
	startGPUTemp += blockOffsetOfOffset.size();
	blockLengthGPU = startGPUTemp;
	CUDA_CALL(cudaMemcpyAsync(startGPUTemp, &blockLength[0], sizeof(int)*blockLength.size(), cudaMemcpyHostToDevice, stream));


	CUDA_CALL(cudaMemcpyAsync(this->stateTableGPU[device_idx], stateTableAllocate, sizeof(RangeQueryStateTable)*this->stateTableLength[device_idx],
		cudaMemcpyHostToDevice, stream));
	//传递完成，开始调用kernel查询
	uint8_t* resultsReturned = (uint8_t*)malloc(sizeof(uint8_t) * (this->trajNum + 1) * rangeNum);

	// timer.stop();
	// cout << "Time 2:" << timer.elapse() << "ms" << endl;

	// timer.start();
	//cudaRangeQueryTestHandler((RangeQueryStateTable*)this->stateTableGPU[device_idx], 
	//	this->stateTableLength[device_idx], resultsReturned, this->trajNum + 1, rangeNum, stream);
	cudaRangeQueryTestHandlerNonMorton((RangeQueryStateTable*)this->stateTableGPU[device_idx],
		this->stateTableLength[device_idx], resultsReturned, this->trajNum + 1, rangeNum, stream,
		blockOffsetInDataGPU, blockLengthGPU, blockOffsetNumGPU, blockOffsetOfOffsetGPU);
	//ofstream fp("queryResult(GTS).txt", ios_base::out);
#ifdef CHECK_CORRECT

	for (int jobID = 0; jobID <= rangeNum - 1; jobID++)
	{
		for (int traID = 0; traID <= this->trajNum; traID++)
		{
			if (resultsReturned[jobID * (this->trajNum + 1) + traID] == 1)
			{
				ResultTable[jobID][traID] = TRUE;
			}
		}
	}
#endif
	//for (vector<uint8_t>::iterator iter = resultsReturned.begin(); iter != resultsReturned.end(); iter++) {
	//	//cout << (*iter) << endl;
	//	//printf("%d\n", *iter);
	//}
	// timer.stop();
	// cout << "Time 3:" << timer.elapse() << "ms" << endl;

	//FILE *fp = fopen("resultQuery.txt", "w+");
	//for (int i = 0; i <= stateTableLength - 1; i++) {
	//	for (int j = 0; j <= stateTableAllocate[i].candidatePointNum - 1; j++) {

	//		if ((resultsReturned[i*maxPointNum + j]) == (uint8_t)(1)) {
	//			fprintf(fp,"%d\n", stateTableAllocate[i].startIdxInAllPoints + j);
	//			fprintf(fp,"%f,%f\n", allPoints[stateTableAllocate[i].startIdxInAllPoints + j].x, allPoints[stateTableAllocate[i].startIdxInAllPoints + j].y);
	//		}

	//	}
	//}
	//查询结束，善后，清空stateTable，清空gpu等
	cudaStreamDestroy(stream);
	CUDA_CALL(cudaFree(GPUTempPtr));
	return 0;
}

int Grid::rangeQueryBatchMultiGPU(MBB* bounds, int rangeNum, CPURangeQueryResult* ResultTable, int* resultSetSize)
{
	MyTimer timer;
	int device_num = 2;
	vector<thread> threads_RQ;
	int rangeNumGPU[2];
	rangeNumGPU[0] = rangeNum / 2;
	rangeNumGPU[1] = rangeNum - rangeNumGPU[0];
	int startIdx[2];
	startIdx[0] = 0;
	startIdx[1] = rangeNumGPU[0];
	void* allocatedGPUMem[2] = { NULL };
	vector<RangeQueryStateTable> stateTableRange[2];
	stateTableRange[0].resize(rangeNum * 50000);
	stateTableRange[1].resize(rangeNum * 50000);

	for (int device_idx=0; device_idx <= device_num - 1; device_idx++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		CUDA_CALL(cudaSetDevice(device_idx));
		CUDA_CALL(cudaMalloc((void**)&this->baseAddrRange[device_idx], (long long int)BIG_MEM * 1024 * 1024));
		CUDA_CALL(cudaMalloc((void**)&this->stateTableGPU[device_idx], (long long int)SMALL_MEM * 1024 * 1024));
		allocatedGPUMem[device_idx] = this->baseAddrRange[device_idx];
		threads_RQ.push_back(thread(std::mem_fn(&Grid::rangeQueryBatchGPU), this, &bounds[startIdx[device_idx
		]], rangeNumGPU[device_idx], &ResultTable[startIdx[device_idx]], resultSetSize, &stateTableRange[device_idx][0], device_idx));
	}

	timer.start();
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Dual GPU Time:" << timer.elapse() << "ms" << endl;
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		CUDA_CALL(cudaFree(allocatedGPUMem[device_idx]));
		CUDA_CALL(cudaFree(this->stateTableGPU[device_idx]));
	}
	return 0;
}


int Grid::rangeQueryBatchMultiGPUWithoutMorton(MBB* bounds, int rangeNum, CPURangeQueryResult* ResultTable, int* resultSetSize)
{
	MyTimer timer;
	int device_num = 2;
	vector<thread> threads_RQ;
	int rangeNumGPU[2];
	rangeNumGPU[0] = rangeNum / 2;
	rangeNumGPU[1] = rangeNum - rangeNumGPU[0];
	int startIdx[2];
	startIdx[0] = 0;
	startIdx[1] = rangeNumGPU[0];
	void* allocatedGPUMem[2] = { NULL };
	vector<RangeQueryStateTable> stateTableRange[2];
	stateTableRange[0].resize(rangeNum * 50000);
	stateTableRange[1].resize(rangeNum * 50000);

	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		CUDA_CALL(cudaSetDevice(device_idx));
		CUDA_CALL(cudaMalloc((void**)&this->baseAddrRange[device_idx], (long long int)BIG_MEM * 1024 * 1024));
		CUDA_CALL(cudaMalloc((void**)&this->stateTableGPU[device_idx], (long long int)SMALL_MEM * 1024 * 1024));
		allocatedGPUMem[device_idx] = this->baseAddrRange[device_idx];
		threads_RQ.push_back(thread(std::mem_fn(&Grid::rangeQueryBatchGPUWithoutMorton), this, &bounds[startIdx[device_idx
		]], rangeNumGPU[device_idx], &ResultTable[startIdx[device_idx]], resultSetSize, &stateTableRange[device_idx][0], device_idx));
	}

	timer.start();
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Dual GPU Time:" << timer.elapse() << "ms" << endl;
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		CUDA_CALL(cudaFree(allocatedGPUMem[device_idx]));
		CUDA_CALL(cudaFree(this->stateTableGPU[device_idx]));
	}
	return 0;
}

int Grid::findMatchNodeInQuadTreeGPU(QuadtreeNode* node, MBB& bound, vector<QuadtreeNode*>* cells, cudaStream_t stream, int queryID, int device_idx)
{
	int totalLevel = int(log2(this->cellnum) / log2(4));
	if (node->isLeaf)
	{
		// printf("level:%d,node:%d\n", node->level, node->NodeID);
		int startCellID = node->NodeID * int(pow(4, (totalLevel - node->level)));
		// printf("startcellID:%d\n", startCellID);
		int startIdx = this->cellPtr[startCellID].pointRangeStart;
		int pointNum = node->numPoints;

		SPoint* dataPtr = NULL;
		//如果gpu内存中没有该node的信息
		if (this->nodeAddrTable[device_idx].find(startCellID) == this->nodeAddrTable[device_idx].end())
		{
			CUDA_CALL(cudaMemcpyAsync(this->baseAddrRange[device_idx], &(this->allPoints[startIdx]), pointNum*sizeof(SPoint), cudaMemcpyHostToDevice, stream));
			dataPtr = (SPoint*)this->baseAddrRange[device_idx];
			this->nodeAddrTable[device_idx].insert(pair<int, void*>(startCellID, this->baseAddrRange[device_idx]));
			this->baseAddrRange[device_idx] = (void*)((char*)this->baseAddrRange[device_idx] + pointNum * sizeof(SPoint));
		}
		//如果有，不再复制，直接用
		else
		{
			//this->stateTableRange[device_idx]->ptr = this->nodeAddrTable[device_idx].find(startCellID)->second;
			dataPtr = (SPoint*)this->nodeAddrTable[device_idx].find(startCellID)->second;
		}

		int pointsInState = 0;
		for (int idx = 0; idx < pointNum; idx+=MAXPOINTINNODE) {
			this->stateTableRange[device_idx]->ptr = dataPtr;
			this->stateTableRange[device_idx]->xmin = bound.xmin;
			this->stateTableRange[device_idx]->xmax = bound.xmax;
			this->stateTableRange[device_idx]->ymin = bound.ymin;
			this->stateTableRange[device_idx]->ymax = bound.ymax;
			if (idx + MAXPOINTINNODE >= pointNum)
				pointsInState = pointNum - idx;
			else
				pointsInState = MAXPOINTINNODE;
			this->stateTableRange[device_idx]->candidatePointNum = pointsInState;
			this->stateTableRange[device_idx]->startIdxInAllPoints = startIdx + idx;
			this->stateTableRange[device_idx]->queryID = queryID;
			this->stateTableRange[device_idx] = this->stateTableRange[device_idx] + 1;
			this->stateTableLength[device_idx] = this->stateTableLength[device_idx] + 1;
			this->testCnt++;
			//printf("%d ", this->testCnt);
			dataPtr += pointsInState;
		}
	}
	else
	{
		if (bound.intersect(node->UL->mbb))
			findMatchNodeInQuadTreeGPU(node->UL, bound, cells, stream, queryID, device_idx);
		if (bound.intersect(node->UR->mbb))
			findMatchNodeInQuadTreeGPU(node->UR, bound, cells, stream, queryID, device_idx);
		if (bound.intersect(node->DL->mbb))
			findMatchNodeInQuadTreeGPU(node->DL, bound, cells, stream, queryID, device_idx);
		if (bound.intersect(node->DR->mbb))
			findMatchNodeInQuadTreeGPU(node->DR, bound, cells, stream, queryID, device_idx);
	}
	return 0;
}

int Grid::findMatchNodeInQuadTreeGPUWithoutMorton(QuadtreeNode * node, MBB & bound, 
	std::vector<QuadtreeNode*>* cells, cudaStream_t stream, int queryID, int device_idx, 
	vector<int>& blockOffsetInData, vector<int>& blockOffsetNum, vector<int>& blockOffsetOfOffset, 
	vector<int>& blockLength)
{
	int totalLevel = int(log2(this->cellnum) / log2(4));
	if (node->isLeaf)
	{
		// printf("level:%d,node:%d\n", node->level, node->NodeID);
		int startCellID = node->NodeID * int(pow(4, (totalLevel - node->level)));
		// printf("startcellID:%d\n", startCellID);
		int startIdx = this->cellPtr[startCellID].pointRangeStart;
		int pointNum = node->numPoints;
		SPoint* dataPtr = NULL;
		// simulate no morton encoding
		int colsOfCells = int(pow(2, (totalLevel - node->level)));


		//如果gpu内存中没有该node的信息
		if (this->nodeAddrTable[device_idx].find(startCellID) == this->nodeAddrTable[device_idx].end())
		{
			CUDA_CALL(cudaMemcpyAsync(this->baseAddrRange[device_idx], &(this->allPoints[startIdx]), pointNum * sizeof(SPoint), cudaMemcpyHostToDevice, stream));
			dataPtr = (SPoint*)this->baseAddrRange[device_idx];
			this->nodeAddrTable[device_idx].insert(pair<int, void*>(startCellID, this->baseAddrRange[device_idx]));
			this->baseAddrRange[device_idx] = (void*)((char*)this->baseAddrRange[device_idx] + pointNum * sizeof(SPoint));
		}
		//如果有，不再复制，直接用
		else
		{
			//this->stateTableRange[device_idx]->ptr = this->nodeAddrTable[device_idx].find(startCellID)->second;
			dataPtr = (SPoint*)this->nodeAddrTable[device_idx].find(startCellID)->second;
		}



		int pointsInState = 0;
		for (int idx = 0; idx < pointNum; idx += MAXPOINTINNODE) {
			this->stateTableRange[device_idx]->ptr = dataPtr;
			this->stateTableRange[device_idx]->xmin = bound.xmin;
			this->stateTableRange[device_idx]->xmax = bound.xmax;
			this->stateTableRange[device_idx]->ymin = bound.ymin;
			this->stateTableRange[device_idx]->ymax = bound.ymax;
			if (idx + MAXPOINTINNODE >= pointNum)
				pointsInState = pointNum - idx;
			else
				pointsInState = MAXPOINTINNODE;
			this->stateTableRange[device_idx]->candidatePointNum = pointsInState;
			this->stateTableRange[device_idx]->startIdxInAllPoints = startIdx + idx;
			this->stateTableRange[device_idx]->queryID = queryID;
			this->stateTableRange[device_idx] = this->stateTableRange[device_idx] + 1;
			this->stateTableLength[device_idx] = this->stateTableLength[device_idx] + 1;
			this->testCnt++;
			
			// simulate no morton encoding
			// each block for a cuda block
			blockOffsetNum.push_back(colsOfCells);
			blockOffsetOfOffset.push_back((int)blockOffsetInData.size());
			int added = 0;
			for (int i = 0; i < colsOfCells; i++) {
				blockOffsetInData.push_back(added);
				if (added + pointsInState / colsOfCells <= pointsInState)
					blockLength.push_back(pointsInState / colsOfCells);
				else
					blockLength.push_back(pointsInState - added);
			}

			//printf("%d ", this->testCnt);
			dataPtr += pointsInState;
		}
	}
	else
	{
		if (bound.intersect(node->UL->mbb))
			findMatchNodeInQuadTreeGPUWithoutMorton(node->UL, bound, cells, stream, queryID, device_idx,
				blockOffsetInData, blockOffsetNum, blockOffsetOfOffset, blockLength);
		if (bound.intersect(node->UR->mbb))
			findMatchNodeInQuadTreeGPUWithoutMorton(node->UR, bound, cells, stream, queryID, device_idx,
				blockOffsetInData, blockOffsetNum, blockOffsetOfOffset, blockLength);
		if (bound.intersect(node->DL->mbb))
			findMatchNodeInQuadTreeGPUWithoutMorton(node->DL, bound, cells, stream, queryID, device_idx,
				blockOffsetInData, blockOffsetNum, blockOffsetOfOffset, blockLength);
		if (bound.intersect(node->DR->mbb))
			findMatchNodeInQuadTreeGPUWithoutMorton(node->DR, bound, cells, stream, queryID, device_idx,
				blockOffsetInData, blockOffsetNum, blockOffsetOfOffset, blockLength);
	}
	return 0;
}









int Grid::FDCalculateParallelHandeler(priority_queue<FDwithID, vector<FDwithID>, cmp>* queue, map<int, int>* freqVectorQ)
{
	this->freqVectors.formPriorityQueue(queue, freqVectorQ);// 关键变量 std::vector<std::map<int, int>> FreqVector; //每个map代表一个轨迹 trajDB中的轨迹
	return 0;
}

// 负责找到top-k的EDR
// thread(std::mem_fn(&Grid::SimilarityMultiThreadHandler), this, queryQueue, qTra, 1, EDRCalculated, kValue, qIDnow)
int Grid::SimilarityMultiThreadHandler(priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue, Trajectory* qTra, int queryTrajNum, priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated, int kValue, int startQueryIdx)
{
	const int k = KSIMILARITY;
	int* numElemInCalculatedQueue = new int[queryTrajNum];
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	for (int qID = startQueryIdx; qID <= startQueryIdx + queryTrajNum - 1; qID++)
	{
		SPoint* queryTra = (SPoint*)malloc(sizeof(SPoint) * qTra[qID].length);
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid; // 888888
		}
		int worstNow = 9999999;
		//timer.start();
		//printf("qID:%d", qID);
		/*MyTimer tt;*/

		while (worstNow > queryQueue[qID].top().FD)// worst now 代表最大EDR
		{
			/*tt.start();*/

			int candidateTrajID[k];// m 保存ID
								   //printf("%d", worstNow);

								   // 第0步：提取topk(topm)FD
			for (int i = 0; i <= k - 1; i++)
			{
				candidateTrajID[i] = queryQueue[qID].top().traID;// FD小顶堆的ID依次pop到candidateTrajID[]中													
				queryQueue[qID].pop();
			}
			//EDR calculate

			//第一步：从AllPoints里提取出来轨迹(m条)
			SPoint** candidateTra = (SPoint**)malloc(sizeof(SPoint*) * k);
			int* candidateTraLength = (int*)malloc(sizeof(int) * k);

			for (int i = 0; i <= k - 1; i++)
			{
				candidateTra[i] = (SPoint*)malloc(sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[i]].trajLength);
				SPoint* tempPtr = candidateTra[i];
				// 遍历每条轨迹 借助Lt
				for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[i]].length - 1; subID++)
				{
					int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[i]].startIdx[subID];
					memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID]);
					//for (int cnt = 0; cnt <= this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID] - 1; cnt++) {
					//	candidateTra[i][cnt] = this->allPoints[idxInAllPoints+cnt];
					//}
					//printf("%d ", this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID]);
					tempPtr += this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID];
				}
				candidateTraLength[i] = this->cellBasedTrajectory[candidateTrajID[i]].trajLength;
			}
			//tt.stop();
			//cout << "Part3.1 time:" << tt.elapse() << endl;
			//tt.start();

			//第二步：计算EDR
			//printf("%d", qID);
			int resultReturned[k];
			this->SimilarityExecuter(queryTra, candidateTra, qTra[qID].length, candidateTraLength, k, resultReturned);
			//tt.stop();
			//cout << "Part3.3 time:" << tt.elapse() << endl;
			for (int i = 0; i <= k - 1; i++)
			{
				if (numElemInCalculatedQueue[qID - startQueryIdx] < kValue)//  qID - startQueryIdxm = 0 比较头元素 第一次计算
				{
					//直接往PQ里加
					FDwithID fd;
					fd.traID = candidateTrajID[i];
					fd.FD = resultReturned[i];
					EDRCalculated[qID].push(fd);
					numElemInCalculatedQueue[qID - startQueryIdx]++;
				}
				else // 非第一次计算
				{
					//看一下是否比PQ里更好，如果是弹出一个差的，换进去一个好的；否则不动优先队列也不更新worstNow。
					int worstInPQ = EDRCalculated[qID].top().FD;
					if (resultReturned[i] < worstInPQ)
					{
						EDRCalculated[qID].pop();
						FDwithID fd;
						fd.traID = candidateTrajID[i];
						fd.FD = resultReturned[i];
						EDRCalculated[qID].push(fd);// 注意 这里fd不一定就是在头部 优先队列会自动排序
					}
				}
			}
			worstNow = EDRCalculated[qID].top().FD;
			//printf("%d,worstNow:%d\t", qID,worstNow);
			//该轮结束，释放内存
			for (int i = 0; i <= k - 1; i++)
				free(candidateTra[i]);
			free(candidateTraLength);
			free(candidateTra); //注意candidateTra二维数组的free
		}
		//timer.stop();
		//cout << "Query Trajectory Length:" << qTra[qID].length << endl;
		//cout << "Part3 time:" << timer.elapse() << endl;
		//timer.start();
		free(queryTra); // free掉临时变量

						//timer.stop();
						//cout << "Part4 time:" << timer.elapse() << endl;
	}
	return 0;
}

int Grid::SimilarityExecuter(SPoint* queryTra, SPoint** candidateTra, int queryLength, int* candidateLength, int candSize, int* resultArray)
{
	for (int i = 0; i <= candSize - 1; i++) // 0 - m-1
	{
		// 计算m个中的每一个与Qt的EDR
		//每个DP问题
		SPoint *CPUqueryTra = queryTra, *CPUCandTra = candidateTra[i];
		int CPUqueryLength = queryLength, CPUCandLength = candidateLength[i];
		int longest = 0;

		const SPoint *tra1, *tra2;
		int len1, len2;
		//printf("%d,%d\t", len1, len2);
		if (CPUCandLength >= CPUqueryLength)
		{
			tra1 = CPUqueryTra;
			tra2 = CPUCandTra;
			len1 = CPUqueryLength;
			len2 = CPUCandLength;
		}
		else
		{
			tra1 = CPUCandTra;
			tra2 = CPUqueryTra;
			len1 = CPUCandLength;
			len2 = CPUqueryLength;
		}

		if (CPUqueryLength >= longest)
		{
			longest = CPUqueryLength;
		}
		else
		{
			longest = CPUCandLength;
		}


		int** stateTable = (int**)malloc(sizeof(int*) * (len1 + 1));
		for (int j = 0; j <= len1; j++)
		{
			stateTable[j] = (int*)malloc(sizeof(int) * (len2 + 1));
		}
		stateTable[0][0] = 0;
		for (int row = 1; row <= len1; row++)
		{
			stateTable[row][0] = row;
		}
		for (int col = 1; col <= len2; col++)
		{
			stateTable[0][col] = col;
		}

		for (int row = 1; row <= len1; row++)
		{
			for (int col = 1; col <= len2; col++)
			{
				SPoint p1 = tra1[row - 1];
				SPoint p2 = tra2[col - 1]; //这样做内存是聚集访问的吗？
				bool subcost;
				if ((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON))
				{
					subcost = 0;
				}
				else
					subcost = 1;

				int myState = 0;
				int state_ismatch = stateTable[row - 1][col - 1] + subcost;
				int state_up = stateTable[row - 1][col] + 1;
				int state_left = stateTable[row][col - 1] + 1;
				//if (state_ismatch < state_up)
				//	myState = state_ismatch;
				//else if (state_left < state_up)
				//	myState = state_left;
				//else
				//	myState = state_up;
				bool c1 = ((state_ismatch < state_up) && (state_ismatch < state_left));
				bool c2 = ((state_left < state_up) && ((state_left < state_ismatch)));
				//去除if的表达方式，是否可以提升性能？ 完全是从GPU优化角度
				myState = c1 * state_ismatch + c2 * state_left + !(c1 || c2) * state_up;

				stateTable[row][col] = myState;
				//	if (row == len1&&col == len2)
				//cout << myState << endl;
			}
		}


		resultArray[i] = stateTable[len1][len2];
		//cout << resultCPU[i] << endl;

		// 再一次二维new数组的delete
		for (int j = 0; j <= len1; j++)
		{
			free(stateTable[j]);
		}
		free(stateTable);
	}
	return 0;
}






// CPU-20-S
int Grid::SimilarityQueryBatchCPUParallel(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue)
{
	MyTimer timer;
	//思路：用不同线程处理不同数量的query
	// 小顶堆找最小的FD 非默认 注意是priority_queue的指针 类似二维
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	// 注意是map的指针 类似二维 因为priority_queue、map都是容器类型
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum];
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++) // 遍历qTra节点
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM;
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "Part1 FV time:" << timer.elapse() << endl;

	timer.start();

#ifdef	FDMULCPU
	// 多CPU线程
	// windows 不要用 太卡
	//cout << "Part2 time: start 阻塞主线程 卡\n";
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join)); // threads_FD已经指明

#else
	// 单CPU线程
	// 计算FD 存在queryQueue数组中 而每个queryQueue是priority_queue 优先队列 with CPU 单线程 多线程跑不动 可能是笔记本原因
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		this->FDCalculateParallelHandeler(&queryQueue[qID], &freqVectors[qID]);
	}

#endif

	timer.stop();
	cout << "Part2 FD time:" << timer.elapse() << endl;

	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];

	int* numElemInCalculatedQueue = new int[queryTrajNum];
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;
	const int k = KSIMILARITY;

	timer.start();
	// way1: 类似range 正而过batch 平分为 16 份
	

	// way2: 让CPU自己决定 虚拟batch个thread并行 实际怎么样呢？  注意参数 1 qID
	
	vector<thread> threads_EDR;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		threads_EDR.push_back(thread(std::mem_fn(&Grid::SimilarityMultiThreadHandler), this, queryQueue, qTra, 1, EDRCalculated, kValue, qID));
		//threads_EDR.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_EDR.begin(), threads_EDR.end(), std::mem_fn(&std::thread::join)); // threads_FD已经指明


	/*
	// 20 threads !!
	const int THREAD_CPU = 20; 
	vector<thread> threads;
	int qIDnow = 0;
	// Tq 串行，每 20 个Tq计算EDR 并行
	// recall that on GPU, Tq*m(40)
	for (int i = qIDnow; qIDnow <= queryTrajNum - 1; i++)// i 没用
	{
		for (int j = 0; j <= THREAD_CPU - 1; j++)
		{
			if (qIDnow>queryTrajNum - 1)
				break;
			threads.push_back(thread(std::mem_fn(&Grid::SimilarityMultiThreadHandler), this, queryQueue, qTra, 1, EDRCalculated, kValue, qIDnow));
			qIDnow++;
		}
		std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
		threads.clear();// 严重的性能下降 根本无法达到 2000% 占有率CPU
	}
	*/

	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}

	timer.stop();
	cout << "Part3 EDR time:" << timer.elapse() << endl;

	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;

	return 0;
}

// CPU-1-S
int Grid::SimilarityQueryBatch(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue)
{
	MyTimer timer;
	//思路：分别处理每一条查询轨迹，用不同stream并行
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum];
	//为查询构造freqVector
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //虚格子
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "Part1 time:" << timer.elapse() << endl;
	timer.start();
	//为剪枝计算Frequency Distance
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Part2 time:" << timer.elapse() << endl;
	//用一个优先队列存储当前最优结果，大顶堆，保证随时可以pop出差的结果
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];
	int* numElemInCalculatedQueue = new int[queryTrajNum]; //保存当前优先队列结果，保证优先队列大小不大于kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//准备好之后，开始做查询
	const int k = KSIMILARITY;
	timer.start();
	// check if the FD is lowerbound for all traj


	// check if the FD is lowerbound for all traj
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		SPoint* queryTra = (SPoint*)malloc(sizeof(SPoint) * qTra[qID].length);
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		int worstNow = 9999999;
		//timer.start();
		//printf("qID:%d", qID);
		/*MyTimer tt;*/
		while (worstNow > queryQueue[qID].top().FD)
		{
			/*tt.start();*/
			int candidateTrajID[k];
			//printf("%d", worstNow);
			//提取topk
			for (int i = 0; i <= k - 1; i++)
			{
				candidateTrajID[i] = queryQueue[qID].top().traID;
				//printf("%d,%d\t", queryQueue[qID].top().traID,queryQueue[qID].top().FD);
				queryQueue[qID].pop();
			}
			//EDR calculate
			//第一步：从AllPoints里提取出来轨迹
			SPoint** candidateTra = (SPoint**)malloc(sizeof(SPoint*) * k);
			int* candidateTraLength = (int*)malloc(sizeof(int) * k);
			for (int i = 0; i <= k - 1; i++)
			{
				candidateTra[i] = (SPoint*)malloc(sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[i]].trajLength);
				SPoint* tempPtr = candidateTra[i];
				for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[i]].length - 1; subID++)
				{
					int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[i]].startIdx[subID];
					memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID]);
					//for (int cnt = 0; cnt <= this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID] - 1; cnt++) {
					//	candidateTra[i][cnt] = this->allPoints[idxInAllPoints+cnt];
					//}
					//printf("%d ", this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID]);
					tempPtr += this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID];
				}
				candidateTraLength[i] = this->cellBasedTrajectory[candidateTrajID[i]].trajLength;
			}
			//tt.stop();
			//cout << "Part3.1 time:" << tt.elapse() << endl;
			//tt.start();
			//第二步：计算EDR
			//printf("%d", qID);
			int resultReturned[k];
			this->SimilarityExecuter(queryTra, candidateTra, qTra[qID].length, candidateTraLength, k, resultReturned);
			//tt.stop();
			//cout << "Part3.3 time:" << tt.elapse() << endl;
			//更新worstNow
			for (int i = 0; i <= k - 1; i++)
			{
				if (numElemInCalculatedQueue[qID] < kValue)
				{
					//直接往PQ里加
					FDwithID fd;
					fd.traID = candidateTrajID[i];
					fd.FD = resultReturned[i];
					EDRCalculated[qID].push(fd);
					numElemInCalculatedQueue[qID]++;
				}
				else
				{
					//看一下是否比PQ里更好，如果是弹出一个差的，换进去一个好的；否则不动优先队列也不更新worstNow。
					int worstInPQ = EDRCalculated[qID].top().FD;
					if (resultReturned[i] < worstInPQ)
					{
						EDRCalculated[qID].pop();
						FDwithID fd;
						fd.traID = candidateTrajID[i];
						fd.FD = resultReturned[i];
						EDRCalculated[qID].push(fd);
					}
				}
			}
			worstNow = EDRCalculated[qID].top().FD;
			//printf("%d,worstNow:%d\t", qID,worstNow);
			//该轮结束，释放内存
			for (int i = 0; i <= k - 1; i++)
				free(candidateTra[i]);
			free(candidateTraLength);
			free(candidateTra);
		}
		//timer.stop();
		//cout << "Query Trajectory Length:" << qTra[qID].length << endl;
		//cout << "Part3 time:" << timer.elapse() << endl;
		//timer.start();
		free(queryTra);

		//timer.stop();
		//cout << "Part4 time:" << timer.elapse() << endl;
	}
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}

	timer.stop();
	cout << "Part3 time:" << timer.elapse() << endl;

	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;

	return 0;
}





// GAT-S-noE
int Grid::SimilarityQueryBatchOnGPU(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue)
//所有（批量）query级别上的GPU并行
//备用思路：分别处理每一条查询轨迹，用不同stream并行
{
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU), (long long int)BIG_MEM * 1024 * 1024));
	void* whileAddrGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024));
	void* whileAddrGPUBase = whileAddrGPU;
	
	//当前分配到的地址
	void* nowAddrGPU = NULL;    
	 
	cudaStream_t defaultStream;
	cudaStreamCreate(&defaultStream);

	MyTimer timer;

	// 小顶堆找最小的FD 非默认 注意是priority_queue的指针 类似二维
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	// 注意是map的指针 类似二维 因为priority_queue、map都是容器类型
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum]; // FVtable

	//为查询构造freqVector 含义：待查询轨迹的freqVector：freqVectors[qID]  = CPU
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //虚格子
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "Part1 time:" << timer.elapse() << endl;


	timer.start();
	//为剪枝计算Frequency Distance
	//CPU并行 parallelly 是的
	vector<thread> threads_FD;
	// 计算FD 存在queryQueue数组中 而每个queryQueue是priority_queue 优先队列 = CPU
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Part2 time:" << timer.elapse() << endl;
	timer.start();
	//MyTimer tt;
	//tt.start();
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];
	
	int* numElemInCalculatedQueue = new int[queryTrajNum]; //保存当前优先队列结果，保证优先队列大小不大于kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;


	//准备好之后，开始做查询

	const int k = KSIMILARITY; //m
	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}
	//查询轨迹信息的准备：

	// Step0：*********** 保存 Tq 查询的轨迹全部信息 CPU->GPU *************
	// CPU
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	SPoint* queryTra = allQueryTra;
	int* allQueryTraOffset = new int[queryTrajNum];
	int* queryTraLength = new int[queryTrajNum];
	// GPU
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU;
	SPoint* queryTraGPUBase = queryTraGPU;
	allQueryTraOffset[0] = 0;
	printf("queryTrajNum:%d", queryTrajNum);
	printf("totalQueryTrajLength:%d\n", totalQueryTrajLength);

	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream));	
		queryTraLength[qID] = qTra[qID].length;
		queryTraGPU = queryTraGPU + qTra[qID].length;
		queryTra += qTra[qID].length;
		if (qID != queryTrajNum - 1)
			allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	nowAddrGPU = queryTraGPU;
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);

	//tt.stop();
	//cout << "Part3.0.1 time:" << tt.elapse() << endl;
	//tt.start();


	// 是最大EDR
	int* worstNow = new int[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		worstNow[qID] = 9999999;
	}
	// 置标志位 每个Tq和整个batch都有一个标志位
	bool* isFinished = new bool[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		isFinished[qID] = FALSE;
	}
	bool isAllFinished = FALSE;


	int* candidateTraLength = (int*)malloc(sizeof(int) * k * queryTrajNum);

	SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);
	int** candidateTrajID = new int*[queryTrajNum];
	for (int i = 0; i <= queryTrajNum - 1; i++)
		candidateTrajID[i] = new int[k]; // m

	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
	
	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum);
	
	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum);

	int candidateTrajNum = 0;

	map<int, void*> traID_baseAddr;

	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU;

	//tt.stop();
	//cout << "Part3.0.2 time:" << tt.elapse() << endl;

	while (!isAllFinished)
	{
		//tt.start();
		int validCandTrajNum = 0; 
		int validQueryTraNum = queryTrajNum; // NeedValidQueryTraNum
		int validQueryIdx = 0;
		SPoint* tempPtr = candidateTra;


		// Step1：*********** CPU->GPU  Tc 相关的拷贝 MAT****************
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
			if (!isFinished[qID])
			{
				if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD))
				{
					validQueryTraNum--;
					isFinished[qID] = TRUE; // valid 结束
					continue;
				}
				else
				{
					for (int i = 0; i <= k - 1; i++)
					{
						candidateTrajID[qID][i] = queryQueue[qID].top().traID;
						queryQueue[qID].pop(); 
						validCandTrajNum++;
					}
					// 存储到CPU 并到GPU
					for (int i = 0; i <= k - 1; i++)
					{
						int CandTrajID = candidateTrajID[qID][i];
						map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr.find(CandTrajID);
						// 如果轨迹还没有保存在GPU中
						if (traID_baseAddr_ITER == traID_baseAddr.end())
						{
							int pointsNumInThisCand = 0;
							SPoint* thisTrajAddr = tempPtr;

							// 遍历 Lt 重建轨迹 因为点的存储是按照cell存的
							for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++) // Lt就是cellBasedTrajectory
							{
								int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
								memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
								tempPtr += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
								pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
							}
							CUDA_CALL(cudaMemcpyAsync(candidateTraGPU, thisTrajAddr, pointsNumInThisCand*sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream));
							traID_baseAddr[candidateTrajID[qID][i]] = candidateTraGPU;
							candidateTraLength[k * validQueryIdx + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
							taskInfoTable[k * validQueryIdx + i].qID = qID;
							taskInfoTable[k * validQueryIdx + i].candTrajID = CandTrajID;
							candidateTrajOffsetTable[k * validQueryIdx + i].objectId = candidateTrajID[qID][i];
							candidateTrajOffsetTable[k * validQueryIdx + i].addr = candidateTraGPU;
							candidateOffsets[k * validQueryIdx + i] = candidateTraGPU;
							candidateTraGPU = (candidateTraGPU + pointsNumInThisCand);
							nowAddrGPU = (void*)candidateTraGPU;
						}
						else
						{
							void* baseAddrGPU = traID_baseAddr_ITER->second;
							candidateTraLength[k * validQueryIdx + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
							taskInfoTable[k * validQueryIdx + i].qID = qID;
							taskInfoTable[k * validQueryIdx + i].candTrajID = CandTrajID;
							candidateTrajOffsetTable[k * validQueryIdx + i].objectId = CandTrajID;
							candidateTrajOffsetTable[k * validQueryIdx + i].addr = baseAddrGPU;
							candidateOffsets[k * validQueryIdx + i] = (SPoint*)baseAddrGPU;
						}
					}
					validQueryIdx++;
				}
			}
		}
		//tt.stop();
		//cout << "Part3.1 time:" << tt.elapse() << endl;
		//tt.start();
		//CUDA_CALL(cudaDeviceSynchronize());

		int* candidateTraLengthGPU = (int*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateTraLengthGPU, candidateTraLength, sizeof(int)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((int*)whileAddrGPU + validCandTrajNum);
		TaskInfoTableForSimilarity* taskInfoTableGPU = (TaskInfoTableForSimilarity*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU, taskInfoTable, sizeof(TaskInfoTableForSimilarity)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU + validCandTrajNum);
		SPoint** candidateOffsetsGPU = (SPoint**)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU, candidateOffsets, sizeof(SPoint*)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((SPoint**)whileAddrGPU + validCandTrajNum);
		int* resultReturned = new int[queryTrajNum * k];
		int* resultReturnedGPU = (int*)whileAddrGPU;
		whileAddrGPU = (void*)((int*)whileAddrGPU + k * queryTrajNum);
		//CUDA_CALL(cudaMalloc((void**)resultReturnedGPU, sizeof(int)*k*queryTrajNum));
		//tt.stop();
		//cout << "Part3.2 time:" << tt.elapse() << endl;
		//tt.start();

		// Step2：***********   如果上面的分配没有错，开始计算EDR  ***************

		if (validQueryTraNum * k == validCandTrajNum)
		{

			// kernel
			EDRDistance_Batch_Handler(validCandTrajNum, taskInfoTableGPU, queryTraGPUBase, queryTraOffsetGPU, candidateOffsetsGPU, queryLengthGPU, candidateTraLengthGPU, resultReturnedGPU, &defaultStream);
			CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost, defaultStream));
			// CUDA_CALL(cudaMemcpy(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost, defaultStream));
			CUDA_CALL(cudaDeviceSynchronize());
		}
		else
		{
			printf("error in line 1007\n");
		}

		//tt.stop();
		//cout << "Part3.3 time:" << tt.elapse() << endl;
		//tt.start();

		//for (int ii = 0; ii < queryTrajNum * k; ii++)
		//	cout << resultReturned[ii] << " ";
		//cout << endl;

		// Step3：***********   并行计算结束后，更新worstNow以及写入结果 *********************
		// 保证最小k个EDR
		for (int idx = 0; idx <= k * validQueryTraNum - 1; idx++)
		{
			int qID = taskInfoTable[idx].qID;
			int i = idx % k;
			if (numElemInCalculatedQueue[qID] < kValue)
			{
				//直接往PQ里加
				FDwithID fd;
				fd.traID = candidateTrajID[qID][i];
				fd.FD = resultReturned[idx];
				// 伪FD 因为EDRCalculated->cmp函数不同
				EDRCalculated[qID].push(fd);

				numElemInCalculatedQueue[qID]++;
			}
			else
			{
				//看一下是否比PQ里更好，如果是弹出一个差的，换进去一个好的；否则不动优先队列也不更新worstNow。
				int worstInPQ = EDRCalculated[qID].top().FD; // 最大EDR
				if (resultReturned[i] < worstInPQ)
				{
					EDRCalculated[qID].pop();
					FDwithID fd;
					fd.traID = candidateTrajID[qID][i];
					fd.FD = resultReturned[idx];
					EDRCalculated[qID].push(fd);
				}
			}
		}

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
			worstNow[qID] = EDRCalculated[qID].top().FD;


		bool temp = TRUE;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			temp = temp && isFinished[qID];
		}
		isAllFinished = temp;
		delete[] resultReturned;

		whileAddrGPU = whileAddrGPUBase;

		//tt.stop();
		//cout << "Part3.4 time:" << tt.elapse() << endl;

	}


	/*
	for (int qID = 0; qID <= queryTrajNum - 1; qID++) {

		timer.start();
		int candidateTrajID[k];
		//printf("qID:%d", qID);
		while (worstNow[qID] > queryQueue[qID].top().FD) {
			//printf("%d", worstNow);
			//提取topk
			for (int i = 0; i <= k - 1; i++) {
				candidateTrajID[i] = queryQueue[qID].top().traID;
				//printf("%d,%d\t", queryQueue[qID].top().traID,queryQueue[qID].top().FD);
				queryQueue[qID].pop();
			}
			//EDR calculate
			//第一步：从AllPoints里提取出来轨迹
			SPoint **candidateTra = (SPoint**)malloc(sizeof(SPoint*)*k);

			for (int i = 0; i <= k - 1; i++) {
				candidateTra[i] = (SPoint*)malloc(sizeof(SPoint)*this->cellBasedTrajectory[candidateTrajID[i]].trajLength);
				SPoint *tempPtr = candidateTra[i];
				for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[i]].length - 1; subID++) {
					int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[i]].startIdx[subID];
					memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint)*this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID]);
					//for (int cnt = 0; cnt <= this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID] - 1; cnt++) {
					//	candidateTra[i][cnt] = this->allPoints[idxInAllPoints+cnt];
					//}
					//printf("%d ", this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID]);
					tempPtr += this->cellBasedTrajectory[candidateTrajID[i]].numOfPointInCell[subID];
				}
				candidateTraLength[i] = this->cellBasedTrajectory[candidateTrajID[i]].trajLength;
			}
			//第二步：计算EDR
			int resultReturned[k];
			this->SimilarityExecuter(queryTra, candidateTra, qTra[qID].length, candidateTraLength, k, resultReturned);
			//更新worstNow
			for (int i = 0; i <= k - 1; i++) {
				if (numElemInCalculatedQueue[qID] < kValue) {
					//直接往PQ里加
					FDwithID fd;
					fd.traID = candidateTrajID[i];
					fd.FD = resultReturned[i];
					EDRCalculated[qID].push(fd);
					numElemInCalculatedQueue[qID]++;
				}
				else {
					//看一下是否比PQ里更好，如果是弹出一个差的，换进去一个好的；否则不动优先队列也不更新worstNow。
					int worstInPQ = EDRCalculated[qID].top().FD;
					if (resultReturned[i] < worstInPQ) {
						EDRCalculated[qID].pop();
						FDwithID fd;
						fd.traID = candidateTrajID[i];
						fd.FD = resultReturned[i];
						EDRCalculated[qID].push(fd);
					}
				}
			}
			worstNow = EDRCalculated[qID].top().FD;
			//printf("%d,worstNow:%d\t", qID,worstNow);
			//该轮结束，释放内存
			for (int i = 0; i <= k - 1; i++)
				free(candidateTra[i]);
			free(candidateTraLength);
			free(candidateTra);

		}
		timer.stop();
		cout << "Query Trajectory Length:" << qTra[qID].length << endl;
		cout << "Part3 time:" << timer.elapse() << endl;
		timer.start();
		free(queryTra);
		for (int i = 0; i <= kValue - 1; i++) {
			topKSimilarityTraj[qID*kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
		timer.stop();
		cout << "Part4 time:" << timer.elapse() << endl;
	}
	*/

	timer.stop();
	cout << "Part3 time:" << timer.elapse() << endl;



	//输出结果
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}


	for (int i = 0; i <= queryTrajNum - 1; i++)
		delete[] candidateTrajID[i];
	free(taskInfoTable);
	free(candidateTrajOffsetTable);
	free(candidateOffsets);
	delete[] candidateTrajID;
	free(candidateTra);
	delete[] isFinished;
	free(candidateTraLength);
	delete[] worstNow;
	free(allQueryTra);
	delete[] allQueryTraOffset;
	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;
	delete[] queryTraLength;
	CUDA_CALL(cudaFree(baseAddrGPU));
	CUDA_CALL(cudaFree(whileAddrGPUBase));
	cudaStreamDestroy(defaultStream);
	return 0;
}

// GAT-S-E
int Grid::SimilarityQueryBatchOnGPUV2(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue)
{
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU), (long long int)BIG_MEM * 1024 * 1024));

	void* whileAddrGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024)); // 256字节 最大1024*1024个
	void* whileAddrGPUBase = whileAddrGPU;
	void* nowAddrGPU = NULL;
	cudaStream_t defaultStream;
	cudaStreamCreate(&defaultStream);

	MyTimer timer;

	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum];
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //虚格子
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "Part1 FV time:" << timer.elapse() << endl;
	timer.start();
#ifdef	FDMULCPU
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join));

#else
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		this->FDCalculateParallelHandeler(&queryQueue[qID], &freqVectors[qID]);
	}

#endif
	timer.stop();
	cout << "Part2 FD time:" << timer.elapse() << endl;
	timer.start();
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];
	

	int* numElemInCalculatedQueue = new int[queryTrajNum];
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}

	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	SPoint* queryTra = allQueryTra;
	int* allQueryTraOffset = new int[queryTrajNum];
	int* queryTraLength = new int[queryTrajNum];
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU;
	SPoint* queryTraGPUBase = queryTraGPU;
	allQueryTraOffset[0] = 0;
	printf("queryTrajNum:%d ", queryTrajNum);
	printf("totalQueryTrajLength:%d\n", totalQueryTrajLength);
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream));
		
		queryTraLength[qID] = qTra[qID].length;
		queryTraGPU = queryTraGPU + qTra[qID].length;
		queryTra += qTra[qID].length;
		if (qID != queryTrajNum - 1)
			allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	nowAddrGPU = queryTraGPU;
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);


	int* worstNow = new int[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		worstNow[qID] = 9999999;
	}

	bool* isFinished = new bool[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		isFinished[qID] = FALSE;
	}
	bool isAllFinished = FALSE;


	// const int k = KSIMILARITY; //m
	
	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU; 
	int candidateTrajNum = 0;
	map<int, void*> traID_baseAddr;

	while (!isAllFinished)
	{
		int validQueryTraNum2 = queryTrajNum;

		int qID = 0;
		for (qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum2--;
			if (!isFinished[qID])
			{

				if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD))
				{
					validQueryTraNum2--;
					isFinished[qID] = TRUE;
					continue;
				}

			}
		}

		// definition of m

		int dynamicm = KSIMILARITY;
		if (validQueryTraNum2 > ALPHA * queryTrajNum)
			dynamicm = KSIMILARITY * queryTrajNum / validQueryTraNum2; // *
		else
			dynamicm = KSIMILARITY; // truncation 

		if (dynamicm > MAXDYNAMICM) dynamicm = KSIMILARITY;

		if (dynamicm % 2) dynamicm += 1;


		SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * dynamicm * queryTrajNum * MAXLENGTH);

		int** candidateTrajID = new int*[queryTrajNum];
		for (int i = 0; i <= queryTrajNum - 1; i++)
			candidateTrajID[i] = new int[dynamicm];
		SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * dynamicm * queryTrajNum);
		int* candidateTraLength = (int*)malloc(sizeof(int) * dynamicm * queryTrajNum);
					
		TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * dynamicm * queryTrajNum);
		OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * dynamicm * queryTrajNum);
		
		int validCandTrajNum = 0; 
		int validQueryTraNum = queryTrajNum; // NeedValidQueryTraNum
		int validQueryIdx = 0;
		SPoint* tempPtr = candidateTra;// CPU

		int ktmp = dynamicm; 
		for (qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
			if (!isFinished[qID])
			{
				
				if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD))
				{
					validQueryTraNum--;
					isFinished[qID] = TRUE;
					continue;
				}
				else
				{
					for (int i = 0; i <= dynamicm - 1; i++)
					{
						if(! queryQueue[qID].empty()){
							candidateTrajID[qID][i] = queryQueue[qID].top().traID;
							queryQueue[qID].pop();
							validCandTrajNum++;
						}
						else {
							validCandTrajNum++;
							candidateTrajID[qID][i] = -1;
							ktmp--;
						}
					}

					for (int i = 0; i <= dynamicm - 1; i++)
					{
						int CandTrajID = candidateTrajID[qID][i];
						if( CandTrajID > 0){
							map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr.find(CandTrajID);
							if (traID_baseAddr_ITER == traID_baseAddr.end())
							{
								int pointsNumInThisCand = 0;
								SPoint* thisTrajAddr = tempPtr;
								for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++)
								{
									int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
									memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
									tempPtr += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
									pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
								}
								CUDA_CALL(cudaMemcpyAsync(candidateTraGPU, thisTrajAddr, pointsNumInThisCand*sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream));
								traID_baseAddr[candidateTrajID[qID][i]] = candidateTraGPU; 
								candidateTraLength[dynamicm * validQueryIdx + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
								taskInfoTable[dynamicm * validQueryIdx + i].qID = qID;
								taskInfoTable[dynamicm * validQueryIdx + i].candTrajID = CandTrajID;
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].objectId = candidateTrajID[qID][i];
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].addr = candidateTraGPU;
								candidateOffsets[dynamicm * validQueryIdx + i] = candidateTraGPU;
								candidateTraGPU = (candidateTraGPU + pointsNumInThisCand);
								nowAddrGPU = (void*)candidateTraGPU;
							}
							else
							{
								void* baseAddrGPU = traID_baseAddr_ITER->second;
								candidateTraLength[dynamicm * validQueryIdx + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
								taskInfoTable[dynamicm * validQueryIdx + i].qID = qID;
								taskInfoTable[dynamicm * validQueryIdx + i].candTrajID = CandTrajID;
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].objectId = CandTrajID;
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].addr = baseAddrGPU;
								candidateOffsets[dynamicm * validQueryIdx + i] = (SPoint*)baseAddrGPU;
							}
						}
						else {
							candidateTraLength[dynamicm * validQueryIdx + i] = -1;
							taskInfoTable[dynamicm * validQueryIdx + i].qID = -1;
							taskInfoTable[dynamicm * validQueryIdx + i].candTrajID = -1;
							candidateTrajOffsetTable[dynamicm * validQueryIdx + i].objectId = -1;
							candidateTrajOffsetTable[dynamicm * validQueryIdx + i].addr = NULL;
							candidateOffsets[dynamicm * validQueryIdx + i] = NULL;

						}
					
					}
					
					validQueryIdx++;
				}
			}
		}
		int* candidateTraLengthGPU = (int*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateTraLengthGPU, candidateTraLength, sizeof(int)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((int*)whileAddrGPU + validCandTrajNum);

		TaskInfoTableForSimilarity* taskInfoTableGPU = (TaskInfoTableForSimilarity*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU, taskInfoTable, sizeof(TaskInfoTableForSimilarity)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU + validCandTrajNum);
		SPoint** candidateOffsetsGPU = (SPoint**)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU, candidateOffsets, sizeof(SPoint*)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((SPoint**)whileAddrGPU + validCandTrajNum);


		int* resultReturned = new int[queryTrajNum * dynamicm];
		int* resultReturnedGPU = (int*)whileAddrGPU;
		whileAddrGPU = (void*)((int*)whileAddrGPU + dynamicm * queryTrajNum);

		if (validQueryTraNum * dynamicm == validCandTrajNum)
		{


			EDRDistance_Batch_Handler(validCandTrajNum, taskInfoTableGPU, queryTraGPUBase, queryTraOffsetGPU, candidateOffsetsGPU, queryLengthGPU, candidateTraLengthGPU, resultReturnedGPU, &defaultStream);
			
			//CUDA_CALL(cudaDeviceSynchronize());
			
			//CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost, defaultStream));
			
			if (resultReturned != NULL && resultReturnedGPU != NULL){
				CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*dynamicm*queryTrajNum, cudaMemcpyDeviceToHost));
			}
			else { 
				cout << "resultReturned = NULLL"<<endl;
			}

			//CUDA_CALL(cudaDeviceSynchronize());
		}
		else
		{
			printf("error in line 1007\n");
		}


		for (int idx = 0; idx <= dynamicm * validQueryTraNum - 1; idx++)
		{
			int qID = taskInfoTable[idx].qID;
			// cout << "qID= " << qID << endl;
			if( qID >= 0 ){
				int i = idx % dynamicm;
				if (numElemInCalculatedQueue[qID] < kValue)
				{
					FDwithID fd;
					fd.traID = candidateTrajID[qID][i];
					fd.FD = resultReturned[idx];
					EDRCalculated[qID].push(fd);

					numElemInCalculatedQueue[qID]++;
				}
				else
				{
					int worstInPQ = EDRCalculated[qID].top().FD;
					if (resultReturned[i] < worstInPQ)
					{
						EDRCalculated[qID].pop();
						FDwithID fd;
						fd.traID = candidateTrajID[qID][i];
						fd.FD = resultReturned[idx];
						EDRCalculated[qID].push(fd);
					}
				}
			}
		}

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
			worstNow[qID] = EDRCalculated[qID].top().FD;


		bool temp = TRUE;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			temp = temp && isFinished[qID];
		}
		isAllFinished = temp;

		delete[] resultReturned;

		whileAddrGPU = whileAddrGPUBase;


		for (int i = 0; i <= queryTrajNum - 1; i++)
			delete[] candidateTrajID[i];
		free(taskInfoTable);
		free(candidateTrajOffsetTable);
		free(candidateOffsets);
		delete[] candidateTrajID;
		free(candidateTra);
		free(candidateTraLength);

	}


	timer.stop();
	cout << "Part3 EDR time:" << timer.elapse() << endl;



	//输出结果
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}

	delete[] isFinished;
	delete[] worstNow;
	free(allQueryTra);
	delete[] allQueryTraOffset;
	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;
	delete[] queryTraLength;
	CUDA_CALL(cudaFree(baseAddrGPU));
	CUDA_CALL(cudaFree(whileAddrGPUBase));
	cudaStreamDestroy(defaultStream);
	return 0;
}






// GAT-S-noE
int Grid::SimilarityQueryBatchOnGPUV3(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue, int device_idx)
{
	MyTimer totaltimer;
	totaltimer.start();

	CUDA_CALL(cudaSetDevice(device_idx));

	void* baseAddrGPU2 = NULL;// 不要使用全局内存地址
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU2), (long long int)BIG_MEM * 1024 * 1024));

	void* whileAddrGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024));
	void* whileAddrGPUBase = whileAddrGPU;

	//当前分配到的地址
	void* nowAddrGPU = NULL;

	cudaStream_t defaultStream;
	cudaStreamCreate(&defaultStream);

	MyTimer timer;

	// 小顶堆找最小的FD 非默认 注意是priority_queue的指针 类似二维
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	// 注意是map的指针 类似二维 因为priority_queue、map都是容器类型
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum]; // FVtable

	//为查询构造freqVector 含义：待查询轨迹的freqVector：freqVectors[qID]  with CPU
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //虚格子
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "device: " << device_idx << "Part1 FV time:" << timer.elapse() << endl;


	timer.start();

#ifdef	FDMULCPU

	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join)); // threads_FD已经指明

#else
	// 单CPU线程
	// 计算FD 存在queryQueue数组中 而每个queryQueue是priority_queue 优先队列 with CPU 单线程 多线程跑不动 可能是笔记本原因
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		this->FDCalculateParallelHandeler(&queryQueue[qID], &freqVectors[qID]);
	}
	// 返回到 queryQueue[]

#endif
	timer.stop();
	cout << "device: " << device_idx << "Part2 FD time:" << timer.elapse() << endl;

	timer.start();

	//MyTimer tt;
	//tt.start();

	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];

	int* numElemInCalculatedQueue = new int[queryTrajNum]; //保存当前优先队列结果，保证优先队列大小不大于kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//准备好之后，开始做查询


	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}
	//查询轨迹信息的准备：

	// Step0：*********** 保存Tq查询的轨迹全部信息 CPU->GPU *************

	// CPU
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	SPoint* queryTra = allQueryTra;
	int* allQueryTraOffset = new int[queryTrajNum];
	int* queryTraLength = new int[queryTrajNum];
	// GPU
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU2;
	SPoint* queryTraGPUBase = queryTraGPU;

	allQueryTraOffset[0] = 0;
	printf("device:%d\tqueryTrajNum:%d\t", device_idx, queryTrajNum);
	printf("totalQueryTrajLength:%d\n", totalQueryTrajLength);
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// 存到queryTra中
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		// 内存拷贝 待查询轨迹 CPU->GPU
		CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream));

		queryTraLength[qID] = qTra[qID].length;
		queryTraGPU = queryTraGPU + qTra[qID].length;
		queryTra += qTra[qID].length;
		if (qID != queryTrajNum - 1)
			allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	nowAddrGPU = queryTraGPU;
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);

	//tt.stop();
	//cout << "Part3.0.1 time:" << tt.elapse() << endl;
	//tt.start();

	//第一步：循环，剪枝
	int* worstNow = new int[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		worstNow[qID] = 9999999;
	}

	bool* isFinished = new bool[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		isFinished[qID] = FALSE;
	}
	bool isAllFinished = FALSE;


	const int k = KSIMILARITY;
	SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);
	int** candidateTrajID = new int*[queryTrajNum];
	for (int i = 0; i <= queryTrajNum - 1; i++)
		candidateTrajID[i] = new int[k];

	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum); // 二维数组 指针

	int* candidateTraLength = (int*)malloc(sizeof(int) * k * queryTrajNum);
	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
	
	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum); // 没用到？ candidateOffsets 重合
	
	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU;
	
	int candidateTrajNum = 0;
	map<int, void*> traID_baseAddr;

	//tt.stop();
	//cout << "Part3.0.2 time:" << tt.elapse() << endl;

	while (!isAllFinished)
	{
		//tt.start();
		int validCandTrajNum = 0;
		int validQueryTraNum = queryTrajNum;
		int validQueryIdx = 0;
		SPoint* tempPtr = candidateTra;
		int qID = 0;


		int ktmp = KSIMILARITY;

		// Step1：*********** CPU->GPU  Tc 相关的拷贝 MAT****************
		for (qID = 0; qID <= queryTrajNum - 1; qID++) 
		{
			if (isFinished[qID])
				validQueryTraNum--;
			if (!isFinished[qID])
			{

				if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD))
				{
					validQueryTraNum--;
					isFinished[qID] = TRUE;
					continue;
				}
				else
				{
					for (int i = 0; i <= k - 1; i++)
					{
						if (!queryQueue[qID].empty()) {
							candidateTrajID[qID][i] = queryQueue[qID].top().traID;

							queryQueue[qID].pop(); 
							validCandTrajNum++;
						}
						else {
							validCandTrajNum++; 
							candidateTrajID[qID][i] = -1; 
							ktmp--;
						}
					}
					for (int i = 0; i <= k - 1; i++)
					{
						int CandTrajID = candidateTrajID[qID][i];
						if (CandTrajID > 0) {
							map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr.find(CandTrajID);
							if (traID_baseAddr_ITER == traID_baseAddr.end())
							{
								int pointsNumInThisCand = 0;
								SPoint* thisTrajAddr = tempPtr;
								for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++) // Lt就是cellBasedTrajectory table
								{
									int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
									memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
									tempPtr += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
									pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
								}
								CUDA_CALL(cudaMemcpyAsync(candidateTraGPU, thisTrajAddr, pointsNumInThisCand*sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream));
								
								traID_baseAddr[candidateTrajID[qID][i]] = candidateTraGPU;// 记录

								candidateTraLength[k * validQueryIdx + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
								taskInfoTable[k * validQueryIdx + i].qID = qID;
								taskInfoTable[k * validQueryIdx + i].candTrajID = CandTrajID;
								candidateTrajOffsetTable[k * validQueryIdx + i].objectId = candidateTrajID[qID][i];// k * validQueryIdx + i
								candidateTrajOffsetTable[k * validQueryIdx + i].addr = candidateTraGPU;
								candidateOffsets[k * validQueryIdx + i] = candidateTraGPU;

								candidateTraGPU = (candidateTraGPU + pointsNumInThisCand);
								nowAddrGPU = (void*)candidateTraGPU;
							}
							else
							{
								void* baseAddrGPU2 = traID_baseAddr_ITER->second;
								candidateTraLength[k * validQueryIdx + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
								taskInfoTable[k * validQueryIdx + i].qID = qID;
								taskInfoTable[k * validQueryIdx + i].candTrajID = CandTrajID;
								candidateTrajOffsetTable[k * validQueryIdx + i].objectId = CandTrajID;
								candidateTrajOffsetTable[k * validQueryIdx + i].addr = baseAddrGPU2;
								candidateOffsets[k * validQueryIdx + i] = (SPoint*)baseAddrGPU2;
							}
						}
						else {
							candidateTraLength[k * validQueryIdx + i] = -1;
							taskInfoTable[k * validQueryIdx + i].qID = -1;
							taskInfoTable[k * validQueryIdx + i].candTrajID = -1;
							candidateTrajOffsetTable[k * validQueryIdx + i].objectId = -1;
							candidateTrajOffsetTable[k * validQueryIdx + i].addr = NULL;
							candidateOffsets[k * validQueryIdx + i] = NULL;

						}

					}
					validQueryIdx++;
				}
			}
		}
		//tt.stop();
		//cout << "Part3.1 time:" << tt.elapse() << endl;
		//tt.start();


		//构建candidateTraj完成，构建candidateTrajLength

		int* candidateTraLengthGPU = (int*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateTraLengthGPU, candidateTraLength, sizeof(int)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((int*)whileAddrGPU + validCandTrajNum);
		TaskInfoTableForSimilarity* taskInfoTableGPU = (TaskInfoTableForSimilarity*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU, taskInfoTable, sizeof(TaskInfoTableForSimilarity)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU + validCandTrajNum);
		SPoint** candidateOffsetsGPU = (SPoint**)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU, candidateOffsets, sizeof(SPoint*)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((SPoint**)whileAddrGPU + validCandTrajNum);


		int* resultReturned = new int[queryTrajNum * k];
		int* resultReturnedGPU = (int*)whileAddrGPU;
		whileAddrGPU = (void*)((int*)whileAddrGPU + k * queryTrajNum);

		//tt.stop();
		//cout << "Part3.2 time:" << tt.elapse() << endl;
		//tt.start();

		//如果上面的分配没有错，开始计算EDR
		if (validQueryTraNum * k == validCandTrajNum)
		{

			// kernel
			EDRDistance_Batch_Handler(validCandTrajNum, taskInfoTableGPU, queryTraGPUBase, queryTraOffsetGPU, candidateOffsetsGPU, queryLengthGPU, candidateTraLengthGPU, resultReturnedGPU, &defaultStream);

			//CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost, defaultStream));

			//CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost));
			if (resultReturned != NULL && resultReturnedGPU != NULL) {
				CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost));
			}
			else {
				cout << "resultReturned = NULL" << endl;
			}
		
			//CUDA_CALL(cudaDeviceSynchronize());	
		}
		else
		{
			printf("error in line 1007\n");
		}

		//tt.stop();
		//cout << "Part3.3 time:" << tt.elapse() << endl;
		//tt.start();

		//并行计算结束后，更新worstNow以及写入结果
		for (int idx = 0; idx <= k * validQueryTraNum - 1; idx++)
		{
			int qID = taskInfoTable[idx].qID;
			if (qID >= 0) {
				int i = idx % k;
				if (numElemInCalculatedQueue[qID] < kValue)
				{
					FDwithID fd;
					fd.traID = candidateTrajID[qID][i];
					fd.FD = resultReturned[idx];
					EDRCalculated[qID].push(fd);

					numElemInCalculatedQueue[qID]++;
				}
				else
				{
					int worstInPQ = EDRCalculated[qID].top().FD;
					if (resultReturned[i] < worstInPQ)
					{
						EDRCalculated[qID].pop();
						FDwithID fd;
						fd.traID = candidateTrajID[qID][i];
						fd.FD = resultReturned[idx];
						EDRCalculated[qID].push(fd);
					}
				}
			}
		}

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
			worstNow[qID] = EDRCalculated[qID].top().FD;


		bool temp = TRUE;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			temp = temp && isFinished[qID];
		}
		isAllFinished = temp;

		delete[] resultReturned;

		whileAddrGPU = whileAddrGPUBase;
		//tt.stop();
		//cout << "Part3.4 time:" << tt.elapse() << endl;
	}


	timer.stop();
	cout << "device: " << device_idx << "Part3 EDR time:" << timer.elapse() << endl;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}
	
	
	for (int i = 0; i <= queryTrajNum - 1; i++)
		delete[] candidateTrajID[i];
	free(taskInfoTable);
	free(candidateTrajOffsetTable);
	free(candidateOffsets);
	delete[] candidateTrajID;
	free(candidateTra);
	delete[] isFinished;
	free(candidateTraLength);
	delete[] worstNow;
	free(allQueryTra);
	delete[] allQueryTraOffset;
	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;
	delete[] queryTraLength;
	CUDA_CALL(cudaFree(baseAddrGPU2));
	CUDA_CALL(cudaFree(whileAddrGPUBase));
	cudaStreamDestroy(defaultStream);

	totaltimer.stop();
	cout << "Single GPU Time:" << totaltimer.elapse() << endl;

	return 0;
}

// GAT-S-E
int Grid::SimilarityQueryBatchOnGPUV4(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue, int device_idx)
{
	
	CUDA_CALL(cudaSetDevice(device_idx));

	void* baseAddrGPU2 = NULL;
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU2), (long long int)BIG_MEM * 1024 * 1024));

	void* whileAddrGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024));
	void* whileAddrGPUBase = whileAddrGPU;

	void* nowAddrGPU = NULL;

	cudaStream_t defaultStream;
	cudaStreamCreate(&defaultStream);

	MyTimer timer;


	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum];

	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //虚格子
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "device: " << device_idx <<" Part1 FV time: " << timer.elapse() << endl;


	timer.start();
#ifdef	FDMULCPU

	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join)); // threads_FD已经指明

#else
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		this->FDCalculateParallelHandeler(&queryQueue[qID], &freqVectors[qID]);
	}

#endif
	timer.stop();
	cout << "device: " << device_idx << " Part2 FD time: " << timer.elapse() << endl;


	timer.start();
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];


	int* numElemInCalculatedQueue = new int[queryTrajNum];
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;


	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}


	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	SPoint* queryTra = allQueryTra;
	int* allQueryTraOffset = new int[queryTrajNum];
	int* queryTraLength = new int[queryTrajNum];
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU2;
	SPoint* queryTraGPUBase = queryTraGPU;
	allQueryTraOffset[0] = 0;
	printf("device:%d\tqueryTrajNum:%d\t", device_idx, queryTrajNum);
	printf("totalQueryTrajLength:%d\n", totalQueryTrajLength);
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// 存到queryTra中
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		// 内存拷贝 待查询轨迹 CPU->GPU
		CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream));

		queryTraLength[qID] = qTra[qID].length;
		queryTraGPU = queryTraGPU + qTra[qID].length;
		queryTra += qTra[qID].length;
		if (qID != queryTrajNum - 1)
			allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	nowAddrGPU = queryTraGPU;
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	int* worstNow = new int[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		worstNow[qID] = 9999999;
	}
	bool* isFinished = new bool[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		isFinished[qID] = FALSE;
	}
	bool isAllFinished = FALSE;


	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU;
	int candidateTrajNum = 0;

	map<int, void*> traID_baseAddr;



	while (!isAllFinished)
	{
		int validQueryTraNum2 = queryTrajNum;

		int qID = 0; 
		for (qID = 0; qID <= queryTrajNum - 1; qID++) 
		{
			if (isFinished[qID])
				validQueryTraNum2--; 
			if (!isFinished[qID])
			{

				if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD)) 
				{
					// 停止pop qID完成搜索
					validQueryTraNum2--;
					isFinished[qID] = TRUE; // valid 结束
					continue;
				}

			}
		}


		int dynamicm = KSIMILARITY;
		if (validQueryTraNum2 > ALPHA * KSIMILARITY)
			dynamicm = KSIMILARITY * queryTrajNum / validQueryTraNum2; // *
		else
			dynamicm = KSIMILARITY;
		if (dynamicm > MAXDYNAMICM) dynamicm = KSIMILARITY;

		if (dynamicm % 2) dynamicm += 1; 

		SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * dynamicm * queryTrajNum * MAXLENGTH);

		int** candidateTrajID = new int*[queryTrajNum];
		for (int i = 0; i <= queryTrajNum - 1; i++)
			candidateTrajID[i] = new int[dynamicm]; 
		SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * dynamicm * queryTrajNum);

		int* candidateTraLength = (int*)malloc(sizeof(int) * dynamicm * queryTrajNum);
		TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * dynamicm * queryTrajNum);
		OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * dynamicm * queryTrajNum);

		int validCandTrajNum = 0;
		int validQueryTraNum = queryTrajNum; // NeedValidQueryTraNum
		int validQueryIdx = 0;
		SPoint* tempPtr = candidateTra;// CPU

		int ktmp = dynamicm;


		// Step1：*********** CPU->GPU  Tc 相关的拷贝 MAT****************

		for (qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
			if (!isFinished[qID])
			{

				if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD)) 
				{
					validQueryTraNum--;
					isFinished[qID] = TRUE;
					continue;
				}
				else
				{
					for (int i = 0; i <= dynamicm - 1; i++)
					{
						if (!queryQueue[qID].empty()) {
							candidateTrajID[qID][i] = queryQueue[qID].top().traID; 
							queryQueue[qID].pop();
							validCandTrajNum++; 
						}
						else {
							// cout << "queryQueue[qID].empty()" << endl;
							validCandTrajNum++;
							candidateTrajID[qID][i] = -1;
							ktmp--;
						}
					}

					// 存储到CPU 并到GPU
					for (int i = 0; i <= dynamicm - 1; i++)
					{
						int CandTrajID = candidateTrajID[qID][i];
						if (CandTrajID > 0) {
							map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr.find(CandTrajID);
							if (traID_baseAddr_ITER == traID_baseAddr.end())
							{
								int pointsNumInThisCand = 0;
								SPoint* thisTrajAddr = tempPtr;
								for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++) // Lt就是cellBasedTrajectory
								{
									int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
									// 从allPoints拷贝点到CPU
									memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
									tempPtr += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
									pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
								}
								// 把这条轨迹提取到candidateTraGPU中
								CUDA_CALL(cudaMemcpyAsync(candidateTraGPU, thisTrajAddr, pointsNumInThisCand*sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream));
								traID_baseAddr[candidateTrajID[qID][i]] = candidateTraGPU;
								candidateTraLength[dynamicm * validQueryIdx + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
								taskInfoTable[dynamicm * validQueryIdx + i].qID = qID;
								taskInfoTable[dynamicm * validQueryIdx + i].candTrajID = CandTrajID;
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].objectId = candidateTrajID[qID][i];// 就是循环开始的CandTrajID
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].addr = candidateTraGPU;
								candidateOffsets[dynamicm * validQueryIdx + i] = candidateTraGPU;
								candidateTraGPU = (candidateTraGPU + pointsNumInThisCand);
								nowAddrGPU = (void*)candidateTraGPU;
							}
							// 如果该轨迹已经复制进了gpu里面，那么只需要按照该轨迹id更新表就行了
							else
							{
								void* baseAddrGPU = traID_baseAddr_ITER->second;
								candidateTraLength[dynamicm * validQueryIdx + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
								taskInfoTable[dynamicm * validQueryIdx + i].qID = qID;
								taskInfoTable[dynamicm * validQueryIdx + i].candTrajID = CandTrajID;
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].objectId = CandTrajID;
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].addr = baseAddrGPU;
								candidateOffsets[dynamicm * validQueryIdx + i] = (SPoint*)baseAddrGPU;
							}
						}
						else {
							candidateTraLength[dynamicm * validQueryIdx + i] = -1;
							taskInfoTable[dynamicm * validQueryIdx + i].qID = -1;
							taskInfoTable[dynamicm * validQueryIdx + i].candTrajID = -1;
							candidateTrajOffsetTable[dynamicm * validQueryIdx + i].objectId = -1;
							candidateTrajOffsetTable[dynamicm * validQueryIdx + i].addr = NULL;
							candidateOffsets[dynamicm * validQueryIdx + i] = NULL;

						}

					}

					validQueryIdx++;
				}
			}
		}

		int* candidateTraLengthGPU = (int*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateTraLengthGPU, candidateTraLength, sizeof(int)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((int*)whileAddrGPU + validCandTrajNum);
		TaskInfoTableForSimilarity* taskInfoTableGPU = (TaskInfoTableForSimilarity*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU, taskInfoTable, sizeof(TaskInfoTableForSimilarity)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU + validCandTrajNum);

		SPoint** candidateOffsetsGPU = (SPoint**)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU, candidateOffsets, sizeof(SPoint*)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((SPoint**)whileAddrGPU + validCandTrajNum);
		int* resultReturned = new int[queryTrajNum * dynamicm]; 
		int* resultReturnedGPU = (int*)whileAddrGPU;
		whileAddrGPU = (void*)((int*)whileAddrGPU + dynamicm * queryTrajNum);

		if (validQueryTraNum * dynamicm == validCandTrajNum)
		{

			EDRDistance_Batch_Handler(validCandTrajNum, taskInfoTableGPU, queryTraGPUBase, queryTraOffsetGPU, candidateOffsetsGPU, queryLengthGPU, candidateTraLengthGPU, resultReturnedGPU, &defaultStream);

			if (resultReturned != NULL && resultReturnedGPU != NULL) {
				CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*dynamicm*queryTrajNum, cudaMemcpyDeviceToHost));
			}
			else {
				cout << "resultReturned = NULL" << endl;
			}

		}
		else
		{
			printf("error in line 1007\n");
		}



		for (int idx = 0; idx <= dynamicm * validQueryTraNum - 1; idx++)
		{
			int qID = taskInfoTable[idx].qID;
			if (qID >= 0) {// 边界条件
				int i = idx % dynamicm;
				if (numElemInCalculatedQueue[qID] < kValue)
				{
					//直接往PQ里加
					FDwithID fd;
					fd.traID = candidateTrajID[qID][i];
					fd.FD = resultReturned[idx];
					EDRCalculated[qID].push(fd);

					numElemInCalculatedQueue[qID]++;
				}
				else
				{
					int worstInPQ = EDRCalculated[qID].top().FD;
					if (resultReturned[i] < worstInPQ)
					{
						EDRCalculated[qID].pop();
						FDwithID fd;
						fd.traID = candidateTrajID[qID][i];
						fd.FD = resultReturned[idx];
						EDRCalculated[qID].push(fd);
					}
				}
			}
		}

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
			worstNow[qID] = EDRCalculated[qID].top().FD;


		bool temp = TRUE;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			temp = temp && isFinished[qID];
		}
		isAllFinished = temp;

		delete[] resultReturned;

		whileAddrGPU = whileAddrGPUBase;

		for (int i = 0; i <= queryTrajNum - 1; i++)
			delete[] candidateTrajID[i];
		free(taskInfoTable);
		free(candidateTrajOffsetTable);
		free(candidateOffsets);
		delete[] candidateTrajID;
		free(candidateTra);
		free(candidateTraLength);

	}


	timer.stop();
	cout << "device: "<< device_idx << " Part3 EDR time: " << timer.elapse() << endl;



	//输出结果
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}


	delete[] isFinished;
	delete[] worstNow;
	free(allQueryTra);
	delete[] allQueryTraOffset;
	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;
	delete[] queryTraLength;
	CUDA_CALL(cudaFree(baseAddrGPU2));
	CUDA_CALL(cudaFree(whileAddrGPUBase));
	cudaStreamDestroy(defaultStream);
	return 0;
}


int Grid::SimilarityQueryBatchOnGPUNoMAT(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue, int device_idx)
{
	MyTimer totaltimer;
	totaltimer.start();

	CUDA_CALL(cudaSetDevice(device_idx));

	void* baseAddrGPU2 = NULL;// 不要使用全局内存地址
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU2), (long long int)BIG_MEM * 1024 * 1024));

	void* whileAddrGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024));
	void* whileAddrGPUBase = whileAddrGPU;

	//当前分配到的地址
	void* nowAddrGPU = NULL;

	cudaStream_t defaultStream;
	cudaStreamCreate(&defaultStream);

	MyTimer timer;

	// 小顶堆找最小的FD 非默认 注意是priority_queue的指针 类似二维
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	// 注意是map的指针 类似二维 因为priority_queue、map都是容器类型
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum]; // FVtable

																  //为查询构造freqVector 含义：待查询轨迹的freqVector：freqVectors[qID]  with CPU
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //虚格子
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "device: " << device_idx << "Part1 FV time:" << timer.elapse() << endl;


	timer.start();

#ifdef	FDMULCPU

	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join)); // threads_FD已经指明

#else
	// 单CPU线程
	// 计算FD 存在queryQueue数组中 而每个queryQueue是priority_queue 优先队列 with CPU 单线程 多线程跑不动 可能是笔记本原因
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		this->FDCalculateParallelHandeler(&queryQueue[qID], &freqVectors[qID]);
	}
	// 返回到 queryQueue[]

#endif
	timer.stop();
	cout << "device: " << device_idx << "Part2 FD time:" << timer.elapse() << endl;

	timer.start();

	//MyTimer tt;
	//tt.start();

	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];

	int* numElemInCalculatedQueue = new int[queryTrajNum]; //保存当前优先队列结果，保证优先队列大小不大于kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//准备好之后，开始做查询


	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}
	//查询轨迹信息的准备：

	// Step0：*********** 保存Tq查询的轨迹全部信息 CPU->GPU *************

	// CPU
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	SPoint* queryTra = allQueryTra;
	int* allQueryTraOffset = new int[queryTrajNum];
	int* queryTraLength = new int[queryTrajNum];
	// GPU
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU2;
	SPoint* queryTraGPUBase = queryTraGPU;

	allQueryTraOffset[0] = 0;
	printf("device:%d\tqueryTrajNum:%d\t", device_idx, queryTrajNum);
	printf("totalQueryTrajLength:%d\n", totalQueryTrajLength);
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// 存到queryTra中
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		// 内存拷贝 待查询轨迹 CPU->GPU
		CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream));

		queryTraLength[qID] = qTra[qID].length;
		queryTraGPU = queryTraGPU + qTra[qID].length;
		queryTra += qTra[qID].length;
		if (qID != queryTrajNum - 1)
			allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	nowAddrGPU = queryTraGPU;
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);

	//tt.stop();
	//cout << "Part3.0.1 time:" << tt.elapse() << endl;
	//tt.start();

	//第一步：循环，剪枝
	int* worstNow = new int[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		worstNow[qID] = 9999999;
	}

	bool* isFinished = new bool[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		isFinished[qID] = FALSE;
	}
	bool isAllFinished = FALSE;


	const int k = KSIMILARITY;
	SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);
	int** candidateTrajID = new int*[queryTrajNum];
	for (int i = 0; i <= queryTrajNum - 1; i++)
		candidateTrajID[i] = new int[k];

	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum); // 二维数组 指针

	int* candidateTraLength = (int*)malloc(sizeof(int) * k * queryTrajNum);
	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);

	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum); // 没用到？ candidateOffsets 重合

	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU;

	int candidateTrajNum = 0;

	//map<int, void*> traID_baseAddr;

	//tt.stop();
	//cout << "Part3.0.2 time:" << tt.elapse() << endl;

	while (!isAllFinished)
	{
		//tt.start();
		int validCandTrajNum = 0;
		int validQueryTraNum = queryTrajNum;
		int validQueryIdx = 0;
		SPoint* tempPtr = candidateTra;
		int qID = 0;


		int ktmp = KSIMILARITY;

		// Step1：*********** CPU->GPU  Tc 相关的拷贝 MAT****************
		for (qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
			if (!isFinished[qID])
			{

				if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD))
				{
					validQueryTraNum--;
					isFinished[qID] = TRUE;
					continue;
				}
				else
				{
					for (int i = 0; i <= k - 1; i++)
					{
						if (!queryQueue[qID].empty()) {
							candidateTrajID[qID][i] = queryQueue[qID].top().traID;

							queryQueue[qID].pop();
							validCandTrajNum++;
						}
						else {
							validCandTrajNum++;
							candidateTrajID[qID][i] = -1;
							ktmp--;
						}
					}
					for (int i = 0; i <= k - 1; i++)
					{
						int CandTrajID = candidateTrajID[qID][i];
						if (CandTrajID > 0) {
							//map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr.find(CandTrajID);
							//if (traID_baseAddr_ITER == traID_baseAddr.end())
							{
								int pointsNumInThisCand = 0;
								SPoint* thisTrajAddr = tempPtr;
								for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++) // Lt就是cellBasedTrajectory table
								{
									int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
									memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
									tempPtr += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
									pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
								}
								CUDA_CALL(cudaMemcpyAsync(candidateTraGPU, thisTrajAddr, pointsNumInThisCand * sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream));

								//traID_baseAddr[candidateTrajID[qID][i]] = candidateTraGPU;// 记录

								candidateTraLength[k * validQueryIdx + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
								taskInfoTable[k * validQueryIdx + i].qID = qID;
								taskInfoTable[k * validQueryIdx + i].candTrajID = CandTrajID;
								candidateTrajOffsetTable[k * validQueryIdx + i].objectId = candidateTrajID[qID][i];// k * validQueryIdx + i
								candidateTrajOffsetTable[k * validQueryIdx + i].addr = candidateTraGPU;
								candidateOffsets[k * validQueryIdx + i] = candidateTraGPU;

								candidateTraGPU = (candidateTraGPU + pointsNumInThisCand);
								nowAddrGPU = (void*)candidateTraGPU;
							}
							//else
							//{
							//	void* baseAddrGPUx = traID_baseAddr_ITER->second;
							//	candidateTraLength[k * validQueryIdx + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
							//	taskInfoTable[k * validQueryIdx + i].qID = qID;
							//	taskInfoTable[k * validQueryIdx + i].candTrajID = CandTrajID;
							//	candidateTrajOffsetTable[k * validQueryIdx + i].objectId = CandTrajID;
							//	candidateTrajOffsetTable[k * validQueryIdx + i].addr = baseAddrGPUx;
							//	candidateOffsets[k * validQueryIdx + i] = (SPoint*)baseAddrGPUx;
							//}
						}
						else {
							candidateTraLength[k * validQueryIdx + i] = -1;
							taskInfoTable[k * validQueryIdx + i].qID = -1;
							taskInfoTable[k * validQueryIdx + i].candTrajID = -1;
							candidateTrajOffsetTable[k * validQueryIdx + i].objectId = -1;
							candidateTrajOffsetTable[k * validQueryIdx + i].addr = NULL;
							candidateOffsets[k * validQueryIdx + i] = NULL;

						}

					}
					validQueryIdx++;
				}
			}
		}
		//tt.stop();
		//cout << "Part3.1 time:" << tt.elapse() << endl;
		//tt.start();


		//构建candidateTraj完成，构建candidateTrajLength

		int* candidateTraLengthGPU = (int*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateTraLengthGPU, candidateTraLength, sizeof(int)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((int*)whileAddrGPU + validCandTrajNum);
		TaskInfoTableForSimilarity* taskInfoTableGPU = (TaskInfoTableForSimilarity*)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU, taskInfoTable, sizeof(TaskInfoTableForSimilarity)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU + validCandTrajNum);
		SPoint** candidateOffsetsGPU = (SPoint**)whileAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU, candidateOffsets, sizeof(SPoint*)*validCandTrajNum, cudaMemcpyHostToDevice, defaultStream));
		whileAddrGPU = (void*)((SPoint**)whileAddrGPU + validCandTrajNum);


		int* resultReturned = new int[queryTrajNum * k];
		int* resultReturnedGPU = (int*)whileAddrGPU;
		whileAddrGPU = (void*)((int*)whileAddrGPU + k * queryTrajNum);

		//tt.stop();
		//cout << "Part3.2 time:" << tt.elapse() << endl;
		//tt.start();

		//如果上面的分配没有错，开始计算EDR
		if (validQueryTraNum * k == validCandTrajNum)
		{

			// kernel
			EDRDistance_Batch_Handler(validCandTrajNum, taskInfoTableGPU, queryTraGPUBase, queryTraOffsetGPU, candidateOffsetsGPU, queryLengthGPU, candidateTraLengthGPU, resultReturnedGPU, &defaultStream);

			//CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost, defaultStream));

			//CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost));
			if (resultReturned != NULL && resultReturnedGPU != NULL) {
				CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU, sizeof(int)*k*queryTrajNum, cudaMemcpyDeviceToHost));
			}
			else {
				cout << "resultReturned = NULL" << endl;
			}

			//CUDA_CALL(cudaDeviceSynchronize());	
		}
		else
		{
			printf("error in line 1007\n");
		}

		//tt.stop();
		//cout << "Part3.3 time:" << tt.elapse() << endl;
		//tt.start();

		//并行计算结束后，更新worstNow以及写入结果
		for (int idx = 0; idx <= k * validQueryTraNum - 1; idx++)
		{
			int qID = taskInfoTable[idx].qID;
			if (qID >= 0) {
				int i = idx % k;
				if (numElemInCalculatedQueue[qID] < kValue)
				{
					FDwithID fd;
					fd.traID = candidateTrajID[qID][i];
					fd.FD = resultReturned[idx];
					EDRCalculated[qID].push(fd);

					numElemInCalculatedQueue[qID]++;
				}
				else
				{
					int worstInPQ = EDRCalculated[qID].top().FD;
					if (resultReturned[i] < worstInPQ)
					{
						EDRCalculated[qID].pop();
						FDwithID fd;
						fd.traID = candidateTrajID[qID][i];
						fd.FD = resultReturned[idx];
						EDRCalculated[qID].push(fd);
					}
				}
			}
		}

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
			worstNow[qID] = EDRCalculated[qID].top().FD;


		bool temp = TRUE;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			temp = temp && isFinished[qID];
		}
		isAllFinished = temp;

		delete[] resultReturned;

		whileAddrGPU = whileAddrGPUBase;
		//tt.stop();
		//cout << "Part3.4 time:" << tt.elapse() << endl;
	}


	timer.stop();
	cout << "device: " << device_idx << "Part3 EDR time:" << timer.elapse() << endl;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}


	for (int i = 0; i <= queryTrajNum - 1; i++)
		delete[] candidateTrajID[i];
	free(taskInfoTable);
	free(candidateTrajOffsetTable);
	free(candidateOffsets);
	delete[] candidateTrajID;
	free(candidateTra);
	delete[] isFinished;
	free(candidateTraLength);
	delete[] worstNow;
	free(allQueryTra);
	delete[] allQueryTraOffset;
	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;
	delete[] queryTraLength;
	CUDA_CALL(cudaFree(baseAddrGPU2));
	CUDA_CALL(cudaFree(whileAddrGPUBase));
	cudaStreamDestroy(defaultStream);

	totaltimer.stop();
	cout << "Single GPU Time:" << totaltimer.elapse() << endl;

	return 0;
}





// FineGrained GAT-S-noE
int Grid::SimilarityQueryBatchOnMultiGPU(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue)

{

	MyTimer totaltimer;
	totaltimer.start();
	// from single GPU to multiple GPUs
	// 分离FV FD 计算

	void* baseAddrSimi[2] = {NULL};
	void* whileAddrGPU[2] = {NULL};
	void* whileAddrGPUBase[2];
	void* nowAddrGPU[2] = {NULL};;
	int num_devices;
	cudaStream_t defaultStream[2];	

	CUDA_CALL(cudaGetDeviceCount(&num_devices));
	//num_devices = 2;
	printf("num of GPU:%d\n",num_devices);

	for (int device_idx = 0; device_idx <= num_devices - 1; device_idx++)
	{
		CUDA_CALL(cudaSetDevice(device_idx));
		CUDA_CALL(cudaMalloc((void**)(&baseAddrSimi[device_idx]), (long long int)BIG_MEM * 1024 * 1024));
		CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU[device_idx]), (long long int)SMALL_MEM * 1024 * 1024));
		whileAddrGPUBase[device_idx] = whileAddrGPU[device_idx];
		nowAddrGPU[device_idx] = baseAddrSimi[device_idx];
		cudaStreamCreate(&defaultStream[device_idx]);
	}

	//当前分配到的地址


	MyTimer timer;
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum]; 
	//为查询构造freqVector
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //虚格子
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "Part1 FV time:" << timer.elapse() << endl;
	
	timer.start();
	//为剪枝计算Frequency Distance
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Part2 FD time:" << timer.elapse() << endl;

	//用一个优先队列存储当前最优结果，大顶堆，保证随时可以pop出差的结果
	timer.start();
	//MyTimer tt;
	//tt.start();
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];
	int* numElemInCalculatedQueue = new int[queryTrajNum]; //保存当前优先队列结果，保证优先队列大小不大于kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//准备好之后，开始做查询
	const int k = KSIMILARITY;

	/*
	int queryNumEachGPU[2];
	queryNumEachGPU[0] = queryTrajNum / 2;
	queryNumEachGPU[1] = queryTrajNum - queryNumEachGPU[0];
	int queryStartIdx[2];
	queryStartIdx[0] = 0;
	queryStartIdx[1] = queryNumEachGPU[0];
	int queryEndIdx[2];
	queryEndIdx[0] = queryNumEachGPU[0] - 1;
	queryEndIdx[1] = queryTrajNum - 1;

	for (int device_idx = 0; device_idx < 2;device_idx++)
	{
		CUDA_CALL(cudaSetDevice(device_idx));
		int totalQueryTrajLength = 0;
		for (int qID = queryStartIdx[device_idx]; qID <= queryEndIdx[device_idx]; qID++)
		{
			totalQueryTrajLength += qTra[qID].length;
		}
		//查询轨迹信息的准备：
		//保存查询的轨迹
		SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
		//保存在allQueryTra中各个轨迹的offset（起始地址）
		int* allQueryTraOffset = new int[queryNumEachGPU[device_idx]];
		SPoint* queryTra = allQueryTra;
		SPoint* queryTraGPU = (SPoint*)baseAddrSimi[device_idx];
		//这个才是保存所有queryTra的基址
		SPoint* queryTraGPUBase = queryTraGPU;
		int* queryTraLength = new int[queryNumEachGPU[device_idx]];
		allQueryTraOffset[0] = 0;
		printf("queryTrajNum:%d", queryNumEachGPU[device_idx]);
		printf("totalQueryTrajLength:%d", totalQueryTrajLength);
		for (int qID = queryStartIdx[device_idx]; qID <= queryEndIdx[device_idx]; qID++)
		{
			int idxInThisGPU = qID - queryStartIdx[device_idx];
			for (int i = 0; i <= qTra[qID].length - 1; i++)
			{
				queryTra[i].x = qTra[qID].points[i].lon;
				queryTra[i].y = qTra[qID].points[i].lat;
				queryTra[i].tID = qTra[qID].tid;
			}
			CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream[device_idx]));
			queryTraLength[idxInThisGPU] = qTra[qID].length;
			queryTraGPU = queryTraGPU + qTra[qID].length;
			queryTra += qTra[qID].length;
			if (qID != queryTrajNum - 1)
				allQueryTraOffset[idxInThisGPU + 1] = allQueryTraOffset[idxInThisGPU] + qTra[qID].length;
		}
		nowAddrGPU[device_idx] = queryTraGPU;
		// queryTraOffsetGPU是保存queryTra中offset的基地址
		int* queryTraOffsetGPU = (int*)nowAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryNumEachGPU[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU + queryNumEachGPU[device_idx]);

		//构建queryLength
		int* queryLengthGPU = (int*)nowAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryNumEachGPU[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU + queryNumEachGPU[device_idx]);
	}

	*/

	//-----------------------------------------------
	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}

	// 查询 Tq 轨迹信息的准备：
	// 对2个GPU拷贝完全一样数据
	// 保存查询的轨迹
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//保存在allQueryTra中各个轨迹的offset（起始地址）
	int* allQueryTraOffset = new int[queryTrajNum];
	SPoint* queryTra = allQueryTra;

	SPoint* queryTraGPU[2];
	SPoint* queryTraGPUBase[2];
	for (int device_idx = 0; device_idx <= num_devices - 1; device_idx++)
	{
		queryTraGPU[device_idx] = (SPoint*)nowAddrGPU[device_idx];
		//这个才是保存所有queryTra的基址
		queryTraGPUBase[device_idx] = queryTraGPU[device_idx];
	}

	int* queryTraLength = new int[queryTrajNum];
	allQueryTraOffset[0] = 0;
	printf("queryTrajNum:%d ", queryTrajNum);
	printf("totalQueryTrajLength:%d\n", totalQueryTrajLength);
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			CUDA_CALL(cudaMemcpyAsync(queryTraGPU[device_idx], queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		queryTraLength[qID] = qTra[qID].length;
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			queryTraGPU[device_idx] = queryTraGPU[device_idx] + qTra[qID].length;
		queryTra += qTra[qID].length;
		if (qID != queryTrajNum - 1)
			allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	int* queryTraOffsetGPU[2];
	int* queryLengthGPU[2];
	for (int device_idx = 0; device_idx <= 1; device_idx++)
	{
		nowAddrGPU[device_idx] = queryTraGPU[device_idx];
		// queryTraOffsetGPU是保存queryTra中offset的基地址
		queryTraOffsetGPU[device_idx] = (int*)nowAddrGPU[device_idx];
		CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU[device_idx], allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU[device_idx] + queryTrajNum);

		//构建queryLength
		queryLengthGPU[device_idx] = (int*)nowAddrGPU[device_idx];
		CUDA_CALL(cudaMemcpyAsync(queryLengthGPU[device_idx], queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU[device_idx] + queryTrajNum);
	}


	//----------------------------------------------------

	/*
	 *
	// 单GPU版本备份
	//-----------------------------------------------
	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}
	//查询轨迹信息的准备：
	//保存查询的轨迹
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//保存在allQueryTra中各个轨迹的offset（起始地址）
	int* allQueryTraOffset = new int[queryTrajNum];
	SPoint* queryTra = allQueryTra;
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU;
	//这个才是保存所有queryTra的基址
	SPoint* queryTraGPUBase = queryTraGPU;
	int* queryTraLength = new int[queryTrajNum];
	allQueryTraOffset[0] = 0;
	printf("queryTrajNum:%d", queryTrajNum);
	printf("totalQueryTrajLength:%d", totalQueryTrajLength);
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream));
		queryTraLength[qID] = qTra[qID].length;
		queryTraGPU = queryTraGPU + qTra[qID].length;
		queryTra += qTra[qID].length;
		if (qID != queryTrajNum - 1)
			allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	nowAddrGPU = queryTraGPU;
	// queryTraOffsetGPU是保存queryTra中offset的基地址
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);

	//构建queryLength
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	//----------------------------------------------------
	//单GPU版本备份
	*/


	//tt.stop();
	//cout << "Part3.0.1 time:" << tt.elapse() << endl;
	//tt.start();
	//第一步：循环，剪枝
	int* worstNow = new int[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		worstNow[qID] = 9999999;
	}
	
	bool* isFinished = new bool[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		isFinished[qID] = FALSE;
	}
	bool isAllFinished = FALSE;

	
	SPoint* candidateTra[2];
	SPoint* candidateTraGPU[2];
	
	map<int, void*> traID_baseAddr[2];

	int candidateTrajNum[2] = { 0 };
	int* candidateTraLength[2];
	TaskInfoTableForSimilarity* taskInfoTable[2];
	OffsetTable* candidateTrajOffsetTable[2];
	SPoint** candidateOffsets[2];
	int** candidateTrajID = new int*[queryTrajNum];

	for (int i = 0; i <= queryTrajNum - 1; i++)
		candidateTrajID[i] = new int[k];

	for (int device_idx = 0; device_idx <= 1; device_idx++)
	{
		CUDA_CALL(cudaSetDevice(device_idx));// no need!!
		candidateTra[device_idx] = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);
		taskInfoTable[device_idx] = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
		candidateTrajOffsetTable[device_idx] = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum);
		candidateOffsets[device_idx] = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum);
		candidateTraGPU[device_idx] = (SPoint*)nowAddrGPU[device_idx];
		candidateTraLength[device_idx] = (int*)malloc(sizeof(int) * k * queryTrajNum);
	}

	/*
	//单GPU版本备份

	SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);

	//保存qid、candID和在candidateTran中的offset的对应关系
	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
	//轨迹唯一，保存在id下轨迹的id和baseAddr （有必要保存轨迹的id吗？）
	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum);
	//保存candidateOffset的基地址们
	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum);
	int candidateTrajNum = 0;
	// traID和在candidateTrajOffsetTable中的idx的对应关系map，主要用于判断轨迹是否已经复制到gpu
	map<int, void*> traID_baseAddr;
	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU;
	
	//单GPU版本备份
	*/


	//tt.stop();
	//cout << "Part3.0.2 time:" << tt.elapse() << endl;

	while (!isAllFinished)
	{
		//tt.start();
		//还有待计算的候选轨迹的总数目，相当于任务数目
		int validCandTrajNum[2] = { 0 };


		// 目前还没有计算完成的查询轨迹的数目
		int validQueryTraNum = queryTrajNum;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
		}


		// 在本轮计算内的第几个轨迹
		int validQueryIdx = 0;
		int validQueryIdxGPU[2] = { 0 };
		int queryEachGPU[2];
		


		// 多GPU加速的关键 每次大while 40*40(m) -> 20*40 20*40 上面的任务个数比之前少一半 （固定k）所以快
		// 任务不可控 数据依赖性太大
		// 剩40 20+20
		// 剩30 15+15 Tq长度完全不可控

		queryEachGPU[0] = validQueryTraNum / 2;
		queryEachGPU[1] = validQueryTraNum - queryEachGPU[0];



		// 线性地址
		SPoint* tempPtr[2];
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			tempPtr[device_idx] = candidateTra[device_idx];

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			// 选择 GPU
			int device_idx = 0;
			if (validQueryIdx < queryEachGPU[0])
				device_idx = 0;
			else
				device_idx = 1;

			// 未完成则拷贝Tc
			if (!isFinished[qID])
			{

				//提取topk作为待查轨迹
				for (int i = 0; i <= k - 1; i++)
				{
					candidateTrajID[qID][i] = queryQueue[qID].top().traID;
					queryQueue[qID].pop();
					validCandTrajNum[device_idx]++;
				}

				for (int i = 0; i <= k - 1; i++)
				{
					int CandTrajID = candidateTrajID[qID][i];
					map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr[device_idx].find(CandTrajID);
					// 如果轨迹还没有保存在GPU中
					if (traID_baseAddr_ITER == traID_baseAddr[device_idx].end())
					{
						int pointsNumInThisCand = 0;
						SPoint* thisTrajAddr = tempPtr[device_idx];
						for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++)
						{
							int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
							memcpy(tempPtr[device_idx], &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
							tempPtr[device_idx] += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
							pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
						}
						// 把这条轨迹提取到candidateTraGPU中
						CUDA_CALL(cudaMemcpyAsync(candidateTraGPU[device_idx], thisTrajAddr, pointsNumInThisCand*sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream[device_idx]));
						traID_baseAddr[device_idx][candidateTrajID[qID][i]] = candidateTraGPU[device_idx];
						//仅保存要计算的query的candidateTraLength
						candidateTraLength[device_idx][k * validQueryIdxGPU[device_idx] + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
						//仅保存要计算的query的offset
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].qID = qID;
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].candTrajID = CandTrajID;
						// 保存轨迹对应的addr
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].objectId = candidateTrajID[qID][i];
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].addr = candidateTraGPU[device_idx];
						candidateOffsets[device_idx][k * validQueryIdxGPU[device_idx] + i] = candidateTraGPU[device_idx];
						//地址往前移动，便于下一次复制
						candidateTraGPU[device_idx] = (candidateTraGPU[device_idx] + pointsNumInThisCand);
						// nowAddrGPU 始终是指向下一个空闲的GPU地址
						nowAddrGPU[device_idx] = (void*)candidateTraGPU[device_idx];
					}
					// 如果该轨迹已经复制进了gpu里面，那么只需要按照该轨迹id更新表就行了
					else
					{
						void* baseAddrGPU = traID_baseAddr_ITER->second;
						//仅保存要计算的query的candidateTraLength
						candidateTraLength[device_idx][k * validQueryIdxGPU[device_idx] + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
						//仅保存要计算的query的offset
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].qID = qID;
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].candTrajID = CandTrajID;
						// 保存轨迹对应的addr
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].objectId = CandTrajID;
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].addr = baseAddrGPU;
						candidateOffsets[device_idx][k * validQueryIdxGPU[device_idx] + i] = (SPoint*)baseAddrGPU;
					}
				}
				validQueryIdx++;
				validQueryIdxGPU[device_idx]++;
			}
		}


		int* candidateTraLengthGPU[2];
		TaskInfoTableForSimilarity* taskInfoTableGPU[2];
		SPoint** candidateOffsetsGPU[2];
		int* resultReturnedGPU[2];
		int* resultReturned = new int[queryTrajNum * k];
		for (int device_idx = 0; device_idx <= 1; device_idx++) {
			candidateTraLengthGPU[device_idx] = (int*)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(candidateTraLengthGPU[device_idx], candidateTraLength[device_idx], sizeof(int)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU 始终是指向下一个空闲的GPU地址
			whileAddrGPU[device_idx] = (void*)((int*)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//构建TaskInfoTable
			taskInfoTableGPU[device_idx] = (TaskInfoTableForSimilarity*)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU[device_idx], taskInfoTable[device_idx], sizeof(TaskInfoTableForSimilarity)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU 始终是指向下一个空闲的GPU地址
			whileAddrGPU[device_idx] = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//复制candidate的地址到gpu中
			candidateOffsetsGPU[device_idx]= (SPoint**)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU[device_idx], candidateOffsets[device_idx], sizeof(SPoint*)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU 始终是指向下一个空闲的GPU地址
			whileAddrGPU[device_idx] = (void*)((SPoint**)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//构建candidateTraj和candidateLength完成，准备并行Similarity search
			//只需要查询isFinished为false的queryTra，至于是哪些可以直接看offsetTableCandidateTra
			
			resultReturnedGPU[device_idx] = (int*)whileAddrGPU[device_idx];
			whileAddrGPU[device_idx] = (void*)((int*)whileAddrGPU[device_idx] + k * queryTrajNum);
		}
		//CUDA_CALL(cudaMalloc((void**)resultReturnedGPU, sizeof(int)*k*queryTrajNum));


		//tt.stop();
		//cout << "Part3.2 time:" << tt.elapse() << endl;
		//tt.start();
		//如果上面的分配没有错，开始计算EDR
		if (validQueryTraNum * k == (validCandTrajNum[0]+ validCandTrajNum[1]))
		{
			// GPU 1 2 工作顺序 由于不同的工作流 stream0 stream1 保证多GPU 同时运行
			// 利用 stream 进行 2 GPU并行
			for (int device_idx = 0; device_idx <= 1; device_idx++){// for 循环进行控制 先GPU-1 再 GPU-2 不过很快 没什么影响！！！
				CUDA_CALL(cudaSetDevice(device_idx));
				if(validCandTrajNum[device_idx]==0) // 边界条件
					continue;
				EDRDistance_Batch_Handler(validCandTrajNum[device_idx], taskInfoTableGPU[device_idx], queryTraGPUBase[device_idx], queryTraOffsetGPU[device_idx], candidateOffsetsGPU[device_idx], queryLengthGPU[device_idx], candidateTraLengthGPU[device_idx], resultReturnedGPU[device_idx], &defaultStream[device_idx]);
			}
			CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU[0], sizeof(int)*k*queryEachGPU[0], cudaMemcpyDeviceToHost, defaultStream[0])); // 注意这里有stream[0]的标志
			CUDA_CALL(cudaMemcpyAsync(resultReturned + k*queryEachGPU[0], resultReturnedGPU[1], sizeof(int)*k*queryEachGPU[1], cudaMemcpyDeviceToHost, defaultStream[1])); // 注意这里有stream[1]的标志
			
			// 这里存在GPU1 GPU2 的隐式同步？？应该是

		}
		else
		{
			printf("error in line 1007\n");
		}

		//tt.stop();
		//cout << "Part3.3 time:" << tt.elapse() << endl;
		//tt.start();

		for (int idx = 0; idx <= k * validQueryTraNum - 1; idx++)
		{
			int device_idx = 0;
			int idxInTaskInfoTable = 0;
			if (idx / k < queryEachGPU[0])
			{
				device_idx = 0;
				idxInTaskInfoTable = idx;
			}
			else
			{
				device_idx = 1;
				idxInTaskInfoTable = idx - k*queryEachGPU[0];
			}
			int qID = taskInfoTable[device_idx][idxInTaskInfoTable].qID;
			
			int i = idx % k;
			if (numElemInCalculatedQueue[qID] < kValue)
			{
				//直接往PQ里加
				FDwithID fd;
				fd.traID = candidateTrajID[qID][i];
				fd.FD = resultReturned[idx];
				EDRCalculated[qID].push(fd);
				numElemInCalculatedQueue[qID]++;
			}
			else
			{
				//看一下是否比PQ里更好，如果是弹出一个差的，换进去一个好的；否则不动优先队列也不更新worstNow。
				int worstInPQ = EDRCalculated[qID].top().FD;
				if (resultReturned[i] < worstInPQ)
				{
					EDRCalculated[qID].pop();
					FDwithID fd;
					fd.traID = candidateTrajID[qID][i];
					fd.FD = resultReturned[idx];
					EDRCalculated[qID].push(fd);
				}
			}
		}

		// 在这里判断 比较清晰 可读性
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			worstNow[qID] = EDRCalculated[qID].top().FD;
			if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD))
			{
				isFinished[qID] = TRUE;
			}
		}

		bool temp = TRUE;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			temp = temp && isFinished[qID];
		}
		isAllFinished = temp;
		delete[] resultReturned;

		//GPU指针回到while开始的地方
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			whileAddrGPU[device_idx] = whileAddrGPUBase[device_idx];

		//tt.stop();
		//cout << "Part3.4 time:" << tt.elapse() << endl;
	} // end of BIG while


	timer.stop();
	cout << "Part3 EDR time:" << timer.elapse() << endl;

	//输出结果
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}

	for (int i = 0; i <= queryTrajNum - 1; i++)
		delete[] candidateTrajID[i];
	for (int device_idx = 0; device_idx <= 1; device_idx++) {
		free(taskInfoTable[device_idx]);
		free(candidateTrajOffsetTable[device_idx]);
		free(candidateOffsets[device_idx]);
		free(candidateTra[device_idx]);
		free(candidateTraLength[device_idx]);
		cudaStreamDestroy(defaultStream[device_idx]);
		CUDA_CALL(cudaFree(baseAddrSimi[device_idx]));
		CUDA_CALL(cudaFree(whileAddrGPU[device_idx]));
	}
	
	delete[] candidateTrajID;
	delete[] isFinished;
	delete[] worstNow;
	free(allQueryTra);
	delete[] allQueryTraOffset;
	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;
	delete[] queryTraLength;
	
	totaltimer.stop();
	cout << "Fined Dual GPU Time:" << totaltimer.elapse() << endl;

	return 0;
}



int Grid::SimilarityQueryBatchOnMultiGPUNoMAT(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue)

{

	MyTimer totaltimer;
	totaltimer.start();
	// from single GPU to multiple GPUs
	// 分离FV FD 计算

	void* baseAddrSimi[2] = { NULL };
	void* whileAddrGPU[2] = { NULL };
	void* whileAddrGPUBase[2];
	void* nowAddrGPU[2] = { NULL };;
	int num_devices;
	cudaStream_t defaultStream[2];

	CUDA_CALL(cudaGetDeviceCount(&num_devices));
	//num_devices = 2;
	printf("num of GPU:%d\n", num_devices);

	for (int device_idx = 0; device_idx <= num_devices - 1; device_idx++)
	{
		CUDA_CALL(cudaSetDevice(device_idx));
		CUDA_CALL(cudaMalloc((void**)(&baseAddrSimi[device_idx]), (long long int)BIG_MEM * 1024 * 1024));
		CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU[device_idx]), (long long int)SMALL_MEM * 1024 * 1024));
		whileAddrGPUBase[device_idx] = whileAddrGPU[device_idx];
		nowAddrGPU[device_idx] = baseAddrSimi[device_idx];
		cudaStreamCreate(&defaultStream[device_idx]);
	}

	//当前分配到的地址


	MyTimer timer;
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum];
	//为查询构造freqVector
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //虚格子
			map<int, int>::iterator iter = freqVectors[qID].find(vituralCellNo);
			if (iter == freqVectors[qID].end())
			{
				freqVectors[qID].insert(pair<int, int>(vituralCellNo, 1));
			}
			else
			{
				freqVectors[qID][vituralCellNo] = freqVectors[qID][vituralCellNo] + 1;
			}
		}
	}
	timer.stop();
	cout << "Part1 FV time:" << timer.elapse() << endl;

	timer.start();
	//为剪枝计算Frequency Distance
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Part2 FD time:" << timer.elapse() << endl;

	//用一个优先队列存储当前最优结果，大顶堆，保证随时可以pop出差的结果
	timer.start();
	//MyTimer tt;
	//tt.start();
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];
	int* numElemInCalculatedQueue = new int[queryTrajNum]; //保存当前优先队列结果，保证优先队列大小不大于kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//准备好之后，开始做查询
	const int k = KSIMILARITY;

	/*
	int queryNumEachGPU[2];
	queryNumEachGPU[0] = queryTrajNum / 2;
	queryNumEachGPU[1] = queryTrajNum - queryNumEachGPU[0];
	int queryStartIdx[2];
	queryStartIdx[0] = 0;
	queryStartIdx[1] = queryNumEachGPU[0];
	int queryEndIdx[2];
	queryEndIdx[0] = queryNumEachGPU[0] - 1;
	queryEndIdx[1] = queryTrajNum - 1;

	for (int device_idx = 0; device_idx < 2;device_idx++)
	{
	CUDA_CALL(cudaSetDevice(device_idx));
	int totalQueryTrajLength = 0;
	for (int qID = queryStartIdx[device_idx]; qID <= queryEndIdx[device_idx]; qID++)
	{
	totalQueryTrajLength += qTra[qID].length;
	}
	//查询轨迹信息的准备：
	//保存查询的轨迹
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//保存在allQueryTra中各个轨迹的offset（起始地址）
	int* allQueryTraOffset = new int[queryNumEachGPU[device_idx]];
	SPoint* queryTra = allQueryTra;
	SPoint* queryTraGPU = (SPoint*)baseAddrSimi[device_idx];
	//这个才是保存所有queryTra的基址
	SPoint* queryTraGPUBase = queryTraGPU;
	int* queryTraLength = new int[queryNumEachGPU[device_idx]];
	allQueryTraOffset[0] = 0;
	printf("queryTrajNum:%d", queryNumEachGPU[device_idx]);
	printf("totalQueryTrajLength:%d", totalQueryTrajLength);
	for (int qID = queryStartIdx[device_idx]; qID <= queryEndIdx[device_idx]; qID++)
	{
	int idxInThisGPU = qID - queryStartIdx[device_idx];
	for (int i = 0; i <= qTra[qID].length - 1; i++)
	{
	queryTra[i].x = qTra[qID].points[i].lon;
	queryTra[i].y = qTra[qID].points[i].lat;
	queryTra[i].tID = qTra[qID].tid;
	}
	CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream[device_idx]));
	queryTraLength[idxInThisGPU] = qTra[qID].length;
	queryTraGPU = queryTraGPU + qTra[qID].length;
	queryTra += qTra[qID].length;
	if (qID != queryTrajNum - 1)
	allQueryTraOffset[idxInThisGPU + 1] = allQueryTraOffset[idxInThisGPU] + qTra[qID].length;
	}
	nowAddrGPU[device_idx] = queryTraGPU;
	// queryTraOffsetGPU是保存queryTra中offset的基地址
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryNumEachGPU[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
	nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU + queryNumEachGPU[device_idx]);

	//构建queryLength
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryNumEachGPU[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
	nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU + queryNumEachGPU[device_idx]);
	}

	*/

	//-----------------------------------------------
	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}

	// 查询 Tq 轨迹信息的准备：
	// 对2个GPU拷贝完全一样数据
	// 保存查询的轨迹
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//保存在allQueryTra中各个轨迹的offset（起始地址）
	int* allQueryTraOffset = new int[queryTrajNum];
	SPoint* queryTra = allQueryTra;

	SPoint* queryTraGPU[2];
	SPoint* queryTraGPUBase[2];
	for (int device_idx = 0; device_idx <= num_devices - 1; device_idx++)
	{
		queryTraGPU[device_idx] = (SPoint*)nowAddrGPU[device_idx];
		//这个才是保存所有queryTra的基址
		queryTraGPUBase[device_idx] = queryTraGPU[device_idx];
	}

	int* queryTraLength = new int[queryTrajNum];
	allQueryTraOffset[0] = 0;
	printf("queryTrajNum:%d ", queryTrajNum);
	printf("totalQueryTrajLength:%d\n", totalQueryTrajLength);
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			CUDA_CALL(cudaMemcpyAsync(queryTraGPU[device_idx], queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		queryTraLength[qID] = qTra[qID].length;
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			queryTraGPU[device_idx] = queryTraGPU[device_idx] + qTra[qID].length;
		queryTra += qTra[qID].length;
		if (qID != queryTrajNum - 1)
			allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	int* queryTraOffsetGPU[2];
	int* queryLengthGPU[2];
	for (int device_idx = 0; device_idx <= 1; device_idx++)
	{
		nowAddrGPU[device_idx] = queryTraGPU[device_idx];
		// queryTraOffsetGPU是保存queryTra中offset的基地址
		queryTraOffsetGPU[device_idx] = (int*)nowAddrGPU[device_idx];
		CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU[device_idx], allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU[device_idx] + queryTrajNum);

		//构建queryLength
		queryLengthGPU[device_idx] = (int*)nowAddrGPU[device_idx];
		CUDA_CALL(cudaMemcpyAsync(queryLengthGPU[device_idx], queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU[device_idx] + queryTrajNum);
	}


	//----------------------------------------------------

	/*
	*
	// 单GPU版本备份
	//-----------------------------------------------
	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
	totalQueryTrajLength += qTra[qID].length;
	}
	//查询轨迹信息的准备：
	//保存查询的轨迹
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//保存在allQueryTra中各个轨迹的offset（起始地址）
	int* allQueryTraOffset = new int[queryTrajNum];
	SPoint* queryTra = allQueryTra;
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU;
	//这个才是保存所有queryTra的基址
	SPoint* queryTraGPUBase = queryTraGPU;
	int* queryTraLength = new int[queryTrajNum];
	allQueryTraOffset[0] = 0;
	printf("queryTrajNum:%d", queryTrajNum);
	printf("totalQueryTrajLength:%d", totalQueryTrajLength);
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
	for (int i = 0; i <= qTra[qID].length - 1; i++)
	{
	queryTra[i].x = qTra[qID].points[i].lon;
	queryTra[i].y = qTra[qID].points[i].lat;
	queryTra[i].tID = qTra[qID].tid;
	}
	CUDA_CALL(cudaMemcpyAsync(queryTraGPU, queryTra, sizeof(SPoint)*qTra[qID].length, cudaMemcpyHostToDevice, defaultStream));
	queryTraLength[qID] = qTra[qID].length;
	queryTraGPU = queryTraGPU + qTra[qID].length;
	queryTra += qTra[qID].length;
	if (qID != queryTrajNum - 1)
	allQueryTraOffset[qID + 1] = allQueryTraOffset[qID] + qTra[qID].length;
	}
	nowAddrGPU = queryTraGPU;
	// queryTraOffsetGPU是保存queryTra中offset的基地址
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);

	//构建queryLength
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	//----------------------------------------------------
	//单GPU版本备份
	*/


	//tt.stop();
	//cout << "Part3.0.1 time:" << tt.elapse() << endl;
	//tt.start();
	//第一步：循环，剪枝
	int* worstNow = new int[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		worstNow[qID] = 9999999;
	}

	bool* isFinished = new bool[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		isFinished[qID] = FALSE;
	}
	bool isAllFinished = FALSE;


	SPoint* candidateTra[2];
	SPoint* candidateTraGPU[2];

	//map<int, void*> traID_baseAddr[2];

	int candidateTrajNum[2] = { 0 };
	int* candidateTraLength[2];
	TaskInfoTableForSimilarity* taskInfoTable[2];
	OffsetTable* candidateTrajOffsetTable[2];
	SPoint** candidateOffsets[2];
	int** candidateTrajID = new int*[queryTrajNum];

	for (int i = 0; i <= queryTrajNum - 1; i++)
		candidateTrajID[i] = new int[k];

	for (int device_idx = 0; device_idx <= 1; device_idx++)
	{
		CUDA_CALL(cudaSetDevice(device_idx));// no need!!
		candidateTra[device_idx] = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);
		taskInfoTable[device_idx] = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
		candidateTrajOffsetTable[device_idx] = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum);
		candidateOffsets[device_idx] = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum);
		candidateTraGPU[device_idx] = (SPoint*)nowAddrGPU[device_idx];
		candidateTraLength[device_idx] = (int*)malloc(sizeof(int) * k * queryTrajNum);
	}

	/*
	//单GPU版本备份

	SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);

	//保存qid、candID和在candidateTran中的offset的对应关系
	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
	//轨迹唯一，保存在id下轨迹的id和baseAddr （有必要保存轨迹的id吗？）
	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum);
	//保存candidateOffset的基地址们
	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum);
	int candidateTrajNum = 0;
	// traID和在candidateTrajOffsetTable中的idx的对应关系map，主要用于判断轨迹是否已经复制到gpu
	map<int, void*> traID_baseAddr;
	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU;

	//单GPU版本备份
	*/


	//tt.stop();
	//cout << "Part3.0.2 time:" << tt.elapse() << endl;

	while (!isAllFinished)
	{
		//tt.start();
		//还有待计算的候选轨迹的总数目，相当于任务数目
		int validCandTrajNum[2] = { 0 };


		// 目前还没有计算完成的查询轨迹的数目
		int validQueryTraNum = queryTrajNum;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
		}


		// 在本轮计算内的第几个轨迹
		int validQueryIdx = 0;
		int validQueryIdxGPU[2] = { 0 };
		int queryEachGPU[2];



		// 多GPU加速的关键 每次大while 40*40(m) -> 20*40 20*40 上面的任务个数比之前少一半 （固定k）所以快
		// 任务不可控 数据依赖性太大
		// 剩40 20+20
		// 剩30 15+15 Tq长度完全不可控

		queryEachGPU[0] = validQueryTraNum / 2;
		queryEachGPU[1] = validQueryTraNum - queryEachGPU[0];



		// 线性地址
		SPoint* tempPtr[2];
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			tempPtr[device_idx] = candidateTra[device_idx];

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			// 选择 GPU
			int device_idx = 0;
			if (validQueryIdx < queryEachGPU[0])
				device_idx = 0;
			else
				device_idx = 1;

			// 未完成则拷贝Tc
			if (!isFinished[qID])
			{

				//提取topk作为待查轨迹
				for (int i = 0; i <= k - 1; i++)
				{
					candidateTrajID[qID][i] = queryQueue[qID].top().traID;
					queryQueue[qID].pop();
					validCandTrajNum[device_idx]++;
				}

				for (int i = 0; i <= k - 1; i++)
				{
					int CandTrajID = candidateTrajID[qID][i];
					//map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr[device_idx].find(CandTrajID);
					// 如果轨迹还没有保存在GPU中
					//if (traID_baseAddr_ITER == traID_baseAddr[device_idx].end())
					{
						int pointsNumInThisCand = 0;
						SPoint* thisTrajAddr = tempPtr[device_idx];
						for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++)
						{
							int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
							memcpy(tempPtr[device_idx], &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
							tempPtr[device_idx] += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
							pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
						}
						// 把这条轨迹提取到candidateTraGPU中
						CUDA_CALL(cudaMemcpyAsync(candidateTraGPU[device_idx], thisTrajAddr, pointsNumInThisCand * sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream[device_idx]));
						
						//traID_baseAddr[device_idx][candidateTrajID[qID][i]] = candidateTraGPU[device_idx];
						
						//仅保存要计算的query的candidateTraLength
						candidateTraLength[device_idx][k * validQueryIdxGPU[device_idx] + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
						//仅保存要计算的query的offset
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].qID = qID;
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].candTrajID = CandTrajID;
						// 保存轨迹对应的addr
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].objectId = candidateTrajID[qID][i];
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].addr = candidateTraGPU[device_idx];
						candidateOffsets[device_idx][k * validQueryIdxGPU[device_idx] + i] = candidateTraGPU[device_idx];
						//地址往前移动，便于下一次复制
						candidateTraGPU[device_idx] = (candidateTraGPU[device_idx] + pointsNumInThisCand);
						// nowAddrGPU 始终是指向下一个空闲的GPU地址
						nowAddrGPU[device_idx] = (void*)candidateTraGPU[device_idx];
					}
					// 如果该轨迹已经复制进了gpu里面，那么只需要按照该轨迹id更新表就行了

					//else
					//{
					//	void* baseAddrGPU = traID_baseAddr_ITER->second;
					//	//仅保存要计算的query的candidateTraLength
					//	candidateTraLength[device_idx][k * validQueryIdxGPU[device_idx] + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
					//	//仅保存要计算的query的offset
					//	taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].qID = qID;
					//	taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].candTrajID = CandTrajID;
					//	// 保存轨迹对应的addr
					//	candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].objectId = CandTrajID;
					//	candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].addr = baseAddrGPU;
					//	candidateOffsets[device_idx][k * validQueryIdxGPU[device_idx] + i] = (SPoint*)baseAddrGPU;
					//}
				}
				validQueryIdx++;
				validQueryIdxGPU[device_idx]++;
			}
		}


		int* candidateTraLengthGPU[2];
		TaskInfoTableForSimilarity* taskInfoTableGPU[2];
		SPoint** candidateOffsetsGPU[2];
		int* resultReturnedGPU[2];
		int* resultReturned = new int[queryTrajNum * k];
		for (int device_idx = 0; device_idx <= 1; device_idx++) {
			candidateTraLengthGPU[device_idx] = (int*)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(candidateTraLengthGPU[device_idx], candidateTraLength[device_idx], sizeof(int)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU 始终是指向下一个空闲的GPU地址
			whileAddrGPU[device_idx] = (void*)((int*)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//构建TaskInfoTable
			taskInfoTableGPU[device_idx] = (TaskInfoTableForSimilarity*)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU[device_idx], taskInfoTable[device_idx], sizeof(TaskInfoTableForSimilarity)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU 始终是指向下一个空闲的GPU地址
			whileAddrGPU[device_idx] = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//复制candidate的地址到gpu中
			candidateOffsetsGPU[device_idx] = (SPoint**)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU[device_idx], candidateOffsets[device_idx], sizeof(SPoint*)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU 始终是指向下一个空闲的GPU地址
			whileAddrGPU[device_idx] = (void*)((SPoint**)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//构建candidateTraj和candidateLength完成，准备并行Similarity search
			//只需要查询isFinished为false的queryTra，至于是哪些可以直接看offsetTableCandidateTra

			resultReturnedGPU[device_idx] = (int*)whileAddrGPU[device_idx];
			whileAddrGPU[device_idx] = (void*)((int*)whileAddrGPU[device_idx] + k * queryTrajNum);
		}
		//CUDA_CALL(cudaMalloc((void**)resultReturnedGPU, sizeof(int)*k*queryTrajNum));


		//tt.stop();
		//cout << "Part3.2 time:" << tt.elapse() << endl;
		//tt.start();
		//如果上面的分配没有错，开始计算EDR
		if (validQueryTraNum * k == (validCandTrajNum[0] + validCandTrajNum[1]))
		{
			// GPU 1 2 工作顺序 由于不同的工作流 stream0 stream1 保证多GPU 同时运行
			// 利用 stream 进行 2 GPU并行
			for (int device_idx = 0; device_idx <= 1; device_idx++) {// for 循环进行控制 先GPU-1 再 GPU-2 不过很快 没什么影响！！！
				CUDA_CALL(cudaSetDevice(device_idx));
				if (validCandTrajNum[device_idx] == 0) // 边界条件
					continue;
				EDRDistance_Batch_Handler(validCandTrajNum[device_idx], taskInfoTableGPU[device_idx], queryTraGPUBase[device_idx], queryTraOffsetGPU[device_idx], candidateOffsetsGPU[device_idx], queryLengthGPU[device_idx], candidateTraLengthGPU[device_idx], resultReturnedGPU[device_idx], &defaultStream[device_idx]);
			}
			CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU[0], sizeof(int)*k*queryEachGPU[0], cudaMemcpyDeviceToHost, defaultStream[0])); // 注意这里有stream[0]的标志
			CUDA_CALL(cudaMemcpyAsync(resultReturned + k*queryEachGPU[0], resultReturnedGPU[1], sizeof(int)*k*queryEachGPU[1], cudaMemcpyDeviceToHost, defaultStream[1])); // 注意这里有stream[1]的标志

																																										   // 这里存在GPU1 GPU2 的隐式同步？？应该是

		}
		else
		{
			printf("error in line 1007\n");
		}

		//tt.stop();
		//cout << "Part3.3 time:" << tt.elapse() << endl;
		//tt.start();

		for (int idx = 0; idx <= k * validQueryTraNum - 1; idx++)
		{
			int device_idx = 0;
			int idxInTaskInfoTable = 0;
			if (idx / k < queryEachGPU[0])
			{
				device_idx = 0;
				idxInTaskInfoTable = idx;
			}
			else
			{
				device_idx = 1;
				idxInTaskInfoTable = idx - k*queryEachGPU[0];
			}
			int qID = taskInfoTable[device_idx][idxInTaskInfoTable].qID;

			int i = idx % k;
			if (numElemInCalculatedQueue[qID] < kValue)
			{
				//直接往PQ里加
				FDwithID fd;
				fd.traID = candidateTrajID[qID][i];
				fd.FD = resultReturned[idx];
				EDRCalculated[qID].push(fd);
				numElemInCalculatedQueue[qID]++;
			}
			else
			{
				//看一下是否比PQ里更好，如果是弹出一个差的，换进去一个好的；否则不动优先队列也不更新worstNow。
				int worstInPQ = EDRCalculated[qID].top().FD;
				if (resultReturned[i] < worstInPQ)
				{
					EDRCalculated[qID].pop();
					FDwithID fd;
					fd.traID = candidateTrajID[qID][i];
					fd.FD = resultReturned[idx];
					EDRCalculated[qID].push(fd);
				}
			}
		}

		// 在这里判断 比较清晰 可读性
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			worstNow[qID] = EDRCalculated[qID].top().FD;
			if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD))
			{
				isFinished[qID] = TRUE;
			}
		}

		bool temp = TRUE;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			temp = temp && isFinished[qID];
		}
		isAllFinished = temp;
		delete[] resultReturned;

		//GPU指针回到while开始的地方
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			whileAddrGPU[device_idx] = whileAddrGPUBase[device_idx];

		//tt.stop();
		//cout << "Part3.4 time:" << tt.elapse() << endl;
	} // end of BIG while


	timer.stop();
	cout << "Part3 EDR time:" << timer.elapse() << endl;

	//输出结果
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int i = 0; i <= kValue - 1; i++)
		{
			topKSimilarityTraj[qID * kValue + i] = EDRCalculated[qID].top().traID;
			EDRCalculated[qID].pop();
		}
	}

	for (int i = 0; i <= queryTrajNum - 1; i++)
		delete[] candidateTrajID[i];
	for (int device_idx = 0; device_idx <= 1; device_idx++) {
		free(taskInfoTable[device_idx]);
		free(candidateTrajOffsetTable[device_idx]);
		free(candidateOffsets[device_idx]);
		free(candidateTra[device_idx]);
		free(candidateTraLength[device_idx]);
		cudaStreamDestroy(defaultStream[device_idx]);
		CUDA_CALL(cudaFree(baseAddrSimi[device_idx]));
		CUDA_CALL(cudaFree(whileAddrGPU[device_idx]));
	}

	delete[] candidateTrajID;
	delete[] isFinished;
	delete[] worstNow;
	free(allQueryTra);
	delete[] allQueryTraOffset;
	delete[] EDRCalculated;
	delete[] numElemInCalculatedQueue;
	delete[] freqVectors;
	delete[] queryQueue;
	delete[] queryTraLength;

	totaltimer.stop();
	cout << "Fined Dual GPU Time:" << totaltimer.elapse() << endl;

	return 0;
}



// CoarseGrainedThreadJoin GAT-S-noE + GAT-S-E 成了串行？ 并没有！
int Grid::SimilarityQueryBatchOnMultiGPUV2(Trajectory* qTra, int queryTrajNum, int* topKSimilarityTraj, int kValue) {

	MyTimer timer;
	int device_num = 2;
	
	int simNumGPU[2];
	// 40*40(m) -> 20*40 + 20*40
	simNumGPU[0] = queryTrajNum / 2;
	simNumGPU[1] = queryTrajNum - simNumGPU[0];
	int startIdx[2];
	startIdx[0] = 0;
	startIdx[1] = simNumGPU[0];
	//void* allocatedGPUMem[2] = { NULL };



	int device_idx = 0;

	cout << "noE\n";
	// noE or E, decide here

	vector<thread> threads_RQ;
	vector<int> stateTableSim[2];
	stateTableSim[0].resize(simNumGPU[0] * kValue);
	stateTableSim[1].resize(simNumGPU[1] * kValue); // 这里真正的结果

	for (device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		//CUDA_CALL(cudaSetDevice(device_idx));
		//CUDA_CALL(cudaMalloc((void**)&this->baseAddrRange[device_idx], (long long int)BIG_MEM * 1024 * 1024));
		//CUDA_CALL(cudaMalloc((void**)&this->stateTableGPU[device_idx], (long long int)SMALL_MEM * 1024 * 1024));
		//allocatedGPUMem[device_idx] = this->baseAddrRange[device_idx];
		threads_RQ.push_back(thread(std::mem_fn(&Grid::SimilarityQueryBatchOnGPUV3), this, &qTra[startIdx[device_idx
		]], simNumGPU[device_idx], &stateTableSim[device_idx][0], kValue, device_idx));
	}
	timer.start();
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));

	// 这里由于thread 会产生默认的同步操作 所以还是balanced比较好

	timer.stop();
	cout << "Coarse Dual GPU Time:" << timer.elapse() << endl; // 包含 FV FD EDR 只看EDR 绝对是减少一半运行时间！！



	cout << "noMAT noE\n";
	// noE or E, decide here

	vector<thread> threads_RQ2;
	vector<int> stateTableSim2[2];
	stateTableSim2[0].resize(simNumGPU[0] * kValue);
	stateTableSim2[1].resize(simNumGPU[1] * kValue); // 这里真正的结果

	for (device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		//CUDA_CALL(cudaSetDevice(device_idx));
		//CUDA_CALL(cudaMalloc((void**)&this->baseAddrRange[device_idx], (long long int)BIG_MEM * 1024 * 1024));
		//CUDA_CALL(cudaMalloc((void**)&this->stateTableGPU[device_idx], (long long int)SMALL_MEM * 1024 * 1024));
		//allocatedGPUMem[device_idx] = this->baseAddrRange[device_idx];
		threads_RQ2.push_back(thread(std::mem_fn(&Grid::SimilarityQueryBatchOnGPUNoMAT), this, &qTra[startIdx[device_idx
		]], simNumGPU[device_idx], &stateTableSim2[device_idx][0], kValue, device_idx));
	}
	timer.start();
	std::for_each(threads_RQ2.begin(), threads_RQ2.end(), std::mem_fn(&std::thread::join));

	// 这里由于thread 会产生默认的同步操作 所以还是balanced比较好
	timer.stop();
	cout << "Coarse Dual GPU Time:" << timer.elapse() << endl; // 包含 FV FD EDR 只看EDR 绝对是减少一半运行时间！！



	//for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	//{
	//	CUDA_CALL(cudaFree(allocatedGPUMem[device_idx]));
	//	CUDA_CALL(cudaFree(this->stateTableGPU[device_idx]));
	//}
	return 0;


}



//int Grid::SimilarityQuery(Trajectory & qTra, Trajectory **candTra, const int candSize, float * EDRdistance)
//{
//	cout << candSize << endl;
//	SPoint *queryTra = (SPoint*)malloc(sizeof(SPoint)*(qTra.length));
//	for (int i = 0; i <= qTra.length - 1; i++) {
//		queryTra[i].x = qTra.points[i].lon;
//		queryTra[i].y = qTra.points[i].lat;
//		queryTra[i].tID = qTra.points[i].tid;
//	}
//
//	SPoint **candidateTra = (SPoint**)malloc(sizeof(SPoint*)*candSize);
//
//	for (int i = 0; i <= candSize - 1; i++) {
//		candidateTra[i] = (SPoint*)malloc(sizeof(SPoint)*(candTra[i]->length)); //调试的时候这一部分总报内存错误，FFFFF
//		for (int j = 0; j <= candTra[i]->length - 1; j++) {
//			candidateTra[i][j].x = candTra[i]->points[j].lon;
//			candidateTra[i][j].y = candTra[i]->points[j].lat;
//			candidateTra[i][j].tID = candTra[i]->points[j].tid;
//		}
//	}
//
//	int queryLength = qTra.length;
//	int *candidateLength = (int*)malloc(sizeof(int)*candSize);
//	for (int i = 0; i <= candSize - 1; i++) {
//		candidateLength[i] = candTra[i]->length;
//	}
//
//	int* result = (int*)malloc(sizeof(int)*candSize);
//
//	MyTimer timer1;
//	timer1.start();
//
//	//CPU
//	int *resultCPU = (int*)malloc(sizeof(int)*candSize);
//	for (int i = 0; i <= candSize - 1; i++) {
//		//每个DP问题
//		SPoint *CPUqueryTra = queryTra, *CPUCandTra = candidateTra[i];
//		int CPUqueryLength = qTra.length, CPUCandLength = candidateLength[i];
//		int longest = 0;
//
//		const SPoint *tra1, *tra2;
//		int len1, len2;
//		if (CPUCandLength >= CPUqueryLength) {
//			tra1 = CPUqueryTra;
//			tra2 = CPUCandTra;
//			len1 = CPUqueryLength;
//			len2 = CPUCandLength;
//		}
//		else
//		{
//			tra1 = CPUCandTra;
//			tra2 = CPUqueryTra;
//			len1 = CPUCandLength;
//			len2 = CPUqueryLength;
//		}
//
//		if (CPUqueryLength >= longest) {
//			longest = CPUqueryLength;
//		}
//		else
//		{
//			longest = CPUCandLength;
//		}
//
//
//		int **stateTable = (int**)malloc(sizeof(int*)*(len1 + 1));
//		for (int j = 0; j <= len1; j++) {
//			stateTable[j] = (int*)malloc(sizeof(int)*(len2 + 1));
//		}
//		stateTable[0][0] = 0;
//		for (int row = 1; row <= len1; row++) {
//			stateTable[row][0] = row;
//		}
//		for (int col = 1; col <= len2; col++) {
//			stateTable[0][col] = col;
//		}
//
//		for (int row = 1; row <= len1; row++) {
//			for (int col = 1; col <= len2; col++) {
//				SPoint p1 = tra1[row - 1];
//				SPoint p2 = tra2[col - 1]; //这样做内存是聚集访问的吗？
//				bool subcost;
//				if((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON)) {
//					subcost = 0;
//				}
//				else
//					subcost = 1;
//				int myState = 0;
//				int state_ismatch = stateTable[row - 1][col - 1] + subcost;
//				int state_up = stateTable[row - 1][col] + 1;
//				int state_left = stateTable[row][col - 1] + 1;
//				if (state_ismatch < state_up)
//					myState = state_ismatch;
//				else if (state_left < state_up)
//					myState = state_left;
//				else
//					myState = state_ismatch;
//
//				stateTable[row][col] = myState;
//				//	if (row == len1&&col == len2)
//						//cout << myState << endl;
//			}
//		}
//
//		resultCPU[i] = stateTable[len1][len2];
//		//cout << resultCPU[i] << endl;
//	}
//	timer1.stop();
//	cout << "CPU Similarity Time:" << timer1.elapse() << "ms" << endl;
//	//GPU
//
//	timer1.start();
//	handleEDRdistance(queryTra, candidateTra, candSize, queryLength, candidateLength, result);
//	timer1.stop();
//	cout << "GPU Similarity Time:" << timer1.elapse() << "ms" << endl;
//
//	for (int i = 0; i <= candSize - 1; i++) {
//		EDRdistance[i] = result[i];
//	}
//	free(queryTra);
//	for (int i = 0; i <= candSize - 1; i++) {
//		free(candidateTra[i]);
//	}
//	free(candidateTra);
//	free(candidateLength);
//	free(result);
//
//	return 0;
//}



