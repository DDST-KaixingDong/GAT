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
		return -1; //̫С����
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
	for (int i = 1; i <= 2 * max(lenx, leny); i++) // ��ʶ�����ƶ�����λ
	{
		if ((i & 1) == 1) // ������ĩβΪ1 ��ȡ�к�
		{
			result += (x >> (xbit - 1) & 1) * (1 << (i - 1)); // x >> (xbit - 1) & 1 ���� ����λ���� ��λ����
			xbit = xbit + 1;
		}
		else // ������ĩβΪ0  ��ȡ�к�
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

//ȷ������
int Grid::WhichCellPointIn(SamplePoint p)
{
	//ע��cell����Ǵ�(xmin,ymax)��ʼ�ģ�������(xmin,ymin)
	int row = (int)((range.ymax - p.lat) / cell_size); //��0��ʼ
	int col = (int)((p.lon - range.xmin) / cell_size); //��0��ʼ
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

//�ѹ켣t������ӹ켣����ӵ�cell����
//��һ�������ǰ��ӹ켣�Ž���cell���棬���һ����item
int Grid::addTrajectoryIntoCell(Trajectory& t)
{
	if (t.length == 0)
		return 1;//�չ켣
	SamplePoint p = t.points[0];
	int lastCellNo = WhichCellPointIn(p);
	int lastCellStartIdx = 0;
	int nowCellNo;
	//cell based traj���ɣ��ǵ�ת����free��
	vector<int>* tempCellBasedTraj = new vector<int>;
	tempCellBasedTraj->reserve(1048577);
	int tempCellNum = 0;
	for (int i = 0; i <= t.length - 1; i++)
	{
		p = t.points[i];
		nowCellNo = WhichCellPointIn(p);

		if (i == t.length - 1)
		{
			//�����һ�����������cellҲ���ϸ�cell�������һ��cell�ˣ����֮
			if (lastCellNo == nowCellNo)
			{
				tempCellNum++;
				tempCellBasedTraj->push_back(nowCellNo);
				cellPtr[nowCellNo].addSubTra(t.tid, lastCellStartIdx, i, i - lastCellStartIdx + 1);
				int vituralCellNo = nowCellNo >> VITURAL_CELL_PARAM; //�����
				this->freqVectors.addPointToFVTable(t.tid, i - lastCellStartIdx + 1, vituralCellNo);
			}
			//������һ�������cell��Ҫ���
			else
			{
				tempCellNum += 2;
				tempCellBasedTraj->push_back(lastCellNo);
				tempCellBasedTraj->push_back(nowCellNo);
				cellPtr[lastCellNo].addSubTra(t.tid, lastCellStartIdx, i - 1, i - 1 - lastCellStartIdx + 1);
				int vituralCellNo = lastCellNo >> VITURAL_CELL_PARAM; //�����
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
				// �ս�һ���ӹ켣����ʼ��һ���ӹ켣
				//cellTra�����һ��
				tempCellNum++;
				tempCellBasedTraj->push_back(lastCellNo);
				//SubTra���
				//printf("cell:%d\n", lastCellNo);
				cellPtr[lastCellNo].addSubTra(t.tid, lastCellStartIdx, i - 1, i - 1 - lastCellStartIdx + 1);
				int vituralCellNo = lastCellNo >> VITURAL_CELL_PARAM; //�����
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
	//����frequency vector
	this->freqVectors.initFVTable(traNum, (this->cellnum)>> this->VITURAL_CELL_PARAM); //ע������������ӵĸ���
	//ע�⣬�켣��Ŵ�1��ʼ
	this->cellBasedTrajectory.resize(traNum + 1); //����cellbasedtraj�Ĺ�ģ���ӹ켣��ʱ�����ֱ����
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
			//choose j; //ѡ���j��subTra�����������滹��û����ʼ�����Ϣ�����ǻ�����point����Ϣ�ӽ�subtraTable�������������Ҫ��subTraTable��free��
			this->cellBasedTrajectory[i].startIdx[cellidx] = this->cellPtr[nowCellID].subTraTable[j].idxInAllPointsArray;
			this->cellBasedTrajectory[i].numOfPointInCell[cellidx] = this->cellPtr[nowCellID].subTraTable[j].numOfPoint;
			tempCntForTraj[nowCellID]++;
		}
		free(tempCntForTraj);
	}
	// Transfer FV to GPU
	//this->freqVectors.transferFVtoGPU();


	////Delta Encoding��cell�����洢
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

	//�����ɺõ�allpoints�ŵ�GPU��
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
	const int device_num = 16; // ��������16������ ˵��MCPU ֪���� ��20������Ū���� ����CPU ������ô
	int nowidx = 0;
	// ����batchƽ��Ϊ16��
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
		threads_RQ.push_back(thread(std::mem_fn(&Grid::rangeQueryBatch), this, &bounds[qID], 1, &ResultTable[qID], resultSetSize));// resultSetSize û�õ�
		//threads_RQ.push_back(thread(std::mem_fn(&Grid::SimilarityMultiThreadHandler), this, queryQueue, qTra, 1, EDRCalculated, kValue, qID));
		//threads_EDR.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}

	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join)); // threads_FD�Ѿ�ָ��

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

	//������ɣ���ʼ����kernel��ѯ
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
	//��ѯ�������ƺ����stateTable�����gpu��
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
	//stateTable�е����Ŀ�����ֵ
	// FILE* testFile = fopen("freq.txt", "w+");
	int maxPointNum = 0;
	for (int i = 0; i <= this->stateTableLength[device_idx] - 1; i++)
	{
		if (stateTableAllocate[i].candidatePointNum > maxPointNum)
			maxPointNum = stateTableAllocate[i].candidatePointNum;
		//fprintf(testFile, "%d ", stateTableAllocate[i].candidatePointNum);
	}
	//����GPU���в��в�ѯ
	//�ȴ���stateTable
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
	//������ɣ���ʼ����kernel��ѯ
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
	//��ѯ�������ƺ����stateTable�����gpu��
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
		//���gpu�ڴ���û�и�node����Ϣ
		if (this->nodeAddrTable[device_idx].find(startCellID) == this->nodeAddrTable[device_idx].end())
		{
			CUDA_CALL(cudaMemcpyAsync(this->baseAddrRange[device_idx], &(this->allPoints[startIdx]), pointNum*sizeof(SPoint), cudaMemcpyHostToDevice, stream));
			dataPtr = (SPoint*)this->baseAddrRange[device_idx];
			this->nodeAddrTable[device_idx].insert(pair<int, void*>(startCellID, this->baseAddrRange[device_idx]));
			this->baseAddrRange[device_idx] = (void*)((char*)this->baseAddrRange[device_idx] + pointNum * sizeof(SPoint));
		}
		//����У����ٸ��ƣ�ֱ����
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


		//���gpu�ڴ���û�и�node����Ϣ
		if (this->nodeAddrTable[device_idx].find(startCellID) == this->nodeAddrTable[device_idx].end())
		{
			CUDA_CALL(cudaMemcpyAsync(this->baseAddrRange[device_idx], &(this->allPoints[startIdx]), pointNum * sizeof(SPoint), cudaMemcpyHostToDevice, stream));
			dataPtr = (SPoint*)this->baseAddrRange[device_idx];
			this->nodeAddrTable[device_idx].insert(pair<int, void*>(startCellID, this->baseAddrRange[device_idx]));
			this->baseAddrRange[device_idx] = (void*)((char*)this->baseAddrRange[device_idx] + pointNum * sizeof(SPoint));
		}
		//����У����ٸ��ƣ�ֱ����
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
	this->freqVectors.formPriorityQueue(queue, freqVectorQ);// �ؼ����� std::vector<std::map<int, int>> FreqVector; //ÿ��map����һ���켣 trajDB�еĹ켣
	return 0;
}

// �����ҵ�top-k��EDR
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

		while (worstNow > queryQueue[qID].top().FD)// worst now �������EDR
		{
			/*tt.start();*/

			int candidateTrajID[k];// m ����ID
								   //printf("%d", worstNow);

								   // ��0������ȡtopk(topm)FD
			for (int i = 0; i <= k - 1; i++)
			{
				candidateTrajID[i] = queryQueue[qID].top().traID;// FDС���ѵ�ID����pop��candidateTrajID[]��													
				queryQueue[qID].pop();
			}
			//EDR calculate

			//��һ������AllPoints����ȡ�����켣(m��)
			SPoint** candidateTra = (SPoint**)malloc(sizeof(SPoint*) * k);
			int* candidateTraLength = (int*)malloc(sizeof(int) * k);

			for (int i = 0; i <= k - 1; i++)
			{
				candidateTra[i] = (SPoint*)malloc(sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[i]].trajLength);
				SPoint* tempPtr = candidateTra[i];
				// ����ÿ���켣 ����Lt
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

			//�ڶ���������EDR
			//printf("%d", qID);
			int resultReturned[k];
			this->SimilarityExecuter(queryTra, candidateTra, qTra[qID].length, candidateTraLength, k, resultReturned);
			//tt.stop();
			//cout << "Part3.3 time:" << tt.elapse() << endl;
			for (int i = 0; i <= k - 1; i++)
			{
				if (numElemInCalculatedQueue[qID - startQueryIdx] < kValue)//  qID - startQueryIdxm = 0 �Ƚ�ͷԪ�� ��һ�μ���
				{
					//ֱ����PQ���
					FDwithID fd;
					fd.traID = candidateTrajID[i];
					fd.FD = resultReturned[i];
					EDRCalculated[qID].push(fd);
					numElemInCalculatedQueue[qID - startQueryIdx]++;
				}
				else // �ǵ�һ�μ���
				{
					//��һ���Ƿ��PQ����ã�����ǵ���һ����ģ�����ȥһ���õģ����򲻶����ȶ���Ҳ������worstNow��
					int worstInPQ = EDRCalculated[qID].top().FD;
					if (resultReturned[i] < worstInPQ)
					{
						EDRCalculated[qID].pop();
						FDwithID fd;
						fd.traID = candidateTrajID[i];
						fd.FD = resultReturned[i];
						EDRCalculated[qID].push(fd);// ע�� ����fd��һ��������ͷ�� ���ȶ��л��Զ�����
					}
				}
			}
			worstNow = EDRCalculated[qID].top().FD;
			//printf("%d,worstNow:%d\t", qID,worstNow);
			//���ֽ������ͷ��ڴ�
			for (int i = 0; i <= k - 1; i++)
				free(candidateTra[i]);
			free(candidateTraLength);
			free(candidateTra); //ע��candidateTra��ά�����free
		}
		//timer.stop();
		//cout << "Query Trajectory Length:" << qTra[qID].length << endl;
		//cout << "Part3 time:" << timer.elapse() << endl;
		//timer.start();
		free(queryTra); // free����ʱ����

						//timer.stop();
						//cout << "Part4 time:" << timer.elapse() << endl;
	}
	return 0;
}

int Grid::SimilarityExecuter(SPoint* queryTra, SPoint** candidateTra, int queryLength, int* candidateLength, int candSize, int* resultArray)
{
	for (int i = 0; i <= candSize - 1; i++) // 0 - m-1
	{
		// ����m���е�ÿһ����Qt��EDR
		//ÿ��DP����
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
				SPoint p2 = tra2[col - 1]; //�������ڴ��Ǿۼ����ʵ���
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
				//ȥ��if�ı�﷽ʽ���Ƿ�����������ܣ� ��ȫ�Ǵ�GPU�Ż��Ƕ�
				myState = c1 * state_ismatch + c2 * state_left + !(c1 || c2) * state_up;

				stateTable[row][col] = myState;
				//	if (row == len1&&col == len2)
				//cout << myState << endl;
			}
		}


		resultArray[i] = stateTable[len1][len2];
		//cout << resultCPU[i] << endl;

		// ��һ�ζ�άnew�����delete
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
	//˼·���ò�ͬ�̴߳���ͬ������query
	// С��������С��FD ��Ĭ�� ע����priority_queue��ָ�� ���ƶ�ά
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	// ע����map��ָ�� ���ƶ�ά ��Ϊpriority_queue��map������������
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum];
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++) // ����qTra�ڵ�
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
	// ��CPU�߳�
	// windows ��Ҫ�� ̫��
	//cout << "Part2 time: start �������߳� ��\n";
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join)); // threads_FD�Ѿ�ָ��

#else
	// ��CPU�߳�
	// ����FD ����queryQueue������ ��ÿ��queryQueue��priority_queue ���ȶ��� with CPU ���߳� ���߳��ܲ��� �����ǱʼǱ�ԭ��
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
	// way1: ����range ������batch ƽ��Ϊ 16 ��
	

	// way2: ��CPU�Լ����� ����batch��thread���� ʵ����ô���أ�  ע����� 1 qID
	
	vector<thread> threads_EDR;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		threads_EDR.push_back(thread(std::mem_fn(&Grid::SimilarityMultiThreadHandler), this, queryQueue, qTra, 1, EDRCalculated, kValue, qID));
		//threads_EDR.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_EDR.begin(), threads_EDR.end(), std::mem_fn(&std::thread::join)); // threads_FD�Ѿ�ָ��


	/*
	// 20 threads !!
	const int THREAD_CPU = 20; 
	vector<thread> threads;
	int qIDnow = 0;
	// Tq ���У�ÿ 20 ��Tq����EDR ����
	// recall that on GPU, Tq*m(40)
	for (int i = qIDnow; qIDnow <= queryTrajNum - 1; i++)// i û��
	{
		for (int j = 0; j <= THREAD_CPU - 1; j++)
		{
			if (qIDnow>queryTrajNum - 1)
				break;
			threads.push_back(thread(std::mem_fn(&Grid::SimilarityMultiThreadHandler), this, queryQueue, qTra, 1, EDRCalculated, kValue, qIDnow));
			qIDnow++;
		}
		std::for_each(threads.begin(), threads.end(), std::mem_fn(&std::thread::join));
		threads.clear();// ���ص������½� �����޷��ﵽ 2000% ռ����CPU
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
	//˼·���ֱ���ÿһ����ѯ�켣���ò�ͬstream����
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum];
	//Ϊ��ѯ����freqVector
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //�����
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
	//Ϊ��֦����Frequency Distance
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Part2 time:" << timer.elapse() << endl;
	//��һ�����ȶ��д洢��ǰ���Ž�����󶥶ѣ���֤��ʱ����pop����Ľ��
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];
	int* numElemInCalculatedQueue = new int[queryTrajNum]; //���浱ǰ���ȶ��н������֤���ȶ��д�С������kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//׼����֮�󣬿�ʼ����ѯ
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
			//��ȡtopk
			for (int i = 0; i <= k - 1; i++)
			{
				candidateTrajID[i] = queryQueue[qID].top().traID;
				//printf("%d,%d\t", queryQueue[qID].top().traID,queryQueue[qID].top().FD);
				queryQueue[qID].pop();
			}
			//EDR calculate
			//��һ������AllPoints����ȡ�����켣
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
			//�ڶ���������EDR
			//printf("%d", qID);
			int resultReturned[k];
			this->SimilarityExecuter(queryTra, candidateTra, qTra[qID].length, candidateTraLength, k, resultReturned);
			//tt.stop();
			//cout << "Part3.3 time:" << tt.elapse() << endl;
			//����worstNow
			for (int i = 0; i <= k - 1; i++)
			{
				if (numElemInCalculatedQueue[qID] < kValue)
				{
					//ֱ����PQ���
					FDwithID fd;
					fd.traID = candidateTrajID[i];
					fd.FD = resultReturned[i];
					EDRCalculated[qID].push(fd);
					numElemInCalculatedQueue[qID]++;
				}
				else
				{
					//��һ���Ƿ��PQ����ã�����ǵ���һ����ģ�����ȥһ���õģ����򲻶����ȶ���Ҳ������worstNow��
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
			//���ֽ������ͷ��ڴ�
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
//���У�������query�����ϵ�GPU����
//����˼·���ֱ���ÿһ����ѯ�켣���ò�ͬstream����
{
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU), (long long int)BIG_MEM * 1024 * 1024));
	void* whileAddrGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024));
	void* whileAddrGPUBase = whileAddrGPU;
	
	//��ǰ���䵽�ĵ�ַ
	void* nowAddrGPU = NULL;    
	 
	cudaStream_t defaultStream;
	cudaStreamCreate(&defaultStream);

	MyTimer timer;

	// С��������С��FD ��Ĭ�� ע����priority_queue��ָ�� ���ƶ�ά
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	// ע����map��ָ�� ���ƶ�ά ��Ϊpriority_queue��map������������
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum]; // FVtable

	//Ϊ��ѯ����freqVector ���壺����ѯ�켣��freqVector��freqVectors[qID]  = CPU
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //�����
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
	//Ϊ��֦����Frequency Distance
	//CPU���� parallelly �ǵ�
	vector<thread> threads_FD;
	// ����FD ����queryQueue������ ��ÿ��queryQueue��priority_queue ���ȶ��� = CPU
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
	
	int* numElemInCalculatedQueue = new int[queryTrajNum]; //���浱ǰ���ȶ��н������֤���ȶ��д�С������kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;


	//׼����֮�󣬿�ʼ����ѯ

	const int k = KSIMILARITY; //m
	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}
	//��ѯ�켣��Ϣ��׼����

	// Step0��*********** ���� Tq ��ѯ�Ĺ켣ȫ����Ϣ CPU->GPU *************
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


	// �����EDR
	int* worstNow = new int[queryTrajNum];
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		worstNow[qID] = 9999999;
	}
	// �ñ�־λ ÿ��Tq������batch����һ����־λ
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


		// Step1��*********** CPU->GPU  Tc ��صĿ��� MAT****************
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
			if (!isFinished[qID])
			{
				if ((queryQueue[qID].empty()) || (worstNow[qID] <= queryQueue[qID].top().FD))
				{
					validQueryTraNum--;
					isFinished[qID] = TRUE; // valid ����
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
					// �洢��CPU ����GPU
					for (int i = 0; i <= k - 1; i++)
					{
						int CandTrajID = candidateTrajID[qID][i];
						map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr.find(CandTrajID);
						// ����켣��û�б�����GPU��
						if (traID_baseAddr_ITER == traID_baseAddr.end())
						{
							int pointsNumInThisCand = 0;
							SPoint* thisTrajAddr = tempPtr;

							// ���� Lt �ؽ��켣 ��Ϊ��Ĵ洢�ǰ���cell���
							for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++) // Lt����cellBasedTrajectory
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

		// Step2��***********   �������ķ���û�д���ʼ����EDR  ***************

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

		// Step3��***********   ���м�������󣬸���worstNow�Լ�д���� *********************
		// ��֤��Сk��EDR
		for (int idx = 0; idx <= k * validQueryTraNum - 1; idx++)
		{
			int qID = taskInfoTable[idx].qID;
			int i = idx % k;
			if (numElemInCalculatedQueue[qID] < kValue)
			{
				//ֱ����PQ���
				FDwithID fd;
				fd.traID = candidateTrajID[qID][i];
				fd.FD = resultReturned[idx];
				// αFD ��ΪEDRCalculated->cmp������ͬ
				EDRCalculated[qID].push(fd);

				numElemInCalculatedQueue[qID]++;
			}
			else
			{
				//��һ���Ƿ��PQ����ã�����ǵ���һ����ģ�����ȥһ���õģ����򲻶����ȶ���Ҳ������worstNow��
				int worstInPQ = EDRCalculated[qID].top().FD; // ���EDR
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
			//��ȡtopk
			for (int i = 0; i <= k - 1; i++) {
				candidateTrajID[i] = queryQueue[qID].top().traID;
				//printf("%d,%d\t", queryQueue[qID].top().traID,queryQueue[qID].top().FD);
				queryQueue[qID].pop();
			}
			//EDR calculate
			//��һ������AllPoints����ȡ�����켣
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
			//�ڶ���������EDR
			int resultReturned[k];
			this->SimilarityExecuter(queryTra, candidateTra, qTra[qID].length, candidateTraLength, k, resultReturned);
			//����worstNow
			for (int i = 0; i <= k - 1; i++) {
				if (numElemInCalculatedQueue[qID] < kValue) {
					//ֱ����PQ���
					FDwithID fd;
					fd.traID = candidateTrajID[i];
					fd.FD = resultReturned[i];
					EDRCalculated[qID].push(fd);
					numElemInCalculatedQueue[qID]++;
				}
				else {
					//��һ���Ƿ��PQ����ã�����ǵ���һ����ģ�����ȥһ���õģ����򲻶����ȶ���Ҳ������worstNow��
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
			//���ֽ������ͷ��ڴ�
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



	//������
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
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024)); // 256�ֽ� ���1024*1024��
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
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //�����
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



	//������
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

	void* baseAddrGPU2 = NULL;// ��Ҫʹ��ȫ���ڴ��ַ
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU2), (long long int)BIG_MEM * 1024 * 1024));

	void* whileAddrGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024));
	void* whileAddrGPUBase = whileAddrGPU;

	//��ǰ���䵽�ĵ�ַ
	void* nowAddrGPU = NULL;

	cudaStream_t defaultStream;
	cudaStreamCreate(&defaultStream);

	MyTimer timer;

	// С��������С��FD ��Ĭ�� ע����priority_queue��ָ�� ���ƶ�ά
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	// ע����map��ָ�� ���ƶ�ά ��Ϊpriority_queue��map������������
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum]; // FVtable

	//Ϊ��ѯ����freqVector ���壺����ѯ�켣��freqVector��freqVectors[qID]  with CPU
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //�����
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
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join)); // threads_FD�Ѿ�ָ��

#else
	// ��CPU�߳�
	// ����FD ����queryQueue������ ��ÿ��queryQueue��priority_queue ���ȶ��� with CPU ���߳� ���߳��ܲ��� �����ǱʼǱ�ԭ��
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		this->FDCalculateParallelHandeler(&queryQueue[qID], &freqVectors[qID]);
	}
	// ���ص� queryQueue[]

#endif
	timer.stop();
	cout << "device: " << device_idx << "Part2 FD time:" << timer.elapse() << endl;

	timer.start();

	//MyTimer tt;
	//tt.start();

	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];

	int* numElemInCalculatedQueue = new int[queryTrajNum]; //���浱ǰ���ȶ��н������֤���ȶ��д�С������kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//׼����֮�󣬿�ʼ����ѯ


	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}
	//��ѯ�켣��Ϣ��׼����

	// Step0��*********** ����Tq��ѯ�Ĺ켣ȫ����Ϣ CPU->GPU *************

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
		// �浽queryTra��
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		// �ڴ濽�� ����ѯ�켣 CPU->GPU
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

	//��һ����ѭ������֦
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

	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum); // ��ά���� ָ��

	int* candidateTraLength = (int*)malloc(sizeof(int) * k * queryTrajNum);
	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
	
	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum); // û�õ��� candidateOffsets �غ�
	
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

		// Step1��*********** CPU->GPU  Tc ��صĿ��� MAT****************
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
								for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++) // Lt����cellBasedTrajectory table
								{
									int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
									memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
									tempPtr += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
									pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
								}
								CUDA_CALL(cudaMemcpyAsync(candidateTraGPU, thisTrajAddr, pointsNumInThisCand*sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream));
								
								traID_baseAddr[candidateTrajID[qID][i]] = candidateTraGPU;// ��¼

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


		//����candidateTraj��ɣ�����candidateTrajLength

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

		//�������ķ���û�д���ʼ����EDR
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

		//���м�������󣬸���worstNow�Լ�д����
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
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //�����
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
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join)); // threads_FD�Ѿ�ָ��

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
		// �浽queryTra��
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		// �ڴ濽�� ����ѯ�켣 CPU->GPU
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
					// ֹͣpop qID�������
					validQueryTraNum2--;
					isFinished[qID] = TRUE; // valid ����
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


		// Step1��*********** CPU->GPU  Tc ��صĿ��� MAT****************

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

					// �洢��CPU ����GPU
					for (int i = 0; i <= dynamicm - 1; i++)
					{
						int CandTrajID = candidateTrajID[qID][i];
						if (CandTrajID > 0) {
							map<int, void*>::iterator traID_baseAddr_ITER = traID_baseAddr.find(CandTrajID);
							if (traID_baseAddr_ITER == traID_baseAddr.end())
							{
								int pointsNumInThisCand = 0;
								SPoint* thisTrajAddr = tempPtr;
								for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++) // Lt����cellBasedTrajectory
								{
									int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
									// ��allPoints�����㵽CPU
									memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
									tempPtr += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
									pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
								}
								// �������켣��ȡ��candidateTraGPU��
								CUDA_CALL(cudaMemcpyAsync(candidateTraGPU, thisTrajAddr, pointsNumInThisCand*sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream));
								traID_baseAddr[candidateTrajID[qID][i]] = candidateTraGPU;
								candidateTraLength[dynamicm * validQueryIdx + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
								taskInfoTable[dynamicm * validQueryIdx + i].qID = qID;
								taskInfoTable[dynamicm * validQueryIdx + i].candTrajID = CandTrajID;
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].objectId = candidateTrajID[qID][i];// ����ѭ����ʼ��CandTrajID
								candidateTrajOffsetTable[dynamicm * validQueryIdx + i].addr = candidateTraGPU;
								candidateOffsets[dynamicm * validQueryIdx + i] = candidateTraGPU;
								candidateTraGPU = (candidateTraGPU + pointsNumInThisCand);
								nowAddrGPU = (void*)candidateTraGPU;
							}
							// ����ù켣�Ѿ����ƽ���gpu���棬��ôֻ��Ҫ���ոù켣id���±������
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
			if (qID >= 0) {// �߽�����
				int i = idx % dynamicm;
				if (numElemInCalculatedQueue[qID] < kValue)
				{
					//ֱ����PQ���
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



	//������
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

	void* baseAddrGPU2 = NULL;// ��Ҫʹ��ȫ���ڴ��ַ
	CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU2), (long long int)BIG_MEM * 1024 * 1024));

	void* whileAddrGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)(&whileAddrGPU), (long long int)SMALL_MEM * 1024 * 1024));
	void* whileAddrGPUBase = whileAddrGPU;

	//��ǰ���䵽�ĵ�ַ
	void* nowAddrGPU = NULL;

	cudaStream_t defaultStream;
	cudaStreamCreate(&defaultStream);

	MyTimer timer;

	// С��������С��FD ��Ĭ�� ע����priority_queue��ָ�� ���ƶ�ά
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	// ע����map��ָ�� ���ƶ�ά ��Ϊpriority_queue��map������������
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum]; // FVtable

																  //Ϊ��ѯ����freqVector ���壺����ѯ�켣��freqVector��freqVectors[qID]  with CPU
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //�����
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
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join)); // threads_FD�Ѿ�ָ��

#else
	// ��CPU�߳�
	// ����FD ����queryQueue������ ��ÿ��queryQueue��priority_queue ���ȶ��� with CPU ���߳� ���߳��ܲ��� �����ǱʼǱ�ԭ��
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		this->FDCalculateParallelHandeler(&queryQueue[qID], &freqVectors[qID]);
	}
	// ���ص� queryQueue[]

#endif
	timer.stop();
	cout << "device: " << device_idx << "Part2 FD time:" << timer.elapse() << endl;

	timer.start();

	//MyTimer tt;
	//tt.start();

	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];

	int* numElemInCalculatedQueue = new int[queryTrajNum]; //���浱ǰ���ȶ��н������֤���ȶ��д�С������kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//׼����֮�󣬿�ʼ����ѯ


	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}
	//��ѯ�켣��Ϣ��׼����

	// Step0��*********** ����Tq��ѯ�Ĺ켣ȫ����Ϣ CPU->GPU *************

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
		// �浽queryTra��
		for (int i = 0; i <= qTra[qID].length - 1; i++)
		{
			queryTra[i].x = qTra[qID].points[i].lon;
			queryTra[i].y = qTra[qID].points[i].lat;
			queryTra[i].tID = qTra[qID].tid;
		}
		// �ڴ濽�� ����ѯ�켣 CPU->GPU
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

	//��һ����ѭ������֦
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

	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum); // ��ά���� ָ��

	int* candidateTraLength = (int*)malloc(sizeof(int) * k * queryTrajNum);
	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);

	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum); // û�õ��� candidateOffsets �غ�

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

		// Step1��*********** CPU->GPU  Tc ��صĿ��� MAT****************
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
								for (int subID = 0; subID <= this->cellBasedTrajectory[candidateTrajID[qID][i]].length - 1; subID++) // Lt����cellBasedTrajectory table
								{
									int idxInAllPoints = this->cellBasedTrajectory[candidateTrajID[qID][i]].startIdx[subID];
									memcpy(tempPtr, &this->allPoints[idxInAllPoints], sizeof(SPoint) * this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID]);
									tempPtr += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
									pointsNumInThisCand += this->cellBasedTrajectory[candidateTrajID[qID][i]].numOfPointInCell[subID];
								}
								CUDA_CALL(cudaMemcpyAsync(candidateTraGPU, thisTrajAddr, pointsNumInThisCand * sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream));

								//traID_baseAddr[candidateTrajID[qID][i]] = candidateTraGPU;// ��¼

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


		//����candidateTraj��ɣ�����candidateTrajLength

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

		//�������ķ���û�д���ʼ����EDR
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

		//���м�������󣬸���worstNow�Լ�д����
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
	// ����FV FD ����

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

	//��ǰ���䵽�ĵ�ַ


	MyTimer timer;
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum]; 
	//Ϊ��ѯ����freqVector
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //�����
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
	//Ϊ��֦����Frequency Distance
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Part2 FD time:" << timer.elapse() << endl;

	//��һ�����ȶ��д洢��ǰ���Ž�����󶥶ѣ���֤��ʱ����pop����Ľ��
	timer.start();
	//MyTimer tt;
	//tt.start();
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];
	int* numElemInCalculatedQueue = new int[queryTrajNum]; //���浱ǰ���ȶ��н������֤���ȶ��д�С������kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//׼����֮�󣬿�ʼ����ѯ
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
		//��ѯ�켣��Ϣ��׼����
		//�����ѯ�Ĺ켣
		SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
		//������allQueryTra�и����켣��offset����ʼ��ַ��
		int* allQueryTraOffset = new int[queryNumEachGPU[device_idx]];
		SPoint* queryTra = allQueryTra;
		SPoint* queryTraGPU = (SPoint*)baseAddrSimi[device_idx];
		//������Ǳ�������queryTra�Ļ�ַ
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
		// queryTraOffsetGPU�Ǳ���queryTra��offset�Ļ���ַ
		int* queryTraOffsetGPU = (int*)nowAddrGPU;
		CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryNumEachGPU[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU + queryNumEachGPU[device_idx]);

		//����queryLength
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

	// ��ѯ Tq �켣��Ϣ��׼����
	// ��2��GPU������ȫһ������
	// �����ѯ�Ĺ켣
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//������allQueryTra�и����켣��offset����ʼ��ַ��
	int* allQueryTraOffset = new int[queryTrajNum];
	SPoint* queryTra = allQueryTra;

	SPoint* queryTraGPU[2];
	SPoint* queryTraGPUBase[2];
	for (int device_idx = 0; device_idx <= num_devices - 1; device_idx++)
	{
		queryTraGPU[device_idx] = (SPoint*)nowAddrGPU[device_idx];
		//������Ǳ�������queryTra�Ļ�ַ
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
		// queryTraOffsetGPU�Ǳ���queryTra��offset�Ļ���ַ
		queryTraOffsetGPU[device_idx] = (int*)nowAddrGPU[device_idx];
		CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU[device_idx], allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU[device_idx] + queryTrajNum);

		//����queryLength
		queryLengthGPU[device_idx] = (int*)nowAddrGPU[device_idx];
		CUDA_CALL(cudaMemcpyAsync(queryLengthGPU[device_idx], queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU[device_idx] + queryTrajNum);
	}


	//----------------------------------------------------

	/*
	 *
	// ��GPU�汾����
	//-----------------------------------------------
	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		totalQueryTrajLength += qTra[qID].length;
	}
	//��ѯ�켣��Ϣ��׼����
	//�����ѯ�Ĺ켣
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//������allQueryTra�и����켣��offset����ʼ��ַ��
	int* allQueryTraOffset = new int[queryTrajNum];
	SPoint* queryTra = allQueryTra;
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU;
	//������Ǳ�������queryTra�Ļ�ַ
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
	// queryTraOffsetGPU�Ǳ���queryTra��offset�Ļ���ַ
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);

	//����queryLength
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	//----------------------------------------------------
	//��GPU�汾����
	*/


	//tt.stop();
	//cout << "Part3.0.1 time:" << tt.elapse() << endl;
	//tt.start();
	//��һ����ѭ������֦
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
	//��GPU�汾����

	SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);

	//����qid��candID����candidateTran�е�offset�Ķ�Ӧ��ϵ
	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
	//�켣Ψһ��������id�¹켣��id��baseAddr ���б�Ҫ����켣��id�𣿣�
	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum);
	//����candidateOffset�Ļ���ַ��
	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum);
	int candidateTrajNum = 0;
	// traID����candidateTrajOffsetTable�е�idx�Ķ�Ӧ��ϵmap����Ҫ�����жϹ켣�Ƿ��Ѿ����Ƶ�gpu
	map<int, void*> traID_baseAddr;
	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU;
	
	//��GPU�汾����
	*/


	//tt.stop();
	//cout << "Part3.0.2 time:" << tt.elapse() << endl;

	while (!isAllFinished)
	{
		//tt.start();
		//���д�����ĺ�ѡ�켣������Ŀ���൱��������Ŀ
		int validCandTrajNum[2] = { 0 };


		// Ŀǰ��û�м�����ɵĲ�ѯ�켣����Ŀ
		int validQueryTraNum = queryTrajNum;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
		}


		// �ڱ��ּ����ڵĵڼ����켣
		int validQueryIdx = 0;
		int validQueryIdxGPU[2] = { 0 };
		int queryEachGPU[2];
		


		// ��GPU���ٵĹؼ� ÿ�δ�while 40*40(m) -> 20*40 20*40 ��������������֮ǰ��һ�� ���̶�k�����Կ�
		// ���񲻿ɿ� ����������̫��
		// ʣ40 20+20
		// ʣ30 15+15 Tq������ȫ���ɿ�

		queryEachGPU[0] = validQueryTraNum / 2;
		queryEachGPU[1] = validQueryTraNum - queryEachGPU[0];



		// ���Ե�ַ
		SPoint* tempPtr[2];
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			tempPtr[device_idx] = candidateTra[device_idx];

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			// ѡ�� GPU
			int device_idx = 0;
			if (validQueryIdx < queryEachGPU[0])
				device_idx = 0;
			else
				device_idx = 1;

			// δ����򿽱�Tc
			if (!isFinished[qID])
			{

				//��ȡtopk��Ϊ����켣
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
					// ����켣��û�б�����GPU��
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
						// �������켣��ȡ��candidateTraGPU��
						CUDA_CALL(cudaMemcpyAsync(candidateTraGPU[device_idx], thisTrajAddr, pointsNumInThisCand*sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream[device_idx]));
						traID_baseAddr[device_idx][candidateTrajID[qID][i]] = candidateTraGPU[device_idx];
						//������Ҫ�����query��candidateTraLength
						candidateTraLength[device_idx][k * validQueryIdxGPU[device_idx] + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
						//������Ҫ�����query��offset
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].qID = qID;
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].candTrajID = CandTrajID;
						// ����켣��Ӧ��addr
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].objectId = candidateTrajID[qID][i];
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].addr = candidateTraGPU[device_idx];
						candidateOffsets[device_idx][k * validQueryIdxGPU[device_idx] + i] = candidateTraGPU[device_idx];
						//��ַ��ǰ�ƶ���������һ�θ���
						candidateTraGPU[device_idx] = (candidateTraGPU[device_idx] + pointsNumInThisCand);
						// nowAddrGPU ʼ����ָ����һ�����е�GPU��ַ
						nowAddrGPU[device_idx] = (void*)candidateTraGPU[device_idx];
					}
					// ����ù켣�Ѿ����ƽ���gpu���棬��ôֻ��Ҫ���ոù켣id���±������
					else
					{
						void* baseAddrGPU = traID_baseAddr_ITER->second;
						//������Ҫ�����query��candidateTraLength
						candidateTraLength[device_idx][k * validQueryIdxGPU[device_idx] + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
						//������Ҫ�����query��offset
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].qID = qID;
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].candTrajID = CandTrajID;
						// ����켣��Ӧ��addr
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
			// nowAddrGPU ʼ����ָ����һ�����е�GPU��ַ
			whileAddrGPU[device_idx] = (void*)((int*)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//����TaskInfoTable
			taskInfoTableGPU[device_idx] = (TaskInfoTableForSimilarity*)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU[device_idx], taskInfoTable[device_idx], sizeof(TaskInfoTableForSimilarity)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU ʼ����ָ����һ�����е�GPU��ַ
			whileAddrGPU[device_idx] = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//����candidate�ĵ�ַ��gpu��
			candidateOffsetsGPU[device_idx]= (SPoint**)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU[device_idx], candidateOffsets[device_idx], sizeof(SPoint*)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU ʼ����ָ����һ�����е�GPU��ַ
			whileAddrGPU[device_idx] = (void*)((SPoint**)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//����candidateTraj��candidateLength��ɣ�׼������Similarity search
			//ֻ��Ҫ��ѯisFinishedΪfalse��queryTra����������Щ����ֱ�ӿ�offsetTableCandidateTra
			
			resultReturnedGPU[device_idx] = (int*)whileAddrGPU[device_idx];
			whileAddrGPU[device_idx] = (void*)((int*)whileAddrGPU[device_idx] + k * queryTrajNum);
		}
		//CUDA_CALL(cudaMalloc((void**)resultReturnedGPU, sizeof(int)*k*queryTrajNum));


		//tt.stop();
		//cout << "Part3.2 time:" << tt.elapse() << endl;
		//tt.start();
		//�������ķ���û�д���ʼ����EDR
		if (validQueryTraNum * k == (validCandTrajNum[0]+ validCandTrajNum[1]))
		{
			// GPU 1 2 ����˳�� ���ڲ�ͬ�Ĺ����� stream0 stream1 ��֤��GPU ͬʱ����
			// ���� stream ���� 2 GPU����
			for (int device_idx = 0; device_idx <= 1; device_idx++){// for ѭ�����п��� ��GPU-1 �� GPU-2 �����ܿ� ûʲôӰ�죡����
				CUDA_CALL(cudaSetDevice(device_idx));
				if(validCandTrajNum[device_idx]==0) // �߽�����
					continue;
				EDRDistance_Batch_Handler(validCandTrajNum[device_idx], taskInfoTableGPU[device_idx], queryTraGPUBase[device_idx], queryTraOffsetGPU[device_idx], candidateOffsetsGPU[device_idx], queryLengthGPU[device_idx], candidateTraLengthGPU[device_idx], resultReturnedGPU[device_idx], &defaultStream[device_idx]);
			}
			CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU[0], sizeof(int)*k*queryEachGPU[0], cudaMemcpyDeviceToHost, defaultStream[0])); // ע��������stream[0]�ı�־
			CUDA_CALL(cudaMemcpyAsync(resultReturned + k*queryEachGPU[0], resultReturnedGPU[1], sizeof(int)*k*queryEachGPU[1], cudaMemcpyDeviceToHost, defaultStream[1])); // ע��������stream[1]�ı�־
			
			// �������GPU1 GPU2 ����ʽͬ������Ӧ����

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
				//ֱ����PQ���
				FDwithID fd;
				fd.traID = candidateTrajID[qID][i];
				fd.FD = resultReturned[idx];
				EDRCalculated[qID].push(fd);
				numElemInCalculatedQueue[qID]++;
			}
			else
			{
				//��һ���Ƿ��PQ����ã�����ǵ���һ����ģ�����ȥһ���õģ����򲻶����ȶ���Ҳ������worstNow��
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

		// �������ж� �Ƚ����� �ɶ���
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

		//GPUָ��ص�while��ʼ�ĵط�
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			whileAddrGPU[device_idx] = whileAddrGPUBase[device_idx];

		//tt.stop();
		//cout << "Part3.4 time:" << tt.elapse() << endl;
	} // end of BIG while


	timer.stop();
	cout << "Part3 EDR time:" << timer.elapse() << endl;

	//������
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
	// ����FV FD ����

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

	//��ǰ���䵽�ĵ�ַ


	MyTimer timer;
	priority_queue<FDwithID, vector<FDwithID>, cmp>* queryQueue = new priority_queue<FDwithID, vector<FDwithID>, cmp>[queryTrajNum];
	map<int, int>* freqVectors = new map<int, int>[queryTrajNum];
	//Ϊ��ѯ����freqVector
	timer.start();
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		for (int pID = 0; pID <= qTra[qID].length - 1; pID++)
		{
			int cellid = WhichCellPointIn(SamplePoint(qTra[qID].points[pID].lon, qTra[qID].points[pID].lat, 1, 1));
			int vituralCellNo = cellid >> VITURAL_CELL_PARAM; //�����
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
	//Ϊ��֦����Frequency Distance
	vector<thread> threads_FD;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		threads_FD.push_back(thread(std::mem_fn(&Grid::FDCalculateParallelHandeler), this, &queryQueue[qID], &freqVectors[qID]));
	}
	std::for_each(threads_FD.begin(), threads_FD.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Part2 FD time:" << timer.elapse() << endl;

	//��һ�����ȶ��д洢��ǰ���Ž�����󶥶ѣ���֤��ʱ����pop����Ľ��
	timer.start();
	//MyTimer tt;
	//tt.start();
	priority_queue<FDwithID, vector<FDwithID>, cmpBig>* EDRCalculated = new priority_queue<FDwithID, vector<FDwithID>, cmpBig>[queryTrajNum];
	int* numElemInCalculatedQueue = new int[queryTrajNum]; //���浱ǰ���ȶ��н������֤���ȶ��д�С������kValue
	for (int i = 0; i <= queryTrajNum - 1; i++)
		numElemInCalculatedQueue[i] = 0;

	//׼����֮�󣬿�ʼ����ѯ
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
	//��ѯ�켣��Ϣ��׼����
	//�����ѯ�Ĺ켣
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//������allQueryTra�и����켣��offset����ʼ��ַ��
	int* allQueryTraOffset = new int[queryNumEachGPU[device_idx]];
	SPoint* queryTra = allQueryTra;
	SPoint* queryTraGPU = (SPoint*)baseAddrSimi[device_idx];
	//������Ǳ�������queryTra�Ļ�ַ
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
	// queryTraOffsetGPU�Ǳ���queryTra��offset�Ļ���ַ
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryNumEachGPU[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
	nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU + queryNumEachGPU[device_idx]);

	//����queryLength
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

	// ��ѯ Tq �켣��Ϣ��׼����
	// ��2��GPU������ȫһ������
	// �����ѯ�Ĺ켣
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//������allQueryTra�и����켣��offset����ʼ��ַ��
	int* allQueryTraOffset = new int[queryTrajNum];
	SPoint* queryTra = allQueryTra;

	SPoint* queryTraGPU[2];
	SPoint* queryTraGPUBase[2];
	for (int device_idx = 0; device_idx <= num_devices - 1; device_idx++)
	{
		queryTraGPU[device_idx] = (SPoint*)nowAddrGPU[device_idx];
		//������Ǳ�������queryTra�Ļ�ַ
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
		// queryTraOffsetGPU�Ǳ���queryTra��offset�Ļ���ַ
		queryTraOffsetGPU[device_idx] = (int*)nowAddrGPU[device_idx];
		CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU[device_idx], allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU[device_idx] + queryTrajNum);

		//����queryLength
		queryLengthGPU[device_idx] = (int*)nowAddrGPU[device_idx];
		CUDA_CALL(cudaMemcpyAsync(queryLengthGPU[device_idx], queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream[device_idx]));
		nowAddrGPU[device_idx] = (void*)((int*)nowAddrGPU[device_idx] + queryTrajNum);
	}


	//----------------------------------------------------

	/*
	*
	// ��GPU�汾����
	//-----------------------------------------------
	int totalQueryTrajLength = 0;
	for (int qID = 0; qID <= queryTrajNum - 1; qID++)
	{
	totalQueryTrajLength += qTra[qID].length;
	}
	//��ѯ�켣��Ϣ��׼����
	//�����ѯ�Ĺ켣
	SPoint* allQueryTra = (SPoint*)malloc(sizeof(SPoint) * totalQueryTrajLength);
	//������allQueryTra�и����켣��offset����ʼ��ַ��
	int* allQueryTraOffset = new int[queryTrajNum];
	SPoint* queryTra = allQueryTra;
	SPoint* queryTraGPU = (SPoint*)baseAddrGPU;
	//������Ǳ�������queryTra�Ļ�ַ
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
	// queryTraOffsetGPU�Ǳ���queryTra��offset�Ļ���ַ
	int* queryTraOffsetGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryTraOffsetGPU, allQueryTraOffset, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);

	//����queryLength
	int* queryLengthGPU = (int*)nowAddrGPU;
	CUDA_CALL(cudaMemcpyAsync(queryLengthGPU, queryTraLength, sizeof(int)*queryTrajNum, cudaMemcpyHostToDevice, defaultStream));
	nowAddrGPU = (void*)((int*)nowAddrGPU + queryTrajNum);
	//----------------------------------------------------
	//��GPU�汾����
	*/


	//tt.stop();
	//cout << "Part3.0.1 time:" << tt.elapse() << endl;
	//tt.start();
	//��һ����ѭ������֦
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
	//��GPU�汾����

	SPoint* candidateTra = (SPoint*)malloc(sizeof(SPoint) * k * queryTrajNum * MAXLENGTH);

	//����qid��candID����candidateTran�е�offset�Ķ�Ӧ��ϵ
	TaskInfoTableForSimilarity* taskInfoTable = (TaskInfoTableForSimilarity *)malloc(sizeof(TaskInfoTableForSimilarity) * k * queryTrajNum);
	//�켣Ψһ��������id�¹켣��id��baseAddr ���б�Ҫ����켣��id�𣿣�
	OffsetTable* candidateTrajOffsetTable = (OffsetTable*)malloc(sizeof(OffsetTable) * k * queryTrajNum);
	//����candidateOffset�Ļ���ַ��
	SPoint** candidateOffsets = (SPoint**)malloc(sizeof(SPoint*) * k * queryTrajNum);
	int candidateTrajNum = 0;
	// traID����candidateTrajOffsetTable�е�idx�Ķ�Ӧ��ϵmap����Ҫ�����жϹ켣�Ƿ��Ѿ����Ƶ�gpu
	map<int, void*> traID_baseAddr;
	SPoint* candidateTraGPU = (SPoint*)nowAddrGPU;

	//��GPU�汾����
	*/


	//tt.stop();
	//cout << "Part3.0.2 time:" << tt.elapse() << endl;

	while (!isAllFinished)
	{
		//tt.start();
		//���д�����ĺ�ѡ�켣������Ŀ���൱��������Ŀ
		int validCandTrajNum[2] = { 0 };


		// Ŀǰ��û�м�����ɵĲ�ѯ�켣����Ŀ
		int validQueryTraNum = queryTrajNum;
		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			if (isFinished[qID])
				validQueryTraNum--;
		}


		// �ڱ��ּ����ڵĵڼ����켣
		int validQueryIdx = 0;
		int validQueryIdxGPU[2] = { 0 };
		int queryEachGPU[2];



		// ��GPU���ٵĹؼ� ÿ�δ�while 40*40(m) -> 20*40 20*40 ��������������֮ǰ��һ�� ���̶�k�����Կ�
		// ���񲻿ɿ� ����������̫��
		// ʣ40 20+20
		// ʣ30 15+15 Tq������ȫ���ɿ�

		queryEachGPU[0] = validQueryTraNum / 2;
		queryEachGPU[1] = validQueryTraNum - queryEachGPU[0];



		// ���Ե�ַ
		SPoint* tempPtr[2];
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			tempPtr[device_idx] = candidateTra[device_idx];

		for (int qID = 0; qID <= queryTrajNum - 1; qID++)
		{
			// ѡ�� GPU
			int device_idx = 0;
			if (validQueryIdx < queryEachGPU[0])
				device_idx = 0;
			else
				device_idx = 1;

			// δ����򿽱�Tc
			if (!isFinished[qID])
			{

				//��ȡtopk��Ϊ����켣
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
					// ����켣��û�б�����GPU��
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
						// �������켣��ȡ��candidateTraGPU��
						CUDA_CALL(cudaMemcpyAsync(candidateTraGPU[device_idx], thisTrajAddr, pointsNumInThisCand * sizeof(SPoint), cudaMemcpyHostToDevice, defaultStream[device_idx]));
						
						//traID_baseAddr[device_idx][candidateTrajID[qID][i]] = candidateTraGPU[device_idx];
						
						//������Ҫ�����query��candidateTraLength
						candidateTraLength[device_idx][k * validQueryIdxGPU[device_idx] + i] = this->cellBasedTrajectory[candidateTrajID[qID][i]].trajLength;
						//������Ҫ�����query��offset
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].qID = qID;
						taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].candTrajID = CandTrajID;
						// ����켣��Ӧ��addr
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].objectId = candidateTrajID[qID][i];
						candidateTrajOffsetTable[device_idx][k * validQueryIdxGPU[device_idx] + i].addr = candidateTraGPU[device_idx];
						candidateOffsets[device_idx][k * validQueryIdxGPU[device_idx] + i] = candidateTraGPU[device_idx];
						//��ַ��ǰ�ƶ���������һ�θ���
						candidateTraGPU[device_idx] = (candidateTraGPU[device_idx] + pointsNumInThisCand);
						// nowAddrGPU ʼ����ָ����һ�����е�GPU��ַ
						nowAddrGPU[device_idx] = (void*)candidateTraGPU[device_idx];
					}
					// ����ù켣�Ѿ����ƽ���gpu���棬��ôֻ��Ҫ���ոù켣id���±������

					//else
					//{
					//	void* baseAddrGPU = traID_baseAddr_ITER->second;
					//	//������Ҫ�����query��candidateTraLength
					//	candidateTraLength[device_idx][k * validQueryIdxGPU[device_idx] + i] = this->cellBasedTrajectory[CandTrajID].trajLength;
					//	//������Ҫ�����query��offset
					//	taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].qID = qID;
					//	taskInfoTable[device_idx][k * validQueryIdxGPU[device_idx] + i].candTrajID = CandTrajID;
					//	// ����켣��Ӧ��addr
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
			// nowAddrGPU ʼ����ָ����һ�����е�GPU��ַ
			whileAddrGPU[device_idx] = (void*)((int*)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//����TaskInfoTable
			taskInfoTableGPU[device_idx] = (TaskInfoTableForSimilarity*)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(taskInfoTableGPU[device_idx], taskInfoTable[device_idx], sizeof(TaskInfoTableForSimilarity)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU ʼ����ָ����һ�����е�GPU��ַ
			whileAddrGPU[device_idx] = (void*)((TaskInfoTableForSimilarity*)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//����candidate�ĵ�ַ��gpu��
			candidateOffsetsGPU[device_idx] = (SPoint**)whileAddrGPU[device_idx];
			CUDA_CALL(cudaMemcpyAsync(candidateOffsetsGPU[device_idx], candidateOffsets[device_idx], sizeof(SPoint*)*validCandTrajNum[device_idx], cudaMemcpyHostToDevice, defaultStream[device_idx]));
			// nowAddrGPU ʼ����ָ����һ�����е�GPU��ַ
			whileAddrGPU[device_idx] = (void*)((SPoint**)whileAddrGPU[device_idx] + validCandTrajNum[device_idx]);

			//����candidateTraj��candidateLength��ɣ�׼������Similarity search
			//ֻ��Ҫ��ѯisFinishedΪfalse��queryTra����������Щ����ֱ�ӿ�offsetTableCandidateTra

			resultReturnedGPU[device_idx] = (int*)whileAddrGPU[device_idx];
			whileAddrGPU[device_idx] = (void*)((int*)whileAddrGPU[device_idx] + k * queryTrajNum);
		}
		//CUDA_CALL(cudaMalloc((void**)resultReturnedGPU, sizeof(int)*k*queryTrajNum));


		//tt.stop();
		//cout << "Part3.2 time:" << tt.elapse() << endl;
		//tt.start();
		//�������ķ���û�д���ʼ����EDR
		if (validQueryTraNum * k == (validCandTrajNum[0] + validCandTrajNum[1]))
		{
			// GPU 1 2 ����˳�� ���ڲ�ͬ�Ĺ����� stream0 stream1 ��֤��GPU ͬʱ����
			// ���� stream ���� 2 GPU����
			for (int device_idx = 0; device_idx <= 1; device_idx++) {// for ѭ�����п��� ��GPU-1 �� GPU-2 �����ܿ� ûʲôӰ�죡����
				CUDA_CALL(cudaSetDevice(device_idx));
				if (validCandTrajNum[device_idx] == 0) // �߽�����
					continue;
				EDRDistance_Batch_Handler(validCandTrajNum[device_idx], taskInfoTableGPU[device_idx], queryTraGPUBase[device_idx], queryTraOffsetGPU[device_idx], candidateOffsetsGPU[device_idx], queryLengthGPU[device_idx], candidateTraLengthGPU[device_idx], resultReturnedGPU[device_idx], &defaultStream[device_idx]);
			}
			CUDA_CALL(cudaMemcpyAsync(resultReturned, resultReturnedGPU[0], sizeof(int)*k*queryEachGPU[0], cudaMemcpyDeviceToHost, defaultStream[0])); // ע��������stream[0]�ı�־
			CUDA_CALL(cudaMemcpyAsync(resultReturned + k*queryEachGPU[0], resultReturnedGPU[1], sizeof(int)*k*queryEachGPU[1], cudaMemcpyDeviceToHost, defaultStream[1])); // ע��������stream[1]�ı�־

																																										   // �������GPU1 GPU2 ����ʽͬ������Ӧ����

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
				//ֱ����PQ���
				FDwithID fd;
				fd.traID = candidateTrajID[qID][i];
				fd.FD = resultReturned[idx];
				EDRCalculated[qID].push(fd);
				numElemInCalculatedQueue[qID]++;
			}
			else
			{
				//��һ���Ƿ��PQ����ã�����ǵ���һ����ģ�����ȥһ���õģ����򲻶����ȶ���Ҳ������worstNow��
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

		// �������ж� �Ƚ����� �ɶ���
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

		//GPUָ��ص�while��ʼ�ĵط�
		for (int device_idx = 0; device_idx <= 1; device_idx++)
			whileAddrGPU[device_idx] = whileAddrGPUBase[device_idx];

		//tt.stop();
		//cout << "Part3.4 time:" << tt.elapse() << endl;
	} // end of BIG while


	timer.stop();
	cout << "Part3 EDR time:" << timer.elapse() << endl;

	//������
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



// CoarseGrainedThreadJoin GAT-S-noE + GAT-S-E ���˴��У� ��û�У�
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
	stateTableSim[1].resize(simNumGPU[1] * kValue); // ���������Ľ��

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

	// ��������thread �����Ĭ�ϵ�ͬ������ ���Ի���balanced�ȽϺ�

	timer.stop();
	cout << "Coarse Dual GPU Time:" << timer.elapse() << endl; // ���� FV FD EDR ֻ��EDR �����Ǽ���һ������ʱ�䣡��



	cout << "noMAT noE\n";
	// noE or E, decide here

	vector<thread> threads_RQ2;
	vector<int> stateTableSim2[2];
	stateTableSim2[0].resize(simNumGPU[0] * kValue);
	stateTableSim2[1].resize(simNumGPU[1] * kValue); // ���������Ľ��

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

	// ��������thread �����Ĭ�ϵ�ͬ������ ���Ի���balanced�ȽϺ�
	timer.stop();
	cout << "Coarse Dual GPU Time:" << timer.elapse() << endl; // ���� FV FD EDR ֻ��EDR �����Ǽ���һ������ʱ�䣡��



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
//		candidateTra[i] = (SPoint*)malloc(sizeof(SPoint)*(candTra[i]->length)); //���Ե�ʱ����һ�����ܱ��ڴ����FFFFF
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
//		//ÿ��DP����
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
//				SPoint p2 = tra2[col - 1]; //�������ڴ��Ǿۼ����ʵ���
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



