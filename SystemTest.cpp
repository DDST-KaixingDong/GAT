#include "SystemTest.h"


using namespace std;
#define sleep(x) Sleep(x*1000)

bool nonmem_cmp2(MBB t1, MBB t2)
{
	if (t1.GetMBBArea() < t2.GetMBBArea())
		return true;
	return false;
}

// only 2 GPUs!!
void ReorderArray(Trajectory* pData, int length)
{
	//if (length % 2) return;
	Trajectory* ouT = new Trajectory[int(length+1) / 2];
	int i = 0;
	Trajectory* jiT = new Trajectory[length - int(length + 1) / 2];
	for (i = 0; i < (length+1) / 2; i++) {
		ouT[i] = pData[2 * i];
		if(length % 2 == 0) jiT[i] = pData[2 * i + 1];
	}

	int p = 0;
	for (i = 0; i < int(length + 1) / 2; i++) {
		pData[p++] = ouT[i];
	}
	for (i = 0; i < length - int(length + 1) / 2; i++) {
		pData[p++] = jiT[i];
	}

	delete[] ouT;
	delete[] jiT;
}

SystemTest::SystemTest()
{
}

SystemTest::~SystemTest()
{
}

SystemTest::SystemTest(Trajectory* tradb, Grid* g, STIG *stig, FSG* fsg, MortonGrid* mgrid)
{
	this->tradb = tradb;
	this->g = g;
	this->stig = stig;
	this->fsg = fsg;
	this->mgrid = mgrid;
}

int SystemTest::STIGrangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)
{
	CUDA_CALL(cudaSetDevice(0));
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	vector<CPURangeQueryResult> resultTable;
	resultTable.resize(rangeQueryNum);
	MBB mbbArray[5000];
	int* resultSize = NULL;
	for (int i = 0; i <= 4999; i++)
		mbbArray[i] = rangeQueryMBB;
	MyTimer timer;

	// single GPU
	printf("********single GPU STIG range query #query=%d:\n", rangeQueryNum);

	void* allocatedGPUMemS = 0;
	CUDA_CALL(cudaMalloc((void**)&this->stig->baseAddrGPU[0], (long long int)BIG_MEM * 1024 * 1024));
	CUDA_CALL(cudaMalloc((void**)&this->stig->stateTableGPU[0], (long long int)SMALL_MEM * 1024 * 1024));
	allocatedGPUMemS = this->stig->baseAddrGPU[0];
	timer.start();
	stig->rangeQueryGPU(mbbArray, rangeQueryNum, &resultTable[0], resultSize, 0);
	timer.stop();
	cout << "single GPU Time of STIG:" << timer.elapse() << "ms" << endl;
	CUDA_CALL(cudaFree(allocatedGPUMemS));
	CUDA_CALL(cudaFree(this->stig->stateTableGPU[0]));

#ifdef USE_MULTIGPU
	// multi-GPU
	printf("********Dual GPU STIG range query #query=%d:\n", rangeQueryNum);
	int device_num = 2;
	vector<thread> threads_RQ;
	int rangeNumGPU[2];
	rangeNumGPU[0] = rangeQueryNum / 2;
	rangeNumGPU[1] = rangeQueryNum - rangeNumGPU[0];
	int startIdx[2];
	startIdx[0] = 0;
	startIdx[1] = rangeNumGPU[0];
	void* allocatedGPUMem[2] = { NULL };
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		// this->freqVectors.formPriorityQueue(&queryQueue[qID], &freqVectors[qID]);
		CUDA_CALL(cudaSetDevice(device_idx));
		CUDA_CALL(cudaMalloc((void**)&this->stig->baseAddrGPU[device_idx], (long long int)BIG_MEM * 1024 * 1024));
		CUDA_CALL(cudaMalloc((void**)&this->stig->stateTableGPU[device_idx], (long long int)SMALL_MEM * 1024 * 1024));
		allocatedGPUMem[device_idx] = this->stig->baseAddrGPU[device_idx];
		threads_RQ.push_back(thread(std::mem_fn(&STIG::rangeQueryGPU), this->stig, &mbbArray[startIdx[device_idx
		]], rangeNumGPU[device_idx], &resultTable[startIdx[1]], resultSize, device_idx));
	}
	timer.start();
	std::for_each(threads_RQ.begin(), threads_RQ.end(), std::mem_fn(&std::thread::join));
	timer.stop();
	cout << "Dual GPU Time of STIG:" << timer.elapse() << "ms" << endl;
	for (int device_idx = 0; device_idx <= device_num - 1; device_idx++)
	{
		CUDA_CALL(cudaFree(allocatedGPUMem[device_idx]));
		CUDA_CALL(cudaFree(this->stig->stateTableGPU[device_idx]));
	}
#else

#endif
#ifdef CHECK_CORRECT
	FILE* fp = fopen("STIGResult.txt", "w+");
	for (int i = 0; i <= rangeQueryNum - 1; i++)
	{
		for (int traID = 1; traID <= this->stig->maxTid; traID++) {
			if (resultTable[i][traID])
				fprintf(fp, "Query %d result: %d\n", i, traID);
		}
	}
	fclose(fp);
#endif
	return 0;
}


int SystemTest::FSGrangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)
{
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	MBB mbbArray[5000];
	int* resultSize = NULL;
	for (int i = 0; i <= 4999; i++)
		//rangeQueryMBB.randomGenerateMBB(mbbArray[i]);
		mbbArray[i] = rangeQueryMBB;
	MyTimer timer;
	vector<CPURangeQueryResult> resultTable;
	resultTable.resize(rangeQueryNum);

	printf("******** single GPU FSG range query #query=%d:\n", rangeQueryNum);
	CUDA_CALL(cudaSetDevice(0));

#ifdef WIN32
	CUDA_CALL(cudaMalloc((void**)(&fsg->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#else
	CUDA_CALL(cudaMalloc((void**)(&fsg->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#endif

	void *allocatedGPUMem = fsg->baseAddrRange[0];
	CUDA_CALL(cudaMalloc((void**)&fsg->stateTableGPU[0], (long long int)SMALL_MEM * 1024 * 1024));
	vector<RangeQueryStateTable> stateTableRange;
	stateTableRange.resize(rangeQueryNum * 50000);

	// why 拆分 有必要么
	const int ONCE_QUERY_NUM = rangeQueryNum;
	timer.start();
	for (int queryIdx = 0; queryIdx < rangeQueryNum; queryIdx += ONCE_QUERY_NUM) {
		int querySize = (queryIdx + ONCE_QUERY_NUM < rangeQueryNum) ? ONCE_QUERY_NUM : rangeQueryNum - queryIdx;
		fsg->baseAddrRange[0] = allocatedGPUMem;
		fsg->rangeQueryBatchGPU(&mbbArray[queryIdx], querySize, &resultTable[queryIdx], resultSize, &stateTableRange[queryIdx], 0);
	}
	timer.stop();
	cout << "Single GPU Time of FSG:" << timer.elapse() << "ms" << endl;
	CUDA_CALL(cudaFree(allocatedGPUMem));
	CUDA_CALL(cudaFree(fsg->stateTableGPU[0]));

#ifdef USE_MULTIGPU
	printf("******** multi-GPU range query FSG #query=%d:\n", rangeQueryNum);
	fsg->rangeQueryBatchMultiGPU(mbbArray, rangeQueryNum, &resultTable[0], resultSize);
#else

#endif
#ifdef CHECK_CORRECT
	FILE* fp = fopen("FSGResult.txt", "w+");
	for (int i = 0; i <= rangeQueryNum - 1; i++)
	{
		for (int traID = 1; traID <= this->fsg->trajNum; traID++) {
			if (resultTable[i][traID])
				fprintf(fp, "Query %d result: %d\n", i, traID);
		}
	}
	fclose(fp);
#endif
	return 0;
}


int SystemTest::rangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)
{
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	MBB mbbArray[5000];
	for (int i = 0; i <= 4999; i++)
		//rangeQueryMBB.randomGenerateMBB(mbbArray[i]);
		mbbArray[i] = rangeQueryMBB;
	MyTimer timer;
	CPURangeQueryResult* resultTable = NULL;
	int* resultSize = NULL;
	
	printf("********single-core CPU range query #query=%d:\n", rangeQueryNum);
	vector<CPURangeQueryResult> rangeQueryResult;
	rangeQueryResult.resize(rangeQueryNum);
	timer.start();
	g->rangeQueryBatch(mbbArray, rangeQueryNum, &rangeQueryResult[0], resultSize);
	timer.stop();
	cout << "single-core CPU Time:" << timer.elapse() << "ms" << endl;

	printf("********multi-core CPU range query #query=%d:\n", rangeQueryNum);
	vector<CPURangeQueryResult> rangeQueryResultMultiCPU;
	rangeQueryResultMultiCPU.resize(rangeQueryNum);
	timer.start();
	g->rangeQueryBatchMultiThread(mbbArray, rangeQueryNum, &rangeQueryResultMultiCPU[0], resultSize);
	timer.stop();
	cout << "multi-core CPU Time:" << timer.elapse() << "ms" << endl;
	
	printf("********single GPU range query #query=%d:\n", rangeQueryNum);
	CUDA_CALL(cudaSetDevice(0));
	vector<CPURangeQueryResult> rangeQueryResultGPU;
	rangeQueryResultGPU.resize(rangeQueryNum);

#ifdef WIN32
	CUDA_CALL(cudaMalloc((void**)(&g->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));// byte为单位
#else
	CUDA_CALL(cudaMalloc((void**)(&g->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#endif

	void *allocatedGPUMem = g->baseAddrRange[0];
	CUDA_CALL(cudaMalloc((void**)&g->stateTableGPU[0], (long long int)SMALL_MEM * 1024 * 1024));
	vector<RangeQueryStateTable> stateTableRange;
	stateTableRange.resize(rangeQueryNum *50000);
	timer.start();

	// 这里没有分开batch做

	g->rangeQueryBatchGPU(mbbArray, rangeQueryNum, &rangeQueryResultGPU[0], resultSize, &stateTableRange[0], 0);
	timer.stop();
	cout << "Single GPU Time:" << timer.elapse() << "ms" << endl;
	CUDA_CALL(cudaFree(allocatedGPUMem));
	CUDA_CALL(cudaFree(g->stateTableGPU[0]));

#ifdef CHECK_CORRECT
	FILE* fp = fopen("GPUResult.txt", "w+");
	for (int i = 0; i <= rangeQueryNum - 1; i++)
	{
		for (int traID = 1; traID <= this->g->trajNum; traID++) {
			if (rangeQueryResultGPU[i][traID])
				fprintf(fp, "Query %d result: %d\n", i, traID);
		}
	}
	fclose(fp);
#endif
#ifdef USE_MULTIGPU
	printf("******** rangeQueryTest:: multi-GPU range query #query=%d:\n", rangeQueryNum);
	vector<CPURangeQueryResult> rangeQueryResultGPUs;
	rangeQueryResultGPUs.resize(rangeQueryNum);
	g->rangeQueryBatchMultiGPU(mbbArray, rangeQueryNum, &rangeQueryResultGPUs[0], resultSize);
#else

#endif
	
	return 0;
}

int SystemTest::rangeQueryTestWithoutMorton(MBB rangeQueryMBB, int rangeQueryNum)
{
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	CPURangeQueryResult* resultTable = NULL;
	MBB mbbArray[5000];
	int* resultSize = NULL;
	for (int i = 0; i <= 4999; i++)
		//rangeQueryMBB.randomGenerateMBB(mbbArray[i]);
		mbbArray[i] = rangeQueryMBB;
	MyTimer timer;



	printf("******** single GPU range query without Morton #query=%d:\n", rangeQueryNum);
	CUDA_CALL(cudaSetDevice(0));
	vector<CPURangeQueryResult> rangeQueryResultGPU;
	rangeQueryResultGPU.resize(rangeQueryNum);

#ifdef WIN32
	CUDA_CALL(cudaMalloc((void**)(&g->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#else
	CUDA_CALL(cudaMalloc((void**)(&g->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#endif

	void *allocatedGPUMem = g->baseAddrRange[0];
	CUDA_CALL(cudaMalloc((void**)&g->stateTableGPU[0], (long long int)SMALL_MEM * 1024 * 1024));
	vector<RangeQueryStateTable> stateTableRange;
	stateTableRange.resize(rangeQueryNum * 50000);
	timer.start();
	g->rangeQueryBatchGPUWithoutMorton(mbbArray, rangeQueryNum, &rangeQueryResultGPU[0], resultSize, &stateTableRange[0], 0);
	timer.stop();
	cout << "Single GPU Time without Morton:" << timer.elapse() << "ms" << endl;
	CUDA_CALL(cudaFree(allocatedGPUMem));
	CUDA_CALL(cudaFree(g->stateTableGPU[0]));
#ifdef CHECK_CORRECT
	FILE* fp = fopen("GPUResult.txt", "w+");
	for (int i = 0; i <= rangeQueryNum - 1; i++)
	{
		for (int traID = 1; traID <= this->g->trajNum; traID++) {
			if (rangeQueryResultGPU[i][traID])
				fprintf(fp, "Query %d result: %d\n", i, traID);
		}
	}
	fclose(fp);
#endif

#ifdef USE_MULTIGPU
	printf("******** multi-GPU range query without Morton #query=%d:\n", rangeQueryNum);
	vector<CPURangeQueryResult> rangeQueryResultGPUs;
	rangeQueryResultGPUs.resize(rangeQueryNum);
	g->rangeQueryBatchMultiGPUWithoutMorton(mbbArray, rangeQueryNum, &rangeQueryResultGPUs[0], resultSize);
#else

#endif

	return 0;
}
 

int SystemTest::MortonGridRangeQueryTest(MBB rangeQueryMBB, int rangeQueryNum)
{
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	MBB mbbArray[5000]; 	
	for (int i = 0; i <= 4999; i++)
		mbbArray[i] = rangeQueryMBB;

	// too early !! can be in mgrid->rangeQueryBatchGPU
	CUDA_CALL(cudaSetDevice(0));
#ifdef WIN32
	CUDA_CALL(cudaMalloc((void**)(&mgrid->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#else
	CUDA_CALL(cudaMalloc((void**)(&mgrid->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#endif
	CUDA_CALL(cudaMalloc((void**)&mgrid->stateTableGPU[0], (long long int)SMALL_MEM * 1024 * 1024));
	void * allocatedGPUMem = mgrid->baseAddrRange[0];


	MyTimer timer;


	const int ONCE_QUERY_NUM = rangeQueryNum;
	int queryIdx;

	printf("\n******** single GPU Morton Grid range query #query=%d:\n", rangeQueryNum);
	int* resultSize = NULL;
	vector<CPURangeQueryResult> resultTable;// bool-vector的vector 二维数组还可以这样定义
	resultTable.resize(rangeQueryNum);
	vector<RangeQueryStateTable> stateTableRange;
	stateTableRange.resize(rangeQueryNum * 50000);


	timer.start(); 
	for (queryIdx = 0; queryIdx < rangeQueryNum; queryIdx += ONCE_QUERY_NUM) {
		int querySize = (queryIdx + ONCE_QUERY_NUM < rangeQueryNum) ? ONCE_QUERY_NUM : rangeQueryNum - queryIdx;
		mgrid->baseAddrRange[0] = allocatedGPUMem;//每次重新开始
		mgrid->rangeQueryBatchGPU(&mbbArray[queryIdx], querySize, &resultTable[queryIdx], resultSize, &stateTableRange[queryIdx], 0);
	}
	timer.stop();
	cout << "Single GPU Time of Morton Grid:" << timer.elapse() << "ms" << endl;


	// seems no need no 2 GPU
	printf("\n******** single GPU Morton Grid range query noMAT #query=%d:\n", rangeQueryNum);
	int* resultSize2 = NULL;
	vector<CPURangeQueryResult> resultTable2;// bool-vector的vector 二维数组还可以这样定义
	resultTable2.resize(rangeQueryNum);
	vector<RangeQueryStateTable> stateTableRange2;
	stateTableRange2.resize(rangeQueryNum * 50000);

	timer.start();
	for (queryIdx = 0; queryIdx < rangeQueryNum; queryIdx += ONCE_QUERY_NUM) {
		int querySize = (queryIdx + ONCE_QUERY_NUM < rangeQueryNum) ? ONCE_QUERY_NUM : rangeQueryNum - queryIdx;
		mgrid->baseAddrRange[0] = allocatedGPUMem;//每次重新开始
		mgrid->rangeQueryBatchGPUNoMAT(&mbbArray[queryIdx], querySize, &resultTable2[queryIdx], resultSize2, &stateTableRange2[queryIdx], 0);
	}
	timer.stop();
	cout << "Single GPU Time of Morton Grid NoMAT:" << timer.elapse() << "ms" << endl;


	CUDA_CALL(cudaFree(allocatedGPUMem)); 
	CUDA_CALL(cudaFree(mgrid->stateTableGPU[0]));

#ifdef USE_MULTIGPU
	printf("\n******** multi-GPU range query Morton Grid #query=%d:\n", rangeQueryNum);
	mgrid->rangeQueryBatchMultiGPU(mbbArray, rangeQueryNum, &resultTable[0], resultSize); // only one resultTable
#else

#endif

#ifdef CHECK_CORRECT
	FILE* fp = fopen("MortonResult.txt", "w+");
	for (int i = 0; i <= rangeQueryNum - 1; i++)
	{
		for (int traID = 1; traID <= this->fsg->trajNum; traID++) {
			if (resultTable[i][traID])
				fprintf(fp, "Query %d result: %d\n", i, traID);
		}
	}
	fclose(fp);
#endif
	return 0;
}


int SystemTest::MortonGridRangeQueryTestV2(MBB rangeQueryMBB, int rangeQueryNum)
{
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	MBB mbbArray[5000];
	for (int i = 0; i <= 4999; i++)
		rangeQueryMBB.randomGenerateMBB(mbbArray[i]);
		//mbbArray[i] = rangeQueryMBB;

	MyTimer timer;
	int* resultSize = NULL;
	vector<CPURangeQueryResult> resultTable;// bool-vector的vector
	resultTable.resize(rangeQueryNum);
	printf("single GPU Morton Grid range query #query=%d:\n", rangeQueryNum);

	MBB* tarray = new MBB[rangeQueryNum];
	vector<MBB> tvec;
	for (int j = 0; j < rangeQueryNum; j++) {
		tarray[j] = mbbArray[j];
		tvec.push_back(mbbArray[j]);
	}
	sort(tvec.begin(), tvec.end(), nonmem_cmp2);

	//for (int j = 0; j < rangeQueryNum; j++) {
	//	tarray[j].printMBB();
	//}
	//cout<<endl;
	//for (int j = 0; j < rangeQueryNum; j++) {
	//	tvec[j].printMBB();
	//}
	
	//CUDA_CALL(cudaSetDevice(0));

	vector<RangeQueryStateTable> stateTableRange; // RangeQueryStateTable结构体的vector
	stateTableRange.resize(rangeQueryNum * 50000);

	const int ONCE_QUERY_NUM = 50; 
	int queryIdx = 0;
	float timerstate = 0;
	for (queryIdx = 0; queryIdx < rangeQueryNum; queryIdx += ONCE_QUERY_NUM) {

#ifdef WIN32
		CUDA_CALL(cudaMalloc((void**)(&mgrid->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#else
		CUDA_CALL(cudaMalloc((void**)(&mgrid->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#endif
		CUDA_CALL(cudaMalloc((void**)&mgrid->stateTableGPU[0], (long long int)SMALL_MEM * 1024 * 1024));

		timer.start();
		void * allocatedGPUMem = mgrid->baseAddrRange[0];
		int querySize = (queryIdx + ONCE_QUERY_NUM < rangeQueryNum) ? ONCE_QUERY_NUM : rangeQueryNum - queryIdx;
		mgrid->baseAddrRange[0] = allocatedGPUMem;	
		mgrid->rangeQueryBatchGPU(&tarray[queryIdx], querySize, &resultTable[queryIdx], resultSize, &stateTableRange[queryIdx], 0);
		timer.stop();
		timerstate += timer.elapse();
		cout << "Single GPU Time of Morton Grid:" << timer.elapse() << "ms" << endl;

		CUDA_CALL(cudaFree(allocatedGPUMem));
		CUDA_CALL(cudaFree(mgrid->stateTableGPU[0]));
	}

	//timer.stop();
	//cout << "Single GPU Time of Morton Grid: GAT-R-noC " << timer.elapse() << "ms" << endl;

	cout << "Single GPU Time of Morton Grid GAT-R-noC: " << timerstate << "ms" << endl;




#ifdef USE_MULTIGPU
	printf("multi-GPU range query Morton Grid #query=%d:\n", rangeQueryNum);
	mgrid->rangeQueryBatchMultiGPU(mbbArray, rangeQueryNum, &resultTable[0], resultSize);
#else

#endif
#ifdef CHECK_CORRECT
	FILE* fp = fopen("MortonResult.txt", "w+");
	for (int i = 0; i <= rangeQueryNum - 1; i++)
	{
		for (int traID = 1; traID <= this->fsg->trajNum; traID++) {
			if (resultTable[i][traID])
				fprintf(fp, "Query %d result: %d\n", i, traID);
		}
	}
	fclose(fp);
#endif
	return 0;

}

int SystemTest::MortonGridRangeQueryTestV3(MBB rangeQueryMBB, int rangeQueryNum)
{
	this->rangeQueryMBB = rangeQueryMBB;
	this->rangeQueryNum = rangeQueryNum;
	MBB mbbArray[5000]; 
	for (int i = 0; i <= 4999; i++)
		rangeQueryMBB.randomGenerateMBB(mbbArray[i]); 
		//mbbArray[i] = rangeQueryMBB;
	MyTimer timer;
	// 查询结果指针
	int* resultSize = NULL;
	// 查询结果 resultTable
	vector<CPURangeQueryResult> resultTable;// bool-vector的vector
	resultTable.resize(rangeQueryNum);
	printf("single GPU Morton Grid range query #query=%d:\n", rangeQueryNum);
	MBB* tarray = new MBB[rangeQueryNum];
	vector<MBB> tvec;
	for (int j = 0; j < rangeQueryNum; j++) {
		tarray[j] = mbbArray[j];
		tvec.push_back(mbbArray[j]);

	}
	sort(tvec.begin(), tvec.end(), nonmem_cmp2);

	vector<RangeQueryStateTable> stateTableRange; // RangeQueryStateTable结构体的vector
	stateTableRange.resize(rangeQueryNum * 50000);

	const int ONCE_QUERY_NUM = 50;

	int queryIdx = 0;

	float timerstate = 0;

	for (queryIdx = 0; queryIdx < rangeQueryNum; queryIdx += ONCE_QUERY_NUM) {
#ifdef WIN32
		CUDA_CALL(cudaMalloc((void**)(&mgrid->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#else
		CUDA_CALL(cudaMalloc((void**)(&mgrid->baseAddrRange[0]), (long long int)BIG_MEM * 1024 * 1024));
#endif
		CUDA_CALL(cudaMalloc((void**)&mgrid->stateTableGPU[0], (long long int)SMALL_MEM * 1024 * 1024));


		timer.start();
		void * allocatedGPUMem = mgrid->baseAddrRange[0];

		int querySize = (queryIdx + ONCE_QUERY_NUM < rangeQueryNum) ? ONCE_QUERY_NUM : rangeQueryNum - queryIdx;
		mgrid->baseAddrRange[0] = allocatedGPUMem; 
		mgrid->rangeQueryBatchGPU(&tvec[queryIdx], querySize, &resultTable[queryIdx], resultSize, &stateTableRange[queryIdx], 0);
		timer.stop();
		timerstate += timer.elapse();
		cout << "Single GPU Time of Morton Grid:" << timer.elapse() << "ms" << endl;

		CUDA_CALL(cudaFree(allocatedGPUMem)); 
		CUDA_CALL(cudaFree(mgrid->stateTableGPU[0]));
	}

	cout << "Single GPU Time of Morton Grid: GAT-R " << timerstate << "ms" << endl;


#ifdef USE_MULTIGPU
	printf("multi-GPU range query Morton Grid #query=%d:\n", rangeQueryNum);
	mgrid->rangeQueryBatchMultiGPU(mbbArray, rangeQueryNum, &resultTable[0], resultSize);
#else

#endif
#ifdef CHECK_CORRECT
	FILE* fp = fopen("MortonResult.txt", "w+");
	for (int i = 0; i <= rangeQueryNum - 1; i++)
	{
		for (int traID = 1; traID <= this->fsg->trajNum; traID++) {
			if (resultTable[i][traID])
				fprintf(fp, "Query %d result: %d\n", i, traID);
		}
	}
	fclose(fp);
#endif
	return 0;
}









int SystemTest::similarityQueryTest(Trajectory t, int similarityScale, int similarityKValue)
{
	baseAddrGPU = NULL;
	Trajectory* qTra = new Trajectory[similarityScale];
	for (int j = 0; j <= similarityScale - 1; j++) {
		qTra[j] = t;
	}
	printf("qTra Length:%d qTra ID:%d\n", qTra[0].length, qTra[0].tid);
	int* simiResult;

	// Similarity on CPU
	// 多CPU版本
	/*
	simiResult = new int[similarityKValue * similarityScale];
	//printf("single-core CPU similarity @ k=%d and #query=%d:\n",similarityKValue,similarityScale);
	//g->SimilarityQueryBatch(qTra, similarityScale, simiResult, similarityKValue);
	printf("multi-core CPU similarity @ k=%d and #query=%d:\n",similarityKValue,similarityScale);
	g->SimilarityQueryBatchCPUParallel(qTra, similarityScale, simiResult, similarityKValue);

	delete[] simiResult;
	*/


	// Similarity on GPU
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("one GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	// grid 调用
	// 40 batch 大小 top-25 m=KSIMILARITY=80
	g->SimilarityQueryBatchOnGPU(qTra, similarityScale, simiResult, similarityKValue); // similarityScale没有分割
	delete[] simiResult;


#ifdef USE_MULTIGPU
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnMultiGPU(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult;
#else

#endif

	delete[] qTra;
	return 0;
}


int SystemTest::similarityQueryTest2(Trajectory* t, int similarityScale, int similarityKValue)
{

	baseAddrGPU = NULL; 
	Trajectory* qTra = new Trajectory[similarityScale];
	for (int j = 0; j <= similarityScale - 1; j++) {
		/*
		qTra[j].points.resize(this->g->cellBasedTrajectory[queryTrajNo].trajLength);
		// form query trajectories
		int cnt = 0;
		qTra[j].length = this->g->cellBasedTrajectory[queryTrajNo].trajLength;
		qTra[j].tid = queryTrajNo;

		for (int subID = 0; subID <= this->g->cellBasedTrajectory[queryTrajNo].length - 1; subID++)
		{
		int idxInAllPoints = this->g->cellBasedTrajectory[queryTrajNo].startIdx[subID];
		for (int pidx = 0; pidx <= this->g->cellBasedTrajectory[queryTrajNo].numOfPointInCell[subID] - 1; pidx++)
		{
		qTra[j].points[cnt + pidx].lat = this->g->allPoints[idxInAllPoints + pidx].y;
		qTra[j].points[cnt + pidx].lon = this->g->allPoints[idxInAllPoints + pidx].x;
		qTra[j].points[cnt + pidx].tid = this->g->allPoints[idxInAllPoints + pidx].tID;
		}
		cnt += this->g->cellBasedTrajectory[queryTrajNo].numOfPointInCell[subID];
		}
		*/
		qTra[j] = t[j];
	}
	int* simiResult;

	// Similarity on CPU
	// GAT-S-CPU
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
															  //printf("single-core CPU similarity @ k=%d and #query=%d:\n",similarityKValue,similarityScale);
															  //g->SimilarityQueryBatch(qTra, similarityScale, simiResult, similarityKValue);
	printf("GAT-S-noC GAT-S-CPU multi-core CPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchCPUParallel(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult;


	//Similarity on GPU
	// GAT-S-noE
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("GAT-S-noC GAT-S-noE one GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	// grid 调用
	// 40 batch 大小 top-25 m=KSIMILARITY=80 
	g->SimilarityQueryBatchOnGPUV3(qTra, similarityScale, simiResult, similarityKValue,0);
	delete[] simiResult; // 释放内存


	// GAT-S-E
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("GAT-S-noC GAT-S-E one GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	// grid 调用
	// 40 batch 大小 top-25 m=KSIMILARITY=80 
	g->SimilarityQueryBatchOnGPUV2(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult; // 释放内存





// 多GPU
#ifdef USE_MULTIGPU
						 /*
						 simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
						 printf("multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
						 g->SimilarityQueryBatchOnMultiGPU(qTra, similarityScale, simiResult, similarityKValue);
						 delete[] simiResult; // 释放内存
						 */
#else

#endif

	delete[] qTra;
	return 0;


}


int SystemTest::similarityQueryTest3(vector<Trajectory> &t, int similarityScale, int similarityKValue)
{
	//baseAddrGPU = NULL;
	Trajectory* qTra = new Trajectory[similarityScale]; // 待查询的轨迹 Tq
	Trajectory* qTra2 = new Trajectory[similarityScale]; // 待查询的轨迹 Tq
	for (int j = 0; j <= similarityScale - 1; j++) {
		/*
		qTra[j].points.resize(this->g->cellBasedTrajectory[queryTrajNo].trajLength);
		// form query trajectories
		int cnt = 0;
		qTra[j].length = this->g->cellBasedTrajectory[queryTrajNo].trajLength;
		qTra[j].tid = queryTrajNo;

		for (int subID = 0; subID <= this->g->cellBasedTrajectory[queryTrajNo].length - 1; subID++)
		{
		int idxInAllPoints = this->g->cellBasedTrajectory[queryTrajNo].startIdx[subID];
		for (int pidx = 0; pidx <= this->g->cellBasedTrajectory[queryTrajNo].numOfPointInCell[subID] - 1; pidx++)
		{
		qTra[j].points[cnt + pidx].lat = this->g->allPoints[idxInAllPoints + pidx].y;
		qTra[j].points[cnt + pidx].lon = this->g->allPoints[idxInAllPoints + pidx].x;
		qTra[j].points[cnt + pidx].tid = this->g->allPoints[idxInAllPoints + pidx].tID;
		}
		cnt += this->g->cellBasedTrajectory[queryTrajNo].numOfPointInCell[subID];
		}
		*/
		qTra[j] = t[j];
		qTra2[j] = t[j];
	}

	ReorderArray(qTra2, similarityScale);

	// printf("qTra Length:%d qTra ID:%d\n", qTra[0].length, qTra[0].tid);

	int* simiResult;

	
	//******** Similarity on CPU 单线程CPU ********
	simiResult = new int[similarityKValue * similarityScale];
	printf("single-core CPU similarity @ k=%d and #query=%d:\n",similarityKValue,similarityScale);
	g->SimilarityQueryBatch(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult;
	

	//******** Similarity on CPU 多线程CPU 比较 ********
	simiResult = new int[similarityKValue * similarityScale];
	//printf("single-core CPU similarity @ k=%d and #query=%d:\n",similarityKValue,similarityScale);
	//g->SimilarityQueryBatch(qTra, similarityScale, simiResult, similarityKValue);
	printf("GAT-S-C GAT-S-CPU multi-core CPU similarity @ k=%d and #query=%d:\n",similarityKValue,similarityScale);
	g->SimilarityQueryBatchCPUParallel(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult;
	

	// fig. 6
	//******** Similarity on GPU 单GPU 比较 ********
	// GAT-S-noE

	printf("One GPU noD + noE\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf(" one GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnGPUV3(qTra, similarityScale, simiResult, similarityKValue,0);
	

	printf("One GPU noMAT noD + noE\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf(" one GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnGPUNoMAT(qTra, similarityScale, simiResult, similarityKValue, 0);


	/*
	// checking result whether need asysc？ donnot know yet??
	for(int qid2 = 0; qid2 <similarityScale;qid2++){
		sort(simiResult+qid2*similarityKValue,simiResult+(qid2+1)*similarityKValue);
	}
	
	for (int qID = 0; qID <= similarityScale - 1; qID++)
	{
		for (int i = 0; i <= similarityKValue - 1; i++)
		{
			cout<<simiResult[qID * similarityKValue + i]<<' ';
		}
		cout << endl;
	}
	cout<<endl;
	*/

	delete[] simiResult; // 释放内存
	

	/*
	printf("D + noE\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf(" one GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnGPUV3(qTra2, similarityScale, simiResult, similarityKValue, 0);
	delete[] simiResult; // 释放内存
	*/


	// 完全没有意义 也就论文扯一扯 不要去掉 没有任何效果提升
	// 20*40 = 800 bloack 
	// 256 thread
	// 完全可以充分利用GPU 
	// this is okay !!
	/*
	//Similarity on GPU
	// GAT-S-E
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("GAT-S-C GAT-S-E one GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	// grid 调用
	g->SimilarityQueryBatchOnGPUV2(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult; // 释放内存
	*/

#ifdef USE_MULTIGPU


	printf("FineGrained noE+noD\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("FineGrained multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnMultiGPU(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult; // 释放内存

	printf("FineGrained noMAT noE+noD\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("FineGrained multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnMultiGPUNoMAT(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult; // 释放内存


	printf("FineGrained noE+D\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("FineGrained multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnMultiGPU(qTra2, similarityScale, simiResult, similarityKValue);
	delete[] simiResult; // 释放内存

	printf("FineGrained noMAT noE+D\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("FineGrained multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnMultiGPUNoMAT(qTra2, similarityScale, simiResult, similarityKValue);
	delete[] simiResult; // 释放内存

	

	printf("CoarseGrainedThreadway noD\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("CoarseGrainedThreadway multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnMultiGPUV2(qTra, similarityScale, simiResult, similarityKValue);
	delete[] simiResult; // 释放内存


	printf("CoarseGrainedThreadway D\n");
	simiResult = new int[similarityKValue * similarityScale]; // 保存轨迹ID
	printf("CoarseGrainedThreadway multi-GPU similarity @ k=%d and #query=%d:\n", similarityKValue, similarityScale);
	g->SimilarityQueryBatchOnMultiGPUV2(qTra2, similarityScale, simiResult, similarityKValue);
	delete[] simiResult; // 释放内存
	
#else

#endif


	delete[] qTra;
	return 0;
}

