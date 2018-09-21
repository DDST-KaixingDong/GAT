#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cudaKernel.h"
#include "thrust/device_ptr.h"
#include "thrust/remove.h"
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include "ConstDefine.h"


#define CUDA_CALL(x) { const cudaError_t a = (x); if (a!= cudaSuccess) { printf("\nCUDA Error: %s(err_num=%d)\n", cudaGetErrorString(a), a); cudaDeviceReset(); assert(0);}}

/*
���м���1�����ģdp
��Ҫ��ǰ����ǰ����dp�Ľ���������ڹ����ڴ���
iter: �ڼ���dp��λ��outputIdx����������ȫ���ڴ�λ�ã�tra1��tra2�������켣����ǰ�����빲���ڴ棻
*/
//__global__ void DPforward(const int iter, const int* outputIdx,const SPoint *tra1,const SPoint *tra2) {
//	SPoint p1 = tra1[threadIdx.x];
//	SPoint p2 = tra2[iter - threadIdx.x - 1]; //�������ڴ��Ǿۼ����ʵ���
//	bool subcost;
//	if((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON)) {
//		subcost = 0;
//	}
//	else
//		subcost = 1;
//
//}

/*
SPoint�汾
case1���켣����С��512
���м���n��DP
��Ҫ��ǰ����ǰ����dp�Ľ���������ڹ����ڴ���
queryTra[],candidateTra[][]:�켣
stateTableGPU[][]:��ÿ��candidate��state��
result[]:����ÿ��candidate��EDR���
�Ż�����
1���켣����share memory����
2��ֱ�Ӵ��ݹ켣����ʹ��ָ��
*/
__global__ void EDRDistance_1(SPoint *queryTra, SPoint **candidateTra, int candidateNum, int queryLength, int *candidateLength, int** stateTableGPU, int *result) {
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	if (blockID >= candidateNum) return;
	if ((threadID >= candidateLength[blockID]) && (threadID >= queryLength)) return;
	const int lenT = candidateLength[blockID];
	//int iterNum = queryLength;
	//if (lenT > queryLength)
	//	iterNum = lenT;
	const int iterNum = queryLength + lenT - 1;
	__shared__ short state[2][MAXTHREAD]; //���ڴ洢ǰ���εĽ��
	state[0][0] = 0;
	state[1][0] = 1;
	state[1][1] = 1;
	//�������켣���򣬱�֤��һ���ȵڶ�����
	//���Ȱѹ켣���ڹ����ڴ���
	__shared__ SPoint queryTraS[MAXTHREAD];
	__shared__ SPoint traData[MAXTHREAD];
	if (threadID < lenT) {
		traData[threadID] = candidateTra[blockID][threadID];
	}
	if (threadID < queryLength) {
		queryTraS[threadID] = queryTra[threadID];
	}
	const SPoint *tra1, *tra2; //��֤tra1��tra2��
	int len1, len2;
	if (lenT >= queryLength) {
		tra1 = queryTraS;
		tra2 = traData;
		len1 = queryLength;
		len2 = lenT;
	}
	else
	{
		tra1 = traData;
		tra2 = queryTraS;
		len1 = lenT;
		len2 = queryLength;
	}

	int myState;
	for (int i = 0; i <= iterNum - 1; i++) {//��i��dp
		if (i < len1 - 1) {
			if (threadID <= i) {
				SPoint p1 = tra1[threadID];
				SPoint p2 = tra2[i - threadID]; //�������ڴ��Ǿۼ����ʵ���
				bool subcost;
				//if((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON)) {
				//	subcost = 0;
				//}
				//else
				//	subcost = 1;
				subcost = !((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON));
				int state_ismatch = state[0][threadID] + subcost;
				int state_up = state[1][threadID] + 1;
				int state_left = state[1][threadID + 1] + 1;
				if (state_ismatch < state_up)
					myState = state_ismatch;
				else if (state_left < state_up)
					myState = state_left;
				else
					myState = state_up;
				//ȥ��if�ı�﷽ʽ���Ƿ�����������ܣ�
				//myState = (state_ismatch < state_up) * state_ismatch + (state_left < state_up) * state_up + (state_left >= state_up) * state_left;

			}
		}
		else if (i > iterNum - len1) {
			if (threadID <= iterNum - i - 1) {
				SPoint p1 = tra1[threadID + len1 - (iterNum - i)];
				SPoint p2 = tra2[len2 - 1 - threadID]; //�������ڴ��Ǿۼ����ʵ���
				bool subcost;
				if ((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON)) {
					subcost = 0;
				}
				else
					subcost = 1;
				int state_ismatch = state[0][threadID + 1] + subcost;
				int state_up = state[1][threadID] + 1;
				int state_left = state[1][threadID + 1] + 1;
				if (state_ismatch < state_up)
					myState = state_ismatch;
				else if (state_left < state_up)
					myState = state_left;
				else
					myState = state_up;
			}
		}
		else
		{
			if (threadID < len1) {
				SPoint p1 = tra1[threadID];
				SPoint p2 = tra2[i - threadID]; //�������ڴ��Ǿۼ����ʵ���
				bool subcost;
				if ((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON)) {
					subcost = 0;
				}
				else
					subcost = 1;
				int state_ismatch = state[0][threadID] + subcost;
				int state_up = state[1][threadID] + 1;
				int state_left = state[1][threadID + 1] + 1;
				if (state_ismatch < state_up)
					myState = state_ismatch;
				else if (state_left < state_up)
					myState = state_left;
				else
					myState = state_up;
			}
		}
		//дmyState��share�ڴ�,ckecked
		int startidx;
		//���Ƚ�������д��ȫ���ڴ棬ȫд
		//startidx�Ǿɵ�����Ӧ����ȫ���ڴ��е�ַ����i-2����
		//����Ӧд��ȫ���ڴ����ʼλ��

		if (i - 2 < len1 - 2) {
			startidx = (i - 2 + 2)*(i - 2 + 3) / 2;
			if (threadID <= i) {
				stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
			}
		}
		else if (i - 2 >= iterNum - len1) {
			startidx = (len1 + 1)*(len2 + 1) - (iterNum - (i - 2))*(iterNum - (i - 2) + 1) / 2;
			if (threadID <= iterNum - i + 1) {
				stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
			}
		}
		else
		{
			startidx = (len1 + 1)*((i - 2) - (len1 - 2)) + len1*(len1 + 1) / 2;
			if (threadID <= len1) {
				stateTableGPU[blockID][threadID + startidx] = state[0][threadID];
			}
		}

		//�ƶ������ݵ�������
		state[0][threadID] = state[1][threadID];
		//д��������
		if (i < len1 - 1) {
			if (threadID <= i)
				state[1][threadID + 1] = myState;
			if (threadID == 0) {
				state[1][0] = i + 2;
				state[1][i + 2] = i + 2;
			}
		}
		else if (i >= iterNum - len1) {
			if (threadID <= iterNum - i - 1)
				state[1][threadID] = myState;
		}
		else
		{
			if (threadID < len1)
				state[1][threadID + 1] = myState;
			if (threadID == 0) {
				state[1][0] = i + 2;
			}
		}
		__syncthreads();
	}
	//�����������һ�μ���һ�����ɽ���0��ɵ�
	if (threadID == 0)
		result[blockID] = myState;
}


//__global__ void testSharedMemory()
//{
//	__shared__ SPoint queryTraS[MAXLENGTH];
//	__shared__ SPoint traData[MAXLENGTH];
//	__shared__ SPoint traData2[MAXLENGTH];
//	SPoint s;
//	s.x = 4;
//	s.y = 5;
//	traData[1535] = s;
//	queryTraS[1535] = s;
//	traData2[1535] = s;
//}

/*
SPoint�汾
ͬʱ�������ɸ�query��EDR��������һ��EDR����Ϊ��λ��ÿ��block����һ��EDR��thread����һ��б����state�Ĳ��м��㡣
case1���켣���ȿɳ���512������ѭ���������512�ġ�
���м���n��DP
��Ҫ��ǰ����ǰ����dp�Ľ���������ڹ����ڴ���
queryTaskNum:�ܹ��м���EDR��������
queryTaskInfo[]��ÿ��task��Ӧ��qID��candidateID��Ϣ����struct�洢
queryTra[],candidateTra[]:�켣���ݣ�candidateTra��֤���ڲ��켣���ظ�
queryTraOffset[],candidateTraOffset[]:ÿ���켣��offset��candidateTra��֤���ڲ��켣���ظ�
queryLength[],candidateLength[]:ÿ���켣�ĳ��ȣ���ʵoffset������ǳ��ȣ�����idx������Ķ�Ӧ
����candidateLength[id]�ǵ�id��candidate Traj�ĳ���
stateTableGPU[][]:��ÿ��candidate��state��
result[]:����ÿ��candidate��EDR���
�Ż�����
1���켣����share memory����
2��ֱ�Ӵ��ݹ켣����ʹ��ָ��
*/

// EDRDistance_Batch_Handler(validCandTrajNum, taskInfoTableGPU, queryTraGPUBase, queryTraOffsetGPU, candidateOffsetsGPU, queryLengthGPU, candidateTraLengthGPU, resultReturnedGPU, &defaultStream);
__global__ void EDRDistance_Batch(int queryTaskNum, TaskInfoTableForSimilarity* taskInfoTable, SPoint *queryTra, int* queryTraOffset, SPoint** candidateTraOffsets, int* queryLength, int *candidateLength, int *result) {

	int blockID = blockIdx.x;
	int threadID = threadIdx.x;

	if (blockID >= queryTaskNum) return;

	__shared__ int thisQueryID;
	__shared__ int thisQueryLength; 
	__shared__ int lenT;
	__shared__ int iterNum;

	thisQueryID = taskInfoTable[blockID].qID;
	if(thisQueryID < 0 ) return; // block���ֿ���

	__shared__ short state[2][MAXLENGTH + 1];

	__shared__ SPoint *queryTraS;
	__shared__ SPoint *traData;
	__shared__ SPoint *tra1, *tra2;
	__shared__ int len1, len2;

	if (threadID == 0) {
		thisQueryID = taskInfoTable[blockID].qID;
		thisQueryLength = queryLength[thisQueryID];
		lenT = candidateLength[blockID];
		iterNum = thisQueryLength + lenT - 1;
		state[0][0] = 0;
		// state[0][1] = 1;
		state[1][0] = 1;
		state[1][1] = 1;
		queryTraS = queryTra + queryTraOffset[thisQueryID]; 
		traData = candidateTraOffsets[blockID]; 
		if (lenT >= thisQueryLength) {
			tra1 = queryTraS;
			tra2 = traData;
			len1 = thisQueryLength;
			len2 = lenT;
		}
		else
		{
			tra1 = traData;
			tra2 = queryTraS;
			len1 = lenT;
			len2 = thisQueryLength;
		}
	}

	__syncthreads(); // ͬ����֤�����߳̿���

	if ((threadID >= lenT) && (threadID >= thisQueryLength)) return;	
	//__shared__ SPoint queryTraS[MAXLENGTH];
	//__shared__ SPoint traData[MAXLENGTH];

	//for (int i = 0; i <= lenT - 1;i+=MAXTHREAD)
	//{
	//	if(threadID+i<lenT)
	//	{
	//		traData[threadID + i] = SPoint(candidateTraOffsets[blockID][threadID + i]);
	//	}
	//}

	//SPoint* queryTraBaseAddr = queryTra + queryTraOffset[thisQueryID];
	//for (int i = 0; i <= thisQueryLength - 1;i+=MAXTHREAD)
	//{
	//	if(threadID+i<thisQueryLength)
	//	{
	//		queryTraS[threadID + i] = *(queryTraBaseAddr + threadID + i);
	//	}
	//}
	
	int myState[5];// 256*4 = 1024 �� __shared ����
	int nodeID;
	SPoint p1;
	SPoint p2;
	bool subcost;

	for (int i = 0; i <= iterNum - 1; i++) { // block ���ѭ��
		if (i < len1 - 1) {
			for (int startIdx = 0; startIdx <= i; startIdx += MAXTHREAD) {
				nodeID = startIdx + threadID; 
				if (nodeID <= i) {

					p1 = tra1[nodeID]; // fetch from global memory
					p2 = tra2[i - nodeID]; // fetch from global memory

					// id1 + id2 = i
					if ((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON)) {
						subcost = 0;
					}
					else
						subcost = 1;
					//subcost = !((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON));
					bool c1 = ((state[0][nodeID] + subcost < (state[1][nodeID] + 1)) && (state[0][nodeID] + subcost < (state[1][nodeID + 1] + 1)));
					bool c2 = (((state[1][nodeID + 1] + 1) < (state[1][nodeID] + 1)) && (((state[1][nodeID + 1] + 1) < state[0][nodeID] + subcost)));
					//ȥ��if�ı�﷽ʽ���Ƿ�����������ܣ�
					myState[nodeID / MAXTHREAD] = c1 * (state[0][nodeID] + subcost) + c2 * (state[1][nodeID + 1] + 1) + !(c1 || c2) * (state[1][nodeID] + 1);
					//if ((state_ismatch < state_up) && (state_ismatch < state_left))
					//	myState[nodeID/MAXTHREAD] = state_ismatch;
					//else if ((state_left < state_up) && ((state_left < state_ismatch)))
					//	myState[nodeID / MAXTHREAD] = state_left;
					//else
					//	myState[nodeID / MAXTHREAD] = state_up;
					////ȥ��if�ı�﷽ʽ���Ƿ�����������ܣ�
					//myState[nodeID / MAXTHREAD] = (state_ismatch < state_up) && (state_ismatch < state_left) * state_ismatch + ((state_left < state_up) && ((state_left < state_ismatch))) * state_left + !(((state_ismatch < state_up) && (state_ismatch < state_left))||(((state_left < state_up) && ((state_left < state_ismatch))))) * state_up;
				}
			}
		}
		else if (i > iterNum - len1) {
			for (int startIdx = 0; startIdx <= iterNum - i - 1; startIdx += MAXTHREAD) {
				nodeID = startIdx + threadID;
				if (nodeID <= iterNum - i - 1) {
					// EDR�����������
					p1 = tra1[nodeID + len1 - (iterNum - i)]; // fetch from global memory
					p2 = tra2[len2 - 1 - nodeID];
					if ((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON)) {
						subcost = 0;
					}
					else
						subcost = 1;
					//if (state_ismatch < state_up)
					//	myState[nodeID / MAXTHREAD] = state_ismatch;
					//else if (state_left < state_up)
					//	myState[nodeID / MAXTHREAD] = state_left;
					//else
					//	myState[nodeID / MAXTHREAD] = state_up;
					bool c1 = (((state[0][nodeID + 1] + subcost) < (state[1][nodeID] + 1)) && ((state[0][nodeID + 1] + subcost) < (state[1][nodeID + 1] + 1)));
					bool c2 = (((state[1][nodeID + 1] + 1) < (state[1][nodeID] + 1)) && (((state[1][nodeID + 1] + 1) < (state[0][nodeID + 1] + subcost))));
					//ȥ��if�ı�﷽ʽ���Ƿ�����������ܣ�
					myState[nodeID / MAXTHREAD] = c1 * (state[0][nodeID + 1] + subcost) + c2 * (state[1][nodeID + 1] + 1) + !(c1 || c2) * (state[1][nodeID] + 1);
				}
			}
		}
		else
		{
			for (int startIdx = 0; startIdx < len1; startIdx += MAXTHREAD) {
				nodeID = startIdx + threadID;
				if (nodeID < len1) { // ע�������ж� ��֤������̸߳������ø��Ǽ��� �����߳̿���
					p1 = tra1[nodeID]; // fetch from global memory
					p2 = tra2[i - nodeID]; //�������ڴ��Ǿۼ����ʵ���?
					if ((fabs(p1.x - p2.x) < EPSILON) && (fabs(p1.y - p2.y)<EPSILON)) {
						subcost = 0;
					}
					else
						subcost = 1;
					//int state_ismatch = (state[0][nodeID] + subcost);
					//int state_up = (state[1][nodeID] + 1);
					//int state_left = (state[1][nodeID + 1] + 1);
					//if (state_ismatch < state_up)
					//	myState[nodeID / MAXTHREAD] = state_ismatch;
					//else if (state_left < state_up)
					//	myState[nodeID / MAXTHREAD] = state_left;
					//else
					//	myState[nodeID / MAXTHREAD] = state_up;
					bool c1 = (((state[0][nodeID] + subcost) < (state[1][nodeID] + 1)) && ((state[0][nodeID] + subcost) < (state[1][nodeID + 1] + 1)));
					bool c2 = (((state[1][nodeID + 1] + 1) < (state[1][nodeID] + 1)) && (((state[1][nodeID + 1] + 1) < (state[0][nodeID] + subcost))));
					//ȥ��if�ı�﷽ʽ���Ƿ�����������ܣ�
					myState[nodeID / MAXTHREAD] = c1 * (state[0][nodeID] + subcost) + c2 * (state[1][nodeID + 1] + 1) + !(c1 || c2) * (state[1][nodeID] + 1);
				}
			}
		}

		// ����state[0]
		//state[1] �� state[0]
		for (int Idx = 0; Idx < MAXLENGTH; Idx += MAXTHREAD)
		{
			if(threadID + Idx < MAXLENGTH)
				state[0][threadID + Idx] = state[1][threadID + Idx];
		}
		//state[0][threadID] = state[1][threadID];

		// ����state[1]
		//д��������
		if (i < len1 - 1) {
			for (int Idx = 0; Idx <= i; Idx += MAXTHREAD) {
				if (threadID + Idx <= i)
					state[1][Idx + threadID + 1] = myState[Idx / MAXTHREAD];
			}
			// �����׶�
			if (threadID == 0) {
				state[1][0] = i + 2; // ����ͷ
				state[1][i + 2] = i + 2; // ����β
			}
		}
		else if (i >= iterNum - len1) {
			//if (threadID <= iterNum - i - 1)
			//	state[1][threadID] = myState;
			for (int Idx = 0; Idx <= iterNum - i - 1; Idx += MAXTHREAD) {
				if (threadID + Idx <= iterNum - i - 1)
					state[1][threadID + Idx] = myState[Idx / MAXTHREAD];
			}
			// �����׶β��ø���
		}
		else
		{
			//if (threadID < len1)
			//	state[1][threadID + 1] = myState;
			//if (threadID == 0) {
			//	state[1][0] = i + 2;
			//}
			for (int Idx = 0; Idx <= len1; Idx += MAXTHREAD) {
				if (threadID + Idx < len1)
					state[1][Idx + threadID + 1] = myState[Idx / MAXTHREAD];
			}
			// ���ֽ׶�
			if (threadID == 0) {
				state[1][0] = i + 2; // ֻ�����ͷ 
			}
		}

		__syncthreads(); // ͬ��һ��block��thread

	}
	// kernel Ӧ��ûʲô������
	// for ѭ������ �����EDR��� 
	// myState[0]; ����EDR
	if (threadID == 0 && blockID < queryTaskNum)
		result[blockID] = myState[0];
	// std::cout << "calc EDR success!\n" ;
}
// EDRDistance_Batch_Handler(validCandTrajNum, taskInfoTableGPU, queryTraGPUBase, queryTraOffsetGPU, candidateOffsetsGPU, queryLengthGPU, candidateTraLengthGPU, resultReturnedGPU, &defaultStream);
int EDRDistance_Batch_Handler(int queryTaskNum, TaskInfoTableForSimilarity* taskInfoTable, SPoint *queryTra, int* queryTraOffset, SPoint** candidateTraOffsets, int* queryLength, int *candidateLength, int *result, cudaStream_t *stream)
{
	//printf("run kernel now\n");
	EDRDistance_Batch <<< queryTaskNum, MAXTHREAD, 0, *stream >> >(queryTaskNum, taskInfoTable, queryTra, queryTraOffset, candidateTraOffsets, queryLength, candidateLength, result);
	//CUDA_CALL(cudaDeviceSynchronize());
	return 0;
}

__device__ inline int binary_search_intPair(intPair* temp, int left, int right, int val)
{
	int mid = (left + right) / 2;
	while (left <= right)
	{
		mid = (left + right) / 2;
		if (temp[mid].int_1 == val)
			return temp[mid].int_2;
		else if (temp[mid].int_1 > val)
		{
			right = mid - 1;
		}
		else
			left = mid + 1;
	}
	return 0;
}

__device__ inline int binary_search_intPair_Neighbor(intPair* temp, int left, int right, int val)
{
	int mid = (left + right) / 2;
	while (left <= right)
	{
		mid = (left + right) / 2;
		if (temp[mid].int_1 == val)
			return mid;
		else if (temp[mid].int_1 > val)
		{
			right = mid - 1;
		}
		else
			left = mid + 1;
	}
	return -1;
}

// -1Ϊû�ҵ�
__device__ inline int binary_search_int(int* temp, int left, int right, int val)
{
	int mid = (left + right) / 2;
	while (left <= right)
	{
		mid = (left + right) / 2;
		if (temp[mid] == val)
			return mid;
		else if (temp[mid] > val)
		{
			right = mid - 1;
		}
		else
			left = mid + 1;
	}
	return -1;
}

__device__ inline int getIdxFromXYGPU(int x, int y)
{
	int lenx, leny;
	if (x == 0)
		lenx = 1;
	else
	{
		lenx = int(log2f(x)) + 1;
	}
	if (y == 0)
		leny = 1;
	else
		leny = int(log2f(y)) + 1;
	int result = 0;
	int xbit = 1, ybit = 1;
	for (int i = 1; i <= 2 * max(lenx, leny); i++)
	{
		if ((i & 1) == 1) //����
		{
			result += (x >> (xbit - 1) & 1) * (1 << (i - 1));
			xbit = xbit + 1;
		}
		else //ż��
		{
			result += (y >> (ybit - 1) & 1) * (1 << (i - 1));
			ybit = ybit + 1;
		}
	}
	return result;
}

__device__ inline int findNeighborGPU(int cellNum, int cellID, int * neighborID)
{
	int x = 0, y = 0;
	for (int bit = 0; bit <= int(log2f(cellNum)) - 1; bit++) {
		if (bit % 2 == 0) {
			//����λ
			x += ((cellID >> bit)&(1))*(1 << (bit / 2));
		}
		else {
			//ż��λ
			y += ((cellID >> bit)&(1))*(1 << (bit / 2));
		}
	}
	int cnt = 0;
	for (int xx = x - 1; xx <= x + 1; xx++) {
		for (int yy = y - 1; yy <= y + 1; yy++) {
			if ((xx != x) || (yy != y))
				neighborID[cnt++] = getIdxFromXYGPU(xx, yy);
			//printf("%d\t", cnt);
		}
	}
	return 0;
}

__device__ inline bool isPositive(short x)
{
	return x >= 0;
}

__global__ void Calculate_FD_Sparse(intPair* queryFVGPU, intPair* FVinfo, intPair* FVTable, intPair* SubbedArray, intPair* SubbedArrayOffset, int SubbedArrayJump, int queryCellLength, int startTrajIdx, int checkNum, int cellNum, int trajNumInDB, int nonZeroFVNumInDB, short* FDistance)
{
	//��һ�׶Σ����м���
	const int MAX_QUERY_CELLNUMBER = 512;
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	int threadIDGlobal = blockDim.x*blockID + threadID;

	__shared__ intPair queryCellTraj[MAX_QUERY_CELLNUMBER];
	__shared__ intPair dbCellTraj[MAX_QUERY_CELLNUMBER];
	//cellchecked��¼��query�г��ֵ�cell��ţ������ڷ��������ʱ�����ǲ����Ѿ������ˡ��Ժ�����ڹ鲢���и��ô˱�����
	__shared__ int cellChecked[MAX_QUERY_CELLNUMBER];
	for (int i = 0; i <= queryCellLength - 1; i += MAXTHREAD) {
		if (threadID + i < queryCellLength)
		{
			queryCellTraj[threadID + i] = queryFVGPU[threadID + i];
		}
	}
	int dbTrajStartIdx = FVinfo[startTrajIdx + blockID].int_2;
	int dbTrajEndIdx;
	if (blockID + startTrajIdx == trajNumInDB - 1)
		dbTrajEndIdx = nonZeroFVNumInDB - 1;
	else
		dbTrajEndIdx = FVinfo[startTrajIdx + blockID + 1].int_2 - 1;

	for (int i = 0; i <= dbTrajEndIdx - dbTrajStartIdx; i += MAXTHREAD)
	{
		if (threadID + i <= dbTrajEndIdx - dbTrajStartIdx)
			dbCellTraj[threadID + i] = FVTable[dbTrajStartIdx + threadID + i];
	}
	//1.1:��query��ȥdb
	for (int i = 0; i < queryCellLength; i += MAXTHREAD)
	{
		if (threadID + i < queryCellLength) {
			int find = binary_search_intPair(dbCellTraj, 0, dbTrajEndIdx - dbTrajStartIdx, queryCellTraj[threadID + i].int_1);
			cellChecked[threadID + i] = queryCellTraj[threadID + i].int_1;
			SubbedArray[SubbedArrayJump * blockID + threadID + i].int_1 = queryCellTraj[threadID + i].int_1;
			SubbedArray[SubbedArrayJump * blockID + threadID + i].int_2 = queryCellTraj[threadID + i].int_2 - find;
		}
		if (threadID == 0) {
			SubbedArrayOffset[blockID].int_1 = queryCellLength - 1;
			SubbedArrayOffset[blockID].int_2 = queryCellLength + dbTrajEndIdx - dbTrajStartIdx;
		}
	}
	//1.2����db��ȥquery��ע��Ӹ���
	for (int i = 0; i <= dbTrajEndIdx - dbTrajStartIdx; i += MAXTHREAD)
	{
		if (threadID + i <= dbTrajEndIdx - dbTrajStartIdx)
		{
			intPair cellNo = dbCellTraj[threadID + i];
			int find = binary_search_int(cellChecked, 0, queryCellLength - 1, cellNo.int_1);
			if (find == -1)
			{
				SubbedArray[SubbedArrayJump * blockID + queryCellLength + threadID + i].int_1 = cellNo.int_1;
				SubbedArray[SubbedArrayJump * blockID + queryCellLength + threadID + i].int_2 = -cellNo.int_2;
			}
			else
				SubbedArray[SubbedArrayJump * blockID + queryCellLength + threadID + i].int_1 = -1;
		}
	}
	__syncthreads();
	//�ڶ��׶Σ��������ڣ�������
	//����׶θ�Ϊÿ��thread����һ��FD
	//2.1���ϲ�ÿ��subbedArray
	if (threadIDGlobal < checkNum) {
		int startMergeIdx = SubbedArrayOffset[threadIDGlobal].int_1 + 1;
		int endMergeIdx = SubbedArrayOffset[threadIDGlobal].int_2;
		int frontPtr = startMergeIdx;
		for (int i = startMergeIdx; i <= endMergeIdx; i++)
		{
			if (SubbedArray[SubbedArrayJump * threadIDGlobal + i].int_1 != -1)
			{
				SubbedArray[SubbedArrayJump * threadIDGlobal + frontPtr] = SubbedArray[SubbedArrayJump * threadIDGlobal + i];
				frontPtr++;
			}
		}
		SubbedArrayOffset[threadIDGlobal].int_2 = frontPtr - 1;
	}
	//2.2 ��������
	int neighborsID[8];
	//cell����ָ�ڼ���Ԫ��
	for (int cell = 0; cell <= SubbedArrayOffset[threadIDGlobal].int_2; cell++)
	{
		findNeighborGPU(cellNum, cell, neighborsID);
		//for (int i = 0; i <= 7; i++)
		//	neighborsID[i] = 11;
		for (int i = 0; i <= 7; i++)
		{
			int find = binary_search_intPair_Neighbor(&SubbedArray[SubbedArrayJump * threadIDGlobal], 0, SubbedArrayOffset[threadIDGlobal].int_1, neighborsID[i]);
			if (find == -1) {
				find = binary_search_intPair_Neighbor(&SubbedArray[SubbedArrayJump * threadIDGlobal], SubbedArrayOffset[threadIDGlobal].int_1 + 1, SubbedArrayOffset[threadIDGlobal].int_2, neighborsID[i]);
			}
			// �����-1��˵�����neighbor��0�����ô���
			if (find != -1)
			{
				if (isPositive(SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2) != isPositive(SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2))
				{
					if (fabsf(SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2) > fabsf(SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2))
					{
						SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2 = SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2 + SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2;
						SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2 = 0;
					}
					else {
						SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2 = SubbedArray[SubbedArrayJump * threadIDGlobal + find].int_2 + SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2;
						SubbedArray[SubbedArrayJump * threadIDGlobal + cell].int_2 = 0;
						break;
					}
				}
			}
		}
	}
	__syncthreads();
	//�����׶Σ�ͳ����������
	//��Ȼ��ÿ��block����һ��FD�ļ���
	if (blockID >= checkNum)
		return;
	int *tempsumPosi = (int*)queryCellTraj;
	int *tempsumNega = (int*)dbCellTraj;
	tempsumPosi[threadID] = 0;
	tempsumNega[threadID] = 0;
	for (int i = 0; i <= SubbedArrayOffset[blockID].int_2; i += MAXTHREAD)
	{
		if (i + threadID <= SubbedArrayOffset[blockID].int_2)
		{
			tempsumPosi[threadID] += (isPositive(SubbedArray[SubbedArrayJump * blockID + i + threadID].int_2)*SubbedArray[SubbedArrayJump * blockID + i + threadID].int_2);
			tempsumNega[threadID] += (-(!isPositive(SubbedArray[SubbedArrayJump * blockID + i + threadID].int_2))*SubbedArray[SubbedArrayJump * blockID + i + threadID].int_2);
		}
	}
	__shared__ int sizeOfTempSum;
	if (threadID == 0)
		sizeOfTempSum = MAXTHREAD;
	__syncthreads();
	while ((sizeOfTempSum>1))
	{
		if (threadID <= (sizeOfTempSum >> 1) - 1)
		{
			tempsumPosi[threadID] = tempsumPosi[threadID] + tempsumPosi[threadID + (sizeOfTempSum >> 1)];
			tempsumNega[threadID] = tempsumNega[threadID] + tempsumNega[threadID + (sizeOfTempSum >> 1)];
		}
		__syncthreads();
		if (threadID == 0)
			sizeOfTempSum = (sizeOfTempSum >> 1);
		__syncthreads();
	}
	if (threadID == 0)
		FDistance[blockID] = (tempsumPosi[0] > tempsumNega[0]) ? tempsumPosi[0] : tempsumNega[0];
}

//ÿ��block����һ��FD�ļ���
__global__ void Calculate_FD_NonColumn(short* queryFVGPU, intPair* FVinfo, intPair* FVTable, int startTrajIdx, int checkNum, int cellNum, int trajNumInDB, int nonZeroFVNumInDB, short* FDistance)
{
	//��һ�׶Σ����м���
	int blockID = blockIdx.x;
	int threadID = threadIdx.x;
	int threadIDGlobal = blockDim.x*blockID + threadID;
	if (blockID >= checkNum)
		return;
	__shared__ intPair taskInfo;
	if (threadID == 0)
		taskInfo = FVinfo[blockID + startTrajIdx];
	int nextCnt;
	if (blockID + startTrajIdx == trajNumInDB - 1)
		nextCnt = nonZeroFVNumInDB;
	else
		nextCnt = FVinfo[blockID + startTrajIdx + 1].int_2;
	__syncthreads();
	for (int i = 0; i <= (cellNum - 1); i += MAXTHREAD)
	{
		int find = binary_search_intPair(FVTable, taskInfo.int_2, (nextCnt - 1), (i + threadID));
		//int find = 1;
		//int k = cellNum*blockID + (i + threadID);
		//queryFVGPU[cellNum*blockID + (i + threadID)] = 2;
		queryFVGPU[cellNum*blockID + (i + threadID)] = queryFVGPU[cellNum*blockID + (i + threadID)] - find;
	}
	//�ڶ��׶Σ��������ڣ�������
	//����׶θ�Ϊÿ��thread����һ��FD
	int neighborsID[8];
	for (int cell = 0; cell <= cellNum - 1; cell++)
	{
		//ֻ��Ҫһ�����߳̾�����
		if (threadIDGlobal >= checkNum)
			break;
		if (queryFVGPU[cellNum*threadIDGlobal + cell] != 0)
		{
			findNeighborGPU(cellNum, cell, neighborsID);
			//for (int i = 0; i <= 7; i++)
			//	neighborsID[i] = 11;
			for (int i = 0; i <= 7; i++)
			{
				if (isPositive(queryFVGPU[cellNum*threadIDGlobal + cell]) != isPositive(queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]])) {
					if (fabsf(queryFVGPU[cellNum*threadIDGlobal + cell]) > fabsf(queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]]))
					{
						queryFVGPU[cellNum*threadIDGlobal + cell] = queryFVGPU[cellNum*threadIDGlobal + cell] + queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]];
						queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]] = 0;
					}
					else
					{
						queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]] = queryFVGPU[cellNum*threadIDGlobal + neighborsID[i]] + queryFVGPU[cellNum*threadIDGlobal + cell];
						queryFVGPU[cellNum*threadIDGlobal + cell] = 0;
						break;
					}
				}
			}
		}
	}
	__syncthreads();
	//�����׶Σ�ͳ����������
	//��Ȼ��ÿ��block����һ��FD�ļ���
	__shared__ int tempsumPosi[MAXTHREAD], tempsumNega[MAXTHREAD];
	tempsumPosi[threadID] = 0;
	tempsumNega[threadID] = 0;
	for (int i = 0; i <= cellNum - 1; i += MAXTHREAD)
	{
		tempsumPosi[threadID] += (isPositive(queryFVGPU[blockID*cellNum + (i + threadID)])*queryFVGPU[blockID*cellNum + (i + threadID)]);
		tempsumNega[threadID] += (-(!isPositive(queryFVGPU[blockID*cellNum + (i + threadID)]))*queryFVGPU[blockID*cellNum + (i + threadID)]);
	}
	__shared__ int sizeOfTempSum;
	if (threadID == 0)
		sizeOfTempSum = MAXTHREAD;
	__syncthreads();
	while ((sizeOfTempSum>1))
	{
		if (threadID <= (sizeOfTempSum >> 1) - 1)
		{
			tempsumPosi[threadID] = tempsumPosi[threadID] + tempsumPosi[threadID + (sizeOfTempSum >> 1)];
			tempsumNega[threadID] = tempsumNega[threadID] + tempsumNega[threadID + (sizeOfTempSum >> 1)];
		}
		__syncthreads();
		if (threadID == 0)
			sizeOfTempSum = (sizeOfTempSum >> 1);
		__syncthreads();
	}
	if (threadID == 0)
		FDistance[blockID] = (tempsumPosi[0] > tempsumNega[0]) ? tempsumPosi[0] : tempsumNega[0];

}

//SubbedArrayJump��SubbedArray��ÿһ���ж��ٸ�Ԫ�أ�������idx��
int Similarity_Pruning_Handler(intPair* queryFVGPU, intPair* FVinfo, intPair* FVTable, intPair* SubbedArray, intPair* SubbedArrayOffset, int SubbedArrayJump, int queryCellLength, int startTrajIdx, int checkNum, int cellNum, int trajNumInDB, int nonZeroFVNumInDB, short* FDistance, cudaStream_t stream)
{
#ifdef NOT_COLUMN_ORIENTED
	Calculate_FD_NonColumn << <checkNum, MAXTHREAD, 0, stream >> >(queryFVGPU, FVinfo, FVTable, startTrajIdx, checkNum, cellNum, trajNumInDB, nonZeroFVNumInDB, FDistance);
#else
	Calculate_FD_Sparse << <checkNum, MAXTHREAD, 0, stream >> >(queryFVGPU, FVinfo, FVTable, SubbedArray, SubbedArrayOffset, SubbedArrayJump, queryCellLength, startTrajIdx, checkNum, cellNum, trajNumInDB, nonZeroFVNumInDB, FDistance);
#endif
	return 0;
}


/*
//�Ȱ����ܷ���һ��SMִ��һ��DP�����������ٷֱ��������kernel
//constructing...
���Ż���
1��queryTra��queryLength����candidateLength����ͨ����ֵ�ķ�ʽֱ�Ӵ��ݵ�SM�ļĴ���������ȫ���ڴ��ʹ��

*/
int handleEDRdistance(SPoint *queryTra, SPoint **candidateTra, int candidateNum, int queryLength, int *candidateLength, int *result) {
	MyTimer time1;
	time1.start();

	int** stateTableGPU = NULL;
	//��GPU��Ϊ״̬������ڴ�
	int** temp = NULL;
	temp = (int**)malloc(sizeof(int*)*candidateNum);
	for (int i = 0; i <= candidateNum - 1; i++) {
		CUDA_CALL(cudaMalloc((void**)&temp[i], sizeof(int)*(candidateLength[i] + 1)*(queryLength + 1)));
	}
	CUDA_CALL(cudaMalloc((void***)&stateTableGPU, sizeof(int*)*candidateNum));
	CUDA_CALL(cudaMemcpy(stateTableGPU, temp, candidateNum*sizeof(int*), cudaMemcpyHostToDevice));

	//Ϊ�洢�Ĺ켣��Ϣ�����ڴ�
	SPoint *queryTraGPU = NULL, **candidateTraGPU = NULL;
	int *candidateLengthGPU = NULL, *resultGPU = NULL;
	CUDA_CALL(cudaMalloc((void**)&queryTraGPU, sizeof(SPoint)*queryLength));
	CUDA_CALL(cudaMalloc((void**)&candidateLengthGPU, sizeof(int)*candidateNum));
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, sizeof(int)*candidateNum));

	SPoint **tempS = (SPoint**)malloc(sizeof(SPoint*)*candidateNum);
	for (int i = 0; i <= candidateNum - 1; i++) {
		CUDA_CALL(cudaMalloc((void**)&tempS[i], sizeof(SPoint)*candidateLength[i]));

	}
	CUDA_CALL(cudaMalloc((void***)&candidateTraGPU, sizeof(SPoint*)*candidateNum));
	CUDA_CALL(cudaMemcpy(candidateTraGPU, tempS, candidateNum*sizeof(SPoint*), cudaMemcpyHostToDevice));
	//
	time1.stop();
	std::cout << time1.elapse() << std::endl;
	time1.start();
	//
	//���ͨ���������ķ������ݹ켣�����Ҫ��켣�����洢
	//��GPU���ݹ켣��Ϣ
	CUDA_CALL(cudaMemcpy(queryTraGPU, queryTra, queryLength*sizeof(SPoint), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(candidateLengthGPU, candidateLength, candidateNum*sizeof(int), cudaMemcpyHostToDevice));

	for (int i = 0; i <= candidateNum - 1; i++) {
		CUDA_CALL(cudaMemcpy(tempS[i], candidateTra[i], candidateLength[i] * sizeof(SPoint), cudaMemcpyHostToDevice));
	}
	//for (int i = 0; i <= candidateNum - 1;i++)
	//	CUDA_CALL(cudaMemcpy(candidateTraGPU[i], candidateTra[i], candidateLength[i]*sizeof(SPoint), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaHostAlloc((void**)&result, candidateNum*sizeof(int), cudaHostAllocWriteCombined | cudaHostAllocMapped));
	CUDA_CALL(cudaHostGetDevicePointer(&resultGPU, result, 0));
	time1.stop();
	std::cout << time1.elapse() << std::endl;
	time1.start();
	//ִ��kernel
	EDRDistance_1 << <candidateNum, MAXTHREAD >> >(queryTraGPU, candidateTraGPU, candidateNum, queryLength, candidateLengthGPU, stateTableGPU, resultGPU);

	//ȡ���
	//result = (int*)malloc(candidateNum*sizeof(int));
	//CUDA_CALL(cudaMemcpy(result, resultGPU, candidateNum*sizeof(int), cudaMemcpyDeviceToHost));
	cudaDeviceSynchronize();
	//	for (int j = 0; j <= candidateNum - 1;j++)
	//		std::cout << result[j] << std::endl;

	//free GPU!!!!!
	time1.stop();
	std::cout << time1.elapse() << std::endl;
	return 0;

}


inline void __getLastCudaError(const char *errorMessage, const char *file, const int line)
{
	cudaError_t err = cudaGetLastError();

	if (cudaSuccess != err)
	{
		fprintf(stderr, "%s(%i) : getLastCudaError() CUDA error : %s : (%d) %s.\n",
			file, line, errorMessage, (int)err, cudaGetErrorString(err));
			exit(EXIT_FAILURE);
	}
}


//using namespace thrust;
//static const int MAXTHREAD = 512; //ÿ��block�߳���

cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);
void CUDAwarmUp() {
	CUDA_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));
	CUDA_CALL(cudaSetDevice(0));
	
}

#ifdef _CELL_BASED_STORAGE
int putCellDataSetIntoGPU(Point* pointsPtr, Point*& pointsPtrGPU, int pointNum) {
	
	CUDA_CALL(cudaMalloc((void**)&pointsPtrGPU, pointNum * sizeof(Point))); //�������ݵ��ڴ�
	//debug
	//std::cout << pointNum << std::endl;
	//debug
	CUDA_CALL(cudaMemcpy(pointsPtrGPU, pointsPtr, pointNum * sizeof(Point), cudaMemcpyHostToDevice));//���ݿ�����gpu��
	return 0;
}
__global__ void cudaRangeQuery(int* rangeStarts, int* rangeEnds, int candidateCellNum, const Point* pointsPtr, const float xmin, const float ymin, const float xmax, const float ymax, const int *resultOffset, Point* resultPtrCuda) {
	int cellNo = blockIdx.x; //candidate����ڼ���cell 0,1,2,....
	if (cellNo >= candidateCellNum) return;
	int tid = threadIdx.x;
	if (tid >= 256) return;
	int pointNum = rangeEnds[cellNo] - rangeStarts[cellNo] + 1;//blockҪ��������cell����ô�����
	const int offset = rangeStarts[cellNo];
	for (int i = tid; i <= pointNum - 1; i += MAXTHREAD) {
		float x = pointsPtr[offset + i].x;
		float y = pointsPtr[offset + i].y;
		uint32_t tid = pointsPtr[offset + i].tID;
		uint32_t time = pointsPtr[offset + i].time;
		if (x <= xmax &&x >= xmin&&y <= ymax&&y >= ymin) {
			resultPtrCuda[resultOffset[cellNo] + i].x = x;
			resultPtrCuda[resultOffset[cellNo] + i].y = y;
			resultPtrCuda[resultOffset[cellNo] + i].tID = tid;
			resultPtrCuda[resultOffset[cellNo] + i].time = time;
		}
		else
			resultPtrCuda[resultOffset[cellNo] + i].tID = -1;
	}
}

__global__ void cudaRangeQueryTest(RangeQueryStateTable* stateTable, int stateTableLength, uint8_t* result, 
	const int maxTrajNum) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;

	__shared__ RangeQueryStateTable sharedStateTable;
	// __shared__ uint8_t resultTemp[10000]; //10K
	
	if(tID == 0)
		sharedStateTable = (stateTable[bID]); // ���岻�� �����������Ƴ����ڴ���

	__syncthreads();//block��threadͬ�� 

	/*
// 4+4*7=32byte
typedef struct RangeQueryStateTable {//  leafnode

	// GPU���
	void* ptr;	// ָ��GPU�ڴ�node��ָ�� �����洢
	int candidatePointNum;	// �����leafnode�еĽڵ���

	float xmin;
	float ymin;
	float xmax;
	float ymax;

	// cpu���
	int queryID;
	int startIdxInAllPoints; //startId in AllPoints����

}RangeQueryStateTable;*/

	int jobID = sharedStateTable.queryID;

	SPoint *baseAddr = (SPoint*)(sharedStateTable.ptr);
	int candidateNum = sharedStateTable.candidatePointNum;

	//int resultOffset = bID*maxPointNumInStateTable; 
	SPoint p;
	/*
	for (int i = 0; i <= candidateNum / MAXTHREAD-1; i++) {
		//p = *(baseAddr + (i*MAXTHREAD + tID));
		p = baseAddr[i*MAXTHREAD + tID];
		//result[i*MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
			//(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(i*MAXTHREAD + tID);//�����֤ͨ�������ֵΪ�����ţ�����Ϊ0
		//result[i*MAXTHREAD + tID + resultOffset].jobID = bID;
		//result[resultOffset + (i*MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && 
		//		(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
		if((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
			result[jobID*maxTrajNum + p.tID] = 1;
		//�����֤ͨ��������Ӧλ����Ϊ1

		//__syncthreads();
	}
	if (tID < candidateNum - candidateNum / MAXTHREAD * MAXTHREAD) {
		p = *(baseAddr + (candidateNum / MAXTHREAD * MAXTHREAD + tID));
		//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
		//	(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(candidateNum / MAXTHREAD * MAXTHREAD + tID);//�����֤ͨ�������ֵΪ�����ţ�����Ϊ0
		//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].jobID = bID;
		//result[resultOffset + (candidateNum / MAXTHREAD * MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) &&
		//	(p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
		if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
			result[jobID*maxTrajNum + p.tID] = 1;
	}
	*/

	/*
	example1:
	int blockSize = 256;
	int numBlocks = (N + blockSize - 1) / blockSize;
	add<<<numBlocks, blockSize>>>(N, x, y);
	__global__
	void add(int n, float *x, float *y)
	{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int stride = blockDim.x * gridDim.x; // ��������stride��ͬ��
	for (int i = index; i < n; i += stride)
	y[i] = x[i] + y[i];
	}

	example2:
	add<<<1, 256>>>(N, x, y);
	__global__
	void add(int n, float *x, float *y)
	{
	  int index = threadIdx.x;
	  int stride = blockDim.x;
	  for (int i = index; i < n; i += stride)
		  y[i] = x[i] + y[i];
	}
	*/

	//  a new version
	// grid-stride loop��ʽ
	//int index = blockIdx.x * blockDim.x + threadIdx.x;
	//int stride = blockDim.x * gridDim.x;
	//for (int i = index; i < n; i += stride) {
	//

	//}

	// stride = MAXTHREAD(blockDim.x)

	// kernel ��ÿ��thread �����첽ִ��
	//new version

	//for (int i = 0; i < candidateNum; i+=MAXTHREAD) {
	//	//p = *(baseAddr + (i*MAXTHREAD + tID));
	//	if (i + tID < candidateNum) {
	//		p = baseAddr[i + tID]; // �����ȫ���ڴ�ȡ��
	//					// p�Ǵ�ȫ���ڴ�ȡ��
	//		//result[i*MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
	//		//(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(i*MAXTHREAD + tID);//�����֤ͨ�������ֵΪ�����ţ�����Ϊ0
	//		//result[i*MAXTHREAD + tID + resultOffset].jobID = bID;
	//		//result[resultOffset + (i*MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && 
	//		//		(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
	//		if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
	//			result[jobID*maxTrajNum + p.tID] = 1;
	//		//�����֤ͨ��������Ӧλ����Ϊ1
	//	}
	//	//__syncthreads();
	//}




	for (int i = 0; i+ tID < candidateNum; i += MAXTHREAD) {
			//p = *(baseAddr + (i*MAXTHREAD + tID));

			p = baseAddr[i + tID];	// fetch from global memory

			if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
				result[jobID*maxTrajNum + p.tID] = 1; // write into global memory
			//jobID*maxTrajNum + p.tID jobID=queryID 1-40 maxTrajNum = this->trajNum + 1  p.tID is not thread, is trajID
			// basic idea
	}


	return;
	


	//else {
	//	//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].idx = 0; //������Ĳ��֣�ֱ����Ϊ��Ч����
	//}
	//__syncthreads();
	//__syncthreads();
	//int globalTID = blockDim.x * blockIdx.x + threadIdx.x;
	//if (globalTID < stateTableLength) {

	//}
}

// offsetNum: how many uncontinuous part in each block
// offsetLen: how many elements in this uncontinuous part
// offset: the offset of all uncontinuous parts
// offsetInoffset: the offset in offset array for each block
__global__ void cudaRangeQueryTestWithoutMorton(RangeQueryStateTable* stateTable, int stateTableLength, uint8_t* result,
	const int maxTrajNum, int* offset, int* offsetLen, int* offsetNum, int* offsetInOffset) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	__shared__ RangeQueryStateTable sharedStateTable;
	// __shared__ uint8_t resultTemp[10000]; //10K
	if (tID == 0)
		sharedStateTable = (stateTable[bID]);
	__syncthreads();
	int jobID = sharedStateTable.queryID;
	SPoint *baseAddr = (SPoint*)(sharedStateTable.ptr);
	int candidateNum = sharedStateTable.candidatePointNum;//��block����Ҫ��ѯ�ĵ�ĸ���
														  //int resultOffset = bID*maxPointNumInStateTable; //��block�Ľ������ʼ��ַ
	SPoint p;
	// all offset of start of array in this block
	__shared__ int offsetLocal[1000];
	__shared__ int offsetLenLocal[1000];
	int continuousNum = offsetNum[bID];
	for (int i = 0; i < continuousNum; i += MAXTHREAD) {
		if (i + tID < continuousNum) {
			offsetLocal[i + tID] = offset[offsetInOffset[bID] + i + tID];
			offsetLenLocal[i + tID] = offsetLen[offsetInOffset[bID] + i + tID];
		}
	}
	__syncthreads();
	/*
	for (int i = 0; i <= candidateNum / MAXTHREAD-1; i++) {
	//p = *(baseAddr + (i*MAXTHREAD + tID));
	p = baseAddr[i*MAXTHREAD + tID];
	//result[i*MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
	//(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(i*MAXTHREAD + tID);//�����֤ͨ�������ֵΪ�����ţ�����Ϊ0
	//result[i*MAXTHREAD + tID + resultOffset].jobID = bID;
	//result[resultOffset + (i*MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
	//		(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
	if((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
	result[jobID*maxTrajNum + p.tID] = 1;
	//�����֤ͨ��������Ӧλ����Ϊ1

	//__syncthreads();
	}
	if (tID < candidateNum - candidateNum / MAXTHREAD * MAXTHREAD) {
	p = *(baseAddr + (candidateNum / MAXTHREAD * MAXTHREAD + tID));
	//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
	//	(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(candidateNum / MAXTHREAD * MAXTHREAD + tID);//�����֤ͨ�������ֵΪ�����ţ�����Ϊ0
	//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].jobID = bID;
	//result[resultOffset + (candidateNum / MAXTHREAD * MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) &&
	//	(p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
	if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
	result[jobID*maxTrajNum + p.tID] = 1;
	}
	*/

	//new version
	for (int i = 0; i < continuousNum; i++) {
		int offsetLength = offsetLenLocal[i];
		int offsetAddr = offsetLocal[i];
		for (int j = 0; j < offsetLength; j += MAXTHREAD) {
			if (j + tID < offsetLength) {
				p = baseAddr[offsetAddr + j + tID];
				//result[i*MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
				//(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(i*MAXTHREAD + tID);//�����֤ͨ�������ֵΪ�����ţ�����Ϊ0
				//result[i*MAXTHREAD + tID + resultOffset].jobID = bID;
				//result[resultOffset + (i*MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && 
				//		(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
				if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
					result[jobID*maxTrajNum + p.tID] = 1;
				//�����֤ͨ��������Ӧλ����Ϊ1
			}
		}
	}
	return;



	//else {
	//	//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].idx = 0; //������Ĳ��֣�ֱ����Ϊ��Ч����
	//}
	//__syncthreads();
	//__syncthreads();
	//int globalTID = blockDim.x * blockIdx.x + threadIdx.x;
	//if (globalTID < stateTableLength) {

	//}
}

__global__ void cudaRangeQuerySTIG(RangeQueryStateTable* stateTable, int stateTableLength, uint8_t* result,
	const int maxTrajNum) {
	int bID = blockIdx.x;
	int tID = threadIdx.x;
	//if (bID > stateTableLength)
	//	return;
	__shared__ RangeQueryStateTable sharedStateTable;
	// __shared__ uint8_t resultTemp[10000]; //10K
	if (tID == 0)
		sharedStateTable = (stateTable[bID]);
	__syncthreads();
	int jobID = sharedStateTable.queryID;
	SPoint *baseAddr = (SPoint*)(sharedStateTable.ptr);
	int candidateNum = sharedStateTable.candidatePointNum;//��block����Ҫ��ѯ�ĵ�ĸ���
														  //int resultOffset = bID*maxPointNumInStateTable; //��block�Ľ������ʼ��ַ
	SPoint p;
	//new version
	for (int i = 0; i < candidateNum; i += MAXTHREAD) {
		//p = *(baseAddr + (i*MAXTHREAD + tID));
		if (i + tID < candidateNum) {
			p = baseAddr[i + tID];
			//result[i*MAXTHREAD + tID + resultOffset].idx = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) &&
			//(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))*(i*MAXTHREAD + tID);//�����֤ͨ�������ֵΪ�����ţ�����Ϊ0
			//result[i*MAXTHREAD + tID + resultOffset].jobID = bID;
			//result[resultOffset + (i*MAXTHREAD + tID)] = ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && 
			//		(p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin));
			if ((p.x<sharedStateTable.xmax) && (p.x>sharedStateTable.xmin) && (p.y<sharedStateTable.ymax) && (p.y>sharedStateTable.ymin))
				result[jobID*maxTrajNum + p.tID] = 1;
			//�����֤ͨ��������Ӧλ����Ϊ1
		}
		//__syncthreads();
	}



	//else {
	//	//result[candidateNum / MAXTHREAD * MAXTHREAD + tID + resultOffset].idx = 0; //������Ĳ��֣�ֱ����Ϊ��Ч����
	//}
	//__syncthreads();
	//__syncthreads();
	//int globalTID = blockDim.x * blockIdx.x + threadIdx.x;
	//if (globalTID < stateTableLength) {

	//}
}

int cudaRangeQuerySTIGHandler(RangeQueryStateTable* stateTableGPU, int stateTableLength, uint8_t *result, int maxTrajNum
	, int maxQueryNum, cudaStream_t stream)
{
	//RangeQueryResultGPU* resultGPU;
	//MyTimer timer;
	uint8_t* resultGPU;
	//int resultByteNum = (maxPointNum)/8+1;//ÿ�������Ҫ�ü���byte���棬���ܰ����أ�ֻ�ܰ��ֽ�
	CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxTrajNum)*maxQueryNum));//selective��һ��
	CUDA_CALL(cudaMemset(resultGPU, 0, (maxTrajNum)*maxQueryNum));
	//timer.start();
	//�����һ���ڴ棬ÿ��stateTable��ռ�ݵ��ڴ������
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxPointNum)*stateTableLength));

	//CUDA_CALL(cudaMalloc((void**)&resultGPU, maxPointNum*stateTableLength*sizeof(RangeQueryResultGPU)));
	//timer.stop();
	//std::cout << "Time 1:" << timer.elapse() << "ms" << std::endl;

	//timer.start();	
	cudaRangeQuerySTIG << <stateTableLength, MAXTHREAD, 0, stream >> >(stateTableGPU, stateTableLength, resultGPU, maxTrajNum);
	CUDA_CALL(cudaDeviceSynchronize());
	//timer.stop();
	//std::cout << "Time 2:" << timer.elapse() << "ms" << std::endl;

	//timer.start();

	CUDA_CALL(cudaMemcpy(result, resultGPU, (maxTrajNum)*maxQueryNum, cudaMemcpyDeviceToHost));

	//timer.stop();
	//std::cout << "Time 3:" << timer.elapse() << "ms" << std::endl;
	CUDA_CALL(cudaFree(resultGPU));
	return 0;
}

int cudaRangeQueryTestHandlerNonMorton(RangeQueryStateTable* stateTableGPU, int stateTableLength, uint8_t *result, int maxTrajNum
	, int maxJobNum, cudaStream_t stream, int* offset, int* offsetLen, int* offsetNum, int* offsetInOffset) {
	/*
	without Morton encoding, each line should be processed seperately.
	divide loops to process each line (continuous points).

	*/
	//RangeQueryResultGPU* resultGPU;
	//MyTimer timer;
	uint8_t* resultGPU;
	//int resultByteNum = (maxPointNum)/8+1;//ÿ�������Ҫ�ü���byte���棬���ܰ����أ�ֻ�ܰ��ֽ�

	CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxTrajNum)*maxJobNum));//selective��һ��
	CUDA_CALL(cudaMemset(resultGPU, 0, (maxTrajNum)*maxJobNum));
	
	//timer.start();
	//�����һ���ڴ棬ÿ��stateTable��ռ�ݵ��ڴ������
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxPointNum)*stateTableLength));

	//CUDA_CALL(cudaMalloc((void**)&resultGPU, maxPointNum*stateTableLength*sizeof(RangeQueryResultGPU)));
	//timer.stop();
	//std::cout << "Time 1:" << timer.elapse() << "ms" << std::endl;

	//timer.start();	
	cudaRangeQueryTestWithoutMorton << <stateTableLength, MAXTHREAD, 0, stream >> >(stateTableGPU, stateTableLength, 
		resultGPU, maxTrajNum, offset, offsetLen, offsetNum, offsetInOffset);
	CUDA_CALL(cudaDeviceSynchronize());
	//timer.stop();
	//std::cout << "Time 2:" << timer.elapse() << "ms" << std::endl;
	//timer.start();

	CUDA_CALL(cudaMemcpy(result, resultGPU, (maxTrajNum)*maxJobNum, cudaMemcpyDeviceToHost));
	CUDA_CALL(cudaFree(resultGPU));

	//timer.stop();
	//std::cout << "Time 3:" << timer.elapse() << "ms" << std::endl;
	return 0;
}

// cudaRangeQueryTestHandler((RangeQueryStateTable*)this->stateTableGPU[device_idx], this->stateTableLength[device_idx], resultsReturned, this->trajNum + 1, rangeNum, stream);

int cudaRangeQueryTestHandler(RangeQueryStateTable* stateTableGPU, int stateTableLength, uint8_t *result, int maxTrajNum
	, int maxJobNum, cudaStream_t stream) {
	//RangeQueryResultGPU* resultGPU;
	//MyTimer timer;

	uint8_t* resultGPU; // ָ��GPU
	//int resultByteNum = (maxPointNum)/8+1;//ÿ�������Ҫ�ü���byte���棬���ܰ����أ�ֻ�ܰ��ֽ�

	CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxTrajNum)*maxJobNum));//selective��һ��
	CUDA_CALL(cudaMemset(resultGPU, 0, (maxTrajNum)*maxJobNum));

	//timer.start();
	//�����һ���ڴ棬ÿ��stateTable��ռ�ݵ��ڴ������
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, (maxPointNum)*stateTableLength));
	
	//CUDA_CALL(cudaMalloc((void**)&resultGPU, maxPointNum*stateTableLength*sizeof(RangeQueryResultGPU)));
	//timer.stop();
	//std::cout << "Time 1:" << timer.elapse() << "ms" << std::endl;
	
	//timer.start();

	//����stateTableLength��filteringʱ���壬1 block ��Ӧһ�� quad-fake-leaf-node 
	cudaRangeQueryTest <<< stateTableLength, MAXTHREAD,0, stream >>>(stateTableGPU, stateTableLength, resultGPU, maxTrajNum);
	
	// // Wait for GPU to finish before accessing on host ��ʽͬ���������õ��� resultGPU

	CUDA_CALL(cudaDeviceSynchronize());//2GPUҲ����ͨ�ţ��� ������


	// CUDA_CALL(cudaDeviceSynchronize(stream)); // Ϊʲôû��stream

	//timer.stop();
	//std::cout << "Time 2:" << timer.elapse() << "ms" << std::endl;
	//timer.start();
	
	// GPU->CPU ���ؽ��
	CUDA_CALL(cudaMemcpyAsync(result, resultGPU, (maxTrajNum)*maxJobNum, cudaMemcpyDeviceToHost));
	//CUDA_CALL(cudaDeviceSynchronize());
	CUDA_CALL(cudaFree(resultGPU)); // not recommended

	//timer.stop();
	//std::cout << "Time 3:" << timer.elapse() << "ms" << std::endl;
	return 0;
}

/*
int cudaRangeQueryHandler(int* candidateCells, int* rangeStarts, int* rangeEnds, int candidateCellNum,float xmin, float ymin, float xmax, float ymax, Point*& resultsGPU, int& resultNum,Point *pointsPtrGPU,Point *&result) {
	//��һ��������ʱû���ã�ע������candidatecells[i]�Ѿ�������cell��id������ֻ��ǿ�
	//���ĸ�������ʾ�ǿյ�cell����
	//����candidate���еĵ��������gpu�ڿ�����ͬ��С�Ŀռ���flag��rangestart��rangeend����Ӧcandidatecell�ڵĲ�������AllPoints����ʼ�±����ֹ�±�
	//�����������͵����ڶ��������������һ���Ǳ�������GPU��ַ���ڶ����ǽ���ĸ���
	//PointsPtrGPU�����ݼ���gpu�ĵ�ַ
	MyTimer timer1;
	timer1.start();
	int counter = 0;
	int *resultOffset = (int*)malloc(sizeof(int)*candidateCellNum);
	//std::cout << candidateCellNum << ":"<<std::endl;
	for (int i = 0; i <= candidateCellNum - 1; i++) {
		resultOffset[i] = counter;
		////debug
		//std::cout << "(" << rangeStarts[i] << "," << rangeEnds[i] << ");"<<"["<<resultOffset[i]<<"]";
		////debug
		counter += rangeEnds[i] - rangeStarts[i] + 1;
	}
	int totalPointNumInCandidate = counter;


	int *rangeStartsCuda = NULL, *rangeEndsCuda = NULL, *resultOffsetCuda = NULL;

	CUDA_CALL(cudaMalloc((void**)&resultsGPU, sizeof(Point)*totalPointNumInCandidate));
	//��range��cell��Ϣд��gpu
	//CUDA_CALL(cudaMalloc((void**)&candidateCellsCuda, sizeof(int)*candidateCellNum));
	CUDA_CALL(cudaMalloc((void**)&rangeStartsCuda, candidateCellNum*sizeof(int)));
	//std::cout << "\n" << candidateCellNum*sizeof(int) << "\n";
	CUDA_CALL(cudaMalloc((void**)&rangeEndsCuda, candidateCellNum*sizeof(int)));
	CUDA_CALL(cudaMalloc((void**)&resultOffsetCuda, candidateCellNum*sizeof(int)));
	//CUDA_CALL(cudaMemcpy(candidateCellsCuda, candidateCells, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(rangeStartsCuda, rangeStarts, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));

	CUDA_CALL(cudaMemcpy(rangeEndsCuda, rangeEnds, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_CALL(cudaMemcpy(resultOffsetCuda, resultOffset, candidateCellNum*sizeof(int), cudaMemcpyHostToDevice));
	////debug
	//CUDA_CALL(cudaMemcpy( rangeStarts, rangeStartsCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//CUDA_CALL(cudaMemcpy( rangeEnds, rangeEndsCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//CUDA_CALL(cudaMemcpy( resultOffset, resultOffsetCuda, candidateCellNum*sizeof(int), cudaMemcpyDeviceToHost));
	//for (int i = 0; i <= candidateCellNum - 1; i++) {
	//	//debug
	//	std::cout << "(" << rangeStarts[i] << "," << rangeEnds[i] << ");" << "[" << resultOffset[i] << "]";
	//	//debug
	//}
	//debug
	timer1.stop();
	std::cout << timer1.ticks() << std::endl;
	timer1.start();
	//����kernel�����ĳ����������������Ӧλ��д����AllPoints�е��±꣬����д��-1
	//ÿ��cell�����һ��block
	cudaRangeQuery <<<candidateCellNum, MAXTHREAD >>>(rangeStartsCuda, rangeEndsCuda, candidateCellNum, pointsPtrGPU, xmin, ymin, xmax, ymax, resultOffsetCuda, resultsGPU);
	//kernel���ý��������������idxsGPU�У����������������ӦԪ����������AllPoints���±꣬��������ϣ�����Ϊ-1
	//CUDA_CALL(cudaFree(candidateCellsCuda));
	CUDA_CALL(cudaFree(rangeStartsCuda));
	CUDA_CALL(cudaFree(rangeEndsCuda));
	CUDA_CALL(cudaFree(resultOffsetCuda));
	//getLastCudaError("Error in Calling 'kernel'");
	//ʹ��Thrustɾ������-1���õ����ս��
	timer1.stop();
	std::cout << timer1.ticks() << std::endl;

	//���н���ϲ�
	//test

	timer1.start();
	Point *resultset = NULL;
	resultset = (Point*)malloc(totalPointNumInCandidate*sizeof(Point));
	CUDA_CALL(cudaMemcpy(resultset, resultsGPU, sizeof(Point)*totalPointNumInCandidate, cudaMemcpyDeviceToHost));
	std::vector<Point> *resultPoint = new std::vector<Point>;
	for (int i = 0; i <= totalPointNumInCandidate - 1; i++) {
		if (resultset[i].tID != -1)
		{
			resultPoint->push_back(resultset[i]);
		}
	}
	result = &resultPoint->at(0);
	free(resultset);
	//test
	timer1.stop();
	std::cout << timer1.ticks() << std::endl;
	
	//���н���ϲ�
	
	//thrust::device_ptr<int> idxsPtr = thrust::device_pointer_cast(idxsGPU);
	//int a;
	//cudaMemcpy(&a, idxsGPU, 1, cudaMemcpyDeviceToHost);
	//size_t num = thrust::remove(idxsPtr, idxsPtr + totalPointNumInCandidate-1, -1) - idxsPtr;
	//int *result = (int*)malloc(sizeof(int)*num);
	//thrust::copy(idxsPtr, idxsPtr + num, result);
	//resultNum = num;
	//resultIdx = result;

	//CUDA_CALL(cudaFree(idxsGPU));


	return 0;
}
*/


#else
int cudaRangeQueryHandler(Point* pointsPtr, int pointNum, float xmin, float ymin, float xmax, float ymax, Point*& resultsPtr, int& resultNum) {
	Point* pointsPtrCuda = NULL;
	Point* resultPtrCuda = NULL;
	CUDA_CALL(cudaMalloc((void**)&pointsPtrCuda, pointNum * sizeof(Point))); //�������ݵ��ڴ�
	CUDA_CALL(cudaMalloc((void**)&resultPtrCuda, pointNum * sizeof(Point))); //gpu�ڴ洢����ĵط�
	CUDA_CALL(cudaMemcpy(pointsPtrCuda, pointsPtr, pointNum * sizeof(Point), cudaMemcpyHostToDevice));//���ݿ�����gpu��

																									  //���ú˺����������ݣ��������gpu��

																									  //ȡ�����ݣ�����
	return 0;
}
#endif




//__global__ void addKernel(int *c, const int *a, const int *b)
//{
//    int i = threadIdx.x;
//    c[i] = a[i] + b[i];
//}
//int main()
//{
//    const int arraySize = 5;
//    const int a[arraySize] = { 1, 2, 3, 4, 5 };
//    const int b[arraySize] = { 10, 20, 30, 40, 50 };
//    int c[arraySize] = { 0 };
//
//    // Add vectors in parallel.
//    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "addWithCuda failed!");
//        return 1;
//    }
//
//    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
//        c[0], c[1], c[2], c[3], c[4]);
//
//    // cudaDeviceReset must be called before exiting in order for profiling and
//    // tracing tools such as Nsight and Visual Profiler to show complete traces.
//    cudaStatus = cudaDeviceReset();
//    if (cudaStatus != cudaSuccess) {
//        fprintf(stderr, "cudaDeviceReset failed!");
//        return 1;
//    }
//
//    return 0;
//}
// Helper function for using CUDA to add vectors in parallel.
