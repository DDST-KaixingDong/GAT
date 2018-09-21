// Data_Structure_Test.cpp
// 20180521 21:01
// test2ST
#ifndef WIN32
#include <unistd.h>
#else
#include <Windows.h>
#define sleep(x) Sleep(x*1000)
#endif

#include <iostream>
#include <fstream>
#include <map>
#include "CalcInAxis.h"
#include "PreProcess.h"
#include "Trajectory.h"
#include <vector>
#include "Grid.h"
#include "cudaKernel.h"
#include "SystemTest.h"
#include "STIG.h"
#include "MortonGrid.h"
#include <stdlib.h>


using namespace std;

//global
map<string, tidLinkTable*> vidTotid;
map<string, tidLinkTable*>::iterator iter;

//global
Trajectory* tradb;
string baseDate = "2014-07-01";
void* baseAddrGPU = NULL; // GPU内存基地址

bool nonmem_cmp(const Trajectory &t1, const Trajectory &t2)
{
	if (t1.length < t2.length)
		return true;
	return false;
}

Trajectory splitTrajectory(const Trajectory* cand, int parts)
{
	Trajectory t;
	int candLen = cand->length;
	t.tid = 888888;
	t.length = 0;
	for (int i = 0; i <= candLen - 2;i++)
	{
		SamplePoint p1 = cand->points[i];
		SamplePoint p2 = cand->points[i + 1];
		switch (parts)
		{
		case 1:
			t.points.push_back(p1);
			t.length += 1;
			break;
		case 2:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.5, p1.lat + (p2.lat - p1.lat)*0.5, 1000, 888888));
			t.length += 2;
			break;
		case 3:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.33, p1.lat + (p2.lat - p1.lat)*0.33, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.66, p1.lat + (p2.lat - p1.lat)*0.66, 1000, 888888));
			t.length += 3;
			break;
		case 4:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.25, p1.lat + (p2.lat - p1.lat)*0.25, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.5, p1.lat + (p2.lat - p1.lat)*0.5, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.75, p1.lat + (p2.lat - p1.lat)*0.75, 1000, 888888));
			t.length += 4;
			break;
		case 5:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.2, p1.lat + (p2.lat - p1.lat)*0.2, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.4, p1.lat + (p2.lat - p1.lat)*0.4, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.6, p1.lat + (p2.lat - p1.lat)*0.6, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.8, p1.lat + (p2.lat - p1.lat)*0.8, 1000, 888888));
			t.length += 5;
			break;
		case 6:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.17, p1.lat + (p2.lat - p1.lat)*0.17, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.34, p1.lat + (p2.lat - p1.lat)*0.34, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.51, p1.lat + (p2.lat - p1.lat)*0.51, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.68, p1.lat + (p2.lat - p1.lat)*0.68, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*0.85, p1.lat + (p2.lat - p1.lat)*0.85, 1000, 888888));
			t.length += 6;
			break;
		case 7:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*1.0 / 7.0, p1.lat + (p2.lat - p1.lat)*1.0 / 7.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*2.0 / 7.0, p1.lat + (p2.lat - p1.lat)*2.0 / 7.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*3.0 / 7.0, p1.lat + (p2.lat - p1.lat)*3.0 / 7.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*4.0 / 7.0, p1.lat + (p2.lat - p1.lat)*4.0 / 7.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*5.0 / 7.0, p1.lat + (p2.lat - p1.lat)*5.0 / 7.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*6.0 / 7.0, p1.lat + (p2.lat - p1.lat)*6.0 / 7.0, 1000, 888888));
			t.length += 7;
			break;
		case 8:
			t.points.push_back(p1);
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*1.0 / 8.0, p1.lat + (p2.lat - p1.lat)*1.0 / 8.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*2.0 / 8.0, p1.lat + (p2.lat - p1.lat)*2.0 / 8.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*3.0 / 8.0, p1.lat + (p2.lat - p1.lat)*3.0 / 8.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*4.0 / 8.0, p1.lat + (p2.lat - p1.lat)*4.0 / 8.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*5.0 / 8.0, p1.lat + (p2.lat - p1.lat)*5.0 / 8.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*6.0 / 8.0, p1.lat + (p2.lat - p1.lat)*6.0 / 8.0, 1000, 888888));
			t.points.push_back(SamplePoint(p1.lon + (p2.lon - p1.lon)*7.0 / 8.0, p1.lat + (p2.lat - p1.lat)*7.0 / 8.0, 1000, 888888));
			t.length += 8;
			break;
		default:
			throw("error number of part");
			break;
		}
	}
	return t;
}


Trajectory reduceTrajectory(const Trajectory* cand, int parts)
{
	Trajectory t;
	int candLen = cand->length;
	t.tid = 888888;
	t.length = 0;

	if (candLen == 1) { t.length = 1; t.points.push_back(cand->points[0]); return t;}

	for (int i = 0; i <= candLen - 1; i+=parts)
	{
		t.points.push_back(cand->points[i]);
		t.length++;
	}
	return t;
}

int main()
{

	int IfSimQuery = 0;

	int datatype =  3;

	int WriteTrajectoryToFile(string outFileName, int numTra);
	cout << "Hello world!" << endl;
	//zero-copy ÉùÃ÷
	//CUDA_CALL(cudaSetDeviceFlags(cudaDeviceMapHost));
	tradb = new Trajectory[MAX_TRAJ_SIZE]; // 50万数量级
	
	PreProcess pp;

	bool first_do = false;
	if (first_do) {
		switch (DATAINDEX) {
		case 0:
			// how to program? 构造函数
			pp = PreProcess("SH_0.txt", "SH_0_OUT.txt");
			break;
		case 1:
			pp = PreProcess("SH_1.txt", "SH_1_OUT.txt");
			break;
		case 2:

			break;
		case 3:

			break;
		case 4:

			break;
		case 5:
			
			break;
		case 6:
			
			break;
		case 7:
			
			break;
		case 8:
			
			break;
		default:
			return 0;
		}
	}
	else{
		switch (DATAINDEX) {
		case 0:
			pp.readTraFromFormatedFile("SH_0_OUT.txt");
			break;
		case 1:
			pp.readTraFromFormatedFile("SH_1_OUT.txt");
			break;
		case 2:
			pp.readTraFromFormatedFile("SH_2_OUT.txt");
			break;
		case 3:
			pp.readTraFromFormatedFile("SH_3_OUT.txt");
			break;
		case 4:
			pp.readTraFromFormatedFile("SH_4_OUT.txt");
			break;
		case 5:
			pp.readTraFromFormatedFile("SH_5_OUT.txt");
			break;
		case 6:
			pp.readTraFromFormatedFile("SH_6_OUT.txt");
			break;
		case 7:
			pp.readTraFromFormatedFile("data_SSmall_SH_OUT.txt");
			break;
		case 8:
			pp.readTraFromFormatedFile("GeoLifeOut.txt");
			break;
		case 9:
			pp.readTraFromFormatedFile("TDATA_OUT.txt");
			break;
		default:
			return 0;
		}
	}
	cout << pp.maxTid << endl;
	sleep(1);
	cout << "read trajectory success!" << endl << "Start building cell index" << endl;
	//for (int i = 1; i <= 10000;i++)
	//{
	//	printf("%d,%d\t", i, tradb[i].length);
	//}

	int cellCV = 2;
	printf("cellCV=%d\n", cellCV);
	MyTimer timer;
	timer.start();
	Grid* g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, cellCV);
	g->addDatasetToGrid(tradb, pp.maxTid);
	cout << "build Grid cell index success!" << endl;
	timer.stop();
	cout << "Grid cell time:" << timer.elapse() << endl;

	//SystemTest test(tradb, g, NULL, NULL, NULL);

	if (!IfSimQuery) {

		timer.start();
		STIG *stig = new STIG();
		stig->initial(10240, 2, tradb, pp.maxTid);
		cout << "build STIG cell index success!" << endl;
		timer.stop();
		cout << "STIG cell time:" << timer.elapse() << endl;



		timer.start();
		MortonGrid *mgrid = new MortonGrid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, cellCV);
		cout << mgrid->cellNum_axis << endl << mgrid->cellnum << endl << mgrid->cell_size << endl;
		mgrid->addDatasetToGrid(tradb, pp.maxTid);
		cout << "build MortonGrid cell index success!" << endl;
		timer.stop();
		cout << "MortonGrid cell time:" << timer.elapse() << endl;


		cout << "build cell index success!" << endl;

		bool IFFSG = 0;
		if (IFFSG) {
			timer.start();
			FSG *fsg = new FSG(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN);
			fsg->addDatasetToGrid(tradb, pp.maxTid);
			cout << "build FSG cell index success!" << endl;
			timer.stop();
			cout << "FSG cell time:" << timer.elapse() << endl;
		}


		SystemTest test(tradb, g, stig, NULL, mgrid);


		printf("\n warm up:--------------------------\n");
		//test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228 ), 80);
		//test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
		//test.MortonGridRangeQueryTestV2(MBB(121.4, 31.128, 121.44, 31.228), 80);
		//test.MortonGridRangeQueryTestV3(MBB(121.4, 31.128, 121.44, 31.228), 80);


		// only one MBB !!
		if(datatype == 1){
			bool warmup = 1;
			if(warmup){
			printf("\n test on performance:--------------------------\n");
			printf("\n test.rangeQueryTest GAT-noG\n");
			test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			printf("\n test.rangeQueryTestWithoutMorton GAT-noM \n");
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.44, 31.228), 80);
			printf("\n test.STIGrangeQueryTest\n");
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			
			if (IFFSG) {
				printf("\n test.FSGrangeQueryTest\n");
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			}

			printf("\n test.MortonGridRangeQueryTest\n");
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			}

			bool fig6 = 0;
			if(fig6){
			printf("\n test on query size-------------------------------------\n");
			printf("\n test.MortonGridRangeQueryTest\n");
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
			printf("\n test.rangeQueryTest GAT-noG\n");
			test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
			test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
			test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
			test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
			printf("\n test.rangeQueryTestWithoutMorton GAT-noM \n");
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.44, 31.228), 40);
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.44, 31.228), 60);
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.44, 31.228), 80);
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.44, 31.228), 100);
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.44, 31.228), 120);
			printf("\n test.STIGrangeQueryTest\n");
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
			if (IFFSG) {
				printf("\n test.FSGrangeQueryTest\n");
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
			}
			printf("\n\n\n\n");
			}

			bool fig9 = 0;
			if(fig9)
			{
			printf("\n fig.9 test on range size-------------------------------------\n");
			printf("\n test.MortonGridRangeQueryTest\n");
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.46, 31.228), 80);
			test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.48, 31.228), 80);
			printf("\n test.rangeQueryTest GAT-noG\n");
			test.rangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
			test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			test.rangeQueryTest(MBB(121.4, 31.128, 121.46, 31.228), 80);
			test.rangeQueryTest(MBB(121.4, 31.128, 121.48, 31.228), 80);
			printf("\n test.rangeQueryTestWithoutMorton GAT-noM \n");
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.42, 31.228), 80);
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.44, 31.228), 80);
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.46, 31.228), 80);
			test.rangeQueryTestWithoutMorton(MBB(121.4, 31.128, 121.48, 31.228), 80);
			printf("\n test.STIGrangeQueryTest\n");
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.46, 31.228), 80);
			test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.48, 31.228), 80);
			if (IFFSG) {
				printf("\n test.FSGrangeQueryTest\n");
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.46, 31.228), 80);
				test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.48, 31.228), 80);
			}
			}

			bool fig10 = 1;
			if (fig10) {
				test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
				test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
				test.MortonGridRangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
			}


			printf("\n\n\n\n");
		}
		else if (datatype == 2) {
			// MBB(116.380, 39.89, 116.480,39.93) //朝阳区 天安门
			// 116.278,39.93,116.378,39.97 // 海淀区 微软北京
			bool warmup = 1;
			if (warmup) {
				printf("\n test on performance:--------------------------\n");
				printf("\n test.rangeQueryTest GAT-noG\n");
				test.rangeQueryTest(MBB(116.278, 39.93, 116.378, 39.97), 80);
				printf("\n test.rangeQueryTestWithoutMorton GAT-noM \n");
				test.rangeQueryTestWithoutMorton(MBB(116.278, 39.93, 116.378, 39.97), 80);
				printf("\n test.STIGrangeQueryTest\n");
				test.STIGrangeQueryTest(MBB(116.278, 39.93, 116.378, 39.97), 80);
				if (IFFSG) {
					printf("\n test.FSGrangeQueryTest\n");
					test.FSGrangeQueryTest(MBB(116.278, 39.93, 116.378, 39.97), 80);
				}
				printf("\n test.MortonGridRangeQueryTest\n");
				test.MortonGridRangeQueryTest(MBB(116.278, 39.93, 116.378, 39.97), 80);
			}
		}

		else if (datatype == 3) {
			//100.00067;20.15757;119.99992;49.22045
			// MBB(116.314, 39.85, 116.414,39.89) //北京车站
			bool warmup = 1;
			if (warmup) {
				printf("\n test on performance:--------------------------\n");
				printf("\n test.rangeQueryTest GAT-noG\n");
				test.rangeQueryTest(MBB(116.314, 39.85, 116.414, 39.89), 80);
				printf("\n test.rangeQueryTestWithoutMorton GAT-noM \n");
				test.rangeQueryTestWithoutMorton(MBB(116.314, 39.85, 116.414, 39.89), 80);
				printf("\n test.STIGrangeQueryTest\n");
				test.STIGrangeQueryTest(MBB(116.314, 39.85, 116.414, 39.89), 80);

				if (IFFSG) {
					printf("\n test.FSGrangeQueryTest\n");
					test.FSGrangeQueryTest(MBB(116.314, 39.85, 116.414, 39.89), 80);
				}
				printf("\n test.MortonGridRangeQueryTest\n");
				test.MortonGridRangeQueryTest(MBB(116.314, 39.85, 116.414, 39.89), 80);
			}


		}
	}



	// printf("warm up:--------------------------\n");
	//	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);
	//	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.42, 31.228), 80);


	  //test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 20);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
	 // test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
	//  test.rangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 140);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 20);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
	//  test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
	//	test.STIGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 140);
	//	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 20);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 40);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 60);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 80);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 100);
	//  test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 120);
	//	test.FSGrangeQueryTest(MBB(121.4, 31.128, 121.44, 31.228), 140);



	//// old version
	//test.similarityQueryTest(47, 20, 5);
	//test.similarityQueryTest(19358, 10, 5);
	//test.similarityQueryTest(1299, 10, 5);
	//test.similarityQueryTest(33708, 10, 5);
	//test.similarityQueryTest(67630, 10, 5);
	//test.similarityQueryTest(99263, 10, 5);
	//test.similarityQueryTest(1309, 10, 5);
	//test.similarityQueryTest(10702, 10, 5);
	//test.similarityQueryTest(74853, 10, 5);
	//test.similarityQueryTest(104728, 10, 5);
	//test.similarityQueryTest(149443, 10, 5);
	//test.similarityQueryTest(66375, 10, 5);	
	//test.similarityQueryTest(149273, 10, 5);
	//test.similarityQueryTest(200276, 10, 5);
	//test.similarityQueryTest(229797, 10, 5);
	//test.similarityQueryTest(249730, 10, 5);
	//test.similarityQueryTest(21032, 10, 5);
	//test.similarityQueryTest(51991, 10, 5);
	//test.similarityQueryTest(69468, 10, 5);
	//test.similarityQueryTest(92757, 10, 5);


	//test.similarityQueryTest(splitTrajectory(24269, 4), 20, 15);

	// printf("test on different length:--------------------------\n"); 不同的split + 和  reduce -
	// test.similarityQueryTest(splitTrajectory(24269, 1), 20, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 2), 20, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 3), 20, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 15);
	//int TRAID = 32519;
	//int TRAID = 84199;
	//test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 1), 40, 15);
	//Trajectory d = tradb[TRAID];
	//test.similarityQueryTest(splitTrajectory(&d, 1), 40, 25);
	//test.similarityQueryTest(splitTrajectory(&d, 2), 40, 25);
	//test.similarityQueryTest(splitTrajectory(&d, 3), 40, 25);
	//test.similarityQueryTest(splitTrajectory(&d, 4), 40, 25);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// test.similarityQueryTest(splitTrajectory(&d, 1), 40, 15);
	// test.similarityQueryTest(splitTrajectory(&d, 2), 40, 15);
	// test.similarityQueryTest(splitTrajectory(&d, 3), 40, 15);
	// test.similarityQueryTest(splitTrajectory(&d, 4), 40, 15);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 40, 15);
	// test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 2), 40, 15);
	// test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 4), 40, 15);
	// test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 6), 40, 15);
	// test.similarityQueryTest(reduceTrajectory(&tradb[TRAID], 8), 40, 15);



	 //test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 10, 5);
	 //test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 20, 5);
	 //test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 30, 5);
	 //test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 40, 5);
	 //test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 50, 5);
	 //test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 60, 5);
	 //test.similarityQueryTest(reduceTrajectory(&tradb[188], 1), 70, 5);



	if(IfSimQuery){

		SystemTest test(tradb, g, NULL, NULL, NULL);
		int TRAID = 2; 
		//int TRAID = 84199; // 995 length good for GPU acceleration

		int stream_size = 60;
		int batch_size0 = stream_size;
		int maxtra = pp.maxTid-10;
		for (int i = 0; i < stream_size;i+= batch_size0) {
			int tmp = stream_size - i;
			int batch_size = tmp > batch_size0 ? batch_size0 : tmp;
			//Trajectory* tarray = new Trajectory[batch_size];
			vector<Trajectory> tvec;
			vector<Trajectory> tvec2;
			for (int j = 0; j < batch_size;j++){
				//srand((unsigned)time(NULL));
				//int TRAID = rand()%400;

				//TRAID += 5000;
				TRAID += int(maxtra/ batch_size0);
				TRAID %= maxtra;
	
				//TRAID = (rand() / double(RAND_MAX)) * maxtra;
				// tradb[TRAID].PrintTraj();
				Trajectory d = reduceTrajectory(&tradb[TRAID], 4);
				//d.PrintTraj();
				Trajectory d2 = splitTrajectory(&d, 5);
				d2.PrintTraj();
				//tarray[j] = d2;
				tvec.push_back(d2);
				tvec2.push_back(d2);
			}
			//for (int i = 0; i < tvec.size(); i++)
			//	tvec.at(i).PrintTraj();
			std::sort(tvec.begin(), tvec.end(), nonmem_cmp);

			//for (int i = 0; i < tvec.size(); i++)
			//	tvec.at(i).PrintTraj();
			//cout << endl;

			//test.similarityQueryTest2(tarray , batch_size, 25); 

			// 25 for fig.12 but there are some extra results.

			test.similarityQueryTest3(tvec, batch_size, 25);    // 最终版本 取代了V1
	
			//delete[] tarray;
			}
		}

/*
	// int TRAID2 = 84199;
	int TRAID2 = 8;
	Trajectory d = reduceTrajectory(&tradb[TRAID2], 4);
	test.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	*/

	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 5);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 10);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 15);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 20);
	//test.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);

	// test.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 60, 25);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 80, 25);
	// test.similarityQueryTest(splitTrajectory(&d, 5), 100, 25);

	// old version
	// test.similarityQueryTest(splitTrajectory(21032, 1), 20, 15);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 30, 15);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 40, 15);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 50, 15);
	// test.similarityQueryTest(splitTrajectory(21032, 1), 60, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 5);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 10);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 15);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 20);
	// test.similarityQueryTest(splitTrajectory(24269, 4), 20, 25);
	 // test.similarityQueryTest(21032, 30, 5);
	 // test.similarityQueryTest(21032, 40, 5);
	 // test.similarityQueryTest(21032, 50, 5);
	 // test.similarityQueryTest(21032, 60, 5);
	 //test.similarityQueryTest(47, 10, 5);
	 //test.similarityQueryTest(47, 20, 5);
	 //test.similarityQueryTest(47, 30, 5);
	 //test.similarityQueryTest(47, 40, 5);
	 //test.similarityQueryTest(47, 50, 5);
	 //test.similarityQueryTest(47, 60, 5);
	 //test.similarityQueryTest(47, 70, 5);
	 //test.similarityQueryTest(47, 20, 5);
	 //test.similarityQueryTest(47, 20, 10);
	 //test.similarityQueryTest(47, 20, 15);
	 //test.similarityQueryTest(47, 20, 20);
	 //test.similarityQueryTest(47, 20, 25);
	 //test.similarityQueryTest(47, 20, 30);
	 //test.similarityQueryTest(47, 20, 35);

	printf("test on different CV:--------------------------\n");
	// printf("cellCV=0\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 0);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test1(tradb, g, stig, fsg);
	// test1.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;
	// printf("cellCV=2\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 2);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test2(tradb, g, stig, fsg);
	// test2.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;
	// printf("cellCV=4\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 4);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test3(tradb, g, stig, fsg);
	// test3.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;
	// printf("cellCV=6\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 6);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test4(tradb, g, stig, fsg);
	// test4.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;
	// printf("cellCV=8\n");
	// g = new Grid(MBB(pp.xmin, pp.ymin, pp.xmax, pp.ymax), CELL_LEN, 8);
	// g->addDatasetToGrid(tradb, pp.maxTid);
	// SystemTest test5(tradb, g, stig, fsg);
	// test5.similarityQueryTest(splitTrajectory(&d, 5), 40, 25);
	// delete g;

	printf("Finished.\n");


	//CPURangeQueryResult* resultTable = NULL;
	//int RangeQueryResultSize = 0;
	//MBB mbbArray[1000];
	//int* resultSize = NULL;
	//for (int i = 0; i <= 999; i++)
	//	mbbArray[i] = MBB(121.1, 31.1, 121.3, 31.3);
	//MyTimer timer;
	//timer.start();
	//g->rangeQueryBatch(mbbArray, 1000, resultTable, resultSize);
	//timer.stop();
	//cout << "CPU Time:" << timer.elapse() << "ms" << endl;


	//CUDA_CALL(cudaMalloc((void**)(&baseAddrGPU), 512 * 1024 * 1024));
	//void* baseAddr = baseAddrGPU;
	//timer.start();
	//g->rangeQueryBatchGPU(mbbArray, 1000, resultTable, resultSize);
	//timer.stop();
	//cout << "GPU Time:" << timer.elapse() << "ms" << endl;
	//CUDA_CALL(cudaFree(baseAddr));
	//baseAddrGPU = NULL;
	//Trajectory* qTra = new Trajectory[100];
	//for (int i = 0; i <= 99; i++)
	//{
	//	qTra[i] = tradb[47]; // length is 1024
	//}
	////for (int i = 1; i <= 9999;i++)
	////{
	////	if (tradb[i].length > 600)
	////		printf("tra:%d,length:%d\n", i, tradb[i].length);
	////}


	////Similarity on CPU
	//int* simiResult = new int[10 * 100];
	//g->SimilarityQueryBatch(qTra, 2, simiResult,50);
	//for (int i = 0; i <= 1; i++) {
	//	cout << "Trajectory:" << i << endl;
	//	for (int j = 40; j <= 49; j++) {
	//		cout << simiResult[i * 50 + j] << "\t" << endl;
	//	}
	//}
	//delete[] simiResult;


	////Similarity on GPU
	//simiResult = new int[10 * 100];
	//g->SimilarityQueryBatchOnGPU(qTra, 2, simiResult, 50);
	//for (int i = 0; i <= 1; i++)
	//{
	//	cout << "Trajectory:" << i << endl;
	//	for (int j = 40; j <= 49; j++)
	//	{
	//		cout << simiResult[i * 50 + j] << "\t" << endl;
	//	}
	//}


	//g->rangeQuery(MBB(121.4, 31.15, 121.6, 31.25), resultTable, &RangeQueryResultSize);
	//g->rangeQueryGPU(MBB(121.4, 31.15, 121.6, 31.25), resultTable, &RangeQueryResultSize);


	//Trajectory **testTra = (Trajectory**)malloc(sizeof(Trajectory*) * 5000);
	//for (int i = 2; i <= 5001; i++) {
	//	testTra[i-2] = &tradb[i];
	//}
	//float *EDRdistance = (float*)malloc(sizeof(float) * 5000);
	//g->SimilarityQuery(tradb[2], testTra, 5000, EDRdistance);


	getchar();
	getchar();
	getchar();

	return 0;
}

int WriteTrajectoryToFile(string outFileName, int numTra)
{
	ofstream fout;
	fout.open(outFileName, ios_base::out);
	for (int i = 1; i <= numTra; i++)
	{
		fout << i << ": ";
		for (int j = 0; j <= tradb[i].length - 1; j++)
		{
			fout << tradb[i].points[j].lon << "," << tradb[i].points[j].lat << ";";
		}
		fout << endl;
	}
	fout.close();
	return 1;
}
