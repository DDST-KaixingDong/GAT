#pragma once
#include "ConstDefine.h"
#include "MBB.h"
#include <string>
#include <iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
// for test



typedef struct subTra {
	int traID;			// tra的ID
	short startpID;		// 在轨迹traID中的startpID
	short endpID;		// 在轨迹traID中的endpID
	short numOfPoint;	// traID的子轨迹的点的个数
	int idxInAllPointsArray;// 起始点在AllPointsArray的id
	subTra* next;
}subTra;

class Cell
{
public:
	Cell();
	~Cell();
	Cell(int x, int y,const MBB& val_mbb);

	bool initial(int x, int y, const MBB& val_mbb);
	int addSubTra(int traID, int startIdx, int endIdx, int numOfPoints);
	int buildSubTraTable();
	int cell_x;	//cell横坐标 行号
	int cell_y;	//cell纵坐标 列号
	MBB mbb;	//cell的范围range

	int anchorPointX; //DeltaEncoding的基准，目前想的是整个cell的中心坐标  增量编码 与monton编码相对
	int anchorPointY;

	int subTraNum; //子轨迹的个数
	int totalPointNum; //cell内点个数

	// 常用变量
	subTra subTraEntry;	//	建立cell过程中链表入口（这一点改掉，太浪费内存）
	subTra* subTraPtr;	//	当前最新数据指针位置
	subTra* subTraTable;//	转化为数组后数组的入口

	std::ofstream fout;//文件接口
	int writeCellToFile(std::string fileName);

#ifdef _CELL_BASED_STORAGE
	// 每个cell在point数组中数据的起始和终止位置
	int pointRangeStart;
	int pointRangeEnd;
#endif // _CELL_BASED_STORAGE

};

