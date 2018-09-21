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
	int traID;			// tra��ID
	short startpID;		// �ڹ켣traID�е�startpID
	short endpID;		// �ڹ켣traID�е�endpID
	short numOfPoint;	// traID���ӹ켣�ĵ�ĸ���
	int idxInAllPointsArray;// ��ʼ����AllPointsArray��id
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
	int cell_x;	//cell������ �к�
	int cell_y;	//cell������ �к�
	MBB mbb;	//cell�ķ�Χrange

	int anchorPointX; //DeltaEncoding�Ļ�׼��Ŀǰ���������cell����������  �������� ��monton�������
	int anchorPointY;

	int subTraNum; //�ӹ켣�ĸ���
	int totalPointNum; //cell�ڵ����

	// ���ñ���
	subTra subTraEntry;	//	����cell������������ڣ���һ��ĵ���̫�˷��ڴ棩
	subTra* subTraPtr;	//	��ǰ��������ָ��λ��
	subTra* subTraTable;//	ת��Ϊ�������������

	std::ofstream fout;//�ļ��ӿ�
	int writeCellToFile(std::string fileName);

#ifdef _CELL_BASED_STORAGE
	// ÿ��cell��point���������ݵ���ʼ����ֹλ��
	int pointRangeStart;
	int pointRangeEnd;
#endif // _CELL_BASED_STORAGE

};

