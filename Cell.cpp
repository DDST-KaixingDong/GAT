#include "Cell.h"

using namespace std;

Cell::Cell()
{
	cell_x = 0;
	cell_y = 0;
	mbb = MBB(0, 0, 0, 0);
	subTraNum = 0;
	subTraPtr = NULL;
}

Cell::Cell(int x, int y, const MBB& val_mbb) {
	cell_x = x;
	cell_y = y;
	mbb = val_mbb;
	subTraNum = 0;
	subTraPtr = NULL;
}

bool Cell::initial(int x, int y, const MBB& val_mbb) {
	cell_x = x;
	cell_y = y;
	mbb = val_mbb;
	anchorPointX = int((mbb.xmin + mbb.xmax) / 2 * 1000000);
	anchorPointY = int((mbb.ymin + mbb.ymax) / 2 * 1000000);
	subTraNum = 0;
	totalPointNum = 0;
	subTraEntry.next = NULL;
	subTraPtr = &subTraEntry;
	subTraTable = NULL;
	return true;
}


// ��������ӹ켣
int Cell::addSubTra(int traID, int startIdx, int endIdx, int numOfPoints)
{
	subTra* newSubTra = (subTra*)malloc(sizeof(subTra));//�ǵ�free������
	if (newSubTra == NULL)
		return 1;
	newSubTra->numOfPoint = numOfPoints;
	newSubTra->traID = traID;
	newSubTra->startpID = startIdx;
	newSubTra->endpID = endIdx;
	newSubTra->next = NULL;
	subTraPtr->next = newSubTra;
	subTraPtr = newSubTra;
	subTraNum++;
	return 0;
}


int Cell::buildSubTraTable()
{
	int countPoints = 0;
	subTraTable = (subTra*)malloc(sizeof(subTra)*subTraNum);
	if (subTraTable == NULL)
		return 1;
	subTra* ptr = subTraEntry.next;
	subTra* nextptr;
	int idx = 0;

	while (ptr!= NULL) {
		nextptr = ptr->next;
		subTraTable[idx] = (*ptr);
		countPoints += (*ptr).numOfPoint;
		free(ptr);
		ptr = nextptr;
		idx++;
	}

	if (idx != subTraNum) //debug info
	{
		cout << idx << "vs" << endl;
		getchar();
	}

	// throw("error in cell buildSubTraTable");

	this->totalPointNum = countPoints;
	return 0;

}

int Cell::writeCellToFile(string fileName)
{
	return 0;
}



Cell::~Cell()
{
}
