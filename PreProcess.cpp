#include "PreProcess.h"
#include "Trajectory.h"
#include "SamplePoint.h"
#include "ConstDefine.h"
#include <map>
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<cstdlib>

using namespace std;


extern int DaysBetween2Date(string date1, string date2);
//�켣��Ŵ�1��ʼ
extern Trajectory* tradb;
extern map<string, tidLinkTable*> vidTotid;
extern map<string, tidLinkTable*>::iterator iter;
extern string baseDate;


PreProcess::PreProcess()
{
    //ctor
}

// split��������
// �򵥵ķָ��
void split(std::string& s, const std::string& delim,std::vector< std::string >* ret)
{
    size_t last = 0;
    size_t index=s.find_first_of(delim,last);
    while (index!=std::string::npos)
    {
        ret->push_back(s.substr(last,index-last));
        last=index+1;
        index=s.find_first_of(delim,last);
    }
    if (index-last>0)
    {
        ret->push_back(s.substr(last,index-last));
    }
}



PreProcess::PreProcess(string fileName,string outFileName)
{
	xmin = 180;
	xmax = 0;
	ymin = 90;
	ymax = 0;
    fin.open(fileName,ios_base::in);
    fout.open(outFileName,ios_base::out); // not used

	// ��ȡtxt�ļ�
    string buffer;
    buffer.assign(istreambuf_iterator<char>(fin),istreambuf_iterator<char>());
    stringstream bufferstream;
    bufferstream.str(buffer);
    string linestr;

	//MyTimer timer;
	//tradb = (Trajectory*)malloc(sizeof(Trajectory) * 100000);
	
    while(getline(bufferstream,linestr))
    {
        string s = linestr;
//		timer.start();
		VLLT vllt = getTraInfoFromString(s);
		if (!validPoint(vllt.lon, vllt.lat))
			continue;
		if (vllt.vid == "0000")
			continue;
		updateMapBound(vllt.lon, vllt.lat);
//		timer.stop();
//		printf("vllt time: %lf\n", timer.elapse());

        iter = vidTotid.find(vllt.vid);
        // maxTid�ǵ�ǰ������Ĺ켣��
        int nowTid = this->maxTid;
        
//		timer.start();
        if(iter!=vidTotid.end())
        {
            tidLinkTable* table = iter->second;
            while(table->next!=NULL)
            {
                table = table->next;
            }
            nowTid = table->tid;

        }
        //�µ�vid
        else
        {
            //�����µĹ켣��ţ���Ϊ�ȷ��䣬���Թ켣��Ŵ�1��ʼ
            this->maxTid++;
			nowTid = this->maxTid;
            tidLinkTable* tidNode = new tidLinkTable();
            tidNode->tid = nowTid;
            tidNode->next = NULL;
            vidTotid.insert(pair<string,tidLinkTable*>(vllt.vid,tidNode));

			//�����¹켣��ͷ���vid��tid
            tradb[nowTid].vid = vllt.vid;
			tradb[nowTid].tid = nowTid;

        }

//		timer.stop();
//		printf("find tid time: %lf\n", timer.elapse());


        //�ҳ��켣�ź���ӹ켣
//		timer.start();
        int ret = tradb[nowTid].addSamplePoints(vllt.lon,vllt.lat,vllt.time);
		
        if(ret == 0)
        {
//			timer.stop();
//			printf("add point time: %lf\n", timer.elapse());
			continue;
        }
        else if((ret == 1) || (ret == 2))
        {
			
            tidLinkTable* node = vidTotid.find(vllt.vid)->second;
            while(node->next!=NULL) node = node->next;
            tidLinkTable* newNode = new tidLinkTable();
            newNode->next = NULL;
			node->next = newNode;
	
			//�������� �µ�vid�в���
            this->maxTid++;
			nowTid = this->maxTid;
            newNode->tid = nowTid;
			if (tradb[nowTid].addSamplePoints(vllt.lon, vllt.lat, vllt.time)!=0)//һ���ܳɹ�
				throw("error");
            tradb[nowTid].vid = vllt.vid;
			tradb[nowTid].tid = nowTid;
        }
        //̫Զ�������õ�
        //��Ϊ���޸���һ������Ǵ�����������ȼ�¼�õ㣬Ȼ����һ�������������������۹�10�����������һ������Ǵ�ģ��Ⱥ�����������
        else if(ret == 3)
        {
//            tradb[nowTid].errPointBuff[tradb[nowTid].errCounter] = SamplePoint(vllt.lon,vllt.lat,vllt.time,nowTid);
//            tradb[nowTid].errCounter++;
//            if(tradb[nowTid].errCounter>=10)
//            {
//
//            }
            continue;
        }

    }
}

inline bool PreProcess::updateMapBound(float lon,float lat)
{
	if (lat < ymin)
		ymin = lat;
	if (lat > ymax)
		ymax = lat;
	if (lon < xmin)
		xmin = lon;
	if (lon > xmax)
		xmax = lon;
	return true;
}

bool PreProcess::validPoint(float lon, float lat)
{
	if (lat < 25)
		return false;
	if (lat > 50)
		return false;
	if (lon < 90)
		return false;
	if (lon > 130)
		return false;
	return true;
}

VLLT PreProcess::getTraInfoFromString(string s)
//get Vid Longitude Latitude TimeStamp from string s
{
    VLLT vllt;
    vector<string> partOfLine;
    string dot = ",";
    split(s,dot,&partOfLine);
	
	// ��������
	if (partOfLine.size() < 12) {
		VLLT a;
		a.vid = "0000";
		return a;
	}

    vllt.vid = partOfLine[2];
    vllt.lon = atof(partOfLine[3].c_str());
    vllt.lat = atof(partOfLine[4].c_str());

	vector<string> datetime;
	split(partOfLine[8], string(" ") , &datetime);
	vector<string> timeValue;
	split(datetime[1], string(":"), &timeValue);
	int days = DaysBetween2Date(baseDate, datetime[0]);
	int hours = atoi(timeValue[0].c_str());
	int minute = atoi(timeValue[1].c_str());
	int second = atoi(timeValue[2].c_str());
	vllt.time = days * 86400 + hours * 3600 + minute * 60 + second;

    return vllt;
}

bool PreProcess::writeTraDataToFile(string outFileName)
{
	FILE* fd = fopen(outFileName.c_str(), "w+");
	fprintf(fd, "%d\n", this->maxTid);
	fprintf(fd, "%f,%f,%f,%f\n", this->xmin, this->ymin, this->xmax, this->ymax);
	for(int tid = 1; tid <= maxTid; tid++)
	{
		fprintf(fd, "%d;", tradb[tid].length);
	}
	fprintf(fd, "\n");
	for (int tid = 1; tid <= maxTid;tid ++)
	{
		for (int idx = 0; idx <= tradb[tid].length - 1;idx++)
		{
			fprintf(fd, "%f,%f,%d,%d;", tradb[tid].points[idx].lat, tradb[tid].points[idx].lon, tradb[tid].points[idx].tid, tradb[tid].points[idx].time);
		}
		fprintf(fd, "\n");
	}
	fclose(fd);
	return true;
}

bool PreProcess::readTraFromFormatedFile(string outFileName)
{
	fin.open(outFileName, ios_base::in);

	// ��ȡtxt�ļ�
	string buffer;
	buffer.assign(istreambuf_iterator<char>(fin), istreambuf_iterator<char>());
	stringstream bufferstream;
	bufferstream.str(buffer);
	string linestr;

	//MyTimer timer;
	//tradb = (Trajectory*)malloc(sizeof(Trajectory) * 100000);

	int cnt = 0;
	char* nextStart = NULL;

	while (getline(bufferstream, linestr))
	{
		string s = linestr;
		if (s.empty()) continue;

		if(cnt==0)
		{
			this->maxTid = std::atoi(s.c_str());
		}
		else if(cnt==1)
		{
			const char* str_c = s.c_str();
			char* thisLineStr = (char*)malloc(s.length());
			memcpy(thisLineStr, str_c, s.length());
			char* num = strtok(thisLineStr, ",");
			this->xmin = std::atof(num);
			num = strtok(NULL, ",");
			this->ymin = std::atof(num);
			num = strtok(NULL, ",");
			this->xmax = std::atof(num);
			num = strtok(NULL, ",");
			this->ymax = std::atof(num);
			free(thisLineStr);
		}
		else if(cnt ==2) 
		{
			const char* str_c = s.c_str();
			char* thisLineStr = (char*)malloc(s.length());
			memcpy(thisLineStr, str_c, s.length());
			char* num = strtok(thisLineStr, ";");
			tradb[0].length = 0;
			tradb[1].length = std::atoi(num);
			tradb[1].points.resize(tradb[1].length);
			for (int i = 2; i <= this->maxTid;i++)
			{
				num = strtok(NULL, ";");
				tradb[i].length = std::atoi(num);
				tradb[i].points.resize(tradb[i].length);
			}
			free(thisLineStr);
		}
		else
		{
			// ���Ϊcnt - 2 �Ĺ켣 ��1 ��ʼ����
			tradb[cnt - 2].tid = cnt - 2;
			const char* str_c = s.c_str();
			char* thisLineStr = (char*)malloc(s.length());
			memcpy(thisLineStr, str_c, s.length());
			if (tradb[cnt - 2].length != 0) {
				nextStart = strchr(thisLineStr, ';');
				char* nums = strtok(thisLineStr, ";");
				char* num = strtok(nums, ",");
				tradb[cnt - 2].points[0].lat = std::atof(num);
				num = strtok(NULL, ",");
				tradb[cnt - 2].points[0].lon = std::atof(num);
				num = strtok(NULL, ",");
				tradb[cnt - 2].points[0].tid = std::atoi(num);
				num = strtok(NULL, ",");
				tradb[cnt - 2].points[0].time = std::atoi(num);
				
				for (int idx = 1; idx <= tradb[cnt - 2].length - 1;idx++)
				{
					nums = nextStart + 1;
					nextStart = strchr(nums, ';');
					num = strtok(nums, ",");
					tradb[cnt - 2].points[idx].lat = std::atof(num);
					num = strtok(NULL, ",");
					tradb[cnt - 2].points[idx].lon = std::atof(num);
					num = strtok(NULL, ",");
					tradb[cnt - 2].points[idx].tid = std::atoi(num);
					num = strtok(NULL, ";");
					tradb[cnt - 2].points[idx].time = std::atoi(num);
				}

			}
			free(thisLineStr);
		}
		cnt++;
	}
	return false;
}

PreProcess::~PreProcess()
{
    //dtor
}

PreProcess& PreProcess::operator=(const PreProcess& pp) {

	if (this == &pp) {
		return *this;
	}
	this->maxTid = pp.maxTid;
	this->xmax = pp.xmax;
	this->xmin = pp.xmax;
	this->ymax = pp.ymax;
	this->ymin = pp.ymin;

}