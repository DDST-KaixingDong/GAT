#include "Trajectory.h"
#include <iostream>
using namespace std;
extern float calculateDistance(float LatA,float LonA,float LatB,float LonB);

Trajectory::Trajectory()
{

    //ctor
}

Trajectory::Trajectory(int tid,std::string vid)
{
    this->tid = tid;
    this->vid = vid;
    this->length = 0;
}


int Trajectory::addSamplePoints(float lon,float lat,int time)
{
    if(this->length>=MAXLENGTH)
    {
        return 1;
    }
    if(this->length>0)
    {
        if(time - this->points[this->length-1].time > MAXGAP)
        {
            return 2;
        }
        if((calculateDistance(lat,lon,this->points[this->length-1].lat,this->points[this->length-1].lon))/(time - this->points[this->length-1].time)>=50)
        {
            return 3;
        }
    }
    //经过检查可以加入这点到轨迹中
	this->points.push_back(SamplePoint(lon, lat, time, this->tid));
    //this->points[this->length].lat = lat;
    //this->points[this->length].lon = lon;
    //this->points[this->length].tid = this->tid;
    //this->points[this->length].time = time;
    this->length++;
    return 0;
}



Trajectory::~Trajectory()
{
    //dtor
}


void Trajectory::PrintTraj() {
	cout << "tid: " << this->tid << "length: " << this->length << endl;
}