#ifndef TRAJECTORY_H
#define TRAJECTORY_H
#include "ConstDefine.h"
#include "SamplePoint.h"



class Trajectory
{
    public:
        Trajectory();
        Trajectory(int tid,std::string vid);
        int addSamplePoints(float lon,float lat,int time);
        virtual ~Trajectory();

		std::vector<SamplePoint> points;


        int tid;
		std::string vid;

        int length = 0;       
        int errCounter = 0;
        SamplePoint errPointBuff[10];

		void PrintTraj();

    protected:

    private:
};

#endif // TRAJECTORY_H
