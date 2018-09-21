#ifndef PREPROCESS_H
#define PREPROCESS_H
#include<fstream>
#include "ConstDefine.h"


typedef struct VLLT
{
	std::string vid;
    float lon;
    float lat;
    int time;
}VLLT;

struct tidLinkTable
{
    int tid;
    tidLinkTable* next = NULL;
};

class PreProcess
{
    public:
		std::ifstream fin;
        std::ofstream fout;
        int maxTid = 0;
        PreProcess();
        PreProcess(std::string fileName, std::string outFileName);
		bool updateMapBound(float lon,float lat);
		bool validPoint(float lon, float lat);
        virtual ~PreProcess();
        VLLT getTraInfoFromString(std::string s);
		bool writeTraDataToFile(std::string outFileName);
		bool readTraFromFormatedFile(std::string outFileName);

		PreProcess & operator=(const PreProcess & pp);

		float xmin, xmax, ymin, ymax; // ±ß½ç


    protected:

    private:
};

#endif // PREPROCESS_H
