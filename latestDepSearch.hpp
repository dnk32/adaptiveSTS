//latestDepSearch class
#ifndef LATESTDEPSEARCH_HPP
#define LATESTDEPSEARCH_HPP

#include <vector>
#include <string>
#include "classDefs.h" // graph node class and heap class

using namespace std;
typedef vector< vector< vector< double > > > dblVec3D; 
class latestDepSearch{
    
    private:

    /*==================
     * member variables
     *==================*/
    // environment boundaries
    int xmax, xmin, ymax, ymin, tmax, tmin;
    // grid spacings
    double xSpGrid, ySpGrid;
    double tSpGrid;
    // number of discrete points in grid
    int nX, nY;

    //spacing and number of points for data vectors
    double xSpData, ySpData, tSpData;
    int nDx, nDy, nDt;
    bool dataParamsSet;

    // pointers to data vectors
    dblVec3D *UxVec;
    dblVec3D *VyVec;
    dblVec3D *OBS;

    //maximum vehicle and flow velocities
    double Vm, Vfm;

    int lim = 3;
    
    // hashbin settings for optimal time search
    int hashBinCap = 10;
    //int hashBinHeight;          // number of seconds inside a hashbin
    int nHashBins;              // number of hash bins
    int hashBinHeight;
    
    // coordinates to start search
    double startx, starty, startt;
    int START_COORD[3];

    /*====================
     * member functions
     *====================*/
    int (*getHashBinNumber)(int, int, double, bool, int, int, int);
    void getFlowVelFromVecs(int x,int y,double t,double &vx,double &vy);
    double getOptTimeCost(graphNode* currNode, int nx, int ny, double Vf, double thF);
    bool isAccessible(int nx,int ny,double t);
    vector< graphNode* > getPlannedPath(graphNode* currNodePtr);
    //void runSearch();
    
    public:
    latestDepSearch():dataParamsSet(false){}
    void setEnvLims(int maxX, int minX, int maxY, int minY, int maxT, int minT);
    void setGridParams(double xGrid, double yGrid, double tGrid, int nXgrid, int nYgrid);
    void setHashBinParams(int nHBins,int hBinHeight, int (*hashFnc)(int, int, double, bool, int, int, int) );
    void setDataVecs(double xData, double yData, double tData, int nXdata, int nYdata, int nTdata , dblVec3D *uVec, dblVec3D *vVec, dblVec3D *obsVec);
    void setVelParams(double vMaxVeh, double vMaxFlow);
    void setStartCoords(double xSt, double ySt);
    vector <double> getEarliestDepTimes();
    void runSearch( vector <double> *depTimeVec);
};

#endif
