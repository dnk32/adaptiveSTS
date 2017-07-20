#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>  //including Eigen library for matrix manipulation
#include <gsl/gsl_spline.h>
#include <chrono>
#include <thread>
#include <mutex>

#include "classDefs.h"
#include "latestDepSearch.hpp"

using namespace std;
using namespace Eigen;

//===========================
// defines for the project
//===========================
#define dxType 2 // if dxType = 1 dxAllowed and dtAllowed computed together
                 // if dxType = 2 dxAllowed and dtAllowed computed seperately

//========================================
// Environment description for flow velocity file
//========================================

const int xSpData = 500;      // x spacing between two grid points
const int ySpData = 500;      // y spacing between two grid points
const int tSpData = 1000;     // time spacing between two grid points in seconds

const int nX = 233; // number of points in the x-direction
const int nY = 120; // number of points in the y-direction
const int nT = 519; // number of points in the t-direction

//========================================
// define parameters related to flow data
//========================================

vector< vector< vector< double > > > vXVec;
vector< vector< vector< double > > > vYVec;

//================================
// obstacle definitions
//================================
vector< vector< vector< double > > > OBS;
// Rectangles: {x1, y1, x2, y2}, x1<x2, y1<y2

//================================
// search space boundaries
//================================

//// limits used for xSt = [10,25], xGl = [40,10];
//int xmax = 42000;
//int xmin = 5000;
//int ymax = 27000;
//int ymin = 8000;
//int tmax = 250000;
//int tmin = 0;

//// limits used for xSt = [20,50], xGl = [25,45];
//int xmax = 28000;
//int xmin = 14000;
//int ymax = 52000;
//int ymin = 43000;
//int tmax = 100000;
//int tmin = 0;

//// limits used for xSt = [20,50], xGl = [25,45];
//int xmax = 22500;
//int xmin = 15500;
//int ymax = 50000;
//int ymin = 45200;
//int tmax = 79000;
//int tmin = 0;

// limits used for xSt = [20,50], xGl = [50,40];
int xmax = 52000;
int xmin = 14000;
int ymax = 52000;
int ymin = 30000;
int tmax = 185000;
int tmin = 0;

//===========================
// Algorithm parameters
//===========================

// maximum flow speed encountered
double Vfm = 0.7288;  // m/s or mm/ms

//max vehicle speed allowed
double Vm = 0.5; // m/s or mm/ms

// discrete time layers considered
int dTLayer = 200;              // all time dimensions will be  multiples of dTlayer.
                                // this is done to speed up the search
// minimum internode distance
double xEpsMin = 50;            // smaller values will lead to large node expansions

// number of runs
const int nTotRuns = 1;

// neighbor selection params for all the runs
int nDivsAll[nTotRuns] = {3};

// neigbbor selection params
int nDivs;

// flow Velocity threshold
double p = 0.1;

//cost function parameters
double k1 = 0.0005;
double k2 = 1.0;
int alpha = 2;

// termination epsilon
int xEndEps =100;

// min vels to compute heuristics
double vMinCost = -Vfm + sqrt(Vfm*Vfm + k1/k2); // vehicle speed that gives min cost heuristic
double vMinH = (Vm<vMinCost)? Vm:vMinCost;       // select the smaller of vMinCost and Vmax

//===========================
// Hash Function Parameters
//===========================
// number of big bins the space is divided to

/* second hash function is used
 * *****************************/

// the dimensions of each hashbin
int dxBin = xEpsMin*8;//250;
int dyBin = xEpsMin*8;//250;
int dtBin = dTLayer*3;//500;

// number of hashBins along each axis
int nXBins = 100;
int nYBins = 100;
int nTBins = (tmax-tmin)/dtBin;//40;

// size of the complete hashBin
int xSpBin = nXBins*dxBin;
int ySpBin = nYBins*dyBin;
int tSpBin = nTBins*dtBin;

// number of total hashBins
int nHashBins = nXBins*nYBins*nTBins;

// number of nodes in each hashbin
int hashBinSize = 200;

//========================
// max path lengths
//========================
int maxPathLength = ceil( (tmax-tmin)/(double)dTLayer );
int neighbPathLength = 0;

//================================
// start and end goal definitions
//================================
// Start and end
int startx = 20000;
int starty = 50000;
int startt = 0;
int endx = 50000;
int endy = 40000;

int START_COORD[3] = {startx, starty, round( (double)startt/dTLayer )*dTLayer};
int GOAL_COORD[2] = {endx, endy};

//===========================
// Include helper functions
//===========================
#include "helperFct.h"
// the helperFct.h file contains the funcitons required by the algorithm.
// The file is inluded at this point because the fucntions use the variables defined above.


//========================================================================
// Parameters to be used with latest departure time search
//========================================================================
// Use same environment boundaries, same start and end coordinates, same 
// velocity and obstacle vectors, same vehicle constraints.
// New definitions required for spatial and temporal discretization levels,
// number of grid points and hashbin function.

    // number of points in the optimal time search grid
    int nXOptT, nYOptT;

// hash function to use with latest departure time search
//**********************************************************************
int getHashBinNumberOptT(int x, int y, double t, bool isVmSmall, int nXT, int nYT, int hashBinHeight){
    if (isVmSmall)
        return (int) ( (x)*nYT + y + nXT*nYT*( (int)t/hashBinHeight ) ); 
    else
        return (int) (x*nYT + y);
}

// function to compute latest departure time for nodes in a grid
//***************************************************************
vector <double>  getDepTimeVect(double xSpGT, double ySpGT, double tSpGT){
    
    // grid spacings
    double xSpGridOptT = xSpGT; 
    double ySpGridOptT = ySpGT; 
    double tSpGridOptT = tSpGT; 

    // number of discrete points in grid
    nXOptT = (int) floor( ( xmax - xmin )/xSpGridOptT + 1 );
    nYOptT = (int) floor( ( ymax - ymin )/ySpGridOptT + 1 );

    // maximum vehicle and flow velocities
    bool isVmSmall = (Vm<=Vfm);
    

    // Hashbin data
    int nHashBinsOptT;
    int hashBinHeight = tSpGridOptT*10;
    if (isVmSmall)
        nHashBinsOptT = nXOptT*nYOptT*( (int)tmax/hashBinHeight + 1);         
    else
        nHashBinsOptT = nXOptT*nYOptT;          

    // vector for latest departure time data
    vector <double> depTimeVec;

    // create latest departure time object and populate parameters
    latestDepSearch dptTimeSearch;
    dptTimeSearch.setEnvLims(xmax,xmin, ymax, ymin, tmax, tmin);
    dptTimeSearch.setGridParams(xSpGridOptT, ySpGridOptT, tSpGridOptT, nXOptT, nYOptT);
    dptTimeSearch.setHashBinParams(nHashBinsOptT, hashBinHeight, &getHashBinNumberOptT);
    dptTimeSearch.setDataVecs(xSpData, ySpData, tSpData, nX, nY, nT, &vXVec, &vYVec, &OBS);
    dptTimeSearch.setVelParams(Vm, Vfm);
    dptTimeSearch.setStartCoords(endx, endy);
    dptTimeSearch.runSearch( &depTimeVec );
    
    return depTimeVec;
}


// function to compute latest departure time from a given coordinate
//*******************************************************************

double getLatestDepTime(int nx, int ny, vector<double> *depTimeVec, double xSpGT, double ySpGT ){
    int xCoord = (nx-xmin)/xSpGT;
    int yCoord = (ny-ymin)/ySpGT;

    double dx = (nx - xCoord*xSpGT)/xSpGT;
    double dy = (ny - yCoord*ySpGT)/ySpGT;

    int hB00 = getHashBinNumberOptT(xCoord, yCoord, 0, false, nXOptT, nYOptT, 0); 
    int hB10 = getHashBinNumberOptT(xCoord+1, yCoord, 0, false, nXOptT, nYOptT, 0); 
    int hB01 = getHashBinNumberOptT(xCoord, yCoord+1, 0, false, nXOptT, nYOptT, 0); 
    int hB11 = getHashBinNumberOptT(xCoord+1, yCoord+1, 0, false, nXOptT, nYOptT, 0); 

    return getInterpVal(dx, dy, (*depTimeVec)[hB00], (*depTimeVec)[hB10], (*depTimeVec)[hB01], (*depTimeVec)[hB11]);    
}

//===========================
// Add Neighbors
//===========================

void addNeighbors(Heap *heap, vector< vector< graphNode* > > *Graph, graphNode *currNodePtr, double vx, double vy,double grid[][2], double pCost[], int nt, int dT,double xEps, int nSt, int nEnd, double xSpGT, double ySpGT, vector <double> *depTimeVector, int *nExcdDepTime){//, mutex *heapGuard, vector<mutex> *graphGuard){

    for (int a=nSt; a<nEnd ; a++){
       int nx = currNodePtr->x + ( vx + grid[a][0] )*dT;     // select x and y positions of the neighb
       int ny = currNodePtr->y + ( vy + grid[a][1] )*dT;
    
       if( !isAccessible(nx,ny,nt) )             // ignore this node if this is not accessible
           continue;
       
       //if ( nt > getLatestDepTime(nx, ny, depTimeVector, xSpGT, ySpGT) ){ // ignore if the time at node is more than the latest time of departure
       //    (*nExcdDepTime)++;
       //    continue; 
       //}
       
       //if( (!isAccessible(nx,ny,nt)) || (!isReachable(nx,ny,nt)) )             // ignore this node if this is not accessible
       //    continue;
    
       double cost = pCost[a]*dT;          // cost = precomputed cost (for unit dT) * dT
       int hashBin = getHashBinNumber(nx,ny,nt);
    
       /*find the node in the graph*/
       graphNode *neighbNode = NULL;
       
       //stSearch = clock();
       {
           //lock_guard<mutex> graphLock( (*graphGuard)[hashBin] );
           for(int i=0; i<(*Graph)[hashBin].size(); i++){
               if( (*Graph)[hashBin][i]->t == nt ){     // proceed if the nodes have the same time cordinate
                   double nDist = sqrt( (double)( (*Graph)[hashBin][i]->x - nx )*(double)( (*Graph)[hashBin][i]->x - nx ) + (double)( (*Graph)[hashBin][i]->y - ny )*(double)( (*Graph)[hashBin][i]->y - ny ) );
                   if ( nDist < xEps ){
                       neighbNode = (*Graph)[hashBin][i];
                       break;
                   }
               }
           }
       }
       //endSearch = clock();
       //searchTime = searchTime + (double)( endSearch - stSearch )/CLOCKS_PER_SEC;
    
       //if neighbor is not in the graph, add it
       if(!neighbNode){
           double heuristic = getHeuristic(nx,ny);
           neighbNode = new graphNode(nx, ny, nt,currNodePtr->g + cost,heuristic,currNodePtr,0);
           neighbNode->pathLength = neighbPathLength;
           {
                //lock_guard<mutex> graphLock( (*graphGuard)[hashBin] );
                (*Graph)[hashBin].push_back(neighbNode);   // add this neighbor to graph
           }
           {
               //lock_guard<mutex> heapLock(*heapGuard);
               heap->pushNode(neighbNode);              // add node to heap
           }
           continue;                               // go to next neighbor
       }
    
       // if Neighbor is not expanded already
       {
           //lock_guard<mutex> graphLock( (*graphGuard)[hashBin] );
           if(!neighbNode->expanded){
               if(neighbNode->g > (currNodePtr->g + cost)){    // if current cost to come to node is higher
                   neighbNode->g = (currNodePtr->g + cost);
                   neighbNode->f = neighbNode->g + neighbNode->h;
                   neighbNode->parent = currNodePtr;
                   neighbNode->pathLength = neighbPathLength;
                   {
                       //lock_guard<mutex> heapLock(*heapGuard);
                       heap->update(neighbNode->heapPos);
                   }
               }
           }
       }
    } // end of looking through neighbors in the grid for a given node currNode
}

//===========================
// Main Code
//===========================
int main(){

    /* Read data files containing velocity information
     * ------------------------------------------------*/
    vXVec = readDataToVecs("../U.txt", nX, nY, nT);
    vYVec = readDataToVecs("../V.txt", nX, nY, nT);

    /* Read in obstacle data
     * -------------------------*/
    OBS = readDataToVecs("../obstacles.txt", nX, nY, 1);
    
    /* Open files to save result
    ----------------------------*/
    ofstream outfProgress;                          // open file to save path parameters
    outfProgress.open("Progress.txt",ofstream::out | ofstream::trunc);

    ofstream outf, thisPathProgress;                // open file to save paths

    // reference path for comparison
    optPath refPath;
    refPath.readPathFromFile("../refPath2.txt");       // read the optimalpath into file
    refPath.createInterpolant();                    // create interpolants using the data

    /* obtain latest departure time for coordinates on the grid
    -------------------------------------------------------------*/
        clock_t depSearchSt = clock();
        
        cout << "Computing latest time of departure..." << endl;
        vector <double> depTimeVec;
        double xSpGridOptT = 200;
        double ySpGridOptT = 200;
        double tSpGridOptT = 200;
    
        depTimeVec = getDepTimeVect(xSpGridOptT, ySpGridOptT, tSpGridOptT);
        cout << "Computing latest departure times complete." << endl;
    
        double depSearchTime = ( clock() - depSearchSt )/CLOCKS_PER_SEC;
        cout << "time for latest departure time search " << depSearchTime << endl;
    /* general variable definitions
    ---------------------------------*/
    int nExpandedNodes;                 // number of Nodes expanded so far

    double cost;                        // cost of transition between nodes

    double vx;                          // flow velocities
    double vy;
    double vF;                          // flow velocity
    double thFlow;                      // angle of the flow
    double thHdg;                       // angle of heading
    double th;                          // heading wrt to flow
    double dvXdx, dvXdy, dvYdx, dvYdy, dvXdt, dvYdt;

    int nx, ny, nt;                     // neighbor spacial cordinates
    graphNode* neighbNode;              // temporary pointer to neighbor node
    double heuristic;                   // cost to go estimate to target
    double dxAllowed, dtAllowed;        // max space-time steps allowed at each expansion
    int dT;                             // selected dt
    double nDist, xEps;

    double vMin, vMax, vNom;            // min, max and nominal speed along a given direction
    double vi;                          // velocity for each node along a direction
    double delV;                        // velocity increment along a direction
    double delX, delT;                  // the t and x neighborhood of a node
    double vStill;                      // velocity selected along an edge
    double vSel;

    bool gFound = false;                // indicator if goal is found

    clock_t stSearch, endSearch;
    double searchTime;

#if dxType==1   // considers 3D hermV to compute dxAllowed, dtAllowed
    Matrix<double, 2, 3> divV;
    Matrix3d hermV, eigenVecs;
    Vector3d eigenVals, maxEigVec;
#elif dxType==2 // considers 2D hermV to compute dxAllowed only
    Matrix<double, 2, 2> divV;
    Matrix2d hermV, eigenVecs;
    Vector2d eigenVals, maxEigVec;
#endif
    double maxEig;
    int maxEigInd;


    double intX1, intX2, intY1, intY2, intT1, intT2;
    double xr, yr, tr;

    vector< vector< double > > pathVect;  // vector to store the computed path
    vector< double > tRow;                     // vector to store a single x y t row of the path

    double meandx, meanSubdx;
    double mindx = 50000.0;
    double maxdx = 0.0;
    double minXeps = 1000.0;
    double maxXeps = 0.0;
    int ndTExceeds;
    int dtAllLessdT = 0;
    int dTmin = 1000;
    int dTmax = 0;
    double dTmean = 0;
    double vFmean = 0; double vFmin = Vfm; double vFmax = 0;
    int nExcdDepTime = 0;
    ofstream tempOut;

    /* Run the graph search for each set of search parameters
     * -------------------------------------------------------*/
    for (int pthNum=0; pthNum<nTotRuns; pthNum++){

        // set the number of directions and velocity levels considered
        nDivs = nDivsAll[pthNum];

        // reset all counters
        auto stTime = chrono::high_resolution_clock::now();
        clock_t searchStTime = clock();
        nExpandedNodes = 0;
        searchTime = 0;

        // create new filestreams to output data
        outf.open("Path"+to_string(pthNum)+".txt", ofstream::out | ofstream::trunc );
        thisPathProgress.open("ProgressPath"+to_string(pthNum)+".txt", ofstream::out | ofstream::trunc );
        /* Neighbor Grid setup
        -----------------------*/
        const int N = ( nDivs + 1 )*( 3*nDivs + 2 ) - ( 2*nDivs + 1 );    // total number of nodes in the hexagonal lattice
        double  grid[N][2];              // velcoities of intermediate grid positions
        double pCost[N];                // cost to reach each intermediate point (pre computed)
        int gridPtr = -1;
        for (int i=0; i<=2*nDivs ; i++){
            for (int j=0; j<=( 2*nDivs - abs(i-nDivs) ); j++){
                gridPtr++;
                grid[gridPtr][0] =  ( j -  ( 2*nDivs - abs(i-nDivs) )/2.0 )*Vm/nDivs;
                grid[gridPtr][1] =  ( i-nDivs )*cos(M_PI/6 )*Vm/nDivs;
                pCost[gridPtr] = k1+k2*pow( sqrt( (double)( grid[gridPtr][0] )*grid[gridPtr][0] + (double)( grid[gridPtr][1] )*grid[gridPtr][1] ), alpha );
            }
        }
        //temporary variables for debugging
        tempOut.open("tempDataOut.txt", ofstream::out | ofstream::trunc);
        meandx = 0;
        dTmean = 0;
        meanSubdx = 0;
        ndTExceeds = 0;
        /* create graph container
        --------------------------*/
        vector< vector<graphNode*> > Graph(nHashBins); // create an empty graph with the specified number of hashbins
        //vector<graphNode*>* Graph = new vector<graphNode*>[nHashBins]; // create an empty graph with the specified number of hashbins
        //std::vector<graphNode*>* Graph = new std::vector<graphNode*>[nHashBins]; // create an empty graph with the specified number of hashbins

        for(int i=0; i<nHashBins; i++)
        {
            Graph[i].clear();
            Graph[i].reserve(hashBinSize);                // allocate enough capacity to each hashBin
        }

        /* create new heap container
        -------------------------------*/
        Heap heap;                                          // new empty heap

        /* Mutexes for guarding heap and graph
         *-------------------------------------*/
        mutex heapGuard;
        vector<mutex> graphGuard(nHashBins);

        /*Add seed node to graph and heap
        ----------------------------------*/
        graphNode* currNodePtr = new graphNode(START_COORD[0],START_COORD[1],START_COORD[2],0,0,NULL,startt);   // create new node pointer with seed node
        int hashBin = getHashBinNumber(currNodePtr->x, currNodePtr->y, currNodePtr->t);                         // find the hash bin of the Graph that the node should go to
        Graph[hashBin].push_back(currNodePtr);                                                                  // add the Node to the correct hashBin in the graph
        //cout << "number of nodes in hashBin " << hashBin << " : " <<  Graph[hashBin].size() << endl;

        heap.pushNode(currNodePtr);     // add pointer to the node in the heap

        /*main graph search loop
        -------------------------*/
        while (!heap.isHeapEmpty()){       // repeat while the heap is not empty
            currNodePtr = heap.popNode();               // pop the node at the root of heap

            // check if current node is the goal location
            nDist = sqrt( (double)( currNodePtr->x-GOAL_COORD[0] )*( currNodePtr->x-GOAL_COORD[0] ) + (double)( currNodePtr->y - GOAL_COORD[1] )*( currNodePtr->y - GOAL_COORD[1] ) );
            if( nDist <= xEndEps ){     // stop if the current node is the end node
                gFound = true;
                break;
            }

            if ( currNodePtr->t > getLatestDepTime(currNodePtr->x, currNodePtr->y, &depTimeVec, xSpGridOptT, ySpGridOptT) ){ // ignore if the time at node is more than the latest time of departure
                nExcdDepTime++;
                continue; 
            }

            if ( currNodePtr->pathLength > maxPathLength )   // ignore node if the path length is more than maximum allowable length
                continue;

            currNodePtr->expanded = true;               // node is expanded
            nExpandedNodes++;
            
            // if not goal, proceed
            getFlowFromVecs(currNodePtr->x,currNodePtr->y,currNodePtr->t,vx,vy, dvXdx, dvXdy, dvYdx, dvYdy, dvXdt, dvYdt);   // get the flow at current node
            thFlow = atan2(vy,vx);                                                  // direction of the flow
            vF = sqrt(vx*vx + vy*vy);
            vSel =(vF<Vm/5.0)? Vm/5.0: vF;//(vF<Vm/1.0)? Vm/1.0: vF;
            
            vFmean = (vFmean*(nExpandedNodes-1) + vF)/nExpandedNodes;
            vFmax = (vF>vFmax)?vF:vFmax;
            vFmin = (vF<vFmin)?vF:vFmin;

            // compute maximum allowable dx and dt
#if dxType==1
            divV << dvXdx, dvXdy, dvXdt, dvYdx, dvYdy, dvYdt;         // divergence of flow
            hermV = divV.transpose()*divV;              // to compute the spectral norm of the divergence
            EigenSolver<Matrix3d> eigSolve(hermV);      // Eigenvalue solver object for the hermV matrix
#elif dxType==2
            divV << dvXdx, dvXdy, dvYdx, dvYdy;         // divergence of flow
            hermV = divV.transpose()*divV;              // to compute the spectral norm of the divergence
            EigenSolver<Matrix2d> eigSolve(hermV);      // Eigenvalue solver object for the hermV matrix
#endif            //divV << dvXdx, dvXdy, dvYdx, dvYdy;         // divergence of flow
            
            eigenVals = eigSolve.eigenvalues().real();  // eigenvalues of the hermV matrix, get only real parts since hermV is symmetric
            eigenVecs = eigSolve.eigenvectors().real(); // get eigenvectors of the matrix

            maxEig = eigenVals.maxCoeff(&maxEigInd);             // get maximum eigenvalue
            maxEigVec = eigenVecs.col(maxEigInd);      // get eigenvector corresponding to max eigenvalue

            dxAllowed = p*vSel/sqrt(maxEig);          // max allowable dx
            intX1 = currNodePtr->x + dxAllowed*maxEigVec(0);
            intY1 = currNodePtr->y + dxAllowed*maxEigVec(1);

#if dxType==1
            intT1 = currNodePtr->t + dxAllowed*maxEigVec(2);
            findDx(intX1,intY1,intT1,maxEigVec(0),maxEigVec(1),maxEigVec(2),vF,vSel,xr,yr,tr);
#elif dxType==2
            dtAllowed = p*vSel/( sqrt( dvXdt*dvXdt + dvYdt*dvYdt ) );   // max allowable dt
            intT1 = currNodePtr->t + dtAllowed;
            findDx(intX1, intY1, currNodePtr->t, maxEigVec(0), maxEigVec(1), vF, vSel, xr, yr);
            findDt(currNodePtr->x, currNodePtr->y, intT1, vF, vSel, tr);
#endif

            dxAllowed = sqrt( (xr-currNodePtr->x)*(xr-currNodePtr->x) + (yr-currNodePtr->y)*(yr-currNodePtr->y) );

            dtAllowed = fabs(tr-currNodePtr->t);

//            if (dxAllowed < dxMax/10.0) dxAllowed = dxMax/10;
//            if (dtAllowed < dtMax/10.0) dtAllowed = dtMax/10;
            
            // select the dt such that, so that all reachable points are withing dxAllowed
            if ( ( (Vm+vF)*dtAllowed ) > dxAllowed ){
                dtAllowed = dxAllowed/(vF+Vm);
                ndTExceeds++;
            }
            
            // round dt to the nearest time layer
            dT = round( (dtAllowed)/dTLayer )*dTLayer;
            if (dT==0)
                dT = dTLayer;
            
            // distance between two nodes in this expansion
            xEps = max( (Vm*dT)/(2*nDivs), xEpsMin);       

            // mean dx and dt limit computation for debugging
            meandx = ( meandx*(nExpandedNodes-1) + dxAllowed )/nExpandedNodes;
            mindx = (mindx>dxAllowed)?dxAllowed:mindx;
            maxdx = (maxdx<dxAllowed)?dxAllowed:maxdx;

            meanSubdx = ( meanSubdx*(nExpandedNodes-1) + xEps )/(nExpandedNodes); 
            maxXeps = (maxXeps<xEps)?xEps:maxXeps;
            minXeps = (minXeps>xEps)?xEps:minXeps;
            //if (dT>10000 && (vF<Vm/10) )
            //    cout<< "highDT"<<endl;
            
            dTmin = (dT < dTmin )?dT : dTmin; 
            dTmax = (dT > dTmax )?dT : dTmax; 
            dTmean = ( dTmean*(nExpandedNodes-1) + dT )/nExpandedNodes;
            // compute temporal coordinate, its the same for all neighbors
            nt = currNodePtr->t + dT;
            neighbPathLength = currNodePtr->pathLength + 1;
            /* add parallelization here
             *-------------------------- */
            addNeighbors(&heap, &Graph, currNodePtr, vx, vy, grid, pCost, nt, dT, xEps, 0, N, xSpGridOptT, ySpGridOptT, &depTimeVec, &nExcdDepTime);//, &heapGuard, &graphGuard);

            if(!(nExpandedNodes%10000)){
                cout << "Number of Nodes expanded : " << nExpandedNodes << ". Number of Nodes in the Heap : " << heap.heapSize << endl;
            }
            if(!(nExpandedNodes%100000)){
                tempOut << "***********************************************************************" << endl;
                tempOut << " nExp           = " << nExpandedNodes << endl;
                tempOut << " dtAll > dxAll/V: " << ndTExceeds << endl;
                tempOut << " dtAll < dT     : " << dtAllLessdT << endl<<endl;

                tempOut << " mindx          = " << mindx << endl;
                tempOut << " maxdx          = " << maxdx << endl;
                tempOut << " mean dxAllowed = " << meandx << endl<<endl;
                ;
                tempOut << " min xEps       = " << minXeps << endl;
                tempOut << " max xEps       = " << maxXeps << endl;
                tempOut << " mean xEps      = " << meanSubdx << endl<<endl;

                tempOut << " dtMin adptv    : " << dTmin << endl;
                tempOut << " dtMax adptv    : " << dTmax << endl;
                tempOut << " dtMean adptv   : " << dTmean << endl <<  endl;
                
                tempOut << " vFmean   : " << vFmean << endl;
                tempOut << " vFmax    : " << vFmax << endl;
                tempOut << " vFmin    : " << vFmin << endl;
                tempOut << "***********************************************************************" << endl;
            }
        } // end of while loop searching through the Heap
        /* End of main graph search loop */

        outfProgress << "dTlayer       : " << dTLayer << endl;
        outfProgress << "xEpsMin       : " << xEpsMin << endl;
        outfProgress << "nDivs         : " << nDivs << endl;
        outfProgress << "errThresh (p) : " << p << endl;
        outfProgress << "nExpanded     : " << nExpandedNodes << endl;

        thisPathProgress << "dTlayer       : " << dTLayer << endl;
        thisPathProgress << "xEpsMin       : " << xEpsMin << endl;
        thisPathProgress << "nDivs         : " << nDivs << endl;
        thisPathProgress << "errThresh (p) : " << p << endl;
        thisPathProgress << "nExpanded     : " << nExpandedNodes << endl;


        if(!heap.isHeapEmpty()){              // if the target is found before the end of the heap
            /* Save path details
            --------------------*/
            vector< graphNode* > path = getPlannedPath(currNodePtr);
            pathVect.clear();
            tRow.clear();

            // print path to file and compute meanError
            cout.precision(8);
            outf.precision(8);
            for (int a=0; a<path.size(); a++)
            {
                outf<<path[a]->x << " " << path[a]->y << " " << path[a]->t << " " << path[a]->g <<endl;
                tRow.push_back( path[a]->x );
                tRow.push_back( path[a]->y );
                tRow.push_back( path[a]->t );
                pathVect.push_back(tRow);
                tRow.clear();
            }

            double mE = refPath.getMeanPathError( pathVect );
            
            cout << "nDivs        : " << nDivs << endl;
            cout << "mE (m)       : " << mE << endl;
            cout << "Path cost    : " << path[0]->g << endl;
            cout << "TIme at Goal : " << path[0]->t << endl;
            cout << "nExpanded    : " << nExpandedNodes << endl;

            outfProgress << "mE (m)        : " << mE << endl;
            outfProgress << "Path cost     : " << path[0]->g << endl;
            outfProgress << "Time at Goal  : " << path[0]->t  << endl;

            thisPathProgress << "mE (m)        : " << mE << endl;
            thisPathProgress << "Path cost     : " << path[0]->g << endl;
            thisPathProgress << "Time at Goal  : " << path[0]->t  << endl;

            }
        else{
            cout << "NO PATH TO TARGET NODE FOUND!!!!!!" << endl;
            cout << "nExpanded    : " << nExpandedNodes << endl;
            outfProgress << "NO PATH TO TARGET NODE FOUND!!!!!!" << endl;
            thisPathProgress << "NO PATH TO TARGET NODE FOUND!!!!!!" << endl;
        }
        clock_t searchEndTime = clock();
        auto endTime = chrono::high_resolution_clock::now();
        //double runTime = chrono::duration_cast< chrono::duration<double> > (endTime-stTime).count();
        double runTime = (searchEndTime - searchStTime )/CLOCKS_PER_SEC;
        cout << "Search running Time : " << runTime << "s" << endl;
        cout << "Total run time with latest departure search " << runTime+depSearchTime << endl;
        outfProgress << "Running Time : " << runTime << "s" << endl;
        outfProgress << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        outfProgress << endl;

        thisPathProgress << "Running Time : " << runTime << "s"  << endl;

        /* delete Graph and Nodes it points to
        -------------------------------------*/
        int meanHashSize = 0;
        int maxHashSize = 0;
        int minHashSize = Graph[0].size();
        for(int i=0; i<nHashBins; i++){
            meanHashSize = meanHashSize + Graph[i].size();
            maxHashSize = (Graph[i].size()>maxHashSize)?Graph[i].size():maxHashSize;
            minHashSize = (Graph[i].size()<minHashSize)?Graph[i].size():minHashSize;
            for(int j=0; j<Graph[i].size(); j++)
                delete Graph[i][j];                // delete the Node pointed to by the graph
        }
        //delete[] Graph;
        cout << "max hash bin size = " << maxHashSize << endl;
        cout << "min hash bin size = " << minHashSize << endl;
        cout << "mean hash bin size = " << (double)meanHashSize/nHashBins << endl;
        cout << endl;
        cout << "mean search time for a node in Graph : " << searchTime/nExpandedNodes << endl;
        cout << endl;
        cout << "number of neighbors that exceeded the latest departure time " << nExcdDepTime << endl;
        cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        outf.close();
        thisPathProgress.close();
        tempOut.close();
    } // end for each pthNumi
    outfProgress.close();
    return 0;
}
