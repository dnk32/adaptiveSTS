#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <Eigen/Dense>  //including Eigen library for matrix manipulation
#include <time.h>
#include <gsl/gsl_spline.h>

#include "classDefs.h"
using namespace std;
using namespace Eigen;

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

std::vector< std::vector< std::vector< double > > > vXVec;
std::vector< std::vector< std::vector< double > > > vYVec;

//================================
// obstacle definitions
//================================
std::vector< std::vector< std::vector< double > > > OBS;
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
int tmax = 180000;
int tmin = 0;

//===========================
// Algorithm parameters
//===========================

// maximum flow speed encountered
double Vfm = 0.7288;  // m/s or mm/ms

//max vehicle speed allowed
double Vm = 0.5; // m/s or mm/ms

// maximum discretization levels
double dxMax = 10000;
double dtMax = 2000;

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
// discrete time layers considered
int dTLayer = 200;              // all time dimensions will be  multiples of dTlayer.
                                // this is done to speed up the search

// the dimensions of each hashbin
int dxBin = 250;
int dyBin = 250;
int dtBin = 500;

// number of hashBins along each axis
int nXBins = 100;
int nYBins = 100;
int nTBins = 20;

// size of the complete hashBin
int xSpBin = nXBins*dxBin;
int ySpBin = nYBins*dyBin;
int tSpBin = nTBins*dtBin;

// number of total hashBins
int nHashBins = nXBins*nYBins*nTBins;

// number of nodes in each hashbin
int hashBinSize = 200;

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


//===========================
// Main Code
//===========================
int main(){
    /* Read data files containing velocity information
     * ------------------------------------------------*/
    vXVec = readDataToVecs("./U.txt", nX, nY, nT);
    vYVec = readDataToVecs("./V.txt", nX, nY, nT);

    /* Read in obstacle data
     * -------------------------*/
    OBS = readDataToVecs("./obstacles.txt", nX, nY, 1);

    /* Open files to save result
    ----------------------------*/
    ofstream outfProgress;                          // open file to save path parameters
    outfProgress.open("Progress.txt",ofstream::out | ofstream::trunc);

    ofstream outf, thisPathProgress;                // open file to save paths

    // reference path for comparison
    optPath refPath;
    refPath.readPathFromFile("./refPath2.txt");       // read the optimalpath into file
    refPath.createInterpolant();                    // create interpolants using the data

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

    clock_t stTime, endTime, stSearch, endSearch;
    double searchTime;

    Matrix2d divV, hermV, eigenVecs;
    Vector2d eigenVals, maxEigVec;
    double maxEig;

    std::vector< std::vector< double > > pathVect;  // vector to store the computed path
    std::vector< double > tRow;                     // vector to store a single x y t row of the path

    double meandx, meandt, meanSubdx;
    int ndTExceeds;
    int dtAllLessdT = 0;
    int dTmin = dtMax;
    int dTmax = 0;
    double dTmean = 0;
    double vFmean = 0; double vFmin = Vfm; double vFmax = 0;
    ofstream tempOut;
    /* Run the graph search for each set of search parameters
     * -------------------------------------------------------*/
    for (int pthNum=0; pthNum<nTotRuns; pthNum++){

        // set the number of directions and velocity levels considered
        nDivs = nDivsAll[pthNum];

        // reset all counters
        stTime = clock();
        nExpandedNodes = 0;
        searchTime = 0;

        // create new filestreams to output data
        outf.open("Path"+to_string(pthNum)+".txt", ofstream::out | ofstream::trunc );
        thisPathProgress.open("ProgressPath"+to_string(pthNum)+".txt", ofstream::out | ofstream::trunc );
        /* Neighbor Grid setup
        -----------------------*/
        const int N = ( nDivs + 1 )*( 3*nDivs + 2 ) - ( 2*nDivs + 1 );    // total number of nodes in the hexagonal lattice
        double grid[N][2];              // velcoities of intermediate grid positions
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
        meandt = 0;
        meanSubdx = 0;
        ndTExceeds = 0;
        /* create graph container
        --------------------------*/
        std::vector<graphNode*>* Graph = new std::vector<graphNode*>[nHashBins]; // create an empty graph with the specified number of hashbins
        //std::vector<graphNode*>* Graph = new std::vector<graphNode*>[nHashBins]; // create an empty graph with the specified number of hashbins

        for(int i=0; i<nHashBins; i++)
        {
            Graph[i].reserve(hashBinSize);                // allocate enough capacity to each hashBin
        }

        /* create new heap container
        -------------------------------*/
        Heap heap;                                          // new empty heap

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
            currNodePtr->expanded = true;               // node is expanded
            nExpandedNodes++;

            // check if current node is the goal location
            nDist = sqrt( (double)( currNodePtr->x-GOAL_COORD[0] )*( currNodePtr->x-GOAL_COORD[0] ) + (double)( currNodePtr->y - GOAL_COORD[1] )*( currNodePtr->y - GOAL_COORD[1] ) );
            if( nDist <= xEndEps ){     // stop if the current node is the end node
                gFound = true;
                break;
            }

            // if not goal, proceed
            getFlowFromVecs(currNodePtr->x,currNodePtr->y,currNodePtr->t,vx,vy, dvXdx, dvXdy, dvYdx, dvYdy, dvXdt, dvYdt);   // get the flow at current node
            thFlow = atan2(vy,vx);                                                  // direction of the flow
            vF = sqrt(vx*vx + vy*vy);
            
            vFmean = (vFmean*(nExpandedNodes-1) + vF)/nExpandedNodes;
            vFmax = (vF>vFmax)?vF:vFmax;
            vFmin = (vF<vFmin)?vF:vFmin;

            // compute maximum allowable dx and dt
            divV << dvXdx, dvXdy, dvYdx, dvYdy;         // divergence of flow
            hermV = divV.transpose()*divV;              // to compute the spectral norm of the divergence

            EigenSolver<Matrix2d> eigSolve(hermV);      // Eigenvalue solver object for the hermV matrix
            Matrix<complex<double>, 2,1> eigVs = eigSolve.eigenvalues();    // eigenvalues of the hermV matrix

            eigenVals = eigVs.real();                   // get only the real parts, since hermV is symmetric
            maxEig = eigenVals(0)>=eigenVals(1)?eigenVals(0):eigenVals(1);   // get maximum eigen value and eigen vector
            
            vSel =(vF<Vm/1.0)? Vm/1.0: vF;
            dxAllowed = min( p*vSel/sqrt(maxEig), dxMax );          // max allowable dx
            dtAllowed = min( p*vSel/( sqrt( dvXdt*dvXdt + dvYdt*dvYdt ) ), dtMax );   // max allowable dt
//            if (dxAllowed < dxMax/10.0) dxAllowed = dxMax/10;
//            if (dtAllowed < dtMax/10.0) dtAllowed = dtMax/10;
            
            // select the dt such that, so that all reachable points are withing dxAllowed
            if ( ( (Vm+vF)*dtAllowed ) > dxAllowed ){
                dtAllowed = dxAllowed/(vF+Vm);
                ndTExceeds++;
            }
            
            // round dt to the nearest time layer
            dT = round( (dtAllowed)/dTLayer )*dTLayer;
            //dT = 1000;
            //if (dtAllowed < dT)
            //    dtAllLessdT++;

            //if (dT<9*dTLayer)
            //    dT = 9*dTLayer;
            // distance between two nodes in this expansion
            xEps = (Vm*dT)/nDivs;       

            // mean dx and dt limit computation for debugging
            meandx = ( meandx*(nExpandedNodes-1) + dxAllowed )/nExpandedNodes;
            meandt = ( meandt*(nExpandedNodes-1) + dtAllowed )/nExpandedNodes;
            meanSubdx = ( meanSubdx*(nExpandedNodes-1) + xEps )/(nExpandedNodes); 

            dTmin = (dT < dTmin )?dT : dTmin; 
            dTmax = (dT > dTmax )?dT : dTmax; 
            dTmean = ( dTmean*(nExpandedNodes-1) + dT )/nExpandedNodes;
            // compute temporal coordinate, its the same for all neighbors
            nt = currNodePtr->t + dT;
            for (int a=0; a<N ; a++){
               nx = currNodePtr->x + ( vx + grid[a][0] )*dT;     // select x and y positions of the neighb
               ny = currNodePtr->y + ( vy + grid[a][1] )*dT;

               if(!isAccessible(nx,ny,nt))             // ignore this node if this is not accessible
                   continue;

               cost = pCost[a]*dT;          // cost = precomputed cost (for unit dT) * dT
               hashBin = getHashBinNumber(nx,ny,nt);

               /*find the node in the graph*/
               neighbNode = NULL;
               stSearch = clock();
               for(int i=0; i<Graph[hashBin].size(); i++){
                   if(Graph[hashBin][i]->t == nt){     // proceed if the nodes have the same time cordinate
                       nDist = sqrt( (double)( Graph[hashBin][i]->x - nx )*(double)( Graph[hashBin][i]->x - nx ) + (double)( Graph[hashBin][i]->y - ny )*(double)( Graph[hashBin][i]->y - ny ) );
                       if ( nDist < xEps/2.0 ){
//                       if(Graph[hashBin][i]->x == nx && Graph[hashBin][i]->y == ny){
                           neighbNode = Graph[hashBin][i];
                           break;
                       }
                   }
               }
               endSearch = clock();
               searchTime = searchTime + (double)( endSearch - stSearch )/CLOCKS_PER_SEC;

               //if neighbor is not in the graph, add it
               if(!neighbNode){
                   heuristic = getHeuristic(nx,ny);
                   neighbNode = new graphNode(nx, ny, nt,currNodePtr->g + cost,heuristic,currNodePtr,0);
                   Graph[hashBin].push_back(neighbNode);   // add this neighbor to graph
                   heap.pushNode(neighbNode);              // add node to heap
                   continue;                               // go to next neighbor
               }

               // if Neighbor is not expanded already
               if(!neighbNode->expanded){
                   if(neighbNode->g > (currNodePtr->g + cost)){    // if current cost to come to node is higher
                       neighbNode->g = (currNodePtr->g + cost);
                       neighbNode->f = neighbNode->g + neighbNode->h;
                       neighbNode->parent = currNodePtr;
                       heap.update(neighbNode->heapPos);
                   }
               }
            } // end of looking through neighbors in the grid for a given node currNode


            if(!(nExpandedNodes%10000)){
                cout << "Number of Nodes expanded : " << nExpandedNodes << ". Number of Nodes in the Heap : " << heap.heapSize << endl;
            }
            if(!(nExpandedNodes%100000)){
                tempOut << "***********************************************************************" << endl;
                tempOut << " mean dxAllowed = " << meandx << endl;
                tempOut << " mean dtAllowed = " << meandt << endl;
                tempOut << " mean xEps      = " << meanSubdx << endl;
                tempOut << " nExp           = " << nExpandedNodes << endl<<endl;
                tempOut << " dtAll > dxAll*V: " << ndTExceeds << endl;
                tempOut << " dtAll < dT     : " << dtAllLessdT << endl;
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

        outfProgress << "dxMax         : " << dxMax << endl;
        outfProgress << "dtMax         : " << dtMax << endl;
        outfProgress << "nDivs         : " << nDivs << endl;
        outfProgress << "errThresh (p) : " << p << endl;
        outfProgress << "nExpanded     : " << nExpandedNodes << endl;

        thisPathProgress << "dxMax         : " << dxMax << endl;
        thisPathProgress << "dtMax         : " << dtMax << endl;
        thisPathProgress << "nDivs         : " << nDivs << endl;
        thisPathProgress << "errThresh (p) : " << p << endl;
        thisPathProgress << "nExpanded     : " << nExpandedNodes << endl;


        if(!heap.isHeapEmpty()){              // if the target is found before the end of the heap
            /* Save path details
            --------------------*/
            std::vector< graphNode* > path = getPlannedPath(currNodePtr);
            pathVect.clear();
            tRow.clear();

            // print path to file and compute meanError
            std::cout.precision(8);
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

        endTime = clock();
        cout << "Running Time : " << double(endTime - stTime)/CLOCKS_PER_SEC << endl;
        outfProgress << "Running Time : " << double(endTime - stTime)/CLOCKS_PER_SEC << endl;
        outfProgress << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        outfProgress << endl;

        thisPathProgress << "Running Time : " << double(endTime - stTime)/CLOCKS_PER_SEC << endl;

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
        delete[] Graph;
        cout << "max hash bin size = " << maxHashSize << endl;
        cout << "min hash bin size = " << minHashSize << endl;
        cout << "mean hash bin size = " << (double)meanHashSize/nHashBins << endl;
        cout << endl;
        cout << "mean search time for a node in Graph : " << searchTime/nExpandedNodes << endl;
        cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
        outf.close();
        thisPathProgress.close();
        tempOut.close();
    } // end for each pthNumi
    outfProgress.close();
    return 0;
}