#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <vector>
#include <time.h>
#include <string.h>
#include <sstream>      // for string stream
#include <algorithm>

//#include "classDefs.h"
#include "latestDepSearch.hpp"


using namespace std;
/*=====================
 * public functions
 *=====================*/

void latestDepSearch::setEnvLims(int maxX, int minX, int maxY, int minY, int maxT, int minT){
    xmax = maxX; xmin = minX; ymax = maxY; ymin = minY; tmax = maxT; tmin = minT;
}

void latestDepSearch::setGridParams(double xGrid, double yGrid, double tGrid, int nXgrid, int nYgrid){
    xSpGrid = xGrid; ySpGrid = yGrid; tSpGrid = tGrid;
    nX = nXgrid;
    nY = nYgrid;
    //hashBinHeight = hashBinCap*tGrid;
}

void latestDepSearch::setHashBinParams(int nHBins, int hBinHeight, int (*hashFnc)(int, int, double, bool, int, int, int) ){
    nHashBins = nHBins;
    getHashBinNumber = hashFnc;    
    hashBinHeight = hBinHeight;
}

void latestDepSearch::setDataVecs(double xData, double yData, double tData, int nXdata, int nYdata, int nTdata , dblVec3D *uVec, dblVec3D *vVec, dblVec3D *obsVec){

    xSpData = xData; ySpData = yData; tSpData = tData;
    nDx = nXdata; nDy = nYdata; nDt = nTdata;
    UxVec = uVec;
    VyVec = vVec;
    OBS = obsVec;
    dataParamsSet = true;
}

void latestDepSearch::setVelParams(double vMaxVeh, double vMaxFlow){
    Vm = vMaxVeh; Vfm = vMaxFlow;
}

void latestDepSearch::setStartCoords(double xSt, double ySt){
    startx = xSt;
    starty = ySt;
    startt = tmax;
}


//vector <double> latestDepSearch::getEarliestDepTimes();

/*=====================
 * private functions
 *=====================*/

vector< graphNode* > latestDepSearch::getPlannedPath(graphNode* currNodePtr){
    vector<graphNode*> path;
    path.clear();
    path.push_back(currNodePtr);
    graphNode* parentNodePtr = currNodePtr->parent;
    while(parentNodePtr)
    {
        path.push_back(parentNodePtr);
        parentNodePtr = parentNodePtr->parent;
    }
    return path;
}


//int latestDepSearch::getHashBinNumber(int x, int y, double t, bool isVmSmall){
//    if (isVmSmall)
//        return (int) ( (x)*nY + y + nX*nY*( (int)t/hashBinHeight ) );    // the x,y cords of any two nodes in the same hashBin will be the same. Only the t would change.
//    else
//        return (int) (x*nY + y);
//}

void latestDepSearch::getFlowVelFromVecs(int x,int y,double t,double &vx,double &vy){
    /* Inputs
        x : x cordinate in the grid
        y : y cordinate in the grid
        t : time at node
        vx : flow velocity in the x direction
        vy : flow velocity in the y direction
    */

    int nx = ( x*xSpGrid + xmin )/xSpData;
    int ny = ( y*ySpGrid + ymin )/ySpData;
    int nt = t/tSpData;

    double dx = ( x*xSpGrid + xmin ) - ( nx*xSpData );
    double dy = ( y*ySpGrid + ymin ) - ( ny*ySpData );
    double dt = t - ( nt*tSpData );

    /* using interpolation
    -----------------------*/
    if(nx == nDx || ny == nDy || nt == nDt){
        vx = (*UxVec)[nt][nx][ny];
        vy = (*VyVec)[nt][nx][ny];
        return;
    }

    double u1 = (double) dx/xSpData * ( (*UxVec)[nt][nx+1][ny] - (*UxVec)[nt][nx][ny] ) + (*UxVec)[nt][nx][ny];
    double u2 = (double) dx/xSpData * ( (*UxVec)[nt][nx+1][ny+1] - (*UxVec)[nt][nx][ny+1] ) + (*UxVec)[nt][nx][ny+1];
    double u3 = (double) dx/xSpData * ( (*UxVec)[nt+1][nx+1][ny] - (*UxVec)[nt+1][nx][ny]) + (*UxVec)[nt+1][nx][ny];
    double u4 = (double) dx/xSpData * ( (*UxVec)[nt+1][nx+1][ny+1] - (*UxVec)[nt+1][nx][ny+1]) + (*UxVec)[nt+1][nx][ny+1];

    double u12 = (double) dy/ySpData * ( u2 - u1 ) + u1;
    double u34 = (double) dy/ySpData * ( u4 - u3 ) + u3;

    vx = (double) dt/tSpData * ( u34 - u12 ) + u12;

    u1 = (double) dx/xSpData * ( (*VyVec)[nt][nx+1][ny] - (*VyVec)[nt][nx][ny] ) + (*VyVec)[nt][nx][ny];
    u2 = (double) dx/xSpData * ( (*VyVec)[nt][nx+1][ny+1] - (*VyVec)[nt][nx][ny+1] ) + (*VyVec)[nt][nx][ny+1];
    u3 = (double) dx/xSpData * ( (*VyVec)[nt+1][nx+1][ny] - (*VyVec)[nt+1][nx][ny]) + (*VyVec)[nt+1][nx][ny];
    u4 = (double) dx/xSpData * ( (*VyVec)[nt+1][nx+1][ny+1] - (*VyVec)[nt+1][nx][ny+1]) + (*VyVec)[nt+1][nx][ny+1];

    u12 = (double) dy/ySpData * ( u2 - u1 ) + u1;
    u34 = (double) dy/ySpData * ( u4 - u3 ) + u3;

    vy = (double) dt/tSpData * ( u34 - u12 ) + u12;

//    vx = (Ux[x][y][t]);
//    vy = (Vy[x][y][t]);
//    vx = (Ux[x][y][0]);
//    vy = (Vy[x][y][0]);
}

double latestDepSearch::getOptTimeCost(graphNode* currNode, int nx, int ny, double Vf, double thF){
    double delY = ( ny-currNode->y )*ySpGrid;
    double delX = ( nx-currNode->x )*xSpGrid;
    double thHdg = atan2( delY, delX);
    
    double thLoc = thHdg - thF;             // heading wrt flow direction
    if ( fabs(thLoc) > M_PI )                // consider the acute angle with flow
        thLoc = 2*M_PI - abs(thLoc);

    if ( Vm<fabs( Vf*sin(thLoc) ) )          // if Vm is not enough to reach this direction
        return -1.0;

    double vTot = Vf*cos(thLoc) + sqrt( Vm*Vm - Vf*sin(thLoc)*Vf*sin(thLoc) );

    return sqrt(delX*delX+delY*delY)/(vTot);// returned value will be negative if Vf>Vm for thLoc>M_PI/2
}

bool latestDepSearch::isAccessible(int nx,int ny,double t){
    int nxorg = ( nx*xSpGrid + xmin )/xSpData;
    int nyorg = ( ny*ySpGrid + ymin )/ySpData;

    bool insideObs = (bool)( (*OBS)[0][nxorg][nyorg] == 1);
    if ( !insideObs && nx>=0 && nx<nX && ny>=0 && ny<nY && t>=tmin && t<=tmax)
		return true;
	else
	{
		return false;
	}

}

void latestDepSearch::runSearch( vector <double> *depTimeVec){
    (*depTimeVec).clear();
    (*depTimeVec).resize(nX*nY);
    fill( (*depTimeVec).begin(), (*depTimeVec).end(), -1.0 );

    bool isVmSmall = (Vm<=Vfm);

    // setting start and end coordinate
    START_COORD[0] = (int)( ( startx - xmin )/xSpGrid ); START_COORD[1] = (int)( ( starty - ymin )/ySpGrid ); START_COORD[2] =  (int)( ( startt - tmin )/tSpGrid );

    clock_t stTime = clock();
    /* general variable definitions
    ---------------------------------*/
    int nExpandedNodes = 0;             // number of Nodes expanded so far

    double cost;                        // cost of transition between nodes

    double vx;                          // flow velocities
    double vy;
    double U;                           // flow speed
    double thFlow;                      // direction of flow

    int nx, ny;                         // neighbor spacial cordinates
    double tAtNeighb;                   // time at neighbor node
    graphNode* neighbNode;              // temporary pointer to neighbor node
    double heuristic;                   // cost to go estimate to target

    clock_t stSearch, endSearch;
    double searchTime = 0;

    int totNodesExplored = 0;           // total number of spatial locations explored
    /* create graph container
    --------------------------*/
    vector< vector<graphNode*> > Graph(nHashBins); // create an empty graph with the specified number of hashbins
    for(int i=0; i<nHashBins; i++)
    {
        Graph[i].clear();
        Graph[i].reserve(hashBinCap);                // allocate enough capacity to each hashBin.
    }

    /* create new heap container
    -------------------------------*/
    Heap heap;                                          // new empty heap

    ///* coverage map
    //-------------------------------*/
    //vector<graphNode*> covMap(nX*nY,NULL);          // coverage map indicating earliest arriving node
    
    /*Add seed node to graph and heap
    ----------------------------------*/
    graphNode* currNodePtr = new graphNode(START_COORD[0],START_COORD[1],START_COORD[2],0,0,NULL,startt);   // create new node pointer with seed node
    int hashBin = getHashBinNumber(currNodePtr->x, currNodePtr->y, currNodePtr->timeAtNode,isVmSmall, nX, nY, hashBinHeight);                         // find the hash bin of the Graph that the node should go to
    Graph[hashBin].push_back(currNodePtr);                                                                  // add the Node to the correct hashBin in the graph

    heap.pushNode(currNodePtr);     // add pointer to the node in the heap

    /*main graph search loop
    -------------------------*/
    while (!heap.isHeapEmpty()){       // repeat while the heap is not empty
        currNodePtr = heap.popNode();                                           // pop the node at the root of heap
        int covMapBin = getHashBinNumber(currNodePtr->x,currNodePtr->y,0,false, nX, nY, hashBinHeight);
        if ( (*depTimeVec)[covMapBin] < 0 ){
            (*depTimeVec)[covMapBin] = currNodePtr->timeAtNode;
            totNodesExplored++;
            if ( !(totNodesExplored % ( (nX*nY)/20 ) ) )
                cout << "Approximately " << ( (double)totNodesExplored )/(nX*nY)*100 << "% complete. " << endl; 
        }

        //if(currNodePtr->x == GOAL_COORD[0] && currNodePtr->y == GOAL_COORD[1])     // stop if the current node is the end node
        //    break;

        if ( currNodePtr->timeAtNode<=tmin || totNodesExplored==(nX*nY) )
            break;

        //getFlowVel(currNodePtr->x,currNodePtr->y,currNodePtr->t,vx,vy);         // get the flow at current node
        getFlowVelFromVecs(currNodePtr->x,currNodePtr->y,currNodePtr->timeAtNode,vx,vy);         // get the flow at current node
        U = sqrt(vx*vx + vy*vy);
        thFlow = atan2(-vy,-vx);                  // flow direction

        /*Find valid neighbors of current node*/
        for (int a=-lim; a<=lim; a++){
            for (int b=-lim; b<=lim; b++){
                if (a==0 && b==0) continue;
                if (lim!=1 && ((abs(a)==abs(b) && abs(b)!=1)||(a==0 && abs(b)!=1)||(b==0 && abs(a)!=1))) continue;      // ignore repeated neighbors in the same direction
                nx = currNodePtr->x + a;
                ny = currNodePtr->y + b;
                cost = getOptTimeCost(currNodePtr, nx, ny, U, thFlow);
                if ( cost<=0 )    // ignore if direction is not accessible
                    continue;
                
                tAtNeighb = currNodePtr->timeAtNode - cost;
                if (!isAccessible(nx, ny, tAtNeighb)) // ignore if the node is obstructed
                        continue;
                
                /*Find the accessible neighbor in the Graph, if not add to graph*/
                hashBin = getHashBinNumber(nx,ny,tAtNeighb, isVmSmall, nX, nY, hashBinHeight);       // if the neighbor is accessible, get its hashbin in the Graph

                //finding if node is in the graph
                stSearch = clock();
                
                neighbNode = NULL;
                if( !isVmSmall ){
                    if( Graph[hashBin].size()!=0 )
                        neighbNode = Graph[hashBin][0];
                }
                else{
                    for(int i=0; i<Graph[hashBin].size(); i++){
                        if ( fabs( Graph[hashBin][i]->timeAtNode - tAtNeighb ) < tSpGrid ){        // only t is considered because in every node in a hashbin, only the t cordinate will differ
                            neighbNode = Graph[hashBin][i];     // if the cordinates match
                            break;                              // no need to check further
                        }
                    }
                }
                endSearch = clock();
                searchTime = searchTime + (double)( endSearch - stSearch )/CLOCKS_PER_SEC;

                //if neighbor is not in the graph, add it
                if(!neighbNode){
                    heuristic = 0;//getHeuristic(nx,ny);
                    neighbNode = new graphNode(nx, ny, (tAtNeighb-tmin)/tSpGrid ,currNodePtr->g + cost,heuristic,currNodePtr,tAtNeighb);
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
                        neighbNode->timeAtNode = tAtNeighb;
                        neighbNode->t = (tAtNeighb-tmin)/tSpGrid;
                        heap.update(neighbNode->heapPos);
                    }
                }

            }// end for b
        }// end for a
        currNodePtr->expanded = true;
        nExpandedNodes++;
        //if(!(nExpandedNodes%10000))
        //        cout << "Number of Nodes expanded : " << nExpandedNodes << ". Number of Nodes in the Heap : " << heap.heapSize << endl;
    }
    
    clock_t endTime = clock();
    
    ///* Getting user input to compute a path
    //----------------------------------------*/
    //char userInput;
    //int pathNumber = 0;
    //double xIn, yIn;
    //int xG, yG;
    //while (true){
    //    //get user input
    //    cout << "Ready to save paths. continue? (y)yes, (n)no" << endl;
    //    cin >> userInput;
    //    while (userInput != 'y' && userInput != 'n'){
    //        cout << "Input not recognize. (y) for more paths, (n) to end" << endl;
    //        cin >> userInput;
    //    }
    //    if (userInput=='n')
    //        break;
    //    cout << "enter x coordinate of desired goal " << endl;
    //    cin >> xIn;
    //    cout << "enter y coordinate of desired goal " << endl;
    //    cin >> yIn;
    //    
    //    xG = (xIn-xmin)/xSpGrid;
    //    yG = (yIn-ymin)/ySpGrid;

    //    if ( !isAccessible(xG,yG,(tmin+tmax)/2) ){
    //        cout << "coordinates are out of bounds!!!" << endl;
    //        continue;
    //    }
    //    pathNumber++;
    //    int covMapBin = getHashBinNumber(xG,yG,0,false);
    //    currNodePtr = covMap[covMapBin];
    //   
    //    if(!currNodePtr){
    //        cout << "no path to goal found" << endl;
    //        continue;
    //    }
    //    /* Save path details
    //    --------------------*/
    //    ofstream outfPath;                          // open file to save paths
    //    outfPath.open("Path"+to_string(pathNumber)+".txt",ofstream::out | ofstream::trunc);
    //    
    //    vector< graphNode* > path = getPlannedPath(currNodePtr);

    //    cout.precision(8);
    //    outfPath.precision(8);
    //    for (int a=0; a<path.size(); a++)
    //    {
    //        outfPath << ( path[a]->x*xSpGrid + xmin ) << " " << ( path[a]->y*ySpGrid + ymin ) << " " << ( path[a]->timeAtNode ) << " " << path[a]->g <<endl;
    //    }

    //    outfPath.close();
    //    
    //    
    //}

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
    cout << "max hash bin size = " << maxHashSize << endl;
    cout << "min hash bin size = " << minHashSize << endl;
    cout << "mean hash bin size = " << (double)meanHashSize/nHashBins << endl;
    cout << endl;
    cout << "mean search time for a node in Graph : " << searchTime/nExpandedNodes << endl;
    cout << "time taken for graph search " << (double)( clock() - stTime )/CLOCKS_PER_SEC << endl; 
}

