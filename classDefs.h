#ifndef CLASSDEFS_H_INCLUDED
#define CLASSDEFS_H_INCLUDED

#include <gsl/gsl_spline.h>

//=============================
// Search direction enumurator
//=============================
enum searchType{
    FREE_END_TIME = 1,          // list of fixed start times provided
    FREE_START_TIME = -1        // list of fixed end times provided
};


//===================
// Class graphNode
//===================

class graphNode{
    public:
    int x;          // node position
    int y;          // node y position
    int t;          // node time
    double g;       // cost to come to node
    double h;       // heuristic to target
    double f;       // cost to go + cost to come, used for sorting heap
    double timeAtNode; //time at node
    int heapPos;    // position of the node in the heap
    bool expanded;  // indicate if the node has been expanded
    int pathLength; // length of path leading up to the node

    graphNode* parent;  //pointer to parent node

    graphNode(int xx, int yy, int tt, double gg, double hh, graphNode* parentNode, double time)       // constructor for a new node
    {
        x = xx; y = yy; t = tt;
        g = gg;
        h = hh;
        f = g+h;
        parent = parentNode;
        expanded = false;
        timeAtNode = time;
        pathLength = 0;
    }
};

//class tempNode{
//    public:
//    int x;          // node position
//    int y;          // node y position
//    int t;          // node time
//};

//=================================
// Class Heap
//=================================

class Heap{
    public:
    std::vector<graphNode*> heapArray;                  //contains the heap
    enum BubbleDirection{UP, DOWN, BOTH} ;              // structure that contains the heap sorting direction, percolating up or percolating down
    int heapSize;                                       // size of the heap

    Heap(){                                             //clear HeapArray on Start
        heapArray.clear();
        heapSize=0;
    }

    void pushNode(graphNode* ptrN)                      // push a graph node pointer to the heap and sort heap so that the minimum cost node is at the root
    {
        heapArray.push_back(ptrN);                      // add node to heap
        heapSize++;
        int newPos = heapArray.size()-1;                // position of new value
        heapArray[newPos]->heapPos = newPos;            // update the heapPos parameter of the Node with the position it is in the heap
        update(newPos,UP);                              // sort heap so that min cost node is at the root
    }

    graphNode* popNode()                                // return the graphNode with the least cost and restore graph to give next least cost node at the root
    {
        graphNode* ret = heapArray[0];                  // return the minimum cost node
        heapArray[0] = heapArray[heapArray.size()-1];   // put last node in the first node position
        heapArray[0]->heapPos = 0;
        heapArray.pop_back();                           // remove last node in the vector
        update(0,DOWN);
        heapSize--;
        return ret;
    }

    void update(int currPos,BubbleDirection dir=BOTH)        // sort the heap. dir is UP for push operations and DOWN for pop operations, BOTH for change of value operations
    {
        if(dir!=DOWN){                                            // a update from a push operation or a change opration
            if(currPos==0 && dir==UP)                             // currently at the root, sorting complete
                return;
            int parentPos = (currPos-1)/2;
            if(heapArray[parentPos]->f > heapArray[currPos]->f){    // if parent cost is higher, swap nodes
                graphNode* tempNode = heapArray[parentPos];
                heapArray[parentPos] = heapArray[currPos];
                heapArray[currPos] = tempNode;
                heapArray[parentPos]->heapPos = parentPos;      // update heapPos parameter of the Node with its current position in the heap
                heapArray[currPos]->heapPos = currPos;
                update(parentPos,UP);                           // recursively go up the heap until the parent has lower cost or until the root is reached
                return;                                         // if parent had smaller value, even in dir == BOTH, no need to check downward direction. Thus return
            }
            if (dir == UP) return;                              // if parent has lower cost, return, but if dir
        }

        // if the sort direction is DOWN (from a pop operation)
        int leftChildPos = 2*currPos + 1;
        int rightChildPos = 2*currPos + 2;
        int minPos;
        if(rightChildPos<heapArray.size())                      //both right and left children are part of the heap
            minPos = (heapArray[leftChildPos]->f > heapArray[rightChildPos]->f)? rightChildPos : leftChildPos;
        else if(leftChildPos<heapArray.size())                 // only the left child is part of the heap
            minPos = leftChildPos;
        else                                                    // neither child is part of the heap, i.e. at the end of the heap, therefore return
            return;

        if(heapArray[minPos]->f < heapArray[currPos]->f){           // if the min cost child has a lower cost than parent, exchange nodes
            graphNode* tempNode = heapArray[currPos];
            heapArray[currPos] = heapArray[minPos];
            heapArray[minPos] = tempNode;
            heapArray[minPos]->heapPos = minPos;      // update heapPos parameter of the Node with its current position in the heap
            heapArray[currPos]->heapPos = currPos;
            update(minPos,DOWN);                                // recursively check if children have lower costs
        }
        return;
    }

    bool isHeapEmpty()
    {
        return(heapSize==0);
    }

    double retHeapRootVal()
    {
        return heapArray[0]->f;
    }
};


//===========================
// Class optPath
// ==========================

// Class describing the actual optimal path
class optPath{
	std::vector< double > pathX;    // x coordinates of the reference path
	std::vector< double > pathY;    // y coordinates of the reference path
	std::vector< double > pathT;    // t coordinates of the reference path

	gsl_interp_accel* accX;     // accelerator objects required by the gsl interpolation
	gsl_interp_accel* accY;
    gsl_spline* splineX;        // the interpolant for X coordinates
    gsl_spline* splineY;        // the interpolant for Y coordinates

    bool pathRead;              // flag to indicate that the path was read from file
    bool interpSet;             // flag to indicate that the interpolants were set

    public:
    int pathLength;             // length of reference path

    optPath(){
	    pathX.clear();
	    pathY.clear();
	    pathT.clear();
        accX = gsl_interp_accel_alloc();
        accY = gsl_interp_accel_alloc();
        pathLength = 0;
        interpSet = false;
        pathRead = false;
	}

    ~optPath(){
        gsl_interp_accel_free(accX);
        gsl_interp_accel_free(accY);
        gsl_spline_free(splineX);
        gsl_spline_free(splineY);
    }

    // function to read the path off file
    void readPathFromFile(std::string fileName){
        std::ifstream inF(fileName);
        double value;
        int k = 0;
        int ptr = 0;

        while( inF>>value){
            ptr = (k++)%3;

            switch( ptr ){
                case 0:
                    pathX.push_back( value );
                    break;
                case 1:
                    pathY.push_back( value );
                    break;
                case 2:
                    pathT.push_back( value );
                    break;
            }
        }
        pathLength = pathT.size();
        pathRead = true;
    }

    // function to set the interpolants using the data
    void createInterpolant(){
        if(!pathRead){
            std::cout << "Paths not read off file. Run readPathFromFile(fileName) function first."<<std::endl;
            return;
        }
        const int pL = pathLength;
        double pathXVec[pL];
        double pathYVec[pL];
        double pathTVec[pL];

        splineX = gsl_spline_alloc( gsl_interp_cspline, pathLength );
        splineY = gsl_spline_alloc( gsl_interp_cspline, pathLength );

        for( int i=0; i<pL; i++){
            pathXVec[i] = pathX[i];
            pathYVec[i] = pathY[i];
            pathTVec[i] = pathT[i];
        }

        gsl_spline_init( splineX, pathTVec, pathXVec, pathLength );
        gsl_spline_init( splineY, pathTVec, pathYVec, pathLength );

        interpSet = true;
    }

    // function to return interpolated values
    void getInterpVal(double &x, double &y, double t){
        if(!interpSet){
            std::cout << "interpolants not set. Run createInterpolant() first" << std::endl;
            return;
        }
        x = gsl_spline_eval( splineX, t, accX );
        y = gsl_spline_eval( splineY, t, accY );

    }

    // function to return the mean error from the reference path
    double getMeanPathError( std::vector< std::vector< double > > path ){
        int pL = path.size();
        double interpX, interpY;
        double accErr = 0;
        for (int i=0; i<pL; i++){
            if ( path[i][2] > pathT[pathLength-1])
                continue;

            getInterpVal( interpX, interpY, path[i][2] );
            accErr = accErr + sqrt( ( interpX-path[i][0] )*( interpX-path[i][0] ) + ( interpY-path[i][1] )*( interpY-path[i][1]) );
        }
        return accErr/pL;
    }

    // returns the specfied coordinate from the reference path
    double getPathVal(int i, int j){
        switch(j){
            case 0:
                return pathX[i];
                break;
            case 1:
                return pathY[i];
                break;
            case 2:
                return pathT[i];
                break;
        }
    }

};

#endif // CLASSDEFS_H_INCLUDED
