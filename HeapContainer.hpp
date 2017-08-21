#ifndef HEAPCONTAINER_HPP 
#define HEAPCONTAINER_HPP

#include <vector>
#include "graphNode.hpp"

using namespace std;
//=================================
// Class Heap
//=================================

class HeapContainer{
    public:
    vector<graphNode*> heapArray;                       //contains the heap
    enum BubbleDirection{UP, DOWN, BOTH} ;              // structure that contains the heap sorting direction, percolating up or percolating down
    int heapSize;                                       // size of the heap

    HeapContainer(){                                    //clear HeapArray on Start
        heapArray.clear();
        heapSize=0;
    }

    void pushNode(graphNode* ptrN);
    graphNode* popNode();
    void update(int currPos,BubbleDirection dir=BOTH);
    bool isHeapEmpty();
    double retHeapRootVal();
};
#endif
