#include "HeapContainer.hpp"

using namespace std;

//push a graph node pointer to the heap and sort heap so that the minimum cost node is at the root
void HeapContainer::pushNode(graphNode* ptrN){
    heapArray.push_back(ptrN);                     
    heapSize++;
    int newPos = heapArray.size()-1;                // position of new value
    heapArray[newPos]->heapPos = newPos;            
    update(newPos,UP);                              
}

//return the graphNode with the least cost
graphNode* HeapContainer::popNode(){
    graphNode* ret = heapArray[0];                  // return the minimum cost node
    heapArray[0] = heapArray[heapArray.size()-1];   // put last node in the first node position
    heapArray[0]->heapPos = 0;
    heapArray.pop_back();                           // remove last node in the vector
    update(0,DOWN);
    heapSize--;
    return ret;
}

//sort the heap. dir is UP for push and DOWN for pop,  BOTH for change of value operations
void HeapContainer::update(int currPos,BubbleDirection dir){
    if(dir!=DOWN){                      // a update from a push operation or a change opration
        if(currPos==0 && dir==UP)       // currently at the root, sorting complete
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

bool HeapContainer::isHeapEmpty(){
    return(heapSize==0);
}

double HeapContainer::retHeapRootVal(){
    return heapArray[0]->f;
}
