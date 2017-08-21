#ifndef GRAPHNODE_HPP
#define GRAPHNODE_HPP
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

    graphNode(int xx, int yy, int tt, double gg, double hh, graphNode* parentNode, double time):
        x(xx), y(yy), t(tt),
        g(gg),
        h(hh),
        f(g+h),
        parent(parentNode),
        expanded(false),
        timeAtNode(time),
        pathLength(0)
    {}
};

#endif
