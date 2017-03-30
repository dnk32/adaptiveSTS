#ifndef ADAPTIVE_DIS_HELPER
#define ADAPTIVE_DIS_HELPER

//===========================
// hashbin function
//===========================
int getHashBinNumber(int x, int y, int t){
    int xBin = ( ( (x-xmin)%xSpBin )/dxBin );
    int yBin = ( ( (y-ymin)%ySpBin )/dyBin );
    int tBin = ( ( (t-tmin)%tSpBin )/dtBin );

    return (xBin + yBin*nXBins + tBin*nXBins*nYBins);
}

//========================================
// function to read data files as vectors
//========================================

std::vector< std::vector< std::vector<double> > > readDataToVecs(string inFname, int nx, int ny, int nt){
    std::vector< std::vector< std::vector<double> > > dataMat;
    std::vector< std::vector< double > > dataRows;
    std::vector<double> data;
    dataMat.clear(); dataRows.clear(); data.clear();

    ifstream inF(inFname);
    cout << "Reading File " << inFname << endl;
    int q = 0;
    double val = 0;
    int N = nx*ny;
    while( inF >> val ){
        q++;
        data.push_back(val);
        if( !(q%ny) ){
            dataRows.push_back(data);
            data.clear();
        }
        if( !(q%N) ){
            dataMat.push_back(dataRows);
            dataRows.clear();
        }
        if( !( q%( ( N*nt )/10 ) ) )
            cout <<  q/(double)( N*nt )*100 << "% of the data read " << endl; 
    }
    return dataMat;
}

//=======================================
// Tri-Linear Interpolation function
//=======================================
double getInterpVal(double lx, double ly, double lt, double val000, double val100, double val010, double val110, double val001, double val101, double val011, double val111){
    double intVal00 = lx*( val100 - val000 ) + val000;
    double intVal01 = lx*( val110 - val010 ) + val010;
    double intVal0 = ly*( intVal01 - intVal00 ) + intVal00;

    intVal00 = lx*( val101-val001 ) + val001;
    intVal01 = lx*( val111-val011 ) + val011;
    double intVal1 = ly*( intVal01 - intVal00 ) + intVal00;

    return lt*( intVal1-intVal0) + intVal0;
}

// interpolate flow velocity using the available data
void getFlowFromVecs(double x,double y,double t,double &vx,double &vy)
{
    /* Inputs
        x : x cordinate in the grid
        y : y cordinate in the grid
        t : t cordinate in the grid
        vx : flow velocity in the x direction
        vy : flow velocity in the y direction
    */

    /* using interpolation
    -----------------------*/
    int nx = x/xSpData;
    int ny = y/ySpData;
    int nt = t/tSpData;

    if( nx == (nX-1) || ny == (nY-1) || nt == (nT-1) ){
        vx = vXVec[nt][nx][ny];
        vy = vYVec[nt][nx][ny];
        return;
    }
    double u1 = (double)( x-nx*xSpData )/xSpData*( vXVec[nt][nx+1][ny]-vXVec[nt][nx][ny] ) + vXVec[nt][nx][ny];
    double u2 = (double)( x-nx*xSpData )/xSpData*( vXVec[nt][nx+1][ny+1]-vXVec[nt][nx][ny+1] ) + vXVec[nt][nx][ny+1];
    double u3 = (double)( x-nx*xSpData )/xSpData*( vXVec[nt+1][nx+1][ny]-vXVec[nt+1][nx][ny] ) + vXVec[nt+1][nx][ny];
    double u4 = (double)( x-nx*xSpData )/xSpData*( vXVec[nt+1][nx+1][ny+1]-vXVec[nt+1][nx][ny+1] ) + vXVec[nt+1][nx][ny+1];

    double u12 = (double)(y-ny*ySpData)/ySpData*(u2-u1) + u1;
    double u34 = (double)(y-ny*ySpData)/ySpData*(u4-u3) + u3;

    vx = (double)(t-nt*tSpData)/tSpData*(u34-u12) + u12;

    u1 = (double)( x-nx*xSpData )/xSpData*( vYVec[nt][nx+1][ny]-vYVec[nt][nx][ny] ) + vYVec[nt][nx][ny];
    u2 = (double)( x-nx*xSpData )/xSpData*( vYVec[nt][nx+1][ny+1]-vYVec[nt][nx][ny+1] ) + vYVec[nt][nx][ny+1];
    u3 = (double)( x-nx*xSpData )/xSpData*( vYVec[nt+1][nx+1][ny]-vYVec[nt+1][nx][ny] ) + vYVec[nt+1][nx][ny];
    u4 = (double)( x-nx*xSpData )/xSpData*( vYVec[nt+1][nx+1][ny+1]-vYVec[nt+1][nx][ny+1] ) + vYVec[nt+1][nx][ny+1];

    u12 = (double)(y-ny*ySpData)/ySpData*(u2-u1) + u1;
    u34 = (double)(y-ny*ySpData)/ySpData*(u4-u3) + u3;

    vy = (double)(t-nt*tSpData)/tSpData*(u34-u12) + u12;

}

void getFlowFromVecs(double x, double y, double t, double &vX, double &vY, double &dvXdx, double &dvXdy, double &dvYdx, double &dvYdy, double &dvXdt, double &dvYdt){
   
    int nx = x/xSpData;
    int ny = y/ySpData;
    int nt = t/tSpData;
    
    // if we are at the boundary in the x-direction
    if ( nx == (nX-1) ){
        dvXdx = ( vXVec[nt][nx][ny] - vXVec[nt][nx-1][ny] )/xSpData;
        dvYdx = ( vYVec[nt][nx][ny] - vYVec[nt][nx-1][ny] )/xSpData;
    }
    else if( nx == 0) {
        dvXdx = ( vXVec[nt][nx+1][ny] - vXVec[nt][nx][ny] )/xSpData;
        dvYdx = ( vYVec[nt][nx+1][ny] - vYVec[nt][nx][ny] )/xSpData;
    }
    else{
        dvXdx = ( vXVec[nt][nx+1][ny] - vXVec[nt][nx-1][ny] )/(2*xSpData);
        dvYdx = ( vYVec[nt][nx+1][ny] - vYVec[nt][nx-1][ny] )/(2*xSpData);
    }
    
    // if at the boundary in the y-direction
    if ( ny == (nY-1) ){
        dvXdy = ( vXVec[nt][nx][ny] - vXVec[nt][nx][ny-1] )/ySpData;
        dvYdy = ( vYVec[nt][nx][ny] - vYVec[nt][nx][ny-1] )/ySpData;
    }
    else if( ny == 0) {
        dvXdy = ( vXVec[nt][nx][ny+1] - vXVec[nt][nx][ny] )/ySpData;
        dvYdy = ( vYVec[nt][nx][ny+1] - vYVec[nt][nx][ny] )/ySpData;
    }
    else{
        dvXdy = ( vXVec[nt][nx][ny+1] - vXVec[nt][nx][ny-1] )/(2*ySpData);
        dvYdy = ( vYVec[nt][nx][ny+1] - vYVec[nt][nx][ny-1] )/(2*ySpData);
    }
    
    // if at the boundary in the t-direction
    if ( nt == (nT-1) ){
        dvXdt = ( vXVec[nt][nx][ny] - vXVec[nt-1][nx][ny] )/tSpData;
        dvYdt = ( vYVec[nt][nx][ny] - vYVec[nt-1][nx][ny] )/tSpData;
    }
    else if( nt == 0) {
        dvXdt = ( vXVec[nt+1][nx][ny] - vXVec[nt][nx][ny] )/tSpData;
        dvYdt = ( vYVec[nt+1][nx][ny] - vYVec[nt][nx][ny] )/tSpData;
    }
    else{
        dvXdt = ( vXVec[nt+1][nx][ny] - vXVec[nt-1][nx][ny] )/(2*tSpData);
        dvYdt = ( vYVec[nt+1][nx][ny] - vYVec[nt-1][nx][ny] )/(2*tSpData);
    }

    // if at the maximum boundary of data
    if( nx == (nX-1) || ny == (nY-1) || nt == (nT-1) ){
        vX = vXVec[nt][nx][ny];
        vY = vYVec[nt][nx][ny];
        return;
    }
    double lx = (double)(x-nx*xSpData)/xSpData;
    double ly = (double)(y-ny*ySpData)/ySpData;
    double lt = (double)(t-nt*tSpData)/tSpData;

    vX = getInterpVal( lx, ly, lt, vXVec[nt][nx][ny], vXVec[nt][nx+1][ny], vXVec[nt][nx][ny+1], vXVec[nt][nx+1][ny+1], vXVec[nt+1][nx][ny], vXVec[nt+1][nx+1][ny], vXVec[nt+1][nx][ny+1], vXVec[nt+1][nx+1][ny+1]); 
    
    vY = getInterpVal( lx, ly, lt, vYVec[nt][nx][ny], vYVec[nt][nx+1][ny], vYVec[nt][nx][ny+1], vYVec[nt][nx+1][ny+1], vYVec[nt+1][nx][ny], vYVec[nt+1][nx+1][ny], vYVec[nt+1][nx][ny+1], vYVec[nt+1][nx+1][ny+1]);

    return;
}

//=====================================================================================
// function to return if a node is on an obstacle or if it is beyond the graph limits
//=====================================================================================
bool isAccessible(int nx,int ny,int nt)
{
    int nxorg = nx/xSpData;
    int nyorg = ny/ySpData;

    if ( nx>=xmin && nx<=xmax && ny>=ymin && ny<=ymax && nt>=tmin && nt<=tmax){
		if( OBS[0][nxorg][nyorg] == 0)
            return true;
    }
	else
	{
		return false;
	}
}
//=============================================
// function to return the heuristic to target
//=============================================
double getHeuristic(int nx,int ny)
{
//    return 0.0;
    double delX = (GOAL_COORD[0] - nx);
	double delY = (GOAL_COORD[1] - ny);

	double dt = sqrt(delX*delX + delY*delY)/(Vfm + vMinH);
    return (k1 + k2*pow(vMinH,alpha) )*dt;
}

//=============================================
// return the planned path in a vector
//=============================================
std::vector< graphNode* > getPlannedPath(graphNode* currNodePtr)
{
    std::vector<graphNode*> path;
    path.clear();
    //path.push_back(currNodePtr);
    //graphNode* parentNodePtr = currNodePtr->parent;
    graphNode* parentNodePtr = currNodePtr;
    while(parentNodePtr)
    {
        path.push_back(parentNodePtr);
        parentNodePtr = parentNodePtr->parent;
    }
    return path;
}

bool isDirAccessible(double vF, double th){
    if ( fabs(th)<=M_PI/2 ){
        if (Vm < vF*fabs( sin(th) ) )
            return false;
    }
    else{
        if ( Vm < vF )
            return false;
    }
    return true;
}

//-----------------------------
// Function to find max dX
//-----------------------------
void findDx(double x1, double y1, double x2, double y2, double t1, double t2, double vF0, double vSel, double &xr, double &yr, double &tr){
    double v1, v2, vx, vy;
    getFlowFromVecs( x1,y1,t1,vx,vy );
    v1 = sqrt( vx*vx+vy*vy );
    getFlowFromVecs( x2,y2,t2,vx,vy );
    v2 = sqrt( vx*vx+vy*vy );
    
    double xt, yt, tt;

    double vDiff1 = fabs(v1-vF0);
    double vDiff2 = fabs(v2-vF0);
    if ( vDiff1 < p*vSel){
        xt = x1; yt = y1; tt = t1;
        x1 = x1 + x1-x2; y1 = 2*y1-y2; t1 = 2*t1-t2;
        x2 = xt; y2 = yt; t2 = tt;
        findDx(x1,y1,x2,y2,t1,t2,vF0,vSel,xr,yr,tr);
        return;
    }
    else if (vDiff2 > p*vSel){
        xt = x2; yt = y2; tt = t2;
        x2 = 2*x2-x1; y2 = 2*y2-y1; t2 = 2*t2-t1;
        x1 = xt; y1 = yt; t1 = tt;
        findDx(x1,y1,x2,y2,t1,t2,vF0,vSel,xr,yr,tr);
        return;
    }
    else {
        xt = (x1+x2)/2; yt = (y1+y2)/2; tt = (t1+t2)/2;
        getFlowFromVecs(xt,yt,tt,vx,vy);
        v1 = sqrt( vx*vx+vy*vy );
        vDiff1 = fabs(v1-vF0);
        if ( ( vDiff1 <= (1.1*p*vSel) ) && ( vDiff1 >= (0.9*p*vSel) ) ){
            xr = xt;
            yr = yt;
            tr = tt;
            return;
        }
        else if ( vDiff1 > p*vSel){
            x1 = xt; y1 = yt; t1 = tt;
            findDx(x1,y1,x2,y2,t1,t2,vF0,vSel,xr,yr,tr);
            return;
        }
        else {
            x2 = xt; y2 = yt; t2 = tt;
            findDx(x1,y1,x2,y2,t1,t2,vF0,vSel,xr,yr,tr);
            return;
        }

    }
}

void findDx(double x1, double y1, double x2, double y2, double t, double vF0, double vSel, double &xr, double &yr){
    if ( x1>=xmax || y1>=ymax || x1<=xmin || y1<=ymin){
        xr = x1; yr = y1;
        return;
    }
    double v1, v2, vx, vy;
    getFlowFromVecs( x1,y1,t,vx,vy );
    v1 = sqrt( vx*vx+vy*vy );
    getFlowFromVecs( x2,y2,t,vx,vy );
    v2 = sqrt( vx*vx+vy*vy );
    
    double xt, yt;

    double vDiff1 = fabs(v1-vF0);
    double vDiff2 = fabs(v2-vF0);
    if ( vDiff1 < p*vSel){
        xt = x1; yt = y1;
        x1 = x1 + x1-x2; y1 = 2*y1-y2;
        x2 = xt; y2 = yt;
        findDx(x1,y1,x2,y2,t,vF0,vSel,xr,yr);
        return;
    }
    else if (vDiff2 > p*vSel){
        xt = x2; yt = y2;
        x2 = 2*x2-x1; y2 = 2*y2-y1;
        x1 = xt; y1 = yt;
        findDx(x1,y1,x2,y2,t,vF0,vSel,xr,yr);
        return;
    }
    else {
        xt = (x1+x2)/2; yt = (y1+y2)/2;
        getFlowFromVecs(xt,yt,t,vx,vy);
        v1 = sqrt( vx*vx+vy*vy );
        vDiff1 = fabs(v1-vF0);
        if ( ( vDiff1 <= (1.1*p*vSel) ) && ( vDiff1 >= (0.9*p*vSel) ) ){
            xr = xt;
            yr = yt;
            return;
        }
        else if ( vDiff1 > p*vSel){
            x1 = xt; y1 = yt;
            findDx(x1,y1,x2,y2,t,vF0,vSel,xr,yr);
            return;
        }
        else {
            x2 = xt; y2 = yt;
            findDx(x1,y1,x2,y2,t,vF0,vSel,xr,yr);
            return;
        }

    }
}

void findDt(double x, double y, double t1, double t2, double vF0, double vSel, double &tr){
    if ( t1>=tmax || t1<=tmin ){
        tr = t1;
        return;
    }
    double v1, v2, vx, vy;
    getFlowFromVecs( x,y,t1,vx,vy );
    v1 = sqrt( vx*vx+vy*vy );
    getFlowFromVecs( x,y,t2,vx,vy );
    v2 = sqrt( vx*vx+vy*vy );
    
    double tt;

    double vDiff1 = fabs(v1-vF0);
    double vDiff2 = fabs(v2-vF0);
    if ( vDiff1 < p*vSel){
        tt = t1;
        t1 = 2*t1-t2;
        t2 = tt;      
        findDt(x,y,t1,t2,vF0,vSel,tr);
        return;
    }
    else if (vDiff2 > p*vSel){
        tt = t2;
        t2 = 2*t2-t1;
        t1 = tt;
        findDt(x,y,t1,t2,vF0,vSel,tr);
        return;
    }
    else {
        tt = (t1+t2)/2;
        getFlowFromVecs(x,y,tt,vx,vy);
        v1 = sqrt( vx*vx+vy*vy );
        vDiff1 = fabs(v1-vF0);
        if ( ( vDiff1 <= (1.1*p*vSel) ) && ( vDiff1 >= (0.9*p*vSel) ) ){
            tr = tt;
            return;
        }
        else if ( vDiff1 > p*vSel){
            t1 = tt;
            findDt(x,y,t1,t2,vF0,vSel,tr);
            return;
        }
        else {
            t2 = tt;
            findDt(x,y,t1,t2,vF0,vSel,tr);
            return;
        }

    }
}

////---------------------------
//// to_string for windows
////---------------------------
//template <typename T>
//std::string to_string(T value)
//{
//	std::ostringstream os ;
//	os << value ;
//	return os.str() ;
//}

/*
******************************
* Unused Functions
******************************

//================================================
// function to get number of neighbors considered
//================================================

int getNumNeighbs(){
    int gridPtr = -1;
    for (int i=0; i<=2*nDivs ; i++){
        for (int j=0; j<=( 2*nDivs - abs(i-nDivs) ); j++){
            gridPtr++;
        }
    }
    return gridPtr;
}
//=============================================
// Recursive Function to compute costs
//=============================================

void recursiveComputeNeighb(int k, double xSt, double ySt, double tSt, double xGl, double yGl, double Est, double dt, int Nt, double delV, vector<double> &E, vector<bool> &reachable, double &dtStatic){
    if ( ( k ) && ( k>Nt ) ){ // return if the k (level iterator) > Nt (max number of levels)
        //cout << "exit condition 1 reached" << endl;
        return;
    }

    // get flow velocities
    double vX, vY;
    getDGFlow(xSt,ySt,tSt,vX,vY);

    double vF = sqrt( vX*vX + vY*vY );
    double thF = atan2(vY,vX);                // flow direction
    double thHdg = atan2( yGl-ySt, xGl-xSt ); // heading direction

    double th = thHdg - thF;                  // heading relative to the flow

    // compute min/max speeds and transition times
    double vMinReq, vMin, vMax, vK;
    double dtMin, dtMax;
    double delX = sqrt( (yGl-ySt)*(yGl-ySt) + (xGl-xSt)*(xGl-xSt)  );
    if (fabs(th)<=M_PI/2){
        vMinReq = vF*sin( fabs(th) );
        vMin = max( 0.0  , vF*cos(th) - sqrt(Vm*Vm - vF*vF*sin(th)*sin(th) )  );
    }
    else{
        vMinReq = vF;
        vMin = 0.0;
    }

    if (vMinReq>Vm){       // if minimum required velocity is greater than the max vehicle speed
        //cout << "exit condition 2 reached" << endl;
        return;           // exit function
    }

    vMax = vF*cos(th) + sqrt( Vm*Vm - vF*vF*sin(th)*sin(th) );

    dtMax = delX/vMin;
    dtMin = delX/vMax;

    if (dtMax < dt){       // if the required direction cannot be reached within dt
        //cout << "exit condition 3 reached" << endl;
        if (!k){
            double dtStaticT = sqrt( k2*(delX)/( k1+k2*vF*vF ) );
            dtStatic = dtStaticT;
            if (dtStaticT>dtMax)
                dtStatic = dtMax;
            if (dtStaticT<dtMin)
                dtStatic = dtMin;

            double vStatic = sqrt( ( delX/dtStatic-vF*cos(th) )*( delX/dtStatic-vF*cos(th) ) + vF*sin(th)*vF*sin(th) );
            double cost = Est + (k1 + k2*pow(vStatic,alpha))*dtStatic/tConvRat;

            reachable.push_back(true);
            E.push_back(cost);
        }
        return;           // exit function
    }

    if (dtMin < dt) {
        double vNet = delX/dt;
        vK = sqrt( ( vNet-vF*cos(th) )*( vNet-vF*cos(th) ) + vF*sin(th)*vF*sin(th) );
        double cost = Est + (k1 + k2*pow(vK,alpha))*dt/tConvRat;

        if ( reachable.size() < (k+1) ){  // resize the cost and reachable vectors to hold upto k+1 vals
            reachable.resize(k+1);
            E.resize(k+1);
            reachable[k] = false;
        }
        if ( !reachable[k] )        // if the current node has not been reached before
            E[k] = cost;
        if ( E[k]>cost )           // if the current cost is less than the previous cost
            E[k] = cost;
        reachable[k] = true;
        vMax = vNet;
    }
    else{
        if (vMax*dt <= delX/5)
            return;
    }

    if (!k){                    // if this is the initial iteration, set Nt and delV
        if (dtMin<=dt)
            delV = (vMax-vMin)/nV;
        else
            delV = (vMax-vMin)/2;

        if ( isinf(dtMax) )
            Nt = NtMax;
        else
            Nt = min( (int)(dtMax/dt), NtMax );
    }
    double vi = (vMin>0)?vMin:delV;     // in vMin=0, set vStart = delV
    while(vi <= vMax){
        double xInt = xSt + vi*cos(thHdg)*dt;
        double yInt = ySt + vi*sin(thHdg)*dt;
        double tInt = tSt + dt;

        vK = sqrt( ( vi-vF*cos(th) )*( vi-vF*cos(th) ) + vF*sin(th)*vF*sin(th) );
        double EInt = Est + (k1 + k2*pow(vK,alpha))*dt/tConvRat;

        recursiveComputeNeighb(k+1,  xInt,  yInt,  tInt,  xGl,  yGl,  EInt, dt, Nt,  delV,  E, reachable, dtStatic);
        vi = vi + delV;
    }
    //cout << "exit condition 4 reached" << endl;
    return;

}
*/
#endif
