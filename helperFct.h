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

vector< vector< vector<double> > > readDataToVecs(string inFname, int nx, int ny, int nt){
    vector< vector< vector<double> > > dataMat;
    vector< vector< double > > dataRows;
    vector<double> data;
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

//=======================================
// Bi-Linear Interpolation function
//=======================================
double getInterpVal(double lx, double ly, double val00, double val10, double val01, double val11){
    double intVal0 = lx*( val10 - val00 ) + val00;
    double intVal1 = lx*( val11 - val01 ) + val01;
    return ly*( intVal1 - intVal0 ) + intVal0;
}

//=====================================================
// Functions to interpolate data from the vectors
//===================================================== 
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
    double delX = (goalx - nx);
	double delY = (goaly - ny);

	double dt = sqrt(delX*delX + delY*delY)/(Vfm + vMinH);
    return (k1 + k2*pow(vMinH,alpha) )*dt;
}

//=============================================
// return the planned path in a vector
//=============================================
vector< graphNode* > getPlannedPath(graphNode* currNodePtr)
{
    vector<graphNode*> path;
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

double getFlowMag(double x, double y, double t){
    double vx, vy;
    getFlowFromVecs(x,y,t,vx,vy);
    return sqrt(vx*vx+vy*vy);
}

//==================================
// Function to find max dX and dT
//==================================
void findDx(double x1, double y1, double t1, double v1, double v2, double v3, double vF0, double vSel, double &xr, double &yr, double &tr){
    if ( x1>=xmax || y1>=ymax || x1<=xmin || y1<=ymin || t1>=tmax || t1<=tmin ){
        xr = x1; yr = y1; tr = t1;
        return;
    }
    double x2, y2, t2;
    
    double vDiff1 = fabs( getFlowMag(x1,y1,t1)-vF0 );
    double vDiffx, vDiffy, vDifft;
    double dfdx, dfdy, dfdt;
    double L, dfdXL;

    while ( fabs(vDiff1-p*vSel)>0.1*p*vSel ){
        x2 = (x1==0)?0.0001:0.9999*x1;
        y2 = (y1==0)?0.0001:0.9999*y1;
        t2 = (t1==0)?0.0001:0.9999*t1;
        
        vDiffx = fabs( getFlowMag(x2,y1,t1)-vF0 );
        vDiffy = fabs( getFlowMag(x1,y2,t1)-vF0 );
        vDifft = fabs( getFlowMag(x1,y1,t2)-vF0 );

        dfdx = (vDiffx - vDiff1)/(x2-x1);
        dfdy = (vDiffy - vDiff1)/(y2-y1);
        dfdt = (vDifft - vDiff1)/(t2-t1);
        
        dfdXL = ( v1*dfdx + v2*dfdy + v3*dfdt );
        if ( fabs(dfdXL) < 0.0001 ){
            xr = x1; yr = y1; tr = t1;
            return;
        }
        L = ( p*vSel-vDiff1 )/dfdXL;

        x1 = v1*L + x1;
        y1 = v2*L + y1;
        t1 = v3*L + t1;

        vDiff1 = fabs( getFlowMag(x1,y1,t1)-vF0 );
    }
    
    xr = x1; yr = y1; tr = t1;
    return;
}

//=========================
// Function to find max dX
//=========================
void findDx(double x1, double y1, double t1, double v1, double v2, double vF0, double vSel, double &xr, double &yr){
    if ( x1>=xmax || y1>=ymax || x1<=xmin || y1<=ymin){ 
        xr = x1; yr = y1;
        return;
    }
    double vDiff1 = fabs( getFlowMag( x1, y1, t1 ) - vF0 );
    double dfdx, dfdy, dfdXL, L;
    double x2, y2;

    while ( fabs(vDiff1-p*vSel) > 0.1*p*vSel ){
        x2 = (x1==0)?0.0001:0.9999*x1;
        y2 = (y1==0)?0.0001:0.9999*y1;

        dfdx = ( fabs( getFlowMag(x2,y1,t1) - vF0  ) - vDiff1 )/(x2-x1); 
        dfdy = ( fabs( getFlowMag(x1,y2,t1) - vF0  ) - vDiff1 )/(y2-y1);
        dfdXL = v1*dfdx + v2*dfdy;
        if ( dfdXL < 0.001 ){
            xr = x1; yr = y1;
            return;
        }

        L = ( p*vSel - vDiff1 )/dfdXL;

        x1 = x1 + v1*L;
        y1 = y1 + v2*L;

        vDiff1 = fabs( getFlowMag( x1, y1, t1 ) - vF0 );
        
    }

    xr = x1; yr = y1;
    return;
}

//=========================
// Function to find max dT
//=========================
void findDt(double x1, double y1, double t1, double vF0, double vSel, double &tr){
    if ( t1>=tmax || t1<=tmin ){
        tr = t1;
        return;
    }
    
    double vDiff1 = fabs( getFlowMag( x1, y1, t1 ) - vF0 );
    double dfdt;
    double t2;

    while ( fabs(vDiff1-p*vSel) >= 0.1*p*vSel ){
        t2 = (t1==0)?0.0001:0.9999*t1;

        dfdt = ( fabs( getFlowMag(x1,y1,t2) - vF0  ) - vDiff1 )/(t2-t1); 
        if ( dfdt < 0.001 ){
            tr = t1;
            return;
        }

        t1 = t1 + ( p*vSel - vDiff1 )/dfdt;
        vDiff1 = fabs( getFlowMag( x1, y1, t1 ) - vF0 );
        
    }

    tr = t1;
    return;
}


//==========================================
// Function to determine node admissibility
//==========================================

bool isReachable(double x, double y, double t){
    double dt = sqrt( (goalx-x)*(goalx-x) + (goaly-y)*(goaly-y) )/( Vfm + Vm);
    if ( (t+dt)<tmax )
        return true;
    else
        return false;
}
////---------------------------
//// to_string for windows
////---------------------------
//template <typename T>
//string to_string(T value)
//{
//	ostringstream os ;
//	os << value ;
//	return os.str() ;
//}

#endif
