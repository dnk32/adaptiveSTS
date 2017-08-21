#include "refPathErr.hpp"
#include <iostream>
#include <sstream>
#include <math.h>


//===================
// Private Functions
//===================

// general function to read data off a file
void refPathErr::readData( string fileName, vector<double> *X, vector<double> *Y, vector<double> *T){
    (*X).clear();
    (*Y).clear();
    (*T).clear();
    
    ifstream inF(fileName);
    double value;
    int k = 0;
    string line;

    while ( !inF.eof() ){
        getline( inF, line );
        stringstream ss(line);
        k = 0;
        while(ss>>value){
            switch (k++){
                case 0:
                    (*X).push_back( value );
                    break;
                case 1:
                    (*Y).push_back( value );
                    break;
                case 2:
                    (*T).push_back( value );
                    break;
            }
        }
    }
    return;
}

// function to set the interpolants using the data
void refPathErr::createInterpolant(){
    if(!pathRead){
        cout << "Paths not read off file. Run readPathFromFile(fileName) function first."<<endl;
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


//===================
// Public Functions
//===================

// function to read the path off file
void refPathErr::readPathFromFile(string fileName){
    readData(fileName, &pathX, &pathY, &pathT);
    pathLength = pathT.size();
    pathRead = true;
    createInterpolant();
}

// function to return the mean error from the reference path
double refPathErr::getMeanPathError( vector< vector< double > > *candPath ){
    int pL = (*candPath).size();
    double interpX, interpY;
    double accErr = 0;
    for (int i=0; i<pL; i++){
        if ( (*candPath)[i][2] > pathT[pathLength-1])
            continue;

        getInterpVal( interpX, interpY, (*candPath)[i][2] );
        accErr = accErr + sqrt( ( interpX-(*candPath)[i][0] )*( interpX-(*candPath)[i][0] ) + ( interpY-(*candPath)[i][1] )*( interpY-(*candPath)[i][1]) );
    }
    return accErr/pL;
}

// function to return the mean error from the reference path
double refPathErr::getMeanPathError( string fileName ){
    vector<double> X,Y,T;
    readData(fileName, &X, &Y, &T);
    int pL = X.size();
    vector <vector <double> > path( pL, vector<double>(3,0.0) );
    for (int i=0; i<pL; i++){
        path[i][0] = X[pL-1-i];
        path[i][1] = Y[pL-1-i];
        path[i][2] = T[pL-1-i];
    }
    return getMeanPathError(&path);
}

// function to return interpolated values
void refPathErr::getInterpVal(double &x, double &y, double t){
    if(!interpSet){
        cout << "interpolants not set. Run createInterpolant() first" << endl;
        return;
    }
    x = gsl_spline_eval( splineX, t, accX );
    y = gsl_spline_eval( splineY, t, accY );

}

// returns the specfied coordinate from the reference path
double refPathErr::getPathVal(int i, int j){
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
