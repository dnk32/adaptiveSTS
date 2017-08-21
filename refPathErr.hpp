#ifndef REFPATHERR_H
#define REFPATHERR_H

#include <gsl/gsl_spline.h>
#include <vector>
#include <string>
#include <fstream>

using namespace std;
//===========================
// Class refPathErr
// ==========================

// Class describing the actual reference path
class refPathErr{
	vector< double > pathX;    // x coordinates of the reference path
	vector< double > pathY;    // y coordinates of the reference path
	vector< double > pathT;    // t coordinates of the reference path

	gsl_interp_accel* accX;     // accelerator objects required by the gsl interpolation
	gsl_interp_accel* accY;
    gsl_spline* splineX;        // the interpolant for X coordinates
    gsl_spline* splineY;        // the interpolant for Y coordinates

    bool pathRead;              // flag to indicate that the path was read from file
    bool interpSet;             // flag to indicate that the interpolants were set

    int pathLength;             // length of reference path
    
    void createInterpolant();
    void readData( string fileName, vector<double> *X, vector<double> *Y, vector<double> *T); 
    public:

    refPathErr(){
	    pathX.clear();
	    pathY.clear();
	    pathT.clear();
        accX = gsl_interp_accel_alloc();
        accY = gsl_interp_accel_alloc();
        pathLength = 0;
        interpSet = false;
        pathRead = false;
	}

    ~refPathErr(){
        gsl_interp_accel_free(accX);
        gsl_interp_accel_free(accY);
        gsl_spline_free(splineX);
        gsl_spline_free(splineY);
    }

    void readPathFromFile(string fileName);
    double getMeanPathError( vector< vector< double > > *path );
    void getInterpVal(double &x, double &y, double t);
    double getPathVal(int i, int j);
    double getMeanPathError( string fileName );
};

#endif
