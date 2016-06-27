/*

Lanscape.h
Header File for environmentArtifact class

(c) Neale Samways 2006-9

*/



#ifndef ENVIRONMENTARTIFACT_H
#define ENVIRONMENTARTIFACT_H

#include<iostream>
#include<cstdlib>

#include "global.h"

using namespace std;

class environmentArtifact
{

public:

	environmentArtifact();					// constructor
	~environmentArtifact();					// destructor

// setters
	void setDimensions( int, int );				// set the dimensions of the artifact
	void setParameter( int, char* );			// set the indexed parameter in the class
	bool init( void );					// initialise the artifact 
	double localUpdate( double, double, double );		// update the concentration at a given locus
	void globalUpdate( int );				// update the entire artifact (i.e. diffuse, possible deployment)
	void deploy( int, int, int, double, double, double );	// make a deployment of artifact according to parameters
	void deposit( int, int, double );			// increase local concentration at the identified locus
	void decay( double );					// globally reduce the concentration of artifact
	void move( int, int );					// move all concentrations in the concentration matrix by the identified amount
	void primaryDeployment( double, double );		// make a deployment at specified location

// getters

	double getSumVolume( void );				// return the aggregate amount of artifact in the environment
	double reportGradient( void );				// return a rating for the amount of "gradient" information in the environment

// printers

	void print( ofstream* );				// output the contents of the concentration matrix (i.e. the entire distribution)
	void displayUpdate( FILE* );				// pipe the contents of the concentration matrix for display
	void printConfiguration( ofstream* );			// write the configuration of the artifact to filestream

// Inlines
	double querySite( double x, double y ){ return( concentrationMatrix[ ( int ) floor( x ) ][ ( int ) floor( y ) ] ); }	// return artifact concentration at a given locus


private:


// functions

	void diffuse( void );						// diffuse the artifact 
	double getDistance( int, int, int, int );			// return the Euclidian distance between given points
	double getPeriodicDistance( int, int, int, int );		// return the Euclidian distance between given points assuming periodic boundaries

//variables

	static const int numParameters;			// number of individual parameters

	char friendlyName[ FILENAME_MAX ];		// human readable label
	
	int worldWidth,			// dimensions of environmentArtifact
		worldHeight,
		worldEpoch,		// elapsed time
		xOffset,		// offset for calculating deployment centriod (x dimension)
		yOffset;		// offset for calculating deployment centriod (y dimension)

	double dConstant;					// diffusion constant (i.e. how quickly the artifact diffuses)
	double **concentrationMatrix;				// pointer to environmentArtifact matrix
	double parameters[ NUMBER_ENVIRONMENT_PARAMS ];		// parameter array (double-casted)
	double deploymentVolume;				// sum volume of nutrient per deployment

	double currentAmp;					// used in calculation of dynamic width deployments
	double currentSigma;					// used in calculation of dynamic width deployments

};

// Environmental parameters
enum
{
	FRIENDLY_NAME = 0,			//	(0)
	DIFFUSIVITY_CONST,			//	(1)
	PERMEATION_FACTOR,			//	(3)
	ENERGY_YIELD,				//	(4)
	BASAL_DIST,				//	(5)
	UPDATE_PERIOD,				//	(6)
	UPDATE_RANGE,				//	(7)
	UPDATE_CONFIG,				//	(8)
	UPDATE_REGIME,				//	(9)
	DEPLOY_RADIUS,				//	(10)
	DEPLOY_MAX,				//	(11)
	DEPLOY_DELTA,				//	(12)
	DECAY_PERIOD,				//	(13)
	DECAY_CONST				//	(14)

};

#endif
