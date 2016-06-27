/*

universe.h
Header File for Universe class


(c) Neale Samways 2006-9

*/


#ifndef UNIVERSE_H
#define UNIVERSE_H

#include<sys/types.h>
#include<sys/stat.h>
#include<sys/unistd.h>

#include<vector>
#include<list>

#include"global.h"
#include"bacteria.h"
#include"environmentArtifact.h"
#include"parameter.h"

using namespace std;

class bacteria;

class universe
{

public:

	universe();
	~universe();

	void init( void );										// prime the virtual universe for simulation
	int runSimulation( void );									// commence simulation 

	// inline setters
	inline void setResultsDirectory( const char* RDName ){ strcpy( resultsDir, RDName ); }		// set the directory in which results will be written
	inline void setParametersDirectory( const char* PDName ){ strcpy( parametersDir, PDName); }	// set the directory containing parameter files
	inline void setVisualisation( bool vis ){ visOn = vis; }					// set the interactive visualisation

	inline void setVerbosity( int verb ){ verbosity = verb; } 					// set the simulation verbosity (level of results reporting)
	inline void setMaxEpochs( int overrideEpochNum ){ maxEpochs = overrideEpochNum; }		// set the maximum simulation time (from command line)

private:

	// functions
	bool setOutStreams( void );									
	void setHeaders( void );
	bool tick( void );
	void refresh( int, int );
	void deathRoutine( int );
	void startGP( void );
	void outputStats( void );
	bool setupEnvironment( char* );
	void removeComments( ifstream* );

	// variables
	vector< bacteria* > populationPtr;		// vector of pointers to every population member
	vector< double* > orientation;			// vector of all population member locations and orientations
	vector< int > deadMembers;			// vector of population members that have died during an epoch
	vector< int > populationStack;			// vector register of all living members - for population randomisation

	char currWD[ FILENAME_MAX ];			// directory in which simulation invoked
	char resultsDir[ FILENAME_MAX ];		// results directory
	char parametersDir[ FILENAME_MAX ];		// parameters directory
	char absResDir[ FILENAME_MAX ];			// absolute reference to results directory
	char absParamDir[ FILENAME_MAX ];		// absolute reference to parameters directory

	environmentArtifact *environmentPtr;		// pointer to all environmental artifacts (nutrients etc.)	

	int populationSize,		// population size at current iteration
		maxEpochs,		// number of epochs the simulation will run for
		numEnvArtifacts,	// number of environmental artifacts (i.e. nutrients, toxins etc.)
		currentEpoch,			// elapsed iterations since initialisation
		numStrains,		// number of separate bacterial strains
		verbosity,		// output verbosity
		deathCount,		// number of deaths per epoch
		memberSelection;	// index of member selected for update

	int behaviour[ TOTAL_PROTEINS ];		// Tally of behaviours undertaken by population per epoch
	int globalEnvParameters[ 3 ];	// global environment details ( see enum below )

	double		nutrientPerEpoch,
			worldWidth,
			worldHeight;

	bool visOn;					// Toggle for interactive population visualisation
	
	parameter	simulationParameters;		// The parameters for the simulation (loaded from parameter file)


	ofstream 	*popSnapshotPtr,		// pointers to OFSTREAMS for outputs
			*environmentSnapshotPtr,
			*individualDataPtr,		// pointer to the individual (BACTERIA) data file
			*individualLocationsPtr,	// pointer to file for individual locations
			*resultsPtr;			// pointer to global results file

	FILE*	gnuPipePtr;			// pipe for population activity
	FILE*	gp2Ptr;				// pipe for nutrient activity

};

enum
{
	WORLD_WIDTH = 0,			//	(0)
	WORLD_HEIGHT,				//	(1)
	NUMBER_FACTORS				//	(2)
};


#endif
