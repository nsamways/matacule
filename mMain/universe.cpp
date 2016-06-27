/*

universe.cpp

member functions and data members for universe class

(c) Neale Samways 2006-9

*/ 

#include<iostream>
#include<fstream>
#include<ctime>
#include<cmath>
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<vector>
#include<algorithm>
#include<string>
#include<sstream>

#include"universe.h"

using namespace std;


universe::universe()							// constructor
{
	// initialise class variables

	// bools
	visOn		= 	false;

	// integers
	populationSize 	= 	0;
	currentEpoch	= 	0;
	numStrains	=	0;
	verbosity 	=	0;
	maxEpochs	=	0;

	// doubles
	worldWidth	=	0.0;
	worldHeight	=	0.0;

	// arrays
	memset( behaviour, 0, sizeof behaviour );
	memset( globalEnvParameters, 0, sizeof globalEnvParameters );

	// pointers
	environmentPtr = NULL;
	popSnapshotPtr = NULL;
	individualDataPtr = NULL;
	resultsPtr = NULL;

	// set default directory names
	strcpy( resultsDir, "./Results/" );
	strcpy( parametersDir, "./Parameters/" );

}


universe::~universe()							// destructor
{					
	// delete agent and associated locations
	for ( int i = 0; i < populationSize; i ++ ){

		delete populationPtr[ i ];
		delete [] orientation[ i ];	

	}

	// delete environment
	delete [] environmentPtr;
	

	// close all output files opened by universe class
	if ( verbosity & O_RESULTS ){
		resultsPtr->close();
		delete resultsPtr;
		resultsPtr = NULL;
	}

	if ( verbosity & O_POP_SNAPSHOT ){
		popSnapshotPtr->close();
		delete popSnapshotPtr;
		popSnapshotPtr = NULL;
	}

	if ( ( verbosity & O_SINGLE_IND_DATA ) || ( verbosity & O_ALL_IND_DATA )  ){
		individualDataPtr->close();
		delete individualDataPtr;
		individualDataPtr = NULL;
	}

	if ( ( verbosity & O_SINGLE_IND_LOCATION ) || ( verbosity & O_ALL_IND_LOCATION )  ){
		individualLocationsPtr->close();
		delete individualLocationsPtr;
		individualLocationsPtr = NULL;
	}

}


void universe::init( void )						// initialise the virtual universe
{	
// get and set necessary parameters, construct virtual environment and prime for simulation

	int initialisedMembers;	// used in initialising population

	double xOffset,			// calcuation of zone of initial deployment
		yOffset,
		xRange,
		yRange;

	double *valPtr;

	char tempFileName[ FILENAME_MAX ];
	char commonName[ FILENAME_MAX ];
	char variableName[ FILENAME_MAX ];
	char convertedFactor[ 5 ];
	char convertedInt[ 5 ];
	char slash[] = "/";
	char rawPopComp[ FILENAME_MAX ];					// holder for population composition string
	char *procPopCompPtr;								// holder for population composition string

	vector< int > initPopComp;								//  vector of population composition values

// setup outputs and file processing info

	// set output streams ( This includes setting the output file for the agent class )
	setOutStreams();

	// check paramters directory ends in /
	if ( parametersDir[ strlen( parametersDir ) - 1 ] != slash[0] )
		strcat( parametersDir, "/" );

// get parameters from external parameter files

	// simulation parameters
	strcpy( tempFileName, parametersDir );

	strcat( tempFileName, "simulation.param" );

	cout << "# Loading parameter file: " << tempFileName << endl;
	if ( !simulationParameters.initialise( NUMBER_SIMULATION_PARAMS, tempFileName ) ){
		cerr << "Unable to initialise simulations parameter file. Terminating" << endl;
		exit ( 1 );	
	}else{
		cout << "-> simulation parameters successfully loaded.\n" << endl;
	}

	// environment parameters
	strcpy( tempFileName, parametersDir );
	strcat( tempFileName, "environment.param" );

	if ( !setupEnvironment( tempFileName ) ){
		cerr << "Unable to initialise environment." << endl;
		exit( 1 );
	}

	// Establish initial population composition

	strcpy( rawPopComp, simulationParameters.cGet( POP_COMPOSITION ));
	procPopCompPtr = strtok( rawPopComp, ":" );

	while( procPopCompPtr != NULL){
		initPopComp.push_back( atoi( procPopCompPtr  ) );
		procPopCompPtr = strtok (NULL, ":");
	}

// now initPopComp contains initial popualation sizes and types

	// set number of strains:
	numStrains = initPopComp.size();

	// determine total population size:
	for ( int i = 0; i < numStrains; i++ )  populationSize += initPopComp[ i ];


	// create the population

	// set class level (static) variables in agent class
	agent::setParameterPath( parametersDir );				// parameters directory
	agent::setBaseOutput( resultsDir );						// base output directory
	agent::setOutputVerbosity( verbosity );					// class verbosity
	agent::setOutputPtr( individualDataPtr );					// proteome output data ptr
	agent::setEnvironmentFactors( globalEnvParameters[ NUMBER_FACTORS ] );	// number of environmental artifacts

	// create initial population  - should have already set: output stream in class, verbosity, number of environmental factors
	for ( int i = 0; i < populationSize; i++ ){

		bacteria *initial = new bacteria;
		populationPtr.push_back( initial );
	}

	// reset individual level flags if necessary (now constructor has been run)
	if ( verbosity & O_SINGLE_IND_DATA )
		populationPtr[ 0 ]->setIndividualVerbosity( verbosity + O_ALL_IND_DATA );

	if ( verbosity & O_SINGLE_IND_PROTEOME )
		populationPtr[ 0 ]->setIndividualVerbosity( verbosity + O_ALL_IND_PROTEOME );

	if ( verbosity & O_SINGLE_IND_LOCATION )
		populationPtr[ 0 ]->setIndividualVerbosity( verbosity + O_ALL_IND_LOCATION );

	// initialise the population based on the strain proportion

	initialisedMembers = 0;

	for ( int i = 1; i <= numStrains ; i++ ){
		cout << "# Loading parameters for Strain: " << i << endl;

		for ( int j = 1; j <= initPopComp[ i - 1 ]; j++ ){
			initialisedMembers++;
			if ( !populationPtr[ initialisedMembers - 1  ]->init( i ) ){
				cerr << "Unable to initialise individual: " << i << ": Terminating..." << endl;
				exit( 1 );
			} 
		}
	}

	// set receptor types manually (i.e. set up pointer to correct member function in bacteria class)
	for ( int i = 0; i < populationSize; i++ ) populationPtr[ i ]->setupReceptor();

	cout <<  "\n->Number of strains: " << numStrains << " \t...Successfully loaded" << endl;

	cout << "->Total population size: " << populationSize << endl;



	// set initial positions and headings randomly within rectangle of dimension provided

	// calculate and check bounds
	xRange = ( simulationParameters.dGet( INIT_POP_LOC ) <= 0.0 ) ? 0.0  : simulationParameters.dGet( INIT_POP_LOC );
	xRange = ( simulationParameters.dGet( INIT_POP_LOC ) >= globalEnvParameters[ WORLD_WIDTH ] ) ? globalEnvParameters[ WORLD_WIDTH ] : simulationParameters.dGet( INIT_POP_LOC );
	yRange = ( simulationParameters.dGet( INIT_POP_LOC ) <= 0.0 ) ? 0.0  : simulationParameters.dGet( INIT_POP_LOC );
	yRange = ( simulationParameters.dGet( INIT_POP_LOC ) >= globalEnvParameters[ WORLD_HEIGHT ] ) ? globalEnvParameters[ WORLD_HEIGHT ] : simulationParameters.dGet( INIT_POP_LOC );

	xOffset = ( ( double ) globalEnvParameters[ WORLD_WIDTH ] - xRange ) / 2.0;
	yOffset = ( ( double ) globalEnvParameters[ WORLD_HEIGHT ] - yRange ) / 2.0;

	// determine location of all initial poulation members
	for ( int i = 0; i < populationSize; i++ ){

		valPtr = new double[ 3 ];

		valPtr[ X_LOCATION ] = xOffset + ( ( ( double ) rand() / RAND_MAX ) * xRange );
		valPtr[ Y_LOCATION ] = yOffset + ( ( ( double ) rand() / RAND_MAX ) * yRange );
		valPtr[ HEADING ] = ( ( double ) rand() / RAND_MAX ) * 2 * PI;

		orientation.push_back( valPtr );	

	}

	// setup primary environmental artifact deployments if necessary
	for ( int i = 0; i < globalEnvParameters[ NUMBER_FACTORS ]; i++ )
		environmentPtr[ i ].primaryDeployment( ( double ) globalEnvParameters[ WORLD_WIDTH ] / 2.0 , ( double ) globalEnvParameters[ WORLD_HEIGHT ] / 2.0 );


	// set max epochs, if not given as command line parameter
	if ( maxEpochs == 0 ) setMaxEpochs( simulationParameters .iGet( MAX_EPOCHS ) ); 

	// start the gnuPlot thread if visualisations active
	if ( visOn ) startGP();

	// record initial states in files if verbosity level correct

	sprintf( convertedInt, "%d", currentEpoch );

	strcpy( commonName, resultsDir );
	strcat( commonName, "/p" );
	strcat( commonName, convertedInt );

	if ( verbosity & O_POP_SNAPSHOT ){
		// record initial population concentrations

		strcpy( tempFileName, commonName );
		strcat( tempFileName, ".location" );

		popSnapshotPtr = new ofstream;				// NOTE: This is where these filestreams are created.
		popSnapshotPtr->open( tempFileName );

		// output all positions to this stream
		
		popSnapshotPtr->close();
	}

	if ( verbosity & O_NUTRIENT_SNAPSHOT ){

		environmentSnapshotPtr = new ofstream;

		for ( int i = 0; i < globalEnvParameters[ NUMBER_FACTORS ]; i++ ){

			sprintf( convertedFactor, "%d", i );

			strcpy( variableName, commonName );
			strcat( variableName, convertedFactor );
			strcat( variableName, ".environment" );

			environmentSnapshotPtr->open( variableName );
			
			environmentPtr[ i ].print( environmentSnapshotPtr );

			environmentSnapshotPtr->close();
		}
	}

// system should now be ready for initial simulation

}


bool universe::setupEnvironment( char *fileNamePtr )			// set up all environmental factors
{

	int currentParam = 0;

	const int numGlobalParams = 3;
	const int numLocalParams = 13;

	char rawGlobalParameter[ FILENAME_MAX ];	// This is the array for global parameters
	char rawLocalParameter[ FILENAME_MAX ];		// This is the array for local parameters (and is re-used for each artifact)

	ifstream inStream;

	cout << "# Loading parameter file: " <<  fileNamePtr << endl;

	// read in the parameter file
	inStream.open( fileNamePtr );

	// check that the file can be opened 	
	if ( inStream.fail() ){
		cerr << "Unable to open parameter file: " << fileNamePtr << endl;
		return ( false );
	}

	// read in the global parameters
	while( currentParam < numGlobalParams ){
		removeComments( &inStream );
		inStream >> rawGlobalParameter;
		// cast global parameters to integer and set in this class
		globalEnvParameters[ currentParam ] = atoi( rawGlobalParameter);
		currentParam++;
	}


	// create environment objects
	environmentPtr = new environmentArtifact[ globalEnvParameters[ NUMBER_FACTORS ] ];		// pointer to environmental artifacts (nutrients etc.)	

	for ( int i = 0; i < globalEnvParameters[ NUMBER_FACTORS ]; i++ ){
		// pass on global information
		environmentPtr[ i ].setDimensions( globalEnvParameters[ WORLD_WIDTH ] , globalEnvParameters[ WORLD_HEIGHT ] );

		// read in local information and pass to new object 

		currentParam = 0;

		while( !inStream.eof() && ( currentParam < numLocalParams ) ){
			removeComments( &inStream );
			inStream >> rawLocalParameter;  					// get parameter
			environmentPtr[ i ].setParameter( currentParam, rawLocalParameter );	// send parameter to environmental object
			currentParam++;
		}

	}

	// initialise all environments
	for ( int i = 0; i < globalEnvParameters[ NUMBER_FACTORS ]; i++ ) environmentPtr[ i ].init();

	

	inStream.close();

	cout << "-> loaded all environmental parameters: " << globalEnvParameters[ NUMBER_FACTORS ] << " environmental artifacts successfully initialised.\n" << endl;

	return ( true );

}


int universe::runSimulation( void )					// start an open-ended simulatory run
{

	cout << "\nRunning Simulation...\n" << endl;

// run the actual simulation -- environment / population / parameters etc should all be set up now
	for ( int i = 1; i <= maxEpochs; i++ )
		if ( ! tick( ) ) break;
	 
	// return number of completed ticks
	return( currentEpoch );
}


bool universe::tick( void )						// Determine behaviour during current epoch
{

// progress the simulation by a single timestep for population
	currentEpoch++;

	double localEnvironment[ globalEnvParameters[ NUMBER_FACTORS ] ];

	// zero counters and pass global information to all population members
	deathCount = 0;

	memset( behaviour, 0, sizeof behaviour );

	agent::setWorldEpoch( currentEpoch );	


// Globally update environment
	
	for ( int i = 0; i < globalEnvParameters[ NUMBER_FACTORS ]; i++ ) environmentPtr[ i ].globalUpdate( currentEpoch );	// globally update all environmental factors


// update population

	for ( int i = 0; i < populationSize; i++ ) populationStack.push_back( i );	// create the holding stack and push individuals back

	for ( int i = 0; i < populationSize; i++ ){					// loop once through the current population

// cout << "energy = " << populationPtr[ memberSelection ]->getEnergy() << endl;

		// select individual to be updated (uniform random selection)
		memberSelection = rand() % populationStack.size();

		// get local environment information and send to selected individual

		for ( int j = 0; j < globalEnvParameters[ NUMBER_FACTORS ]; j++ ){
			localEnvironment[ j ] = environmentPtr[ j ].querySite( orientation[ populationStack[ memberSelection ] ][ X_LOCATION ], orientation[ populationStack[ memberSelection ]][ Y_LOCATION ] );
		}
		
		// determine behaviour of selected individual and update environment accordingly

		// pass environment info, cycle RN
		populationPtr[ populationStack[ memberSelection ] ]->internalUpdate( localEnvironment );

		refresh( populationStack[ memberSelection ] , populationPtr[ populationStack[ memberSelection ] ]->updateBehaviour( ) );		// once behavioural period

		// update recorded locations if necessary
		if ( verbosity & O_SINGLE_IND_LOCATION )
			*individualLocationsPtr << orientation[ populationStack[ memberSelection ] ][ X_LOCATION ] << "\t" << orientation[ populationStack[ memberSelection ]][ Y_LOCATION ] << endl;


// 			*individualLocationsPtr << populationPtr[ populationStack[ memberSelection ] ]->getUniqueID() << "\t" << orientation[ populationStack[ memberSelection ] ][ X_LOCATION ] << "\t" << orientation[ populationStack[ memberSelection ]][ Y_LOCATION ] << endl;

		// expunge from the vector of non-selected population members
		populationStack.erase( populationStack.begin() + memberSelection );

	}

	// we have iterated through the current population once randomly and updated behaviour
	// dead members are stored in deadMembers, and have been deleted from memory, but their pointers remain in the population vector
	
	// update the population size with births, and expunge dead from population vector
		// at this point, the population vector size reflects births but not deaths, so can only have increased
	populationSize += behaviour[ B_REPRODUCTION ];		// add new population members

	// now cycle has finished, pointers to dead individuals can be removed from the population vector 

	// sort if more than one death ( to ensure it is only dead members that are expunged from the vector)
	if ( deathCount > 1 ) sort( deadMembers.begin() , deadMembers.end() );
	
	for ( int i = deathCount -1; i >= 0; i-- ){

		populationPtr.erase( populationPtr.begin() + deadMembers[ i ] );
		orientation.erase( orientation.begin() + deadMembers[ i ] );
		deadMembers.pop_back();
		populationSize -= 1;	
		// now dead individuals have actually been removed, we expunge vector populationPtr 
	}

//	output population, birth and death statistics to file

// 	for ( int i = 0; i < populationSize; i++ ){
// 		if ( !populationPtr[ i ]){
// 
// cerr << "Dead members present in current population" << endl;
// exit (1);	
// 
// // 			populationPtr.erase( populationPtr.begin() + i );
// // 			orientation.erase( orientation.begin() + i  );
// 		}
// 	}


	// NOTE: This is where any reporting should be done, as the population will have been incremented / deleted and expunged
	// check for snapshot output condition

	if ( (  currentEpoch % simulationParameters .iGet( RES_FREQ ) == 0 ) && ( simulationParameters .iGet( RES_FREQ ) > 0 ) ){

		if ( ( currentEpoch <= simulationParameters .iGet( R_END ) ) && ( currentEpoch >= simulationParameters .iGet( R_START ) ) ){

			char commonName[ FILENAME_MAX ];
			char variableName[ FILENAME_MAX ];
			char convertedInt[ 5 ];
	
			sprintf( convertedInt, "%d", currentEpoch );
	
			strcpy( commonName, resultsDir );
			strcat( commonName, "/p" );
			strcat( commonName, convertedInt ); 
		
			// go through all snapshot options
	
			if ( verbosity & O_POP_SNAPSHOT ){
	
				if ( populationSize > 0 ){			// this condition must be true
					strcpy( variableName, commonName );
					strcat( variableName, ".location" );
		
					popSnapshotPtr->open( variableName );
					
					for ( int i = 0; i < populationSize; i++ ){
						*popSnapshotPtr << orientation[ i ][ X_LOCATION ] << "\t";
						*popSnapshotPtr << orientation[ i ][ Y_LOCATION ] << "\n";
					}
	
					popSnapshotPtr->close();
				}
			}
	
			if ( verbosity & O_NUTRIENT_SNAPSHOT ){

				for ( int i = 0; i < globalEnvParameters[ NUMBER_FACTORS ]; i++ ){

					char convertedFactor[ 5 ];
	
					sprintf( convertedFactor, "_%d", i );

					strcpy( variableName, commonName );
					strcat( variableName, convertedFactor );
					strcat( variableName, ".environment" );
		
					environmentSnapshotPtr->open( variableName );
					
					environmentPtr[ i ].print( environmentSnapshotPtr );
		
					environmentSnapshotPtr->close();
				}
			}

		}
	}

	// optionally output positions to GNUPlot

	if ( populationSize > 0 ){
	
		if ( visOn  ){
			fprintf(gnuPipePtr, "plot '-' using 1:2 with dots\n");
			for ( int q = 0; q < populationSize; q++ ){
	
				fprintf( gnuPipePtr,"%f\t",orientation[ q ][ X_LOCATION ] );
				fprintf( gnuPipePtr,"%f\n",orientation[ q ][ Y_LOCATION ] );
			}
		
			fprintf( gnuPipePtr, "e\n" );	// signal end of data
			fflush( gnuPipePtr );		// flush file buffers	
		}

		outputStats();		// print stats for current tick, in to results file

		return ( true );

	}

	return ( false );	// population size has reached zero, and will terminate
}


void universe::refresh( int member, int action )			// update the virtual universe
{	

// automatic actions

	// add random (thermal) noise to the orientation
	orientation[member][ HEADING ] += (  ( 2.0 * ( 0.5 - ( rand() / ( double )  RAND_MAX ) ) ) *  ( PI / 20.0 ) );

	// consume nutrient 
	// calculate amount, subtract nutrient from environment and add energy to individual

	// subtract requested amount from environment, and add to individual
	for ( int i = 0; i < globalEnvParameters[ NUMBER_FACTORS ]; i++ ){

		populationPtr[ member ]->modifyEnergy( environmentPtr[ i ].localUpdate( orientation[ member ][ X_LOCATION ], orientation[ member ][ Y_LOCATION ], populationPtr[ member ]->getParameter( MAX_INTAKE )  ) );
	}
	
// optional actions

	switch( action ){

		case B_NONE:{		// run
		
			double nextLocX = orientation[ member ][ X_LOCATION ] + ( populationPtr[ member ]->getParameter( MAX_MOVEMENT ) ) * sin( orientation[ member ][ HEADING ] );
			double nextLocY = orientation[ member ][ Y_LOCATION ] + ( populationPtr[ member ]->getParameter( MAX_MOVEMENT ) ) * cos( orientation[ member ][ HEADING ] );

			
			if ( nextLocX  < 0.0 ){
				// set flag, reset heading
				nextLocX = 0.0;
				// randomly re--orientate between (270, 90)
				orientation[member][ HEADING ] = ( rand() / ( double )  RAND_MAX )  * PI;

				if ( nextLocY < 0 ){
					nextLocY = 0.0;
					orientation[member][ HEADING ] =   ( rand() / ( double )  RAND_MAX )  * PI / 2.0;
				} else {

					if  ( nextLocY > globalEnvParameters[ WORLD_HEIGHT ] - 1 ){
						nextLocY = globalEnvParameters[ WORLD_HEIGHT ] - 1;
 						orientation[member][ HEADING ] =   PI / 2.0 + ( ( rand() / ( double )  RAND_MAX )  * PI / 2.0 );
					}
				}
			}

			if ( nextLocX  >  ( globalEnvParameters[ WORLD_WIDTH ] - 1 ) ){
			
				nextLocX = globalEnvParameters[ WORLD_WIDTH ] - 1;

				orientation[member][ HEADING ] =  PI + (  ( rand() / ( double )  RAND_MAX )  * PI );

				if ( nextLocY < 0 ){
					nextLocY = 0.0;
					orientation[member][ HEADING ] = ( 3 * PI )  / 2.0 +  ( rand() / ( double )  RAND_MAX )  * PI / 2.0;
				} else {

					if  ( nextLocY > globalEnvParameters[ WORLD_HEIGHT ] - 1 ){
						nextLocY = globalEnvParameters[ WORLD_HEIGHT ] - 1;
 						orientation[member][ HEADING ] =  ( PI / 2.0 ) + ( ( rand() / ( double )  RAND_MAX )  * PI / 2.0 );
					}
				}
			}


//  Y co-ordinates:


			if ( nextLocY  < 0.0 ){
			// this means that nextLocX is valid wrap value between -3pi/2 and pi / 2
				nextLocY = 0.0;
				orientation[member][ HEADING ] = ( 3 * PI )  / 2.0 +  ( rand() / ( double )  RAND_MAX )  * PI;
			}else{
				if ( nextLocY > globalEnvParameters[ WORLD_HEIGHT ] - 1 ){
					nextLocY  = globalEnvParameters[ WORLD_HEIGHT ] - 1;
					orientation[member][ HEADING ] = ( 3 * PI )  / 2.0 +  ( rand() / ( double )  RAND_MAX )  * PI;
				}
			}


// NOTE: old method

// 			nextLocX =  ( nextLocX < 0.0 ) ? 0.0 : nextLocX;
// 			nextLocX =  ( nextLocX > ( globalEnvParameters[ WORLD_WIDTH ] - 1 ) ) ? globalEnvParameters[ WORLD_WIDTH ] - 1 : nextLocX;
// 			
// 			nextLocY =  ( nextLocY < 0.0 ) ? 0.0 : nextLocY;
// 			nextLocY =  ( nextLocY > globalEnvParameters[ WORLD_HEIGHT ] - 1  ) ? globalEnvParameters[ WORLD_HEIGHT ] - 1 : nextLocY;	 

			orientation[ member ][ X_LOCATION ] = nextLocX;
			orientation[ member ][ Y_LOCATION ] = nextLocY;

			
			break;
		}

		case B_DEATH:{			// Starvation

			deathRoutine( member );
			
			break;	
		}

		case B_REPRODUCTION:{			// reproduce

			// create child instance
			bacteria *childPtr = new bacteria;

			// initialise and pop on to population stack
			childPtr->init( populationPtr[ member ]);
			childPtr->setupReceptor();
			populationPtr.push_back( childPtr );

			// set coordinates and orientation
			double* valPtr;
			valPtr = new double[ 3 ];
	
			valPtr[ X_LOCATION ] = orientation[ member ][X_LOCATION ];
			valPtr[ Y_LOCATION ] = orientation[ member ][Y_LOCATION ];
			valPtr[ HEADING ] = (  rand() / (double) RAND_MAX ) * 2 * PI;
	
			orientation.push_back( valPtr );	

			break;
		}

		case B_TUMBLE:{			// tumble
			
			// change heading
			orientation[member][ HEADING ] = ( rand() / ( double ) RAND_MAX ) * 2 * PI;

			break;
		}

		default:

		cerr << "Unspecified behaviour requested ( " << action << " ). Terminating." << endl;
		exit( 1 );
	}
	
	behaviour[ action ]++;

}


void universe::deathRoutine( int memberIdx )				// remove an individual from the simulation
{

	delete populationPtr[ memberIdx ];	// delete individual and set to NULL
	populationPtr[ memberIdx ] = NULL;	

	delete [] orientation[ memberIdx ];	// delete dead individuals' orientation vector

	deadMembers.push_back( memberIdx );
	deathCount++;

}


void universe::outputStats( void )					// write per-epoch statistics to the global output file
{
	// zero necessary statistics counters

	double avgGenomeSize = 0.0;
	double avgTotGenes = 0.0;
	double avgDistinctGenes = 0.0;
 	double totalConnections = 0.0;
 	double connectionWeight = 0.0;

// 	double minGenSize = populationPtr[ 0 ]->getGenomeLength();
// 	double maxGenSize = minGenSize;

	double currGenSize;
	double currTotGenes;
	double currDistinctGenes;
	int strainCount[ numStrains ];

	for ( int i = 0; i < numStrains; i++ ) strainCount[ i ] = 0;

	// get stats from population

// 	// calculate population averages
	for ( int j = 0; j < populationSize; j++ ){
 
 		// get average length of genome
 		currGenSize = populationPtr[ j ]->getGenomeLength();
 		currTotGenes = populationPtr[ j ]->getSumGenes();
 		currDistinctGenes = populationPtr[ j ]->getDistinctGenes();
 		totalConnections += populationPtr[ j ]->getTotalConnections();
 		connectionWeight += populationPtr[ j ]->getSumConnectivity();
 			
		avgGenomeSize += currGenSize;
		avgTotGenes += currTotGenes;
		avgDistinctGenes += currDistinctGenes;
		
		// count strains
		strainCount[ populationPtr[ j ]->getStrain() -1 ]++;


	}

	totalConnections /= populationSize;
	connectionWeight /= populationSize;
	avgGenomeSize /= populationSize;
 	avgTotGenes /= populationSize;
 	avgDistinctGenes /= populationSize;

	// dump statistics from the last tick to the output file
	if ( verbosity & O_RESULTS ){
		 *resultsPtr
		<< currentEpoch << "\t"				// epoch number					(1)
		<< populationSize << "\t"			// population size at end of epoch		(2)
		<< behaviour[ B_NONE ] << "\t" 			// run						(3)
		<< behaviour[ B_TUMBLE ] << "\t"		// tumble					(4)
		<< behaviour[ B_DEATH ] << "\t"			// death					(5)
		<< behaviour[ B_REPRODUCTION ] << "\t"		// birth					(6)
		<< avgGenomeSize  << "\t"			// average genome size				(7)
		<< avgTotGenes << "\t"				// average number of genes			(8)
 		<< avgDistinctGenes << "\t"			// average number of proteins			(9)
 		<< connectionWeight << "\t"			// average aggregate connection weight		(10)
		<< totalConnections << "\t"			// complete number of GRN connections		(11)
		<< environmentPtr->getSumVolume() << "\t";		// nutrient sum					(12)

	for ( int i = 0; i < numStrains; i++ ){
		*resultsPtr
		<< strainCount[ i ] << "\t";
	}

 		*resultsPtr << "\n";
	}
}


bool universe::setOutStreams( void )					// set up the file streams for recording data
{
	bool fileError = false;

	char filePath[ FILENAME_MAX ];

	if ( verbosity > 0 ){
		// create results directory
				
		mkdir( resultsDir, S_IRWXU | S_IRWXG );
	
		// create output streams
	
		if ( verbosity & O_RESULTS ){
			resultsPtr = new ofstream;
			strcpy( filePath, resultsDir );
			strcat( filePath,"/results.csv" );

			resultsPtr->open( filePath );

			setHeaders( );	
			
		}

		// set ouput data stream for agent class
		if ( ( verbosity & O_SINGLE_IND_DATA ) || ( verbosity & O_ALL_IND_DATA )  ){
			individualDataPtr = new ofstream;
			strcpy( filePath, resultsDir );
			strcat( filePath,"/individualData.csv" );

			individualDataPtr->open( filePath );
		}

		// set ouput location stream for agent class
		if ( ( verbosity & O_SINGLE_IND_LOCATION ) || ( verbosity & O_ALL_IND_LOCATION )  ){
			individualLocationsPtr = new ofstream;
			strcpy( filePath, resultsDir );
			strcat( filePath,"/individualLocations.csv" );

			individualLocationsPtr->open( filePath );
		}


	}

	return ( fileError );

}


void universe::setHeaders( void )					// write human-readable headers in the global output file
{

	*resultsPtr << "# Epoch\tPopulation Size\tSum Nutrient\tD: AGE\tD: internal poisoning\tD: starvation\tD: external poisoning\tB:Reproduction\tB:Tumble\tB:Run\tAverage genome length\tMin genome legth\tMax genome length\tAverage Gene Number\tProteome Size\tSimple Connection\tAA Connection Weight\tSum of non-null connections\n";
	
}


void universe::removeComments( ifstream* inputFile )			// remove extraneous data from the file stream
{
	bool	comments = true;
	char	c = 0;

	while ( comments ){

		// ignore any line feeds left in the stream
		while ( inputFile->peek() == '\n' || inputFile->peek() == ' ' || inputFile->peek() == '\t' || inputFile->peek() == ',') 
			inputFile->get();	

		while ( inputFile->peek() == '#' ) inputFile->ignore( 512, '\n' );

		inputFile->get( c );

		if ( c == '\n' || c == '\t' || c == ' ' || c == '#' || c == ',' )
			comments = true;
		else{
			comments = false;
			inputFile->putback( c );
		}

	}
}


void universe::startGP( void )						// start a GNUplot thread for real-time display
{

// setup gnuplot
		gnuPipePtr = popen( "/usr/bin/gnuplot", "w" );
		fprintf( gnuPipePtr, "set title \"Population Snapshot\"\n" );
		fprintf( gnuPipePtr, "unset key\n" );
		fprintf( gnuPipePtr, "unset xtics\n" );
		fprintf( gnuPipePtr, "unset ytics\n" );
		fprintf( gnuPipePtr, "set xrange[0:%i ] \n", globalEnvParameters[ WORLD_WIDTH ] );
		fprintf( gnuPipePtr, "set yrange[0:%i ] \n", globalEnvParameters[ WORLD_HEIGHT ] );
		fprintf( gnuPipePtr, "set size square\n" );


}
