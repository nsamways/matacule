 /*

agent.cpp

member functions and data members for agent class

(c) Neale Samways 2006-9

*/ 

// test to see if this is committed

#include<iostream>
#include<fstream>

#include"agent.h"

using namespace std;

	char agent::baseParameterPath[ FILENAME_MAX ];
	char agent::baseOutputPath[ FILENAME_MAX ];

	int agent::iCount 			=	0;
	int agent::worldEpoch 			=	0;
	int agent::recordLevel			=	0;
	int agent::totalEnvironmentFactors 	=	0;

	ofstream* agent::dataOutPtr		= 	NULL;
	

agent::agent()									// constructor
{

	// create all new abstract objects
	dnaPtr		=	new genome;		// genome
	genotypePtr 	=	new genotype;		// functional aspects of proteins
	chemistryPtr 	=	new chemistry;		// pointer to chemistry
	proteomeOutPtr 	=	NULL;			// pointer to output for proteome recording

	// initialise integer values	
	strain			=	0;
	age 			=	0;
	numOffspring		=	0;
	creationEpoch		=	0;
	destructionEpoch	=	0;
	parentID			=	0;

	energy 			=	0.0;

	receptorType		=	-1;

	alpha			=	0;

	setIndividualVerbosity( agent::recordLevel );

	uniqueID = iCount++;

	paramsPtr = NULL;

	// doubleParameters array
	individualParametersPtr = new double[ NUMBER_INDIVIDUAL_PARAMS ];
	for ( int i = 0; i < NUMBER_INDIVIDUAL_PARAMS; i++ ) individualParametersPtr[ i ] = 0.0;

}


agent::~agent()									// destructor
{
	// if verbosity set, then output data to file	
	// notify of destruction epoch (in case simulation terminates whilst individuals still alive)
	destructionEpoch = worldEpoch;

	if ( individualVerbosity & O_ALL_IND_DATA )
		outputBaseData();


	// delete all created objects
	
	delete dnaPtr;
	dnaPtr = NULL;

	delete genotypePtr;
	genotypePtr = NULL;

	delete chemistryPtr;
	chemistryPtr = NULL;

	// close proteome recording stream if not null

	if ( proteomeOutPtr != NULL ){
		proteomeOutPtr->close();
		delete proteomeOutPtr;
	}

	proteomeOutPtr = NULL;

	delete paramsPtr;
	paramsPtr = NULL;

	delete [] individualParametersPtr;
	individualParametersPtr = NULL;


}


bool agent::init( int sType )							// initialiser for initial population members
{
	strain = sType;	
	bool loadProteomeConcs = false;

	// determine base path to strain files
	char strainPathPtr[ FILENAME_MAX ];
	char individualBasePathPtr[ FILENAME_MAX ];

	char convertedInt[ 5 ];

	sprintf( convertedInt, "%d", strain );

	strcpy( strainPathPtr, baseParameterPath );
	strcat( strainPathPtr,"/strain" );
	strcat( strainPathPtr, convertedInt );

	// read individual parameters
	strcpy( individualBasePathPtr, strainPathPtr );
	strcat( individualBasePathPtr, "/individual" );

	if ( !loadParameters( individualBasePathPtr ) ){
		cerr << "Unable to load parameters for individual: " << uniqueID << " , strain " << strain << endl;
		return( false );
	}

	// set up chemistry & check initialisation conditions

	// set the random deletion probability
	chemistry::setInteractionDeletionProb( paramsPtr->dGet( INTERACTION_DEL_PROB ) );

	if ( paramsPtr->dGet( CHEMISTRY_SETUP ) > 0.0 ){
		chemistryPtr->randomizeInteractions( paramsPtr->dGet( CHEMISTRY_SETUP ) );	// randomize
	}else{
		
		if ( paramsPtr->dGet( CHEMISTRY_SETUP ) < 0 ){
			if ( !chemistryPtr->load( strainPathPtr ) ){				// load from file
				cerr << "Unable to load chemistry for individual: " << uniqueID << " , strain " << strain << endl;
				return( false );
			}
		}
	}	// default to zero'd interaction matrix

	// set up genome
	if ( paramsPtr->iGet( GENOME_INITIALISATION ) < 0 ){ 
	// signal to load proteomic Concs
	loadProteomeConcs = true;

	// load genome from file
		if ( !dnaPtr->initialise( strainPathPtr ) ){
			cerr << "Unable to initialise genome from file: " << strainPathPtr << endl;
			return( false );
		}

	}else{
	// create random genome
		dnaPtr->initialise( paramsPtr->iGet( GENOME_INITIALISATION ) ); 
	}

	// set up recordng if activated
	if ( individualVerbosity & O_ALL_IND_PROTEOME )	proteomeRecordInit();

	// Genome constructed; setup Genotype
	genotypePtr->initialise( dnaPtr );
	genotypePtr->setupProteome( NULL );				// set up the proteome
	genotypePtr->setChemistry( chemistryPtr, 0.0 );		// send pointer to chemistry

	if ( loadProteomeConcs ){ 								// load preset values if possible. NOTE: this must be done *after* setting up the proteome
		if ( !genotypePtr->loadProteomeConcs( strainPathPtr ) ) cout << " ! Using default values !" << endl;
	}
			    
	// set final variables
	creationEpoch = worldEpoch;

	energy = ( INITIAL_ENERGY);	

	return ( true );	// The routine was able to initialise the DNA
}


void agent::init( agent* parentPtr )						// overwritten initialiser for offspring
{
	// initialise an individual based on *parentPtr ( with probabilistic mutation )
	int dupLength = 0;
	int maxDeletion = 0;
	int maxDuplication = 0;
	int parDNALength = parentPtr->getGenomeLength();

	double dice = 0.0;

// update parent parameters

	// increase offspring count
	parentPtr->numOffspring++;

	// subtract energy for reproduction and divide between parent and offspring
	parentPtr->modifyEnergy( parentPtr->energy / -2.0 );
	energy = parentPtr->energy;

// initialise new offspring

	//copy parent parameters
	for ( int i = 0; i < NUMBER_INDIVIDUAL_PARAMS; i++ ) individualParametersPtr[ i ] = parentPtr->individualParametersPtr[ i ];

	creationEpoch = worldEpoch;

	// set type and lineage
	strain = parentPtr->strain;
	parentID = parentPtr->uniqueID;

	// set up chemistry based on parent and pass pointer and mutation probability to genotype
	*chemistryPtr = parentPtr->chemistryPtr;
	genotypePtr->setChemistry( chemistryPtr, individualParametersPtr[ CHEMISTRY_MUT_PROB ] );

	// set up genome
	// set max mutation change NOTE: a value of -1 indicates entire genomic reproduction

	// Evolution routine -- determine mutations
	dice = rand() / ( double )( RAND_MAX );

	// determine genomic deletion
	if ( ( dice < individualParametersPtr[ DELETION_PROB ] ) && parDNALength > 11 ){

		// calculate largest allowable deletion
		if ( individualParametersPtr[ MAX_MUT_CHANGE ] < 0 ){		// complete deletion possible
			maxDeletion = ( parDNALength - 10 );
		} else {
			maxDeletion = min( parDNALength - 10, ( int ) individualParametersPtr[ MAX_MUT_CHANGE ] );
		}
		// nb, this is uniform
		dupLength += ( ( 1 + rand()% maxDeletion ) * -1 );

	}

	dice = rand() / ( double ) RAND_MAX;

	if ( dice < individualParametersPtr[ DUPLICATION_PROB ] ) {
	// calculate largest allowable duplication
	
		if ( individualParametersPtr[ MAX_MUT_CHANGE ] < 0 ){
			maxDuplication = parDNALength;
		} else {
			maxDuplication = min( parDNALength, ( int ) individualParametersPtr[ MAX_MUT_CHANGE ] );
		}

		dupLength += ( 1 + ( rand()% maxDuplication ));
	}

	// initialise the DNA
	dnaPtr->initialise( parentPtr->dnaPtr, individualParametersPtr[ SNP_RATE ], dupLength );

	// initialise genotype
	genotypePtr->initialise( dnaPtr );

	// set up recordng if activated
	if ( individualVerbosity & O_ALL_IND_PROTEOME )	proteomeRecordInit();

	// set up the proteome, by passing parent
	genotypePtr->setupProteome( parentPtr->getProteomePointer() );	


	// NOTE: This is where individual parameters could be evolved in future implimentations

}


void agent::proteomeRecordInit( void )						// open outputfile in the appropriate location
{
 
	strcpy( individualOutputPath, baseOutputPath );
	strcat( individualOutputPath, "/member" );

	char convertedInt[ 5 ];
	sprintf( convertedInt, "%i", uniqueID );
		
	strcat( individualOutputPath, convertedInt );
	strcat( individualOutputPath, ".proteome" );

	proteomeOutPtr = new ofstream;

	proteomeOutPtr->open( individualOutputPath );

	genotypePtr->setProteomeStream( proteomeOutPtr );

	genotypePtr->setRecording( true );

	
}


void agent::internalUpdate( double *environmentPtr )				// Cycle GRN / protein network etc.
{

	// update receptor activity, cycle regulatory network

	for ( int i = 0; i < 1; i++ ){

		// update proteome according to environmental factors
		for ( int j = 0; j < totalEnvironmentFactors; j++ )
			genotypePtr->setProteinConcentration( j , environmentPtr[ j ]  );

		// cycle regulatory network
		genotypePtr->cycleRN( );
	}

	// dump new proteomic concentrations if necessary
	if ( individualVerbosity & O_ALL_IND_PROTEOME ) genotypePtr->printProteinConcs();

}


int agent::updateBehaviour( void )				// determine behaviour for current epoch
{

	int behaviour = 0;


	return( behaviour );


}


void agent::printConfig( ostream* oConfPtr )					// print individual parameters to file
{

	*oConfPtr << "\nIndividual: " << uniqueID << "\n";

	// loop through all parameters
	
	for ( int i = 0; i < NUMBER_INDIVIDUAL_PARAMS ; i++ ){
	
	*oConfPtr << paramsPtr->dGet( i ) << "\n";

	}

	*oConfPtr << endl;


}


void agent::outputBaseData( void )						// write all individual data to file
{

	*dataOutPtr << ":::BEGIN MEMBER INFO:::"
		<< "\nPID:\t" 
		<< parentID 
		<< "\nUID:\t"
		<< uniqueID
		<< "\nSTRAIN:\t"
		<< strain
		<< "\nCE:\t"
		<< creationEpoch
		<< "\nDE:\t"
		<< destructionEpoch
		<< "\nAGE:\t"
		<< age
		<< "\nOS:\t" << numOffspring << "\n";

 	dnaPtr->print( dataOutPtr ); 		// print DNA sequence to file
	genotypePtr->print( dataOutPtr );	// print proteome information to file
	chemistryPtr->print( dataOutPtr );	// print chemistry to file

	*dataOutPtr << ":::END MEMBER INFO:::" << endl;

}


bool agent::loadParameters( char *bPtr )					// load parameters from file - used for initial population
{


	bool success = false;

	char filePathPtr[ FILENAME_MAX ];				// path to parameter file

	strcpy( filePathPtr, bPtr );
	strcat( filePathPtr,".param" );


	paramsPtr = new parameter;
	
	success = paramsPtr->initialise( NUMBER_INDIVIDUAL_PARAMS, filePathPtr );

	// copy parameters to data members - treat all as doubles	
	// get parameters in to local array

	for ( int i = 0; i < NUMBER_INDIVIDUAL_PARAMS; i ++ )
		individualParametersPtr[ i ] = paramsPtr->dGet( i );

	return( success );

}
