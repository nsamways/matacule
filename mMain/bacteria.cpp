/*

bacteria_C.cpp

member functions and data members for bacteria, derrived class
RN update routine determined by receptorType variable

(c) Neale Samways 2006-9

*/ 

#include<iostream>
#include<fstream>

#include"bacteria.h"


bacteria::bacteria()									// constructor
{

// create vector of prior environmental factors
	priorEnvironmentPtr = new double [ totalEnvironmentFactors ];

	for ( int i = 0; i < totalEnvironmentFactors; i++ ) priorEnvironmentPtr[ i ] = 0.0;

}

bacteria::~bacteria()									// destructor
{

	delete [] priorEnvironmentPtr;

}


void bacteria::internalUpdateP( double *environmentPtr )		// Cycle GRN / protein network etc.
{

	// NOTE: update routine for PTD individuals

	// determine factor deltas and pass to genotype
	for ( int i = 0; i < totalEnvironmentFactors; i++ ){

		genotypePtr->modifyProteinConcentration( i, alpha * ( environmentPtr[ i ] - priorEnvironmentPtr[ i ] ) );

		priorEnvironmentPtr[ i ] = environmentPtr[ i ];
	}

	// cycle regulatory network
	genotypePtr->cycleRN( );

	// dump new proteomic concentrations if necessary
	if ( individualVerbosity & O_ALL_IND_PROTEOME ) genotypePtr->printProteinConcs();
}

void bacteria::internalUpdateI( double *environmentPtr )		// Cycle GRN / protein network etc.  (new routine)
{
	//  NOTE: update routine for IV individuals

	// pass input values to network
	genotypePtr->updateInputConc(  environmentPtr  );

	// loop for specified number of epochs
	for ( int i = 0;  i < 1 ; i++ ){

		// cycle regulatory network now that inputs have been calculated 
		genotypePtr->cycleRN( );

	}	// end update loops

	// dump new proteomic concentrations if necessary
	if ( individualVerbosity & O_ALL_IND_PROTEOME ) genotypePtr->printProteinConcs();

}

void bacteria::internalUpdate( double *environmentPtr )
{

	(this->*internalUpdatePtr)( environmentPtr );

}

void bacteria::setupReceptor( void )						// set interalUpdatePtr function
{

	// set the update pointer
	if ( individualParametersPtr[ RECEPTOR_TYPE ] == 0 ){
		internalUpdatePtr = &bacteria::internalUpdateI;	// IV model
	}else{
		internalUpdatePtr = &bacteria::internalUpdateP;	//PTD model
	}

	// set alpha
	alpha = individualParametersPtr[ GAIN_PARAMETER ];

}

int bacteria::updateBehaviour( void )						// constructor
{

	int behaviour = 0;

	double dice = QUICKRAND;

	// expend energy
	modifyEnergy (- individualParametersPtr[ BASAL_ENERGY_EXPEND ]);

	// check for death or reproduction
	if ( energy <= 0.0 ){
		return( B_DEATH );			// died of starvation
	}else{
		if ( energy >= MAX_ENERGY ) return( B_REPRODUCTION );
	}

	// Reproduction or death have not occurred, thus check for protein modulated behaviours
	if ( ( genotypePtr->getProteinConcentration( B_TUMBLE ) ) > dice ) behaviour = B_TUMBLE;

	return( behaviour );
}
