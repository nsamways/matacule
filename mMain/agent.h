/*

agent.h
Header File for agent class

(c) Neale Samways 2006-9

*/

#ifndef AGENT_H
#define AGENT_H

#include<iostream>
#include<ostream>
#include<iomanip>

#include"global.h"
#include"chemistry.h"
#include"parameter.h"
#include"genotype.h"
#include"genome.h"

using namespace std;

class agent
{

public:
	agent();						// constructor
	virtual ~agent();					// destructor

	// Initialisers	
	bool init( int );					// initialiser for initial population genomes
	void init( agent* );					// initialiser for offspring

	// Updaters
	virtual void internalUpdate( double* );			// Cycle GRN / protein network etc.
	virtual int updateBehaviour( void );			// per epoch behaviour routine

	//printers
	void printConfig( ostream* oConfPtr );			// print individual parameters to file
	void outputBaseData( void );				// write all individual data to file

	// Inline Getters
	inline int getUniqueID( void ){ return uniqueID; }									// return UID
	inline int getStrain( void ){ return strain; }										// return strain
	inline int getAge( void ){ return age; }										// return age
	inline double getEnergy( void ){ return ( energy );}									// return energy
	inline int getGenomeLength( void ){ return( dnaPtr->getLength() ); }							// return genome length	
	inline int getSumGenes( void ){ return( genotypePtr->getTotalGenes() );}						// return sum genes in genome
	inline int getDistinctGenes( void ){ return( genotypePtr->getProteomeSize() ); }					// return distinct genes in genome
	inline double getParameter( int index ){ return individualParametersPtr[ index ]; } 					// return selected parameter
	inline int getTotalConnections( void ){ return genotypePtr->getConnectionCount(); }					// return sum of RN connections
	inline double getSumConnectivity( void ){ return genotypePtr->getConnectionDensity(); }					// return sum of non-zero connections
	inline double getProteomicConcentration( int index ){ return genotypePtr->getProteinConcentration( index ); }		// return protein concentration	
	inline protein* getProteomePointer( void ){ return genotypePtr->getProteomePtr(); }					// return pointer to proteome	

	// Inline Setters
	inline void setIndividualVerbosity( int verbosityOverride ){ individualVerbosity = verbosityOverride; }			// override individual verbosity
	inline void modifyEnergy( double eDelta ){ energy += eDelta; }								// increment / decrement energy
	inline static void setWorldEpoch( int wEpoc ){ agent::worldEpoch = wEpoc; }						// set class world epoch
	inline static void setOutputPtr( ofstream* oPtr ){ agent::dataOutPtr = oPtr; }						// set class pointer to output file
	inline static void setBaseOutput( char *baseName ){ strcpy( agent::baseOutputPath, baseName ); }			// set class base filename
	inline static void setOutputVerbosity( int verbosity ) { agent::recordLevel = verbosity; }				// set class verbosity
	inline static void setParameterPath( char *baseName ){ strcpy( agent::baseParameterPath, baseName ); }			// set class base parameter file 
	inline static void setEnvironmentFactors( int numEnvFactors ){ agent::totalEnvironmentFactors = numEnvFactors; genotype::setTotalReceptors( numEnvFactors );}	// set num environmental artifact 
	inline void resetProteome( void ){ genotypePtr->resetProteomeConc(); }							// set all protein values to initial

	// Inline others


protected:

	// variables
	genome *dnaPtr;					// pointer to genome
	genotype *genotypePtr;				// pointer to genotype	
	chemistry *chemistryPtr;			// pointer to individuals' chemistry interaction matrix
	parameter *paramsPtr;				// pointer to the initial indivduals' (common) parameters


	double *individualParametersPtr;		// pointer to the array of double-casted parameters
	double energy;					// internal energy level

	int individualVerbosity;			// effective verbosity

	static int totalEnvironmentFactors;		// number of environmental factors (nutrients etc.)



private:

	// functions
	void proteomeRecordInit( void );		// set up stream for proteome recording
	bool loadParameters( char* );			// load individual parameters into parameter array

	// variables


	int uniqueID,					// unique identifier
		strain,					// simple tag, inherited from parent for data analysis
		receptorType,			// simple tag identifying PTD or IV model
		creationEpoch,				// world epoch at birth
		destructionEpoch,			// world epoch ay death
		age,						// number of live epochs
		parentID,					// for calculating lineage
		numOffspring;				// Number of offspring spawned from individual

	double alpha;					// Multiplicatory constant for calculating environmental artifact sensing

	ofstream* proteomeOutPtr;			// pointer to output stream for proteome ( not always used )

	char individualOutputPath[ FILENAME_MAX ];	// individual level output path 

	// static (class level) variables
	static int iCount,				// cumulative population count
		recordLevel,				// class level verbosity
		worldEpoch;				// current epoch

	static ofstream* dataOutPtr;			// pointer to output file stream	( not always used )

	static char baseParameterPath[ FILENAME_MAX ],	// abs route to parameter file 		
		baseOutputPath[ FILENAME_MAX ];		// abs route to output file 		( not always used )

	
};

#endif
