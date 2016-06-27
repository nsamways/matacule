/*

genotype.h
Header File for genotype class

(c) Neale Samways 2006-9

*/


#ifndef GENOTYPE_H
#define GENOTYPE_H

#include<vector>
#include<iostream>
#include<ostream>

#include"global.h"
#include"genome.h"
#include"operon.h"
#include"chemistry.h"
#include"protein.h"

using namespace std;

class genotype
{

public:

	genotype();					// default constructor
	~genotype();					// destructor

	// setters
	void initialise( genome* );				// determine functional genome content
	void setChemistry( chemistry* , double );	// set the pointer to the chemistry and set mutation probability
	void setupProteome( protein* );			// construct the proteome
	void resetProteomeConc( void );			// reset all proteins to initial state
	bool loadProteomeConcs( char* );		// load the protein concentrations from a file

	// others
	void cycleRN( void );					// update protein concentrations based on regulatory network
	void updateInputConc( double* );		// update the concentration of input proteins based on current internal and external state

	// printers
	void print( ofstream* );				// print out information
	void printProteinConcs( void );			// dump proteomic concentrations to file

	// reporters
	int getConnectionCount( void );			// return the number of active connections in the GRN
	double getConnectionDensity( void );		// returns summed absolute weighting of all connections
	bool checkConnectionPath( int, int ); 		// check if a path exists between two nodes of the GRN

	// Inline getters
	inline int getTotalGenes( void ){ return( totalGenes ); }								// returns the number of genes in the genome
	inline int getProteomeSize( void ){ return( proteomeSize ); }							// return distinct number of proteins
	inline double getProteinConcentration( int index ){ return( proteomePtr[ index ].getConcentration() ); }	// return protein concentration
	inline protein* getProteomePtr( void ){ return proteomePtr; }							// returns pointer to proteome	

	// inline setters
	inline static void setTotalReceptors( int nReceptors ){ totalReceptors = nReceptors; }			// set the number of env receptors
	inline void setRecording( bool toggle ){ record = toggle; }							// switch on recording
	inline void setProteomeStream( ofstream *strPtr ){ proteomeRegister = strPtr; }					// set the file stream
	inline void setProteinConcentration( int index, double concentrationDelta ){ proteomePtr[ index ].setConcentration( concentrationDelta );}	// set concentration of indicated protein
	inline void modifyProteinConcentration( int index, double concentrationDelta ){ proteomePtr[ index ].modifyConcentration( concentrationDelta );}	// modify concentration of indicated protein

private:

	// functions
	void printProteome( ofstream* );	// record current protein concentrations to output file
	void removeComments( ifstream* );	// clean input stream and tokenize

	// variables
	bool record;				// toggle for recording all proteome activity

	int 	proteomeSize,			// total number of *distinct* proteins in the proteome
		totalGenes,				// total number of genes in the genome
		totalOperons,			// total number of operons in the genome
		inputProteinNumber;		// total number of receptor proteins

	double 	chemistryMutationProb,		// probability of mutating chemistry NB. set to NULL if not offspring
		phiDeltaP[ TOTAL_PROTEINS ],	// array of changes to protein concentrations
		phiDeltaN[ TOTAL_PROTEINS ];	// array of changes to protein concentrations

	static int totalReceptors;		// number of receptors (i.e. number of environmental artifacts)

	vector< int > activeProteins;		// vector of active elements of the proteome

	ofstream *proteomeRegister;		// pointer to the file stream for proteome recording

	operon *operonPtr;			// array of operons
	protein *proteomePtr;			// the proteome
	chemistry *chemIntPtr;			// pointer to individuals' chemistry object

};

#endif
