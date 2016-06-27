/*

operon.cpp

member functions and data members for operon class

(c) Neale Samways 2006-9

*/ 

#include<iostream>
#include<iomanip>
#include<fstream>

#include"operon.h"

using namespace std;

operon::operon()						// constructor
{
	
}


operon::~operon()						// destructor
{

}


void operon::init( int locus, genome *currGenome )		// initialise operon by searching genome
{
	int shift = PROMOTOR_LENGTH + CIS_REGION_LENGTH; 		// makes it less cumbersome for functional region calculation

// 	// copy the promotor signature	NOTE: This is unlikely to be needed
// 	for ( int i = 0; i < PROMOTOR_LENGTH; i++ )
// 		promotor[ i ] = currGenome->getBase( locus +i );

	for ( int i = 0; i < INIT_GENE_NUM; i++ ) ptn[ i ] = currGenome->getBase( locus + shift + i );

}


void operon::transcribe( protein* proteomePtr )			// increase concentration of proteins based on genes in operon
{
	// nb. takes the proteome and updates accordingly
	for ( int i = 0; i < INIT_GENE_NUM; i++ )				// go through genes in order ( one entry per gene )
			proteomePtr[ ptn[ i ] ].modifyConcentration( BASAL_INCREMENT_LEVEL );

	// All genes in this operon will now have been transcribed; only present genes are transcribed.
}


void operon::print( ofstream *outputFilePtr )			// write encoded proteins to file
{
	*outputFilePtr << "\nGenes: ";

	for ( int i = 0; i < INIT_GENE_NUM; i++ ){
		*outputFilePtr << ptn[ i ];
		if ( i != INIT_GENE_NUM -1){ *outputFilePtr << " , "; 
			}else{
			*outputFilePtr << endl; }
	}
}
