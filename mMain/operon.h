/*

operon.h
Header File for operon class

(c) Neale Samways 2006-9

*/


#ifndef OPERON_H
#define OPERON_H

#include"genome.h"
#include"protein.h"

using namespace std;


class operon
{

public:

	operon();				// constructor
	~operon();				// destructor
	void init( int, genome* );		// initialiser for pre-specified genome
	void transcribe( protein* );		// express all genes within operon by basal amount

	void print( ofstream* );		// print operon to file

	// Inlines
	int getProteinIdent( int index ){ return ptn[ index ]; }	// return identitiy of requested gene
	int getNumGenes( void ){ return( INIT_GENE_NUM ); }		// return number of genes in operon
private:
	
// 	int cisRegion[ CIS_REGION_LENGTH ];
	int ptn[ INIT_GENE_NUM ];				// proteins encoded by gene i.e. when transcribed, these proteins get more concentrated 

};

#endif
