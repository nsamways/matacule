/*

genome.h
Header File for genome class

(c) Neale Samways 2006-9

*/

#ifndef GENOME_H
#define GENOME_H

#include<fstream>

using namespace std;

class genome
{

public:

	genome();					// default constructor
	~genome();					// destructor
	
	// setters
	void initialise( int );				// initialiser for random genomes
	bool initialise( char* );			// initialiser for pre-specified genome
	bool initialise( genome*, double, int );	// initialiser for child genome

	// printers
	void print( ostream* );				// write genome to file

	// inline setters
	inline void setLength( int length ){ currentLength = length; }		// change the length of the genome

	// inline getters
	inline int getLength( void ){ return currentLength; }			// return length of the genome
	inline int* getStart( void ){ return sequencePtr; }			// return the first element of the DNA
	inline int getBase( int locus ){ return sequencePtr[ locus ]; }		// return the base at arbitraty position

private:

	int *sequencePtr;	// pointer to integer array of values (genome sequence)
	int currentLength;	// length of genome

};

#endif
