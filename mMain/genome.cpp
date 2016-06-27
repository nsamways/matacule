/*

genome.cpp

member functions and data members for genome class

(c) Neale Samways 2006-9

*/ 


#include<iostream>
#include<fstream>
#include<cstring>

#include"genotype.h"
#include"global.h"

using namespace std;

genome::genome()									// constructor
{
	currentLength = 0;
	sequencePtr = NULL;	
}


genome::~genome()									// destructor
{
	delete [] sequencePtr;
	sequencePtr = NULL;
}


void genome::initialise( int baseLength )						// initialise with random values
{

	currentLength = baseLength;
	sequencePtr = new int[ baseLength ];
	
	// construct random sequence
	for ( int i = 0; i < baseLength; i++ )
		sequencePtr[ i ] = rand() % 16;

}


bool genome::initialise( char* basePtr )						// initialise from pre-determined sequence: load genome from file
{
	int locus = 0;

	vector< int > loadedGenome;

	char filePathPtr[ FILENAME_MAX ];

	ifstream inStream;


	strcpy( filePathPtr, basePtr );
	strcat( filePathPtr, "/genome.csv" );

// load genome from file
	
	// check the the file can be opened 	
	inStream.open( filePathPtr );

	if ( inStream.fail() ){
		cerr << "Unable to open genome file: " << filePathPtr << endl;
		return ( false );
	}

	// read in the genome
	while( !inStream.eof() ){

		// remove any comments etc
	
		bool	comments = true;
		char	c = 0;
		
		while ( comments ){
	
			// ignore any line feeds left in the stream
			while ( inStream.peek() == '\n' || inStream.peek() == ' ' || inStream.peek() == '\t' || inStream.peek() == ',') 
				inStream.get();	
	
			inStream.get( c );
	
			if ( c == '\n' || c == '\t' || c == ' ' || c == '#' || c == ',' )
				comments = true;
			else{
				comments = false;
				inStream.putback( c );
			}
	
		}
			int base;
			inStream >> base;
			loadedGenome.push_back( base );
			locus++;
	}

	inStream.close();

	currentLength = loadedGenome.size();
 
	sequencePtr = new int[ currentLength ];
	
	for ( int i = 0; i < currentLength; i++ )
		sequencePtr[ i ] = loadedGenome[ i ];

	loadedGenome.clear();

	return( true );
}


bool genome::initialise( genome* parentPtr, double SNPrate, int changeSize )		// create offspring genome
{

	int spliceIndex = -1;

	if ( changeSize > 0 ){
		if ( changeSize >= parentPtr->currentLength )
			changeSize = parentPtr->currentLength -1;

		spliceIndex = changeSize + ( rand()%( parentPtr->currentLength - changeSize ) );
	}

	if ( changeSize < 0 )
		spliceIndex = rand()% ( parentPtr->currentLength + changeSize );

	// set the current length for the object (dont delete this statement!)
	currentLength = parentPtr->currentLength + changeSize;

	sequencePtr = new int[ currentLength ];

	int j = 0;
	bool duplicated = false;

	for ( int i = 0; i < parentPtr->currentLength; i++ ){
		
		if ( ( ( double)rand() / RAND_MAX ) < SNPrate ){

			// mutate
			sequencePtr[ j ] = rand()%16;
		}else{

			sequencePtr[ j ] = parentPtr->sequencePtr[ i ];
			}
	
		j++;

		if ( ( i == ( spliceIndex ) ) && ( duplicated == false ) ){
			duplicated = true;
			i -= changeSize;
		}
	}
	
	return( true );
}


void genome::print( ostream *outputStreamPtr )						// write genome to file
{

	*outputStreamPtr << "GL:\t" << currentLength << "\n";
	*outputStreamPtr << "GS:\t";
	for ( int i = 0; i < ( currentLength - 1 ); i++ )
		*outputStreamPtr << sequencePtr[ i ] << ",";
	
	*outputStreamPtr << sequencePtr[ currentLength - 1 ] << endl;	
}
