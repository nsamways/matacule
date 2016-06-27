/*

parameter.cpp

member functions and data members for parameter class

(c) Neale Samways 2007-9

*/ 

#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdlib>

#include"parameter.h"

using namespace std;

parameter::parameter()							// constructor
{

// set up the integer,double amd char parameter arrays

	dParameterPtr = NULL;
	iParameterPtr = NULL;
	cParameterPtr = NULL;

}


parameter::~parameter()							// destructor
{

		delete [] dParameterPtr;
		dParameterPtr = NULL;

		delete [] iParameterPtr;
		iParameterPtr = NULL;

		for ( int i = 0; i < parameterCount; i++ )
			delete [] cParameterPtr[ i ];

		delete [] cParameterPtr;
		cParameterPtr = NULL;
}


bool parameter::initialise( int pTot, char *fileLocation )		// set up parameter object based on number of entries
{
	bool success = false;

	// set number of parameters
	parameterCount = pTot;

	// set up arrays
	dParameterPtr = new double[ pTot ];
	iParameterPtr = new int[ pTot ];
	cParameterPtr = new char*[ pTot ];

	for ( int i = 0; i < pTot; i++ )
		cParameterPtr[ i ] = new char[ FILENAME_MAX ];
	
	// load in parameters
	
	success = loadParameters( fileLocation );	

	return( success );
}


int parameter::iGet( int paramNum )					// return requested integer casted parameter
{

	if ( paramNum < parameterCount ){
		return ( iParameterPtr[ paramNum ] );
	}else{
		cerr << "Request for undefined parameter (integer) : " << paramNum << " ( of: " << parameterCount << " ) Exiting." << endl;
		exit( 1 );
	}

}


double parameter::dGet( int paramNum )					// return requested double casted parameter
{

	if ( paramNum < parameterCount ){
		return ( dParameterPtr[ paramNum ] );
	}else{
		cerr << "Request for undefined parameter: (double): " << paramNum << " ( of: " << parameterCount << " ) Exiting." << endl;
		exit( 1 );
	}
}


char* parameter::cGet( int paramNum )					// return requested char casted parameter
{

	if ( paramNum < parameterCount ){
		return ( cParameterPtr[ paramNum ] );
	}else{
		cerr << "Request for undefined parameter: (char): " << paramNum << " ( of: " << parameterCount << " ) Exiting." << endl;
		exit( 1 );
	}
}


void parameter::print( void )						// write parameters to file
{
	// print parameters

	cout << "integer:\t\tDouble\n";

	for ( int i = 0; i < parameterCount; i ++ ){
		cout << iParameterPtr[ i ] << "\t\t" << dParameterPtr[ i ] << "\n";
	}
	
	cout << endl;

}


bool parameter::loadParameters( char *fileNamePtr )			// load parameters from file
{
	char rawParameter[ parameterCount ][ FILENAME_MAX ];

	int currentParam = 0;
	ifstream inStream;

	// read in the parameter file

	inStream.open( fileNamePtr );

	// check the the file can be opened 	
	if ( inStream.fail() ){
		cerr << "Unable to open parameter file: " << fileNamePtr << endl;
		return ( false );
	}

	// read in the parameters
	while( !inStream.eof() && ( currentParam < parameterCount ) ){
		removeComments( &inStream );
		inStream >> rawParameter[ currentParam ];
		currentParam++;

	}

	inStream.close();

	// cast the parameter to the correct type and copy to array
	for ( int i = 0; i < parameterCount; i++ ){
		// integer array:
		iParameterPtr[ i ] = atoi(  rawParameter[ i ] );

		// double array;
		dParameterPtr[ i ] = atof( rawParameter[ i ] );

		// char array
		strcpy( cParameterPtr[ i ], rawParameter[ i ]);

	}

	return ( true );

}


void parameter::removeComments( ifstream* inputFile )			// remove extraneous characters from input stream
{
	bool	comments = true;
	char	c = 0;

	while ( comments ){

		// ignore any line feeds left in the stream
		while ( inputFile->peek() == '\n' || inputFile->peek() == ' ' || inputFile->peek() == '\t' ) 
			inputFile->get();	

		while ( inputFile->peek() == '#' ) inputFile->ignore( 512, '\n' );

		inputFile->get( c );

		if ( c == '\n' || c == '\t' || c == ' ' || c == '#' )
			comments = true;
		else{
			comments = false;
			inputFile->putback( c );
		}

	}
}
