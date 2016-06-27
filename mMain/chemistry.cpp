/*

chemistry.cpp

member functions and data members for chemistry class

(c) Neale Samways 2006-9

*/ 

#include<iostream>
#include<cmath>

#include"chemistry.h"

using namespace std;

	double chemistry::interactionDeletionProb = 0.0;

chemistry::chemistry()									// constructor
{

// set up the reactions matrix to no interactions

	reactions = new double*[ TOTAL_PROTEINS ];

	for ( int i = 0; i < TOTAL_PROTEINS; i++ ){
		reactions[ i ] = new double[ TOTAL_PROTEINS ];
		
		for ( int j = 0; j < TOTAL_PROTEINS; j++ )
			reactions[ i ][ j ] = 0.0;
	}

}


chemistry::~chemistry()									// destructor
{

// destruct the matrix
	for ( int i = 0; i < TOTAL_PROTEINS; i++ )
		delete[] reactions[ i ];
	
	delete [] reactions;
	
	reactions = NULL;	
	
}


void chemistry::operator= ( chemistry *initialChem )					// copy constructor
{

//	copy the chemical reactions matrix

	for ( int i = 0; i < TOTAL_PROTEINS ; i++ ){
		for ( int j = 0; j < TOTAL_PROTEINS; j++ ){

			this->reactions[ i ][ j ] = initialChem->reactions[ i ][ j ];
		}
	}
}


void chemistry::mutate( int product, int reactant, double mProb )			// add Gaussian noise to specified interaction
{

// apply mutation only to active proteins
// get passed only active individual reactions and mutation probability
	
// gaussian mutation

			if ( ( ( double) rand() / RAND_MAX ) < mProb ){
				// create random number with Gaussian probability

				reactions[ product ][ reactant ] += generateGaussian( 0.0, 0.2 );
				
				// apply bounds by clamping
				reactions[ product ][ reactant ] = ( reactions[ product ][ reactant ] <= 1.0 ) ? reactions[ product ][ reactant ] : 1.0;
				reactions[ product ][ reactant ] = ( reactions[ product ][ reactant ] >= -1.0 ) ? reactions[ product ][ reactant ] : -1.0;
			}

			// stochastically delete interaction
			if ( ( ( double) rand() / RAND_MAX ) < interactionDeletionProb ) reactions[ product ][ reactant ] = 0.0;

}


void chemistry::print( ofstream *outStreamPtr )						// Write matrix to file stream	
{
	// display the interaction matrix

	*outStreamPtr << "\n#Chemistry Interaction Matrix\n";

	for ( int i = 0; i < TOTAL_PROTEINS; i++ ){

		for ( int j = 0; j < TOTAL_PROTEINS; j++ ) *outStreamPtr << reactions[ i ][ j ] << ",";

		*outStreamPtr << endl;
	}
	
	*outStreamPtr << "#End Chemistry\n" << endl;

}


void chemistry::randomizeInteractions( double threshold )				// Randomize the interaction matrix
{

	double index;

	for ( int i = 0; i < TOTAL_PROTEINS - 1 ; i++ ){

		for ( int j = 0; j < TOTAL_PROTEINS ; j++ ){

			if ( ( ( double ) rand() / RAND_MAX ) < threshold ){
				// give a random interaction
				index = ( ( ( double ) rand() / RAND_MAX ) * 2.0 ) - 1.0;
			}else{
				index = 0.0;
			}	

			reactions[ i ][ j ] = index;
		}
	}

}


bool chemistry::load( char *basePathPtr )						// load a pre-specified chemistry
{
	// load chemistry interaction matrix from file

	double index = -1.0;
	ifstream inStream;

	char filePathPtr[ FILENAME_MAX ];	// absolute name of strain directory

	strcpy( filePathPtr, basePathPtr );
	strcat( filePathPtr,"/chemistry.csv" );

	// open the chemistry file
	
	inStream.open( filePathPtr );
	
	// check the the file can be opened 	
	if ( inStream.fail() ){
		cerr << "Unable to open chemistry file:" << filePathPtr << endl;
		return( false );
	}

	for ( int i = 0; i < TOTAL_PROTEINS; i++ ){

		for ( int j = 0; j < TOTAL_PROTEINS; j++ ){

			removeComments( &inStream);
			inStream >> index;
			setInteraction( i, j, index );
		}
	}
	
	inStream.close();

	return( true );

}


void chemistry::removeComments( ifstream* inputFile )					// remove extraneous data from the stream
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


double chemistry::generateGaussian( double m, double s )				// generate a Gaussian number
{
	// normal distribution with mean m and standard deviation s
	static double normal_x2;
	static bool normal_x2_valid;

	double normal_x1;

	double w;                           // radius
	if ( normal_x2_valid ) {              // we have a valid result from last call
		normal_x2_valid = false;
		return( normal_x2 * s + m );
	}

	// make two normally distributed variates by Box-Muller transformation
	do {
		double rand1 = ( double ) rand()/RAND_MAX;
		double rand2 = ( double ) rand()/RAND_MAX;
		normal_x1 = 2. * rand1 - 1.;
		normal_x2 = 2. * rand2 - 1.;
		w = normal_x1*normal_x1 + normal_x2*normal_x2;
	} while (w >= 1. || w < 1E-30);

	w = sqrt( log( w ) * ( -2.0 / w ) );
	normal_x1 *= w;  normal_x2 *= w;    // normal_x1 and normal_x2 are independent normally distributed variates
	normal_x2_valid = true;                // save normal_x2 for next call


	return( normal_x1 * s + m );
}
