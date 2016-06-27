// test instance for bacteria simulation


#include"universe.h"
#include<iostream>
#include<iomanip>
#include<cstring>

using namespace::std;

void usage( void );
int convertSeed( char* );
char* strrev( char* );

int main( int argc, char *argv[] )
{
	int arg = 1;
	int verbosity = 0;
	int completedEpochs = 0;
	int sOffset = 0;

	bool veryVerbose = false;
	bool visualisation = false;
	char resultsDir[ FILENAME_MAX ];	
	char parametersDir[ FILENAME_MAX ];
	char delCommand[ FILENAME_MAX ];

	universe *worldPtr;
	worldPtr = new universe;

	// set default results and parameter directories
	strcpy( resultsDir, "./Results");
	strcpy( parametersDir, "./Parameters");

	if ( argc > 1 ){

		// read in parameters
		arg = 1;
		while( arg < argc )
		{

			//  -g	( number of iterations)			
			if( strcmp( argv[arg], "-g") == 0 )
			{
				arg++;
				if ( argv[arg] == NULL ){ delete worldPtr; usage();}
				worldPtr->setMaxEpochs( atoi( argv[arg] ) );
				goto loop;
			}

			// -v (verbosity)
			if( strcmp( argv[arg], "-v") == 0 )
			{
				arg++;
				if ( argv[arg] == NULL ){delete worldPtr; usage();}
				verbosity = atoi( argv[ arg ] );
				worldPtr->setVerbosity( verbosity );
				goto loop;
			}

			// -z (interactive visualisations)
			if( strcmp( argv[arg], "-z") == 0 )
			{
				arg++;
				if ( argv[arg] == NULL ){delete worldPtr; usage();}
				visualisation = static_cast<bool>(atoi( argv[arg] ));
				worldPtr->setVisualisation( visualisation );
				goto loop;
			}

			// -y (show very verbose data)
			if( strcmp( argv[arg], "-y") == 0 )
			{
				arg++;
				if ( argv[arg] == NULL ){delete worldPtr; usage();}
				veryVerbose = true;
				goto loop;
			}

			// -r (results directory)
			if( strcmp( argv[arg], "-r" ) == 0 )
			{
				arg++;
				if ( argv[arg] == NULL ){ delete worldPtr; usage();}
				strcpy( resultsDir, argv[arg] );
				worldPtr->setResultsDirectory( resultsDir );
				goto loop;
			}

			// -p (parameters directory)
			if( strcmp( argv[arg], "-p" ) == 0 )
			{
				arg++;
				if ( argv[arg] == NULL ){ delete worldPtr; usage();}
				strcpy( parametersDir, argv[arg] );
				worldPtr->setParametersDirectory( parametersDir );
				goto loop;
			}

			// -s (random seed offset)
			if( strcmp( argv[arg], "-s" ) == 0 )
			{
				arg++;
				if ( argv[arg] == NULL ){ delete worldPtr; usage();}
				sOffset = convertSeed( argv[arg] );
				goto loop;
			}

loop:
			arg++;
		}
	}

	// seed the random number generator	

 	srand( sOffset );

	// remove old results directory

	strcpy( delCommand, "rm -rf " );
	strcat( delCommand, resultsDir );

	if ( veryVerbose ){
		cout << "\nSeeding random number generator..."; 

		cout << "\nSeed = " << sOffset << "\nRandom number: " << rand() %1000 << endl;

		cout << "\nDeleting Results directory" << endl;
	}

	system( delCommand );

	// initialise the universe object
	worldPtr->init();

	// start simulation
	completedEpochs = worldPtr->runSimulation();

	// give output messages
	cout << "Run ended after " << completedEpochs << " epochs.\n" << endl;

	delete worldPtr;

	return ( 0 );

}


void usage( void )
{
        // print out error messages
        cout << "Incorrect usage: terminating." << endl;
        exit( 1 );
}


int convertSeed( char* inputString )
{
	char c,
		asciiNumber [ 256 ];

	int seed = 0,
		sLength = 0;

	// loop through string, determining whether character is digit or non digit
	for ( unsigned int i = 0; i < strlen( inputString ); i++ ){

		c = inputString[ i ];

		if ( isdigit ( inputString[ i ] ) ){
			int dgt = atoi( &c );
			seed += (int) ( dgt * ( pow( 10.0, ( double ) sLength ) ) );	
			sLength++;
		}
	}
	
	// convert the int to character string
	sprintf( asciiNumber, "%i", seed );
	// reverse the string
	strrev(	asciiNumber );
	// convert back to int
	seed = atoi( asciiNumber );

	return ( seed );

}


char *strrev( char *s )
{
	char *u, *v;
	char c;

	for ( u = s, v = s + strlen( s ) -1; u < v; ++u, --v ) {
		c = *u;
		*u = *v;
		*v = c;
	}

  return ( s );

}
