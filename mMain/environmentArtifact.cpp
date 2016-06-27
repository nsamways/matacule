/*

environmentArtifact.cpp

member functions and data members for environmentArtifact class

(c) Neale Samways 2006-9

*/

#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>
#include<cstring>

using std::setprecision;
using std::setiosflags;

#include"environmentArtifact.h"

using namespace std;

// static variables
const int environmentArtifact::numParameters = 13;

//functions
environmentArtifact::environmentArtifact()									// constructor
{

// integers
	worldWidth 	= 0;
	worldHeight 	= 0;
	worldEpoch 	= 0;
	xOffset		= 0;
	yOffset		= 0;

// doubles 
	dConstant = 0.0;
	deploymentVolume = 0.0;

// arrays
	memset( parameters, 0, sizeof parameters );
	

}


environmentArtifact::~environmentArtifact()									// destructor
{

	// dispose ladscape matrix
	for ( int i = 0; i < worldWidth; i++ )
		delete[] concentrationMatrix[ i ];
		
		delete[] concentrationMatrix;

	concentrationMatrix = NULL;

}


void environmentArtifact::setDimensions( int width , int height )						// set the dimensions of the artifact
{

	worldWidth = width;
	worldHeight = height;

}


void environmentArtifact::setParameter( int index, char *rawParameter )						// set the indexed parameter in the class
{
	// cast and set the parameter accordingly
	// nb. the first parameter is a string, thus not cast to double
	if ( index == 0 ){
		strcpy( friendlyName, rawParameter );			// set human readable label
	}else{	parameters[ index ] = atof( rawParameter );}

}


bool environmentArtifact::init( void )										// initialise the artifact 
{

	// create the matrix of concentration values
	concentrationMatrix = new double*[ worldWidth ];

	for ( int i = 0; i < worldWidth; i++ ){
		concentrationMatrix[ i ] = new double[ worldHeight ];
		
		for ( int j = 0; j < worldHeight; j++ )
			concentrationMatrix[ i ][ j ] = parameters[ BASAL_DIST ];
	}


	// calculate deployment offsets
	xOffset = ( int ) ( ( ( double ) worldWidth - parameters[ UPDATE_RANGE ] ) / 2.0 );
	yOffset = ( int ) ( ( ( double ) worldHeight - parameters[ UPDATE_RANGE ] ) / 2.0 );

	// deploy basal nutrient amount
	deploy( 0, 0, 0, 0, parameters[ BASAL_DIST ], 0.0 );

	// set the diffusion constant
	dConstant = parameters[ DIFFUSIVITY_CONST ];

	// determine periodic volume of deployment
	// 1 - linear
	if ( parameters[ UPDATE_CONFIG ] == 1 ){
		deploymentVolume = PI * parameters[ DEPLOY_RADIUS ] * parameters[ DEPLOY_RADIUS ] * parameters[ DEPLOY_MAX ] / 3.0;
	}

	// 2 - Gaussian
	if ( parameters[ UPDATE_CONFIG ] == 2 ){
		deploymentVolume = 2 * PI * parameters[ DEPLOY_RADIUS ] * parameters[ DEPLOY_RADIUS ] * parameters[ DEPLOY_MAX ];
	}

	return( true );
}


double environmentArtifact::localUpdate( double xLocus, double yLocus, double value )				// update the concentration at a given locus
{
	// change local concentration according to individual's consumption
	double energy,
		consumed;

	// convert double to integer: 
	int Xc = ( int ) floor( xLocus );
	int Yc = ( int ) floor( yLocus );

	// calculate amount of artifact consumed; lower boundary of zero automatically enforced
	consumed = min( concentrationMatrix[ Xc ][ Yc ] , value );

	concentrationMatrix[ Xc ][ Yc ] -= consumed;

	// return the amount of energy gained from the transaction
	energy = consumed * parameters[	ENERGY_YIELD ];

	return( energy );
}


void environmentArtifact::globalUpdate( int worldEpoch )							// update the entire artifact (i.e. diffuse, possible deployment)
{
// change complete environment according to epoch
	double cX,
		cY;

	// check for new deployment
	if ( worldEpoch % ( int )parameters[ UPDATE_PERIOD ] == 0 ){
		
		if ( parameters[ UPDATE_REGIME ] > 0.0 ){
			// plot centroid coordinates according to simple reciprocating function	
			// use NUTRIENT REGIME value to determine periodicity
			cX = ( parameters[ UPDATE_RANGE ] / 2.0 ) * sin( PI * ( worldEpoch / parameters[ UPDATE_PERIOD ] ) / parameters[ UPDATE_REGIME ] ) + ( worldWidth / 2.0 );
			cY = ( parameters[ UPDATE_RANGE ] / 2.0 ) * sin( 2.0 * PI * (  worldEpoch / parameters[ UPDATE_PERIOD ] ) / parameters[ UPDATE_REGIME ] ) + ( worldHeight / 2.0 ); 
		}else{
			if ( parameters[ UPDATE_REGIME ] == 0.0 ){
			// plot centroid randomly within range
				cX = xOffset + ( ( double ) rand() /  RAND_MAX ) * parameters[ UPDATE_RANGE ];
				cY = yOffset + ( ( double ) rand() /  RAND_MAX ) * parameters[ UPDATE_RANGE ];
			}else{
			// plot circular function
			cX = ( parameters[ UPDATE_RANGE ] / 2.0 ) * cos( (  2 * PI * ( worldEpoch /  parameters[ UPDATE_PERIOD ] )  / -parameters[ UPDATE_REGIME ]  ) ) + ( worldWidth / 2.0 );
			cY = ( parameters[ UPDATE_RANGE ] / 2.0 ) * sin( ( 2* PI * ( worldEpoch  /  parameters[ UPDATE_PERIOD ] )  /  -parameters[ UPDATE_REGIME ]  ) ) + ( worldHeight / 2.0 ); 

			}

		}

		// make any uniform deployment
		if (  parameters[ BASAL_DIST ] ) deploy( 0, 0, 0, parameters[ BASAL_DIST ], 0, 0);
//NOTE:
		// make gradient deployment
		deploy( ( int ) cX, ( int ) cY, ( int ) parameters[ UPDATE_CONFIG ], parameters[ DEPLOY_MAX ], parameters[ DEPLOY_RADIUS ], parameters[ DEPLOY_DELTA ] ) ;
	}
	// diffuse artifact ( according to diffusion constant ) 
	diffuse();

}


void environmentArtifact::primaryDeployment( double cX, double cY )	// place deployment at specified location
{
	deploy( ( int ) cX, ( int ) cY, ( int ) parameters[ UPDATE_CONFIG ], parameters[ DEPLOY_MAX ], parameters[ DEPLOY_RADIUS ], parameters[ DEPLOY_DELTA ] ) ;


}


void environmentArtifact::deploy( int centerX, int centerY, int type, double maxAmp, double sigma, double delta )	// make a deployment of artifact according to parameters
{

// put a deployment of chemical in to the environment

	double distance;

	switch( type ){
		
		case 0:{
			// uniform (constant) 
			
			for ( int i = 0; i < worldWidth; i++ ){
				for ( int j = 0; j < worldHeight; j++ )
					concentrationMatrix[i][j] += maxAmp;
			}

			break;
		}
		
		case 1:{	
			// linear increment

			// calculate parameters depending on delta values

			// determine current Amp / sigma
			currentSigma = sigma + ( ( 2.0 * ( ( ( double ) rand() / ( RAND_MAX ) ) -0.5 ) ) * parameters[ DEPLOY_DELTA ] );
			currentAmp = 3.0 * deploymentVolume / ( PI * currentSigma * currentSigma );
 			

			double increment = 0;
			double grad = - ( currentAmp / currentSigma );

			for ( int i = 0; i < worldWidth; i++ ){
				for ( int j = 0; j < worldHeight; j++ ){

				distance = getDistance( i, j, centerX, centerY );

				increment = grad * distance + currentAmp;
				increment = ( increment > 0 ) ? increment : 0;		

				concentrationMatrix[ i ][ j ] += increment;	

				}
			}

			break;
		}			

		case 2:{
			// indexed Gaussian

			currentSigma = sigma + ( ( 2.0 * ( ( ( double ) rand() / ( RAND_MAX ) ) -0.5 ) ) * parameters[ DEPLOY_DELTA ] );
			currentAmp = deploymentVolume / ( 2.0 * PI * currentSigma * currentSigma );

			for ( int i = 0; i < worldWidth; i++ ){
			
				for ( int j = 0; j < worldHeight; j++ ){
		
					distance = getDistance( i, j, centerX, centerY );  
					concentrationMatrix[ i ][ j ] += currentAmp * exp(-(( distance * distance) / ( 2 * ( currentSigma * currentSigma ))));

				}
			}	

			break;
		}		

		case 3:{
			// uniform probabilistic distribution
			// sigma is the probability per site, max amp is the size of the deployment

			for ( int i = 0; i < worldWidth; i++ ){
				for (int j = 0; j < worldHeight; j++){
					if ( ( ( double ) rand() / RAND_MAX ) > sigma ) 
						concentrationMatrix[i][j] += maxAmp;
				}
			}
			
			break;
		}
	}

}


void environmentArtifact::deposit( int xc, int yc, double amount )						// increase local concentration at the identified locus
{	
	concentrationMatrix[ xc ][ yc ] += abs( amount );
}


void environmentArtifact::decay( double volume )								// globally reduce the concentration of artifact
{

	for ( int i = 0; i < worldWidth; i++ ){
	
		for ( int j = 0; j < worldHeight; j++ ){
			concentrationMatrix[ i ][ j ] -= volume;
			concentrationMatrix[ i ][ j ] = ( concentrationMatrix[ i ][ j ] < 0 ) ? 0 : concentrationMatrix[ i ][ j ];
		}
	}	

}


void environmentArtifact::move( int xShift, int yShift )							// move all concentrations in the concentration matrix by the identified amount
{
	// shift the peak of chemical by x,y
	
	// create holding matrix
	double **holdingMatrixPtr = new double*[ worldHeight ];

	int k,
		l;

	for ( int i = 0; i < worldWidth; i++ ){

		holdingMatrixPtr[ i ] = new double[ worldHeight ];

		for ( int j = 0; j < worldHeight; j++ ){
			
			k = i - xShift; 

			l = j - yShift;

			if ( k < 0 ) k += worldWidth;
			if ( l < 0 ) l += worldHeight;

			holdingMatrixPtr[ i ][ j ] = concentrationMatrix[ k % ( worldWidth ) ][ l % ( worldHeight ) ];
		}
	}

	// delete old concentrationMatrix and reassign pointer to new matrix

	for ( int i = 0; i < worldWidth; i++ ) delete[] concentrationMatrix[i];
		
		delete[] concentrationMatrix;

		concentrationMatrix = holdingMatrixPtr;

}


double environmentArtifact::getSumVolume( void )								// return the aggregate amount of artifact in the environment
{
	double total = 0.0;

	for ( int i = 0; i < ( worldWidth  ); i++ ){
		for ( int j = 0; j < ( worldHeight ); j++ )
			total += concentrationMatrix[ i ][ j ];		
	}

	return( total );
}


double environmentArtifact::reportGradient( void )								// return a rating for the amount of "gradient" information in the environment
{

	double gradientScore = 0.0;
	int windowSize	=	5;
	bool gradient = false;
	bool partGradient = true;
	int l;
	int k;

	// function to show how much "gradient info" is in the environmentArtifact

/*
	solution 1: 

	for all squares, move in the cardinal directions as far as possible, until the siqn of the gradient changes, adding one for every
	square moved along.

	solution 2:

	for all squares, calculate if there is any possible gradient within the window size.If there is a gradient (in any direction), increment
	the score by one. This is essentially a count of all positions in which it would not be clear where to move.
*/


	// start from "windowSize" units in

	for ( int i = windowSize ; i < worldWidth - windowSize; i++ ){
		for ( int j = windowSize; j < worldHeight - windowSize; j++ ){

				gradient = false;

				// 'UP '

				partGradient = true;
				l = 1;
				while( ( partGradient == true ) && ( windowSize >= l ) ){
					k  = l - 1;
					if ( concentrationMatrix[ i ][ j + k ] >= concentrationMatrix[ i ][ j + l ] ) partGradient = false;
					l++;	
				}
				
				if ( partGradient ) gradient = true;

				// 'DOWN '

				partGradient = true;
				l = 1;
				while( ( partGradient == true ) && ( windowSize >= l ) ){
					k  = l - 1;
					if ( concentrationMatrix[ i ][ j - k ] >= concentrationMatrix[ i ][ j - l ] ) partGradient = false;
					l++;	
				}
				
				if ( partGradient ) gradient = true;

				// 'LEFT'

				partGradient = true;
				l = 1;
				while( ( partGradient == true ) && ( windowSize >= l ) ){
					k  = l - 1;
					if ( concentrationMatrix[ i - k ][ j ] >= concentrationMatrix[ i - l ][ j ] ) partGradient = false;
					l++;	
				}
				
				if ( partGradient ) gradient = true;


				// 'RIGHT'

				partGradient = true;
				l = 1;
				while( ( partGradient == true ) && ( windowSize >= l ) ){
					k  = l - 1;
					if ( concentrationMatrix[ i + k ][ j ] >= concentrationMatrix[ i + l ][ j ] ) partGradient = false;
					l++;	
				}
				
				if ( partGradient ) gradient = true;


				// now have looped through all cardinalities
				if ( !gradient ) gradientScore++;
		}
	}
	
	return( gradientScore );

}


void environmentArtifact::print( ofstream* outfile )								// output the contents of the concentration matrix (i.e. the entire distribution)
{

	for (int i = 0; i < worldWidth; i++){
		for (int j = 0; j < worldHeight; j++){
			*outfile << setprecision(3) << setiosflags(ios::fixed | ios::showpoint) << concentrationMatrix[i][j] << "\t";
		}
		*outfile << "\n";
	}
	
	*outfile << endl;
}


void environmentArtifact::displayUpdate( FILE* displayPtr )							// pipe the contents of the concentration matrix for display					
{
	fprintf( displayPtr,"set pm3d \n" );		
	
	fprintf( displayPtr, "splot '-' matrix\n" );
	for ( int i = 0; i < worldHeight; i++ ){
		for ( int j = 0; j < worldWidth; j++ ) fprintf( displayPtr,"%f\t",concentrationMatrix[ j ][ i ]);
		
		fprintf( displayPtr,"\n" );
	}
	
	fprintf( displayPtr, "e\n" );
	fprintf( displayPtr,"set autoscale z \n" );
	fflush( displayPtr );
}


void environmentArtifact::printConfiguration( ofstream* reportFile )						// write the configuration of the artifact to filestream
{

	*reportFile << "---\n"
		<< "Configuration for environmental artifact: " << friendlyName
		<< "\nDiffusivity constant:\t"
		<< parameters[ DIFFUSIVITY_CONST ]
		<< "\nPermeation factor: \t"
		<< parameters[ PERMEATION_FACTOR ]
		<< "\nEnergy yield \t"
		<< parameters[ ENERGY_YIELD ]
		<< endl;


}


// Private functions


void environmentArtifact::diffuse( void )									// diffuse the artifact 
{

	// diffuse chemical in the environment by the Euler forward method
	
	// create holding matrix + variables

	double cU, cD, cL, cR;

	double **holdingMatrix = new double*[ worldWidth ];

	for ( int i = 0; i < worldWidth; i++ ){

		holdingMatrix[ i ] = new double[ worldHeight ];

		for ( int j = 0; j < worldHeight; j++ )
			holdingMatrix[ i ][ j ] = 0.0;
	
	}
	
		// create the diffusion in the new matrix
	// inner sqaure first 
/*	cout << "World width  = " << worldWidth << endl;*/
	for ( int i = 1; i < ( worldWidth - 1 ); i++ ){
		for ( int j = 1; j < ( worldHeight - 1 ); j++ ){
/*cout << "i = " << i << " : j = " << j << endl;*/
			cU = concentrationMatrix[ i ][ j + 1 ];
			cD = concentrationMatrix[ i ][ j - 1 ];	
			cL = concentrationMatrix[ i - 1 ][ j ];	
			cR = concentrationMatrix[ i + 1 ][ j ];	

			holdingMatrix[ i ][ j ] = concentrationMatrix[ i ][ j ] + ( dConstant * ( -4 *  concentrationMatrix[ i ][ j ] + cU + cD + cL + cR ));
		}
	}

	// corners 
		// bottom left

	holdingMatrix[ 0 ][ 0 ] = concentrationMatrix[ 0 ][ 0 ] + ( dConstant * ( -4 * concentrationMatrix[ 0 ][ 0 ] + ( 2 * concentrationMatrix[ 1 ][ 0 ] ) + ( 2* concentrationMatrix[ 0 ][ 1 ] ) ) );

		// top left

	holdingMatrix[ 0 ][ worldHeight - 1 ] = concentrationMatrix[ 0 ][ worldHeight -1 ] + ( dConstant * ( -4 * concentrationMatrix[ 0 ][ worldHeight -1 ] + ( 2 * concentrationMatrix[ 0 ][ worldHeight - 2 ]) + ( 2 * concentrationMatrix[ 1 ][ worldHeight -1 ]) ) );

		// top right

	holdingMatrix[ worldWidth -1 ][ worldHeight -1 ] = concentrationMatrix[ worldWidth -1 ][ worldHeight -1 ] + ( dConstant * ( -4 * concentrationMatrix[ worldWidth -1 ][ worldHeight -1 ] + ( 2 * concentrationMatrix[ worldWidth - 2 ][ worldHeight -1 ]) + ( 2 * concentrationMatrix[ worldWidth -1 ][ worldHeight - 2] ) ) );

		// bottom right

	holdingMatrix[ worldWidth -1 ][ 0 ] = concentrationMatrix[ worldWidth -1 ][ 0 ] + ( dConstant * ( -4 * concentrationMatrix[ worldWidth -1 ][ 0 ] + ( 2 * concentrationMatrix[ worldWidth  -2 ][ 0 ] ) + ( 2 * concentrationMatrix[ worldWidth -1 ][ 1 ] ) ) );


	// edges
		// top edge

	for ( int i = 1; i < ( worldWidth - 1 ); i++ ){

			cD = 2 * concentrationMatrix[ i ][ 1 ];	
			cL = concentrationMatrix[ i - 1 ][ 0 ];	
			cR = concentrationMatrix[ i + 1 ][ 0 ];	

			holdingMatrix[ i ][ 0 ] = concentrationMatrix[ i ][ 0 ] + ( dConstant * ( -4 *  concentrationMatrix[ i ][ 0 ] + cD + cL + cR ));
	}
		// bottom edge

	for ( int i = 1; i < ( worldWidth - 1 ); i++ ){

			cU = 2 * concentrationMatrix[ i ][ worldHeight - 2 ];	
			cL = concentrationMatrix[ i - 1 ][ worldHeight - 1];	
			cR = concentrationMatrix[ i + 1 ][ worldHeight - 1];	

			holdingMatrix[ i ][ worldHeight - 1] = concentrationMatrix[ i ][ worldHeight - 1] + ( dConstant * ( -4 *  concentrationMatrix[ i ][ worldHeight -1 ] + cU + cL + cR ));
	}

		// left edge

	for ( int i = 1; i < ( worldHeight - 1 ); i++ ){

			cU = concentrationMatrix[ 0 ][ i - 1];
			cD = concentrationMatrix[ 0 ][ i + 1 ];	
			cR = 2 * concentrationMatrix[ 1 ][ i ];	

			holdingMatrix[ 0 ][ i ] = concentrationMatrix[ 0 ][ i ] + ( dConstant * ( -4 *  concentrationMatrix[ 0 ][ i ] + cU + cD + cR ));
	}
		// right edge

	for ( int i = 1; i < ( worldHeight - 1 ); i++ ){

			cU = concentrationMatrix[ worldWidth -1 ][ i - 1 ];
			cD = concentrationMatrix[ worldWidth - 1][ i + 1 ];	
			cL = 2 * concentrationMatrix[ worldWidth - 2 ][ i ];	

			holdingMatrix[ worldWidth -1 ][ i ] = concentrationMatrix[ worldWidth - 1][ i ] + ( dConstant * ( -4 *  concentrationMatrix[ worldWidth -1 ][ i ] + cU + cD + cR ));
	}
 
	// delete old concentrationMatrix and reassign pointer to new matrix

	for ( int i = 0; i < worldWidth; i++ ) delete[] concentrationMatrix[i];
		
		delete[] concentrationMatrix;

		concentrationMatrix = holdingMatrix;
	
}


double environmentArtifact::getPeriodicDistance( int currentX, int currentY, int targetX, int targetY)			// return the distance between given points
{

	double x1 = 0.0;
	double x2 = 0.0;
	double x = 0.0;

	double y1 = 0.0;
	double y2 = 0.0;
	double y = 0.0;

	double dist = 0.0;
	
	x1 = ( double ) abs(( currentX - targetX ));
	y1 = ( double ) abs(( currentY - targetY ));

	if ( currentX > targetX ){
		x2 = ( double )( targetX + worldWidth - currentX );
	}else{
		x2 = ( double )( currentX + worldWidth - targetX );
	}
	
	if ( currentY > targetY ){
		y2 = ( double )( targetY + worldHeight - currentY );
	}else{
		y2 = ( double )( currentY + worldHeight - targetY );
	}
	
	x = min( x1, x2 );
	y = min( y1, y2 );

	dist = sqrt( ( x * x ) + ( y * y ) );

	return ( dist );
}

double environmentArtifact::getDistance( int currentX, int currentY, int targetX, int targetY)			// return the distance between given points
{

	double dist =  sqrt( ( ( targetX - currentX ) * ( targetX - currentX ) ) + ( targetY - currentY ) * ( targetY - currentY ) );

	return ( dist );
}
