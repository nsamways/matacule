/*

genotype.cpp

member functions and data members for genotype class

(c) Neale Samways 2006-9

*/ 

#include<iostream>
#include<iomanip>
#include<fstream>

#include"genotype.h"

using namespace std;

	int genotype::totalReceptors = 0;

genotype::genotype()								// constructor
{
	proteomeSize	=	0;
	totalGenes 	=	0;
	totalOperons 	=	0,
	inputProteinNumber  =	0;

	chemistryMutationProb	=	0.0;

	operonPtr 	=	NULL;
	chemIntPtr	=	NULL;
	proteomePtr	=	NULL;

	record 		=	false;

}


genotype::~genotype()								// destructor
{
	
	// termination housekeeping	

	delete [] operonPtr;
	operonPtr = NULL;
	
	// also deallocate the proteome

	delete [] proteomePtr;
	proteomePtr = NULL;

}


void genotype::initialise( genome *genomePtr )					// set up the genotype object
{

	// determine all present operons - promotor sequence is: 0,15

	int genomeLength = genomePtr->getLength();

	int positionMarker[ genomeLength ];			

	for ( int i = 0; i < genomeLength - ( OPERON_LENGTH - 1 ); i++ ){ 		// check to the last possible operon position
		if ( genomePtr->getBase( i ) == 0 ){
			if ( genomePtr->getBase( i + 1 ) == 15 ){
				// we have located a valid operon

				positionMarker[ totalOperons ] = i;			// notify that an operon is presnt starting at this loci
				++totalOperons;
				// we can skip the next set of base sequences
				i += ( OPERON_LENGTH - 1 ); // This cannot exceed the array bounds
			}
		}
	}

// all operon positions stored in positionMarker[]

	// initialise the operons, and determine total number of genes in genome
	operonPtr = new operon[ totalOperons ];

	for ( int i = 0; i < totalOperons; i++ ){
		operonPtr[ i ].init( positionMarker[ i ], genomePtr );	// initialise the operons by passing values to operon
		totalGenes += operonPtr[ i ].getNumGenes();		// sum the number of genes in the genome
	}

}


void genotype::setChemistry( chemistry* cPtr, double cmProb )		// set the pointer to the chemistry and set mutation probability
{

	chemIntPtr = cPtr;
	chemistryMutationProb = cmProb;

}


void genotype::setupProteome( protein *parentProteomePtr )			// set array of producible proteins and mutate chemistry for offspring
{
// set up the proteome

	proteomePtr = new protein[ TOTAL_PROTEINS ];			// proteome object 

	bool proteinRegister[ TOTAL_PROTEINS ];					// TRUE at specific index indicates existing genes 
	
	// initialise indicator register to zero - i.e. no chemicals present
	for ( int j = 0; j < TOTAL_PROTEINS; j++ ) proteinRegister[ j ] = false;

	// if the parent argument is not NULL, inherit parent proteomic concentrations
	if ( parentProteomePtr ){

		for ( int i = 0; i < TOTAL_PROTEINS; i++ ){
			// inherit to new proteome
			proteomePtr[ i ].setConcentration( parentProteomePtr[ i ].getConcentration() );
		}
	}

	// initialisation for initial proteomes begins here (passed pointer is not NULL)

	// ensure protein register entries for all genes ( from all operons ) are made 
	for ( int i = 0; i < totalOperons; i++ ){
		for ( int j = 0; j < INIT_GENE_NUM; j++ )
			proteinRegister[ operonPtr[ i ].getProteinIdent( j ) ] = true;
	}

	// include any necessary proteins here:
// 	proteinRegister[ INT_ENERGY_PROTEIN ] = true;

	// set the active proteins vector
	for ( int i = 0; i < TOTAL_PROTEINS; i++ )
		if ( proteinRegister[ i ] == true ) activeProteins.push_back( i );

	proteomeSize = activeProteins.size();

	// set the initial concentrations if from starting population,else mutate chemistry according to active proteins
	if ( !parentProteomePtr ){
		for ( int i = 0; i < TOTAL_PROTEINS; i++ ){
			if ( proteinRegister[ i ] == true ) proteomePtr[ i ].setConcentration( INITIAL_PROTEIN_CONC );
		}
	}else{	// mutate the chemistry
		for ( unsigned int i = 0; i < activeProteins.size(); i++ ){
			for ( unsigned int j = 0; j < activeProteins.size() ; j++ ) chemIntPtr->mutate( activeProteins[ i ] , activeProteins[ j ], chemistryMutationProb );
		}
	}

	// set the instances	
	for ( int i = 0; i < totalOperons; i++ ){
		for ( int j = 0; j < INIT_GENE_NUM; j++ ) proteomePtr[ operonPtr[ i ].getProteinIdent( j ) ].increaseInstances();
	}	

	// identify proteins in external file
	if ( record ){

		*proteomeRegister << "# ";
		for ( int j = 0; j < proteomeSize; j++ ) *proteomeRegister <<  activeProteins[ j ] << "\t";
			
		*proteomeRegister  << "\n";
	}

}


void genotype::resetProteomeConc( void )					// set the proteome concentrations back to base level
{

	for ( int i = 0; i < TOTAL_PROTEINS; i++ )
		proteomePtr[ i ].setConcentration( INITIAL_PROTEIN_CONC );

}


bool genotype::loadProteomeConcs( char* basePathPtr )
{

	// load protein concentrations from file

	double index = -1.0;
	ifstream inStream;

	char filePathPtr[ FILENAME_MAX ];	// absolute name of strain directory

	strcpy( filePathPtr, basePathPtr );
	strcat( filePathPtr,"/proteomeConc.csv" );

	// open the concentration file
	
	inStream.open( filePathPtr );
	
	// check the the file can be opened 	
	if ( inStream.fail() ){
		cerr << "Unable to open initial protein concentration file:" << filePathPtr << "  ";
		return( false );
	}

	for ( int i = 0; i < proteomeSize; i++ ){

		removeComments( &inStream);
		inStream >> index;
		proteomePtr[ activeProteins[ i ] ].setConcentration( index );	
	}

	inStream.close();

	return( true );

}


void genotype::updateInputConc( double *inputPtr )
{

      	double interaction = 0.0;						// Interaction value	

	memset( phiDeltaP, 0, sizeof phiDeltaP );				// set the array of updates to zero
	memset( phiDeltaN, 0, sizeof phiDeltaN );				// set the array of updates to zero

	// determine epsilon values from chemical reactions matrix
	// calculate changes for input proteins i, by all proteins j
	   for ( int i = 0; i < totalReceptors; i++ ){
		for ( int j = 0; j < proteomeSize; j++ ){

			interaction = chemIntPtr->getInteraction( i , activeProteins[ j ] );

			if ( interaction > 0.0 ){
				phiDeltaP[ i ] += interaction * ( proteomePtr[ activeProteins[ j ] ].getConcentration() );
			}else{
				phiDeltaN[ i ] += interaction * ( proteomePtr[ activeProteins[ j ] ].getConcentration() );
			}

		}
	}

	   for ( int i = 0; i < totalReceptors; i++ ){
  		// determine proper values of phiDelta
		phiDeltaP[  i  ] *= ( 1.0 - proteomePtr[ i ].getConcentration()  ) ; 
		phiDeltaN[  i  ] *= proteomePtr[ i ].getConcentration() ; 

		// set concentration
		proteomePtr[ i ].setConcentration(  inputPtr[ i ] );

		// modify concentration based on interactions
		proteomePtr[  i  ].modifyConcentration( phiDeltaP[ i ] + phiDeltaN[ i ]  );
	    }

}


void genotype::cycleRN( void )							// cycle through network, updating concentrations
{

	double interaction = 0.0;							// Interaction value	

	memset( phiDeltaP, 0, sizeof phiDeltaP );				// set the array of updates to zero
	memset( phiDeltaN, 0, sizeof phiDeltaN );				// set the array of updates to zero

	// determine epsilon values from chemical reactions matrix -- use the activeProteins values to avoid running through entire proteome
	// calculate changes for protein i, by all proteins j
	for ( int i = totalReceptors; i <  proteomeSize; i++ ){
		for ( int j = 0; j < proteomeSize; j++ ){

			interaction = chemIntPtr->getInteraction( activeProteins[ i ], activeProteins[ j ] );

			if ( interaction > 0.0 ){
				phiDeltaP[ activeProteins[ i ] ] += interaction * ( proteomePtr[ activeProteins[ j ] ].getConcentration() );
			}else{
				phiDeltaN[ activeProteins[ i ] ] += interaction * ( proteomePtr[ activeProteins[ j ] ].getConcentration() );
			}

		}
	}
	
	// modify levels by ammount determined by chemical reactions matrix
	for ( int i = 0; i < proteomeSize; i++ ){
		// determine proper values of phiDelta
		phiDeltaP[ activeProteins[ i ] ] *= ( 1.0 - proteomePtr[ activeProteins[ i ] ].getConcentration() ); 
		phiDeltaN[ activeProteins[ i ] ] *= proteomePtr[ activeProteins[ i ] ].getConcentration() ; 
		// modify concentration
		proteomePtr[ activeProteins[ i ] ].modifyConcentration( phiDeltaP[ activeProteins[ i ] ] + phiDeltaN[ activeProteins[ i ] ] );

	}

}


void genotype::printProteinConcs( void )
{

	for ( int j = 0; j < proteomeSize; j++ ){ *proteomeRegister << proteomePtr[ activeProteins[ j ] ].getConcentration() << "\t"; }
	*proteomeRegister  << "\n";

}


void genotype::print( ofstream *outFilePtr  )					// print out operon details
{

	// display the contents of the genotype - I.e. show functionality
	*outFilePtr << "\nTotal Operons: " << totalOperons << endl;
	
	for (int i = 0; i < totalOperons; i ++){
		*outFilePtr << "\nOperon " << (i + 1) << endl;
		operonPtr[ i ].print( outFilePtr );
	}

	printProteome( outFilePtr );

}


void genotype::printProteome( ofstream *outFilePtr )				// print out proteome details
{
	*outFilePtr << "\nProteome:\n"
	<< "\nNumber of genes = " << totalGenes << endl;
	*outFilePtr << "Number of distinct genes = " << proteomeSize << endl;
	*outFilePtr << "\n\n#Gene Index:\t\tIdentity:\t\tConcentraion:\t\tCount:" << endl;
	for ( int i = 0; i < proteomeSize; i++ )
		*outFilePtr << " " << i << "\t" <<  activeProteins[ i ] << "\t" << proteomePtr[ activeProteins[ i ] ].getConcentration() << "\t" << proteomePtr[ activeProteins[ i ] ].getInstances() << endl;
	
	*outFilePtr << "#End Index" << endl;
}


int genotype::getConnectionCount( void )					// count number of non--zero connections in matrix
{

	int numConnections = 0;
	
	for ( int i = 0; i < proteomeSize; i++ ){
		for ( int j = 0; j < proteomeSize; j++ )
			if ( chemIntPtr->getInteraction( i, j ) != 0.0 ) numConnections++;
	}

	return( numConnections );

}


double genotype::getConnectionDensity( void )					// return the summed absolute value of weighted connections
{

	double summedWeight = 0.0;
	
	for ( int i = 0; i < proteomeSize; i++ ){
		for ( int j = 0; j < proteomeSize; j++ )
			summedWeight += abs( chemIntPtr->getInteraction( i, j ) );
	}

	return( summedWeight );


}


bool genotype::checkConnectionPath( int source, int target )			// return true if direct *or* indirect connection exists between source and target
{

	vector< int > validVector;			// list of valid nodes left to to try
	vector< int > visitedVector;			// list of already visited nodes

	bool notConnected = true;			// flag showing whether the target node has not been found
	bool visited = false;
	bool existing = false;

	for ( int i = 0; i < proteomeSize; i++ )									// start from source node
		if ( chemIntPtr->getInteraction( activeProteins[ i ], source ) != 0.0 ) validVector.push_back( activeProteins[ i ] );
	
	visitedVector.push_back( source );					// target is the first visited node so push on to visited stack

	// now go through all nodes in validVector 
	while( ( validVector.size() > 0 ) && notConnected ){
		int k = validVector.back();		// makes it easier -- k is current targeted node
		visitedVector.push_back( k );		// we're visiting this now
		validVector.pop_back();			// take it off the stack too

		visited = false;			// reset the visited flag; by default assume a node has not been visited

		if ( k == target )			// check if it is the target node
			notConnected = false;
		
			
		// check if already visited
		for ( unsigned int m = 0; m < visitedVector.size(); m++ )
			if ( visitedVector[ m ] ==  k ) visited = true; 
		
		if ( !visited ){						// not checked this node before, so push unvisited connections on to stack
			
			for ( unsigned int j = 0; j < activeProteins.size(); j++ ){				// add its unvisited connections to the valid vector
				if ( chemIntPtr->getInteraction( activeProteins[ j ], k ) != 0.0 ){  	// there is a connection
					existing = false;						// reset flag
					for ( unsigned int m = 0; m < visitedVector.size(); m++ )		// check if it has been visited
						if ( visitedVector[ m ] ==  activeProteins[ j ] ) existing = true; 			

					for ( unsigned int m = 0; m < validVector.size(); m++ ) // also check if already added!
						if ( validVector[ m ] == activeProteins[ j ] ) existing = true;  	

					if ( !existing )
						validVector.push_back( activeProteins[ j ] );	// not visited node previously, so push on to stack
				}
			}
		}// we've checked for a connection to the target node, reduced validVector by one node, and possibly added more 
	} 
	
	// either everywhere is visited, or we have found a connection

	return( !notConnected );

}


void genotype::removeComments( ifstream* inputFile )					// remove extraneous data from the stream
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
