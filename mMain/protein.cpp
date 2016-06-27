/*

protein.cpp

member functions and data members for protein class

(c) Neale Samways 2007-9

*/

#include<iostream>

#include"protein.h"

using namespace std;

protein::protein()							// constructor
{

	concentration = 0.0;
	geneticInstances = 0;

}


protein::~protein( void )						// destructor
{


}


void protein::modifyConcentration( double deltaConcentration )		// set concentration
{

	concentration += deltaConcentration;

	// apply bounds to concentration
	// volumes cannot be negative, or above the preset threshold

	concentration = ( concentration < 0.0 ) ? 0.0 : concentration;
	concentration = ( concentration > MAX_PROTEIN_VOLUME ) ? MAX_PROTEIN_VOLUME : concentration;

}


void protein::setConcentration( double newConcentration )		// change concentration
{
	concentration = newConcentration;

	// apply bounds to concentration
	// volumes cannot be negative, or above the preset threshold

	concentration = ( concentration < 0.0 ) ? 0.0 : concentration;
	concentration = ( concentration > MAX_PROTEIN_VOLUME ) ? MAX_PROTEIN_VOLUME : concentration;

}
