/*

parameter.h
Header File for parameter class

(c) Neale Samways 2007-9

*/

#ifndef PARAMETER_H
#define PARAMETER_H

#include<iostream>
#include<cmath>
#include<cstring>

#include"global.h"

using namespace std;

class parameter
{

public:

	parameter();				// constructor
	~parameter();				// destructor

	bool initialise( int, char* );		// set up parameter object based on number of entries
	int iGet( int );			// return requested integer casted parameter
	double dGet ( int );			// return requested double casted parameter
	char* cGet( int );			// return requested char casted parameter

	void print( void );			// write parameter values (of all casts) to filestream

private:

	// functions
	bool loadParameters( char* );		// load parameters from file
	void removeComments( ifstream* );	// remove extraneous characters from input stream

	// variables
	int parameterCount;			// number of parameters contained in object

	double *dParameterPtr;			// pointer to double array of parameters
	int *iParameterPtr;			// pointer to integer array of parameters
	char **cParameterPtr;			// pointer to char array of parameters



};

#endif




