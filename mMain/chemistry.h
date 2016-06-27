/*

chemistry.h
Header File for chemistry class

(c) Neale Samways 2006-9

*/

#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<fstream>

#include"global.h"

using namespace std;


class chemistry
{

public:

	chemistry();					// constructor
	~chemistry();					// destructor

	bool load( char* );				// load interaction matrix from file

	// initialisers
	void randomizeInteractions( double );		// Randomize the interaction matrix
	void mutate( int , int , double );		// add Gaussian noise to specified interaction

	// printers
	void print( ofstream* );			// Write matrix to file stream

	// copy operator
	void operator= ( chemistry* );			// Copy given interaction matrix

	// Inline getters
	inline double getInteraction( int p, int q ){ return( reactions[ p ][ q ] ); }					// return specified interaction

	// inline setters
	inline void setInteraction( int p, int q, double rRate ){ reactions[ p ][ q ] = rRate; }			// set interaction according to arguments
	inline static void setInteractionDeletionProb( double dProb ){ chemistry::interactionDeletionProb = dProb; }	// set interactino deletion probability

private:

	// functions
	void removeComments( ifstream* );					// discard extraneous data from input stream
	double generateGaussian( double, double );				// return a Gaussian distributed number

	// variables
	double** reactions;							// interaction matrix

	static double interactionDeletionProb;					// interaction deletion probability

};

#endif




