/*

bacteria.h
Header File for derrived bacteria class -- new prototype

adaptation selectable

(c) Neale Samways 2006-9

*/

#ifndef bacteria_H
#define bacteria_H

#include<iostream>
#include<ostream>
#include<iomanip>

#include"agent.h"

using namespace std;

class bacteria : public agent
{

public:


	bacteria();							// constructor
	~bacteria();						// destructor

	void setupReceptor( void );						// set interalUpdatePtr function
	void internalUpdate( double* );					// update function to be called

	int updateBehaviour( void );					// overwritten behavioural update routine


	inline int getReceptorType( void ){ return(  ( int )individualParametersPtr[ RECEPTOR_TYPE ] ); };				// return the agent type


private:

	void (bacteria::*internalUpdatePtr)( double* ) ;	// instance dependent RN update routine (RT assigned function pointer)
	void internalUpdateP( double* );				// Cycle RN / protein network : PTD model
	void internalUpdateI( double* );				// Cycle RN /protein network : IV model

	double *priorEnvironmentPtr;

	double alpha;


};



#endif
