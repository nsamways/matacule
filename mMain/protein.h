/*

protein.h
Header File for protein class

(c) Neale Samways 2006-9

*/

#ifndef PROTEIN_H
#define PROTEIN_H

#include<iostream>

#include "global.h"

using namespace std;

class protein
{

public:

	protein();							// constructor
	~protein();							// destructor


	// setters	
	void	setConcentration ( double );				// set concentration
	void	modifyConcentration( double );				// change concentration


	// Inlines
	inline int getInstances(){ return geneticInstances; }		// return number of genetic instances of specified protein
	inline double getConcentration(){ return concentration; }	// return concentration of protein
	inline void increaseInstances(){ ++geneticInstances; }		// increment the number of genetic instances of protein

private:

	int geneticInstances;						// number of genes coding for this protein
	double concentration;						// concentration of protein

};

#endif
