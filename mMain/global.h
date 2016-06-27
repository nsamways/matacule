/*

global.h
global parameters, variables and enumerations (Most classes depend on this file)

(c) Neale Samways 2006-9

*/

#ifndef GLOBAL_H
#define GLOBAL_H

#define PI 3.1415927

// Hard-coded variables
#define PROMOTOR_LENGTH	2
#define CIS_REGION_LENGTH	0
#define INIT_GENE_NUM		3
#define TOTAL_PROTEINS		16
#define BASAL_INCREMENT_LEVEL	1
#define OPERON_LENGTH ( PROMOTOR_LENGTH + CIS_REGION_LENGTH + INIT_GENE_NUM )
#define MAX_PROTEIN_VOLUME	1.0
#define INITIAL_PROTEIN_CONC	0.5

#define NUMBER_INDIVIDUAL_PARAMS	13
#define NUMBER_SIMULATION_PARAMS	6
#define NUMBER_ENVIRONMENT_PARAMS	14

#define NUMBER_BEHAVIOURAL_PROTEINS	1
#define INITIAL_ENERGY	24
#define MAX_ENERGY	25
// Hard-coded functions
#define QUICKRAND		( rand() / ( double )  RAND_MAX )

// Simulation parameters
enum
{
	POP_COMPOSITION = 0,		//	(0)
	RES_FREQ,					//	(1)
	R_START,					//	(2)
	R_END,						//	(3)
	INIT_POP_LOC,				//	(4)
	MAX_EPOCHS				//	(5)
	
};


// Individual parameters
enum
{
	RECEPTOR_TYPE=0,			// 	(0)
	GAIN_PARAMETER,				// 	(1)
	GENOME_INITIALISATION,		//	(2)
	CHEMISTRY_SETUP,			//	(3)
	SNP_RATE,						//	(4)
	DELETION_PROB,				//	(5)
	DUPLICATION_PROB, 			//	(6)
	MAX_MUT_CHANGE,			//	(7)
	CHEMISTRY_MUT_PROB,		//	(8)
	INTERACTION_DEL_PROB,		//	(9)
	BASAL_ENERGY_EXPEND,		//	(10)
	MAX_MOVEMENT,				//	(11)
	MAX_INTAKE					//	(12)

};

// Verbosity options
enum 
{
	O_NONE 			= 0x000,		// silent mode
	O_RESULTS		= 0x001,		// global statistics
	O_POP_SNAPSHOT		= 0x002,		// instantaneous snapshot of location of all population members
	O_NUTRIENT_SNAPSHOT	= 0x004,		// instantaneous snapshot of food landscape
	O_TOXIN_SNAPSHOT 	= 0x008,		// instantaneous snapshot of toxin landscape
	O_SINGLE_IND_DATA	= 0x010,		// individual data for first member
	O_SINGLE_IND_LOCATION	= 0x020,		// behaviour data for first member
	O_SINGLE_IND_PROTEOME	= 0x040,		// proteome data for first member
	O_ALL_IND_DATA		= 0x080,		// individual data for all members
	O_ALL_IND_LOCATION	= 0x100,		// behavioural data for all members
	O_ALL_IND_PROTEOME	= 0x200			// proteome behaviour for all members
};


enum
{
	B_NONE = 0,	// 0
	B_DEATH,	// 1
	B_REPRODUCTION,	// 2
	B_V3,		// 3
	B_V4,		// 4
	B_V5,		// 5
	B_V6,		// 6
	B_V7,		// 7
	B_V8,		// 8
	B_V9,		// 9
	B_V10,		// 10
	B_V11,		// 11
	B_V12,		// 12
	B_V13,		// 13
	B_TUMBLE,	// 14
	B_ENERGY	// 15
};

enum
{

	UNIFORM_DET = 0,			// 	(0)	uniform deterministic deployment of chemical
	LINEAR_DET,				//	(1) 	linear increment in gradient
	GAUSSIAN_DET,				//	(2)
	UNIFORM_PROB				//	(3)

};

enum
{

	X_LOCATION = 0,			// 	(0)	x coordinate of individual
	Y_LOCATION,			//	(1) 	y coordinate of individual
	HEADING				//	(2)	heading angle (in radians)

};

#endif
