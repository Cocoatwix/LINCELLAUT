
/* A simple program for calculating various things about linear
 *  cellular automata.
 *
 * Mar 27, 2022
 *
 */
 
/* The following resources were used as a reference:
docs.microsoft.com/en-us/cpp/c-language/cpp-integer-limits?view=msvc-170
stackoverflow.com/questions/3219393

stackoverflow.com/questions/1190870
*/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h> //So I can get the maximum integer

#include <time.h> //For seeding random numbers 
#include <math.h>

#include "headers/helper.h"
#include "headers/bigint.h" //Arbitrary precision

#include "headers/linalg.h" 
#include "headers/cycles.h"  //Allows us to use Floyd's Algorithm
#include "headers/modular.h" //Modular square roots and inverses
#include "headers/factors.h" //For LCM()

#include "headers/algebra.h"

#define FREE(v) free(v); v = NULL

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define SET_BIG_NUM(str, bigint, msg) \
if (strtoBIT((str), &(bigint)) != 1) \
{ \
	fprintf(stderr, "%s\n", (msg)); \
	returnvalue = EXIT_FAILURE; \
	goto FREE_VARIABLES; \
}

#define ADD_FILE_PREFIX(filepath) \
if (haveOutputPrefix) \
	strcat(filepath, outputprefix)
	
#define PROHIBIT_UNIXFLAGS \
if (unixflagsCount > 0) \
{ \
	fprintf(stderr, "UNIX-like flags are currently not supported with this tool.\n"); \
	returnvalue = EXIT_FAILURE; \
	goto FREE_VARIABLES; \
}


int main(int argc, char* argv[])
{
	srand(time(NULL)); //Seeding random number generator
	
	//Check to see if MAXBUNCH for BigIntT structs is a power of ten
	//If it isn't, warn the user
	for (int i = 10; i < MAXBUNCH; i *= 10)
		if (MAXBUNCH % i != 0)
		{
			printf("MAXBUNCH is not a power of ten. Some functionality involving the use" \
			" of BigIntT structs may not work correctly. Change the MAXBUNCH constant within" \
			" bigint.c to fix this issue.\n");
			break;
		}
	
	//Used to tell the program what to return when jumping to FREE_VARIABLES
	int returnvalue = EXIT_SUCCESS;
	
	//The maximum length for a string in the .config file
	const int MAXSTRLEN = 101;
	
	VectorTypeT vectorType = row;
	
	//Whether to suppress certain output (decided via "quiet" key in .config file)
	bool quietmode = FALSE;
	
	int iterations;
	int modulus;
	
	//For holding str representation of mod, used for BigIntT number (if needed)
	char* bigintmodstring = malloc(MAXSTRLEN*sizeof(char));
	char* resumemodstring = malloc(MAXSTRLEN*sizeof(char));
	
	char* updatefilepath     = malloc(MAXSTRLEN*sizeof(char));
	char* initialfilepath    = malloc(MAXSTRLEN*sizeof(char));
	char* resumefilepath     = malloc(MAXSTRLEN*sizeof(char));
	char* sentinelfilepath   = malloc(MAXSTRLEN*sizeof(char));
	char* polynomialfilepath = malloc(MAXSTRLEN*sizeof(char));
	char* iterfilepath       = malloc(MAXSTRLEN*sizeof(char));
	
	bool haveOutputPrefix = FALSE; //Says whether we have a prefix to add to any output file paths
	char* outputprefix = malloc(MAXSTRLEN*sizeof(char));
	outputprefix[0] = '\0';
	
	//Creating some default BigIntTPs that are used often
	int oneArr[1] = {1};
	BigIntTP zero = empty_BigIntT(1);
	BigIntTP one  = new_BigIntT(oneArr, 1);
	
	BigIntMatrixTP UPDATEMATRIX  = NULL;
	BigIntMatrixTP INITIALMATRIX = NULL;
	
	//Temporarily holds config data
	char* systemData  = malloc(MAXSTRLEN*sizeof(char));
	
	char* tempStr; //Holds temporary info when using string functions
	
	//Read .config file to get appropriate data loaded
	FILE* system = fopen("config/system.config", "r");
	
	if (system == NULL)
	{
		fprintf(stderr, "Unable to open .config file. Check the config folder to ensure the file exists!\n");
		returnvalue = EXIT_FAILURE;
		goto FREE_VARIABLES;
	}
	
	while (fscanf(system, "%10s", systemData) == 1)
	{
		//printf("%s\n", systemData);
		if (! strcmp(systemData, "mod"))
		{
			//Getting modulus as string first in case we need to convert to BigIntT
			if (fscanf(system, "%101s", bigintmodstring) != 1)
			{
				fprintf(stderr, "Unable to read modulus from .config file. " \
				"Defaulting to a modulus of 3.\n");
				bigintmodstring[0] = '3';
				bigintmodstring[1] = '\0';
				modulus = 3;
				continue;
			}
			
			if (strlen(bigintmodstring) < 10)
			{
				modulus = (int)strtol(bigintmodstring, &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Invalid modulus provided in .config file. " \
					"Defaulting to a modulus of 3.\n");
					modulus = 3;
				}
			}
		}
		
		//If the user wants to change how vectors are displayed to the console
		else if (! strcmp(systemData, "vectorType"))
		{
			if (fscanf(system, "%s", systemData) != 1)
			{
				fprintf(stderr, "Unable to read vector type from .config file. " \
				"Defaulting to \"row\" vector type.\n");
			}
			
			else
				if (!strcmp(systemData, "col"))
					vectorType = col;
		}
		
		else if (! strcmp(systemData, "iterations"))
		{
			if (fscanf(system, "%d", &iterations) != 1)
			{
				fprintf(stderr, "Unable to read number of iterations from .config file. " \
				"Defaulting to 1 iteration.\n");
				iterations = 1;
			}
		}
		
		else if (! strcmp(systemData, "update"))
		{
			if (fscanf(system, "%s", updatefilepath) != 1)
			{
				fprintf(stderr, "Unable to read update matrix path from .config file. " \
				"Continuing without setting an update matrix...\n");
				updatefilepath = NULL;
				UPDATEMATRIX = NULL;
				continue;
			}
			
			UPDATEMATRIX = read_BigIntMatrixT(updatefilepath);
			if (UPDATEMATRIX == NULL)
			{
				fprintf(stderr, "Unable to set update matrix from .config file. " \
				"Continuing without setting an update matrix...\n");
			}
		}
		
		else if (! strcmp(systemData, "initial"))
		{
			if (fscanf(system, "%s", initialfilepath) != 1)
			{
				fprintf(stderr, "Unable to read initial matrix path from .config file. " \
				"Continuing without setting an initial matrix...\n");
				initialfilepath = NULL;
				INITIALMATRIX = NULL;
				continue;
			}
			
			INITIALMATRIX = read_BigIntMatrixT(initialfilepath);
			if (INITIALMATRIX == NULL)
			{
				fprintf(stderr, "Unable to set initial matrix from .config file. " \
				"Continuing without setting an initial matrix...\n");
			}
		}
		
		else if (! strcmp(systemData, "resumeMat"))
		{
			if (fscanf(system, "%s", resumefilepath) != 1)
			{
				fprintf(stderr, "Unable to read resume matrix path from .config file. " \
				"Continuing without setting a resume matrix...\n");
				resumefilepath = NULL;
			}
		}
		
		else if (! strcmp(systemData, "sentinel"))
		{
			if (fscanf(system, "%s", sentinelfilepath) != 1)
			{
				fprintf(stderr, "Unable to read sentinel matrix path from .config file. " \
				"Continuing without setting a sentinel matrix...\n");
				sentinelfilepath = NULL;
			}
		}
		
		else if (! strcmp(systemData, "polynomial"))
		{
			if (fscanf(system, "%s", polynomialfilepath) != 1)
			{
				fprintf(stderr, "Unable to read polynomial path from .config file. " \
				"Continuing without setting a polynomial...\n");
				polynomialfilepath = NULL;
			}
		}
		
		else if (! strcmp(systemData, "resumeMod"))
		{
			if (fscanf(system, "%101s", resumemodstring) != 1)
			{
				fprintf(stderr, "Unable to read resume modulus from .config file. " \
				"Defaulting to a resume modulus of 3.\n");
				resumemodstring[0] = '3';
				resumemodstring[1] = '\0';
			}
		}
		
		else if (! strcmp(systemData, "itername"))
		{
			if (fscanf(system, "%101s", iterfilepath) != 1)
			{
				fprintf(stderr, "Unable to read iteration file path from .config file. " \
				"Continuing without setting an iteration file path...\n");
				iterfilepath = NULL;
			}
		}
		
		else if (!strcmp(systemData, "quiet"))
		{
			if (fscanf(system, "%101s", systemData) != 1)
			{
				fprintf(stderr, "Unable to read value of \"quiet\" key in .config file. " \
				"Quiet mode will be turned off.\n");
			}
			
			else
				if (!strcmp(systemData, "TRUE"))
					quietmode = TRUE;
		}
		
		else if (!strcmp(systemData, "outpath"))
		{
			if (fscanf(system, "%101s", systemData) != 1)
			{
				fprintf(stderr, "Unable to read default output file path from .config file. " \
				"Output files will be created in the working directory.\n");
			}
			else
			{
				haveOutputPrefix = TRUE;
				strcpy(outputprefix, systemData);
				if (outputprefix[strlen(outputprefix)-1] != '/')
					strcat(outputprefix, "/");
			}
		}
	}
	
	FREE(systemData);
	if (fclose(system) == EOF)
	{
		fprintf(stderr, "Unable to properly close .config file.\n");
		returnvalue = EXIT_FAILURE;
		goto FREE_VARIABLES;
	}
	
	/* $$$ */
	
	//This holds the names of all the UNIX-type-compatiable CLI tools
	const int UNIXTYPECOUNT = 24;
	
	const char* CLINAMES[] = {"iterate",    "inverse",          "rowreduce", "factor",        "order",
	                          "evalpoly",   "chara",            "orbits",    "splitorbits",   "orbitreps",
														"branchreps", "orbitspaces",      "floyd",     "cycconvmat",    "ccmzerosearch",
														"vectprops",  "vectpolys",        "matprops",  "fibmultsearch", "dynamics",
														"orbitmaps2", "oddminpolysearch", "stablelift", "allstablelifts"};

	//Lists of argument names for the CLI tools above
	//The empty strings are so that I don't have to allocate memory for
	// an irregular array shape.
	const char* ARGNAMES[][6] = {{"iterations", "", "", "", "", ""}, //iterate
		                           {"modulus", "", "", "", "", ""},    //inverse
															 {"modulus", "", "", "", "", ""},    //rowreduce
															 {"modulus", "", "", "", "", ""},    //factor
															 {"modulus", "", "", "", "", ""},    //order
															 {"modulus", "multiplyByInitial", "", "", "", ""},  //evalpoly
															 {"modulus", "", "", "", "", ""},                   //chara
															 {"modulus", "fileoutput", "", "", "", ""},         //orbits
															 {"modulus", "fileoutput", "", "", "", ""},         //splitorbits
															 {"modulus", "fileoutput", "", "", "", ""},         //orbitreps
															 {"modulus", "", "", "", "", ""},                   //branchreps
															 {"modulus", "minpolys", "fileoutput", "", "", ""}, //orbitspaces
															 {"modulus", "", "", "", "", ""},                   //floyd
															 {"from", "to", "mod", "", "", ""},                 //cycconvmat
															 {"resume", "size", "mod", "", "", ""},             //ccmzerosearch
															 {"baseMod", "modPower", "resume", "fileoutput", "", ""},     //vectprops
															 {"baseMod", "modPower", "", "", "", ""},                     //vectpolys
															 {"maxpower", "modulus", "", "", "", ""},                     //matprops
															 {"bound", "", "", "", "", ""},                               //fibmultsearch
															 {"maxPower", "modulus", "allConfigs", "fileoutput", "", ""}, //dynamics
		                           {"maxpower", "modulus", "fileoutput", "belowBound", "aboveBound", ""}, //orbitmaps2
															 {"maxpower", "size", "polysize", "modulus", "resume", "fileoutput"},   //oddminpolysearch 
															 {"baseMod", "modPower", "", "", "", ""},                               //stablelift
															 {"baseMod", "modPower", "", "", "", ""}};                              //allstablelifts

	//The maximum number of command line arguments each tool takes, +1
	const int MAXARGC[] = {2, 2, 2, 2, 2, 
	                       3, 2, 3, 3, 3,
												 2, 4, 2, 4, 4,
												 5, 3, 3, 2, 5,
												 6, 7, 3, 3};

	//Will hold key-value pairs if user prefers unix flags
	//Size of unixflags will be argc/2
	//Each key-value pair contains two strings
	char*** unixflags = NULL;
	char** oargv = NULL; //Ordered CLI arguments
	char** kvpair; //Points to a key-value pair; basically a temp variable
	int maxargc = 0; //Max number of CLI arguments our selected tool uses
	int toolindex = -1; //Where in the UNIX-type data our tool's data is stored
	
	//If we have command line arguments
	if (argc > 1)
	{
		/* Preparing UNIX-type arguments, if possible. */
		
		//First, create an array of any unix-type flags that were passed via the CLI
		unixflags = malloc((argc/2)*sizeof(char**));
		int unixflagsCount = 0;
		
		for (int i = 0; i < argc/2; i += 1)
		{
			unixflags[i] = malloc(2*sizeof(char*));
			unixflags[i][0] = malloc(MAXSTRLEN*sizeof(char));
			unixflags[i][1] = malloc(MAXSTRLEN*sizeof(char));
			unixflags[i][0][0] = '\0';
			unixflags[i][1][0] = '\0';
		}
		
		for (int i = 2; i < argc; i += 1)
		{
			if (argv[i][0] == '-') //Found a unix flag
			{
				//First, check to make sure that the next argument is even valid...
				// or if we even have another argument...
				if ((i+1 >= argc) || (argv[i+1][0] == '-'))
				{
					fprintf(stderr, "Argument flag '%s' has no given value.\n", argv[i]+1);
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
				
				strcpy(unixflags[unixflagsCount][0], argv[i]+1);
				i += 1;
				strcpy(unixflags[unixflagsCount][1], argv[i]);
				unixflagsCount += 1;
			}
		}
		
		//Let's see how the unix flags look
		/*
		printf("Flags:\n");
		for (int i = 0; i < argc/2; i += 1)
		{
			if (unixflags[i][0][0] != '\0')
				printf("%s\t%s\n", unixflags[i][0], unixflags[i][1]);
		}
		printf("\n");
		*/
		
		//Get the index of our tool's UNIX-type data
		for (int i = 0; i < UNIXTYPECOUNT; i += 1)
		{
			//If we find a tool with the given name that has UNIX-type arguments
			if (!strcmp(argv[1], CLINAMES[i]))
			{
				toolindex = i;
				
				oargv = malloc((MAXARGC[toolindex]+1)*sizeof(char*)); //Ordered arguments
				for (int i = 0; i <= MAXARGC[toolindex]; i += 1)
				{
					oargv[i] = malloc(MAXSTRLEN*sizeof(char));
					oargv[i][0] = '\0';
				}
				
				//The duplicates
				strcpy(oargv[0], "./lincellaut");
				strcpy(oargv[1], argv[1]);
				
				//Prep oargv for use below
				if (unixflagsCount > 0)
				{
					for (int key = 0; key < unixflagsCount; key += 1)
					{
						kvpair = unixflags[key]; //Our current key-value pair
						
						//Compare our key to the arguments our tool takes
						for (int toolArg = 0; toolArg < MAXARGC[toolindex]-1; toolArg += 1)
							if (! strcmp(kvpair[0], ARGNAMES[toolindex][toolArg]))
							{
								strcpy(oargv[toolArg+2], kvpair[1]);
								break;
							}
					}
				}
				
				//Positional arguments
				else
					for (int val = 2; val < argc; val += 1)
						strcpy(oargv[val], argv[val]);
				
				break;
			}
		}
		
		//If we didn't find a UNIX-type argument tool, or there were no
		// UNIX-type arguments passed
		if (toolindex == -1)
		{
			oargv = malloc(argc*sizeof(char*)); //Ordered arguments
			for (int i = 0; i < argc; i += 1)
			{
				oargv[i] = malloc(MAXSTRLEN*sizeof(char));
				oargv[i][0] = '\0';
			}
			
			//The duplicates
			strcpy(oargv[0], "./lincellaut");
			strcpy(oargv[1], argv[1]);
			
			//Defaulting to regular positional arguments
			for (int i = 2; i < argc; i += 1)
				strcpy(oargv[i], argv[i]);
		}
		
		/*
		printf("oargv:\n");
		for (int t = 0; t <= MAXARGC[toolindex]; t += 1)
			printf("%s\n", oargv[t]);
		*/
		
		/* -- COMMAND LINE TOOLS -- */
		
		//If we want to iterate the given update matrix
		if (! strcmp(argv[1], "iterate"))
		{
			BigIntMatrixTP F;
			BigIntMatrixTP F_2;
			BigIntMatrixTP F_result = NULL;
			BigIntMatrixTP power = NULL;
			
			int tempArr[1];
			BigIntTP bigIters = NULL;
			
			BigIntTP bigMod;
			SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
		
			//If the user provided a custom number of iterations
			if (oargv[2][0] != '\0')
			{
				iterations = (int)strtol(oargv[2], &tempStr, 10);
				
				//If we didn't read any digits for iterations
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Invalid number of iterations provided at command line.\n");
					bigMod = free_BigIntT(bigMod);
					
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
			}
			
			printf("Iterations: %d\n", iterations);
			if (iterations < 0)
			{
				F_2 = read_BigIntMatrixT(updatefilepath);
				F   = big_inverse(F_2, bigMod);
				
				if (F == NULL)
				{
					printf("An inverse for the update matrix modulo ");
					printi(bigMod);
					printf(" doesn't exist.\n");
				}
				
				else
				{
					iterations *= -1;
					F_2 = free_BigIntMatrixT(F_2);
					F_2 = read_BigIntMatrixT(initialfilepath);
				}
			}
			
			else
			{
				F   = read_BigIntMatrixT(updatefilepath);
				F_2 = read_BigIntMatrixT(initialfilepath);
			}
			
			if ((F == NULL) || (F_2 == NULL))
			{
				fprintf(stderr, "Unable to read matrices from config file.\n");
				F = free_BigIntMatrixT(F);
				F_2 = free_BigIntMatrixT(F_2);
				bigMod = free_BigIntT(bigMod);
				
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (big_cols(F) != big_rows(F_2))
			{
				fprintf(stderr, "Given matrices cannot be multiplied.\n");
				F = free_BigIntMatrixT(F);
				F_2 = free_BigIntMatrixT(F_2);
				bigMod = free_BigIntT(bigMod);
				
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Iterate F a few times
			//The minus 1 is for easier conversion between ORBITVIS results
			if (iterations == 1)
			{
				F_result = new_BigIntMatrixT(big_rows(F), big_cols(F_2));
				big_mat_mul(F, F_2, F_result);
				modbm(F_result, bigMod);
			}
			
			else if (iterations > 0)
			{
				tempArr[0] = iterations;
				bigIters = new_BigIntT(tempArr, 1);
				F_result = new_BigIntMatrixT(big_rows(F_2), big_cols(F_2));
				
				power = new_BigIntMatrixT(big_rows(F), big_cols(F));
				powbm(F, power, bigIters, bigMod);
				big_mat_mul(power, F_2, F_result);
				modbm(F_result, bigMod);
			}
			
			//Prevents the matrix from being printed when
			// an inverse doesn't exist and the iterations
			// provided was negative
			if (iterations > 0)
			{
				if (vectorType == row)
					printbm_row(F_result);
				else
					printbm(F_result);
			}
			
			else if (iterations == 0)
			{
				if (vectorType == row)
					printbm_row(F_2);
				else
					printbm(F_2);
			}
			
			if ((vectorType == row) && (big_rows(F_2) != big_cols(F_2)))
				printf("\n");

			F        = free_BigIntMatrixT(F);
			F_2      = free_BigIntMatrixT(F_2);
			F_result = free_BigIntMatrixT(F_result);
			power    = free_BigIntMatrixT(power);
			
			bigIters = free_BigIntT(bigIters);
		}
		
		
		//Find the inverse of the update matrix
		else if (!strcmp(argv[1], "inverse"))
		{
			BigIntMatrixTP F = UPDATEMATRIX;
			BigIntMatrixTP Finv;
			BigIntTP bigMod;
			
			if (F == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from config file.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			printf("Update matrix:\n");
			printbm(F);
			
			//If the user specified a modulus at the command line
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from .config file.");
			}
			
			Finv = big_inverse(F, bigMod);
			
			if (Finv == NULL)
			{
				printf("An inverse for the update matrix modulo ");
				printi(bigMod);
				printf(" doesn't exist.\n");
			}
			
			else
			{
				printf("Inverse modulo ");
				printi(bigMod);
				printf(":\n");
				printbm(Finv);
			}
			
			Finv = free_BigIntMatrixT(Finv);
			bigMod = free_BigIntT(bigMod);
		}
		
		
		//If we want to row reduce the update matrix
		// by the given modulus.
		//./lincellaut rowreduce [modulus]
		else if (!strcmp(argv[1], "rowreduce"))
		{
			BigIntTP bigMod;
			
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from .config file.");
			}
			
			BigIntMatrixTP reducedMatrix = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_cols(UPDATEMATRIX));
			
			printf("Note that this tool is meant to be used only with prime and prime-power moduli. Results for general composite moduli are undefined.\n");
			
			big_row_reduce(UPDATEMATRIX, bigMod, reducedMatrix);
			printf("Modulus: ");
			printi(bigMod);
			printf("\nMatrix:\n");
			printbm(UPDATEMATRIX);
			printf("\nRow-reduced matrix:\n");
			printbm(reducedMatrix);
			printf("\n");
			
			bigMod = free_BigIntT(bigMod);
			reducedMatrix = free_BigIntMatrixT(reducedMatrix);
		}
		
		
		//Factor the given polynomial from the .config file
		else if (!strcmp(argv[1], "factor"))
		{
			BigIntTP bigMod;
			
			//If the user specified a modulus at the command line
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from .config file.");
			}
			
			BigPolyTP P = read_BigPolyT(polynomialfilepath);
			if (P == NULL)
			{
				fprintf(stderr, "Unable to read polynomial at %s.\n", polynomialfilepath);
				bigMod = free_BigIntT(bigMod);
				
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Now, let's factor the polynomial
			BigFactorsTP factoredP = factor_BigPolyT(P, bigMod);
			printf("Original: ");
			printp(P);
			printf("\nFactored mod ");
			printi(bigMod);
			printf(": ");
			printpf(factoredP);
			printf("\n");
			
			P = free_BigPolyT(P);
			factoredP = free_BigFactorsT(factoredP);
			bigMod = free_BigIntT(bigMod);
		}
		
		
		//Find the order of the given polynomial
		else if (!strcmp(argv[1], "order"))
		{
			BigIntTP bigMod, negOne, temp, temp2, bigOrder;
			BigIntTP orderMustDivide; //The number which our extensions' orders must divide
			BigIntTP* lambdaCoeffs;
			
			BigPolyTP P, unityPoly, negOnePoly, zeroPoly, lambdaPoly;
			BigPolyTP tempPoly, quotientPoly, remainderPoly;
			
			BigFactorsTP factoredP;
			BigPolyTP* coprimeFactors;
			int* mults; //Holds the multiplicities of the factors in coprimeFactors
			
			//For making sure we check the correct possible order for multiple roots
			BigIntTP* primePowersForOrder;
			
			int tArray[1] = {1};        //For initialising value of iterationExt
			BigIntTP* polyCoeffs;       //Holds the coeffs for a MultiVarExtT's min poly
			MultiVarExtTP oneExt;       //Holds the number 1 for comparison purposes
			MultiVarExtTP tempExt;      //Holds temporary results of computations
			MultiVarExtTP factorExt;    //Extension representing a factor of the given polynomial
			MultiVarExtTP iterationExt; //Extension for iterating factorExt (holds 1*<extension>)
			
			//These are for getting all the possible orders
			// for a given extension
			typedef struct ordstruct
			{
				BigIntTP possibleOrder;
				BigIntTP* factors;
				int numOfFactors;
			} OrderStructT, *OrderStructTP;
			
			typedef struct linkedlist
			{
				OrderStructTP anOrder;
				void* next;
			} OrderListT, *OrderListTP;
			
			int positionOfInterest;
			int* factorPositions; //This will help up create divisors of orderMustDivide
			BigIntTP* primeFactors;
			
			OrderStructTP temporaryOrder;
			OrderListTP possibleOrders;
			
			//For helping to navigate the linked list
			OrderListTP prevListPointer;
			OrderListTP listPointer;
			
			//If the user specified a modulus at the command line
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from .config file.");
			}
			
			P = read_BigPolyT(polynomialfilepath);
			if (P == NULL)
			{
				fprintf(stderr, "Unable to read polynomial at %s.\n", polynomialfilepath);
				bigMod = free_BigIntT(bigMod);
				
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			bigOrder = empty_BigIntT(1);
			copy_BigIntT(one, bigOrder);
			
			//Initialise helper variables
			temp = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			orderMustDivide = empty_BigIntT(1);
			tempPoly = empty_BigPolyT();
			zeroPoly = constant_BigPolyT(zero);
			quotientPoly = empty_BigPolyT();
			remainderPoly = empty_BigPolyT();
			
			//These are for creating polynomials of the form x^c - 1
			negOne = empty_BigIntT(1);
			subtract_BigIntT(bigMod, one, negOne);
			negOnePoly = constant_BigPolyT(negOne);
			
			lambdaCoeffs = malloc(2*sizeof(BigIntTP));
			lambdaCoeffs[0] = zero;
			lambdaCoeffs[1] = one;
			lambdaPoly = new_BigPolyT(lambdaCoeffs, 2);
			
			//Create x^1 - c
			lambdaCoeffs[0] = negOne;
			unityPoly = new_BigPolyT(lambdaCoeffs, 2);
			
			//The order of a polynomial will be the lcm of the multiplicative orders
			// of the algebraic objects defined by its factors... I think.
			//So, let's factor the polynomial and get the multiplicative orders.
			factoredP = factor_BigPolyT(P, bigMod);
			mults = multiplicities(factoredP);
			coprimeFactors = extract_coprime_factors(factoredP);
			
			primePowersForOrder = malloc(count_unique_factors(factoredP)*sizeof(BigIntTP));

			printf("P: ");
			printp(P);
			printf("\nFactored P: ");
			printpf(factoredP);
			printf("\nCoprime factors of P: ");
			for (int i = 0; i < count_unique_factors(factoredP); i += 1)
			{
				if (i != 0)
					printf(", ");
				printp(coprimeFactors[i]);
				
				//Initialising each element to one
				primePowersForOrder[i] = new_BigIntT(tArray, 1);
				
				//THIS WILL NOT WORK FOR VERY LARGE POLYNOMIALS
				clear_BigIntT(temp);
				set_bunch(temp, 0, mults[i]);
				
				//Find the first power of p that's greater than the
				// multiplicity of this particular factor
				while (compare_BigIntT(temp, primePowersForOrder[i]) > 0)
				{
					multiply_BigIntT(primePowersForOrder[i], bigMod, temp2);
					copy_BigIntT(temp2, primePowersForOrder[i]);
				}
				
			}
			printf("\n");
			
			//Now, we need to create MultiVarExtTs for each coprimeFactor
			// and find out what their multiplicative orders are
			for (int f = 0; f < count_unique_factors(factoredP); f += 1)
			{
				//Making sure we're not considering factors of the form x^k
				divide_BigPolyT(coprimeFactors[f], lambdaPoly, quotientPoly, remainderPoly, bigMod);
				if (compare_BigPolyT(remainderPoly, zeroPoly) != 0)
				{
					//Creating the field extension defined by the given factor of
					// our polynomial so we can find its order
					oneExt       = new_MultiVarExtT(1);
					tempExt      = new_MultiVarExtT(1);
					factorExt    = new_MultiVarExtT(1);
					iterationExt = new_MultiVarExtT(1);
					set_MultiVarExtT_mod(oneExt, bigMod);
					set_MultiVarExtT_mod(factorExt, bigMod);
					set_MultiVarExtT_mod(iterationExt, bigMod);
					
					polyCoeffs = extract_coefficients(coprimeFactors[f]);
					add_extension(oneExt, polyCoeffs, degree(coprimeFactors[f])+1, "t");
					add_extension(tempExt, polyCoeffs, degree(coprimeFactors[f])+1, "t");
					add_extension(factorExt, polyCoeffs, degree(coprimeFactors[f])+1, "t");
					add_extension(iterationExt, polyCoeffs, degree(coprimeFactors[f])+1, "t");
					
					for (int i = 0; i < degree(coprimeFactors[f])+1; i += 1)
						polyCoeffs[i] = free_BigIntT(polyCoeffs[i]);
					FREE(polyCoeffs);
					
					//Setting the values of our extensions to prepare for iteration
					increment_MultiVarExtT(oneExt);
					increment_MultiVarExtT(factorExt);
					set_MultiVarExtT_coefficient(iterationExt, tArray, one);
					
					//Now, compute the value which any extension's order must divide (p^d - 1)
					//THIS WILL NOT WORK FOR EXTREMELY HIGH VALUES
					set_bunch(temp2, 0, degree(coprimeFactors[f]));
					pow_BigIntT(bigMod, temp2, temp);
					subtract_BigIntT(temp, one, orderMustDivide);
					
					/*
					REMOVE THIS ONCE TESTING IS DONE
					resize_BigIntT(orderMustDivide, 1);
					set_bunch(orderMustDivide, 0, 360);
					*/
					
					//So, whatever the order of our extension is, it must divide orderMustDivide
					//Let's go through all the possible divisors and find the one that works
					primeFactors = prime_factors_of_BigIntT(orderMustDivide);
					
					
					#ifdef VERBOSE
						printi(bigMod);
						printf("^%d - 1 = ", degree(coprimeFactors[f]));
						printi(orderMustDivide);
						
						printf("\nPrime factors: ");
						for (int i = 1; i <= extract_bunch(primeFactors[0], 0); i += 1)
						{
							if (i != 1)
								printf(", ");
							printi(primeFactors[i]);
						}
						printf("\n");
					#endif
					
					possibleOrders = calloc(1, sizeof(OrderListT));
					
					//Now, we must create every possible divisor of orderMustDivide 
					// and sort them.
					factorPositions = calloc(extract_bunch(primeFactors[0], 0), sizeof(int));
					for (int numOfFactors = 1; numOfFactors <= extract_bunch(primeFactors[0], 0); numOfFactors += 1)
					{
						//Initialise factorPositions to get ready for next round of factors
						for (int i = 0; i < numOfFactors; i += 1)
							factorPositions[i] = i;
						
						//Create all possible factors using the given number of factors
						while (factorPositions[0] <= extract_bunch(primeFactors[0], 0) - numOfFactors)
						{
							//Create OrderStructTs here
							temporaryOrder = malloc(sizeof(OrderStructT));
							temporaryOrder->numOfFactors = numOfFactors;
							temporaryOrder->factors = malloc(numOfFactors*sizeof(BigIntTP));
							temporaryOrder->possibleOrder = empty_BigIntT(1);
							
							//Printing out the current factorisation of our factor for testing
							//Also storing the factorisation within our temporary OrderStructT
							//ALSO multiplying the facorisation togetehr to get the possibleOrder
							copy_BigIntT(one, temp);
							for (int i = 0; i < numOfFactors; i += 1)
							{
								/*
								if (i != 0)
									printf("*");
								printi(primeFactors[factorPositions[i]+1]);
								*/
								
								temporaryOrder->factors[i] = primeFactors[factorPositions[i]+1];
								
								multiply_BigIntT(temp, primeFactors[factorPositions[i]+1], temp2);
								copy_BigIntT(temp2, temp);
							}
							
							copy_BigIntT(temp, temporaryOrder->possibleOrder);
							
							//Now, we have to insert temporaryOrder in its correct position
							// in the sorted list of OrderStructTs.
							//I'm just using Bubble Sort. Sue me.
							listPointer = possibleOrders;
							prevListPointer = NULL;
							if (possibleOrders->anOrder == NULL) //If we haven't started the list yet
							{
								possibleOrders->anOrder = temporaryOrder;
							}
								
							else
							{
								while (TRUE)
								{
									//Check to see if current order is greater than the one at listPointer
									if (compare_BigIntT(temporaryOrder->possibleOrder, listPointer->anOrder->possibleOrder) > 0)
									{
										//Add temporaryOrder to the end of the list
										if (listPointer->next == NULL)
										{
											listPointer->next = calloc(1, sizeof(OrderListT));
											((OrderListTP)(listPointer->next))->anOrder = temporaryOrder;
											break;
										}
										
										//Navigate further through the list
										else
										{
											prevListPointer = listPointer;
											listPointer = listPointer->next;
										}
									}
									
									//If our current order is less than the one at listPointer
									else
									{
										//if there is no previous element to place our current order in front of
										if (prevListPointer == NULL)
										{
											prevListPointer = malloc(sizeof(OrderListT));
											prevListPointer->anOrder = temporaryOrder;
											prevListPointer->next = listPointer;
											possibleOrders = prevListPointer;
										}
										
										else
										{
											prevListPointer->next = malloc(sizeof(OrderListT));
											((OrderListTP)(prevListPointer->next))->anOrder = temporaryOrder;
											((OrderListTP)(prevListPointer->next))->next = listPointer;
										}
										
										break;
									}
								}
							}
							
							//Now, for my own sanity, let's print out the linked list each pass to
							// make sure it's ordering correctly
							#ifdef VERBOSE
								listPointer = possibleOrders;
								printf("[");
								while (listPointer != NULL)
								{
									printi(listPointer->anOrder->possibleOrder);
									printf(",");
									
									listPointer = listPointer->next;
								}
								printf("]\n");
							#endif
							
							//Increment factorPositions to next unique factor
							positionOfInterest = numOfFactors - 1;
							while (positionOfInterest >= 0)
							{
								#ifdef VERBOSE
									//Print some stuff out here, showing what's actually going on under the hood
									/*
									printf("[");
									for (int i = 0; i < numOfFactors; i += 1)
									{
										if (i != 0)
											printf(", ");
										printf("%d", factorPositions[i]);
									}
									printf("]\n");
									getchar();
									*/
								#endif
								
								factorPositions[positionOfInterest] += 1;
								
								//If our positionOfInterest goes too far
								if (factorPositions[positionOfInterest] > extract_bunch(primeFactors[0], 0) - (numOfFactors - positionOfInterest))
								{
									//Push all relevant positions forward... except if our position of interest is 0.
									// In that case, we don't need to push anything forward since we're done
									if (positionOfInterest > 0)
										for (int i = positionOfInterest; i < numOfFactors; i += 1)
											factorPositions[i] = factorPositions[positionOfInterest-1] + 2 + (i - positionOfInterest);
									
									//Roll over to next positionOfInterest
									positionOfInterest -= 1;
								}
								
								//If we increment a position to the same prime factor it was pointing to before
								else if (compare_BigIntT(primeFactors[factorPositions[positionOfInterest]+1],
								                         primeFactors[factorPositions[positionOfInterest]]) == 0)
							  {
									//Increment current position again
									//Push all relevant positions forward
									for (int i = positionOfInterest + 1; i < numOfFactors; i += 1)
										factorPositions[i] = factorPositions[positionOfInterest] + 1 + (i - positionOfInterest);
							  }
								
								//If we have a good set of factorPositions that'll give us a new factor
								//It should always be a new factor, so we don't need to check for duplicates
								else
									break;
							}
						}
					}
					FREE(factorPositions);
					
					//NOW, we need to iterate through our entire list of possible orders
					// and find the first one where our extension becomes 1 when raised to
					// that power.
					
					//Before we check all the prime factors, we should check whether 1
					// is a valid order for the extension.
					mult_sim_MultiVarExtT(factorExt, iterationExt, tempExt);
					reduce_MultiVarExtT(tempExt);
					copy_MultiVarExtT(tempExt, factorExt);
					
					//If the order is NOT one
					if (compare_MultiVarExtT(factorExt, oneExt) != 1)
					{
						copy_MultiVarExtT(oneExt, factorExt);
						listPointer = possibleOrders;
						
						while (listPointer != NULL)
						{
							//Go through all prime factors of our order and iterate
							// the extension accordingly
							for (int i = 0; i < listPointer->anOrder->numOfFactors; i += 1)
							{
								//Add the factor of p we need to check "embedded orders"
								if (i == 0)
									multiply_BigIntT(listPointer->anOrder->factors[i], primePowersForOrder[f], temp);
								else
									copy_BigIntT(listPointer->anOrder->factors[i], temp);
								
								if (i != 0)
								{
									//This makes sure we don't raise our extension to
									// too high a power.
									subtract_BigIntT(temp, one, temp2);
									copy_BigIntT(temp2, temp);
								}
								
								while (! is_zero(temp))
								{
									//Do multiplication/iteration
									mult_sim_MultiVarExtT(factorExt, iterationExt, tempExt);
									
									reduce_MultiVarExtT(tempExt);
									copy_MultiVarExtT(tempExt, factorExt);
									
									subtract_BigIntT(temp, one, temp2);
									copy_BigIntT(temp2, temp);
								}
								
								copy_MultiVarExtT(factorExt, iterationExt);
							}
							
							/*
							printf("extension iterated to power ");
							printi(listPointer->anOrder->possibleOrder);
							printf(" * ");
							printi(primePowersForOrder[f]);
							printf(": ");
							printmve_row(factorExt);
							printf("\n");
							*/
							
							//Check to see if the extension is one.
							//If it is, we've found the order
							if (compare_MultiVarExtT(factorExt, oneExt) == 1)
							{
								multiply_BigIntT(listPointer->anOrder->possibleOrder, primePowersForOrder[f], temp2);
								big_lcm(bigOrder, temp2, temp);
								copy_BigIntT(temp, bigOrder);
								
								break;
							}
							
							else
							{
								copy_MultiVarExtT(oneExt, factorExt);
								clear_MultiVarExtT(iterationExt);
								set_MultiVarExtT_coefficient(iterationExt, tArray, one);
								
								//Try the next order
								listPointer = listPointer->next;
							}
						}
						
						//Build up order of the entire polynomial using LCM
						printf("Order of ");
						printp(coprimeFactors[f]);
						printf(" mod ");
						printi(bigMod);
						printf(": ");
						
						//temp2 = listPointer->anOrder->possibleOrder * primePowersForOrder[f],
						// so temp2 is some factor of p^d - 1, multiplied by some number of
						// factors of p.
						printi(temp2);
						
						printf("\n");
					}
					
					else
					{
						//Build up order of the entire polynomial using LCM
						printf("Order of ");
						printp(coprimeFactors[f]);
						printf(" mod ");
						printi(bigMod);
						printf(": 1\n");
					}
					
					oneExt = free_MultiVarExtT(oneExt);
					tempExt = free_MultiVarExtT(tempExt);
					factorExt = free_MultiVarExtT(factorExt);
					iterationExt = free_MultiVarExtT(iterationExt);
					
					for (int i = 1; i <= extract_bunch(primeFactors[0], 0); i += 1)
						primeFactors[i] = free_BigIntT(primeFactors[i]);
					primeFactors[0] = free_BigIntT(primeFactors[0]);
					FREE(primeFactors);
					
					//What a pain in the ass
					prevListPointer = NULL;
					listPointer = possibleOrders;
					do
					{
						prevListPointer = listPointer;
						listPointer = listPointer->next;
						
						free_BigIntT(prevListPointer->anOrder->possibleOrder);
						FREE(prevListPointer->anOrder->factors);
						FREE(prevListPointer->anOrder);
						FREE(prevListPointer);
					} 
					while (listPointer != NULL);
				}
			}
			
			printf("Order of ");
			printp(P);
			printf(": ");
			printi(bigOrder);
			printf("\n");
			
			bigOrder = free_BigIntT(bigOrder);
			
			FREE(lambdaCoeffs);
			
			temp = free_BigIntT(temp);
			temp2 = free_BigIntT(temp2);
			negOne = free_BigIntT(negOne);
			bigMod = free_BigIntT(bigMod);
			orderMustDivide = free_BigIntT(orderMustDivide);
			
			P = free_BigPolyT(P);
			tempPoly = free_BigPolyT(tempPoly);
			zeroPoly = free_BigPolyT(zeroPoly);
			unityPoly = free_BigPolyT(unityPoly);
			negOnePoly = free_BigPolyT(negOnePoly);
			lambdaPoly = free_BigPolyT(lambdaPoly);
			quotientPoly = free_BigPolyT(quotientPoly);
			remainderPoly = free_BigPolyT(remainderPoly);
			
			for (int i = 0; i < count_unique_factors(factoredP); i += 1)
			{
				coprimeFactors[i] = free_BigPolyT(coprimeFactors[i]);
				primePowersForOrder[i] = free_BigIntT(primePowersForOrder[i]);
			}
			FREE(coprimeFactors);
			FREE(primePowersForOrder);
			factoredP = free_BigFactorsT(factoredP);
		}
		
		
		//If the user wants to plug the update matrix into a polynomial expression
		else if (! strcmp(argv[1], "evalpoly"))
		{
			bool multiplyByInitial = FALSE;
			
			BigIntTP bigMod;
			BigPolyTP thePoly;
			
			BigIntMatrixTP evaluatedMatrix;
			
			//For if the user wants to multiply by INITIALMATRIX
			BigIntMatrixTP resultingMatrix = NULL;
			
			if (big_rows(UPDATEMATRIX) != big_cols(UPDATEMATRIX))
			{
				fprintf(stderr, "Given update matrix is not square.");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from .config file.");
			}
			
			if (oargv[3][0] != '\0')
			{
				if (!strcmp(oargv[3], "TRUE"))
				{
					//Check to make sure the given INITIALMATRIX is of the correct dimensions
					if (big_rows(INITIALMATRIX) == big_cols(UPDATEMATRIX))
						multiplyByInitial = TRUE;
					else
					{
						fprintf(stderr, "Given dimensions for initial matrix don't match with the given update matrix.\n");
						bigMod = free_BigIntT(bigMod);
						returnvalue = EXIT_SUCCESS;
						goto FREE_VARIABLES;
					}
				}
			}
			
			thePoly = read_BigPolyT(polynomialfilepath);
			
			if (thePoly == NULL)
			{
				fprintf(stderr, "Unable to read polynomial at %s.\n", polynomialfilepath);
				bigMod = free_BigIntT(bigMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			evaluatedMatrix = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_cols(UPDATEMATRIX));
			eval_BigPolyT(thePoly, UPDATEMATRIX, evaluatedMatrix, bigMod);
			
			printf("m (modulus): ");
			printi(bigMod);
			printf("\nf(x) (polynomial): ");
			printp(thePoly);
			printf("\nA (matrix):\n");
			printbm(UPDATEMATRIX);
			
			if (multiplyByInitial)
			{
				resultingMatrix = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_cols(INITIALMATRIX));
				big_mat_mul(evaluatedMatrix, INITIALMATRIX, resultingMatrix);
				modbm(resultingMatrix, bigMod);
				
				printf("\nB (initial matrix):\n");
				printbm(INITIALMATRIX);
				printf("\nf(A)B mod m:\n");
				printbm(resultingMatrix);
			}
			
			else
			{
				printf("\nf(A) mod m:\n");
				printbm(evaluatedMatrix);				
			}
			
			bigMod = free_BigIntT(bigMod);
			thePoly = free_BigPolyT(thePoly);
			evaluatedMatrix = free_BigIntMatrixT(evaluatedMatrix);
			resultingMatrix = free_BigIntMatrixT(resultingMatrix);
		}
		
		
		//Find the characteristic equation of the update matrix
		else if (! strcmp(argv[1], "chara"))
		{
			BigIntTP bigMod = NULL;
			
			BigIntMatrixTP bigMatrix;
			
			BigPolyTP bigEqn;
			BigFactorsTP bigEqnFactors = NULL;
			BigPolyTP* minEqn = NULL;
			
			//If the user provided a custom modulus
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			
			//Extract the modulus in the .config file otherwise
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			//Now, read the update matrix
			bigMatrix = read_BigIntMatrixT(updatefilepath);
			
			if (bigMatrix == NULL)
			{
				fprintf(stderr, "Unable to read matrix from %s.\n", updatefilepath);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Actually calculate the characteristic equation here
			printf("Matrix:\n");
			printbm(bigMatrix);
			printf("Modulus: ");
			printi(bigMod);
			printf("\n");
			bigEqn = chara_poly(bigMatrix, bigMod);
			printf("\nCharacteristic polynomial: ");
			printp(bigEqn);
			
			printf("\nFactored characteristic polynomial: ");
			bigEqnFactors = factor_BigPolyT(bigEqn, bigMod);
			printpf(bigEqnFactors);
			printf("\nMinimal polynomial: ");
			minEqn = min_poly(bigMatrix, NULL, bigMod, bigEqn);
			old_printpf(minEqn);
			printf("\n\n");
			
			bigMod  = free_BigIntT(bigMod);
			
			bigMatrix     = free_BigIntMatrixT(bigMatrix);
			bigEqn        = free_BigPolyT(bigEqn);
			bigEqnFactors = free_BigFactorsT(bigEqnFactors);
			minEqn        = free_BigPolyT_factors(minEqn);
		}
		
		
		//If the user wants to find all matrices with a particular characteristic polynomial
		else if (! strcmp(argv[1], "allcharas"))
		{
			//This should be rewritten to use .polynomial files to
			// define the characteristic polynomial, not command line
			// arguments. This would allow UNIX-type arguments to
			// be used.
			PROHIBIT_UNIXFLAGS
			
			BigIntTP bigMod;
			BigIntTP negOne;
			
			BigIntTP** currMatElements;
			BigIntMatrixTP currMat;
			
			BigIntTP* coeffs;
			BigPolyTP targetEqn;
			BigPolyTP tempPoly;
			BigPolyTP negTargetEqn;
			BigPolyTP tempEqn;
			BigPolyTP negOnePoly;
			
			bool newCycle;
			int tempCurrCycle;
			int numOfCycleLengths = 0;
			int maxArraySize = 0;
			int* cycleLengths = NULL;
			int* cycleLengthCounts = NULL;
			
			BigIntTP** currVectElements;
			BigIntMatrixTP currVect;
			CycleInfoTP theCycle = NULL;
			
			if (argc == 2)
			{
				printf("No characteristic polynomial was given on command line.\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			
			negOne = empty_BigIntT(1);
			subtract_BigIntT(bigMod, one, negOne);
			negOnePoly = new_BigPolyT(&negOne, 1);
			
			//Get coefficients for targetEqn
			coeffs = malloc((argc-2)*sizeof(BigIntTP));
			for (int i = 2; i < argc; i += 1)
				strtoBIT(argv[i], &coeffs[i-2]);
			
			targetEqn = new_BigPolyT(coeffs, argc-2);
			reduce_BigPolyT(targetEqn);
			
			negTargetEqn = new_BigPolyT(&one, 1);
			tempPoly = new_BigPolyT(&one, 1);
			
			multiply_BigPolyT(targetEqn, negOnePoly, tempPoly);
			mod_BigPolyT(tempPoly, bigMod, negTargetEqn);
			
			currMatElements = new_BigIntT_array(degree(targetEqn), degree(targetEqn));
			currMat = new_BigIntMatrixT(degree(targetEqn), degree(targetEqn));
			
			currVectElements = new_BigIntT_array(degree(targetEqn), 1);
			currVect = new_BigIntMatrixT(degree(targetEqn), 1);
			
			printf("Target characteristic eqn: ");
			printp(targetEqn);
			printf(" or ");
			printp(negTargetEqn);
			printf("\nModulus: ");
			printi(bigMod);
			printf("\n~~~\n");
			
			//Loop through all matrices, find the ones with the correct characteristic polynomial
			do
			{
				set_big_matrix(currMat, currMatElements);
				tempEqn = chara_poly(currMat, bigMod);
				
				//If we found a matrix with the correct characteristic polynomial
				if ((compare_BigPolyT(tempEqn, targetEqn) == 0) ||
				    (compare_BigPolyT(tempEqn, negTargetEqn) == 0))
				{
					printbm(currMat);
					
					//Now, let's find cycle length counts
					do
					{
						set_big_matrix(currVect, currVectElements);
						big_floyd(currMat, currVect, bigMod, &theCycle);
						newCycle = TRUE;
						
						//Try and see if we've found this cycle length before
						for (int i = 0; i < numOfCycleLengths; i += 1)
						{
							if (cycleLengths[i] == omega(theCycle))
							{
								cycleLengthCounts[i] += 1;
								newCycle = FALSE;
								break;
							}
						}
						
						//If we need to add a new cycle length to our list
						if (newCycle)
						{
							numOfCycleLengths += 1;
							
							if (maxArraySize < numOfCycleLengths)
							{
								cycleLengths = realloc(cycleLengths, numOfCycleLengths*sizeof(int));
								cycleLengthCounts = realloc(cycleLengthCounts, numOfCycleLengths*sizeof(int));
								
								maxArraySize = numOfCycleLengths;
							}
							
							cycleLengths[numOfCycleLengths-1] = omega(theCycle);
							cycleLengthCounts[numOfCycleLengths-1] = 1;
						}
					}
					while (!increment_BigIntT_array(currVectElements, degree(targetEqn), 1, one, bigMod));
					
					//Print out our found cycle lengths IN ORDER
					newCycle = TRUE;
					tempCurrCycle = 1;
					printf("1 (%d)", cycleLengthCounts[0]);
					while (newCycle)
					{
						newCycle = FALSE;
						for (int i = 0; i < numOfCycleLengths; i += 1)
							if (cycleLengths[i] > tempCurrCycle)
							{
								newCycle = TRUE;
								tempCurrCycle = cycleLengths[i];
								
								printf(", %d (%d)", tempCurrCycle, cycleLengthCounts[i]);
							}
					}
					printf("\n\n");

					numOfCycleLengths = 0;
				}
				
				tempEqn = free_BigPolyT(tempEqn);
			}
			while (!increment_BigIntT_array(currMatElements, degree(targetEqn), degree(targetEqn), one, bigMod));
			
			
			bigMod = free_BigIntT(bigMod);
			negOne = free_BigIntT(negOne);
			
			currMatElements  = free_BigIntT_array(currMatElements, big_rows(currMat), big_rows(currMat));
			currVectElements = free_BigIntT_array(currVectElements, big_rows(currMat), 1);
			currMat  = free_BigIntMatrixT(currMat);
			currVect = free_BigIntMatrixT(currVect);
			theCycle = free_CycleInfoT(theCycle);
			FREE(cycleLengths);
			FREE(cycleLengthCounts);
			
			for (int i = 0; i < argc-2; i += 1)
				coeffs[i] = free_BigIntT(coeffs[i]);
			FREE(coeffs);
			
			targetEqn    = free_BigPolyT(targetEqn);
			negTargetEqn = free_BigPolyT(negTargetEqn);
			negOnePoly   = free_BigPolyT(negOnePoly);
			tempPoly     = free_BigPolyT(tempPoly);
		}
		
		
		//If the user wants to calculate the specific vectors in a
		// matrix's orbits
		else if (! strcmp(argv[1], "orbits"))
		{
			BigIntTP bigMod;
			BigIntTP overflowCheck;
			int smallMod;
			int totalVects = 1;
			int foundVectors = 0;
			int tempVectCount = 0;
			BigIntTP tempIndex1;
			BigIntTP tempIndex2;
			
			BigIntMatrixTP  A;
			BigIntMatrixTP  currVect;
			BigIntTP** currVectElements;
			
			CycleInfoTP theCycle = NULL;
			
			int numOfOrbits = 0;
			BigIntMatrixTP* orbitReps = NULL;
			bool foundNewOrbit = FALSE;
			
			int oneArr[1] = {1};
			BigIntTP temp;
			BigIntTP temp2;
			BigIntTP temp3;
			
			BigIntMatrixTP tempVect;
			BigIntMatrixTP tempVect2;
			
			FILE* outputFile = NULL;
			char* outputFileName = NULL;
			bool fileoutput = FALSE;
			
			FILE* graphFile = NULL;
			char* graphFileName = NULL;
			
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Invalid modulus passed on command line.");
				smallMod = (int)strtol(oargv[2], &tempStr, 10);
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Invalid modulus passed in config file.");
				smallMod = (int)strtol(bigintmodstring, &tempStr, 10);
			}
			
			if (oargv[3][0] != '\0')
				if (! strcmp(oargv[3], "TRUE"))
					fileoutput = TRUE;
				
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from config file.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (big_rows(A) != big_cols(A))
			{
				fprintf(stderr, "Given update matrix isn't square.\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
				
			//Creating filename and file
			if (fileoutput)
			{
				outputFileName = malloc(MAXSTRLEN*sizeof(char));
				outputFileName[0] = '\0';
				
				graphFileName = malloc(MAXSTRLEN*sizeof(char));
				graphFileName[0] = '\0';
				
				ADD_FILE_PREFIX(outputFileName);
				ADD_FILE_PREFIX(graphFileName);
				
				strcat(outputFileName, "orbits ");
				strcat(graphFileName, "graph ");
				append_BigIntT(outputFileName, bigMod);
				append_BigIntT(graphFileName, bigMod);
				strcat(outputFileName, " F");
				strcat(graphFileName, " F");
				for (int row = 0; row < big_rows(A); row += 1)
					for (int col = 0; col < big_rows(A); col += 1)
					{
						append_BigIntT(outputFileName, big_element(A, row, col));
						append_BigIntT(graphFileName, big_element(A, row, col));
					}
				strcat(outputFileName, ".txt");
				strcat(graphFileName, ".graph");
				
				outputFile = fopen(outputFileName, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
				
				graphFile = fopen(graphFileName, "w");
				if (graphFile == NULL)
					fprintf(stderr, "Unable to create graph file.\n");
			}
			
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			currVect = new_BigIntMatrixT(big_rows(A), 1);
			
			temp = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			temp3 = empty_BigIntT(1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			
			tempIndex1 = empty_BigIntT(1);
			tempIndex2 = empty_BigIntT(1);
			
			overflowCheck = new_BigIntT(oneArr, 1);
			for (int a = 0; a < big_rows(A); a += 1)
			{
				multiply_BigIntT(bigMod, overflowCheck, temp);
				copy_BigIntT(temp, overflowCheck);
				
				totalVects *= smallMod;
			}
			
			//Iterate through all vectors, find all orbits
			do
			{
				foundNewOrbit = FALSE;
				
				set_big_matrix(currVect, currVectElements);
				big_floyd(A, currVect, bigMod, &theCycle);
				
				if (numOfOrbits == 0)
				{
					foundNewOrbit = TRUE;
					
					numOfOrbits += 1;
					orbitReps = malloc(sizeof(BigIntMatrixTP));
					orbitReps[0] = new_BigIntMatrixT(big_rows(A), 1);
					
					copy_BigIntMatrixT(rep(theCycle), orbitReps[0]);
				}
				
				//Sift through our currently-found orbits to see if our currVect is in them
				else
				{
					foundNewOrbit = TRUE;
					for (int ourRep = 0; ourRep < numOfOrbits; ourRep += 1)
					{
						copy_BigIntMatrixT(orbitReps[ourRep], tempVect);
						do
						{
							//If we found our rep of currVect's cycle in a cycle
							if (compare_BigIntMatrixT(tempVect, rep(theCycle)))
							{
								foundNewOrbit = FALSE;
								break;
							}
							
							big_mat_mul(A, tempVect, tempVect2);
							modbm(tempVect2, bigMod);
							copy_BigIntMatrixT(tempVect2, tempVect);
						}
						while (!compare_BigIntMatrixT(tempVect, orbitReps[ourRep]));
						
						if (!foundNewOrbit)
							break;
					}
					
					if (foundNewOrbit)
					{
						numOfOrbits += 1;
						orbitReps = realloc(orbitReps, numOfOrbits*sizeof(BigIntMatrixTP));
						orbitReps[numOfOrbits-1] = new_BigIntMatrixT(big_rows(A), 1);
						copy_BigIntMatrixT(rep(theCycle), orbitReps[numOfOrbits-1]);
					}
				}
				
				if (foundNewOrbit)
				{
					printf("Orbit #%d:\n", numOfOrbits-1);
					copy_BigIntMatrixT(orbitReps[numOfOrbits-1], tempVect);
					
					if (outputFile != NULL)
						fprintf(outputFile, "Orbit #%d (ω = %d)\n", numOfOrbits-1, omega(theCycle));
					
					for (int v = 0; v < omega(theCycle); v += 1)
					{
						foundVectors += 1;
						
						printbm_row(tempVect);
						printf("\n");
						
						if (outputFile != NULL)
						{
							fprintf(outputFile, "%d: ", v);
							fprintbm_row(outputFile, tempVect);
							fprintf(outputFile, "\n");
						}
						
						big_mat_mul(A, tempVect, tempVect2);
						modbm(tempVect2, bigMod);
						copy_BigIntMatrixT(tempVect2, tempVect);
					}
					
					printf("Vectors in cycles: %d\n", foundVectors);
					printf("\n");
					if (outputFile != NULL)
						fprintf(outputFile, "\n");
				}
				
				//If we've accounted for all the vectors in the module
				//The first condition prevents this from triggering if 
				if ((size(overflowCheck) == 1) && (foundVectors == totalVects))
					break;
			}
			while (!increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
			
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
				
			if (graphFile != NULL)
			{
				fprintf(graphFile, "~b:False\n");
				fprintf(graphFile, "~n:%d\n~l:above\n", foundVectors);
				for (int i = 0; i < big_rows(A); i += 1)
					copy_BigIntT(zero, currVectElements[i][0]);
				
				do
				{
					set_big_matrix(currVect, currVectElements);
					fprintf(graphFile, "%d:", tempVectCount);
					fprintbm_row(graphFile, currVect);
					fprintf(graphFile, "\n");
					tempVectCount += 1;
				}
				while (!increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
				
				fprintf(graphFile, "~c:direction\n");
				
				//Go through every cycle, add the cycles to the file
				for (int r = 0; r < numOfOrbits; r += 1)
				{
					copy_BigIntMatrixT(orbitReps[r], tempVect);
					big_mat_mul(A, tempVect, tempVect2);
					modbm(tempVect2, bigMod);
					
					//Loop around until we've recorded every mapping in the orbit
					do
					{
						//Get numbers for each vector
						copy_BigIntT(zero, tempIndex1);
						copy_BigIntT(zero, tempIndex2);
						copy_BigIntT(one, temp);
						
						for (int i = 0; i < big_rows(A); i += 1)
						{
							//Get correct powers of our modulus
							for (int j = 0; j < i; j += 1)
							{
								multiply_BigIntT(temp, bigMod, temp2);
								copy_BigIntT(temp2, temp);
							}
							
							//temp holds correct power
							multiply_BigIntT(temp, big_element(tempVect, i, 0), temp2);
							add_BigIntT(temp2, tempIndex1, temp3);
							copy_BigIntT(temp3, tempIndex1);
							
							multiply_BigIntT(temp, big_element(tempVect2, i, 0), temp2);
							add_BigIntT(temp2, tempIndex2, temp3);
							copy_BigIntT(temp3, tempIndex2);
						}
						
						//Now our tempIndices hold the correct numbers representing each vector
						fprinti(graphFile, tempIndex1);
						fprintf(graphFile, ",");
						fprinti(graphFile, tempIndex2);
						fprintf(graphFile, "\n");
						
						//Update our vectors to the next mapping in the orbit
						copy_BigIntMatrixT(tempVect2, tempVect);
						big_mat_mul(A, tempVect, tempVect2);
						modbm(tempVect2, bigMod);
					}
					while (!compare_BigIntMatrixT(tempVect, orbitReps[r]));
				}
				
				if (fclose(graphFile) == EOF)
					fprintf(stderr, "Unable to save graph file.\n");
			}
			
			FREE(outputFileName);
			FREE(graphFileName);
			bigMod = free_BigIntT(bigMod);
			
			for (int i = 0; i < numOfOrbits; i += 1)
				orbitReps[i] = free_BigIntMatrixT(orbitReps[i]);
			FREE(orbitReps);
			
			currVectElements = free_BigIntT_array(currVectElements, big_rows(A), 1);
			currVect = free_BigIntMatrixT(currVect);
			
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect);
			
			tempIndex1 = free_BigIntT(tempIndex1);
			tempIndex2 = free_BigIntT(tempIndex2);
			temp = free_BigIntT(temp);
			temp2 = free_BigIntT(temp);
			overflowCheck = free_BigIntT(overflowCheck);
			
			theCycle = free_CycleInfoT(theCycle);
			
			A = free_BigIntMatrixT(A);
		}
		
		
		//If the user wants to calculate orbits over a splitting field for
		// the given update matrix's minimal polynomial
		else if (! strcmp(argv[1], "splitorbits"))
		{
			//Let's just hope we don't have more extensions than this
			char* extensionNames[5] = {"α", "β", "δ", "ε", "θ"};
			int uniqueExtensionCombos = 1; //For helping determine how many vectors are possible
			
			BigIntTP bigMod;
			BigIntMatrixTP bigF;
			
			BigPolyTP* minPolyFactors; 
			BigPolyTP* extDefns = NULL; //For extracting coeffs from min poly factors
			BigIntTP*  tempFactor;
			bool isSameDefn;
			
			//How many extensions we need to add to our field
			int  numOfExtensionsNeeded = 0;
			int* zeroPointer = NULL; //For placing coefficients into MultiVarExtTs
			
			//Will hold all our extensions
			MultiVarExtTP repExtension = NULL;
			
			void*** genericFelements = NULL;
			GenericMatrixTP genericF = NULL;
			
			void*** currVectElements = NULL;
			GenericMatrixTP currVect = NULL;
			
			GenericMatrixTP tempVect = NULL; //Facilitates calculations
			
			CycleInfoTP theCycle = NULL;
			bool newCycleLength  = FALSE;
			bool newCycle        = FALSE;
			int foundCycleIndex; //Used for holding a position in the arrays below
			
			int     numOfCycleLengthsFound = 0;
			int*    numOfCyclesFound = NULL;  //Holds how many cycles of each length we've found
			int*    foundCycleLengths = NULL; //Holds all the cycle lengths we've come across so far
			void*** foundReps;                //Holds representatives for all the cycles we've found thus far
			
			bool moreVectors = TRUE;
			
			BigIntTP temp, temp2;
			
			int cycleCounter = 0; //How many cycles have we found?
			BigIntTP vectorCounter = NULL;
			BigIntTP totalVectors = NULL;
			
			char* fileOutputName = NULL;
			FILE* fileOutput = NULL;
			
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from .config file.");
			}
			
			bigF = read_BigIntMatrixT(updatefilepath);
			if (bigF == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from %s.\n", updatefilepath);
				bigMod = free_BigIntT(bigMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (big_rows(bigF) != big_cols(bigF))
			{
				fprintf(stderr, "Given matrix isn't square.\n");
				bigMod = free_BigIntT(bigMod);
				bigF = free_BigIntMatrixT(bigF);
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			if (oargv[3][0] != '\0')
				if (!strcmp(oargv[3], "TRUE"))
				{
					fileOutputName = malloc(MAXSTRLEN*sizeof(char));
					fileOutputName[0] = '\0';

					ADD_FILE_PREFIX(fileOutputName);
					
					strcat(fileOutputName, "splitorbits F");
					
					//Get elements of matrix
					for (int row = 0; row < big_rows(bigF); row += 1)
						for (int col = 0; col < big_cols(bigF); col += 1)
							append_BigIntT(fileOutputName, big_element(bigF, row, col));
						
					strcat(fileOutputName, " ");
					append_BigIntT(fileOutputName, bigMod);
					strcat(fileOutputName, ".txt");
					
					fileOutput = fopen(fileOutputName, "a");
					if (fileOutput == NULL)
						fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
				}
				
			vectorCounter = empty_BigIntT(1);
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			
			printf("Update matrix:\n");
			printbm(bigF);
			printf("Modulus: ");
			printi(bigMod);
			minPolyFactors = min_poly(bigF, NULL, bigMod, NULL);
			printf("\nMin poly: ");
			old_printpf(minPolyFactors);
			printf("\n");
			
			//Now, we check to see how many extensions we need for
			// our splitting field
			copy_BigIntT(one, temp);
			for (int i = 1; compare_BigIntT(constant(minPolyFactors[0]), temp) >= 0; i += 1)
			{
				isSameDefn = FALSE;
				if (degree(minPolyFactors[i]) >= 2)
				{
					//We need to ensure that this factor is unique.
					for (int extDefn = 0; extDefn < numOfExtensionsNeeded; extDefn += 1)
					{
						isSameDefn = TRUE;
						if (degree(extDefns[extDefn]) == degree(minPolyFactors[i]))
						{
							if (compare_BigPolyT(extDefns[extDefn], minPolyFactors[i]) != 0)
								isSameDefn = FALSE;
						}
						
						else
							isSameDefn = FALSE;
						
						if (isSameDefn)
							break;
					}
					
					//Only add new extension if we don't already have it
					if (!isSameDefn)
					{
						numOfExtensionsNeeded += 1;
						uniqueExtensionCombos *= degree(minPolyFactors[i]);
						
						extDefns = realloc(extDefns, numOfExtensionsNeeded*sizeof(BigPolyTP));
						extDefns[numOfExtensionsNeeded-1] = empty_BigPolyT();
						copy_BigPolyT(minPolyFactors[i], extDefns[numOfExtensionsNeeded-1]);
					}
				}
				
				add_BigIntT(temp, one, temp2);
				copy_BigIntT(temp2, temp);
			}
			
			if (numOfExtensionsNeeded == 0)
			{
				printf("\nThe splitting field for this system contains no new extensions.\n" \
			         "Use another tool to analyse this system.\n");
			}
			
			else
			{
				//First, count how many vectors are in our space
				totalVectors = empty_BigIntT(1);
				copy_BigIntT(one, totalVectors);
				for (int i = 0; i < uniqueExtensionCombos*big_rows(bigF); i += 1)
				{
					multiply_BigIntT(bigMod, totalVectors, temp);
					copy_BigIntT(temp, totalVectors);
				}
				
				//We have the number of extensions we need, now create a 
				// representative MultiVarExtT that has all the extensions
				// we need.
				
				repExtension = new_MultiVarExtT(numOfExtensionsNeeded);
				set_MultiVarExtT_mod(repExtension, bigMod);
				
				zeroPointer = calloc(numOfExtensionsNeeded, sizeof(int));
				
				//Adding the extensions we need
				for (int i = 0; i < numOfExtensionsNeeded; i += 1)
				{
					tempFactor = extract_coefficients(extDefns[i]);
					add_extension(repExtension, tempFactor, degree(extDefns[i])+1, extensionNames[i]);
					
					for (int tempCoeff = 0; tempCoeff < degree(extDefns[i])+1; tempCoeff += 1)
						tempFactor[tempCoeff] = free_BigIntT(tempFactor[tempCoeff]);
					FREE(tempFactor);
				}
				
				//Now, initialise our matrix and vector
				genericFelements = malloc(big_rows(bigF)*sizeof(void**));
				currVectElements = malloc(big_rows(bigF)*sizeof(void**));
				for (int i = 0; i < big_rows(bigF); i += 1)
				{
					genericFelements[i] = malloc(big_cols(bigF)*sizeof(void*));
					currVectElements[i] = malloc(sizeof(void*));
					for (int j = 0; j < big_cols(bigF); j += 1)
					{
						genericFelements[i][j] = new_MultiVarExtT(numOfExtensionsNeeded);
						copy_MultiVarExtT(repExtension, genericFelements[i][j]);
						set_MultiVarExtT_coefficient(genericFelements[i][j], zeroPointer, big_element(bigF, i, j));
					}
					currVectElements[i][0] = new_MultiVarExtT(numOfExtensionsNeeded);
					copy_MultiVarExtT(repExtension, currVectElements[i][0]);
				}
				
				genericF = new_MultiVarExtMatrixT(big_rows(bigF), big_cols(bigF), numOfExtensionsNeeded);
				currVect = new_MultiVarExtMatrixT(big_rows(bigF), 1, numOfExtensionsNeeded);
				tempVect = new_MultiVarExtMatrixT(big_rows(bigF), 1, numOfExtensionsNeeded);
				set_GenericMatrixT(genericF, genericFelements);
				set_GenericMatrixT(currVect, currVectElements);
				
				//Setting tempVect's extensions
				set_GenericMatrixT(tempVect, currVectElements);
				
				foundCycleLengths = calloc(1, sizeof(int));
				numOfCyclesFound  = calloc(1, sizeof(int));
				foundReps = NULL;
				
				//Write MultiVarExt definitions to output file
				if (fileOutput != NULL)
				{
					fprintmve(fileOutput, repExtension);
					fprintf(fileOutput, "\n");
					if (fclose(fileOutput) == EOF)
						fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
					else
					{
						fileOutput = fopen(fileOutputName, "a");
						if (fileOutput == NULL)
							fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
					}
				}
				
				//Iterate through all possible vectors over our splitting field
				while (moreVectors)
				{
					//Now, let's actually do the calculations
					generic_floyd(genericF, currVect, &theCycle);
					
					//Check to see if our found cycle length has been found before
					newCycleLength = TRUE;
					for (int cyclen = 0; (cyclen < numOfCycleLengthsFound) && (newCycleLength); cyclen += 1)
						if (foundCycleLengths[cyclen] == omega(theCycle))
						{
							newCycleLength = FALSE;
							foundCycleIndex = cyclen;
						}
						
					//Record data for new cycle found
					if (newCycleLength)
					{
						foundCycleIndex = numOfCycleLengthsFound;
						
						numOfCycleLengthsFound += 1;
						foundCycleLengths = realloc(foundCycleLengths, numOfCycleLengthsFound*sizeof(int));
						numOfCyclesFound  = realloc(numOfCyclesFound, numOfCycleLengthsFound*sizeof(int));
						
						numOfCyclesFound[foundCycleIndex] = 0;
						foundCycleLengths[foundCycleIndex] = omega(theCycle);
						
						foundReps = realloc(foundReps, numOfCycleLengthsFound*sizeof(void**));
						foundReps[numOfCycleLengthsFound-1] = NULL;
					}
					
					//Search to make sure that the cycle we calculated doesn't already have a rep
					newCycle = TRUE;
					for (int repr = 0; (repr < numOfCyclesFound[foundCycleIndex]) && newCycle; repr += 1)
					{
						//Iterate the current representative, see if our new cycle's rep is somewhere
						// in the cycle.
						copy_sim_GenericMatrixT(foundReps[foundCycleIndex][repr], currVect);
						for (int iter = 0; (iter < foundCycleLengths[foundCycleIndex]-1) && newCycle; iter += 1)
						{
							gen_mat_mul(genericF, currVect, tempVect);
							reduce_GenericMatrixT(tempVect);

							if (compare_GenericMatrixT(tempVect, rep(theCycle)))
								newCycle = FALSE;
							
							copy_sim_GenericMatrixT(tempVect, currVect);
						}
					}
					
					//Now, add the new cycle representative if needed
					if (newCycle)
					{
						numOfCyclesFound[foundCycleIndex] += 1;
						foundReps[foundCycleIndex] = realloc(foundReps[foundCycleIndex], numOfCyclesFound[foundCycleIndex]*sizeof(void*));
						foundReps[foundCycleIndex][numOfCyclesFound[foundCycleIndex]-1] = new_MultiVarExtMatrixT(big_rows(bigF), 
																																																		 1, 
																																																		 numOfExtensionsNeeded);
						//Making sure our new foundReps vector has the correct extension definitions and stuff
						set_GenericMatrixT(foundReps[foundCycleIndex][numOfCyclesFound[foundCycleIndex]-1], currVectElements);
						copy_sim_GenericMatrixT(rep(theCycle), foundReps[foundCycleIndex][numOfCyclesFound[foundCycleIndex]-1]);
						
						//Now, we print out the new cycle
						printf("Cycle #%d:\n", cycleCounter++);
						if (fileOutput != NULL)
						{
							fprintf(fileOutput, "Cycle #%d:\n", cycleCounter-1);
							if (fclose(fileOutput) == EOF)
								fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							else
							{
								fileOutput = fopen(fileOutputName, "a");
								if (fileOutput == NULL)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							}
						}
						
						set_GenericMatrixT(currVect, currVectElements);
						
						printgm_row(currVect);
						printf("\n");
						if (fileOutput != NULL)
						{
							fprintgm_row(fileOutput, currVect);
							fprintf(fileOutput, "\n");
							if (fclose(fileOutput) == EOF)
								fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							else
							{
								fileOutput = fopen(fileOutputName, "a");
								if (fileOutput == NULL)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							}
						}
						
						add_BigIntT(vectorCounter, one, temp);
						copy_BigIntT(temp, vectorCounter);
							
						copy_sim_GenericMatrixT(rep(theCycle), currVect);
						for (int iter = 0; iter < foundCycleLengths[foundCycleIndex]-1; iter += 1)
						{
							add_BigIntT(vectorCounter, one, temp);
							copy_BigIntT(temp, vectorCounter);
							
							gen_mat_mul(genericF, currVect, tempVect);
							reduce_GenericMatrixT(tempVect);
							
							printgm_row(tempVect);
							printf("\n");
							if (fileOutput != NULL)
							{
								fprintgm_row(fileOutput, tempVect);
								fprintf(fileOutput, "\n");
								if (fclose(fileOutput) == EOF)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
								else
								{
									fileOutput = fopen(fileOutputName, "a");
									if (fileOutput == NULL)
										fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
								}
							}
							
							copy_sim_GenericMatrixT(tempVect, currVect);
						}
						printf("Vectors in cycles: ");
						printi(vectorCounter);
						printf("\n\n");
						if (fileOutput != NULL)
						{
							fprintf(fileOutput, "\n");
							if (fclose(fileOutput) == EOF)
								fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							else
							{
								fileOutput = fopen(fileOutputName, "a");
								if (fileOutput == NULL)
									fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
							}
						}
					}
					
					//Iterate our vector
					moreVectors = FALSE;
					if (compare_BigIntT(vectorCounter, totalVectors) != 0)
					{
						for (int component = 0; component < big_rows(bigF); component += 1)
						{
							if (! increment_MultiVarExtT(currVectElements[component][0]))
							{
								moreVectors = TRUE;
								set_GenericMatrixT(currVect, currVectElements);
								break;
							}
						}
					}
				}
			}
			
			repExtension = free_MultiVarExtT(repExtension);
			
			copy_BigIntT(one, temp);
			for (int i = 1; compare_BigIntT(constant(minPolyFactors[0]), temp) >= 0; i += 1)
			{
				minPolyFactors[i] = free_BigPolyT(minPolyFactors[i]);
				add_BigIntT(temp, one, temp2);
				copy_BigIntT(temp2, temp);
			}
			minPolyFactors[0] = free_BigPolyT(minPolyFactors[0]);
			FREE(minPolyFactors);
				
			if (fileOutput != NULL)
				if (fclose(fileOutput) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
			
			if (fileOutputName != NULL)
			{
				FREE(fileOutputName);
			}
			
			for (int i = 0; i < numOfExtensionsNeeded; i += 1)
				extDefns[i] = free_BigPolyT(extDefns[i]);

			FREE(extDefns);
			
			if (genericFelements != NULL)
			{
				for (int i = 0; i < big_rows(bigF); i += 1)
				{
					for (int j = 0; j < big_cols(bigF); j += 1)
						genericFelements[i][j] = free_MultiVarExtT(genericFelements[i][j]);
					
					FREE(genericFelements[i]);
				}
				FREE(genericFelements);
			}
			
			if (currVectElements != NULL)
			{
				for (int i = 0; i < big_rows(bigF); i += 1)
				{
					currVectElements[i][0] = free_MultiVarExtT(currVectElements[i][0]);
					FREE(currVectElements[i]);
				}
				FREE(currVectElements);
			}
		
			bigF = free_BigIntMatrixT(bigF);
			genericF = free_GenericMatrixT(genericF);
			currVect = free_GenericMatrixT(currVect);
			tempVect = free_GenericMatrixT(tempVect);
			
			if (zeroPointer != NULL)
			{
				FREE(zeroPointer);
			}
			
			theCycle = free_CycleInfoT(theCycle);
			FREE(foundCycleLengths);
			for (int cyclen = 0; cyclen < numOfCycleLengthsFound; cyclen += 1)
			{
				for (int cyc = 0; cyc < numOfCyclesFound[cyclen]; cyc += 1)
				{
					foundReps[cyclen][cyc] = free_GenericMatrixT(foundReps[cyclen][cyc]);
				}
				FREE(foundReps[cyclen]);
			}
			FREE(foundReps);
			FREE(numOfCyclesFound);
			
			temp   = free_BigIntT(temp);
			temp2  = free_BigIntT(temp2);
			bigMod = free_BigIntT(bigMod);
			
			totalVectors = free_BigIntT(totalVectors);
			vectorCounter = free_BigIntT(vectorCounter);
		}
		
		
		//If we want to find a representative from each of the update matrix's orbits
		else if (! strcmp(argv[1], "orbitreps"))
		{
			BigIntTP bigMod;
			
			BigIntMatrixTP A;
			BigIntMatrixTP currVect;
			BigIntMatrixTP tempVect;
			BigIntMatrixTP tempVect2;
			
			BigIntMatrixTP** reps = NULL;
			
			int** foundCycleLengths; //[[cycle length, how many reps have cycle length]]
			BigIntTP** currVectElements;
			
			CycleInfoTP theCycle = NULL;
			
			int index;
			bool isNewCycle;
			bool hasRep;
			
			//Says whether to create the .txt output
			bool fileoutput = FALSE;
			
			char* outputfilename;
			FILE* outputFile;
			
			//If user provided modulus at command line
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			//If the user wanted to change the default file output behaviour
			if (oargv[3][0] != '\0')
				if (!strcmp(oargv[3], "TRUE"))
					fileoutput = TRUE;
			
			A = UPDATEMATRIX;
			if ((A == NULL))
			{
				fprintf(stderr, "Unable to read update matrix file.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (big_rows(A) != big_cols(A))
			{
				fprintf(stderr, "Given update matrix is not square.\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			//First entry says how many cycle lengths we've found so far
			foundCycleLengths = malloc(sizeof(int*));
			foundCycleLengths[0] = calloc(2, sizeof(int));
			
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			
			//Construct output file name
			outputfilename = malloc(MAXSTRLEN*sizeof(char));
			outputfilename[0] = '\0';
			
			//If we have a default path to put output files into
			ADD_FILE_PREFIX(outputfilename);
			
			strcat(outputfilename, "orbitreps ");
			append_BigIntT(outputfilename, bigMod);
			strcat(outputfilename, " F");
			
			//Matrix
			for (int row = 0; row < big_rows(A); row += 1)
				for (int col = 0; col < big_rows(A); col += 1)
					append_BigIntT(outputfilename, big_element(A, row, col));
				
			strcat(outputfilename, ".txt");
			
			currVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			
			//Increment through all possible vectors
			printf("Matrix:\n");
			printbm(A);
			printf("Mod ");
			printi(bigMod);
			printf("\n\n");
			do
			{
				//Getting cycle length, adding to list if needed
				set_big_matrix(currVect, currVectElements);
				big_floyd(A, currVect, bigMod, &theCycle);
				
				isNewCycle = TRUE;
				for (int i = 1; i <= foundCycleLengths[0][0]; i += 1)
					if (foundCycleLengths[i][0] == omega(theCycle))
					{
						isNewCycle = FALSE;
						index = i;
						break;
					}
					
				if (isNewCycle)
				{
					//Allocate space for new cycle length info
					foundCycleLengths[0][0] += 1;
					foundCycleLengths = realloc(foundCycleLengths, (foundCycleLengths[0][0]+1)*sizeof(int*));
					foundCycleLengths[foundCycleLengths[0][0]] = calloc(2, sizeof(int));
					foundCycleLengths[foundCycleLengths[0][0]][0] = omega(theCycle);
					index = foundCycleLengths[0][0];
					
					//Also allocate space for new representatives
					reps = realloc(reps, foundCycleLengths[0][0]*sizeof(BigIntMatrixTP*));
					reps[index-1] = NULL;
				}
				
				//Now that we have our cycle length, we have to check every vector
				// in our currVect's orbit to see if we already have a representative
				if (foundCycleLengths[index][1] != 0)
				{
					copy_BigIntMatrixT(rep(theCycle), tempVect);
					hasRep = FALSE;
					
					//Now, check tempVect's orbit for reps
					for (int i = 0; i < omega(theCycle); i += 1)
					{
						for (int r = 0; r < foundCycleLengths[index][1]; r += 1)
							if (compare_BigIntMatrixT(tempVect, reps[index-1][r]))
							{
								hasRep = TRUE;
								break;
							}
							
						if (hasRep)
							break;
						
						big_mat_mul(A, tempVect, tempVect2);
						modbm(tempVect2, bigMod);
						copy_BigIntMatrixT(tempVect2, tempVect);
					}
					
					//If we found a new representative, add it to the list
					if (!hasRep)
					{
						foundCycleLengths[index][1] += 1;
						reps[index-1] = realloc(reps[index-1], foundCycleLengths[index][1]*sizeof(BigIntMatrixTP));
						
						reps[index-1][foundCycleLengths[index][1]-1] = new_BigIntMatrixT(big_rows(A), 1);
						copy_BigIntMatrixT(rep(theCycle), reps[index-1][foundCycleLengths[index][1]-1]);
					}
				}
			
				//If we don't have any representatives yet
				else
				{
					//Increment number of reps we've found
					foundCycleLengths[index][1] += 1;
					reps[index-1] = realloc(reps[index-1], foundCycleLengths[index][1]*sizeof(BigIntMatrixTP));
					
					reps[index-1][0] = new_BigIntMatrixT(big_rows(A), 1);
					copy_BigIntMatrixT(rep(theCycle), reps[index-1][0]);
				}
			}
			while (! increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
			
			//Create save file
			if (fileoutput)
			{
				outputFile = fopen(outputfilename, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create save file. Results will still print to stdout.\n");
			}
			else
				outputFile = NULL;
			
			//Now, print out the representatives!
			for (int cyclenindex = 1; cyclenindex <= foundCycleLengths[0][0]; cyclenindex += 1)
			{
				printf("Cycle length: %d\n", foundCycleLengths[cyclenindex][0]);
				printf("Representatives:\n");
				
				if (outputFile != NULL)
				{
					fprintf(outputFile, "Cycle length: %d\n", foundCycleLengths[cyclenindex][0]);
					fprintf(outputFile, "Reps: \n");
				}
				
				for (int cycrep = 0; cycrep < foundCycleLengths[cyclenindex][1]; cycrep += 1)
				{
					printbm_row(reps[cyclenindex-1][cycrep]);
					printf("\n");
					
					if (outputFile != NULL)
					{						
						fprintbm_row(outputFile, reps[cyclenindex-1][cycrep]);
						fprintf(outputFile, "\n");
					}
				}
				
				if (outputFile != NULL)
					fprintf(outputFile, "\n");
			}
			
			if (fileoutput)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save data.\n");
			
			
			currVectElements = free_BigIntT_array(currVectElements, big_rows(A), 1);
			
			for (int i = 1; i <= foundCycleLengths[0][0]; i += 1)
			{
				for (int j = 0; j < foundCycleLengths[i][1]; j += 1)
					reps[i-1][j] = free_BigIntMatrixT(reps[i-1][j]);
				
				FREE(reps[i-1]);
				FREE(foundCycleLengths[i]);
			}
			FREE(reps);

			FREE(foundCycleLengths[0]);
			FREE(foundCycleLengths);
			
			bigMod = free_BigIntT(bigMod);
			
			theCycle = free_CycleInfoT(theCycle);

			currVect  = free_BigIntMatrixT(currVect);
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect2);
			
			FREE(outputfilename);
		}
		
		
		//If we want to generate transient length representatives
		else if (!strcmp(argv[1], "branchreps"))
		{
			BigIntMatrixTP A;
			BigIntMatrixTP I;
			BigIntMatrixTP currVect  = NULL;
			BigIntMatrixTP tempVect  = NULL;
			BigIntMatrixTP tempVect2 = NULL;
			BigIntMatrixTP tempVect3 = NULL;
			BigIntMatrixTP tempRoot  = NULL;
			BigIntTP**     currVectElements;
			
			CycleInfoTP theCycle = NULL;
			
			/* First vector is the "leaf" of the branch, the vector with more steps to get into a cycle
			 * Second vector is the "root" of the branch, the vector closest to a cycle (and which 
			 *  possibly leads to another branch).
			 * Third vector is the "base root" of the branch, the root you'd get to if you follow the branch down
			 */
			BigIntMatrixTP** branches; //Holds the transient branch representations
			int numOfBranches = 0;
			
			//Holds references to the vectors we're using
			//Helps to save on memory usage since there's a lot of redundancy in our tree representation
			BigIntMatrixTP* vectCatalogue = NULL;
			int sizeOfCatalogue = 0;
			
			BigIntTP bigMod;
			
			bool isNewRoot;
			bool foundGrowth; //If we've found the place where a vector grows from
			bool vectInsideVine;
			bool hasConnections = TRUE;
			bool isTopLeaf;   //For determining which leaves are strictly necessary for transient length representatives
			
			char* fileName = NULL;
			FILE* outputFile = NULL;
			
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus fom command line.");
			}
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read matrix in %s.\n", updatefilepath);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Checking to ensure a transient region exists
			I = identity_BigIntMatrixT(big_rows(A));
			big_floyd(A, I, bigMod, &theCycle);
			if (tau(theCycle) == 0)
			{
				printf("Given matrix is invertible under given modulus, so no transient region exists.\n");
				
				bigMod = free_BigIntT(bigMod);
				A = free_BigIntMatrixT(A);
				I = free_BigIntMatrixT(I);
				theCycle = free_CycleInfoT(theCycle);
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			fileName = malloc(MAXSTRLEN*sizeof(char));
			fileName[0] = '\0';
			ADD_FILE_PREFIX(fileName);
			
			strcat(fileName, "branchreps ");
			append_BigIntT(fileName, bigMod);
			strcat(fileName, " F");
			for (int row = 0; row < big_rows(A); row += 1)
				for (int col = 0; col < big_cols(A); col += 1)
					append_BigIntT(fileName, big_element(A, row, col));
			strcat(fileName, ".txt");
			
			outputFile = fopen(fileName, "a");
			if (outputFile == NULL)
				fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
			
			branches    = malloc(sizeof(BigIntMatrixTP*));
			branches[0] = malloc(3*sizeof(BigIntMatrixTP));
			
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			
			currVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			tempVect3 = new_BigIntMatrixT(big_rows(A), 1);
			tempRoot  = new_BigIntMatrixT(big_rows(A), 1);
			
			printf("Modulus: ");
			printi(bigMod);
			printf("\nMatrix:\n");
			printbm(A);
			printf("~~~~~~~~~~\n");
			
			//Loop over each vector, check if it's in a cycle or not
			do
			{
				//This is here to prevent the program from counting new branches after we blocked off a single leaf
				// from adding to the branch count
				hasConnections = TRUE; //Says whether the leaf we're looking at has other connections feeding into it
				
				set_big_matrix(currVect, currVectElements);
				big_floyd(A, currVect, bigMod, &theCycle);
				
				//If we found a vector that isn't a part of a cycle
				if (tau(theCycle) != 0)
				{
					//If we found our first branch
					if (numOfBranches == 0)
					{
						//printf("First branch.\n");
						
						//Insert vector into leaf
						branches[0][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
						
						//Get root (which we already know since this is the only branch and we have the transient length)
						copy_BigIntMatrixT(currVect, tempVect);
						for (int t = 0; t < tau(theCycle); t += 1)
						{
							big_mat_mul(A, tempVect, tempVect2);
							modbm(tempVect2, bigMod);
							copy_BigIntMatrixT(tempVect2, tempVect);
						}
						branches[0][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, tempVect);
						
						//Get base root
						branches[0][2] = branches[0][1];
					}
					
					else
					{
						//Find the root of our current vector (acting as a leaf)
						copy_BigIntMatrixT(currVect, tempRoot);
						for (int t = 0; t < tau(theCycle); t += 1)
						{
							big_mat_mul(A, tempRoot, tempVect);
							modbm(tempVect, bigMod);
							copy_BigIntMatrixT(tempVect, tempRoot);
						}
						
						//Check to see if we've already found this root
						isNewRoot = TRUE;
						for (int r = 0; r < numOfBranches; r += 1)
							if (compare_BigIntMatrixT(branches[r][2], tempRoot))
							{
								isNewRoot = FALSE;
								break;
							}
							
						//Easy case. We just create a new branch like before
						if (isNewRoot)
						{
							/*printf("Has a new root:\n");
							printbm(tempRoot);
							printf("\n");*/
							
							branches[numOfBranches][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
							branches[numOfBranches][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, tempRoot);
							branches[numOfBranches][2] = branches[numOfBranches][1];
						}
						
						//Now, we search through the existing branches to find where currVect grows from
						else
						{
							foundGrowth = FALSE;   //Says whether we found a place to attach currVect onto
							vectInsideVine = TRUE; //Says whether the vect is intrinsically inside another branch, meaning we can ignore it
							
							copy_BigIntMatrixT(currVect, tempVect);
							while (!foundGrowth)
							{
								for (int twig = 0; twig < numOfBranches; twig += 1)
								{
									//If we found a branch with the same base root
									if (compare_BigIntMatrixT(branches[twig][2], tempRoot))
									{
										//Is tempVect the root of our current branch?
										if (compare_BigIntMatrixT(tempVect, branches[twig][1]))
										{
											foundGrowth = TRUE;
											if (!vectInsideVine)
											{
												/*printf("Vector connects to base of branch:\n");
												printbm(branches[twig][1]);
												printf("\n");*/
												
												//Create new branch which feeds into the root of our current one
												branches[numOfBranches][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
												branches[numOfBranches][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][1]);
												branches[numOfBranches][2] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][2]);
											}
											else
											{
												//printf("Vector is inside the branch.\n");
												hasConnections = FALSE;
											}
											
											break; //Don't need to do anything if our vector is inside the branch
										}
										
										//Search through entire branch, checking if our tempVect can be found
										copy_BigIntMatrixT(branches[twig][0], tempVect2);
										while (!compare_BigIntMatrixT(tempVect2, branches[twig][1]))
										{
											//If we found where our vector feeds into on the current branch
											if (compare_BigIntMatrixT(tempVect2, tempVect))
											{
												foundGrowth = TRUE;
												
												//If our currVect connects to the leaf
												if (compare_BigIntMatrixT(tempVect2, branches[twig][0]))
												{
													hasConnections = FALSE;
													
													//Check to see if anything else connects to the leaf
													for (int br = 0; br < numOfBranches; br += 1)
													{
														if (compare_BigIntMatrixT(branches[br][1], tempVect2))
														{
															hasConnections = TRUE;
															break;
														}
													}
													
													//If nothing else connects to the leaf, extend the branch so that
													// currVect becomes the leaf
													if (!hasConnections)
													{
														/*printf("Vector connects to a single leaf:\n");
														printbm(branches[twig][0]);
														printf("\n");*/
														branches[twig][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
														break;
													}
												}
												
												//If currVect connects somehwere in the middle of the branch,
												// and isn't INSIDE the branch
												if (!vectInsideVine)
												{
													/*printf("Vector connects to branch at:\n");
													printbm(tempVect2);
													printf("\n");*/
													
													//New branch with currVect
													branches[numOfBranches][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, currVect);
													branches[numOfBranches][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, tempVect);
													branches[numOfBranches][2] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][2]);
													
													numOfBranches += 1;
													branches = realloc(branches, (numOfBranches+1)*sizeof(BigIntMatrixTP*));
													branches[numOfBranches] = malloc(3*sizeof(BigIntMatrixTP));
													
													//Bottom split of old branch
													branches[numOfBranches][0] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[numOfBranches-1][1]);
													branches[numOfBranches][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][1]);
													branches[numOfBranches][2] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[twig][2]);
													
													//Top split
													branches[twig][1] = BigIntMatrixT_catalogue_get(&vectCatalogue, &sizeOfCatalogue, branches[numOfBranches][0]);
													
													//Split old branch into two branches
													//We want our new branch's root to be a leaf for one of the split parts
												}
												else
												{
													//printf("Vector is inside the branch.\n");
													hasConnections = FALSE;
												}
												
												break;
											}
											
											//Iterate the vector in the branch we're checking
											big_mat_mul(A, tempVect2, tempVect3);
											modbm(tempVect3, bigMod);
											copy_BigIntMatrixT(tempVect3, tempVect2);
										}
									}
									
									if (foundGrowth)
										break;
								}
								
								//Iterate our vector, check again
								big_mat_mul(A, tempVect, tempVect2);
								modbm(tempVect2, bigMod);
								copy_BigIntMatrixT(tempVect2, tempVect);
								vectInsideVine = FALSE;
							}
						}
					}
					
					//Increase our branch count
					if (hasConnections)
					{
						numOfBranches += 1;
						branches = realloc(branches, (numOfBranches+1)*sizeof(BigIntMatrixTP*));
						branches[numOfBranches] = malloc(3*sizeof(BigIntMatrixTP));
					}
				}
			}
			while (!increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
			
			//Print out our branches
			/*
			printf("\n~~~~~Branches:~~~~~\n");
			for (int br = 0; br < numOfBranches; br += 1)
			{
				printf("Leaf:\n");
				printbm_row(branches[br][0]);
				printf("\nRoot:\n");
				printbm_row(branches[br][1]);
				printf("\nBase root:\n");
				printbm_row(branches[br][2]);
				printf("\n~~~~~~~~~~~~~~~~~~~\n");
				getchar();
			}
			*/
			
			//Now, we simply figure out which leaves are "true leaves",
			// that is, leaves that aren't the roots of any branch
			printf("Leaves of the transient tree:\n");
			for (int maybeLeaf = 0; maybeLeaf < numOfBranches; maybeLeaf += 1)
			{
				isTopLeaf = TRUE;
				for (int rootCheck = 0; rootCheck < numOfBranches; rootCheck += 1)
					if (compare_BigIntMatrixT(branches[maybeLeaf][0], branches[rootCheck][1]))
					{
						isTopLeaf = FALSE;
						break;
					}
					
				if (isTopLeaf)
				{
					if (vectorType == row)
					{
						printbm_row(branches[maybeLeaf][0]);
						if (outputFile != NULL)
							fprintbm_row(outputFile, branches[maybeLeaf][0]);
					}
					else
					{
						printbm(branches[maybeLeaf][0]);
						if (outputFile != NULL)
							fprintbm(outputFile, branches[maybeLeaf][0]);
					}
					
					printf("\n");
					if (outputFile != NULL)
						fprintf(outputFile, "\n");
				}
			}
			
			
			for (int i = 0; i <= numOfBranches; i += 1)
			{				
				FREE(branches[i]);
			}
			FREE(branches);
			
			for (int i = 0; i < sizeOfCatalogue; i += 1)
				vectCatalogue[i] = free_BigIntMatrixT(vectCatalogue[i]);
			FREE(vectCatalogue);
			
			free_BigIntT_array(currVectElements, big_rows(A), 1);
			
			bigMod = free_BigIntT(bigMod);
			
			A         = free_BigIntMatrixT(A);
			I         = free_BigIntMatrixT(I);
			currVect  = free_BigIntMatrixT(currVect);
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect);
			tempVect3 = free_BigIntMatrixT(tempVect);
			tempRoot  = free_BigIntMatrixT(tempRoot);
			theCycle  = free_CycleInfoT(theCycle);
			
			FREE(fileName);
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
		}
		
		
		//If the user wants to calculate each unique cycle space
		// within a given LCA system
		else if (! strcmp(argv[1], "orbitspaces"))
		{
			bool calcMinPoly = FALSE; //Do we calculate min polys for each cyclespace?
			bool newCyclespace;
			
			BigIntTP bigMod;
			BigIntTP** currVectElements;
			
			BigIntMatrixTP A;
			BigIntMatrixTP currVect;
			BigIntMatrixTP tempVect;
			BigIntMatrixTP tempVect2;
			
			BigIntTP** tempSetCyclespace;
			BigIntMatrixTP tempCyclespace;
			BigIntMatrixTP tempCyclespace2;
			
			int numOfCyclespaces = 0;
			BigIntMatrixTP* cyclespaces = NULL;

			BigPolyTP   charaPoly   = NULL;
			BigPolyTP** minpolys    = NULL;
			
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			//If the user wants to turn on minimal polynomial calculation
			if (oargv[3][0] != '\0')
				if (! strcmp(oargv[3], "TRUE"))
					calcMinPoly = TRUE;
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from config file.\n");
				bigMod = free_BigIntT(bigMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (big_rows(A) != big_cols(A))
			{
				fprintf(stderr, "Given update matrix is not square.\n");
				A = free_BigIntMatrixT(A);
				bigMod = free_BigIntT(bigMod);
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			if (calcMinPoly)
				charaPoly = chara_poly(A, bigMod);
			
			printf("Matrix:\n");
			printbm(A);
			printf("Modulus: ");
			printi(bigMod);
			printf("\n\n");
			
			currVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			
			tempCyclespace  = new_BigIntMatrixT(big_rows(A), big_cols(A));
			tempCyclespace2 = new_BigIntMatrixT(big_rows(A), big_cols(A));
			tempSetCyclespace = new_BigIntT_array(big_rows(A), big_cols(A));
			
			do
			{
				set_big_matrix(currVect, currVectElements);
				copy_BigIntMatrixT(currVect, tempVect);
				clear_BigIntT_array(tempSetCyclespace, big_rows(A), big_cols(A));
				
				//Calculate the max number of possible linearly-independent vectors
				// for the cyclespace, and put them in tempSetCyclespace
				for (int entry = 0; entry < big_rows(A); entry += 1)
					copy_BigIntT(currVectElements[entry][0], tempSetCyclespace[0][entry]);
				
				for (int v = 1; v < big_rows(A); v += 1)
				{
					big_mat_mul(A, tempVect, tempVect2);
					copy_BigIntMatrixT(tempVect2, tempVect);
					modbm(tempVect, bigMod);
					
					for (int entry = 0; entry < big_rows(A); entry += 1)
						copy_BigIntT(big_element(tempVect, entry, 0), tempSetCyclespace[v][entry]);
				}
				
				/*printbm_row(currVect);
				printf("\n"); */
				set_big_matrix(tempCyclespace, tempSetCyclespace);
				/*printbm(tempCyclespace);
				printf("~~~\n");*/
				big_row_echelon(tempCyclespace, bigMod, tempCyclespace2, NULL);
				big_reduced_row_echelon(tempCyclespace2, bigMod, tempCyclespace, NULL);
				//printbm(tempCyclespace);
				
				//Now, search to see if we've already accounted for the current cyclespace
				newCyclespace = TRUE;
				for (int m = 0; m < numOfCyclespaces; m += 1)
					if (compare_BigIntMatrixT(tempCyclespace, cyclespaces[m]))
					{
						newCyclespace = FALSE;
						break;
					}
				
				if (newCyclespace)
				{
					numOfCyclespaces += 1;
					cyclespaces = realloc(cyclespaces, numOfCyclespaces*sizeof(BigIntMatrixTP));
					cyclespaces[numOfCyclespaces-1] = new_BigIntMatrixT(big_rows(A), big_cols(A));
					copy_BigIntMatrixT(tempCyclespace, cyclespaces[numOfCyclespaces-1]);
					
					//Calculate min poly if needed
					if (calcMinPoly)
					{
						minpolys = realloc(minpolys, numOfCyclespaces*sizeof(BigPolyTP*));
						minpolys[numOfCyclespaces-1] = min_poly(A, currVect, bigMod, charaPoly);
					}
				}
			}
			while (! increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigMod));
			
			
			//Now, we print out all the cyclespaces we've found
			printf("Unique cyclespaces found:\n");
			for (int m = 0; m < numOfCyclespaces; m += 1)
			{
				big_rowsp(cyclespaces[m]);
				
				//If needed, also print out the min poly associated with the cyclespace
				if (calcMinPoly)
				{
					printf(" : ");
					old_printpf(minpolys[m]);
					minpolys[m] = free_BigPolyT_factors(minpolys[m]);
				}
				
				printf("\n");
				cyclespaces[m] = free_BigIntMatrixT(cyclespaces[m]);
			}
			FREE(cyclespaces);
			if (minpolys != NULL)
			{
				FREE(minpolys);
			}
			
			bigMod = free_BigIntT(bigMod);
			
			charaPoly   = free_BigPolyT(charaPoly);
			
			currVectElements  = free_BigIntT_array(currVectElements, big_rows(A), 1);
			tempSetCyclespace = free_BigIntT_array(tempSetCyclespace, big_rows(A), big_cols(A));
			
			A         = free_BigIntMatrixT(A);
			currVect  = free_BigIntMatrixT(currVect);
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect2);
			
			tempCyclespace  = free_BigIntMatrixT(tempCyclespace);
			tempCyclespace2 = free_BigIntMatrixT(tempCyclespace2);
		}
		
		
		//If the user wants to use Floyd's Cycle Detection Algorithm
		else if (! strcmp(argv[1], "floyd"))
		{
			BigIntTP bigModulus;
			BigIntMatrixTP initial;
			BigIntMatrixTP update;
			CycleInfoTP coolCycle = NULL;
			
			//If the user specified a custom modulus
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], bigModulus, "Unable to read modulus from command line.");
			}
			
			//Use the modulus provided in the .config file
			else
			{
				SET_BIG_NUM(bigintmodstring, bigModulus, "Unable to read modulus from config file. Try giving the modulus as a CLI argument.");
			}
			
			//Now, get the matrices from the provided files
			initial = read_BigIntMatrixT(initialfilepath);
			update  = read_BigIntMatrixT(updatefilepath);
			
			if ((initial == NULL) || (update == NULL))
			{
				fprintf(stderr, "Unable to read .matrix files.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			else if ((big_rows(initial) != big_cols(update)) ||
			         (big_rows(initial) != big_rows(update)))
		  {
			  fprintf(stderr, "Given matrices have inappropriate dimensions for iterating.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
		  }
			
			//Get cycle info
			big_floyd(update, initial, bigModulus, &coolCycle);
			
			if (coolCycle == NULL)
				printf("Update matrix is not a square matrix.\n");
			
			else
			{
				printcycle(coolCycle, vectorType);
				if ((vectorType == row) && (big_rows(initial) != big_cols(initial)))
					printf("\n");
			}
			
			bigModulus = free_BigIntT(bigModulus);
			initial    = free_BigIntMatrixT(initial);
			update     = free_BigIntMatrixT(update);
			coolCycle  = free_CycleInfoT(coolCycle);
		}
		
		
		//Generate some basic rotation matrices for the given modulus
		else if (! strcmp(argv[1], "rots"))
		{
			PROHIBIT_UNIXFLAGS
			
			bool has30 = TRUE;
			bool has45 = TRUE;
			
			int inv2;
			int sqrt3;
			
			//If user provided modulus on CLI
			if (argc > 2)
			{
				modulus = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read modulus from command line.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
			}
			
			if ((num_inverse(2, modulus) == -1))
			{
				printf("45degrot and 30degrot don't exist mod %d.\n", modulus);
				has45 = FALSE;
				has30 = FALSE;
			}
			
			else
			{
				if (square_root(2, modulus) == -1)
				{
					printf("45degrot doesn't exist mod %d.\n", modulus);
					has45 = FALSE;
				}
				
				if (square_root(3, modulus) == -1)
				{
					printf("30degrot does not exist mod %d.\n", modulus);
					has30 = FALSE;
				}
			}
			
			
			if (has45)
			{
				printf("45degrot mod %d:\n", modulus);
				
				printf("%2d %2d\n",
				(num_inverse(2, modulus)*square_root(2, modulus)) % modulus,
				modulus - ((num_inverse(2, modulus)*square_root(2, modulus)) % modulus));
				
				printf("%2d %2d\n",
				(num_inverse(2, modulus)*square_root(2, modulus)) % modulus,
				(num_inverse(2, modulus)*square_root(2, modulus)) % modulus);
			}
			
			if (has30)
			{
				inv2  = (modulus - num_inverse(2, modulus)) % modulus;
				sqrt3 = square_root(3, modulus);
				
				printf("30degrot mod %d:\n", modulus);
				
				printf("%2d %2d\n",
				(inv2*sqrt3) % modulus,
				((modulus - inv2) % modulus));
				
				printf("%2d %2d\n",
				inv2, 
				(inv2*sqrt3) % modulus);
			}
		}
		
		
		//If we want to look for matrices that have column vectors
		// that cycle with a cycle length less than the matrix itself
		else if (! strcmp(argv[1], "cycmatsearch"))
		{
			PROHIBIT_UNIXFLAGS
			
			//If user didn't provide enough arguments for the tool, explain how to use it
			if (argc < 6)
			{
				printf(ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "resume size maxmod cycles..." ANSI_COLOR_RESET \
				": Searches for a matrix whose column vectors cycle with cycle lengths less than the matrix itself.\n");
				printf(" - " ANSI_COLOR_CYAN "resume" ANSI_COLOR_RESET \
				": Says whether to resume computation at a particular matrix.\n");
				printf(" - " ANSI_COLOR_CYAN "size" ANSI_COLOR_RESET \
				": Tells what size matrix to use.\n");
				printf(" - " ANSI_COLOR_CYAN "maxmod" ANSI_COLOR_RESET \
				": Tells which modulus to stop searching at.\n");
				printf(" - " ANSI_COLOR_CYAN "cycles..." ANSI_COLOR_RESET \
				": A list of cycle lengths for their respective vectors.\n\n");
				
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			
			FILE* textOutput; 
			char* textOutputName;
			
			BigIntTP maxMod;  //What's the largest modulus we should search to?
			BigIntTP currMod; //What modulus are we currently checking?
			
			BigIntTP temp;
			BigIntTP temp2;
			BigIntTP tempModCounter;
			BigIntTP tempPercentCounter;
			
			BigIntTP progressUpdateElement; //For keeping track of progress through a modulus
			int progressRow = 1;
			int progressCol = 1; //Used to decide which element in the matrix to look at for progress updates
			bool printProgress = TRUE;
			
			BigIntTP** currMatElements; //Holds matrix numbers so we can set the matrix easily
			
			CycleInfoTP theCycle; //Used to get matrices in a cycle
			
			BigIntMatrixTP zeroMat; //The zero matrix
			BigIntMatrixTP currMat; //Holds the matrix we're currently testing for cyclic vectors
			BigIntMatrixTP tempMat; //Holds results of matrix multiplication
			BigIntMatrixTP tempMat2;
			
			bool checkedAllMatrices; //Says whether we've exhausted all matrices under a given modulus
			bool matrixIsCyclic;     //Says whether the matrix we tested has cyclic column vectors or not
			int i, j; //Used for for-loops so we don't have to keep declaring new variables
			int size; //The size of our matrix to use
			int currIteration;
			int nextIteration; //Used to keep track of which iteration of the matrix to check next
			int cycleLCM = 1;
			
			int start[] = {2}; //Modulus to start with
			
			int* colVectCycles; //Holds the different cycle lengths 
			
			size = (int)strtol(argv[3], &tempStr, 10);
			
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read matrix size from command line.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			currMatElements = new_BigIntT_array(size, size);
			
			//If we need to resume computation from a specific matrix
			if (! strcmp(argv[2], "TRUE"))
			{
				currMat = read_BigIntMatrixT(resumefilepath);
				SET_BIG_NUM(bigintmodstring, currMod, "Unable to read resume modulus from config file.");
				
				if (currMat == NULL)
				{
					fprintf(stderr, "Unable to read resume matrix from config file.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
				
				if ((big_rows(currMat) != size) || (big_cols(currMat) != size))
				{
					fprintf(stderr, "Given resume matrix is of incorrect size (%dx%d instead of %dx%d).\n",
					big_rows(currMat), big_cols(currMat), size, size);
					returnvalue = EXIT_SUCCESS;
					goto FREE_VARIABLES;
				}
				
				//Now, set the values for currMatElements
				for (i = 0; i < size; i += 1)
					for (j = 0; j < size; j += 1)
						copy_BigIntT(big_element(currMat, i, j), currMatElements[i][j]);
			}
			
			else
			{
				//Start looking for matrices at mod 2 from the zero matrix
				currMod = new_BigIntT(start, 1);
				currMat = new_BigIntMatrixT(size, size);
			}
			
			//If we couldn't read the modulus
			SET_BIG_NUM(argv[4], maxMod, "Unable to read modulus from command line.");
			
			//If the user didn't provide enough cycles for the given matrix size
			if (argc < 5 + size)
			{
				printf("Too few cycles passed for given matrix size.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Get column vector cycles
			colVectCycles = malloc(size*sizeof(int));
			for (i = 0; i < size; i += 1)
			{
				colVectCycles[i] = (int)strtol(argv[5+i], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read cycle lengths from command line.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
			}
			
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			
			tempModCounter        = empty_BigIntT(1);
			tempPercentCounter    = empty_BigIntT(1);
			progressUpdateElement = empty_BigIntT(1);
			
			//Find the LCM of our cycles
			for (i = 0; i < size; i += 1)
				cycleLCM = LCM(cycleLCM, colVectCycles[i]);
			
			zeroMat  = new_BigIntMatrixT(size, size);
			tempMat  = new_BigIntMatrixT(size, size);
			tempMat2 = new_BigIntMatrixT(size, size);
			
			//Initialising a reusable CycleInfoTP to save time on
			// memory allocations.
			theCycle = new_CycleInfoT();
			
			//Forming name for text file output.
			textOutputName = malloc(MAXSTRLEN*sizeof(char));
			textOutputName[0] = '\0';
			ADD_FILE_PREFIX(textOutputName);
				
			strcat(textOutputName, argv[1]); //cycmatsearch
			strcat(textOutputName, " ");
			strcat(textOutputName, argv[3]); //size
			strcat(textOutputName, " ");
			strcat(textOutputName, argv[4]); //maxmod
			strcat(textOutputName, " ");
			for (i = 0; i < size-1; i += 1)
			{
				strcat(textOutputName, argv[5+i]);
				strcat(textOutputName, " ");
			}
			strcat(textOutputName, argv[size+4]);
			strcat(textOutputName, ".txt");
			
			textOutput = fopen(textOutputName, "w");
			if (textOutput == NULL)
			{
				fprintf(stderr, "Unable to create/open %s for writing.\n", textOutputName);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Prepare to monitor whether the progress element in our matrix increased or not
			copy_BigIntT(currMatElements[progressRow][progressCol], progressUpdateElement);

			
			//Search all moduli until we get to the specified limit
			while (compare_BigIntT(currMod, maxMod) <= 0)
			{
				printf("Currently searching mod ");
				printi(currMod);
				printf("...\n");
				
				fprintf(textOutput, "~~~Mod ");
				fprinti(textOutput, currMod);
				fprintf(textOutput, "~~~\n");

				//Search through all matrices under the current modulus
				checkedAllMatrices = FALSE;
				while (!checkedAllMatrices)
				{
					//Getting all the relevant information we need to check for vector cycles
					set_big_matrix(currMat, currMatElements);
					big_floyd(currMat, currMat, currMod, &theCycle);
					copy_BigIntMatrixT(rep(theCycle), tempMat);
					
					if ((printProgress) && 
					    (compare_BigIntT(currMatElements[progressRow][progressCol], progressUpdateElement) != 0))
					{
						//Calculate proper fraction to show for percentage
						copy_BigIntT(one, tempModCounter);
						copy_BigIntT(zero, tempPercentCounter);
						for (int r = progressRow; r < size; r += 1)
						{
							for (int c = 0; c < size; c += 1)
							{
								//Making sure we start from the correct element in the matrix
								if ((r == progressRow) && (c == 0))
									c = progressCol;
								
								multiply_BigIntT(currMatElements[r][c], tempModCounter, temp);
								add_BigIntT(tempPercentCounter, temp, temp2);
								copy_BigIntT(temp2, tempPercentCounter);
								
								//Increment modCounter
								multiply_BigIntT(tempModCounter, currMod, temp);
								copy_BigIntT(temp, tempModCounter);
							}
						}
						
						printi(tempPercentCounter);
						printf(" / ");
						printi(tempModCounter);
						printf(" searched...\n");
						
						copy_BigIntT(currMatElements[progressRow][progressCol], progressUpdateElement);
					}
					
					currIteration = 0;
					
					//Prevent matrices that don't have the correct cycle length
					// from being printed
					matrixIsCyclic = FALSE; 
				
					//Check to see if the vectors actually cycle how we want them to
					//The matrix's cycle must be the LCM of the cycle lengths if
					// the vectors are to behave how we want them to.
					while ((omega(theCycle) == cycleLCM))
					{
						//Assume the matrix is cyclic until proven otherwise
						matrixIsCyclic = TRUE;
						
						//Get a value of nextIteration that's greater than before,
						// then work from there
						for (i = 0; i < size; i += 1)
							if (colVectCycles[i] > currIteration)
							{
								nextIteration = colVectCycles[i];
								break;
							}
						
						//Looking for the next smallest iteration we need to check
						for (i = 1; i < size; i += 1)
							if ((colVectCycles[i] < nextIteration) && (colVectCycles[i] > currIteration))
								nextIteration = colVectCycles[i];
							
						//If we've checked all the iterations we need to check
						if (nextIteration == currIteration)
							break;
								
						//Iterate matrix to the correct iteration
						for (i = 0; i < nextIteration - currIteration; i += 1)
						{
							big_mat_mul(currMat, tempMat, tempMat2);
							modbm(tempMat2, currMod);
							copy_BigIntMatrixT(tempMat2, tempMat);
						}
						
						//Check for repeated vectors
						for (i = 0; i < size; i += 1)
							if (colVectCycles[i] == nextIteration)
								if (compare_BigIntMatrixT_cols(tempMat, rep(theCycle), i) != 1)
								{
									matrixIsCyclic = FALSE;
									break;
								}
								
							
						currIteration = nextIteration;
						
						//If we've determined the matrix doesn't work
						if (!matrixIsCyclic)
							break;
					}
					
					//If we found a matrix that cycles how we want, print it
					if ((matrixIsCyclic))
					{
						printf("Mod ");
						printi(currMod);
						printf("\n");
						printbm(currMat);
						printf("\n");
						
						fprintbm_nopad(textOutput, currMat);
						fprintf(textOutput, "\n");
						
						if (fclose(textOutput) == EOF)
						{
							fprintf(stderr, "Error saving file %s.\n", textOutputName);
							checkedAllMatrices = FALSE;
							break;
						}
						
						textOutput = fopen(textOutputName, "a");
						if (textOutput == NULL)
						{
							fprintf(stderr, "Error reopening file %s after saving.\n", textOutputName);
							checkedAllMatrices = FALSE;
							break;
						}
					}

					//Iterate to the next matrix
					checkedAllMatrices = increment_BigIntT_array(currMatElements, size, size, one, currMod);
				}
				
				//Increment modulus
				add_BigIntT(one, currMod, temp);
				copy_BigIntT(temp, currMod);
			}
			
			if (fclose(textOutput) == EOF)
				fprintf(stderr, "Unable to close %s.\n", textOutputName);
			
			//Freeing memory
			maxMod  = free_BigIntT(maxMod);
			currMod = free_BigIntT(currMod);
			temp    = free_BigIntT(temp);
			temp2   = free_BigIntT(temp2);
			
			tempModCounter = free_BigIntT(tempModCounter);
			tempPercentCounter = free_BigIntT(tempPercentCounter);
			progressUpdateElement = free_BigIntT(progressUpdateElement);
			
			currMat  = free_BigIntMatrixT(currMat);
			zeroMat  = free_BigIntMatrixT(zeroMat);
			tempMat  = free_BigIntMatrixT(tempMat);
			tempMat2 = free_BigIntMatrixT(tempMat2);
			
			theCycle = free_CycleInfoT(theCycle);
			
			FREE(colVectCycles);
			FREE(textOutputName);
			currMatElements = free_BigIntT_array(currMatElements, size, size);
			
			printf("Finished searching.\n");
		}
		
		
		//If we want to compute some "cycle converting matrices"
		else if (!strcmp(argv[1], "cycconvmat"))
		{
			BigIntTP bigMod, from, to;
			BigIntMatrixTP sumMat, A;
			
			//If the user didn't supply the required arguments "from" and "to"
			if ((oargv[2][0] == '\0') || (oargv[3][0] == '\0'))
			{
				printf("Usage: " ANSI_COLOR_YELLOW "cycconvmat " ANSI_COLOR_CYAN "from to [mod]" ANSI_COLOR_RESET "\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Check to see if user provided a custom modulus
			if (oargv[4][0] != '\0')
			{
				SET_BIG_NUM(oargv[4], bigMod, "Unable to read modulus from command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			SET_BIG_NUM(oargv[2], from, "Unable to read starting cycle length from command line.");
			SET_BIG_NUM(oargv[3], to, "Unable to read conversion cycle length from command line.");
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix provided in config file.\n");
				bigMod = free_BigIntT(bigMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}

			sumMat = new_BigIntMatrixT(big_rows(A), big_rows(A));
			ccm(A, sumMat, from, to, bigMod);
			
			//Now, print out our resulting "cycle converting matrix"
			printf("Matrix:\n");
			printbm(A);
			printf("Modulus: ");
			printi(bigMod);
			printf("\n\nCCM from ");
			printi(from);
			printf(" to ");
			printi(to);
			printf(":\n");
			printbm(sumMat);
			
			from   = free_BigIntT(from);
			to     = free_BigIntT(to);
			bigMod = free_BigIntT(bigMod);
			
			A      = free_BigIntMatrixT(A);
			sumMat = free_BigIntMatrixT(sumMat);
		}
		
		
		//If we want to look for matrices that have CCMs equal to zero
		else if (! strcmp(argv[1], "ccmzerosearch"))
		{
			//Checking to see if the user provided all the required arguments
			if (oargv[3][0] == '\0')
			{
				printf("Usage: " ANSI_COLOR_YELLOW "ccmzerosearch " ANSI_COLOR_CYAN "resume size [mod]" ANSI_COLOR_RESET "\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			BigIntTP   bigMod, temp, temp2;
			BigIntTP   prevLastElement; //Holds the last element in the matrix; for debugging
			BigIntTP*  cycleLengthFactors;
			BigIntTP** currMatElements;
			BigIntTP** currVectElements;
			
			BigIntMatrixTP currMat, currVect, tempCCM, zeroMat;
			BigIntMatrixTP resumeMat = NULL; //The matrix to resume at if specified
			
			CycleInfoTP theCycle = NULL;
			
			int matSize;
			int indexCounter; //For freeing and printing
			
			//Holds the cycle length of each matrix and vector we check
			int cycleLengthArray[1] = {0};
			BigIntTP bigOmega, bigVectOmega;
			
			//Holds the number of vectors with each possible cycle length
			int* cycleLengthCounts;
			
			bool checkedAllMatrices = FALSE;
			bool checkedAllVectors  = FALSE;
			bool hasAllZeros        = FALSE;
			
			bool hasInterestingCycles = FALSE;
			
			bool debug = TRUE;
			int matrixCounter = 0;
			
			char* outputfilename;
			FILE* outputFile;
			BigPolyTP  charaPoly;
			BigPolyTP* charaPolyFactors = NULL;
			
			//Checking to see if the user provided a modulus on the command line
			if (oargv[4][0] != '\0')
			{
				SET_BIG_NUM(oargv[4], bigMod, "Unable to read modulus on command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			matSize = (int)strtol(oargv[3], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Invalid matrix size provided on command line.\n");
				bigMod = free_BigIntT(bigMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//If the user wants to use the resume matrix
			if (!strcmp(oargv[2], "TRUE"))
			{
				resumeMat = read_BigIntMatrixT(resumefilepath);
				if (resumeMat == NULL)
				{
					fprintf(stderr, "Unable to read resume matrix specified in config file. Computation " \
					"will start from the zero matrix.\n");
				}
			}
			
			//Create output name
			outputfilename = malloc(MAXSTRLEN*sizeof(char));
			outputfilename[0] = '\0';
			
			//Maybe I should make this into a macro
			ADD_FILE_PREFIX(outputfilename);
			
			strcat(outputfilename, "ccmzerosearch "); //ccmzerosearch
			strcat(outputfilename, oargv[3]); //matrix size
			strcat(outputfilename, " ");
			append_BigIntT(outputfilename, bigMod);
			
			strcat(outputfilename, ".txt");
			outputFile = fopen(outputfilename, "w");
			if (outputFile == NULL)
				fprintf(stderr, "Unable to save found matrices. Continuing without saving...\n");
			else
				printf("Found matrices will be saved at %s\n", outputfilename);
			
			//Initialising elements for matrices and vectors
			currMatElements  = new_BigIntT_array(matSize, matSize);
			currVectElements = new_BigIntT_array(matSize, 1);
			
			currMat  = new_BigIntMatrixT(matSize, matSize);
			zeroMat  = new_BigIntMatrixT(matSize, matSize);
			currVect = new_BigIntMatrixT(matSize, 1);
			
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			prevLastElement = empty_BigIntT(1);
			
			//Allowing user to resume computation
			if (resumeMat != NULL)
				for (int i = 0; i < matSize; i += 1)
					for (int j = 0; j < matSize; j += 1)
						copy_BigIntT(big_element(resumeMat, i, j), currMatElements[i][j]);
			
			//Loop until we're checked every matrix under the given modulus
			while (!checkedAllMatrices)
			{
				set_big_matrix(currMat, currMatElements);
				big_floyd(currMat, currMat, bigMod, &theCycle);
				
				if (debug)
				{
					if (compare_BigIntT(currMatElements[matSize-1][matSize-2], prevLastElement) != 0)
					{
						printf("Current matrix:\n");
						printbm(currMat);
					}
				}
				
				//This program currently assumes that the cycle length of the
				// matrix will be less than MAXBUNCH
				cycleLengthArray[0] = omega(theCycle);
				bigOmega = new_BigIntT(cycleLengthArray, 1);
				cycleLengthFactors = divisors_of_BigIntT(bigOmega);
				
				//Check each factor to see what its CCM looks like
				indexCounter = 1;
				copy_BigIntT(zero, temp);
				hasAllZeros = TRUE;
				while (compare_BigIntT(temp, cycleLengthFactors[0]) < 0)
				{
					tempCCM = new_BigIntMatrixT(matSize, matSize);
					ccm(currMat, tempCCM, bigOmega, cycleLengthFactors[indexCounter], bigMod);
					
					//If the CCM is not the zero matrix, we can stop looking
					if (!compare_BigIntMatrixT(tempCCM, zeroMat))
					{
						tempCCM = free_BigIntMatrixT(tempCCM);
						hasAllZeros = FALSE;
						break;
					}
					
					indexCounter += 1;
					add_BigIntT(temp, one, temp2);
					copy_BigIntT(temp2, temp);
					tempCCM = free_BigIntMatrixT(tempCCM);
				}
				
				//Calculate cycle length counts so we can determine if the CCM
				// zeros are justified or not
				cycleLengthCounts = calloc(indexCounter, sizeof(int));
				
				checkedAllVectors = FALSE;
				while (!checkedAllVectors)
				{
					set_big_matrix(currVect, currVectElements);
					big_floyd(currMat, currVect, bigMod, &theCycle);

					cycleLengthArray[0] = omega(theCycle);
					bigVectOmega = new_BigIntT(cycleLengthArray, 1);
					
					//First, check if vector cycle length is equal to matrix
					//This prevents possible index errors below
					if (compare_BigIntT(bigVectOmega, bigOmega) == 0)
						cycleLengthCounts[indexCounter-1] += 1;
					
					else
					{
						//Find which cycle length group this vector belongs to
						for (int w = 1; w <= indexCounter-1; w += 1)
							if (compare_BigIntT(cycleLengthFactors[w], bigVectOmega) == 0)
							{
								cycleLengthCounts[w-1] += 1;
								break;
							}
					}

					//Iterate to next vector
					checkedAllVectors = increment_BigIntT_array(currVectElements, matSize, 1, one, bigMod);
				}
				
				
				//Now, check to see if any "interesting cycle lengths" exist for the matrix
				// This is, cycle lengths that aren't maximal or zero
				hasInterestingCycles = FALSE;
				for (int c = 0; c < indexCounter-1; c += 1)
				{
					//We only care about fixed points if there's more than just the zero vector
					if (c == 0)
					{
						if (cycleLengthCounts[0] > 1)
						{
							hasInterestingCycles = TRUE;
							break;
						}
					}
					
					else
					{
						if (cycleLengthCounts[c] > 0)
						{
							hasInterestingCycles = TRUE;
							break;
						}
					}
				}
				
				//Print out the matrix we found, as well as any relevant information regarding it
				if ((hasAllZeros) && (compare_BigIntT(one, bigOmega) != 0) && (hasInterestingCycles))
				{
					matrixCounter += 1;
					printbm(currMat);
					printf("Cycle length counts: ");
					
					if (outputFile != NULL)
					{
						fprintbm_nopad(outputFile, currMat);
						fprintf(outputFile, "Cycle length counts: ");
					}
					
					indexCounter = 1;
					copy_BigIntT(zero, temp);
					while (compare_BigIntT(temp, cycleLengthFactors[0]) < 0)
					{
						printi(cycleLengthFactors[indexCounter]);
						printf(" (%d), ", cycleLengthCounts[indexCounter-1]);
						
						if (outputFile != NULL)
						{
							fprinti(outputFile, cycleLengthFactors[indexCounter]);
							fprintf(outputFile, " (%d), ", cycleLengthCounts[indexCounter-1]);
						}
						
						indexCounter += 1;
						add_BigIntT(temp, one, temp2);
						copy_BigIntT(temp2, temp);
					}
					printi(bigOmega);
					printf(" (%d)\n", cycleLengthCounts[indexCounter-1]);
					
					charaPoly = chara_poly(currMat, bigMod);
					charaPolyFactors = old_factor_BigPolyT(charaPoly, bigMod);
					printf("Chara poly: ");
					printp(charaPoly);
					printf("\nFactored chara poly: ");
					old_printpf(charaPolyFactors);
					printf("\n\n");
					
					if (outputFile != NULL)
					{
						fprinti(outputFile, bigOmega);
						fprintf(outputFile, " (%d)\n", cycleLengthCounts[indexCounter-1]);
						
						fprintf(outputFile, "Chara poly: ");
						old_fprintpf(outputFile, charaPolyFactors);
						fprintf(outputFile, "\n\n");
						
						//Saving file
						if (fclose(outputFile) == EOF)
						{
							fprintf(stderr, "Unable to save data. Continuing without saving...\n");
							outputFile = NULL;
						}
						
						else
						{
							outputFile = fopen(outputfilename, "a");
							if (outputFile == NULL)
								fprintf(stderr, "Unable to save data. Continuing without saving...\n");
						}
					}
					
					charaPoly = free_BigPolyT(charaPoly);
					charaPolyFactors = free_BigPolyT_factors(charaPolyFactors);
				}
				
				bigOmega = free_BigIntT(bigOmega);
				bigVectOmega = free_BigIntT(bigVectOmega);
				FREE(cycleLengthCounts);
				
				//Freeing this array is a massive pain
				//I have to do it each loop since the factorisation matrix allocates new
				// arrays each time.
				//I should probably change this in the future for speed.
				indexCounter = 1;
				copy_BigIntT(zero, temp);
				while (compare_BigIntT(temp, cycleLengthFactors[0]) < 0)
				{
					cycleLengthFactors[indexCounter] = free_BigIntT(cycleLengthFactors[indexCounter]);
					indexCounter += 1;
					add_BigIntT(temp, one, temp2);
					copy_BigIntT(temp2, temp);
				}
				cycleLengthFactors[0] = free_BigIntT(cycleLengthFactors[0]);
				FREE(cycleLengthFactors);
				
				//Iterate to next matrix
				copy_BigIntT(currMatElements[matSize-1][matSize-2], prevLastElement);
				checkedAllMatrices = increment_BigIntT_array(currMatElements, matSize, matSize, one, bigMod);
			}
			
			if (debug)
				printf("Number of matrices with unjustified CCMs: %d\n", matrixCounter);
			
			
			currMatElements  = free_BigIntT_array(currMatElements, matSize, matSize);
			currVectElements = free_BigIntT_array(currVectElements, matSize, 1);
			
			currMat   = free_BigIntMatrixT(currMat);
			zeroMat   = free_BigIntMatrixT(zeroMat);
			currVect  = free_BigIntMatrixT(currVect);
			resumeMat = free_BigIntMatrixT(resumeMat);
			
			theCycle = free_CycleInfoT(theCycle);
			
			bigMod = free_BigIntT(bigMod);
			temp   = free_BigIntT(temp);
			temp2  = free_BigIntT(temp2);
			
			prevLastElement = free_BigIntT(prevLastElement);
			
			FREE(outputfilename);
			if (fclose(outputFile) == EOF)
				fprintf(stderr, "Unable to save data.\n");
		}
		
		
		//If the user wants to calculate every vector's 
		// "minimal" annihilating polynomial.
		else if (!strcmp(argv[1], "vectprops"))
		{
			if ((oargv[2][0] == '\0') || (oargv[3][0] == '\0'))
			{
				printf("Usage: " ANSI_COLOR_YELLOW "vectprops " ANSI_COLOR_CYAN "baseMod modPower [resume] [fileoutput]" ANSI_COLOR_RESET "\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			bool checkedAllVects = FALSE;
			
			int modPower, numOfAnnihPolys;
			BigIntTP baseMod, bigMod, temp;
			
			BigIntTP** currVectElements;
			
			BigIntMatrixTP currVect;
			BigIntMatrixTP resumeVect = NULL;
			
			BigPolyTP* generators;
			
			char* fileName   = NULL;
			FILE* outputFile = NULL;
			
			SET_BIG_NUM(oargv[2], baseMod, "Unable to read base of prime/prime-power modulus from command line.");
			
			modPower = (int)strtol(oargv[3], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read modulus exponent from command line.\n");
				baseMod = free_BigIntT(baseMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (big_rows(UPDATEMATRIX) != big_cols(UPDATEMATRIX))
			{
				fprintf(stderr, "Given update matrix is not square.\n");
				baseMod = free_BigIntT(baseMod);
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			currVectElements = new_BigIntT_array(big_rows(UPDATEMATRIX), 1);
			
			//If the user wants to resume computation from a
			// particular vector
			if (!strcmp(oargv[4], "TRUE"))
			{
				resumeVect = read_BigIntMatrixT(resumefilepath);
				if (resumeVect == NULL)
				{
					fprintf(stderr, "Unable to read resume vector from .config file. ");
					fprintf(stderr, "Beginning computation from the zero vector instead.\n");
				}
				
				else if ((big_rows(resumeVect) != big_rows(UPDATEMATRIX)) || 
								 (big_cols(resumeVect) != 1))
				{
					fprintf(stderr, "Resume vector is of incorrect dimensions (must be %d by 1). ", big_rows(UPDATEMATRIX));
					fprintf(stderr, "Beginning computation from the zero vector instead.\n");
				}
				
				else //Set currVect to resumeVect
					for (int elem = 0; elem < big_rows(UPDATEMATRIX); elem += 1)
						copy_BigIntT(big_element(resumeVect, elem, 0), currVectElements[elem][0]);
					
				resumeVect = free_BigIntMatrixT(resumeVect);
			}
	
			currVect = new_BigIntMatrixT(big_rows(UPDATEMATRIX), 1);
			
			temp = empty_BigIntT(1);
			bigMod = empty_BigIntT(1);
			copy_BigIntT(one, bigMod);
			
			//Compute bigMod
			for (int i = 0; i < modPower; i += 1)
			{
				multiply_BigIntT(baseMod, bigMod, temp);
				copy_BigIntT(temp, bigMod);
			}
			
			//If user wants file output
			if (!strcmp(oargv[5], "TRUE"))
			{
				fileName = malloc(MAXSTRLEN*sizeof(char));
				fileName[0] = '\0';
				
				ADD_FILE_PREFIX(fileName);
				
				strcat(fileName, "vectprops F");
				for (int x = 0; x < big_rows(UPDATEMATRIX); x += 1)
					for (int y = 0; y < big_cols(UPDATEMATRIX); y += 1)
						append_BigIntT(fileName, big_element(UPDATEMATRIX, x, y));
				strcat(fileName, " ");
				append_BigIntT(fileName, bigMod);
				strcat(fileName, ".txt");
				
				outputFile = fopen(fileName, "a");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to open file for saving output. Continuing without saving...\n");
			}
			
			//Print out the update matrix for my sanity
			printf("\nNote that degree 0 annihilating polynomials are not currently calculated.\n\n");
			printf("Modulus: ");
			printi(bigMod);
			printf("\nUpdate matrix: \n");
			printbm(UPDATEMATRIX);
			printf("\n");
			
			//Iterate through all vectors in the module
			while (!checkedAllVects)
			{
				set_big_matrix(currVect, currVectElements);
				
				//Get the generators for currVect's ideal of annihilating polynomials
				generators = ann_generators(UPDATEMATRIX, currVect, baseMod, modPower);
				numOfAnnihPolys = extract_bunch(constant(generators[0]), 0);
				
				//Basically, we're ignoring the zero vector, as the zero vector will be the
				// only one without non-degree 0 generators.
				if (numOfAnnihPolys > 0)
				{					
					printbm_row(currVect);
					printf("'s annideal generators: ");
					for (int i = 1; i <= numOfAnnihPolys; i += 1)
					{
						if (i != 1)
							printf(", ");
						printp(generators[i]);
					}
					printf("\n");
					
					if (outputFile != NULL)
					{
						fprintbm_row(outputFile, currVect);
						fprintf(outputFile, " : ");
						fprintp(outputFile, generators[1]);
						fprintf(outputFile, "\n");
						
						if (fclose(outputFile) == EOF)
							fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
						else
						{
							outputFile = fopen(fileName, "a");
							if (outputFile == NULL)
								fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
						}
					}
				}
				
				for (int i = 0; i <= numOfAnnihPolys; i += 1)
					generators[i] = free_BigPolyT(generators[i]);
				FREE(generators);
				
				//Increment currVect
				checkedAllVects = increment_BigIntT_array(currVectElements, big_rows(UPDATEMATRIX), 1, one, bigMod);
			}
			
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
				
			if (fileName != NULL)
			{
				FREE(fileName);
			}
			
			currVectElements = free_BigIntT_array(currVectElements, big_rows(UPDATEMATRIX), 1);
			
			currVect = free_BigIntMatrixT(currVect);
			
			temp    = free_BigIntT(temp);
			bigMod  = free_BigIntT(bigMod);
			baseMod = free_BigIntT(baseMod);
		}
		
		
		//If the user wants to calculate every annihilating polynomial of a vector under
		// an update matrix up to a certain degree
		else if (!strcmp(argv[1], "vectpolys"))
		{
			int modPower, numOfAnnihPolys;
			BigIntTP baseMod;
			
			BigPolyTP* generators;
			
			if ((oargv[2][0] == '\0') || (oargv[3][0] == '\0'))
			{
				printf("Usage: " ANSI_COLOR_YELLOW "vectpolys " ANSI_COLOR_CYAN "baseMod modPower " ANSI_COLOR_RESET "\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			SET_BIG_NUM(oargv[2], baseMod, "Unable to read base of modulus from command line.");
			
			modPower = (int)strtol(oargv[3], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read modulus exponent from command line.\n");
				baseMod = free_BigIntT(baseMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (big_rows(UPDATEMATRIX) != big_cols(UPDATEMATRIX))
			{
				fprintf(stderr, "Given update matrix is not square.\n");
				baseMod = free_BigIntT(baseMod);
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			if ((big_rows(INITIALMATRIX) != big_rows(UPDATEMATRIX)) || (big_cols(INITIALMATRIX) != 1))
			{
				fprintf(stderr, "Given vector's dimensions are invalid. The number of rows must be equal to the matrix's rows " \
				"and the number of columns must be 1.\n");
				baseMod = free_BigIntT(baseMod);
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			//Let's hope we don't have like a billion generators...
			generators = ann_generators(UPDATEMATRIX, INITIALMATRIX, baseMod, modPower);
			numOfAnnihPolys = extract_bunch(constant(generators[0]), 0);
			
			printf("\nNote that degree 0 annihilating polynomials are not currently calculated.\n\n");
			printf("Modulus = ");
			printi(baseMod);
			printf("^%d\nA = \n", modPower);
			printbm(UPDATEMATRIX);
			printf("\nGenerators for Ann(");
			printbm_row(INITIALMATRIX);
			printf(") under A:\n");
			for (int i = 1; i <= numOfAnnihPolys; i += 1)
			{
				printp(generators[i]);
				printf("\n");
			}
			
			for (int i = 0; i <= numOfAnnihPolys; i += 1)
				generators[i] = free_BigPolyT(generators[i]);
			FREE(generators);
			
			baseMod = free_BigIntT(baseMod);
		}
		
		
		//If the user wants to calculate the algebraic relations a given matrix satisfies
		// (simplify higher powers)
		else if (! strcmp(argv[1], "matprops"))
		{
			int maxpower;
			
			BigIntTP bigMod;
			BigIntTP temp  = NULL;
			BigIntTP temp2 = NULL;
			BigIntTP temp3 = NULL;
			
			BigIntMatrixTP A;
			BigPolyTP charaPoly;
			BigIntTP* OGcoeffs = NULL; //Holds the coefficients of charaPoly
			BigIntTP* newcoeffs = NULL; //Holds the coefficients for higher powers of A
			
			if (oargv[2][0] != '\0')
			{
				maxpower = (int)strtol(oargv[2], &tempStr, 10);
				
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read maxpower from command line.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
			}
			
			else
			{
				printf("Usage: " ANSI_COLOR_YELLOW "matprops " ANSI_COLOR_CYAN "maxpower [modulus]" ANSI_COLOR_RESET \
				": Calculates some useful matrix properties.\n");
				printf(" - " ANSI_COLOR_CYAN "maxpower" ANSI_COLOR_RESET \
				": The highest power of the update matrix to calculate an equivalent expression for.\n");
				printf(" - " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
				": Overrides the modulus specified in the config file.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (oargv[3][0] != '\0')
			{
				SET_BIG_NUM(oargv[3], bigMod, "Invalid modulus passed on command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Invalid modulus in config file.");
			}
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix from config file.\n");
				bigMod = free_BigIntT(bigMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			charaPoly = chara_poly(A, bigMod);
			printf("Matrix's characteristic polynomial modulo ");
			printi(bigMod);
			printf(": ");
			printp(charaPoly);
			printf("\n");
			
			//Check to make sure there's actually some computation to do
			if (degree(charaPoly) <= maxpower)
			{
				temp  = empty_BigIntT(1);
				temp2 = empty_BigIntT(1);
				temp3 = empty_BigIntT(1);
				OGcoeffs = extract_coefficients(charaPoly);
				
				//Now we need to divide by the leading coefficient, if it exists
				//If not, then don't do anything
				if (compare_BigIntT(OGcoeffs[degree(charaPoly)], one) != 0)
					big_num_inverse(OGcoeffs[degree(charaPoly)], bigMod, temp2);
				else
					copy_BigIntT(one, temp2);
				
				//Okay, we have the coefficients. Now to negate them.
				for (int i = 0; i < degree(charaPoly); i += 1)
				{
					subtract_BigIntT(bigMod, OGcoeffs[i], temp);
					mod_BigIntT(temp, bigMod, OGcoeffs[i]);
					multiply_BigIntT(temp2, OGcoeffs[i], temp);
					mod_BigIntT(temp, bigMod, OGcoeffs[i]);
				}
				
				//WHEN A NUMBER IS ZERO, ITS TERM SHOULDN'T BE PRINTED
				
				//This gives us the first equivalent expression
				printf("A^%d = ", degree(charaPoly));
				printi(OGcoeffs[0]);
				printf(" + ");
				for (int i = 1; i < degree(charaPoly)-1; i += 1)
				{
					printi(OGcoeffs[i]);
					printf("A^%d + ", i);
				}
				printi(OGcoeffs[degree(charaPoly)-1]);
				printf("A^%d\n", degree(charaPoly)-1);
				
				newcoeffs = malloc((degree(charaPoly)-1)*sizeof(BigIntTP));
				for (int i = 0; i < degree(charaPoly); i += 1)
				{
					newcoeffs[i] = empty_BigIntT(1);
					copy_BigIntT(OGcoeffs[i], newcoeffs[i]);
				}
				
				//Now, we compute the rest of the requested powers
				for (int i = degree(charaPoly)+1; i <= maxpower; i += 1)
				{
					//Shift all coefficients up a power
					copy_BigIntT(newcoeffs[degree(charaPoly)-1], temp);
					
					for (int n = degree(charaPoly)-2; n >= 0; n -= 1)
						copy_BigIntT(newcoeffs[n], newcoeffs[n+1]);
					copy_BigIntT(zero, newcoeffs[0]);
					
					//Now, simplify the big term into smaller exponent terms
					for (int n = 0; n < degree(charaPoly); n += 1)
					{
						multiply_BigIntT(temp, OGcoeffs[n], temp2);
						add_BigIntT(temp2, newcoeffs[n], temp3);
						mod_BigIntT(temp3, bigMod, newcoeffs[n]);
					}
					
					//Now, we print out the new expression
					printf("A^%d = ", i);
					printi(newcoeffs[0]);
					printf(" + ");
					for (int n = 1; n < degree(charaPoly)-1; n += 1)
					{
						printi(newcoeffs[n]);
						printf("A^%d + ", n);
					}
					printi(newcoeffs[degree(charaPoly)-1]);
					printf("A^%d\n", degree(charaPoly)-1);
				}
			}
			
			else
				printf("Degree of characteristic polynomial is less than maxpower. No computation to do.\n");
			
			for (int i = 0; i <= degree(charaPoly); i += 1)
			{
				if (i < degree(charaPoly))
					newcoeffs[i] = free_BigIntT(newcoeffs[i]);
				
				OGcoeffs[i] = free_BigIntT(OGcoeffs[i]);
			}
			FREE(OGcoeffs);
			FREE(newcoeffs);
			
			bigMod = free_BigIntT(bigMod);
			temp   = free_BigIntT(temp);
			temp2  = free_BigIntT(temp2);
			temp3  = free_BigIntT(temp3);
			
			A = free_BigIntMatrixT(A);
			
			charaPoly = free_BigPolyT(charaPoly);
		}
		
		
		//If we want to step over a matrix space to see what the
		// characteristic polynomials look like over specific
		// step sizes.
		else if (! strcmp(argv[1], "charawalk"))
		{
			PROHIBIT_UNIXFLAGS
			
			int counter;
			
			BigIntTP mod  = NULL;
			BigIntTP step = NULL; //How much to step in each direction in the matrix space
			BigIntTP temp;
			BigIntTP temp2;
			
			//Holds the BigIntTPs describing what matrix we're
			// currently testing
			BigIntTP** currentMatrixElements;
			
			BigPolyTP  charaPoly;
			BigPolyTP* factorList = NULL;
			
			BigIntMatrixTP startingMatrix;
			BigIntMatrixTP currentMatrix;
			BigIntMatrixTP identity;
			
			CycleInfoTP theCycle;
			
			startingMatrix = read_BigIntMatrixT(updatefilepath);
			if (startingMatrix == NULL)
			{
				fprintf(stderr, "Unable to read .matrix file at %s.\n", updatefilepath);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Getting user provided step size
			SET_BIG_NUM(argv[2], step, "Invalid step size passed on command line.");
			
			//If user provided a specific modulus on CLI
			if (argc > 3)
			{
				SET_BIG_NUM(argv[3], mod, "Invalid modulus passed on command line.");
			}
			
			else
			{
				SET_BIG_NUM(bigintmodstring, mod, "Unable to read modulus from config file.");
			}
			
			//Initialise currentMatrix
			currentMatrixElements = malloc(big_rows(startingMatrix)*sizeof(BigIntTP*));
			for (int i = 0; i < big_rows(startingMatrix); i += 1)
			{
				currentMatrixElements[i] = malloc(big_rows(startingMatrix)*sizeof(BigIntTP));
				for (int j = 0; j < big_rows(startingMatrix); j += 1)
				{
					currentMatrixElements[i][j] = empty_BigIntT(1);
					copy_BigIntT(big_element(startingMatrix, i, j), currentMatrixElements[i][j]);
				}
			}
			
			identity = identity_BigIntMatrixT(big_rows(startingMatrix));
			currentMatrix = new_BigIntMatrixT(big_rows(startingMatrix), big_rows(startingMatrix));
			set_big_matrix(currentMatrix, currentMatrixElements);
			
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			
			theCycle = new_CycleInfoT();
			
			//Check matrices at given step sizes apart in all directions
			for (int incRow = 0; incRow < big_rows(startingMatrix); incRow += 1)
			{
				for (int incCol = 0; incCol < big_rows(startingMatrix); incCol += 1)
				{
					//Step in one direction in the matrix space until we look back
					do
					{
						add_BigIntT(currentMatrixElements[incRow][incCol], step, temp);
						mod_BigIntT(temp, mod, currentMatrixElements[incRow][incCol]);
						set_big_matrix(currentMatrix, currentMatrixElements);
						
						big_floyd(currentMatrix, identity, mod, &theCycle);
						
						charaPoly  = chara_poly(currentMatrix, mod);
						factorList = old_factor_BigPolyT(charaPoly, mod);
						
						printbm(currentMatrix);
						printf("Cycle length: %d\n", omega(theCycle));
						printf("Transient length: %d\n", tau(theCycle));
						
						//It would be real beneficial to make a "clear_CycleInfoT() function"
						theCycle = free_CycleInfoT(theCycle);
						theCycle = new_CycleInfoT();
						
						//Print all found factors
						printp(charaPoly);
						printf("\n");
						copy_BigIntT(zero, temp);
						counter = 1;
						while (compare_BigIntT(temp, constant(factorList[0])) < 0)
						{
							printf("(");
							printp(factorList[counter]);
							printf(")");
							
							counter += 1;
							add_BigIntT(temp, one, temp2);
							copy_BigIntT(temp2, temp);
						}
						printf("\n");
					}
					while (!compare_BigIntMatrixT(startingMatrix, currentMatrix));
					
					printf("------------------------------\n");
				}
			}
			
			
			currentMatrixElements = free_BigIntT_array(currentMatrixElements, 
			                                           big_rows(startingMatrix), 
																								 big_rows(startingMatrix));
			
			mod   = free_BigIntT(mod);
			step  = free_BigIntT(step);
			temp  = free_BigIntT(temp);
			temp2 = free_BigIntT(temp2);
			
			charaPoly = free_BigPolyT(charaPoly);
			
			FREE(factorList);
			
			startingMatrix = free_BigIntMatrixT(startingMatrix);
			currentMatrix  = free_BigIntMatrixT(currentMatrix);
			identity       = free_BigIntMatrixT(identity);
		}
		
		
		//If we want to see all possible cycle lengths for a particular modulus
		else if (! strcmp(argv[1], "fibcyclelens"))
		{
			PROHIBIT_UNIXFLAGS
			
			BigIntTP bigModulus;
			
			int* cycleLengths = malloc(30*sizeof(int));
			int  indexCounter = 0;
			int  counter;
			
			//Fix hardcoding later
			BigIntTP startingX = empty_BigIntT(1);
			BigIntTP startingY = empty_BigIntT(1);
			
			BigIntTP currX = empty_BigIntT(1);
			BigIntTP currY = empty_BigIntT(1);
			
			BigIntTP tempInt = empty_BigIntT(1);
			
			int ok;
			
			//Use modulus provided on CLI
			if (argc == 3)
			{
				SET_BIG_NUM(argv[2], bigModulus, "Unable to read modulus provided on command line.");
			}
			
			//Use config modulus instead
			else
			{
				SET_BIG_NUM(bigintmodstring, bigModulus, "Unable to read modulus from config file.");
			}
			
			printf("Modulus: ");
			printi(bigModulus);
			printf("\n");
			
			//Iterate through all possible starting vectors
			while (compare_BigIntT(startingX, bigModulus) < 0)
			{
				while (compare_BigIntT(startingY, bigModulus) < 0)
				{
					//Prepping our starting vector
					copy_BigIntT(startingX, currX);
					copy_BigIntT(startingY, currY);
					counter = 0;
					
					//Increment until we get back to the start
					do
					{
						//Perform the iteration
						add_BigIntT(currX, currY, tempInt);
						mod_BigIntT(tempInt, bigModulus, currY);
						copy_BigIntT(currX, tempInt);
						copy_BigIntT(currY, currX);
						copy_BigIntT(tempInt, currY);
						
						counter += 1;
					}
					while ((compare_BigIntT(currX, startingX) != 0) ||
								 (compare_BigIntT(currY, startingY) != 0));
								 
					//Now, check to see if we've already gotten this cycle length
					ok = 1;
					for (int x = 0; x < indexCounter; x += 1)
					{
						if (cycleLengths[x] == counter)
						{
							ok = 0;
							break;
						}
					}
						
					//If we didn't find it, add the new cycle length
					if (ok)
					{
						cycleLengths[indexCounter] = counter;
						indexCounter += 1;
						printf("Cycle length: %d\n", counter);
						printf("Starting vector: [");
						printi(startingX);
						printf(", ");
						printi(startingY);
						printf("]\n");
					}
					
					//Incrementing Y
					add_BigIntT(startingY, one, tempInt);
					copy_BigIntT(tempInt, startingY);
				}
				//Incrementing X, resetting Y
				copy_BigIntT(zero, startingY);
				add_BigIntT(startingX, one, tempInt);
				copy_BigIntT(tempInt, startingX);
			}
			
			bigModulus = (bigModulus == NULL) ? NULL : free_BigIntT(bigModulus);
			startingX = free_BigIntT(startingX);
			startingY = free_BigIntT(startingY);
			currX     = free_BigIntT(currX);
			currY     = free_BigIntT(currY);
			
			tempInt  = free_BigIntT(tempInt);
			
			FREE(cycleLengths);
		}
		
		
		//If we want to check the Fibonacci numbers to see if multiples of numbers
		// appear before multiples of powers of those numbers
		else if (! strcmp(argv[1], "fibmultsearch"))
		{
			int hundred[] = {100};
			int oneArr[]  = {001};
			int start[]   = {001}; //Where the program starts counting
			
			//int debugCounter = 0; //For debugging
			BigIntTP upperbound;
			BigIntTP currNum = new_BigIntT(start, 1);
			BigIntTP counter = new_BigIntT(oneArr, 1);
			
			BigIntTP fibA    = new_BigIntT(oneArr, 1);
			BigIntTP fibB    = empty_BigIntT(1);
			BigIntTP fibTemp = empty_BigIntT(1);
			
			if (oargv[2][0] != '\0')
			{
				SET_BIG_NUM(oargv[2], upperbound, "Invalid upper bound passed.");
			}
			else
				upperbound = new_BigIntT(hundred, 1);
			
			printf("pos\tnum\t\tfib\n");
			while (compare_BigIntT(upperbound, currNum) >= 0)
			{	
				//Testing to see if our current Fib number is a multiple of currNum
				mod_BigIntT(fibA, currNum, fibTemp);
				
				if (compare_BigIntT(zero, fibTemp) == 0)
				{
					printi(counter);
					printf("\t");
					
					printi(currNum);
					printf("\tdivides\t");
					printi(fibA);
					
					//Now we check to see if it's a nultiple of some power of currNum
					divide_BigIntT(fibA, currNum, fibTemp);
					copy_BigIntT(fibTemp, fibA);
					mod_BigIntT(fibA, currNum, fibTemp);
					
					if (compare_BigIntT(zero, fibTemp) == 0)
					{
						printf(", along with ");
						printi(currNum);
						printf("^2");
					}
					
					printf("\n");
					
					//Prepare for the next number to check
					add_BigIntT(currNum, one, fibTemp);
					copy_BigIntT(fibTemp, currNum);
					copy_BigIntT(one, fibA);
					copy_BigIntT(zero, fibB);
					copy_BigIntT(one, counter);
				}
				
				else
				{
					//Calculate next Fibonacci number
					add_BigIntT(fibA, fibB, fibTemp);
					copy_BigIntT(fibA, fibB);
					copy_BigIntT(fibTemp, fibA);
					
					add_BigIntT(one, counter, fibTemp);
					copy_BigIntT(fibTemp, counter);
				}
			}
			
			upperbound = free_BigIntT(upperbound);
			currNum    = free_BigIntT(currNum);
			counter    = free_BigIntT(counter); 
			fibA       = free_BigIntT(fibA);
			fibB       = free_BigIntT(fibB);
			fibTemp    = free_BigIntT(fibTemp);
		}
		
		
		//If we want to iterate each vector in a space and record their cycle 
		// and transient lengths, then repeat for higher-powered moduli.
		else if (! strcmp(argv[1], "dynamics"))
		{
			int highpower = 2;
			int highmod;   //Holds the modulus raised to the correct power
			int maxmod;    //Holds the highest mod we obtain; used for iterating through matrix lifts

			int* possibleCycleLengths = NULL;
			int* currVectElements;
			
			//Each vector's lower-powered cycle lengths,
			// used for finding correct previous tuple
			int* currCycleLengths;
			
			int firstIndex; //Used to help create our lookup table for cycle lengths
			int tempModulus;
			
			//Holds counts for all possible cycle length pairs
			//For instance, if a vector % highmod has a cycle length
			// of 20, while % modulus it has a cycle length of 2,
			// then cycleTable[w=2][w=20] += 1;
			int** cycleTable = NULL;
			
			//Each time we go up a modulus, we have to keep track of the pevious
			// pairs of cycles we've found.
			int** prevCycleTuples = NULL;
			
			bool checkedAllVects = FALSE;
			bool matchedTuple;
			
			IntMatrixTP A; //Holds a lift to our update matrix
			IntMatrixTP originalA; //Holds our original update matrix
			IntMatrixTP currVect;
			
			int** matrixLiftOffsets; //Used to iterate over all possible lifts of the given matrix A
			int matrixLiftIncrement; //The step size when creating new matrix lifts
			
			CycleInfoTP theCycle = NULL;
			
			char* outputfilename = NULL;
			FILE* outputFile = NULL;
			
			//These variables are for recording all possible dynamic configurations
			// [[[1, 1, 1, 10], [1, 1, 4, 20], ...], [[1, 1, 1, 5], [1, 1, 2, 15], ...], ...]
			int*** allConfigs = NULL;      //Holds cleaner representations of all found dynamics configurations
			int* allConfigsLengths = NULL; //Holds how many cycle tuples each config in allConfigs has
			int numOfConfigs = 0;
			
			//Holds an example matrix for each different dynamic behaviour found
			IntMatrixTP* allConfigsMatrices = NULL;
			
			bool findAllConfigs = FALSE;
			bool fileoutput = FALSE;
			bool isUniqueConfig = TRUE; //Used to determine whether we've found a new configuration once we've computed a new one
			
			//User-provided upper bound for moduli power
			if (oargv[2][0] != '\0')
			{
				highpower = (int)strtol(oargv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read upper bound from command line.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
			}
			
			//If user provided a modulus for us
			if (oargv[3][0] != '\0')
			{
				modulus = (int)strtol(oargv[3], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read modulus from command line.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
			}
			
			//If user wants to find all possible dynamic configurations
			if (!strcmp(oargv[4], "TRUE"))
				findAllConfigs = TRUE;
				
			//If user wants to store an output file
			if (!strcmp(oargv[5], "TRUE"))
				fileoutput = TRUE;
			
			originalA = read_IntMatrixT(updatefilepath);
			if (originalA == NULL)
			{
				fprintf(stderr, "Unable to read .matrix file at %s.\n", updatefilepath);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			else if (rows(originalA) != cols(originalA))
			{
				printf("Update matrix provided is not square.\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			printf("Matrix:\n");
			printm(originalA);
			printf("~~~~~~~\n");
			A = new_IntMatrixT(rows(originalA), rows(originalA));
			copy_IntMatrixT(originalA, A);
			
			//Creating filename
			if ((!findAllConfigs) && (fileoutput))
			{
				outputfilename = malloc(MAXSTRLEN*sizeof(char));
				outputfilename[0] = '\0';
				
				ADD_FILE_PREFIX(outputfilename);
				
				strcat(outputfilename, "dynamics "); //dynamics
				strcat(outputfilename, oargv[2]); //maxPower
				strcat(outputfilename, " ");
				
				//Add modulus
				append_int(outputfilename, modulus);
				strcat(outputfilename, " ");
				
				//Now, include the update matrix elements
				strcat(outputfilename, "F");
				for (int x = 0; x < rows(A); x += 1)
					for (int y = 0; y < rows(A); y += 1)
						append_int(outputfilename, element(A, x, y));
					
				strcat(outputfilename, ".txt");
				printf("Vectors will be saved to %s\n", outputfilename);
				
				//Now, open the output file
				outputFile = fopen(outputfilename, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
			}
			
			//Initialise things for matrix lifts
			matrixLiftIncrement = modulus;
			
			matrixLiftOffsets = malloc(rows(A)*sizeof(int*));
			for (int i = 0; i < rows(A); i += 1)
				matrixLiftOffsets[i] = calloc(rows(A), sizeof(int));
			
			maxmod = modulus;
			for (int m = 1; m < highpower; m += 1)
				maxmod *= modulus;
			
			currVectElements = calloc(rows(A), sizeof(int));
			currVect         = new_IntMatrixT(rows(A), 1);
			
			//Dealing with exiting when highpower == 1 is done near the bottom of the tool
			while ((increment_int_array(matrixLiftOffsets, rows(A), rows(A), 1, maxmod/matrixLiftIncrement) != TRUE) ||
			       (highpower == 1))
			{
				highmod = 1;
				for (int modulusCounter = 1; modulusCounter <= highpower; modulusCounter += 1)
				{
					checkedAllVects = FALSE;
					
					highmod *= modulus;
					
					//Now that we have the matrix, calculate the cycle length
					// and use it to reason all possible cycle lengths
					// for alloting memory
					theCycle = floyd(A, A, highmod);
					
					if (!findAllConfigs)
					{
						printf("Current modulus: %d\n", highmod);
						printf("Matrix's multiplicative order mod %d: %d\n", highmod, omega(theCycle));
					}
					
					//First number in the pointer says how many cycle lengths
					// are stored in the vector
					if (possibleCycleLengths != NULL)
					{
						FREE(possibleCycleLengths);
					}
					possibleCycleLengths = malloc(2*sizeof(int));
					possibleCycleLengths[0] = 0;
					
					//Adding all divisors of the matrix's cycle length to our list
					//Reallocating the array each time is slow, but the numbers
					// we're dealing with should be low enough so it isn't a problem
					for (int i = 1; i <= omega(theCycle)/2; i += 1)
						if (omega(theCycle) % i == 0)
						{
							possibleCycleLengths[0] += 1;
							possibleCycleLengths[possibleCycleLengths[0]] = omega(theCycle) / i;
							possibleCycleLengths = realloc(possibleCycleLengths, 
																						 (possibleCycleLengths[0]+2)*sizeof(int));
						}
						
					//Add 1 as a cycle length
					possibleCycleLengths[0] += 1;
					possibleCycleLengths[possibleCycleLengths[0]] = 1;
					
					theCycle = free_CycleInfoT(theCycle);
					
					//Now, construct our cycleTable to hold cycle counts
					//This gets reallocated for each power we go up
					if (modulusCounter == 1)
					{
						//We only need one row on the first pass; we have no tuples yet
						cycleTable    = malloc(sizeof(int*));
						cycleTable[0] = calloc(possibleCycleLengths[0], sizeof(int));
					}
					
					else
					{
						cycleTable = malloc(prevCycleTuples[0][0]*sizeof(int*));
						for (int i = 0; i < prevCycleTuples[0][0]; i += 1)
						{
							//This line can be improved by only allocating enough
							// ints for those equal to or greater than the cycles
							// within the tuple
							cycleTable[i] = calloc(possibleCycleLengths[0], sizeof(int));
						}
					}
					
					//Helps find which previous tuple this vector corresponds to
					//Contains the cycle lengths of our vector under all moduli
					// less than the current one. The last element is the
					// current modulus' cycle length
					currCycleLengths = malloc(modulusCounter*sizeof(int));
				
					//Iterate until we've tested every vector
					while (!checkedAllVects)
					{
						set_column(currVect, currVectElements);
						
						//Start with highest modulus, work our way down
						//Allows taking vectors mod tempModulus to be easier
						tempModulus = modulus;
						for (int i = 1; i < modulusCounter; tempModulus *= modulus, i += 1);
						
						for (int i = modulusCounter-1; i >= 0; i -= 1)
						{
							theCycle = floyd(A, currVect, tempModulus);
							currCycleLengths[i] = omega(theCycle);
							
							theCycle = free_CycleInfoT(theCycle);
							
							tempModulus /= modulus;
							modm(currVect, tempModulus);
						}
						
						/*
						printf("currTuple: (");
						for (int cc = 0; cc < modulusCounter; cc += 1)
							printf("%d,", currCycleLengths[cc]);
						printf(")\n");*/

						//Now, output to our file the current vector's group only if:
						// we're on the last modulus to check and if the file exists
						if ((modulusCounter == highpower) && (outputFile != NULL))
						{
							fprintf(outputFile, "(");
							for (int i = 0; i < modulusCounter-1; i += 1)
								fprintf(outputFile, "%d, ", currCycleLengths[i]);
							fprintf(outputFile, "%d) : ", currCycleLengths[modulusCounter-1]);
							
							//Now, output the vector
							set_column(currVect, currVectElements);							
							fprintm_row(outputFile, currVect);
							fprintf(outputFile, "\n");
							
							//Now, save the file
							if (fclose(outputFile) == EOF)
							{
								fprintf(stderr, "Unable to properly save data. Continuing without saving...\n");
								outputFile = NULL;
							}
							
							else
								outputFile = fopen(outputfilename, "a");
						}
						
						//Now, find the correct firstIndex which matches our cycle lengths
						if (modulusCounter == 1)
							firstIndex = 0;
						
						else
						{
							for (int tuple = 1; tuple <= prevCycleTuples[0][0]; tuple += 1)
							{
								//Now, check the specific elements of each tuple to see if they
								// match with ours
								matchedTuple = TRUE;
								
								for (int tupleElem = 0; tupleElem < modulusCounter-1; tupleElem += 1)
								{
									if (prevCycleTuples[tuple][tupleElem] != currCycleLengths[tupleElem])
									{
										matchedTuple = FALSE;
										break;
									}
								}
								
								//The for loop will reach this at some point since the vectors must fall into
								// one of the defined cycle reductions (tuples)
								if (matchedTuple)
								{
									/*
									printf("Matched prev tuple: (");
									for (int pc = 0; pc < modulusCounter-1; pc += 1)
										printf("%d,", prevCycleTuples[tuple][pc]);
									printf(")\n");*/
									
									firstIndex = tuple-1;
									break;
								}
							}
						}
						
						//Now, we simply add our count to the cycleTable entry
						//Find our secondIndex
						for (int secondIndex = 1; 
						secondIndex <= possibleCycleLengths[0]; 
						secondIndex += 1)
						{
							if (possibleCycleLengths[secondIndex] == currCycleLengths[modulusCounter-1])
							{
								/*
								printf("Cycle length found: %d\n", possibleCycleLengths[secondIndex]);*/
								cycleTable[firstIndex][secondIndex-1] += 1;
								break;
							}
						}

						//Increment through all possible vectors (% modulus)
						//When I switch this program to use BigIntMatrixTPs, I can use my
						// new function here
						checkedAllVects = TRUE;
						for (int i = rows(A)-1; i >= 0; i -= 1)
						{
							if (currVectElements[i] != highmod - 1)
							{
								checkedAllVects = FALSE;
								currVectElements[i] += 1;
								break;
							}
							
							//Carry to next vector component
							else
								currVectElements[i] = 0;
						}
					}
					
					FREE(currCycleLengths);
					
					//Now that we have a completed cycleTable up to a particular modulus,
					// we need to get our list of tuples from it.
					if (prevCycleTuples == NULL)
					{
						prevCycleTuples    = malloc(sizeof(int*));
						prevCycleTuples[0] = malloc(sizeof(int));
						prevCycleTuples[0][0] = 0;
						
						//Add cycle lengths as tuples
						for (int i = 0; i < possibleCycleLengths[0]; i += 1)
						{
							if (cycleTable[0][i] != 0)
							{
								//Increasing tuple count
								prevCycleTuples[0][0] += 1;
								prevCycleTuples = realloc(prevCycleTuples, (prevCycleTuples[0][0]+1)*sizeof(int*));
								prevCycleTuples[prevCycleTuples[0][0]] = malloc(sizeof(int));
								prevCycleTuples[prevCycleTuples[0][0]][0] = possibleCycleLengths[i+1];
							}
						}
						
						//Prevents a segmentation fault if highpower == 1
						if (modulusCounter != highpower)
						{
							FREE(cycleTable[0]);
							FREE(cycleTable);
						}
						
						if (highpower == 1)
						{
							//If we're only checking one iteration, we don't need to
							// prepare this value for a next iteration. It does, however,
							// need to say how many tuples we had last iteration, which 
							// in this case was 1 (the null tuple)
							prevCycleTuples[0][0] = 1;
						}
					}
					
					else
					{
						//Get cycleTable ready for next iteration
						if (modulusCounter != highpower)
						{
							int prevPrevTupleCount = prevCycleTuples[0][0];
							
							//New tuples get put here, then prevCycleTuples points to this after
							int** newPrevCycleTuples = malloc(sizeof(int*));
							newPrevCycleTuples[0] = calloc(1, sizeof(int));
							
							for (int tuple = 0; tuple < prevPrevTupleCount; tuple += 1)
							{
								for (int newCycle = 0; newCycle < possibleCycleLengths[0]; newCycle += 1)
								{
									if (cycleTable[tuple][newCycle] != 0)
									{
										//Add new tuple
										newPrevCycleTuples[0][0] += 1;
										newPrevCycleTuples = realloc(newPrevCycleTuples, (newPrevCycleTuples[0][0]+1)*sizeof(int*));
										newPrevCycleTuples[newPrevCycleTuples[0][0]] = malloc(modulusCounter*sizeof(int));
										
										//Re-add elements from previous tuple, add new cycle length
										for (int tupleElem = 0; tupleElem < modulusCounter-1; tupleElem += 1)
											newPrevCycleTuples[newPrevCycleTuples[0][0]][tupleElem] = prevCycleTuples[tuple+1][tupleElem];
										newPrevCycleTuples[newPrevCycleTuples[0][0]][modulusCounter-1] = possibleCycleLengths[newCycle+1];
									}
								}
							}

							for (int i = 0; i < prevPrevTupleCount; i += 1)
							{
								FREE(cycleTable[i]);
							}
							FREE(cycleTable);
							
							//Now, deallocate prevCycleLengths, have it point to new tuples
							for (int i = 1; i < prevCycleTuples[0][0]; i += 1)
							{
								FREE(prevCycleTuples[i]);
							}
							FREE(prevCycleTuples[0]);
							FREE(prevCycleTuples);
							prevCycleTuples = newPrevCycleTuples;
						}
					}
				}
				
				//Now, we print all the numbers we've collected, or prepare to check for a unique tuple
				if (!findAllConfigs)
					printf("\nCycle tuples:\n");
				else
				{
					//Allocating space to hold possible new config and its length
					numOfConfigs += 1;
					allConfigsLengths = realloc(allConfigsLengths, numOfConfigs*sizeof(int));
					allConfigsLengths[numOfConfigs-1] = 0;
					
					allConfigs = realloc(allConfigs, numOfConfigs*sizeof(int**));
					allConfigs[numOfConfigs-1] = NULL;
				}
				
				//Loop through all our arrays and extract the tuple data
				for (int i = 0; i < prevCycleTuples[0][0]; i += 1)
				{
					for (int j = 0; j < possibleCycleLengths[0]; j += 1)
					{
						//Only print out info if vectors have that tuple of cycle lengths
						if (cycleTable[i][j] != 0)
						{
							if (findAllConfigs)
							{
								//Create proper space for storing current tuple
								allConfigsLengths[numOfConfigs-1] += 1;
								allConfigs[numOfConfigs-1] = realloc(allConfigs[numOfConfigs-1], 
																										 allConfigsLengths[numOfConfigs-1]*sizeof(int*));
								allConfigs[numOfConfigs-1][allConfigsLengths[numOfConfigs-1]-1] = malloc((highpower+1)*sizeof(int));
							}
							
							else
								printf("(");
							
							//Here, we either print the tuples or store them in allConfigs, checking whether they're unique
							//printf("(");
							for (int tupleElem = 0; tupleElem < highpower-1; tupleElem += 1)
							{
								if (!findAllConfigs)
									printf("%d, ", prevCycleTuples[i+1][tupleElem]);
								
								//Store tuple in our record of tuple configurations, even if it's not a new one
								//We'll worry about whether it's a new one later on
								// If it's not new, we'll just reallocate our memory back to how it was before
								else
									allConfigs[numOfConfigs-1][allConfigsLengths[numOfConfigs-1]-1][tupleElem] = prevCycleTuples[i+1][tupleElem];
								
								//printf("%d,", prevCycleTuples[i+1][tupleElem]);
							}

							//printf("%d) = %d\n", possibleCycleLengths[j+1], cycleTable[i][j]);
							
							if (!findAllConfigs)
								printf("%d) = %d\n", possibleCycleLengths[j+1], cycleTable[i][j]);
							
							else
							{
								allConfigs[numOfConfigs-1][allConfigsLengths[numOfConfigs-1]-1][highpower-1] = possibleCycleLengths[j+1];
								allConfigs[numOfConfigs-1][allConfigsLengths[numOfConfigs-1]-1][highpower] = cycleTable[i][j];
							}
						}
					}
				}
				
				//Now, we check to see if we've seen this specific tuple before
				//Essentially, we're trying to prove that it's not unique
				isUniqueConfig = TRUE;
				
				//Print out configs 
				/*
				for (int prevConfig = 0; prevConfig < numOfConfigs-1; prevConfig += 1)
				{
					for (int tuple = 0; tuple < allConfigsLengths[prevConfig]; tuple += 1)
					{
						printf("(");
						for (int entry = 0; entry < highpower+1; entry += 1)
							printf("%d,", allConfigs[prevConfig][tuple][entry]);
						printf(")\n");
					}
					printf("~~~\n");
				}
				*/
				
				//Determining whether configurations are equal or not hinges on the
				// assumptipon that configurations with equivalent tuples will appear
				// in exactly the same order.
				//This SHOULD always be the case due to how configs are copied into
				// allConfigs, but I'm leaving a note here just in case.
				for (int prevConfig = 0; prevConfig < numOfConfigs-1; prevConfig += 1)
				{					
					//If the two configs have the same number of tuples
					if (allConfigsLengths[prevConfig] == allConfigsLengths[numOfConfigs-1])
					{
						isUniqueConfig = FALSE;
						for (int currTuple = 0; currTuple < allConfigsLengths[numOfConfigs-1]; currTuple += 1)
						{
							for (int entry = 0; entry < highpower+1; entry += 1)
							{
								//If at least one number in the tuple is different
								if (allConfigs[prevConfig][currTuple][entry] != allConfigs[numOfConfigs-1][currTuple][entry])
								{
									isUniqueConfig = TRUE;
									break;
								}
							}
								
							if (isUniqueConfig)
								break;
						}
					}
					
					//If we've guaranteed that this configuration isn't unique
					if (!isUniqueConfig)
					{
						//printf("#%d is not unique!\n", numOfConfigs);
						//Deallocate new tuples since they're not actually unique and we don't need them
						//This is definitely inefficient, but I don't want a potential memory leak
						for (int n = 0; n < allConfigsLengths[numOfConfigs-1]; n += 1)
						{
							FREE(allConfigs[numOfConfigs-1][n]);
						}
						FREE(allConfigs[numOfConfigs-1]);
						allConfigsLengths[numOfConfigs-1] = 0;
						numOfConfigs -= 1;
						break;
					}
				}
				
				//If it is unique, add the matrix to the array 
				//The allConfigsMatrices == NULL prevents already stored matrices from being overriden when there's only one unique config
				if ((findAllConfigs) && ((isUniqueConfig) || ((numOfConfigs == 1) && (allConfigsMatrices == NULL))))
				{
					allConfigsMatrices = realloc(allConfigsMatrices, numOfConfigs*sizeof(IntMatrixTP));
					allConfigsMatrices[numOfConfigs-1] = new_IntMatrixT(rows(A), rows(A));
					copy_IntMatrixT(A, allConfigsMatrices[numOfConfigs-1]);
				}
				
				if ((outputFile != NULL) && (fclose(outputFile) == EOF))
					fprintf(stderr, "Unable to save data.\n");
				
				//Don't need to free theCycle since it gets freed above
				theCycle = NULL;
				
				for (int i = 0; i < prevCycleTuples[0][0]; i += 1)
				{
					FREE(cycleTable[i]);
				}
				FREE(cycleTable);
				FREE(possibleCycleLengths);
				
				for (int i = 1; i < prevCycleTuples[0][0]; i += 1)
				{
					FREE(prevCycleTuples[i]);
				}
				FREE(prevCycleTuples[0]);
				FREE(prevCycleTuples);
				
				if (outputfilename != NULL)
				{
					FREE(outputfilename);
				}
				
				//Ensuring we don't have to loop a bunch of times when only computing one dynamics configuration
				if ((!findAllConfigs) || (highpower == 1))
					break;
				
				//Set A to the next lift
				else
				{
					for (int r = 0; r < rows(originalA); r += 1)
						for (int c = 0; c < rows(originalA); c += 1)
						{
							matrixLiftOffsets[r][c] *= matrixLiftIncrement;
							matrixLiftOffsets[r][c] += element(originalA, r, c);
						}
						
					set_matrix(A, matrixLiftOffsets);
					modm(A, maxmod);
					
					for (int r = 0; r < rows(originalA); r += 1)
						for (int c = 0; c < rows(originalA); c += 1)
						{
							matrixLiftOffsets[r][c] -= element(originalA, r, c);
							matrixLiftOffsets[r][c] /= matrixLiftIncrement;
						}
				}
			}
			
			//Now, let's try and print out our record of unique dynamics configurations
			for (int config = 0; config < numOfConfigs; config += 1)
			{
				printf("Config #%d:\n", config);
				for (int tuple = 0; tuple < allConfigsLengths[config]; tuple += 1)
				{
					printf("(");
					for (int num = 0; num < highpower-1; num += 1)
						printf("%d, ", allConfigs[config][tuple][num]);
					printf("%d) -> %d\n", allConfigs[config][tuple][highpower-1], allConfigs[config][tuple][highpower]);
				}
				printf("Representative matrix:\n");
				printm(allConfigsMatrices[config]);
				printf("\n");
			}
			
			if (findAllConfigs)
			{
				for (int config = 0; config < numOfConfigs; config += 1)
				{
					allConfigsMatrices[config] = free_IntMatrixT(allConfigsMatrices[config]);
					for (int tuple = 0; tuple < allConfigsLengths[config]; tuple += 1)
					{
						FREE(allConfigs[config][tuple]);
					}
					FREE(allConfigs[config]);
				}
				FREE(allConfigs);
				
				FREE(allConfigsLengths);
				
				FREE(allConfigsMatrices);
			}
			
			for (int i = 0; i < rows(originalA); i += 1)
			{
				FREE(matrixLiftOffsets[i]);
			}
			FREE(matrixLiftOffsets);
			
			FREE(currVectElements);
			currVect = free_IntMatrixT(currVect);
				
			originalA = free_IntMatrixT(originalA);
			A = free_IntMatrixT(A);
		}
		
		
		//If the user wants to see how orbits from higher moduli map to lower moduli orbits
		else if (! strcmp(argv[1], "orbitmaps"))
		{
			PROHIBIT_UNIXFLAGS
			
			BigIntTP bigMod;
			BigIntTP bigModPower; //Holds bigMod^maxpower
			BigIntTP temp;
			BigIntTP temp2;
			
			int oneArr[1] = {1};
			
			BigIntMatrixTP A;
			BigIntMatrixTP currVect;
			BigIntMatrixTP tempVect;
			BigIntMatrixTP tempVect2;
			BigIntTP** currVectElements;
			
			CycleInfoTP theCycle = NULL;
			
			BigIntMatrixTP** orbitReps;
			int** orbitLengths;
			int** orbitMaps; //Keeps track of how the orbits map onto each other
			int*  numOfOrbits;
			
			int   groupCount = 0;
			int   currCycleLength;
			int   prevCycleLength; //Used to sort mappings by cycle length
			int*  currOrbitMapMatch = NULL; //For organising output into groupings of reductions
			int** currOrbitMapGroup = NULL; //Holds the vectors that match our match
			
			int maxpower;     //Holds the max power as a regular int
			int reducedOrbit; //Holds which orbit our current orbit feeds into
			
			bool isNewOrbit; //For checking whether particular orbits are already in our list
			bool matchingMap; //For finding maps that match our grouping for output
			bool printedAllMaps = FALSE;
			bool newCycleLength; //For determining whether we found a new cycle length when printing results
			
			FILE* outputFile = NULL;
			char* outputFileName = NULL;
			bool fileoutput = FALSE;
			
			FILE* graphFile = NULL; //GRAPHVIS output file
			char* graphFileName = NULL;
			int tempCounter = 0;
			
			if (argc > 2)
			{
				maxpower = (int)strtol(argv[2], &tempStr, 10);
				if (tempStr[0] != '\0')
				{
					fprintf(stderr, "Unable to read maxpower from command line.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
			}
			else
			{
				printf(ANSI_COLOR_YELLOW "orbitmaps " ANSI_COLOR_CYAN "maxpower [modulus] [fileoutput]" ANSI_COLOR_RESET \
				": Calculates orbit representatives for a higher-powered modulus and sees how they map down to lower-powered moduli.\n");
				printf(" - " ANSI_COLOR_CYAN "maxpower" ANSI_COLOR_RESET \
				": The highest (initial) power to use for the modulus.\n");
				printf(" - " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET \
				": Says what (base) modulus to use.\n\n");
				printf(" - " ANSI_COLOR_CYAN "[fileoutput]" ANSI_COLOR_RESET \
				": Says whether to create an output file.\n\n");
				
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			if (argc > 3)
			{
				SET_BIG_NUM(argv[3], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			if (argc > 4)
				if (!strcmp(argv[4], "TRUE"))
					fileoutput = TRUE;
			
			A = read_BigIntMatrixT(updatefilepath);
			if (A == NULL)
			{
				fprintf(stderr, "Unable to read update matrix provided in config file.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Create filename
			if (fileoutput)
			{
				outputFileName = malloc(MAXSTRLEN*sizeof(char));
				outputFileName[0] = '\0';
				graphFileName = malloc(MAXSTRLEN*sizeof(char));
				graphFileName[0] = '\0';
				
				ADD_FILE_PREFIX(outputFileName);
				ADD_FILE_PREFIX(graphFileName);
				
				strcat(outputFileName, "orbitmaps ");
				strcat(graphFileName, "graph ");
				append_int(outputFileName, maxpower);
				append_int(graphFileName, maxpower);
				strcat(outputFileName, " ");
				strcat(graphFileName, " ");
				append_BigIntT(outputFileName, bigMod);
				append_BigIntT(graphFileName, bigMod);
				strcat(outputFileName, " F");
				strcat(graphFileName, " F");
				for (int x = 0; x < big_rows(A); x += 1)
					for (int y = 0; y < big_cols(A); y += 1)
					{
						append_BigIntT(outputFileName, big_element(A, x, y));
						append_BigIntT(graphFileName, big_element(A, x, y));
					}
					
				strcat(outputFileName, ".txt");
				strcat(graphFileName, ".graph");
				
				outputFile = fopen(outputFileName, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create output file. Continuing without saving...\n");
			}
			
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			
			bigModPower = new_BigIntT(oneArr, 1);
			
			orbitReps    = calloc(maxpower, sizeof(BigIntMatrixTP*));
			orbitLengths = calloc(maxpower, sizeof(int*));
			numOfOrbits  = calloc(maxpower, sizeof(int));
			
			currVectElements = new_BigIntT_array(big_rows(A), 1);
			currVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect  = new_BigIntMatrixT(big_rows(A), 1);
			tempVect2 = new_BigIntMatrixT(big_rows(A), 1);
			
			//Calculate bigMod^maxpower
			for (int i = 0; i < maxpower; i += 1)
			{
				multiply_BigIntT(bigMod, bigModPower, temp);
				copy_BigIntT(temp, bigModPower);
			}
			
			printf("Base modulus: ");
			printi(bigMod);
			printf("\nMax power: %d\n", maxpower);
			printf("Matrix:\n");
			printbm(A);
			printf("~~~~~~~~~\n");
			
			//Iterate over each vector in our highest modulus to find all orbits
			do
			{
				set_big_matrix(currVect, currVectElements);
				big_floyd(A, currVect, bigModPower, &theCycle);
				
				if (tau(theCycle) != 0)
					continue;
				
				//If we don't have any orbits yet
				if (numOfOrbits[maxpower-1] == 0)
				{
					numOfOrbits[maxpower-1] += 1;
					orbitReps[maxpower-1] = realloc(orbitReps[maxpower-1], numOfOrbits[maxpower-1]*sizeof(BigIntMatrixTP));
					orbitLengths[maxpower-1] = realloc(orbitLengths[maxpower-1], numOfOrbits[maxpower-1]*sizeof(int));
					
					//Saving orbit data
					orbitReps[maxpower-1][numOfOrbits[maxpower-1]-1] = new_BigIntMatrixT(big_rows(A), 1);
					copy_BigIntMatrixT(rep(theCycle), orbitReps[maxpower-1][numOfOrbits[maxpower-1]-1]); //Does it implicitly cast?
					orbitLengths[maxpower-1][numOfOrbits[maxpower-1]-1] = omega(theCycle);
				}
				
				//If we already have orbits, we have to check and see if 
				// we've already accounted for our current orbit
				else
				{
					isNewOrbit = TRUE;
					
					//Iterate over all the orbits we already have
					for (int i = 0; i < numOfOrbits[maxpower-1]; i += 1)
					{
						copy_BigIntMatrixT(orbitReps[maxpower-1][i], tempVect);
						
						//Check to see if any vector in the cycle matches currVect
						do
						{
							if (compare_BigIntMatrixT(tempVect, currVect))
							{
								isNewOrbit = FALSE;
								break;
							}
							
							big_mat_mul(A, tempVect, tempVect2);
							modbm(tempVect2, bigModPower);
							copy_BigIntMatrixT(tempVect2, tempVect);
						}
						while (!compare_BigIntMatrixT(tempVect, orbitReps[maxpower-1][i])); //Will this implicitly cast?
						
						if (!isNewOrbit)
							break;
					}
					
					if (isNewOrbit)
					{
						numOfOrbits[maxpower-1] += 1;
						orbitReps[maxpower-1] = realloc(orbitReps[maxpower-1], numOfOrbits[maxpower-1]*sizeof(BigIntMatrixTP));
						orbitLengths[maxpower-1] = realloc(orbitLengths[maxpower-1], numOfOrbits[maxpower-1]*sizeof(int));
						
						orbitLengths[maxpower-1][numOfOrbits[maxpower-1]-1] = omega(theCycle);
						orbitReps[maxpower-1][numOfOrbits[maxpower-1]-1] = new_BigIntMatrixT(big_rows(A), 1);
						copy_BigIntMatrixT(rep(theCycle), orbitReps[maxpower-1][numOfOrbits[maxpower-1]-1]);
					}
				}
			}
			while (!increment_BigIntT_array(currVectElements, big_rows(A), 1, one, bigModPower));
			
			orbitMaps = malloc(numOfOrbits[maxpower-1]*sizeof(int*));
			for (int i = 0; i < numOfOrbits[maxpower-1]; i += 1)
			{
				orbitMaps[i] = calloc(maxpower, sizeof(int));
				orbitMaps[i][maxpower-1] = i;
			}
			
			//Now, we have to find the orbitreps for all the lower moduli
			for (int i = maxpower-2; i >= 0; i -= 1)
			{
				//Calculate new modulus
				divide_BigIntT(bigModPower, bigMod, temp);
				copy_BigIntT(temp, bigModPower);
				
				//Go through each rep from the last modulus we calculated
				for (int cyc = 0; cyc < numOfOrbits[i+1]; cyc += 1)
				{
					copy_BigIntMatrixT(orbitReps[i+1][cyc], currVect);
					modbm(currVect, bigModPower);
					
					if (numOfOrbits[i] == 0)
					{
						reducedOrbit = 0;
						numOfOrbits[i] += 1;
						orbitReps[i] = realloc(orbitReps[i], numOfOrbits[i]*sizeof(BigIntMatrixTP));
						orbitLengths[i] = realloc(orbitLengths[i], numOfOrbits[i]*sizeof(int));
						
						big_floyd(A, currVect, bigModPower, &theCycle);
						orbitLengths[i][numOfOrbits[i]-1] = omega(theCycle);
						orbitReps[i][numOfOrbits[i]-1] = new_BigIntMatrixT(big_rows(A), 1);
						copy_BigIntMatrixT(rep(theCycle), orbitReps[i][numOfOrbits[i]-1]);
					}
					
					//Check to see if current higher mod orbit maps to a 
					// lower mod orbit we already have
					else
					{
						isNewOrbit = TRUE;
						
						//Iterate over all our current reps
						for (int newOrbs = 0; newOrbs < numOfOrbits[i]; newOrbs += 1)
						{
							copy_BigIntMatrixT(orbitReps[i][newOrbs], tempVect);
							do
							{
								//If our rep is already accounted for with a cycle
								if (compare_BigIntMatrixT(tempVect, currVect))
								{
									isNewOrbit = FALSE;
									reducedOrbit = newOrbs;
									break;
								}
								
								big_mat_mul(A, tempVect, tempVect2);
								modbm(tempVect2, bigModPower);
								copy_BigIntMatrixT(tempVect2, tempVect);
							}
							while (!compare_BigIntMatrixT(tempVect, orbitReps[i][newOrbs]));
							
							if (!isNewOrbit)
								break;
						}
						
						//Add new orbit to list if necessary
						if (isNewOrbit)
						{
							reducedOrbit = numOfOrbits[i];
							numOfOrbits[i] += 1;
							orbitReps[i] = realloc(orbitReps[i], numOfOrbits[i]*sizeof(BigIntMatrixTP));
							orbitLengths[i] = realloc(orbitLengths[i], numOfOrbits[i]*sizeof(int));
							
							big_floyd(A, currVect, bigModPower, &theCycle);
							orbitLengths[i][numOfOrbits[i]-1] = omega(theCycle);
							orbitReps[i][numOfOrbits[i]-1] = new_BigIntMatrixT(big_rows(A), 1);
							copy_BigIntMatrixT(rep(theCycle), orbitReps[i][numOfOrbits[i]-1]);
						}
					}
					
					//Make sure to keep track of mappings
					for (int map = 0; map < numOfOrbits[maxpower-1]; map += 1)
						if (orbitMaps[map][i+1] == cyc)
							orbitMaps[map][i] = reducedOrbit;
				}
			}
			
			//Now, let's check to see if I got the orbitreps correctly
			for (int e = maxpower-1; e >= 0; e -= 1)
			{
				printf("For mod ");
				printi(bigMod);
				printf("^%d:\n", e+1);
				
				for (int i = 0; i < numOfOrbits[e]; i += 1)
				{
					printf("Rep: ");
					if (vectorType == row)
					{
						printbm_row(orbitReps[e][i]);
						printf(",\tOrbit length: %d\n", orbitLengths[e][i]);
					}
					else if (vectorType == col)
					{
						printf("\n");
						printbm(orbitReps[e][i]);
						printf("Orbit length: %d\n", orbitLengths[e][i]);
					}
					
				}
				
				printf("\n");
			}
			
			currOrbitMapMatch = calloc((maxpower-1), sizeof(int));
			
			//Go through the mappings, organise them based on common mappings
			while (!printedAllMaps)
			{
				groupCount = 0;
				
				for (int map = 0; map < numOfOrbits[maxpower-1]; map += 1)
				{
					matchingMap = TRUE;
					
					//Checking each individual reduction to see if they match
					for (int modLevel = 0; modLevel < maxpower-1; modLevel += 1)
					{
						//If the map we checked isn't equal somewhere
						if (orbitMaps[map][modLevel] != currOrbitMapMatch[modLevel])
						{
							matchingMap = FALSE;
							break;
						}
					}
					
					//Add matching map to the group
					if (matchingMap)
					{
						groupCount += 1;
						currOrbitMapGroup = realloc(currOrbitMapGroup, groupCount*sizeof(int*));
						currOrbitMapGroup[groupCount-1] = orbitMaps[map];
					}
				}
				
				//Now, we should have every matching mapping, sorted by cycle length
				//Print them out, and save any relevant files
				prevCycleLength = 0;
				while (TRUE)
				{
					newCycleLength = FALSE;
					currCycleLength = -1;
					
					//Go through mappings and find the next lowest cycle length
					for (int map = 0; map < groupCount; map += 1)
					{
						if ((orbitLengths[maxpower-1][currOrbitMapGroup[map][maxpower-1]] > prevCycleLength) &&
						    ((currCycleLength == -1) ||
								(orbitLengths[maxpower-1][currOrbitMapGroup[map][maxpower-1]] < currCycleLength)))
						{
							currCycleLength = orbitLengths[maxpower-1][currOrbitMapGroup[map][maxpower-1]];
							newCycleLength = TRUE;
						}
					}
					
					if (newCycleLength)
					{
						for (int map = 0; map < groupCount; map += 1)
						{
							if (orbitLengths[maxpower-1][currOrbitMapGroup[map][maxpower-1]] == currCycleLength)
							{
								for (int modLevel = maxpower-1; modLevel > 0; modLevel -= 1)
								{
									if (vectorType == row)
									{
										printbm_row(orbitReps[modLevel][currOrbitMapGroup[map][modLevel]]);
										printf(" (ω = %d) -> ", orbitLengths[modLevel][currOrbitMapGroup[map][modLevel]]);
									}
									else if (vectorType == col)
									{
										printbm(orbitReps[modLevel][currOrbitMapGroup[map][modLevel]]);
										printf("↓ (ω = %d)\n", orbitLengths[modLevel][currOrbitMapGroup[map][modLevel]]);
									}
									if (outputFile != NULL)
									{
										fprintbm_row(outputFile, orbitReps[modLevel][currOrbitMapGroup[map][modLevel]]);
										fprintf(outputFile, " (ω = %d) -> ", orbitLengths[modLevel][currOrbitMapGroup[map][modLevel]]);
									}
								}
								if (vectorType == row)
								{
									printbm_row(orbitReps[0][currOrbitMapGroup[map][0]]);
									printf(" (ω = %d)\n", orbitLengths[0][currOrbitMapGroup[map][0]]);
								}
								else if (vectorType == col)
								{
									printbm(orbitReps[0][currOrbitMapGroup[map][0]]);
									printf("(ω = %d)\n\n", orbitLengths[0][currOrbitMapGroup[map][0]]);
								}
								if (outputFile != NULL)
								{
									fprintbm_row(outputFile, orbitReps[0][currOrbitMapGroup[map][0]]);
									fprintf(outputFile, " (ω = %d)\n", orbitLengths[0][currOrbitMapGroup[map][0]]);
								}
							}
						}
						
						prevCycleLength = currCycleLength;
					}
					
					else
						break;
				}
				
				//Preventing extra newlines from cluttering the screen
				if (groupCount > 0)
				{
					printf("\n");
					if (outputFile != NULL)
						fprintf(outputFile, "\n");
				}
				
				//Now, we change which mapping we're looking for and repeat
				for (int i = maxpower-2; i >= 0; i -= 1)
				{
					printedAllMaps = TRUE;
					
					currOrbitMapMatch[i] += 1;
					if (i == 0)
					{
						printf("==========\n\n");
						if (outputFile != NULL)
							fprintf(outputFile, "=========\n\n");
					}
					
					if (currOrbitMapMatch[i] == numOfOrbits[i])
						currOrbitMapMatch[i] = 0;
					else
					{
						printedAllMaps = FALSE;
						break;
					}
				}
			}
			
			if (fileoutput)
			{
				graphFile = fopen(graphFileName, "w");
				if (graphFile == NULL)
					fprintf(stderr, "Unable to save .graph file.\n");
			}
			
			if (graphFile != NULL)
			{
				fprintf(graphFile, "~b:False\n"); //Turning on the bounding box
				
				fprintf(graphFile, "~n:");
				for (int i = 0; i < maxpower; i += 1)
					tempCounter += numOfOrbits[i];
				
				fprintf(graphFile, "%d\n~l:above\n", tempCounter);
				tempCounter = 0;
				
				//Give labels to every vector we've saved
				for (int p = 0; p < maxpower; p += 1)
					for (int v = 0; v < numOfOrbits[p]; v += 1)
					{
						fprintf(graphFile, "%d:", tempCounter);
						fprintbm_row(graphFile, orbitReps[p][v]);
						fprintf(graphFile, "\n");
						tempCounter += 1;
					}
					
				fprintf(graphFile, "~c:direction\n");
				
				//Add our reductions to the .graph file
				for (int map = 0; map < numOfOrbits[maxpower-1]; map += 1)
				{
					for (int pow = maxpower-1; pow > 0; pow -= 1)
					{
						//ID of first vector
						tempCounter = 0;
						for (int prevs = 0; prevs < pow; tempCounter += numOfOrbits[prevs], prevs += 1);
						tempCounter += orbitMaps[map][pow];
						fprintf(graphFile, "%d,", tempCounter);
						
						//ID of second vector
						tempCounter = 0;
						for (int prevs = 0; prevs < pow-1;  tempCounter += numOfOrbits[prevs], prevs += 1);
						tempCounter += orbitMaps[map][pow-1];
						fprintf(graphFile, "%d\n", tempCounter);
					}
				}
			}
			
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to close output file. Unable to save data.\n");
				
			if (graphFile != NULL)
				if (fclose(graphFile) == EOF)
					fprintf(stderr, "Unable to close .graph file. Unable to save graph.\n");
			
			bigMod      = free_BigIntT(bigMod);
			bigModPower = free_BigIntT(bigModPower);
			temp        = free_BigIntT(temp);
			temp2       = free_BigIntT(temp2);
			
			for (int i = 0; i < maxpower; i += 1)
			{
				for (int j = 0; j < numOfOrbits[i]; j += 1)
					orbitReps[i][j] = free_BigIntMatrixT(orbitReps[i][j]);
				
				FREE(orbitReps[i]);
				FREE(orbitLengths[i]);
			}
			FREE(orbitReps);
			FREE(orbitLengths);
			
			for (int i = 0; i < numOfOrbits[maxpower-1]; i += 1)
			{
				FREE(orbitMaps[i]);
			}
			FREE(orbitMaps);
			
			FREE(numOfOrbits);
			FREE(currOrbitMapMatch);
			FREE(currOrbitMapGroup);
			
			currVectElements = free_BigIntT_array(currVectElements, big_rows(A), 1);
			currVect  = free_BigIntMatrixT(currVect);
			tempVect  = free_BigIntMatrixT(tempVect);
			tempVect2 = free_BigIntMatrixT(tempVect2);
			theCycle  = free_CycleInfoT(theCycle);
			
			FREE(outputFileName);
			FREE(graphFileName);
			
			A = free_BigIntMatrixT(A);
		}
		
		
		//If the user wants to calculate the orbit maps for ALL lifts of a particular matrix
		else if (! strcmp(argv[1], "orbitmaps2"))
		{
			bool isNewVector;
			bool isAInvertible;
			bool useTreeSpaces = FALSE; //Says whether we should use spaces or numbers for tree indentation
			
			BigIntTP baseMod;  //The base prime
			int maxPower;
			int currPower = 1; //For keeping track of which layer of moduli we're on
			
			int numArr[1] = {1};
			const int numOfTemps = 3;
			BigIntTP temps[numOfTemps];
			
			//For keeping track of how far we are in the computations
			BigIntTP progressDenominator, progressNumerator;
			
			const BigIntMatrixTP initialA = UPDATEMATRIX;
			BigIntMatrixTP Ainv = NULL;
			
			BigIntMatrixTP** orbitReps;   //Holds representatives for all the orbits, from smallest modulus to biggest
			int*  orbitRepsCount;         //How many vectors are in each layer of orbitReps?
			int** orbitRepsCycleLengths;  //What are the cycle lengths of the orbit reps? 
			
			BigIntTP**  vectorLiftElements; //Holds the raw elements responsible for "lifting" each vector
			BigIntTP*** matrixLiftElements; //Holds the raw elements for each matrix lift at each level
			
			BigIntMatrixTP currMat;
			BigIntMatrixTP currPowMat; //For holding powers of currMat to speed up computations
			
			const int numOfTempVects = 3;
			BigIntMatrixTP tempVects[numOfTempVects];
			CycleInfoTP theCycle = NULL;

			int stemStart; //For simplifying tree-printing code
			
			FILE* outputFile     = NULL;
			char* outputFileName = NULL;
			char* outputFileBaseName = NULL;
			int   fileCounter = 0; //For numbering the different files created by the tool
			
			//File for holding the unique mapping structures found while computing
			FILE* mapOutputFile     = NULL;
			char* mapOutputFileName = NULL;
			int   currentLayer; //For keeping track of which layer we're printing from orbitMapsCatalogue
			
			//Matrices for starting and stopping computation, if requested
			BigIntMatrixTP resumeMatrix   = NULL;
			BigIntMatrixTP sentinelMatrix = NULL;
			
			//For keeping track of the structures of the orbit maps
			typedef struct orbitmapnode
			{
				int subNodes; //Says how many subnodes this particular node has
				int cyclen;   //The cycle length of this particular node
			}
			OMN;
			
			typedef struct orbitmaptree
			{
				int occurrences;  //How many times has this particular tree structure occurred?
				BigIntMatrixTP matrixRep; //A matrix that causes this particular tree structure
				
				int* layerCount; //How many nodes are on each layer?
				OMN** nodes;     //Array of nodes
			}
			OMT;
			
			bool isNewTree; //For saying whether we've found a new tree to add to our collection
			int orbitMappingID;             //For printing out which mapping configuration a matrix's LCA belongs to
			int orbitMapsCount      = 0;
			OMT* orbitMapsCatalogue = NULL; //Holds all the different structures we encounter when computing orbit maps
			
			//For determining if different tree structures are actually equal
			bool foundMatchingNode;
			OMN* nodeTracker = NULL;
			
			if (oargv[2][0] == '\0') //If the user hasn't given the mandatory argument
			{
				fprintf(stderr, "The argument 'maxpower' must be specified.\n");
				printf("Usage: " ANSI_COLOR_YELLOW "orbitmaps2 " ANSI_COLOR_CYAN "maxpower [modulus] [fileoutput] [belowBound] [aboveBound]" ANSI_COLOR_RESET "\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			if (big_rows(initialA) != big_cols(initialA))
			{
				fprintf(stderr, "Given update matrix isn't square.\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			maxPower = (int)strtol(oargv[2], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read maxpower from command line.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Has to go here since we don't read maxPower until now
			BigIntTP modList[maxPower]; //Holds every power of baseMod we'll need, from baseMod to baseMod^maxPower
			
			int stemCollectionCycleLengths[maxPower];
			int stemCollectionIndices[maxPower];
			BigIntMatrixTP stemCollection[maxPower]; //For holding stems in the correct order
			
			//For outputting unique tree structures to an output file
			int subCount[maxPower-1]; //Counts how many children nodes each node has in the respective layer
			int mapIndexCount[maxPower]; //Where we are in the currently-being-printed tree
			
			if (oargv[3][0] != '\0')
			{
				SET_BIG_NUM(oargv[3], baseMod, "Unable to read modulus from command line.");
				printf("modulus: %s\n", oargv[3]);
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, baseMod, "Unable to read modulus from config file.");
			}
			
			if (!strcmp(oargv[4], "TRUE"))
			{
				outputFileName     = malloc(MAXSTRLEN*sizeof(char));
				mapOutputFileName  = malloc(MAXSTRLEN*sizeof(char));
				outputFileBaseName = malloc(MAXSTRLEN*sizeof(char));
				outputFileBaseName[0] = '\0';
				mapOutputFileName[0]  = '\0';
				ADD_FILE_PREFIX(outputFileBaseName);
				ADD_FILE_PREFIX(mapOutputFileName);
				
				strcat(outputFileBaseName, "orbitmaps2 ");
				strcat(mapOutputFileName,  "orbitmaps2 maps ");
				strcat(outputFileBaseName, oargv[2]);
				strcat(mapOutputFileName,  oargv[2]);
				strcat(outputFileBaseName, " ");
				strcat(mapOutputFileName,  " ");
				append_BigIntT(outputFileBaseName, baseMod);
				append_BigIntT(mapOutputFileName,  baseMod);
				strcat(outputFileBaseName, " F");
				strcat(mapOutputFileName,  " F");
				for (int row = 0; row < big_rows(initialA); row += 1)
					for (int col = 0; col < big_rows(initialA); col += 1)
					{
						append_BigIntT(outputFileBaseName, big_element(initialA, row, col));
						append_BigIntT(mapOutputFileName, big_element(initialA, row, col));
						
						if ((row != big_rows(initialA)-1) || (col != big_rows(initialA)-1))
						{
							strcat(outputFileBaseName, " ");
							strcat(mapOutputFileName,  " ");
						}
						else
							strcat(outputFileBaseName, " f");
					}
				strcpy(outputFileName, outputFileBaseName);
				append_int(outputFileName, fileCounter);
				strcat(outputFileName, ".txt");
				strcat(mapOutputFileName, ".txt");
				
				//Don't store any matrix data in quiet mode
				if (!quietmode)
				{
					outputFile = fopen(outputFileName, "a");
					if (outputFile == NULL)
						fprintf(stderr, "Unable to open output file. Continuing without saving...\n");
				}
			}
			
			//If the user wants to resume computation at a particular matrix
			if (!strcmp(oargv[5], "TRUE"))
			{
				resumeMatrix = read_BigIntMatrixT(resumefilepath);
				if (resumeMatrix == NULL)
				{
					fprintf(stderr, "Unable to read resume matrix from config file.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
				
				else if ((big_rows(resumeMatrix) != big_rows(initialA)) ||
				         (big_cols(resumeMatrix) != big_cols(initialA)))
				{
					fprintf(stderr, "Given resume matrix doesn't match dimensions of update matrix.\n");
					returnvalue = EXIT_SUCCESS;
					goto FREE_VARIABLES;
				}
			}
			
			//If the user wants to stop computation at a particular matrix
			if (!strcmp(oargv[6], "TRUE"))
			{
				sentinelMatrix = read_BigIntMatrixT(sentinelfilepath);
				if (sentinelMatrix == NULL)
				{
					fprintf(stderr, "Unable to read sentinel matrix from config file.\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
				
				else if ((big_rows(sentinelMatrix) != big_rows(initialA)) ||
				         (big_cols(sentinelMatrix) != big_cols(initialA)))
				{
					fprintf(stderr, "Given sentinel matrix doesn't match dimensions of update matrix.\n");
					returnvalue = EXIT_SUCCESS;
					goto FREE_VARIABLES;
				}
			}
	
			//Initialise list of all moduli we'll need to use
			modList[0] = empty_BigIntT(1);
			copy_BigIntT(baseMod, modList[0]);
			for (int i = 1; i < maxPower; i += 1)
			{
				modList[i] = empty_BigIntT(1);
				multiply_BigIntT(modList[i-1], baseMod, modList[i]);
			}
			
			matrixLiftElements = malloc((maxPower-1)*sizeof(BigIntTP**));
			for (int i = 0; i < maxPower-1; i += 1)
			{
				matrixLiftElements[i] = new_BigIntT_array(big_rows(initialA), big_rows(initialA));
				for (int row = 0; row < big_rows(initialA); row += 1)
					for (int col = 0; col < big_cols(initialA); col += 1)
					{
						if (resumeMatrix == NULL)
							copy_BigIntT(big_element(initialA, row, col), matrixLiftElements[i][row][col]);
						else
							mod_BigIntT(big_element(resumeMatrix, row, col), modList[i+1], matrixLiftElements[i][row][col]);
					}
			}
			vectorLiftElements = new_BigIntT_array(big_rows(initialA), 1);
			currMat    = new_BigIntMatrixT(big_rows(initialA), big_rows(initialA));
			currPowMat = new_BigIntMatrixT(big_rows(initialA), big_rows(initialA));
			
			orbitReps      = calloc(maxPower, sizeof(BigIntMatrixTP*));
			orbitRepsCount = calloc(maxPower, sizeof(int));
			orbitRepsCycleLengths = calloc(maxPower, sizeof(int*));
			
			for (int i = 0; i < numOfTempVects; i += 1)
				tempVects[i] = new_BigIntMatrixT(big_rows(initialA), 1);
			
			for (int t = 0; t < numOfTemps; t += 1)
				temps[t] = empty_BigIntT(1);
			
			//Macro for printing out header information to the output file
			#define CREATEHEADER(OF, printBase) \
			if (printBase) \
			{ \
				fprintf(OF, "base:\n"); \
				for (int row = 0; row < big_rows(initialA); row += 1) \
				{ \
					for (int col = 0; col < big_rows(initialA); col += 1) \
					{ \
						if (col != 0) \
							fprintf(OF, " "); \
						fprinti(OF, matrixLiftElements[0][row][col]); \
					} \
					fprintf(OF, "\n"); \
				} \
			} \
			if (resumeMatrix != NULL) \
			{ \
				fprintf(OF, "start:\n");      \
				fprintbm(OF, resumeMatrix);   \
			} \
			if (sentinelMatrix != NULL) \
			{ \
				fprintf(OF, "end:\n");        \
				fprintbm(OF, sentinelMatrix); \
			} \
			fprintf(OF, "~\n")
			
			if (outputFile != NULL)
			{
				CREATEHEADER(outputFile, TRUE);
			}
			
			//Check to see if our given matrix is invertible
			//If it is, there's a whole section of code we can skip
			Ainv = big_inverse(initialA, baseMod);
			if (Ainv != NULL)
			{
				printf("Given matrix is invertible. Certain computations will be skipped to improve performance!\n");
				isAInvertible = TRUE;
			}
			else
				isAInvertible = FALSE;
			
			Ainv = free_BigIntMatrixT(Ainv);
			
			if (quietmode)
				printf("Quiet mode is on. Only select data will be output.\n");
			
			printf("Mapping down from mod ");
			printi(baseMod);
			printf("^%d...\n", maxPower);
			
			//Getting the progress variables ready
			progressNumerator   = empty_BigIntT(1);
			progressDenominator = new_BigIntT(numArr, 1);
			
			//progressDenominator = p^{L*L}
			for (int i = 0; i < big_rows(initialA)*big_rows(initialA); i += 1)
			{
				multiply_BigIntT(progressDenominator, baseMod, temps[0]);
				copy_BigIntT(temps[0], progressDenominator);
			}
			
			//Let's calculate our base orbit reps, since those won't change throughout
			// the course of the program running
			orbitRepsCount[0] = big_orbit_reps(initialA, baseMod, orbitReps, orbitRepsCycleLengths);
			
			//Now, let's start generating some orbit maps
			currPower += 1;
			while ((currPower > 1) && (!compare_BigIntMatrixT(sentinelMatrix, currMat)))
			{
				set_big_matrix(currMat, matrixLiftElements[maxPower-2]);
				if (!quietmode)
				{
					printf("M\n");
					printbm(currMat);
					if (outputFile != NULL)
					{
						fprintf(outputFile, "M\n");
						fprintbm(outputFile, currMat);
					}
				}
				
				//Clear out all vectors from mod baseMod^currPower upward
				for (int L = currPower-1; L < maxPower; L += 1)
				{
					if (orbitReps[L] != NULL)
					{
						for (int v = 0; v < orbitRepsCount[L]; v += 1)
						{
							orbitReps[L][v] = free_BigIntMatrixT(orbitReps[L][v]);
							orbitRepsCycleLengths[L][v] = 0;
						}
						FREE(orbitReps[L]);
						FREE(orbitRepsCycleLengths[L]);
						orbitRepsCount[L] = 0;
					}
				}
				
				//Get our orbit reps
				do
				{
					//Iterate through all orbit reps from lower modulus
					for (int lowRep = 0; lowRep < orbitRepsCount[currPower-2]; lowRep += 1)
					{
						//Getting ready to check every possible lift vector
						for (int i = 0; i < big_rows(initialA); i += 1)
							copy_BigIntT(big_element(orbitReps[currPower-2][lowRep], i, 0), vectorLiftElements[i][0]);
						
						//Check every lift vector
						do
						{
							set_big_matrix(tempVects[0], vectorLiftElements);
							big_floyd(currMat, tempVects[0], modList[currPower-1], &theCycle);
							copy_BigIntMatrixT(rep(theCycle), tempVects[0]);
							
							//Now, iterate until tempVects[0] is a lift of orbitReps[currPower-2][lowRep]
							//This is guaranteed to happen as orbitReps[currPower-2][lowRep] is in a cycle
							
							//This part of the code isn't necessary if our matrix is invertible
							// since the rep(theCycle) will be the same as the vector we started with
							while (!isAInvertible)
							{
								copy_BigIntMatrixT(tempVects[0], tempVects[1]);
								modbm(tempVects[1], modList[currPower-2]); //Reduce to see if we have a lift of orbitReps[currPower-2][lowRep]
								
								if (! compare_BigIntMatrixT(tempVects[1], orbitReps[currPower-2][lowRep]))
								{
									big_mat_mul(currMat, tempVects[0], tempVects[1]);
									copy_BigIntMatrixT(tempVects[1], tempVects[0]);
									modbm(tempVects[0], modList[currPower-1]);
								}
								
								else
									break;
							}
							
							//Okay, so tempVects[0] is in a cycle and is a lift of orbitReps[currPower-2][lowReps]
							//Now, if we need to, we check to see if it's already accounted for by the other cycles.
							if (orbitRepsCount[currPower-1] == 0)
							{
								orbitRepsCount[currPower-1] = 1;
								orbitReps[currPower-1] = realloc(orbitReps[currPower-1], sizeof(BigIntMatrixTP));
								orbitReps[currPower-1][0] = new_BigIntMatrixT(big_rows(initialA), 1);
								copy_BigIntMatrixT(tempVects[0], orbitReps[currPower-1][0]);
								
								orbitRepsCycleLengths[currPower-1] = realloc(orbitRepsCycleLengths[currPower-1], sizeof(int));
								orbitRepsCycleLengths[currPower-1][0] = omega(theCycle);
							}
							
							else
							{
								isNewVector = TRUE;
								copy_BigIntMatrixT(tempVects[0], tempVects[2]);
								
								//Get matrix ready for doing quick iterations
								//Let's hope the cycle length never gets too big...
								numArr[0] = orbitRepsCycleLengths[currPower-2][lowRep];
								temps[0] = free_BigIntT(temps[0]); //Why didn't I just make a set function for BigIntTs?
								temps[0] = new_BigIntT(numArr, 1);
								powbm(currMat, currPowMat, temps[0], modList[currPower-1]);
								
								do
								{
									//Look at all the orbit reps we have so far, see if we have any duplicates
									//We only need to check specific vectors in the cycle of tempVects[0] since
									// we're only keeping track of lifts from specific vectors, so if tempVects[0]
									// is to be in another cycle, it'll have to cycle back to one of the lifts we
									// already have
									for (int otherRep = 0; otherRep < orbitRepsCount[currPower-1]; otherRep += 1)
									{
										if (compare_BigIntMatrixT(orbitReps[currPower-1][otherRep], tempVects[0]))
										{
											isNewVector = FALSE;
											break;
										}
									}
									
									//Iterate tempVects[0] to next lift of orbitReps[currPower-2][lowRep]
									big_mat_mul(currPowMat, tempVects[0], tempVects[1]);
									modbm(tempVects[1], modList[currPower-1]);
									copy_BigIntMatrixT(tempVects[1], tempVects[0]);
								}
								while ((! compare_BigIntMatrixT(tempVects[0], tempVects[2])) && (isNewVector));
								
								//Yay, we found a new rep!
								if (isNewVector)
								{
									orbitRepsCount[currPower-1] += 1;
									orbitReps[currPower-1] = realloc(orbitReps[currPower-1], orbitRepsCount[currPower-1]*sizeof(BigIntMatrixTP));
									orbitReps[currPower-1][orbitRepsCount[currPower-1]-1] = new_BigIntMatrixT(big_rows(initialA), 1);
									copy_BigIntMatrixT(tempVects[0], orbitReps[currPower-1][orbitRepsCount[currPower-1]-1]);
									
									orbitRepsCycleLengths[currPower-1] = realloc(orbitRepsCycleLengths[currPower-1], orbitRepsCount[currPower-1]*sizeof(int));
									orbitRepsCycleLengths[currPower-1][orbitRepsCount[currPower-1]-1] = omega(theCycle);
								}
							}
						}
						while (! step_BigIntT_array(vectorLiftElements, big_rows(initialA), 1, modList[currPower-2], modList[currPower-1]));
					  //I'm suspicious of the step_BigIntT_array() above me, but it should be fine...
					}
					
					//Find the orbit reps for the next modulus up
					currPower += 1;
				}
				while (currPower <= maxPower);
				currPower -= 1; //Accounting for the weird way I set up this loop
				
				//Now, print out the orbit reps (or at least, start preparing to)
				
				for (int i = 0; i < maxPower; i += 1)
				{
					stemCollection[i] = NULL;
					stemCollectionIndices[i] = 0;
				}
				
				//Template:
				//End of headers are marked with a ~
				/* M
				 * <matrix>
				 * <indent level> <vector> <cycle length>
				 * M
				 * <matrix>
				 * <indent level> <vector> <cycle length>
				 * etc.
				 */
				 
				//Create a new orbitmaptree structure to mimic the tree we're printing
				//Then, compare this structure with the ones we already have to see if it's unique
				OMT tempTree = {
												 1, 
												 new_BigIntMatrixT(big_rows(initialA), big_rows(initialA)), 
												 calloc(maxPower, sizeof(int)), 
												 calloc(maxPower, sizeof(OMN*))
											 };
				//Put the current matrix into the mapping to act as a representative
				copy_BigIntMatrixT(currMat, tempTree.matrixRep);

				while (stemCollectionIndices[maxPower-1] < orbitRepsCount[maxPower-1])
				{
					stemCollection[maxPower-1] = orbitReps[maxPower-1][stemCollectionIndices[maxPower-1]];
					stemCollectionCycleLengths[maxPower-1] = orbitRepsCycleLengths[maxPower-1][stemCollectionIndices[maxPower-1]];
					stemCollectionIndices[maxPower-1] += 1;
					
					//Reduce this vector down to see if it matches with a lower vector on the tree
					// If not, we'll have to start a new branch
					for (int stem = maxPower-2; stem >= 0; stem -= 1)
					{
						stemStart = -1;
						copy_BigIntMatrixT(stemCollection[maxPower-1], tempVects[0]);
						modbm(tempVects[0], modList[stem]);
						
						if (! compare_BigIntMatrixT(tempVects[0], stemCollection[stem]))
						{
							stemCollection[stem] = orbitReps[stem][stemCollectionIndices[stem]];
							stemCollectionCycleLengths[stem] = orbitRepsCycleLengths[stem][stemCollectionIndices[stem]];
							stemCollectionIndices[stem] += 1;
								
							if (stem == 0)
								stemStart = 0;
						}
						
						else
							stemStart = stem+1;

						//Now, actually print out this tree structure if there's something to print
						if (stemStart != -1)
						{
							for (int v = stemStart; v < maxPower; v += 1)
							{
								if (!quietmode)
								{
									if (useTreeSpaces)
										for (int i = 0; i < v; i += 1)
											printf(" ");
										
									else
										printf("%d ", v);
								
									for (int i = 0; i < big_rows(initialA); i += 1)
									{
										if (i != 0)
											printf(",");
										printi(big_element(stemCollection[v], i, 0));
									}
									printf(" %d\n", stemCollectionCycleLengths[v]);
								}
								
								//Make space in orbitmaptree for vector
								tempTree.layerCount[v] += 1;
								tempTree.nodes[v] = realloc(tempTree.nodes[v], tempTree.layerCount[v]*sizeof(OMN));
								
								//Add vector
								OMN tempNode = {0, stemCollectionCycleLengths[v]};
								tempTree.nodes[v][tempTree.layerCount[v]-1] = tempNode;
								
								//Increase number of children for parent node, if it exists
								if (v > 0)
									tempTree.nodes[v-1][tempTree.layerCount[v-1]-1].subNodes += 1;
								
								if (outputFile != NULL)
								{
									if (useTreeSpaces)
										for (int i = 0; i < v; i += 1)
											fprintf(outputFile, " ");
										
									else
										fprintf(outputFile, "%d ", v);
									
									for (int i = 0; i < big_rows(initialA); i += 1)
									{
										if (i != 0)
											fprintf(outputFile, ",");
										fprinti(outputFile, big_element(stemCollection[v], i, 0));
									}
									fprintf(outputFile, " %d\n", stemCollectionCycleLengths[v]);
								}
							}
							break;
						}
					}
				}
				
				//Here's the place where we check whether the tree we just built is unique.
				//Thanks, past me, for adding this comment. It is actually very useful
				if (orbitMapsCount == 0)
				{
					isNewTree = TRUE;
					orbitMapsCount += 1;
					orbitMapsCatalogue = malloc(sizeof(OMT));
					orbitMapsCatalogue[0] = tempTree;
					tempTree.occurrences += 1;
					
					if (!quietmode)
						printf("Belongs to mapping 0.\n");

					if (outputFile != NULL)
						fprintf(outputFile, "Belongs to mapping 0.\n");
				}
				
				//Compare our most recent tree with all the other trees we have
				//If it's unique, add it
				else
				{
					isNewTree = TRUE;
					for (int tree = 0; tree < orbitMapsCount; tree += 1)
					{
						isNewTree = FALSE;
						
						//Check to see if the layer counts are equal
						for (int LC = 0; LC < maxPower; LC += 1)
						{
							if (orbitMapsCatalogue[tree].layerCount[LC] != tempTree.layerCount[LC])
							{
								isNewTree = TRUE;
								break;
							}
						}
						
						//If the layer counts were the same
						if (!isNewTree)
						{
							//isNewTree = FALSE;
							for (int LC = 0; LC < maxPower; LC += 1)
							{
								//Here, nodes are copied into nodeTracker BY VALUE
								//If I wanted to speed this code up, I could make nodeTracker
								// hold references to those nodes
								nodeTracker = realloc(nodeTracker, orbitMapsCatalogue[tree].layerCount[LC]*sizeof(OMT));
								for (int N = 0; N < orbitMapsCatalogue[tree].layerCount[LC]; N += 1)
									nodeTracker[N] = orbitMapsCatalogue[tree].nodes[LC][N];
								
								//Now, we make sure that each node in tempTree has a representative in nodeTracker
								for (int N = 0; N < tempTree.layerCount[LC]; N += 1)
								{
									foundMatchingNode = FALSE;
									for (int matchN = 0; matchN < orbitMapsCatalogue[tree].layerCount[LC]; matchN += 1)
									{
										if ((nodeTracker[matchN].subNodes == tempTree.nodes[LC][N].subNodes) &&
												(nodeTracker[matchN].cyclen   == tempTree.nodes[LC][N].cyclen))
										{
											foundMatchingNode = TRUE;
											nodeTracker[matchN].cyclen = -1;
											break;
										}
									}
									
									//If there was at least one node in our new tree that didn't match
									// the tree we're looking at in the catalogue
									if (!foundMatchingNode)
									{
										isNewTree = TRUE;
										break;
									}
								}
								
								//If we know tempTree isn't this tree, check the next
								if (isNewTree)
									break;
							}
						}
						
						//If we found a matching tree to tempTree in orbitMapsCatalogue, we don't need to look any further
						if (!isNewTree)
						{
							orbitMappingID = tree;
							//Add +1 occurrence to the tree structure that's the same as tempTree
							orbitMapsCatalogue[tree].occurrences += 1;
							break;
						}
					}
					
					//If we found a new tree to add
					if (isNewTree)
					{
						if (!quietmode)
							printf("Belongs to mapping %d.\n", orbitMapsCount);
						
						if (outputFile != NULL)
							fprintf(outputFile, "Belongs to mapping %d.\n", orbitMapsCount);
						
						orbitMapsCount += 1;
						orbitMapsCatalogue = realloc(orbitMapsCatalogue, orbitMapsCount*sizeof(OMT));
						orbitMapsCatalogue[orbitMapsCount-1] = tempTree;
						tempTree.occurrences += 1;
					}
					
					//Discard tree
					else
					{
						if (!quietmode)
							printf("Belongs to mapping %d.\n", orbitMappingID);
						
						if (outputFile != NULL)
							fprintf(outputFile, "Belongs to mapping %d.\n", orbitMappingID);
						
						for (int layer = 0; layer < maxPower; layer += 1)
						{
							FREE(tempTree.nodes[layer]);
						}
						tempTree.matrixRep = free_BigIntMatrixT(tempTree.matrixRep);
						FREE(tempTree.nodes);
						FREE(tempTree.layerCount);
					}
				}
				
				//Save datafile
				if (outputFile != NULL)
				{
					if (fclose(outputFile) == EOF)
						fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
					
					outputFile = fopen(outputFileName, "a");
					if (outputFile == NULL)
						fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
				}
				
				//Increment through our matrix lifts in such a way as to reduce
				// the amount of times we need to recompute lifts
				for (int level = maxPower-2; level >= 0; level -= 1)
				{
					//Free these reps since they'll be replaced for the next lift we check
					for (int v = 0; v < orbitRepsCount[level+1]; v += 1)
					{
						orbitReps[level+1][v] = free_BigIntMatrixT(orbitReps[level+1][v]);
						orbitRepsCycleLengths[level+1][v] = 0;
					}
					FREE(orbitReps[level+1]);
					FREE(orbitRepsCycleLengths[level+1]);
					orbitRepsCount[level+1] = 0;
					
					//Print progress report to stdout
					if (level == 0)
					{
						add_BigIntT(one, progressNumerator, temps[0]);
						copy_BigIntT(temps[0], progressNumerator);
						printi(progressNumerator);
						printf(" / ");
						printi(progressDenominator);
						printf(" computed...\n");
					}
					
					//If we're done looking at all the matrix lifts that keep the
					// lower vector lifts the same
					if (step_BigIntT_array(matrixLiftElements[level], 
																 big_rows(initialA), 
																 big_rows(initialA), 
																 modList[level], 
																 modList[level+1]))
					{
						currPower -= 1;
					}
					
					else
						break;
				}
				
				//Now, make sure that if we changed a lower-moduli matrix lift, we carry that
				// new lift up to the higher matrix lifts so the reduction process that keeps
				// lift vectors constant works properly.
				if (currPower > 1)
					for (int rewrite = 1; rewrite + currPower <= maxPower; rewrite += 1)
						for (int row = 0; row < big_rows(initialA); row += 1)
							for (int col = 0; col < big_rows(initialA); col += 1)
								copy_BigIntT(matrixLiftElements[currPower-2][row][col], 
							               matrixLiftElements[currPower-2+rewrite][row][col]);
				
				//If it's time to create a new output file (split when the lowest
				// lift gets incremented)
				if ((currPower == 2) && (outputFile != NULL))
				{
					if (fclose(outputFile) == EOF)
						fprintf(stderr, "Unable to save output file. Will attempt to save next file created...\n");
					
					fileCounter += 1;
					strcpy(outputFileName, outputFileBaseName);
					append_int(outputFileName, fileCounter);
					strcat(outputFileName, ".txt");
					
					outputFile = fopen(outputFileName, "a");
					if (outputFile == NULL)
						fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
					
					//If file created successfully, make sure to add the file header to the top
					else
					{
						CREATEHEADER(outputFile, TRUE);
					}
				}
			}
			
			//Now, output the tree structures to a file, if desired
			if (mapOutputFileName != NULL)
			{
				mapOutputFile = fopen(mapOutputFileName, "a");
				if (mapOutputFile == NULL)
					fprintf(stderr, "Unable to create output file for saving unique orbit map structures.\n");
				
				else
				{
					printf("Unique mapping structures:\n");
					CREATEHEADER(mapOutputFile, FALSE);
					
					for (int T = 0; T < orbitMapsCount; T += 1)
					{
						for (int i = 0; i < maxPower; i += 1)
						{
							if (i < maxPower-1)
								subCount[i] = 0;
							mapIndexCount[i] = 0;
						}
						
						fprintf(mapOutputFile, ":%d (%d)\n", T, orbitMapsCatalogue[T].occurrences);
						printf(":%d (%d)\n", T, orbitMapsCatalogue[T].occurrences);
						
						//Print the representative matrix for the mapping
						fprintf(mapOutputFile, "M\n");
						printf("M\n");
						fprintbm(mapOutputFile, orbitMapsCatalogue[T].matrixRep);
						printbm(orbitMapsCatalogue[T].matrixRep);
						printf("~\n");
						fprintf(mapOutputFile, "~\n");
						do
						{
							//Generate the number of nodes we need to print in each layer for the current layer 0 node
							for (int L = 0; L < maxPower-1; L += 1)
								subCount[L] = orbitMapsCatalogue[T].nodes[L][mapIndexCount[L]].subNodes;
							
							//Now, iterate through the tree
							fprintf(mapOutputFile, "0 %d\n", orbitMapsCatalogue[T].nodes[0][mapIndexCount[0]].cyclen);
							printf("0 %d\n", orbitMapsCatalogue[T].nodes[0][mapIndexCount[0]].cyclen);
							currentLayer = 1;
							do
							{
								if (subCount[currentLayer-1] > 0)
								{
									fprintf(mapOutputFile, "%d %d\n", currentLayer, 
									                                  orbitMapsCatalogue[T].nodes[currentLayer][mapIndexCount[currentLayer]].cyclen);
									printf("%d %d\n", currentLayer, orbitMapsCatalogue[T].nodes[currentLayer][mapIndexCount[currentLayer]].cyclen);
									subCount[currentLayer-1] -= 1; //currentLayer should never be zero here due to the while condition
									mapIndexCount[currentLayer] += 1;
									
									//Keep going down until we find the lowest layer of the tree to consecutively print
									if (currentLayer < maxPower-1)
										currentLayer += 1;
								}
								
								//Recalculate subCount, if needed
								else
								{
									currentLayer -= 1;
									if ((currentLayer > 0) && (subCount[currentLayer-1] > 0))
										for (int L = currentLayer; L < maxPower-1; L += 1)
											subCount[L] = orbitMapsCatalogue[T].nodes[L][mapIndexCount[L]].subNodes;
										
									//We don't need to increment mapIndexCount[currentLayer] here as 
									// it'll be incremented above in the next loop
								}
							}
							while (currentLayer > 0);
							mapIndexCount[0] += 1;
						}
						while (mapIndexCount[maxPower-1] < orbitMapsCatalogue[T].layerCount[maxPower-1]);
					}
				}
				if (fclose(mapOutputFile) == EOF)
					fprintf(stderr, "Unable to save unique orbit map structures to output file.\n");
			}
			
			else
			{
				//Print out mapping structures, just in a different way
				printf("orbitMapsCatalogue:\n");
				for (int T = 0; T < orbitMapsCount; T += 1)
				{
					printf("%d (%d):\n", T, orbitMapsCatalogue[T].occurrences);
					printf("M\n");
					printbm(orbitMapsCatalogue[T].matrixRep);
					for (int L = 0; L < maxPower; L += 1)
					{
						printf("~Nodes for mod ");
						printi(baseMod);
						printf("^%d (%d)~", L+1, orbitMapsCatalogue[T].layerCount[L]);
						for (int N = 0; N < orbitMapsCatalogue[T].layerCount[L]; N += 1)
						{
							if (N != 0)
								printf(", ");
							if (N % 8 == 0) //This breaks up each line into more manageable rows
								printf("\n");
							
							printf("{%d, %d}", orbitMapsCatalogue[T].nodes[L][N].subNodes, 
																 orbitMapsCatalogue[T].nodes[L][N].cyclen);
						}
						printf("\n");
					}
					printf("\n");
				}
			}
			
			//FREEEVERYTHING:
			if (outputFile != NULL)
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Unable to save output file.\n");
			FREE(outputFileName);
			FREE(mapOutputFileName);
			FREE(outputFileBaseName);
			
			resumeMatrix   = free_BigIntMatrixT(resumeMatrix);
			sentinelMatrix = free_BigIntMatrixT(sentinelMatrix);
			
			for (int i = 0; i <= maxargc; i += 1)
			{
				FREE(oargv[i]);
			}
			FREE(oargv);
			
			baseMod = free_BigIntT(baseMod);
			
			progressNumerator   = free_BigIntT(progressNumerator);
			progressDenominator = free_BigIntT(progressDenominator);
			
			for (int t = 0; t < numOfTemps; t += 1)
				temps[t] = free_BigIntT(temps[t]);
			
			for (int m = 0; m < maxPower; m += 1)
				modList[m] = free_BigIntT(modList[m]);
			
			for (int i = 0; i < maxPower-1; i += 1)
				matrixLiftElements[i] = free_BigIntT_array(matrixLiftElements[i], big_rows(initialA), big_rows(initialA));
			FREE(matrixLiftElements);
			
			vectorLiftElements = free_BigIntT_array(vectorLiftElements, big_rows(initialA), 1);
			currMat    = free_BigIntMatrixT(currMat);
			currPowMat = free_BigIntMatrixT(currPowMat);
			
			if (orbitMapsCatalogue != NULL)
			{
				for (int om = 0; om < orbitMapsCount; om += 1)
				{
					for (int layer = 0; layer < maxPower; layer += 1)
					{
						FREE(orbitMapsCatalogue[om].nodes[layer]);
					}
					orbitMapsCatalogue[om].matrixRep = free_BigIntMatrixT(orbitMapsCatalogue[om].matrixRep);
					FREE(orbitMapsCatalogue[om].nodes);
					FREE(orbitMapsCatalogue[om].layerCount);
				}
				FREE(orbitMapsCatalogue);
			}
			
			for (int i = 0; i < maxPower; i += 1)
			{
				for (int j = 0; j < orbitRepsCount[i]; j += 1)
					orbitReps[i][j] = free_BigIntMatrixT(orbitReps[i][j]);
				FREE(orbitReps[i]);
				FREE(orbitRepsCycleLengths[i]);
			}
			FREE(orbitReps);
			FREE(orbitRepsCount);
			FREE(orbitRepsCycleLengths);
			
			theCycle  = free_CycleInfoT(theCycle);
			for (int i = 0; i < numOfTempVects; i += 1)
				tempVects[i] = free_BigIntMatrixT(tempVects[i]);
			
			FREE(nodeTracker);
		}
		
		
		//If the user wants to find an annihilating polynomial which
		// isn't a multiple of the minimal polynomial
		else if (! strcmp(argv[1], "oddminpolysearch"))
		{
			//Checking to make sure user entered all required arguments
			for (int i = 2; i <= 4; i += 1)
				if (oargv[i][0] == '\0')
				{
					printf("Usage: " ANSI_COLOR_YELLOW "oddminpolysearch " ANSI_COLOR_CYAN \
					"maxpower size polysize [modulus] [resume] [fileoutput]" ANSI_COLOR_RESET "\n");
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
			
			int matSize, maxPower, currPower, polySize;
			BigIntTP bigMod = NULL;
			BigIntTP currMod = NULL;
			BigIntTP** currMatElements = NULL;
			BigIntMatrixTP currMat = NULL;
			
			bool foundMin, foundMonic;
			BigIntTP* minPolyCoeffs   = NULL;
			BigPolyTP minPoly         = NULL;
			BigPolyTP minMonicPoly    = NULL; //Non-monic annihilating polys can be smaller than minPoly
			BigIntTP* annilPolyCoeffs = NULL;
			BigPolyTP annilPoly       = NULL;
			
			BigPolyTP remainderPoly = NULL;
			BigPolyTP zeroPoly      = NULL;
			
			BigIntTP temp, coeffCompare;
			
			BigIntMatrixTP zeroMat = NULL;
			BigIntMatrixTP tempMat = NULL;
			
			char* outputFileName = NULL;
			FILE* outputFile = NULL;
			
			maxPower = (int)strtol(oargv[2], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read maxPower from command line.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			matSize = (int)strtol(oargv[3], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read matrix size from command line.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			polySize = (int)strtol(oargv[4], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read polynomial size from command line.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (oargv[5][0] != '\0')
			{
				SET_BIG_NUM(oargv[5], bigMod, "Unable to read modulus from command line.");
			}
			else
			{
				SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			}
			
			//Check to see if the user wants to resume computation
			if (!strcmp(oargv[6], "TRUE"))
			{
				currMat = read_BigIntMatrixT(resumefilepath);
				if (currMat == NULL)
				{
					fprintf(stderr, "Unable to read resume matrix from %s.\n", resumefilepath);
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
				
				if ((big_rows(currMat) != big_cols(currMat)) || (big_rows(currMat) != matSize))
				{
					fprintf(stderr, "Resume matrix given isn't the correct dimensions specified.\n");
					currMat = free_BigIntMatrixT(currMat);
					bigMod = free_BigIntT(bigMod);
					
					returnvalue = EXIT_FAILURE;
					goto FREE_VARIABLES;
				}
				
				currMatElements = malloc(matSize*sizeof(BigIntTP*));
				for (int i = 0; i < matSize; i += 1)
				{
					currMatElements[i] = malloc(matSize*sizeof(BigIntTP));
					for (int j = 0; j < matSize; j += 1)
					{
						currMatElements[i][j] = empty_BigIntT(1);
						copy_BigIntT(big_element(currMat, i, j), currMatElements[i][j]);
					}
				}
			}
			
			else
			{
				currMat = new_BigIntMatrixT(matSize, matSize);
				
				currMatElements = malloc(matSize*sizeof(BigIntTP*));
				for (int i = 0; i < matSize; i += 1)
				{
					currMatElements[i] = malloc(matSize*sizeof(BigIntTP));
					for (int j = 0; j < matSize; j += 1)
						currMatElements[i][j] = empty_BigIntT(1);
				}
			}
			
			//If the user wants a file output
			if (!strcmp(oargv[7], "TRUE"))
			{
				outputFileName = malloc(MAXSTRLEN*sizeof(char));
				outputFileName[0] = '\0';
				
				ADD_FILE_PREFIX(outputFileName);
				
				//Let's hope there's enough space in the string for all this!
				strcat(outputFileName, "oddminpolysearch ");
				strcat(outputFileName, oargv[2]);
				strcat(outputFileName, " ");
				strcat(outputFileName, oargv[3]);
				strcat(outputFileName, " ");
				strcat(outputFileName, oargv[4]);
				strcat(outputFileName, " ");
				append_BigIntT(outputFileName, bigMod);
				strcat(outputFileName, ".txt");	
				
				outputFile = fopen(outputFileName, "w");
				if (outputFile == NULL)
					fprintf(stderr, "Unable to create file %s. Continuing without saving...\n", outputFileName);
			}
			
			temp = empty_BigIntT(1);
			
			zeroPoly = constant_BigPolyT(temp);
			
			currMod = empty_BigIntT(1);
			copy_BigIntT(bigMod, currMod);
			
			minPolyCoeffs   = malloc(polySize*sizeof(BigIntTP));
			annilPolyCoeffs = malloc(polySize*sizeof(BigIntTP));
			for (int i = 0; i < polySize; i += 1)
			{
				minPolyCoeffs[i]   = empty_BigIntT(1);
				annilPolyCoeffs[i] = empty_BigIntT(1);
			}
			
			minPoly      = empty_BigPolyT();
			minMonicPoly = empty_BigPolyT();
			annilPoly    = empty_BigPolyT();
			
			tempMat = new_BigIntMatrixT(matSize, matSize);
			zeroMat = new_BigIntMatrixT(matSize, matSize);
			
			//Loops over all matrices under current modulus
			for (currPower = 1; currPower <= maxPower; currPower += 1)
			{
				if (outputFile != NULL)
				{
					fprintf(outputFile, "~~Mod ");
					fprinti(outputFile, currMod);
					fprintf(outputFile, "~~\n");
					if (fclose(outputFile) == EOF)
					{
						fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
						outputFile = NULL;
					}
					else
					{
						outputFile = fopen(outputFileName, "a");
						if (outputFile == NULL)
							fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
					}
				}
				
				do
				{
					set_big_matrix(currMat, currMatElements);
					resize_BigPolyT(minPoly, polySize);
					resize_BigPolyT(minMonicPoly, polySize);
					
					//Find minimal polynomial
					
					for (int i = 0; i < polySize; i += 1)
					{
						clear_BigIntT(minPolyCoeffs[i]);
						clear_BigIntT(annilPolyCoeffs[i]);
					}
					
					foundMin   = FALSE;
					foundMonic = FALSE;
					while ((!foundMin) || (!foundMonic))
					{
						for (int i = 0; i < polySize; i += 1)
						{
							add_BigIntT(minPolyCoeffs[i], one, temp);
							if (compare_BigIntT(temp, currMod) == 0)
								clear_BigIntT(minPolyCoeffs[i]);
							else
							{
								copy_BigIntT(temp, minPolyCoeffs[i]);
								break;
							}
						}
						
						if (!foundMin)
						{
							set_BigPolyT(minPoly, minPolyCoeffs);
							set_BigPolyT(minMonicPoly, minPolyCoeffs);
							reduce_BigPolyT(minPoly);
							eval_BigPolyT(minPoly, currMat, tempMat, currMod);
							
							resize_BigPolyT(minPoly, polySize);
						}
						else if (!foundMonic)
						{
							set_BigPolyT(minMonicPoly, minPolyCoeffs);
							reduce_BigPolyT(minMonicPoly);
							eval_BigPolyT(minMonicPoly, currMat, tempMat, currMod);
							
							resize_BigPolyT(minMonicPoly, polySize);
						}
						
						//If we've found our minimum polynomial for currMat
						if (compare_BigIntMatrixT(tempMat, zeroMat))
						{
							if (!foundMin)
							{
								reduce_BigPolyT(minPoly);
								foundMin = TRUE;
							}
							
							//Check to see if polynomial is monic
							if (!foundMonic)
							{
								for (int i = polySize-1; i >= 0; i -= 1)
								{
									if (!is_zero(minPolyCoeffs[i]))
									{
										if (compare_BigIntT(one, minPolyCoeffs[i]) == 0)
										{
											reduce_BigPolyT(minMonicPoly);
											foundMonic = TRUE;
										}
										
										break;
									}
								}
							}
						}
					}
					
					//For checking annihilating polynomials, we'll only check
					// monic polynomials of size polySize
					for (int c = 0; c < polySize-1; c += 1)
						copy_BigIntT(currMod, annilPolyCoeffs[c]);
					
					//We've found our minimal polynomial(s)
					//Now to see if all other annihilating polynomials are
					// multiples of them.
					
					//foundMin and foundMonic are now gonna be used to keep track of
					// whether we've found a polynomial that's annihilating
					foundMin = FALSE;
					foundMonic = FALSE;
					
					remainderPoly = empty_BigPolyT();
					coeffCompare = empty_BigIntT(1);
					
					//First, we should check whether our first monic polynomial divides the first annihilating polynomial
					copy_BigPolyT(minMonicPoly, annilPoly);
					while ((!foundMin) || (!foundMonic))
					{
						//Check to see if annilPoly is an annihilating polynomial
						eval_BigPolyT(annilPoly, currMat, tempMat, currMod);
						if (compare_BigIntMatrixT(tempMat, zeroMat))
						{
							if (!foundMin)
							{
								//Check if annilPoly is a multiple of minPoly
								divide_BigPolyT(annilPoly, minPoly, NULL, remainderPoly, currMod);
								if (compare_BigPolyT(remainderPoly, zeroPoly) != 0)
									foundMin = TRUE;
							}
							
							if (!foundMonic)
							{
								//Check if annilPoly is a multiple of minMonicPoly
								divide_BigPolyT(annilPoly, minMonicPoly, NULL, remainderPoly, currMod);
								if (compare_BigPolyT(remainderPoly, zeroPoly) != 0)
									foundMonic = TRUE;
							}
							
							if ((foundMin) && (foundMonic))
							{
								printf("\ncurrMat:\n");
								printbm(currMat);
								printf("Minimal poly: ");
								printp(minPoly);
								printf("\nMonic poly: ");
								printp(minMonicPoly);
								printf("\n");
								printp(annilPoly);
								printf(" is not a multiple of ");
								printp(minPoly);
								printf("\n");
								printp(annilPoly);
								printf(" is not a multiple of ");
								printp(minMonicPoly);
								printf("\n");
								
								if (outputFile != NULL)
								{
									fprintbm(outputFile, currMat);
									fprintf(outputFile, "Minimal poly: ");
									fprintp(outputFile, minPoly);
									fprintf(outputFile, "\nMonic poly: ");
									fprintp(outputFile, minMonicPoly);
									fprintf(outputFile, "\n");
									fprintp(outputFile, minPoly);
									fprintf(outputFile, " !| ");
									fprintp(outputFile, annilPoly);
									fprintf(outputFile, "\n");
									fprintp(outputFile, minMonicPoly);
									fprintf(outputFile, " !| ");
									fprintp(outputFile, annilPoly);
									fprintf(outputFile, "\n\n");
									
									if (fclose(outputFile) == EOF)
									{
										fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
										outputFile = NULL;
									}
									else
									{
										outputFile = fopen(outputFileName, "a");
										if (outputFile == NULL)
											fprintf(stderr, "Unable to save output file. Continuing without saving...\n");
									}
								}
							}
						}
						
						//Increment annilPolyCoeffs
						for (int i = 0; i < polySize; i += 1)
						{
							add_BigIntT(annilPolyCoeffs[i], one, temp);
							if (i != polySize-1)
								copy_BigIntT(currMod, coeffCompare);
							else
								add_BigIntT(one, one, coeffCompare);
							
							if (compare_BigIntT(temp, coeffCompare) >= 0)
							{
								clear_BigIntT(annilPolyCoeffs[i]);
								if (i == polySize-1)
								{
									foundMin = TRUE;
									foundMonic = TRUE;
								}
							}
							else
							{
								copy_BigIntT(temp, annilPolyCoeffs[i]);
								break;
							}
						}
						
						//Geez, I love spending CPU cycles on this :)
						resize_BigPolyT(annilPoly, polySize);
						set_BigPolyT(annilPoly, annilPolyCoeffs);
						reduce_BigPolyT(annilPoly);
					}
					
					remainderPoly = free_BigPolyT(remainderPoly);
					coeffCompare  = free_BigIntT(coeffCompare);
				}
				while (!increment_BigIntT_array(currMatElements, matSize, matSize, one, currMod));
				
				//Get next power of the modulus
				multiply_BigIntT(currMod, bigMod, temp);
				copy_BigIntT(temp, currMod);
			}
			
			if (outputFile != NULL)
			{
				if (fclose(outputFile) == EOF)
					fprintf(stderr, "Output file may not have saved correctly.\n");
			}
			
			bigMod  = free_BigIntT(bigMod);
			currMatElements = free_BigIntT_array(currMatElements, matSize, matSize);
			currMat = free_BigIntMatrixT(currMat);
			
			temp         = free_BigIntT(temp);
			currMod      = free_BigIntT(currMod);
			//coeffCompare = free_BigIntT(coeffCompare);
			
			tempMat = free_BigIntMatrixT(tempMat);
			zeroMat = free_BigIntMatrixT(zeroMat);
			
			for (int i = 0; i < polySize; i += 1)
			{
				minPolyCoeffs[i]   = free_BigIntT(minPolyCoeffs[i]);
				annilPolyCoeffs[i] = free_BigIntT(annilPolyCoeffs[i]);
			}
			FREE(minPolyCoeffs);
			FREE(annilPolyCoeffs);
			
			minPoly      = free_BigPolyT(minPoly);
			minMonicPoly = free_BigPolyT(minMonicPoly);
			annilPoly    = free_BigPolyT(annilPoly);
			
			//remainderPoly = free_BigPolyT(remainderPoly);
			zeroPoly      = free_BigPolyT(zeroPoly);
			
			if (outputFileName != NULL)
			{
				FREE(outputFileName);
			}
		}
		
		
		//If we want to find a lift of a matrix where the cycle length stays the same
		else if (!strcmp(argv[1], "stablelift"))
		{
			bool debugMode = FALSE;
			
			/* This CLI tool might be able to be used to find Wall-Sun-Sun primes.
			 * A Wall-Sun-Sun prime is one where, when finding a stable lift of the
			 * 2x2 Fibonacci matrix, the stable lift is the same as the 2x2 Fibonacci
			 * matrix...
			 *
			 */
			if ((oargv[2][0] == '\0') || (oargv[3][0] == '\0'))
			{
				printf("Usage: " ANSI_COLOR_YELLOW "stablelift " ANSI_COLOR_CYAN "baseMod modPower" ANSI_COLOR_RESET "\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (big_rows(UPDATEMATRIX) != big_cols(UPDATEMATRIX))
			{
				fprintf(stderr, "Update matrix must be square.\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			//Says whether the program needs to be run again to find a proper stable lift
			bool isIncompleteSolution = FALSE;
			
			int numArr[1] = {1};
			int modPower;
			BigIntTP baseMod = NULL;
			BigIntTP bigMod = NULL;
			BigIntTP biggerMod = NULL;
			
			int cycLen;
			CycleInfoTP theCycle = NULL;
			BigIntMatrixTP I;
			
			BigIntTP** matrixElements = NULL;
			BigIntTP** matrixDeltaRows = NULL; //Holds the changes to the matrix in row form for row reduction
			int deltaRowCounter = 0;
			int deltaColCounter = 0;
			
			BigIntTP bigCycLen;
			BigIntMatrixTP basePDigitMatrix; //UPDATE^cycLen mod baseMod^(modPower + 1)
			BigIntMatrixTP tempMat = NULL;
			BigIntMatrixTP tempMat2 = NULL;
			BigIntMatrixTP tempPDigitMatrix = NULL; //<modified UPDATE>^cycLen mod baseMod^(modPower + 1)
			
			BigIntTP* baseDeltaRow; //The row we're trying to eliminate
			
			//Keeps track of what shifts we need to perform on UPDATEMATRIX to get a stable lift
			BigIntMatrixTP matrixShiftWatcher = NULL; 
			BigIntTP* matrixShiftWatcherRow = NULL;
			
			BigIntMatrixTP deltaRowsToReduce = NULL;
			BigIntMatrixTP deltaRowsToReduceTemp = NULL;
			
			//Holds the combination of matrix terms additions we need to get a stable lift
			BigPolyTP deltaCombo = NULL;
			BigIntTP* deltaComboCoeffs;
			int matrixIndexCounter = 0; //Helps to index matrices as 1D arrays
			
			BigIntTP temp, temp2, temp3, temp4; //Insanity
			BigIntTP negOne = NULL;
			
			//Tells the algorithm how many levels deep UPDATEMATRIX is as a 
			// stable lift (is it the product of a previous stable lift computation?)
			int stableLiftLevel = -1;
			BigIntTP deltaMod = NULL;
			BigIntTP liftCompMod = NULL;
			//My usage of deltaMod and liftCompMod may be completely wrong
			// We're gonna find out
			
			SET_BIG_NUM(oargv[2], baseMod, "Unable to read base of modulus from command line.");
			
			modPower = (int)strtol(oargv[3], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read modulus power from command line.\n");
				baseMod = free_BigIntT(baseMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			//Now, let's automatically compute what level of stable lift the given matrix is
			I     = identity_BigIntMatrixT(big_rows(UPDATEMATRIX));
			temp  = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			temp3 = empty_BigIntT(1);
			temp4 = empty_BigIntT(1);
			
			copy_BigIntT(baseMod, temp);
			big_floyd(UPDATEMATRIX, I, temp, &theCycle);
			
			do
			{
				//How many powers of the base modulus does the matrix's cycle length stay the same?
				stableLiftLevel += 1;
				
				//Keeps track of the cycle length from the previous modulus
				cycLen = omega(theCycle);
				
				//Increase modulus power
				multiply_BigIntT(temp, baseMod, temp2);
				copy_BigIntT(temp2, temp);
				
				//Calculate cycle length of next modulus
				big_floyd(UPDATEMATRIX, I, temp, &theCycle);
			}
			while ((omega(theCycle) == cycLen) && (stableLiftLevel < modPower));
			
			//I think if stableLiftLevel == modPower, then the matrix we're dealing with
			// is already a stable lift.
			if (stableLiftLevel == modPower)
			{
				printbm(UPDATEMATRIX);
				printf("is already a stable lift modulo ");
				printi(baseMod);
				printf("^%d. No computation needed.\n", modPower+1);
				goto FOUNDSTABLELIFT;
			}
			
			//Now, stableLiftLevel should hold the correct level
			if (debugMode)
				printf("stableLiftLevel: %d\n", stableLiftLevel);
			
			negOne = empty_BigIntT(1);
			
			//Construct moduli values
			bigMod     = new_BigIntT(numArr, 1);
			deltaMod   = new_BigIntT(numArr, 1);
			biggerMod  = empty_BigIntT(1);
			liftCompMod = new_BigIntT(numArr, 1);
			for (int i = 0; i < modPower; i += 1)
			{
				multiply_BigIntT(bigMod, baseMod, temp);
				copy_BigIntT(temp, bigMod);
				
				if (i == stableLiftLevel)
					copy_BigIntT(temp, liftCompMod);
				
				if (i == modPower - stableLiftLevel - 1)
					copy_BigIntT(temp, deltaMod);
			}
			copy_BigIntT(bigMod, biggerMod);
			multiply_BigIntT(biggerMod, baseMod, temp);
			copy_BigIntT(temp, biggerMod);
			
			subtract_BigIntT(biggerMod, one, negOne);
			
			//Get cycle length of matrix mod bigMod
			big_floyd(UPDATEMATRIX, I, bigMod, &theCycle);
			cycLen = omega(theCycle);
			
			//Extract matrix entries into a 2D array so we can manipulate them
			matrixElements = new_BigIntT_array(big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
			for (int row = 0; row < big_rows(UPDATEMATRIX); row += 1)
				for (int col = 0; col < big_rows(UPDATEMATRIX); col += 1)
					copy_BigIntT(big_element(UPDATEMATRIX, row, col), matrixElements[row][col]);
				
			//Get ready to find matrix deltas
			basePDigitMatrix = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
			matrixDeltaRows  = new_BigIntT_array(big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX),
			                                     big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX));
			tempMat = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
			tempMat2 = new_BigIntMatrixT(big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX), 
			                             big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX));
			
			deltaCombo            = empty_BigPolyT();
			deltaRowsToReduce     = new_BigIntMatrixT(1, big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX));
			deltaRowsToReduceTemp = new_BigIntMatrixT(1, big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX));
			matrixShiftWatcher    = new_BigIntMatrixT(1, big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX));
			
			matrixShiftWatcherRow = malloc(big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX)*sizeof(BigIntTP));
			for (int i = 0; i < big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX); i += 1)
				matrixShiftWatcherRow[i] = empty_BigIntT(1);
			
			//THIS WILL ONLY WORK FOR SMALL ENOUGH CYCLELENGTHS
			//FIX THIS LATER
			numArr[0] = cycLen;
			bigCycLen = new_BigIntT(numArr, 1);
			
			powbm(UPDATEMATRIX, basePDigitMatrix, bigCycLen, biggerMod);
			baseDeltaRow = malloc(big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX)*sizeof(BigIntTP));
			for (int r = 0; r < big_rows(basePDigitMatrix); r += 1)
				for (int c = 0; c < big_rows(basePDigitMatrix); c += 1)
				{
					baseDeltaRow[matrixIndexCounter] = empty_BigIntT(1);
					
					if (r == c)
						subtract_BigIntT(big_element(basePDigitMatrix, r, c), one, temp);
					else
						copy_BigIntT(big_element(basePDigitMatrix, r, c), temp);
					
					divide_BigIntT(temp, deltaMod, baseDeltaRow[matrixIndexCounter]);
					matrixIndexCounter += 1;
				}
			matrixIndexCounter = 0;
			
			//The row of matrix entry changes we're trying to eliminate
			// Once we find a way to eliminate it, we've found a way to make a stablelift
			if (debugMode)
			{
				printf("At iteration %d, baseDeltaRow represents how different our matrix iterate \
is from I. We want this row to be zero.\n", cycLen);
				printf("baseDeltaRow = [");
				for (int i = 0; i < big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX); i += 1)
				{
					if (i != 0)
						printf(", ");
					printi(baseDeltaRow[i]);
				}
				printf("]\n");
			}
			
			tempPDigitMatrix = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
			
			//Go through all component places in the matrix,
			// find how shifting the elements shifts the iterates
			for (int rowToShift = 0; rowToShift < big_rows(UPDATEMATRIX); rowToShift += 1)
			{
				for (int colToShift = 0; colToShift < big_rows(UPDATEMATRIX); colToShift += 1)
				{
					//Since the rows of our working matrix get shifted as we work,
					// matrixShiftWatcher keeps track of what components correspond to
					// which changes
					if (debugMode)
					{
						printf("~~~Beginning of loop~~~\n");
						/* printf("matrixShiftWatcher keeps track of the row-echelon moves. \
Each row says what each row in deltaRowsToReduce represents in terms of raw component changes.\n"); */
						printf("matrixShiftWatcher:\n");
						printbm(matrixShiftWatcher);
						/* printf("deltaRowsToReduce represents the row-echeloned changes that each matrix \
component causes to UPDATEMATRIX^%d mod p^{k+1}, written mod p^k.\n", cycLen); */
						printf("deltaRowsToReduce:\n");
						printbm(deltaRowsToReduce);
						printf("\n");
					}
					
					//Shift specific element
					add_BigIntT(baseMod, matrixElements[rowToShift][colToShift], temp);
					mod_BigIntT(temp, biggerMod, matrixElements[rowToShift][colToShift]);
					set_big_matrix(tempMat, matrixElements);
					
					//Revert matrixElements
					copy_BigIntT(big_element(UPDATEMATRIX, rowToShift, colToShift), matrixElements[rowToShift][colToShift]);
					
					//Calculate new deltas
					
					/* When trying to calculate deltas using an UPDATEMATRIX that's
					 * already a stable lift, things get weird...
					 * it seems like tempPDigitMatrix changes by increments of the
					 * baseMod, not the bigMod like usual...
					 * Maybe the amount we need to add for getting new shift matrices changes
					 * when we're working with a stable lift, or maybe we need to do the
					 * calculations under a different modulus
					 *
					 */
					powbm(tempMat, tempPDigitMatrix, bigCycLen, biggerMod);
					
					if (debugMode)
					{
						printf("basePDigitMatrix is simply UPDATEMATRIX^%d mod p^{k+1}.\n", cycLen);
						printf("basePDigitMatrix = \n");
						printbm(basePDigitMatrix);
						printf("tempPDigitMatrix keeps track of how much each component of our \
matrix changes UPDATEMATRIX^%d mod p^{k+1}.\n", cycLen);
						printf("tempPDigitMatrix = \n");
						printbm(tempPDigitMatrix);
					}
					
					deltaColCounter = 0;
					for (int deltaRow = 0; deltaRow < big_rows(UPDATEMATRIX); deltaRow += 1)
						for (int deltaCol = 0; deltaCol < big_rows(UPDATEMATRIX); deltaCol += 1)
						{
							multiply_BigIntT(big_element(basePDigitMatrix, deltaRow, deltaCol), negOne, temp);
							add_BigIntT(big_element(tempPDigitMatrix, deltaRow, deltaCol), temp, temp2);
							mod_BigIntT(temp2, biggerMod, temp);
							
							//Change the modulus used here depending on how tempPDigitMatrix comes out
							divide_BigIntT(temp, deltaMod, matrixDeltaRows[deltaRowCounter][deltaColCounter]);
							deltaColCounter += 1;
						}
					
					//So which modulus should I use? I think baseMod, but I'm not sure...
					if (deltaRowCounter != 0)
						resize_BigIntMatrixT(matrixShiftWatcher, big_rows(matrixShiftWatcher)+1, big_cols(matrixShiftWatcher));
					
					set_big_row(deltaRowsToReduce, matrixDeltaRows[deltaRowCounter], deltaRowCounter);
					
					for (int i = 0; i < big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX); i += 1)
					{
						if (i != deltaRowCounter)
							copy_BigIntT(zero, matrixShiftWatcherRow[i]);
						else
							copy_BigIntT(one, matrixShiftWatcherRow[i]);
					}
					set_big_row(matrixShiftWatcher, matrixShiftWatcherRow, deltaRowCounter);
					
					if (debugMode)
					{
						printf("After matrixShiftWatcher has been set\n");
						printf("matrixShiftWatcher:\n");
						printbm(matrixShiftWatcher);
						printf("deltaRowsToReduce:\n");
						printbm(deltaRowsToReduce);
						printf("\n");
					}
					
					//Prepare matrixShiftWatcher and deltaRowsToReduce for new rows
					if (debugMode)
					{
						printf("\n--BEFORE big_row_echelon()--\n");
						printf("deltaRowsToReduce = \n");
						printbm(deltaRowsToReduce);
						printf("\nliftCompMod = ");
						printi(liftCompMod);
						printf("\ndeltaRowsToReduceTemp = \n");
						printbm(deltaRowsToReduceTemp);
						printf("\nmatrixShiftWatcher = \n");
						printbm(matrixShiftWatcher);
					}
					
					big_row_echelon(deltaRowsToReduce, liftCompMod, deltaRowsToReduceTemp, matrixShiftWatcher);
					copy_BigIntMatrixT(deltaRowsToReduceTemp, deltaRowsToReduce);
					
					resize_BigIntMatrixT(deltaRowsToReduce, big_rows(deltaRowsToReduce)+1, big_cols(deltaRowsToReduce));
					resize_BigIntMatrixT(deltaRowsToReduceTemp, big_rows(deltaRowsToReduceTemp)+1, big_cols(deltaRowsToReduceTemp));
					
					if (debugMode)
					{
						printf("After row echelon and resize\n");
						printf("matrixShiftWatcher:\n");
						printbm(matrixShiftWatcher);
						printf("deltaRowsToReduce:\n");
						printbm(deltaRowsToReduce);
						printf("\n");
					}
					
					//Now, try to eliminate baseDeltaRow to see if we found a stable lift
					set_big_row(deltaRowsToReduce, baseDeltaRow, deltaRowCounter+1);
					
					//If we successfully eliminate the bottom row
					if (debugMode)
					{
						printf("After baseDeltaRow is set\n");
						printf("matrixShiftWatcher:\n");
						printbm(matrixShiftWatcher);
						printf("deltaRowsToReduce:\n");
						printbm(deltaRowsToReduce);
						printf("\n");
					}
					
					if (big_eliminate_bottom(deltaRowsToReduce, liftCompMod, deltaCombo))
					{
						//Used to keep track of all the different additions we need to do
						// to get a stablelift
						if (debugMode)
						{
							printf("deltaCombo = ");
							printp(deltaCombo);
							printf("\nmatrixDeltaRows = \n");
						}
						set_big_matrix(tempMat2, matrixDeltaRows);
						if (debugMode)
						{
							printbm(tempMat2);
							printf("\n");
						}
						
						//We'll reuse matrixShiftWatcherRow for convenience
						//It'll hold a representation of which components of the matrix need to be
						// modified to be a stable lift
						copy_BigIntT(zero, matrixShiftWatcherRow[deltaRowCounter]);
						
						//We need one of the bottom row plus whatever the polynomial specifies
						deltaComboCoeffs = extract_coefficients(deltaCombo);
						for (int dRow = 0; dRow < degree(deltaCombo)+1; dRow += 1)
						{
							if (compare_BigIntT(zero, deltaComboCoeffs[dRow]) == 0)
								continue;
							
							for (int dElement = 0; dElement < big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX); dElement += 1)
							{
								multiply_BigIntT(big_element(matrixShiftWatcher, dRow, dElement), deltaComboCoeffs[dRow], temp2);
								add_BigIntT(matrixShiftWatcherRow[dElement], temp2, temp);
								
								mod_BigIntT(temp, liftCompMod, matrixShiftWatcherRow[dElement]);
							}
						}
						
						if (debugMode)
						{
							printf("matrixShiftWatcherRow = [");
							for (int i = 0; i < big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX); i += 1)
							{
								if (i != 0)
									printf(", ");
								printi(matrixShiftWatcherRow[i]);
							}
							printf("]\n");
						}
						
						//Now, we create the matrix
						for (int mRow = 0; mRow < big_rows(UPDATEMATRIX); mRow += 1)
						{
							for (int mCol = 0; mCol < big_rows(UPDATEMATRIX); mCol += 1)
							{
								multiply_BigIntT(baseMod, matrixShiftWatcherRow[matrixIndexCounter], temp2);
								
								//Now, we check whether the found solution is valid. If not, inform the user
								mod_BigIntT(temp2, bigMod, temp3);
								if ((compare_BigIntT(matrixShiftWatcherRow[matrixIndexCounter], zero) != 0) && 
								    (compare_BigIntT(temp3, zero) != 0))
									isIncompleteSolution = TRUE;
								
								add_BigIntT(matrixElements[mRow][mCol], temp2, temp3);
								mod_BigIntT(temp3, biggerMod, matrixElements[mRow][mCol]);
								
								matrixIndexCounter += 1;
							}
						}
						
						if (isIncompleteSolution)
							printf("The matrix found is likely not a stable lift. Run the tool again on the given matrix.\n");
						
						set_big_matrix(tempMat, matrixElements);
						printf("Stable lift = \n");
						printbm(tempMat);
						
						//Break out of loop
						goto FOUNDSTABLELIFT;
					}
					
					if (debugMode)
					{
						printf("End of loop\n");
						printf("Now, deltaRowsToReduce should be reduced. If the bottom-most row \
is not zero, then we haven't found a stable lift yet.\n");
						printf("matrixShiftWatcher:\n");
						printbm(matrixShiftWatcher);
						printf("deltaRowsToReduce:\n");
						printbm(deltaRowsToReduce);
						printf("\n");
					}

					deltaRowCounter += 1;
				}
			}
			
			printf("No matrix found. A stable lift likely doesn't exist for the given matrix.\n");
			
			FOUNDSTABLELIFT:
			
			temp        = free_BigIntT(temp);
			temp2       = free_BigIntT(temp2);
			temp3       = free_BigIntT(temp3);
			temp4       = free_BigIntT(temp4);
			negOne      = free_BigIntT(negOne);
			bigMod      = free_BigIntT(bigMod);
			baseMod     = free_BigIntT(baseMod);
			deltaMod    = free_BigIntT(deltaMod);
			biggerMod   = free_BigIntT(biggerMod);
			liftCompMod = free_BigIntT(liftCompMod);
			
			theCycle = free_CycleInfoT(theCycle);
			I = free_BigIntMatrixT(I);
			
			matrixElements  = free_BigIntT_array(matrixElements, big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
			matrixDeltaRows = free_BigIntT_array(matrixDeltaRows, big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX),
			                                                      big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX));
			
			bigCycLen = free_BigIntT(bigCycLen);
			basePDigitMatrix = free_BigIntMatrixT(basePDigitMatrix);
			tempPDigitMatrix = free_BigIntMatrixT(tempPDigitMatrix);
			
			tempMat            = free_BigIntMatrixT(tempMat);
			tempMat2           = free_BigIntMatrixT(tempMat2);
			matrixShiftWatcher = free_BigIntMatrixT(matrixShiftWatcher);
			if (matrixShiftWatcherRow != NULL)
			{
				for (int i = 0; i < big_rows(UPDATEMATRIX)*big_rows(UPDATEMATRIX); i += 1)
				{
					matrixShiftWatcherRow[i] = free_BigIntT(matrixShiftWatcherRow[i]);
					baseDeltaRow[i]          = free_BigIntT(baseDeltaRow[i]);
				}
				FREE(matrixShiftWatcherRow);
				FREE(baseDeltaRow);
			}
			
			deltaRowsToReduce     = free_BigIntMatrixT(deltaRowsToReduce);
			deltaRowsToReduceTemp = free_BigIntMatrixT(deltaRowsToReduceTemp);
			
			deltaCombo = free_BigPolyT(deltaCombo);
		}
		
		
		//Half-test tool to find all stable lifts of a given matrix under the given modulus
		//This is by no means coded to be optimised. In fact, this tool will likely be
		// horrendously slow.
		else if (!strcmp(argv[1], "allstablelifts"))
		{
			//allstablelifts baseMod modPower
			//baseMod^modPower is the modulus of the thing we want to lift,
			// not the lift modulus itself (which is baseMod^{modPower+1}).
			if ((oargv[2][0] == '\0') || (oargv[3][0] == '\0'))
			{
				printf("Usage: " ANSI_COLOR_YELLOW "allstablelifts " ANSI_COLOR_CYAN "baseMod modPower" ANSI_COLOR_RESET "\n");
				goto FREE_VARIABLES;
			}
			
			BigIntTP baseMod, liftMultiple, bigMod, temp;
			int modPower;
			
			CycleInfoTP OGcycle = NULL;
			CycleInfoTP liftCycle = NULL;
			
			BigIntMatrixTP I;
			
			int numOfStableLifts = 0;
			BigIntTP** tempStableLiftElements;
			BigIntMatrixTP tempStableLift, tempMatrix;
			BigIntMatrixTP* listOfStableLifts = NULL;
			
			if (big_rows(UPDATEMATRIX) != big_cols(UPDATEMATRIX))
			{
				fprintf(stderr, "Given update matrix isn't square.\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			SET_BIG_NUM(oargv[2], baseMod, "Unable to read base mod from command line.");
			
			modPower = (int)strtol(oargv[3], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read modPower from command line.\n");
				baseMod = free_BigIntT(baseMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			temp = empty_BigIntT(1);
			bigMod = empty_BigIntT(1);
			liftMultiple = empty_BigIntT(1);
			copy_BigIntT(one, liftMultiple);
			
			//Calculate our powers of p we need
			for (int i = 0; i < modPower; i += 1)
			{
				multiply_BigIntT(liftMultiple, baseMod, temp);
				copy_BigIntT(temp, liftMultiple);
			}
			copy_BigIntT(liftMultiple, bigMod);
			multiply_BigIntT(bigMod, baseMod, temp);
			copy_BigIntT(temp, bigMod);
			
			//Now, get our cycle info for the thing we're lifting
			I = identity_BigIntMatrixT(big_rows(UPDATEMATRIX));
			big_floyd(UPDATEMATRIX, I, liftMultiple, &OGcycle);
			
			tempStableLiftElements = malloc(big_rows(UPDATEMATRIX)*sizeof(BigIntTP*));
			for (int r = 0; r < big_rows(UPDATEMATRIX); r += 1)
			{
				tempStableLiftElements[r] = malloc(big_rows(UPDATEMATRIX)*sizeof(BigIntTP));
				for (int c = 0; c < big_rows(UPDATEMATRIX); c += 1)
					tempStableLiftElements[r][c] = empty_BigIntT(1);
			}
		
			tempStableLift = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
			tempMatrix = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
			
			printf("Update =\n");
			printbm(UPDATEMATRIX);
			printf("\nStable lifts of Update modulo ");
			printi(bigMod);
			printf(":\n");
			
			//Loop over all possible lifts of UPDATEMATRIX
			do
			{
				//Make lift
				set_big_matrix(tempMatrix, tempStableLiftElements);
				big_mat_add(tempMatrix, UPDATEMATRIX, tempStableLift);
				
				//Check if cycle length equals omega(OGcycle)
				big_floyd(tempStableLift, I, bigMod, &liftCycle);
				if (omega(liftCycle) == omega(OGcycle))
				{
					//Add lift to list
					numOfStableLifts += 1;
					listOfStableLifts = realloc(listOfStableLifts, numOfStableLifts*sizeof(BigIntMatrixTP));
					listOfStableLifts[numOfStableLifts-1] = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
					copy_BigIntMatrixT(tempStableLift, listOfStableLifts[numOfStableLifts-1]);
					
					printbm(tempStableLift);
					printf("\n");
				}
			}
			while (!step_BigIntT_array(tempStableLiftElements, 
			                           big_rows(UPDATEMATRIX), 
																 big_rows(UPDATEMATRIX), 
																 liftMultiple, 
																 bigMod));

			for (int i = 0; i < numOfStableLifts; i += 1)
				listOfStableLifts[i] = free_BigIntMatrixT(listOfStableLifts[i]);
			FREE(listOfStableLifts);
			
			for (int r = 0; r < big_rows(UPDATEMATRIX); r += 1)
			{
				for (int c = 0; c < big_rows(UPDATEMATRIX); c += 1)
					tempStableLiftElements[r][c] = free_BigIntT(tempStableLiftElements[r][c]);
				FREE(tempStableLiftElements[r]);
			}
			FREE(tempStableLiftElements);
			
			temp = free_BigIntT(temp);
			bigMod = free_BigIntT(bigMod);
			baseMod = free_BigIntT(baseMod);
			liftMultiple = free_BigIntT(liftMultiple);
			
			OGcycle = free_CycleInfoT(OGcycle);
			liftCycle = free_CycleInfoT(liftCycle);
			
			I = free_BigIntMatrixT(I);
			tempMatrix = free_BigIntMatrixT(tempMatrix);
			tempStableLift = free_BigIntMatrixT(tempStableLift);
		}
		
		
		//Half-test tool to go through ALL matrices under a given modulus and 
		// count how many stable lifts each one has.
		else if (!strcmp(argv[1], "allstablelifts2"))
		{
			//allstablelifts2 size baseMod modPower [isPrime]
			//baseMod^modPower is the modulus of the thing we want to lift,
			// not the lift modulus itself (which is baseMod^{modPower+1}).
			
			PROHIBIT_UNIXFLAGS
			
			if (argc < 5)
			{
				printf("allstablelifts2 size baseMod modPower\n");
				goto FREE_VARIABLES;
			}
			
			//If the user entered a prime modulus, compute minimal
			// polynomials for each matrix we check.
			bool isPrime = FALSE;
			BigPolyTP* minPoly;
			char* typeOfMinPoly;
			int numOfLambdas;
			
			//This keeps track of how many stable lifts a matrix should have,
			// given its minimal polynomial.
			typedef struct mpti
			{
				char* type;
				int count;
			}
			MinPolyTypeIndexT;
			
			int numOfMinPolyTypes = 7;
			int numOfOutliers = 0; //How many matrices don't follow the established patterns.
			MinPolyTypeIndexT stableLiftCounts[numOfMinPolyTypes];
			stableLiftCounts[0].type = "Irreducible Quadratic";
			stableLiftCounts[1].type = "Single Linear";
			stableLiftCounts[2].type = "Single Nilpotent";
			stableLiftCounts[3].type = "Multiple Nilpotent";
			stableLiftCounts[4].type = "Single Linear, Single Lambda";
			stableLiftCounts[5].type = "Multiple Root";
			stableLiftCounts[6].type = "Distinct Roots";
			for (int i = 0; i < numOfMinPolyTypes; i += 1)
				stableLiftCounts[i].count = -1;
			
			BigIntTP baseMod, liftMultiple, bigMod;
			BigIntTP temp, temp2;
			int modPower, size;
			
			CycleInfoTP OGcycle = NULL;
			CycleInfoTP liftCycle = NULL;
			
			BigIntMatrixTP I;
			
			BigIntMatrixTP currentMatrix;
			BigIntTP** currentMatrixElements;
			
			int numOfStableLifts = 0;
			BigIntTP** tempStableLiftElements;
			BigIntMatrixTP tempStableLift, tempMatrix;
			BigIntMatrixTP* listOfStableLifts = NULL;
			
			size = (int)strtol(argv[2], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read matrix size from command line.\n");
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			SET_BIG_NUM(argv[3], baseMod, "Unable to read base mod from command line.");
			
			modPower = (int)strtol(argv[4], &tempStr, 10);
			if (tempStr[0] != '\0')
			{
				fprintf(stderr, "Unable to read modPower from command line.\n");
				baseMod = free_BigIntT(baseMod);
				returnvalue = EXIT_FAILURE;
				goto FREE_VARIABLES;
			}
			
			if (argc > 5)
				if (!strcmp(argv[5], "TRUE"))
					isPrime = TRUE;
			
			temp = empty_BigIntT(1);
			temp2 = empty_BigIntT(1);
			bigMod = empty_BigIntT(1);
			liftMultiple = empty_BigIntT(1);
			copy_BigIntT(one, liftMultiple);
			
			//Calculate our powers of p we need
			for (int i = 0; i < modPower; i += 1)
			{
				multiply_BigIntT(liftMultiple, baseMod, temp);
				copy_BigIntT(temp, liftMultiple);
			}
			copy_BigIntT(liftMultiple, bigMod);
			multiply_BigIntT(bigMod, baseMod, temp);
			copy_BigIntT(temp, bigMod);
			
			//Now, get our cycle info for the thing we're lifting
			I = identity_BigIntMatrixT(size);
			
			tempStableLiftElements = malloc(size*sizeof(BigIntTP*));
			currentMatrixElements = malloc(size*sizeof(BigIntTP*));
			for (int r = 0; r < size; r += 1)
			{
				tempStableLiftElements[r] = malloc(size*sizeof(BigIntTP));
				currentMatrixElements[r] = malloc(size*sizeof(BigIntTP));
				for (int c = 0; c < size; c += 1)
				{
					tempStableLiftElements[r][c] = empty_BigIntT(1);
					currentMatrixElements[r][c] = empty_BigIntT(1);
				}
			}
				
			currentMatrix = new_BigIntMatrixT(size, size);
			tempStableLift = new_BigIntMatrixT(size, size);
			tempMatrix = new_BigIntMatrixT(size, size);
			
			//Loop over all matrices in our matrix space
			do
			{
				set_big_matrix(currentMatrix, currentMatrixElements);
				big_floyd(currentMatrix, I, liftMultiple, &OGcycle);
				
				//Loop over all possible lifts of currentMatrix
				do
				{
					//Make lift
					set_big_matrix(tempMatrix, tempStableLiftElements);
					big_mat_add(tempMatrix, currentMatrix, tempStableLift);
					
					//Check if cycle length equals omega(OGcycle)
					big_floyd(tempStableLift, I, bigMod, &liftCycle);
					if (omega(liftCycle) == omega(OGcycle))
					{
						//Add lift to list
						numOfStableLifts += 1;
						listOfStableLifts = realloc(listOfStableLifts, numOfStableLifts*sizeof(BigIntMatrixTP));
						listOfStableLifts[numOfStableLifts-1] = new_BigIntMatrixT(size, size);
						copy_BigIntMatrixT(tempStableLift, listOfStableLifts[numOfStableLifts-1]);
					}
				}
				while (!step_BigIntT_array(tempStableLiftElements, size, size, liftMultiple, bigMod));

				//printbm(currentMatrix);
				if (isPrime) //Compute minimal polynomial
				{
					minPoly = min_poly(currentMatrix, NULL, liftMultiple, NULL);
					//printf("with minimal polynomial ");
					//old_printpf(minPoly);
					//printf(" ");
					
					//Now, let's try and classify the minimal polynomial (for 2x2 matrices).
					if (compare_BigIntT(constant(minPoly[0]), one) == 0)
					{
						if (degree(minPoly[1]) == 2)
							typeOfMinPoly = "Irreducible Quadratic";
						else if (compare_BigIntT(constant(minPoly[1]), zero) != 0)
							typeOfMinPoly = "Single Linear";
						else
							typeOfMinPoly = "Single Nilpotent";
					}
					
					else
					{
						clear_BigIntT(temp);
						numOfLambdas = 0;
						for (int i = 1; compare_BigIntT(temp, constant(minPoly[0])) < 0; i += 1)
						{
							if (compare_BigIntT(constant(minPoly[i]), zero) == 0)
								numOfLambdas += 1;
							
							add_BigIntT(temp, one, temp2);
							copy_BigIntT(temp2, temp);
						}
						
						if (numOfLambdas == 2)
							typeOfMinPoly = "Multiple Nilpotent";
						else if (numOfLambdas == 1)
							typeOfMinPoly = "Single Linear, Single Lambda";
						else
						{
							if (compare_BigIntT(constant(minPoly[1]), constant(minPoly[2])) == 0)
								typeOfMinPoly = "Multiple Root";
							else
								typeOfMinPoly = "Distinct Roots";
						}
					}
					
					//With the minimal polynomial type determined, let's check and see
					// if its number of stable lifts matches with its type of minimal
					// polynomial.
					for (int i = 0; i < numOfMinPolyTypes; i += 1)
						if (!strcmp(typeOfMinPoly, stableLiftCounts[i].type))
						{
							if (stableLiftCounts[i].count == -1)
								stableLiftCounts[i].count = numOfStableLifts;
							else if (stableLiftCounts[i].count != numOfStableLifts)
							{
								numOfOutliers += 1;
								//printf("does not follow the established pattern for its minimal polynomial. It ");
							}
							
							break;
						}
				}
				/*
				if (numOfStableLifts != 1)
					printf("has %d stable lifts modulo ", numOfStableLifts);
				else
					printf("has %d stable lift modulo ", numOfStableLifts);
				printi(bigMod);
				printf(".\n\n");
				*/

				//Reset stable lift counts for next matrix
				for (int i = 0; i < numOfStableLifts; i += 1)
					listOfStableLifts[i] = free_BigIntMatrixT(listOfStableLifts[i]);
				FREE(listOfStableLifts);
				numOfStableLifts = 0;
			}
			while (!increment_BigIntT_array(currentMatrixElements, size, size, one, liftMultiple));
			
			if (isPrime)
			{
				printf("Minimal polynomial type counts:\n");
				for (int i = 0; i < numOfMinPolyTypes; i += 1)
					printf("%s : %d\n", stableLiftCounts[i].type, stableLiftCounts[i].count);
				
				printf("\nNumber of matrix outliers: %d\n", numOfOutliers);
			}
			
			for (int r = 0; r < size; r += 1)
			{
				for (int c = 0; c < size; c += 1)
				{
					currentMatrixElements[r][c] = free_BigIntT(currentMatrixElements[r][c]);
					tempStableLiftElements[r][c] = free_BigIntT(tempStableLiftElements[r][c]);
				}
				FREE(tempStableLiftElements[r]);
				FREE(currentMatrixElements[r]);
			}
			FREE(tempStableLiftElements);
			FREE(currentMatrixElements);
			
			temp = free_BigIntT(temp);
			temp2 = free_BigIntT(temp2);
			bigMod = free_BigIntT(bigMod);
			baseMod = free_BigIntT(baseMod);
			liftMultiple = free_BigIntT(liftMultiple);
			
			OGcycle = free_CycleInfoT(OGcycle);
			liftCycle = free_CycleInfoT(liftCycle);
			
			I = free_BigIntMatrixT(I);
			tempMatrix = free_BigIntMatrixT(tempMatrix);
			currentMatrix = free_BigIntMatrixT(currentMatrix);
			tempStableLift = free_BigIntMatrixT(tempStableLift);
		}
		
		
		//If we want to get a sense for what the cyclespace of a particular vector looks like
		//This might become an extension of vectPolys, honestly
		else if (!strcmp(argv[1], "2023"))
		{
			PROHIBIT_UNIXFLAGS
			
			if (big_rows(UPDATEMATRIX) != big_cols(UPDATEMATRIX))
			{
				fprintf(stderr, "Given update matrix isn't square.\n");
				returnvalue = EXIT_SUCCESS;
				goto FREE_VARIABLES;
			}
			
			//Number of row vectors to have in matrixToReduce
			const int NUMOFCYCLESPACEGENERATORS = big_rows(INITIALMATRIX)+1;
			
			BigPolyTP coolPoly = NULL;
			
			BigIntTP bigMod;
			SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from config file.");
			
			//This'll be used to get a sense for what the cyclespace of our initial matrix looks like
			BigIntTP**     matrixElemsToReduce = malloc(NUMOFCYCLESPACEGENERATORS*sizeof(BigIntTP*));
			BigIntMatrixTP matrixToReduce;
			BigIntMatrixTP reducedMatrix;
			
			BigIntMatrixTP iteratedVect = new_BigIntMatrixT(big_rows(INITIALMATRIX), 1);
			BigIntMatrixTP tempVect     = new_BigIntMatrixT(big_rows(INITIALMATRIX), 1);
			
			
			for (int i = 0; i < NUMOFCYCLESPACEGENERATORS; i += 1)
			{
				matrixElemsToReduce[i] = malloc(big_rows(INITIALMATRIX)*sizeof(BigIntTP));
				for (int j = 0; j < big_rows(INITIALMATRIX); j += 1)
					matrixElemsToReduce[i][j] = empty_BigIntT(1);
			}
			
			copy_BigIntMatrixT(INITIALMATRIX, iteratedVect);
			
			//Iterate vector
			for (int iteration = 0; iteration < NUMOFCYCLESPACEGENERATORS; iteration += 1)
			{
				//Put iterates into a matrix
				for (int entry = 0; entry < big_rows(INITIALMATRIX); entry += 1)
					copy_BigIntT(big_element(iteratedVect, entry, 0), matrixElemsToReduce[iteration][entry]);
				
				//Iterate vector
				big_mat_mul(UPDATEMATRIX, iteratedVect, tempVect);
				modbm(tempVect, bigMod);
				copy_BigIntMatrixT(tempVect, iteratedVect);
			}
			
			//Reduce matrix
			//matrixToReduce = new_BigIntMatrixT(NUMOFCYCLESPACEGENERATORS, big_rows(UPDATEMATRIX));
			matrixToReduce = new_BigIntMatrixT(big_rows(UPDATEMATRIX), big_rows(UPDATEMATRIX));
			reducedMatrix  = new_BigIntMatrixT(NUMOFCYCLESPACEGENERATORS, big_rows(UPDATEMATRIX));
			//set_big_matrix(matrixToReduce, matrixElemsToReduce);
			copy_BigIntMatrixT(UPDATEMATRIX, matrixToReduce);
			
			printf("matrixToReduce:\n");
			printbm(matrixToReduce);
			
			//big_row_echelon(matrixToReduce, bigMod, reducedMatrix, NULL);
			coolPoly = empty_BigPolyT();
			big_eliminate_bottom(matrixToReduce, bigMod, coolPoly);
			
			printf("reducedMatrix after elimination:\n");
			printbm(matrixToReduce);
			printf("coolPoly: ");
			printp(coolPoly);
			printf("\n");
			
			for (int i = 0; i < NUMOFCYCLESPACEGENERATORS; i += 1)
			{
				for (int j = 0; j < big_rows(INITIALMATRIX); j += 1)
					matrixElemsToReduce[i][j] = free_BigIntT(matrixElemsToReduce[i][j]);
				FREE(matrixElemsToReduce[i]);
			}
			FREE(matrixElemsToReduce);
			
			matrixToReduce = free_BigIntMatrixT(matrixToReduce);
			reducedMatrix  = free_BigIntMatrixT(reducedMatrix);
			iteratedVect   = free_BigIntMatrixT(iteratedVect);
			tempVect       = free_BigIntMatrixT(tempVect);
			
			bigMod = free_BigIntT(bigMod);
			
			coolPoly = free_BigPolyT(coolPoly);
		}
		
		
		else if (! strcmp(argv[1], "multidebug"))
		{
			//2 + 2λ^2 + 4λ^3 + 1λ^4
			/*
			Extension definitions:
			2 + 2(t)^2 + 4(t)^3 + 1(t)^4 = 0
			1(t)^3
			Extension definitions:
			2 + 2(t)^2 + 4(t)^3 + 1(t)^4 = 0
			1(t)^3
			Extension definitions:
			2 + 2(t)^2 + 4(t)^3 + 1(t)^4 = 0
			1(t)^3
			*/
			
			/*
			printf("a : \n");
			printmve(a);
			printf("\na's BigIntDirectorT: ");
			display_BigIntDirectorT(a->coeffs, NULL, 0);
			printf("\n");
			exit(0);
			*/
			
			/*					
					polyCoeffs = extract_coefficients(coprimeFactors[f]);
					add_extension(oneExt, polyCoeffs, degree(coprimeFactors[f])+1, "t");
					add_extension(tempExt, polyCoeffs, degree(coprimeFactors[f])+1, "t");
					add_extension(factorExt, polyCoeffs, degree(coprimeFactors[f])+1, "t");
					add_extension(iterationExt, polyCoeffs, degree(coprimeFactors[f])+1, "t");
					
					for (int i = 0; i < degree(coprimeFactors[f])+1; i += 1)
						polyCoeffs[i] = free_BigIntT(polyCoeffs[i]);
					FREE(polyCoeffs);
					
					//Setting the values of our extensions to prepare for iteration
					increment_MultiVarExtT(oneExt);
					increment_MultiVarExtT(factorExt);
					set_MultiVarExtT_coefficient(iterationExt, tArray, one);
			*/
			
			BigIntTP bigMod;
			SET_BIG_NUM(bigintmodstring, bigMod, "Unable to read modulus from .config file.");
			
			BigPolyTP coolPoly = read_BigPolyT("data/test.polynomial");
			BigFactorsTP coolFactors = factor_BigPolyT(coolPoly, bigMod);
			BigPolyTP* coolFactorsSequel = extract_factors(coolFactors);
			
			BigIntTP** polyCoeffs = malloc(count_factors(coolFactors)*sizeof(BigIntTP*));
			for (int i = 0; i < count_factors(coolFactors); i += 1)
				polyCoeffs[i] = extract_coefficients(coolFactorsSequel[i]);
			
			int coeffsCount = 3;
			char* extNames[3] = {"α", "β", "γ"};
			int someCoeffsLocs[3][3] = {{1, 1, 1},  //αβγ
				                          {0, 2, 3},  //(β^2)(γ^3)
															    {0, 0, 4}}; //γ^4
			
			//Initialising our MultiVarExtT
			MultiVarExtTP a = new_MultiVarExtT(3);
			set_MultiVarExtT_mod(a, bigMod);
			
			for (int i = 0; i < count_factors(coolFactors); i += 1)
				add_extension(a, polyCoeffs[i], degree(coolFactorsSequel[i])+1, extNames[i]);
			
			for (int i = 0; i < coeffsCount; i += 1)
				set_MultiVarExtT_coefficient(a, someCoeffsLocs[i], one);
			
			printf("Modulus: ");
			printi(bigMod);
			printf("\ncoolPoly: ");
			printp(coolPoly);
			printf("\ncoolFactors: ");
			printpf(coolFactors);
			printf("\na: \n");
			printmve(a);
			
			printf("\nTesting the display_MultiVarExtT_internals() function...\n");
			display_MultiVarExtT_internals(a);
			
			
			for (int i = 0; i < count_factors(coolFactors); i += 1)
			{
				for (int j = 0; j < degree(coolFactorsSequel[i])+1; j += 1)
				{
					polyCoeffs[i][j] = free_BigIntT(polyCoeffs[i][j]);
				}
				FREE(polyCoeffs[i]);
				
			coolFactorsSequel[i] = free_BigPolyT(coolFactorsSequel[i]);
			}
			FREE(polyCoeffs);
			FREE(coolFactorsSequel);
			
			a = free_MultiVarExtT(a);
			bigMod = free_BigIntT(bigMod);
			coolPoly = free_BigPolyT(coolPoly);
			coolFactors = free_BigFactorsT(coolFactors);
		}
	}

	
	else
	{		
		printf(ANSI_COLOR_GREEN "LINCELLAUT by Zach Strong.\n" ANSI_COLOR_RESET);
		printf("Usage: lincellaut <tool> [options]\n\n");
		printf("Tools:\n");
		
		printf(" * " ANSI_COLOR_YELLOW "iterate " ANSI_COLOR_CYAN "[iterations]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "inverse " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "rowreduce " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "factor " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "order " ANSI_COLOR_RED "(UNFINISHED) " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "evalpoly " ANSI_COLOR_CYAN "[modulus] [multiplyByInitial]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "chara " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "allcharas " ANSI_COLOR_CYAN "coeffs..." ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "orbits " ANSI_COLOR_CYAN "[modulus] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "splitorbits " ANSI_COLOR_CYAN "[modulus] [fileoutput] " ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "orbitreps " ANSI_COLOR_CYAN "[modulus] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "branchreps " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "orbitspaces " ANSI_COLOR_RED "(UNFINISHED) " ANSI_COLOR_CYAN "[modulus] [minpolys] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "floyd " ANSI_COLOR_CYAN "[modulus]" ANSI_COLOR_RESET "\n");
		printf(" - " ANSI_COLOR_YELLOW "cycmatsearch " ANSI_COLOR_CYAN "resume size maxmod cycles..." ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "cycconvmat " ANSI_COLOR_CYAN "from to [mod]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "ccmzerosearch " ANSI_COLOR_CYAN "resume size [mod]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "vectprops " ANSI_COLOR_CYAN "baseMod modPower [resume] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "vectpolys " ANSI_COLOR_CYAN "baseMod modPower" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "matprops " ANSI_COLOR_CYAN "maxpower [modulus]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "fibmultsearch " ANSI_COLOR_CYAN "[bound]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "dynamics " ANSI_COLOR_CYAN "[maxPower] [modulus] [allConfigs] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "orbitmaps2 " ANSI_COLOR_CYAN "maxpower [modulus] [fileoutput] [belowBound] [aboveBound]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "oddminpolysearch " ANSI_COLOR_CYAN "maxpower size polysize [modulus] [resume] [fileoutput]" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "stablelift " ANSI_COLOR_CYAN "baseMod modPower" ANSI_COLOR_RESET "\n");
		printf(" * " ANSI_COLOR_YELLOW "allstablelifts " ANSI_COLOR_CYAN "baseMod modPower" ANSI_COLOR_RESET "\n");
		
		printf("Tools marked with a \"*\" in the list can use \"UNIX-like\" arguments.\n");
		printf("\nFor a more complete description of LINCELLAUT's usage, " \
		"refer to documentation/'LINCELLAUT CLI Usage.txt'.\n");
	}
	
	FREE_VARIABLES:
	FREE(updatefilepath);
	FREE(initialfilepath);
	FREE(resumefilepath);
	FREE(sentinelfilepath);
	FREE(polynomialfilepath);
	FREE(iterfilepath);
	FREE(bigintmodstring);
	FREE(resumemodstring);
	FREE(outputprefix);
	
	one  = free_BigIntT(one);
	zero = free_BigIntT(zero);
	
	UPDATEMATRIX  = free_BigIntMatrixT(UPDATEMATRIX);
	INITIALMATRIX = free_BigIntMatrixT(INITIALMATRIX);
	
	if (unixflags != NULL)
	{
		for (int r = 0; r < argc/2; r += 1)
		{
			for (int c = 0; c < 2; c += 1)
			{
				FREE(unixflags[r][c]);
			}
			FREE(unixflags[r]);
		}
		FREE(unixflags);
	}
	
	if (oargv != NULL)
	{
		if (toolindex == -1)
		{
			for (int i = 0; i < argc; i += 1)
			{
				FREE(oargv[i]);
			}
		}
		else
		{
			for (int i = 0; i <= MAXARGC[toolindex]; i += 1)
			{
				FREE(oargv[i]);
			}
		}
		FREE(oargv);
	}
	
	return returnvalue;
}
