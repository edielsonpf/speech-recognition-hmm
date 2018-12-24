#include <stdio.h>
#include <stdlib.h>

#include "estruturas.h"
#include "loadConfig.h"
#include "loadSubunits.h"
#include "HMMFiles.h"
#include "forback.h"

int main(int argc, char *argv[])
{
	// Local variables
	struct config trainConfig; // training configuration settings 
	struct modelo_HMM HMM;     // stores HMMs
	char **subunits;           // list of phonetic subunits      
	//int i;
	
	// Testing arguments
	if (argc != 2)
	{
		printf("Usage: treina config.txt\n");
		printf("where config.txt is a textfile with configurations settings (see manual for details).\n");
		return (1);
	}
	
	printf("Processing file %s\n",argv[1]);
	
	// Loading training configurations
	loadConfig(argv[1],&trainConfig,&HMM);
	
	// Showing training configurations
	//showConfig(trainConfig,HMM);
	
	// Loading list of phonetic subunits
	loadSubunits(trainConfig.subunitsFile,&HMM.n_fones,&subunits);
	/*
	printf("Number of phonetic subunits: %d\n",HMM.n_fones);
	printf("Phonetic subunits: \n");
	for(i=0;i<HMM.n_fones;i++)
		printf("%s\n",subunits[i]);
	*/
	
	// Training HMMs
	forwardBackward(subunits,trainConfig,&HMM);
	
	// Freeing memory
	freeSubunits(&subunits,HMM.n_fones);
    freeHMM(&HMM);	
	return EXIT_SUCCESS;
}
