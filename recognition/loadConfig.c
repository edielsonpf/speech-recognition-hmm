/******************************************************************************
 * Function for config files reading
 * Developed by Carlos Alberto Ynoguti
 *
 * Usage: 
 *
 * 
 * 
 * The organization of the config.cfg is described in the file CONFIGURE.TXT
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "estruturas.h"	
#include "loadConfig.h"

int loadConfig
(
  char *configFile,
  //~ int *searchAlgorithm, // 0: Level Building, 1: One Step, 2: Herrmann Ney
  //~ int *useBeamSearch, // 0 yes, 1 no
	//~ float *beamSearchThreshold, // threshold for node pruning
	//~ int *numberOfLevels, // maximum number of words per utterance
	//~ int *useWordDurationModel, // use (0) or not (1) the word duration model
	//~ int *useGrammar, // use (0) or not (1) the language model
	struct configBusca *searchConfig,
	char *grammarFile, // file with the language model
	char *vocabFile, // file with the vocabulary
	char *HMMFile, // file with the HMM parameters
	int *nUtterances, // number of utterances to be recognized
	char *utterancesList // file with the list of utterances to be recognized
)
{
  // Local variables
  FILE *ptr; // file handlers
  char data[512],aux[512];
	
	// Allocating memory	
	//data = malloc(512*sizeof(char));
	//aux = malloc(512*sizeof(char));
  
	// Opening configuration file
  ptr = fopen(configFile,"rt");
	
	if (ptr == NULL)
	{
		printf("Error opening file %s.\n",configFile);
		return 1;
	}
	
  // Reading informations
	fgets(data, 500, ptr);
  
	// removing control characters (\r \n) frm the end of string
	while (data[strlen(data)-1]=='\r' || data[strlen(data)-1]=='\n')
		data[strlen(data)-1] = '\0';

  while (!feof(ptr))
  {
		// Search algorithm
    if (strstr(data,"Algoritmo=") != NULL)
		{
	  	sscanf(data,"Algoritmo=%s",aux);
			if (strcmp(aux,"LevelBuilding") == 0)
	    	searchConfig->searchAlgorithm = 0;
	  	else if (strcmp(aux,"OneStep") == 0)
	    	searchConfig->searchAlgorithm = 1;
	    else if (strcmp(aux,"HerrmannNey") == 0)
	    	searchConfig->searchAlgorithm = 2;	
		}
		// Beam Search
    if (strstr(data,"Beam search=") != NULL)
		{
	  	sscanf(data,"Beam search=%s",aux);
			if (strcmp(aux,"sim") == 0)
	    	searchConfig->useBeamSearch = 1;
	  	else 
  	  	searchConfig->useBeamSearch = 0;
		}
		// Beam Search threshold
    if (strstr(data,"Limiar de poda=") != NULL)
	  	sscanf(data,"Limiar de poda=%f",&searchConfig->beamSearchThreshold);
		// Number of levels
    if (strstr(data,"Numero de niveis=") != NULL)
	  	sscanf(data,"Numero de niveis=%d",&searchConfig->numberOfLevels);
		// Word duration model
    if (strstr(data,"Modelo de duracao=") != NULL)
		{
	  	sscanf(data,"Modelo de duracao=%s",aux);
			if (strcmp(aux,"sim") == 0)
	    	searchConfig->useWordDurationModel = 1;
	  	else 
  	  	searchConfig->useWordDurationModel = 0;
    }
		// Grammar
    if (strstr(data,"Usa gramatica=") != NULL)
		{
	  	sscanf(data,"Usa gramatica=%s",aux);
			if (strcmp(aux,"sim") == 0)
	    	searchConfig->useGrammar = 1;
	  	else 
  	  	searchConfig->useGrammar = 0;
		}
		// Grammar File
    if (strstr(data,"Gramatica=") != NULL)
	  	sscanf(data,"Gramatica=%s",grammarFile);
		// Vocabulary File
    if (strstr(data,"Vocabulario=") != NULL)
	  	sscanf(data,"Vocabulario=%s",vocabFile);
		// HMM File
    if (strstr(data,"Modelos HMM=") != NULL)
	  	sscanf(data,"Modelos HMM=%s",HMMFile);
		// Number of utterances to be recognized
    if (strstr(data,"Numero de locucoes=") != NULL)
	  	sscanf(data,"Numero de locucoes=%d",nUtterances);
		// File with the list of utterances to be recognized
    if (strstr(data,"Locucoes=") != NULL)
	  	sscanf(data,"Locucoes=%s",utterancesList);
		
		// Reading new entry from file
		if (!feof(ptr))
			fgets(data, 500, ptr);
		
  	// removing control characters (\r \n) frm the end of string
		while (data[strlen(data)-1]=='\r' || data[strlen(data)-1]=='\n')
			data[strlen(data)-1] = '\0';


  } // end of configuration file reading

  fclose(ptr);
	
	// Deallocating memory
	//free(data);
	//free(aux);
	
  return (0);
}
