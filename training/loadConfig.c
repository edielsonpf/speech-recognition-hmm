#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "estruturas.h"	
#include "loadConfig.h"

int loadConfig
(
  char *configFile,
  struct config *trainConfig,
  struct modelo_HMM *HMM
)
{
  // Local variables
  FILE *ptr; // file handlers
  char data[512],aux[512]; // aux variables for file reading
  int i; // counter

  // Initializing some variables (for valgrind to stop complaining)
  HMM->n_par = 1;	
  
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
  	// Initialization mode
    if (strstr(data,"Inicializacao=") != NULL)
	{
	  sscanf(data,"Inicializacao=%s",aux);
	  if (strcmp(aux,"Padrao") == 0)
        trainConfig->initializationMode = 0;
	  else if (strcmp(aux,"SegmentalKMeans") == 0)
	   	trainConfig->initializationMode = 1;
	  else if (strcmp(aux,"ModelosPreTreinados") == 0)
	   	trainConfig->initializationMode = 2;	
	}
	// HMM file
    if (strstr(data,"Arquivo de saida=") != NULL)
	  	sscanf(data,"Arquivo de saida=%s",trainConfig->HMMFile);
	// Phonetic subunits file
    if (strstr(data,"Nome subunidades=") != NULL)
	  	sscanf(data,"Nome subunidades=%s",trainConfig->subunitsFile);
	// File with transcription files listing
    if (strstr(data,"Arquivo transcricoes") != NULL)
	  	sscanf(data,"Arquivo transcricoes=%s",trainConfig->transcriptionsFiles);
	// File with training utterances files listing
    if (strstr(data,"Arquivo locucoes") != NULL)
	  	sscanf(data,"Arquivo locucoes=%s",trainConfig->utterancesFiles);
	// Number of training utterances
    if (strstr(data,"Numero de arquivos=") != NULL)
	  	sscanf(data,"Numero de arquivos=%d",&trainConfig->nUtterances);
	// Name of file with pre-trained models
    if (strstr(data,"Nome pre-treinado=") != NULL)
	  	sscanf(data,"Nome pre-treinado=%s",trainConfig->preTrained);
	// Number of states for each HMM
	if (strstr(data,"Numero de estados=") != NULL)
	  	sscanf(data,"Numero de estados=%d",&HMM->n_estados);
    // Maximun allowed skip beteen states 
	if (strstr(data,"Salto maximo=") != NULL)
	  	sscanf(data,"Salto maximo=%d",&HMM->salto_maximo);	  	  	
	// Number of acoustic parameters
    if (strstr(data,"Numero de parametros=") != NULL)
    {
	  	sscanf(data,"Numero de parametros=%d",&HMM->n_par);
	  	HMM->parametros = malloc(sizeof(char *)*HMM->n_par);
	  	for (i=0;i<HMM->n_par;i++)
	  		HMM->parametros[i] = malloc(sizeof(char)*20);
	  	HMM->delta = malloc(sizeof(int)*HMM->n_par);
	  	HMM->n_gaussianas = malloc(sizeof(int)*HMM->n_par);	
	  	HMM->ordem = malloc(sizeof(int)*HMM->n_par);
    }
	// Acoustic parameter information
	for (i=0;i<HMM->n_par;i++)
	{
		// Tipos de parametros
		sprintf(aux,"Parametro %d=",i);
		if (strstr(data,aux) != NULL)
		{
			strcat(aux,"%s");
			sscanf(data,aux,HMM->parametros[i]);
		}
		// Numero de janelas adjacentes para o calculo dos parametros delta
		sprintf(aux,"Delta %d=",i);
		if (strstr(data,aux) != NULL)
		{
			strcat(aux,"%d");
			sscanf(data,aux,&HMM->delta[i]);
		}
	    // Numero de gaussianas
		sprintf(aux,"Gaussianas %d=",i);
		if (strstr(data,aux) != NULL)
		{
			strcat(aux,"%d");
			sscanf(data,aux,&HMM->n_gaussianas[i]);
		}
	    // Ordem dos parametros
		sprintf(aux,"Ordem %d=",i);
		if (strstr(data,aux) != NULL)
		{
			strcat(aux,"%d");
			sscanf(data,aux,&HMM->ordem[i]);
		}
	}
	// Uso ou nao de PCA p/ reducao da dimensao do vetor de parametros
	if (strstr(data,"Usa PCA=") != NULL)
	{
		sscanf(data,"Usa PCA=%s",aux);
		if(strcmp(aux,"sim")==0)
			HMM->usaPCA = 1;
		else
			HMM->usaPCA = 0;
	}

	// Nome do arquivo com matriz de transformacao
	if (strstr(data,"Nome PCA=") != NULL)
		sscanf(data,"Nome PCA=%s",HMM->nomePCA);

	// Ordem do vetor de parametros reduzido
	if (strstr(data,"Ordem reduzida=") != NULL)
		sscanf(data,"Ordem reduzida=%d",&HMM->ordemRed);
	
	
	// Reading new entry from file
	if (!feof(ptr))
		fgets(data, 500, ptr);
		//fscanf(ptr,"%s",data);
		
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

// Function that shows the training configuration in the screen (just for debug)
void showConfig
(
  struct config trainConfig,
  struct modelo_HMM HMM
)
{
	int i; // counter
	
	printf("Initialization mode: %d\n",trainConfig.initializationMode);
	printf("Output file: %s\n",trainConfig.HMMFile);
	printf("Phonetic subunits file: %s\n",trainConfig.subunitsFile);
	printf("File with transcription files listing: %s\n",trainConfig.transcriptionsFiles);
	printf("File with training utterances listing: %s\n",trainConfig.utterancesFiles);
	printf("Number of training utterances: %d\n",trainConfig.nUtterances);
	printf("File with pre-trained HMMs: %s\n",trainConfig.preTrained);
	printf("States per HMM: %d\n",HMM.n_estados);
	printf("Maximun allowed skip between states: %d\n",HMM.salto_maximo);
	printf("Number of acoustic parameters: %d\n",HMM.n_par);
	for (i=0;i<HMM.n_par;i++)
	{
		printf("Parametro %d: %s\n",i,HMM.parametros[i]);
		printf("Delta: %d\n",HMM.delta[i]);
		printf("Gaussianas: %d\n",HMM.n_gaussianas[i]);
		printf("Ordem: %d\n",HMM.ordem[i]);
	}
	printf("PCA: ");
	if(HMM.usaPCA == 1)
	  	printf("sim\n");
	else
	  	printf("nao\n");
	if(HMM.usaPCA == 1)
	{
		printf("PCA file: %s\n",HMM.nomePCA);
		printf("Ordem reduzida: %d\n",HMM.ordemRed);
	}  	
}

