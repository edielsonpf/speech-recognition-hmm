#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "estruturas.h"
#include "par.h"
#include "oneStep.h"
#include "herrmanNey.h"

void reconheceFrases
(
  struct configBusca searchConfig, // search algorithm configurations	
	char *utterancesList, // list of utterances to be recognized
	struct modelo_HMM HMM, // acoustic models
	struct bigram *grammar, // stores the bigram information
	long nTerms,	// number of entries in the grammar
	struct vocab vocabulary, // stores the vocabulary
	char *name_out
)
{	
	FILE *ptr; // handle for file
	char utteranceFile[512]; // name of file with the utterance to be recognized
	struct locucao x; // stores the acoustic parameters of the utterance to be recognized
	int i,j; // counters	
	char frase[512]; // recognized string
	time_t t_inic,t_final; // instantes inicial e final
	int decorrido;
	int n=0; // counts the  number of recognized utterances
	FILE *arq_out;
	
	ptr = fopen(utterancesList,"rt");
	if (ptr == NULL)
	{
		printf("Error opening file %s.\n",utterancesList);
		exit(1);
	}
	// Reading informations
	fgets(utteranceFile,512, ptr);
	
	// removing control characters (\r \n) frm the end of string
	while (utteranceFile[strlen(utteranceFile)-1]=='\r' || utteranceFile[strlen(utteranceFile)-1]=='\n')
			utteranceFile[strlen(utteranceFile)-1] = '\0';

    // Opening output file
    arq_out = fopen(name_out,"wt");
    if (arq_out == NULL)
    {
      printf("Error opening output file %s\n",name_out);	
      exit(1);
    }
	
  while (!feof(ptr))
  {
    printf("Processing file %s...\n",utteranceFile);
		
	// Iniciando contagem de tempo de processamento
    time(&t_inic);
	
    // Parametrizando sinal de voz (aloca ponteiro 'x')
    // CalcPar e definido em Par.h
    // bUsaPCA,NomePCA,OrdemRed definidos em HMM.h
	calcPar(utteranceFile,HMM,&x.n_quadros,&x.par);

    // Algoritmo de busca
    switch (searchConfig.searchAlgorithm)
    {
      case 0:
        puts("Level Building sucks!\n");
        break;
      case 1:
        oneStep(x,frase,searchConfig,HMM,grammar,nTerms,vocabulary);
        break;
      default:
        HerrmanNey(x,frase,searchConfig,HMM,grammar,nTerms,vocabulary);
    }

    // removing control characters (\r \n) frm the end of string
	while (frase[strlen(frase)-1]=='\r' || frase[strlen(frase)-1]=='\n')
		frase[strlen(frase)-1] = '\0';
		
    // Parando cronometragem
    time(&t_final);

    // Calculando tempo decorrido
    decorrido = (int) difftime(t_final,t_inic);
	printf("Tempo: %d segundos\n",decorrido);
    
    // Apresentando frase reconhecida
    printf("%s\n",frase);

    // Salvando resultados em arquivo
    for (i=0;i<strlen(frase);i++) // eliminando virgulas para sclite
      if(frase[i]==',')
        frase[i] = ' ';
    fprintf(arq_out,"%s (%d)\n",frase,n++);

    // Desalocando ponteiros
    for (i=0;i<HMM.n_par;i++)
      for (j=0;j<x.n_quadros;j++)
        free(x.par[i][j]);
    for (i=0;i<HMM.n_par;i++)
      free(x.par[i]);
    free(x.par);

		// Reading new entry from file
		if (!feof(ptr))
			fgets(utteranceFile, 500, ptr);
		
  	// removing control characters (\r \n) frm the end of string
		while (utteranceFile[strlen(utteranceFile)-1]=='\r' || utteranceFile[strlen(utteranceFile)-1]=='\n')
			utteranceFile[strlen(utteranceFile)-1] = '\0';
  }
	fclose(ptr);
	fclose(arq_out);
	
}
