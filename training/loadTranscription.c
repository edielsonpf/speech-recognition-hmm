#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "mystrfun.h"

void loadTranscription
(
	char *filename,
	char **fones,
	int n_fonemas,
	int **transcription1,
	int *length1
)

{
	FILE *ptr; // file handler
	char transc[512];
	int achou;
	int n_caracteres;
	char **caracteres; // palavras da frase
	const int MaxCrtPlv = 10;
    int i,j; // counters
    int length;
    int *transcription;
    	
	// Reading transcription from file
	ptr=fopen(filename,"rt");
	fgets(transc,512,ptr);
	fclose(ptr);
	
	// removing control characters (\r \n) frm the end of string
    while (transc[strlen(transc)-1]=='\r' || transc[strlen(transc)-1]=='\n')
  	  transc[strlen(transc)-1] = '\0';
	
	// Length of string
	n_caracteres = strlen(transc);
	
	// Allocating memory
	caracteres = malloc(n_caracteres*sizeof(char*));
    for (i=0;i<n_caracteres;i++)
      caracteres[i] = malloc(MaxCrtPlv*sizeof(char));
	
	// Separating phones
    codpal(transc,caracteres,&length," ","/",10,300);

    transcription = malloc(sizeof(int)*length);

    // Atribuindo indice aos fonemas
    for (i=0;i < length;i++)
    {
      achou=0;
      j=0;
      while ((achou == 0) && (i < length) && (j < n_fonemas))
      {
        if ((strcmp(caracteres[i],fones[j]) == 0) && (j < n_fonemas))
        {
          achou = 1;
          transcription[i] = j;
        }
        j++;
      }
      if (achou == 0)
      {
        printf("Fone %s desconhecido.\n",caracteres[i]);
        exit(1);
      }
    }
  
    // Desalocando ponteiros
    for (i=0;i<n_caracteres;i++)
      free(caracteres[i]);
    free(caracteres);

	*length1 = length;
	*transcription1 = transcription;
}
