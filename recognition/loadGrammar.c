#include <stdlib.h>
#include <stdio.h>
#include "estruturas.h"
#include "loadGrammar.h"

// Funcao que carrega modelo de gramatica bigram de palavras, gerado pelo programa Gramatica.exe

void loadGrammar
(
  char *nome, // nome do arquivo com os dados da gramatica
  long *n_termos, // numero de pares na gramatica
  struct bigram **grammar1 // pares da gramatica
)
{
  FILE *arquivo; // handle p/ leitura de arquivo
  long i; // contador
  struct bigram *grammar; // armazena gramatica bigram
  grammar = *grammar1;

  // Abrindo e lendo arquivo com gramatica
	arquivo = fopen(nome,"rb");
	if(arquivo==NULL)
	{
		printf("Error opening file %s.\n",nome);
		exit(1);
	}

  // Lendo numero de pares da gramatica
	fread(n_termos,sizeof(long),1,arquivo);

  // Alocando memoria p/ gramatica
	grammar = malloc(*n_termos*sizeof(struct bigram));

  // Carregando gramatica a partir de arquivo
	for (i=0;i<*n_termos;i++)
  {
		fread(&grammar[i].primeira,sizeof(int),1,arquivo);
		fread(&grammar[i].segunda,sizeof(int),1,arquivo);
		fread(&grammar[i].prob,sizeof(float),1,arquivo);
  }
  
  fclose(arquivo);

  *grammar1 = grammar;
}
