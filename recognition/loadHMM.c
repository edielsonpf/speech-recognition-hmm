#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "estruturas.h"
#include "loadHMM.h"
#include "constantes.h"

// Le informacoes do arquivo .inf referente aos modelos HMM
int loadHMMInfo
(
  char *nome_HMM, // nome do arquivo com os modelos HMM das subunidades
  struct modelo_HMM *HMM// modelos HMM das subunidades foneticas
)
{
  char data[512],buffer[512]; // var aux p/ leitura do arquivo de informacao
  int i; // contador
  char NomeConfig[256]; // nome do arquivo com informacoes sobre o arquivo '.HMM'
  FILE *ptr; // file handlers
			
	// Inicializando algumas variaveis
	HMM->n_par = 0;
	HMM->n_estados = 3;
  HMM->salto_maximo = 2;

	
  // Mudando a extensao do arquivo de .hmm para .inf
	strcpy(NomeConfig,nome_HMM);
  char *p = strchr(NomeConfig,'.');
	strcpy(p+1,"inf");
  
		// Opening configuration file
  ptr = fopen(NomeConfig,"rt");
	
	if (ptr == NULL)
	{
		printf("Error opening file %s.\n",NomeConfig);
		return 1;
	}
	
  // Reading informations
	fgets(data, 500, ptr);
 
	// removing control characters (\r \n) frm the end of string
	while (data[strlen(data)-1]=='\r' || data[strlen(data)-1]=='\n')
		data[strlen(data)-1] = '\0';


  while (!feof(ptr))
  {
		// Numero de subunidades
		if (strstr(data,"Numero de subunidades=") != NULL)
	  	sscanf(data,"Numero de subunidades=%d",&HMM->n_fones);
  	// Numero de parametros
		if (strstr(data,"Numero de parametros=") != NULL)
		{
	  	sscanf(data,"Numero de parametros=%d",&HMM->n_par);
			// Allocating memory for HMM.parametros
			HMM->parametros = malloc(HMM->n_par*sizeof(char *));
			for (i=0;i<HMM->n_par;i++)
				HMM->parametros[i] = malloc(128*sizeof(char));
			// Allocating memory for HMM->delta
			HMM->delta = malloc(HMM->n_par*sizeof(int));
			// Allocating memory for HMM->n_gaussianas
			HMM->n_gaussianas = malloc(HMM->n_par*sizeof(int));
			// Allocating memory for HMM->ordem
			HMM->ordem = malloc(HMM->n_par*sizeof(int));
		}
		
		for (i=0;i<HMM->n_par;i++)
		{
			// Tipos de parametros
			sprintf(buffer,"Parametro %d=",i);
			if (strstr(data,buffer) != NULL)
			{
				strcat(buffer,"%s");
				sscanf(data,buffer,HMM->parametros[i]);
			}
			// Numero de janelas adjacentes para o calculo dos parametros delta
			sprintf(buffer,"Delta %d=",i);
			if (strstr(data,buffer) != NULL)
			{
				strcat(buffer,"%d");
				sscanf(data,buffer,&HMM->delta[i]);
			}
	    // Numero de gaussianas
			sprintf(buffer,"Gaussianas %d=",i);
			if (strstr(data,buffer) != NULL)
			{
				strcat(buffer,"%d");
				sscanf(data,buffer,&HMM->n_gaussianas[i]);
			}
	    // Ordem dos parametros
			sprintf(buffer,"Ordem %d=",i);
			if (strstr(data,buffer) != NULL)
			{
				strcat(buffer,"%d");
				sscanf(data,buffer,&HMM->ordem[i]);
			}
		}
		
		// Uso ou nao de PCA p/ reducao da dimensao do vetor de parametros
		if (strstr(data,"Usa PCA=") != NULL)
		{
			sscanf(data,"Usa PCA=%s",buffer);
			if(strcmp(buffer,"sim")==0)
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
		
  	// removing control characters (\r \n) frm the end of string
		while (data[strlen(data)-1]=='\r' || data[strlen(data)-1]=='\n')
			data[strlen(data)-1] = '\0';
	}	

	fclose(ptr);	
	return 0;	
}
//-------------------------------------------------------------------------

// Carrega modelos HMM das subunidades foneticas
void loadHMM
(
  char *nome_HMM, // nome do arquivo com os modelos HMM das subunidades
  struct modelo_HMM *HMM// modelos HMM das subunidades foneticas
)
{
  // Declaracao das variaveis locais
  FILE *arquivo; // handle para leitura de arquivo
  int register i,j,k,l,m; // contadores

	// Reading information from file .inf
	loadHMMInfo(nome_HMM,HMM);

  // Alocando memoria para os ponteiros A e B de HMM
  // (Matrizes de transicao e de emissao)

	// Matriz de transicao para os modelos HMM das subunidades
	HMM->A = malloc(HMM->n_fones*sizeof(double **));
	for (i=0;i<HMM->n_fones;i++)
		HMM->A[i] = malloc(HMM->n_estados*sizeof(double *));
	for (i=0;i<HMM->n_fones;i++)
		for (j=0;j<HMM->n_estados;j++)
			HMM->A[i][j] = malloc((HMM->salto_maximo+1)*sizeof(double));

	// Matriz de emissao para os modelos HMM das subunidades
	HMM->B = malloc(HMM->n_par*sizeof(struct mistura **));
	for (i=0;i<HMM->n_par;i++)
		HMM->B[i] = malloc(HMM->n_fones*sizeof(struct mistura *));
	for (i=0;i<HMM->n_par;i++)
		for (j=0;j<HMM->n_fones;j++)
			HMM->B[i][j] = malloc(HMM->n_estados*sizeof(struct mistura));
	for (i=0;i<HMM->n_par;i++)
		for (j=0;j<HMM->n_fones;j++)
			for (k=0;k<HMM->n_estados;k++)
			{
				HMM->B[i][j][k].c = malloc(HMM->n_gaussianas[i]*sizeof(double));
				HMM->B[i][j][k].m = malloc(HMM->n_gaussianas[i]*sizeof(double *));
				HMM->B[i][j][k].v = malloc(HMM->n_gaussianas[i]*sizeof(double *));
			}
	for (i=0;i<HMM->n_par;i++)
		for (j=0;j<HMM->n_fones;j++)
			for (k=0;k<HMM->n_estados;k++)
				for (l=0;l<HMM->n_gaussianas[i];l++)
				{
					HMM->B[i][j][k].m[l] = malloc(HMM->ordem[i]*sizeof(double));
					HMM->B[i][j][k].v[l] = malloc(HMM->ordem[i]*sizeof(double));
				}

  // Lendo arquivo .HMM
	arquivo = fopen(nome_HMM,"rb");			
  if (arquivo == NULL)
	{
		printf("Error opening file %s\n",nome_HMM);
		exit(1);
	}

  // Carregando matrizes de transicao
  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      for (k=0;k<(HMM->salto_maximo+1);k++)
				fread(&HMM->A[i][j][k], sizeof(double), 1, arquivo);

  // Carregando matrizes de emissao
  for (i=0;i<HMM->n_par;i++)
    for (j=0;j<HMM->n_fones;j++)
      for (k=0;k<HMM->n_estados;k++)
        for (l=0;l<HMM->n_gaussianas[i];l++)
        {
          fread(&HMM->B[i][j][k].c[l], sizeof(double), 1, arquivo);
          for (m=0;m<HMM->ordem[i];m++)
          {
            fread(&HMM->B[i][j][k].m[l][m], sizeof(double), 1, arquivo);
            fread(&HMM->B[i][j][k].v[l][m], sizeof(double), 1, arquivo);
          }
        }

  fclose(arquivo);

  // Calculando log10 das probabilidades de transicao
  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      for (k=0;k<(HMM->salto_maximo+1);k++)
        if (HMM->A[i][j][k] != 0.0)
          HMM->A[i][j][k] = log10(HMM->A[i][j][k]);
        else
          HMM->A[i][j][k] = -Inf;

}

//-------------------------------------------------------------------------
// Funcao que desaloca os ponteiros de HMM
void freeHMM
(
  struct modelo_HMM *HMM// modelos HMM das subunidades foneticas
)
{
  int register i,j,k,l; // contador

  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      free(HMM->A[i][j]);
  for (i=0;i<HMM->n_fones;i++)
    free(HMM->A[i]);
  free(HMM->A);
  HMM->A = NULL;

  for (i=0;i<HMM->n_par;i++)
    for (j=0;j<HMM->n_fones;j++)
      for (k=0;k<HMM->n_estados;k++)
        for (l=0;l<HMM->n_gaussianas[i];l++)
        {
          free(HMM->B[i][j][k].m[l]);
          free(HMM->B[i][j][k].v[l]);
        }
  for (i=0;i<HMM->n_par;i++)
    for (j=0;j<HMM->n_fones;j++)
      for (k=0;k<HMM->n_estados;k++)
      {
        free(HMM->B[i][j][k].c);
        free(HMM->B[i][j][k].m);
        free(HMM->B[i][j][k].v);
      }
  for (i=0;i<HMM->n_par;i++)
    for (j=0;j<HMM->n_fones;j++)
      free(HMM->B[i][j]);
  for (i=0;i<HMM->n_par;i++)
    free(HMM->B[i]);
  free(HMM->B);
  HMM->B = NULL;

  free(HMM->delta);
  HMM->delta = NULL;

  free(HMM->n_gaussianas);
  HMM->n_gaussianas = NULL;

  free(HMM->ordem);
  HMM->ordem = NULL;

  for (i=0;i<HMM->n_par;i++)
    free(HMM->parametros[i]);
  free(HMM->parametros);
  HMM->parametros = NULL;
}
