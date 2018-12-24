#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "estruturas.h"
#include "constantes.h"
#include "HMMFiles.h"


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
        free(HMM->B[i][j][k].m);
        free(HMM->B[i][j][k].c);
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

// Funcao que guarda em arquivo modelos HMM treinados
//
// - numero de subunidades foneticas (int)
// - numero de parametros (int)
// - Para cada parametro: - numero de gaussianas na mistura (int)
//                        - ordem (int)
//                        - extensao dos arquivos de parametros (char [MAXEXT])
void SalvaHMM
(
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  char nome_saida[256], // nome do arquivo com os parâmetros HMM das sub-unidades
  int epoca // numero da epoca de treinamento
)
{
  // Declaracao das variaveis locais
  FILE *arquivo; // handle para arquivo de saida
  int register i,j,k,l,m; // contadores
  
  // Abrindo arquivo com os parametros pre-treinados
  arquivo = fopen(nome_saida,"wb");
  if (arquivo == NULL)
  {
  	printf("Error creating file %s. Disk full? Verify permissions.\n",nome_saida);
    exit(1);
  }
  // Salvando matrizes de transicao
  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      for (k=0;k<(HMM->salto_maximo+1);k++)
        fwrite(&HMM->A[i][j][k],sizeof(double),1,arquivo);

  // Salvando matrizes de emissao
  for (i=0;i<HMM->n_par;i++)
    for (j=0;j<HMM->n_fones;j++)
      for (k=0;k<HMM->n_estados;k++)
        for (l=0;l<HMM->n_gaussianas[i];l++)
        {
          fwrite(&HMM->B[i][j][k].c[l],sizeof(double),1,arquivo);
          for (m=0;m<HMM->ordem[i];m++)
          {
            fwrite(&HMM->B[i][j][k].m[l][m],sizeof(double),1,arquivo);
            fwrite(&HMM->B[i][j][k].v[l][m],sizeof(double),1,arquivo);
          }
        }
  fclose(arquivo);
  /*
  // Salvando dados de treinamento em arquivo texto
  strcpy(nome_saida,"epoca");
  strcat(nome_saida,IntToStr(epoca).c_str());
  strcat(nome_saida,".txt");
  if ((arquivo = FormArquivo->AbreArquivo(nome_saida,"wb")) == NULL)
    exit(1);
  else
  {
    for (i=0;i<HMM->n_par;i++)
    {
      fprintf(arquivo,"Parametro: %d\n",i);
      for (j=0;j<HMM->n_fones;j++)
      {
        fprintf(arquivo,"Fone: %d\n",j);
        for (k=0;k<HMM->n_estados;k++)
        {
          fprintf(arquivo,"Estado: %d\n",k);
          for (l=0;l<HMM->n_gaussianas[i];l++)
          {
            fprintf(arquivo,"Gaussiana: %d\n",l);
            fprintf(arquivo,"Coeficiente de ponderação: %f\n",FloatToStr(HMM->B[i][j][k].c[l]).c_str());
            fprintf(arquivo,"Média      Variância\n");
            for (m=0;m<HMM->ordem[i];m++)
              fprintf(arquivo,"%f %f\n",FloatToStr(HMM->B[i][j][k].m[l][m]).c_str(),FloatToStr(HMM->B[i][j][k].v[l][m]).c_str());
          }
        }
      }
    }
    fclose(arquivo);
  }
  */
}

/////////////////////////////////////////////////////////////////////////////////////
// Gera arquivo de informacao sobre o arquivo HMM sendo treinado
void GeraArquivoInformacao
(
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  char* nome_saida // nome do arquivo de saida
)
{
  int register i; // contador
  char NomeInf[500]; // nome do arquivo de informacoes sobre o arquivo HMM
  FILE *arquivo; // handle for output file
  
  // Gerando nome do arquivo de saida
  // Verificando arquivo de informacoes dos modelos HMM das subunidades
  strcpy(NomeInf,nome_saida);
  char *p = strchr(NomeInf,'.');
  strcpy(p+1,"inf");
  remove(NomeInf); // apagando versoes anteriores, se existirem

  // Opening output file
  arquivo = fopen(NomeInf,"wt");
  if(arquivo==NULL)
  {
  	printf("Erro ao criar arquivo %s. Disco cheio? Permissoes?\n",NomeInf);
  	exit(1);
  }
  
  // Dados das subunidades foneticas
  fprintf(arquivo,"[Subunidades]\n");
  
  // Numero de subunidades foneticas
  fprintf(arquivo,"Numero de subunidades=%d\n",HMM->n_fones);

  // Dados dos parametros acusticos
  fprintf(arquivo,"[Parametros]\n");
  
  // Numero de parametros
  fprintf(arquivo,"Numero de parametros=%d\n",HMM->n_par);

  for (i=0;i<HMM->n_par;i++)
  {
    // Parametros a serem utilizados
    fprintf(arquivo,"Parametro %d=%s\n",i,HMM->parametros[i]);

    // Numero de janelas adjacentes utilizadas no calculo dos parametros delta
    fprintf(arquivo,"Delta %d=%d\n",i,HMM->delta[i]);

    // Numero de gaussianas nas misturas
    fprintf(arquivo,"Gaussianas %d=%d\n",i,HMM->n_gaussianas[i]);

    // Ordem dos vetores de parametros
    fprintf(arquivo,"Ordem %d=%d\n",i,HMM->ordem[i]);
  }

  // Uso ou nao de PCA p/ reducao da dimensao do vetor de parametros
  if (HMM->usaPCA == 1)
  {
    // Usa PCA
    fprintf(arquivo,"PCA\n");
    
    fprintf(arquivo,"Usa PCA=sim\n");

    // Nome do arquivo com matriz de transformacao
    fprintf(arquivo,"Nome PCA=%s\n",HMM->nomePCA);

    // Ordem do vetor de parametros reduzido
    fprintf(arquivo,"Ordem reduzida=%d\n",HMM->ordemRed);
  }
  fclose(arquivo); 
}
