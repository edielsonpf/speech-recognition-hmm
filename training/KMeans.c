//---------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "estruturas.h"
#include "lbg.h"
#include "par.h"
#include "loadTranscription.h"
#include "KMeans.h"

//---------------------------------------------------------------------------

// Funcao para inicializacao dos modelos das subunidades foneticas via
// Segmental K-Means
void K_Means
(
	char **subunits,        // list of phonetic subunits
	struct config trainConfig, // training configuration settings  
	struct modelo_HMM *HMM // modelos HMM das subunidades foneticas
)
{
  // Declaracao das variaveis locais
  FILE *arquivo; // handle para arquivo
  double aux; // var aux p/ o calculo das variancias das gaussianas
  double distancia; // var aux p/ o calculo das variancias das gaussianas
  int final; // quadro final da selecao
  int register fone,frase,i,j,k,l,t; // contadores
  int inicio; // quadro inicial da selecao
  double **medias; // var aux p/ o calculo das medias das gaussianas
  double minimo; // var aux p/ o calculo das variancias das gaussianas
  char NomeArquivo[256]; // var aux p/ formacao de nomes de arquivos
  int n_simbolos; // numero de simbolos para cada estado
  int **numero_vetores; // armazena as contagens
  int qual; // var aux p/ o calculo das variancias das gaussianas
  int *quantos; // numero de vetores selecionados p/ cada estado
  int *quantosg; // numero de vetores selecionados para cada gaussiana
  int resto; // numero de estados que irao ter um simbolo a mais no calculo da media
  double **variancias; // var aux p/ o calculo das variancias das gaussianas
  double ****vetores; // lista de vetores para incializacao dos modelos
  struct locucao x; // vetores de parametros de treinamento
  FILE *ptrUtterances;
  FILE *ptrTranscriptions;
  char wavFile[512]; // name of file with the utterance
  char transcFile[512]; // name of file with phonetic transcription of the utterance
  int *transcription; // phonetic transcription of the training utterance
  int length; // number of phonetic subunits in phonetic transcription

  // Gerando arquivos para armazenar os parametros de treinamento de cada subunidade
  for (fone=0;fone<HMM->n_fones;fone++)
  {
    for (i=0;i<HMM->n_estados;i++)
    {
      sprintf(NomeArquivo,"fone%d.s%d",fone,i);

      arquivo = fopen(NomeArquivo,"wb");
      if (arquivo == NULL)
      {
      	printf("Error creating file %s. Verify permissions or if disk is full.\n",NomeArquivo);
        exit(1);
      }
      else
        fclose(arquivo);
    }
  }

  // Alocando memoria
  numero_vetores = malloc(sizeof(int *)*HMM->n_fones);
  for (i=0;i<HMM->n_fones;i++)
    numero_vetores[i] = malloc(sizeof(int)*HMM->n_estados);

  quantos = malloc(sizeof(int)*HMM->n_estados);

  // Inicializando variavel que armazena o numero de vetores para cada estado de cada subunidade
  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      numero_vetores[i][j] = 0;

  ptrUtterances = fopen(trainConfig.utterancesFiles,"rt");
  if(ptrUtterances==NULL)
  {
    puts("Error opening utterances list file\n");
    exit(1);
  } 
  
  ptrTranscriptions = fopen(trainConfig.transcriptionsFiles,"rt");
  if(ptrTranscriptions==NULL)
  {
    puts("Error opening transcriptions list file\n");
    exit(1);
  } 
  
  
  for (frase=0;frase<trainConfig.nUtterances;frase++)
  {

    // Mensagem de processamento
    printf("%cSeparando vetores: arquivo %d de %d",13,frase+1,trainConfig.nUtterances);
    fflush(stdout);
	
    // Loading transcription
	fscanf(ptrTranscriptions,"%s",transcFile); // reading filename
	loadTranscription(transcFile,subunits,HMM->n_fones,&transcription,&length); // reading phonetic transcription
	
    // Parametrizando sinal de voz (aloca ponteiro 'x')
	fscanf(ptrUtterances,"%s",wavFile);
    calcPar(wavFile,*HMM,&x.n_quadros,&x.par);

    // Verificando quantos simbolos para cada estado
    n_simbolos = x.n_quadros / (HMM->n_estados*length);
    resto = x.n_quadros % (HMM->n_estados*length);

    // Separando vetores
    inicio = 0; // vetor inicial da selecao
    final = 0; // vetor final da selecao
    for (j=0;j<length;j++)
    {
      for (k=0;k<HMM->n_estados;k++)
      {
        // Numero de vetores para este estado
        if ((k+j*HMM->n_estados) < resto)
        {
          final = inicio + n_simbolos + 1;
          numero_vetores[transcription[j]][k] += n_simbolos+1;
        }
        else
        {
          final = inicio + n_simbolos;
          numero_vetores[transcription[j]][k] += n_simbolos;
        }

        // Abrindo arquivo do fone e estado correspondentes
        sprintf(NomeArquivo,"fone%d.s%d",transcription[j],k);

		arquivo = fopen(NomeArquivo,"a+b");
		if (arquivo == NULL)
		{
		  printf("Error opening file %s.\n",NomeArquivo);
          exit(1);
		}
		
        // Armazenando vetores a lista
        for (t=inicio;t<final;t++)
          for (i=0;i<HMM->n_par;i++)
            for (l=0;l<HMM->ordem[i];l++)
              fwrite(&x.par[i][t][l],sizeof(double),1,arquivo);
        fclose(arquivo);

        // Vetor inicial para o proximo estado
        inicio = final;
      }
    }

    // Desalocando ponteiro com os parametros
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<x.n_quadros;j++)
        free(x.par[i][j]);
    for (i=0;i<HMM->n_par;i++)
      free(x.par[i]);
    free(x.par);
    free(transcription);
  }  // end for (frase)
  printf("\n");
  
  fclose(ptrUtterances);
  fclose(ptrTranscriptions);
  
  // Inicializando medias e variancias das misturas
  fone = 0;
  while (fone < HMM->n_fones) 
  {
    // Mensagem de processamento
    printf("%cCalculado médias e variâncias: Fone %d de %d",13,fone+1,HMM->n_fones);
	fflush(stdout);
	
    // Alocando memoria para matriz de parametros correspondentes a este fone
    vetores = malloc(sizeof(double ***)*HMM->n_par);
    for (i=0;i<HMM->n_par;i++)
      vetores[i] = malloc(sizeof(double **)*HMM->n_estados);
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_estados;j++)
        vetores[i][j] = malloc(sizeof(double *)*numero_vetores[fone][j]);
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_estados;j++)
        for (k=0;k<(numero_vetores[fone][j]);k++)
          vetores[i][j][k] = malloc(sizeof(double)*HMM->ordem[i]);

    for (i=0;i<HMM->n_estados;i++)
    {
      sprintf(NomeArquivo,"fone%d.s%d",fone,i);	

	  arquivo = fopen(NomeArquivo,"rb");
      if (arquivo == NULL)
      {
        printf("Error opening file %s.\n",NomeArquivo);	
        exit(1);
      }
      
      for (j=0;j<numero_vetores[fone][i];j++)
        for (k=0;k<HMM->n_par;k++)
          for (l=0;l<HMM->ordem[k];l++)
            fread(&vetores[k][i][j][l],sizeof(double),1,arquivo);
      fclose(arquivo);
    }

    // Inicializando misturas
    for (i=0;i<HMM->n_par;i++)
    {
      for (j=0;j<HMM->n_estados;j++)
      {
        // Medias
        medias = malloc(sizeof(double *)*HMM->n_gaussianas[i]);
        for (k=0;k<HMM->n_gaussianas[i];k++)
          medias[k] = malloc(sizeof(double)*HMM->ordem[i]);

        variancias = malloc(sizeof(double *)*HMM->n_gaussianas[i]);
        for (k=0;k<HMM->n_gaussianas[i];k++)
          variancias[k] = malloc(sizeof(double)*HMM->ordem[i]);

        quantosg = malloc(sizeof(int)*HMM->n_gaussianas[i]);

        lbg(&medias,HMM->n_gaussianas[i],numero_vetores[fone][j],HMM->ordem[i],vetores[i][j]);

        for (k=0;k<HMM->n_gaussianas[i];k++)
          for (l=0;l<HMM->ordem[i];l++)
            HMM->B[i][fone][j].m[k][l] = medias[k][l];

        // Inicialziando variancias
        for (k=0;k<HMM->n_gaussianas[i];k++)
          for (l=0;l<HMM->ordem[i];l++)
            variancias[k][l] = 0.0;

        for (k=0;k<HMM->n_gaussianas[i];k++)
          quantosg[k] = 0;

        for (k=0;k<numero_vetores[fone][j];k++)
        {
          // Verificando a qual gaussiana pertence o vetor sob analise
          minimo = 1.0e20;
          for (l=0;l<HMM->n_gaussianas[i];l++)
          {
            distancia = deucl(vetores[i][j][k],medias[l],HMM->ordem[i]);
            if (distancia < minimo)
            {
              minimo = distancia;
              qual = l;
            }
          }

          // Atualizando contagem dos vetores p/ cada gaussiana
          quantosg[qual]++;

          // Somando valores
          for (l=0;l<HMM->ordem[i];l++)
          {
            aux = vetores[i][j][k][l] - medias[qual][l];
            aux *= aux;
            variancias[qual][l] += aux;
          }
        }
        for (k=0;k<HMM->n_gaussianas[i];k++)
          for (l=0;l<HMM->ordem[i];l++)
            variancias[k][l] /= quantosg[k];

        for (k=0;k<HMM->n_gaussianas[i];k++)
          for (l=0;l<HMM->ordem[i];l++)
            HMM->B[i][fone][j].v[k][l] = variancias[k][l];

        // Coeficientes das gaussianas
        aux = 0.0;
        for (k=0;k<HMM->n_gaussianas[i];k++)
          aux += quantosg[k];
        for (k=0;k<HMM->n_gaussianas[i];k++)
          HMM->B[i][fone][j].c[k] = quantosg[k] / aux;

        // Desalocando ponteiros
        for (k=0;k<HMM->n_gaussianas[i];k++)
          free(medias[k]);
        free(medias);

        for (k=0;k<HMM->n_gaussianas[i];k++)
          free(variancias[k]);
        free(variancias);

        free(quantosg);
      }
    }

    // Desalocando 'vetores'
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_estados;j++)
        for (k=0;k<numero_vetores[fone][j];k++)
          free(vetores[i][j][k]);
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_estados;j++)
        free(vetores[i][j]);
    for (i=0;i<HMM->n_par;i++)
      free(vetores[i]);
    free(vetores);

    fone++; // proximo fone
  } // end while(fones)
  printf("\n");

  // Desalocando ponteiro n_vetores
  for (i=0;i<HMM->n_fones;i++)
    free(numero_vetores[i]);
  free(numero_vetores);

  free(quantos);

  // Apagando arquivos temporarios
  for (fone=0;fone<HMM->n_fones;fone++)
    for (i=0;i<HMM->n_estados;i++)
    {
      sprintf(NomeArquivo,"fone%d.s%d",fone,i);	
      remove(NomeArquivo);
    }
}
