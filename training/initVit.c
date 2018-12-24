//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "estruturas.h"
#include "par.h"
#include "gausFunc.h"
#include "lbg.h"
#include "loadTranscription.h"
#include "forback.h"

// Funcao para inicializacao dos modelos das subunidades foneticas via
// Algoritmo de Viterbi
void Inicializa_Viterbi
(
	char **subunits,        // list of phonetic subunits
	struct config trainConfig, // training configuration settings  
	struct modelo_HMM *HMM // modelos HMM das subunidades foneticas
)
{
  // Declaracao das variaveis locais
  double **a; // matriz de transicao p/ os modelos das locucoes
  int ***A2; // matriz aux p/ o calculo da matriz de transicao p/ os fonemas
  FILE *arquivo; // handle para arquivo
  double aux; // var aux p/ o calculo das variancias das gaussianas
  int* caminho_otimo; // seq de estados determinada pelo algoritmo de Viterbi
  struct mistura **b; // matriz de emissao para os modelos das locucoes
  double ***bj; // probabilidade de cada simbolo da locucao, dado o modelo
  double distancia; // var aux p/ o calculo das variancias das gaussianas
  int estado; // identifica a que estado corresponde um dado vetor
  int register fone,frase,i,j,k,l,m,t; // contadores
  int fonema; // identifica a que fone corresponde um dado vetor
  double **medias; // var aux p/ o calculo das medias das gaussianas
  double minimo; // var aux p/ o calculo das variancias das gaussianas
  double ****Nm; // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
  char NomeArquivo[256]; // var aux p/ formacao de nomes de arquivos 
  int **numero_vetores; // armazena as contagens
  int qual; // var aux p/ o calculo das variancias das gaussianas
  int *quantos; // numero de vetores selecionados p/ cada estado
  int *quantosg; // numero de vetores selecionados para cada gaussiana
  int soma; // var aux p/ atualizacao das matrizes de transicao
  double **variancias; // var aux p/ o calculo das variancias das gaussianas
  double ****vetores; // lista de vetores para incializacao dos modelos
  struct v_psi vpsi; // var aux p/ determinacao do caminho otimo (Viterbi)
  struct locucao x; // vetores de parametros de treinamento
  int *transcription; // phonetic transcription of the training utterance
  int length; // number of phonetic subunits in phonetic transcription
  FILE *ptrUtterances;
  FILE *ptrTranscriptions;
  char wavFile[512]; // name of file with the utterance
  char transcFile[512]; // name of file with phonetic transcription of the utterance
  
  
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
  A2 = malloc(sizeof(int **)*HMM->n_fones);
  for (i=0;i<HMM->n_fones;i++)
    A2[i] = malloc(sizeof(int*)*HMM->n_estados);
  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      A2[i][j] = malloc(sizeof(int)*(HMM->salto_maximo+1));

  numero_vetores = malloc(sizeof(int *)*HMM->n_fones);
  for (i=0;i<HMM->n_fones;i++)
    numero_vetores[i] = malloc(sizeof(int)*HMM->n_estados);

  quantos = malloc(sizeof(int)*HMM->n_estados);

  // Zerando contagens para as matrizes de transicao
  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      for (k=0;k<(HMM->salto_maximo+1);k++)
        A2[i][j][k] = 0;

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
    printf("%cSeparando vetores: arquivo %d de %d.",13,frase+1,trainConfig.nUtterances);
	fflush(stdout);
	
    // Loading transcription
	fgets(transcFile,512,ptrTranscriptions); // reading filename
	
	// removing control characters (\r \n) frm the end of string
    while (transcFile[strlen(transcFile)-1]=='\r' || transcFile[strlen(transcFile)-1]=='\n')
  	  transcFile[strlen(transcFile)-1] = '\0';
	
	loadTranscription(transcFile,subunits,HMM->n_fones,&transcription,&length); // reading phonetic transcription

    // Criando modelos das frases
    // Matriz de transicao de estados
    a = malloc(sizeof(double*)*(HMM->n_estados*length));
    for (i=0;i<(HMM->n_estados*length);i++)
      a[i] = malloc(sizeof(double)*(HMM->salto_maximo+1));

    // Matriz de emissao
    b = malloc(sizeof(struct mistura *)*HMM->n_par);
    for (i=0;i<HMM->n_par;i++)
      b[i] = malloc(sizeof(struct mistura)*(HMM->n_estados*length));
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_estados*length;j++)
      {
        b[i][j].c = malloc(sizeof(double)*HMM->n_gaussianas[i]);
        b[i][j].m = malloc(sizeof(double *)*HMM->n_gaussianas[i]);
        b[i][j].v = malloc(sizeof(double *)*HMM->n_gaussianas[i]);
      }
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_estados*length;j++)
        for (k=0;k<HMM->n_gaussianas[i];k++)
        {
          b[i][j].m[k] = malloc(sizeof(double)*HMM->ordem[i]);
          b[i][j].v[k] = malloc(sizeof(double)*HMM->ordem[i]);
        }

    //
    for (i=0;i<length;i++)
      for (j=0;j<HMM->n_estados;j++)
        for (k=0;k<(HMM->salto_maximo+1);k++)
          a[j+i*HMM->n_estados][k] = HMM->A[transcription[i]][j][k];

    for (l=0;l<HMM->n_par;l++)
      for (i=0;i<length;i++)
        for (j=0;j<HMM->n_estados;j++)
          for (k=0;k<HMM->n_gaussianas[l];k++)
          {
            b[l][j+i*HMM->n_estados].c[k] = HMM->B[l][transcription[i]][j].c[k];
            for (m=0;m<HMM->ordem[l];m++)
            {
              b[l][j+i*HMM->n_estados].m[k][m] = HMM->B[l][transcription[i]][j].m[k][m];
              b[l][j+i*HMM->n_estados].v[k][m] = HMM->B[l][transcription[i]][j].v[k][m];
            }
          }

    // Parametrizando sinal de voz (aloca ponteiro 'x')
    fgets(wavFile, 500, ptrUtterances);
    
    // removing control characters (\r \n) frm the end of string
    while (wavFile[strlen(wavFile)-1]=='\r' || wavFile[strlen(wavFile)-1]=='\n')
  	  wavFile[strlen(wavFile)-1] = '\0';
    
    calcPar(wavFile,*HMM,&x.n_quadros,&x.par);

    // Estruturas para armazenamento da probabilidade de simbolo e da contribuicao
    // de cada gaussiana p/ a prob. simb. da locucao.
    // Tambem aloca ponteiros bj e Nm.
    ProbSimbolo(b,&bj,length,*HMM,&Nm,x.n_quadros,(double ***)x.par);

    // Alocando memoria
    vpsi.v = malloc(sizeof(double*)*HMM->n_estados*length);
    vpsi.psi = malloc(sizeof(int*)*HMM->n_estados*length);
    for (i=0;i<HMM->n_estados*length;i++)
    {
      vpsi.v[i] = malloc(sizeof(double)*x.n_quadros);
      vpsi.psi[i] = malloc(sizeof(int)*x.n_quadros);
    }

    caminho_otimo = malloc(sizeof(int)*(x.n_quadros+1));

    // Verificando caminho otimo (funcao Viterbi definida em 'ForBack.cpp')
    Viterbi(a,bj,length,*HMM,x.n_quadros,vpsi);

    // Recuperando caminho otimo
    caminho_otimo[x.n_quadros] = HMM->n_estados*length-1;
    for (i=(x.n_quadros-1);i>0;i--)
      caminho_otimo[i] = vpsi.psi[caminho_otimo[i+1]][i];
    caminho_otimo[0] = 0;

    // Atualizando matrizes de transicao
    for (t=0;t<x.n_quadros;t++)
    {
      l = caminho_otimo[t];
      j = caminho_otimo[t+1];
      A2[transcription[(j/3)]][l - (l/3)*HMM->n_estados][j - l]++;
    }

    // Separando vetores p/ atualizacao das medias e variancias das matrizes de emissao
    for (t=0;t<x.n_quadros;t++)
    {
      // Determinando a que fone corresponde o vetor
      fonema = caminho_otimo[t+1] / HMM->n_estados;

      // Determinando a que estado do fone corresponde o vetor
      estado = caminho_otimo[t+1] % HMM->n_estados;

      // Atualizando contagem dos vetores para cada estado de cada fonema
      numero_vetores[transcription[fonema]][estado]++;

      // Abrindo arquivo do fone e estado correspondentes
      sprintf(NomeArquivo,"fone%d.s%d",transcription[fonema],estado);
	  arquivo = fopen(NomeArquivo,"a+b"); 
      if (arquivo == NULL)
      {
      	printf("Error opening file %s.\n",NomeArquivo);
        exit(1);
      }
      
      // Armazenando vetor
      for (i=0;i<HMM->n_par;i++)
        for (l=0;l<HMM->ordem[i];l++)
          fwrite(&x.par[i][t][l],sizeof(double),1,arquivo);
      fclose(arquivo);
    } // end for (t)

    // Desalocando ponteiros
    for (i=0;i<(HMM->n_estados*length);i++)
      free(a[i]);
    free(a);

    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<(HMM->n_estados*length);j++)
        for (k=0;k<HMM->n_gaussianas[i];k++)
        {
          free(b[i][j].m[k]);
          free(b[i][j].v[k]);
        }
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<(HMM->n_estados*length);j++)
    {
      free(b[i][j].c);
      free(b[i][j].m);
      free(b[i][j].v);
    }
    for (i=0;i<HMM->n_par;i++)
      free(b[i]);
    free(b);

    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<x.n_quadros;j++)
        free(x.par[i][j]);
    for (i=0;i<HMM->n_par;i++)
      free(x.par[i]);
    free(x.par);

	free(transcription);
	
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<x.n_quadros;j++)
        for (k=0;k<HMM->n_estados*length;k++)
          free(Nm[i][j][k]);
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<x.n_quadros;j++)
        free(Nm[i][j]);
    for (i=0;i<HMM->n_par;i++)
      free(Nm[i]);
    free(Nm);

    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_estados*length;j++)
        free(bj[i][j]);
    for (i=0;i<HMM->n_par;i++)
      free(bj[i]);
    free(bj);

    for (i=0;i<HMM->n_estados*length;i++)
    {
      free(vpsi.v[i]);
      free(vpsi.psi[i]);
    }
    free(vpsi.v);
    free(vpsi.psi);

    free(caminho_otimo);
  } // end for(frase)
  printf("\n");

  fclose(ptrUtterances);
  fclose(ptrTranscriptions);

  // Atualizando matrizes de transicao
  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
    {
      soma = 0;
      for (k=0;k<(HMM->salto_maximo+1);k++)
        soma += A2[i][j][k];
      if (soma != 0)
        for (k=0;k<(HMM->salto_maximo+1);k++)
          HMM->A[i][j][k] = (double)A2[i][j][k] / (double)soma;
      else
        for (k=0;k<(HMM->salto_maximo+1);k++)
          HMM->A[i][j][k] = 0.0;
    }

  // Inicializando medias e variancias das misturas

  fone = 0;
  while (fone < HMM->n_fones)
  //for (fone=0;fone<HMM->n_fones;fone++)
  {
    // Mensagem de processamento
    printf("%cCalculando médias e variâncias: Fone %d de %d",13,fone+1,HMM->n_fones);
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
      	printf("Error creating file %s. Verify permissions or if disk is full.\n",NomeArquivo);
        exit(1);
      }

      for (j=0;j<numero_vetores[fone][i];j++)
        for (k=0;k<HMM->n_par;k++)
          for (l=0;l<HMM->ordem[k];l++)
            fread(&vetores[k][i][j][l],sizeof(double),1,arquivo);
      fclose(arquivo);
    } // end for (i)

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
      } // end for (j)
    } // end for (i)

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
  } // end while (fone)
  printf("\n");
  
  // Desalocando ponteiros
  for (i=0;i<HMM->n_fones;i++)
    free(numero_vetores[i]);
  free(numero_vetores);

  free(quantos);

  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      free(A2[i][j]);
  for (i=0;i<HMM->n_fones;i++)
    free(A2[i]);
  free(A2);

  // Apagando arquivos temporarios
  for (fone=0;fone<HMM->n_fones;fone++)
    for (i=0;i<HMM->n_estados;i++)
    {
      sprintf(NomeArquivo,"fone%d.s%d",fone,i);
      remove(NomeArquivo);
    }
}
