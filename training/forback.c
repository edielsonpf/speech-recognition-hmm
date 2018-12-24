#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//---------------------------------------------------------------------------
#include "estruturas.h"
#include "constantes.h"
#include "KMeans.h"
#include "HMMFiles.h"
#include "par.h"
#include "gausFunc.h"
#include "loadTranscription.h"
#include "initVit.h"
#include "forback.h"

//---------------------------------------------------------------------------

void forwardBackward
(
	char **subunits,        // list of phonetic subunits
	struct config trainConfig, // training configuration settings  
	struct modelo_HMM *HMM // modelos HMM das subunidades foneticas
)
{
	// Declaracao das variaveis locais
    double **a; // matriz de transicao p/ os modelos das locucoes
    double **alfa; // variavel FORWARD
    double **beta; // variavel BACKWARD
    double ***bj; // probabilidade de cada simbolo da locucao, dado o modelo
    double **den_A; // var aux para atualizacao das probs de transicao dos modelos das subunidades
    double ****den_B_c; // var aux para atualizacao das probs de emissao dos modelos das subunidades
    double distorcao; // derivada da variacao de p_media
    int frase; // conta as locucoes de treinamento
    struct mistura **b; // matriz de emissao para os modelos das locucoes
    int epoca; // conta as epocas de treinamento
    double ***num_A; // var aux para atualizacao das probs de transicao dos modelos das subunidades
    struct mistura ***num_B; // var aux para atualizacao das probs de emissao dos modelos das subunidades
    double *c; // fator de escala
    int register i,j,k,l,m; // contadores
    int nArquivosConv; // numero de arquivos utilizados para a verificacao da convergencia
    double ****Nm; // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
    int nSubstituicoes; // numero de vezes que o sistema teve que atribuir um valor minimo para coeficientes de ponderacao e variancias das gaussianas
    double p_media[2]; // P(O/M) media para todas as locucoes(p_media[1] epoca atual, p_media[0] epoca anterior)
    struct locucao x; // vetores de parametros de treinamento
    FILE *ptrUtterances;
    FILE *ptrTranscriptions;
    char wavFile[512]; // name of file with the utterance
    char transcFile[512]; // name of file with phonetic transcription of the utterance
    int *transcription; // phonetic transcription of the training utterance
    int length; // number of phonetic subunits in phonetic transcription
     
    // Generating information file
    GeraArquivoInformacao(HMM,trainConfig.HMMFile);

    // Numero de arquivos utilizados para verificacao da convergencia do treinamento
    if (trainConfig.nUtterances > 100)
      nArquivosConv = (int) (trainConfig.nUtterances/10);
    else
      nArquivosConv = trainConfig.nUtterances;

    // Allocating memory

    // Matriz de transicao para os modelos HMM das subunidades
    HMM->A = malloc(sizeof(double **)*HMM->n_fones);
    for (i=0;i<HMM->n_fones;i++)
      HMM->A[i] = malloc(sizeof(double *)*HMM->n_estados);
    for (i=0;i<HMM->n_fones;i++)
      for (j=0;j<HMM->n_estados;j++)
        HMM->A[i][j] = malloc(sizeof(double)*(HMM->salto_maximo+1));

    // Matrizes aux p/ atualizacao da matriz A
    num_A = malloc(sizeof(double **)*HMM->n_fones);
    den_A = malloc(sizeof(double *)*HMM->n_fones);
    for (i=0;i<HMM->n_fones;i++)
    {
      num_A[i] = malloc(sizeof(double *)*HMM->n_estados);
      den_A[i] = malloc(sizeof(double)*HMM->n_estados);
    }
    for (i=0;i<HMM->n_fones;i++)
      for (j=0;j<HMM->n_estados;j++)
        num_A[i][j] = malloc(sizeof(double)*(HMM->salto_maximo+1));

    // Matriz de emissao para os modelos HMM das subunidades
    HMM->B = malloc(sizeof(struct mistura **)*HMM->n_par);
    for (i=0;i<HMM->n_par;i++)
      HMM->B[i] = malloc(sizeof(struct mistura *)*HMM->n_fones);
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_fones;j++)
        HMM->B[i][j] = malloc(sizeof(struct mistura)*HMM->n_estados);
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_fones;j++)
        for (k=0;k<HMM->n_estados;k++)
        {
          HMM->B[i][j][k].c = malloc(sizeof(double)*HMM->n_gaussianas[i]);
          HMM->B[i][j][k].m = malloc(sizeof(double *)*HMM->n_gaussianas[i]);
          HMM->B[i][j][k].v = malloc(sizeof(double *)*HMM->n_gaussianas[i]);
        }
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_fones;j++)
        for (k=0;k<HMM->n_estados;k++)
          for (l=0;l<HMM->n_gaussianas[i];l++)
          {
            HMM->B[i][j][k].m[l] = malloc(sizeof(double)*HMM->ordem[i]);
            HMM->B[i][j][k].v[l] = malloc(sizeof(double)*HMM->ordem[i]);
          }

    // Matrizes aux p/ atualizacao da matriz B
    num_B = malloc(sizeof(struct mistura **)*HMM->n_par);
    den_B_c = malloc(sizeof(double ***)*HMM->n_par);
    for (i=0;i<HMM->n_par;i++)
    {
      num_B[i] = malloc(sizeof(struct mistura *)*HMM->n_fones);
      den_B_c[i] = malloc(sizeof(double **)*HMM->n_fones);
    }
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_fones;j++)
      {
        num_B[i][j] = malloc(sizeof(struct mistura)*HMM->n_estados);
        den_B_c[i][j] = malloc(sizeof(double *)*HMM->n_estados);
      }
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_fones;j++)
        for (k=0;k<HMM->n_estados;k++)
        {
          num_B[i][j][k].c = malloc(sizeof(double)*HMM->n_gaussianas[i]);
          num_B[i][j][k].m = malloc(sizeof(double *)*HMM->n_gaussianas[i]);
          num_B[i][j][k].v = malloc(sizeof(double *)*HMM->n_gaussianas[i]);
          den_B_c[i][j][k] = malloc(sizeof(double)*HMM->n_gaussianas[i]);
        }
    for (i=0;i<HMM->n_par;i++)
      for (j=0;j<HMM->n_fones;j++)
        for (k=0;k<HMM->n_estados;k++)
          for (l=0;l<HMM->n_gaussianas[i];l++)
          {
            num_B[i][j][k].m[l] = malloc(sizeof(double)*HMM->ordem[i]);
            num_B[i][j][k].v[l] = malloc(sizeof(double)*HMM->ordem[i]);
          }

    // Inicializacao dos modelos HMM das subunidades

    // Matriz de transicao de estados (matriz 'A')
    for (i=0;i<HMM->n_fones;i++)
    {
      HMM->A[i][0][0] = 1.0/3.0; HMM->A[i][0][1] = 1.0/3.0; HMM->A[i][0][2] = 1.0/3.0;
      HMM->A[i][1][0] = 0.5;     HMM->A[i][1][1] = 0.5;     HMM->A[i][1][2] = 0.0;
      HMM->A[i][2][0] = 0.5;     HMM->A[i][2][1] = 0.5;     HMM->A[i][2][2] = 0.0;
    }

    // Matriz de emissao (matriz 'B')
    switch (trainConfig.initializationMode)
    {
    	case 0: // Inicializacao por distribuicao uniforme
      		for (i=0;i<HMM->n_par;i++)
		        for (j=0;j<HMM->n_fones;j++)
        			for (k=0;k<HMM->n_estados;k++)
			            for (l=0;l<HMM->n_gaussianas[i];l++)
             			{
                			HMM->B[i][j][k].c[l] = 1.0 / (double)HMM->n_gaussianas[i];
                			for (m=0;m<HMM->ordem[i];m++)
                			{
                  				HMM->B[i][j][k].m[l][m] = 10.0 + (double)l/10.0;
                  				HMM->B[i][j][k].v[l][m] = 100.0 + (double)l/10.0;
                			}
              			}
      	break;

    	case 1: // Inicializacao pelo algoritmo Segmental K-means

	    // Passo 1: Segmentacao uniforme
        printf("Inicializando modelos via Segmental K-Means\n");
        printf("Passo1: segmentação uniforme\n");

        // Inicializando modelos com segmentacao uniforme
        K_Means(subunits,trainConfig,HMM);

        // Substituindo valores muito pequenos dos coeficientes de ponderacao e das
        // variancias das gaussianas por um valor minimo (variavel MINIMO em 'HMM.cpp')
        nSubstituicoes = ValoresMinimos(HMM);
        if (nSubstituicoes != 0)
        	printf("%d substituições realizadas.\n",nSubstituicoes);

    	// Verificando probabilidades de transicao e coeficientes das gaussianas
    	TestaHMM(HMM);

	    // Salvando modelos em arquivo
    	SalvaHMM(HMM,trainConfig.HMMFile,-2);

    	// Calculando verossimilhanca media com estes modelos
    	p_media[1] = CalculaProbMedia(subunits,trainConfig,HMM,nArquivosConv);

	    // Apresentando mensagens
    	printf("P(O/M): %f\n",p_media[1]);

    	// Passo 2: Segmentacao via algoritmo de Viterbi
    	printf("Passo 2: segmentação via algoritmo de Viterbi\n");

	    distorcao = 1.0;
    	epoca = 0;
	    while (distorcao > 0.001)
    	{
        	// Mensagem de processamento
	        printf("Época: %d\n",epoca+1);

    	    // Inicializando modelos via segmentacao por Viterbi
        	Inicializa_Viterbi(subunits,trainConfig,HMM);

	        // Substituindo valores muito pequenos dos coeficientes de ponderacao e das
    	    // variancias das gaussianas por um valor minimo (variavel MINIMO em 'HMM.cpp')
        	nSubstituicoes = ValoresMinimos(HMM);
	        if (nSubstituicoes != 0)
    	        printf("%d substituições realizadas.\n",nSubstituicoes);

        	// Verificando probabilidades de transicao e coeficientes das gaussianas
	        TestaHMM(HMM);

    	    // Salvando modelos
        	SalvaHMM(HMM,trainConfig.HMMFile,-1);

	        // Calculando verossimilhanca media com estes modelos
    	    p_media[0] = p_media[1];
        	p_media[1] = CalculaProbMedia(subunits,trainConfig,HMM,nArquivosConv);
	        // Calculando distorcao
    	    distorcao = (p_media[1]-p_media[0])/fabs(p_media[1]);
        	if (distorcao < 0.0)
        		distorcao = 1.0;

	        // Apresentando mensagens
    	    printf("P(O/M): %f\n",p_media[1]);
        	printf("Distorção: %f\n",distorcao);

	        // Nova epoca
    	    epoca++;

        	//distorcao = 0.0; // usa-se esta linha para usar o Viterbi uma epoca apenas
		} // end while (distorcao)
    	break;
        
    	default: // Carregando modelo pre-treinado
    		loadHMM(trainConfig.preTrained,HMM);
    
    		// Calculando verossimilhanca media com estes modelos
    		p_media[1] = CalculaProbMedia(subunits,trainConfig,HMM,nArquivosConv);
    		distorcao = 1.0;

    		// Apresentando mensagens
    		printf("P(O/M): %f\n",p_media[1]);
	} // end switch (initializationMode)


  	// Treinando o sistema
  	printf("Treinando modelos das subunidades via Baum-Welch\n");

	if (trainConfig.initializationMode == 0) // Inicializando p_media e distorcao p/ inicializacao padrao
  	{
    	p_media[0] = -Inf;
    	p_media[1] = -Inf;
    	distorcao = 1.0;
  	}
  	else
    	p_media[0] = p_media[1]; // K-Means e modelos pre-treinados

  	distorcao = 1.0;
  	//distorcao = 0.0; // linha usada para treinamento usando apenas o algoritmo de Viterbi
  	epoca = 0;
  	while (distorcao > 0.001) 
  	{
    	// Mensagem de processamento
    	printf("Época: %d\n",epoca+1);

    	// Inicializacao das variaveis auxiliares num_A, den_A, num_B, den_B
    	for (i=0;i<HMM->n_fones;i++)
      		for (j=0;j<HMM->n_estados;j++)
      		{
        		den_A[i][j] = 0.0;
        		for (k=0;k<(HMM->salto_maximo+1);k++)
          			num_A[i][j][k] = 0.0;
      		}

    	for (i=0;i<HMM->n_par;i++)
      		for (j=0;j<HMM->n_fones;j++)
        		for (k=0;k<HMM->n_estados;k++)
         			for (l=0;l<HMM->n_gaussianas[i];l++)
          			{
            			num_B[i][j][k].c[l] = 0.0;
            			den_B_c[i][j][k][l] = 0.0;

            			for (m=0;m<HMM->ordem[i];m++)
            			{
              				num_B[i][j][k].m[l][m] = 0.0;
              				num_B[i][j][k].v[l][m] = 0.0;
            			}
          			}

		// Processando as locucoes
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
      		printf("%cTreinando modelos: arquivo %d de %d",13,frase+1,trainConfig.nUtterances);
      		fflush(stdout);

			// Loading transcription
			fgets(transcFile,512,ptrTranscriptions); // reading filename
			
			// removing control characters (\r \n) frm the end of string
    		while (transcFile[strlen(transcFile)-1]=='\r' || transcFile[strlen(transcFile)-1]=='\n')
  	  			transcFile[strlen(transcFile)-1] = '\0';
			
			loadTranscription(transcFile,subunits,HMM->n_fones,&transcription,&length); // reading phonetic transcription
			
      		// Allocating memory
      		// Matriz de transicao de estados
	        a = malloc(sizeof(double*)*(HMM->n_estados*length));
    	    for (i=0;i<(HMM->n_estados*length);i++)
        		a[i] = malloc(sizeof(double)*(HMM->salto_maximo+1));

        	// Matriz de emissao
        	b = malloc(sizeof(struct mistura *)*HMM->n_par);
        	for (i=0;i<HMM->n_par;i++)
          		b[i] = malloc(sizeof(struct mistura)*HMM->n_estados*length);
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

			
			// Criando modelos das frases
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
			// Variavel FORWARD
        	alfa = malloc(sizeof(double *)*length*HMM->n_estados);
        	for (i=0;i<(length*HMM->n_estados);i++)
          		alfa[i] = malloc(sizeof(double)*x.n_quadros);

        	// Veriavel BACKWARD
        	beta = malloc(sizeof(double *)*length*HMM->n_estados);
        	for (i=0;i<(length*HMM->n_estados);i++)
          		beta[i] = malloc(sizeof(double)*x.n_quadros);

        	// Fator de escala
        	c = malloc(sizeof(double)*x.n_quadros);

	      	// Calculando variaveis FORWARD
      		Foward(a,alfa,b,bj,c,length,HMM,x.n_quadros);

      		// Calculando variaveis BACKWARD
      		Backward(a,b,beta,bj,c,length,HMM,x.n_quadros);

      		// Somando contagens p/ atualizacao dos modelos das subunidades
      		SomaContagens(a,alfa,b,beta,bj,c,length,den_A,den_B_c,transcription,HMM,Nm,num_A,num_B,x);

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
      		for (i=0;i<(length*HMM->n_estados);i++)
        		free(alfa[i]);
      		free(alfa);

      		for (i=0;i<(length*HMM->n_estados);i++)
        		free(beta[i]);
      		free(beta);

      		free(c);

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
		} // fim loop for frase
		printf("\n");
		
    	// Atualizando modelos das subunidades apos uma epoca de treinamento
    	// Atualizando matrizes de transicao
      	Atualiza_A(den_A,HMM,num_A);

      	// Atualizando matrizes de emissao
      	Atualiza_B(den_B_c,HMM,num_B);

      	// Substituindo valores muito pequenos dos coeficientes de ponderacao e das
      	// variancias das gaussianas por um valor minimo (variavel MINIMO em 'HMM.cpp')
      	nSubstituicoes = ValoresMinimos(HMM);
      	if (nSubstituicoes != 0)
        	printf("%d substituições realizadas.\n",nSubstituicoes);

      	// Verificando probabilidades de transicao e coeficientes das gaussianas
      	TestaHMM(HMM);

      	// Salvando parametros treinados em arquivo
      	SalvaHMM(HMM,trainConfig.HMMFile,epoca);

      	// Calculando P(O/M)
      	p_media[0] = p_media[1];
      	p_media[1] = CalculaProbMedia(subunits,trainConfig,HMM,nArquivosConv);

      	// Apresentando mensagens
        printf("P(O/M): %f\n",p_media[1]);

      	// Calculando distorcao
      	distorcao = (p_media[1]-p_media[0])/fabs(p_media[1]);
      	if (distorcao < 0.0)
        	distorcao = 1.0;

      	// Apresentando mensagens
      	printf("Distorção: %f\n",distorcao);

      	// Salvando informacoes em arquivo
      	//FormMain->RichEdit1->Lines->SaveToFile("relat.doc");
    
    	// Nova epoca
    	epoca++;
  	} // end while (distorcao > 0.001)

  	// Desalocando ponteiros
  	for (i=0;i<HMM->n_fones;i++)
    	for (j=0;j<HMM->n_estados;j++)
	      	free(num_A[i][j]);
  	for (i=0;i<HMM->n_fones;i++)
  	{	
    	free(num_A[i]);
    	free(den_A[i]);
  	}
  	free(num_A);
  	free(den_A);

  	for (i=0;i<HMM->n_par;i++)
    	for (j=0;j<HMM->n_fones;j++)
      		for (k=0;k<HMM->n_estados;k++)
        		for (l=0;l<HMM->n_gaussianas[i];l++)
        		{
          			free(num_B[i][j][k].m[l]);
          			free(num_B[i][j][k].v[l]);
        		}
  	for (i=0;i<HMM->n_par;i++)
    	for (j=0;j<HMM->n_fones;j++)
      		for (k=0;k<HMM->n_estados;k++)
      		{
        		free(num_B[i][j][k].c);
        		free(num_B[i][j][k].m);
        		free(num_B[i][j][k].v);
        		free(den_B_c[i][j][k]);
      		}
  	for (i=0;i<HMM->n_par;i++)
    	for (j=0;j<HMM->n_fones;j++)
    	{
      		free(num_B[i][j]);
      		free(den_B_c[i][j]);
    	}
  	for (i=0;i<HMM->n_par;i++)
  	{
    	free(num_B[i]);
    	free(den_B_c[i]);
  	}
  	free(num_B);
  	free(den_B_c);
} // end  forwardBackward

//------------------------------------------------------------------------------
// Funcao que implementa o passo FORWARD do algoritmo Baum-Welch
double Foward
(
  double **a,         // matriz de transicao do modelo HMM da locucao
  double **alfa,    // variavel FORWARD
  struct mistura **b, // matriz de emissao do modelo HMM da locucao
  double ***bj,       // probabilidade de cada simbolo da locucao, dado o modelo
  double *c,        // fator de escala
  int comp,           // numero de subunidades na locucao
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  long tamanho        // comprimento da sequencia de entrada
)
{
  // Declaracao das variaveis locais
  double aux;         // variavel auxiliar
  int register i,j,k; // Indicador de estado.
  int register t;     // Indicador de tempo.

  // Inicializando o fator de escala
  for(t=0;t<tamanho;t++)
    c[t] = 0.0;

  // Inicializando a matriz alfa
  for(i=0;i<(HMM->n_estados*comp);i++)
    for(t=0;t<tamanho;t++)
      alfa[i][t] = 0.0;

  // Inicializacao
  alfa[0][0] = 1.0;
  for (i=0;i<HMM->n_par;i++)
    alfa[0][0] *= bj[i][0][0];
  c[0] = 1.0 / alfa[0][0];
  alfa[0][0] *= c[0];

  // Inducao
  for(t=0;t<(tamanho-1);t++)
  {
    for(i=0;i<(HMM->n_estados*comp);i++)
    {
      alfa[i][(t+1)] = 0.0;
      for (j=0;j<(HMM->salto_maximo+1);j++)
        if ((i-j) >= 0)
       	  alfa[i][t+1] += alfa[(i-j)][t]*a[(i-j)][j];

      for (k=0;k<HMM->n_par;k++)
        alfa[i][t+1] *= bj[k][i][t+1];

      c[(t+1)] += alfa[i][(t+1)];
    }
    c[(t+1)] = 1.0/c[(t+1)];

    for(i=0;i<(HMM->n_estados*comp);i++)
      alfa[i][(t+1)] *= c[(t+1)];
  }

  // Calculando P(O|M)
  aux = 0.0;
  for (t=0;t<tamanho;t++)
    aux += log10(c[t]);

  return(-aux);

}

//------------------------------------------------------------------------------
// Funcao que implementa o passo BACKWARD do algoritmo Baum-Welch
void Backward
(
  double **a,         // matriz de transicao do modelo HMM da locucao
  struct mistura **b, // matriz de emissao do modelo HMM da locucao
  double **beta,    // variavel BACKWARD
  double ***bj,       // probabilidade de cada simbolo da locucao, dado o modelo
  double *c,          // fator de escala
  int comp,           // numero de subunidades na locucao
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  long tamanho        // numero de quadros na locucao
)
{
  // Declaracao das variaveis locais
  double aux; // var aux p/ calcular as variaveis BACKWARD
  int register i,j,k; // Indicador de estado.
  int register t;   // Indicador de simbolo.

  // Inicializa a matriz beta
  for(i=0;i<(HMM->n_estados*comp);i++)
    for(t=0;t<(tamanho-1);t++)
      beta[i][t] = 0.0;

  // Calculo de beta
  // Inicializacao
  for(i=0;i<(HMM->n_estados*comp);i++)
    beta[i][tamanho-1] = 1.0 * c[tamanho-1];

  // Inducao
  for(t=(int)tamanho-2;t>=0;t--)
  {
    for (i=0;i<HMM->n_estados*comp;i++)
    {
      for (j=0;j<(HMM->salto_maximo+1);j++)
        if ((j+i) < HMM->n_estados*comp)
        {
      	  aux = beta[(i+j)][t+1] * a[i][j];
    	  for (k=0;k<HMM->n_par;k++)
            aux *= bj[k][i+j][t+1];
     	  beta[i][t] += aux;
        }
    }
    for (i=0;i<HMM->n_estados*comp;i++)
      beta[i][t] *= c[t];
  }
}

//------------------------------------------------------------------------------
// Funcao que soma as contagens de cada locucao para atualizacao das matrizes de
// transicao e emissao dos modelos HMM das subunidades
void SomaContagens
(
  double **a,                // matriz de transicao do modelo HMM da locucao
  double **alfa,             // variavel FORWARD
  struct mistura **b,        // matriz de emissao do modelo HMM da locucao
  double **beta,             // variavel BACKWARD
  double ***bj,              // probabilidade de cada simbolo da locucao, dado o modelo
  double *c,                 // fator de escala
  int comp,                  // numero de subunidades na locucao
  double **den_A,            // var aux p/ atualizacao da matriz A
  double ****den_B_c,        // var aux para atualizacao das probs de emissao dos modelos das subunidades
  int *modelo,               // modelo fonetico da locucao de treinamento
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  double ****Nm,             // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
  double ***num_A,           // var aux p/ atualizacao da matriz A
  struct mistura ***num_B,   // var aux para atualizacao das probs de emissao dos modelos das subunidades
  struct locucao x           // locucao parametrizada
)
{
  // Declaracao das variaveis locais
  double aux,aux1; // variaveis auxiliares
  int register i,j,k,l,m,t; // contadores

  // Matrizes de transicao
  // Atualizando num_A
  for (i=0;i<comp;i++)
    for (j=0;j<HMM->n_estados;j++)
      for (k=0;k<(HMM->salto_maximo+1);k++)
        if ((i*HMM->n_estados+j+k) < (comp*HMM->n_estados))
          for (t=0;t<(x.n_quadros-1);t++)
          {
            aux = 1.0;
            for (l=0;l<HMM->n_par;l++)
              aux *= bj[l][i*HMM->n_estados+j+k][t+1];
            aux *= alfa[i*HMM->n_estados+j][t]*a[i*HMM->n_estados+j][k]*beta[i*HMM->n_estados+j+k][t+1];

            num_A[modelo[i]][j][k] += aux;
          }

  // Atualizando den_A
  for (i=0;i<comp;i++)
    for (j=0;j<HMM->n_estados;j++)
      for (t=0;t<(x.n_quadros-1);t++)
        den_A[modelo[i]][j] += alfa[i*HMM->n_estados+j][t]*beta[i*HMM->n_estados+j][t]/c[t];

  // Matrizes de emissao
  // Atualizando num_B
  for (i=0;i<comp;i++)
    for (j=0;j<HMM->n_par;j++)
      for (k=0;k<HMM->n_estados;k++)
        for (l=0;l<HMM->n_gaussianas[j];l++)
          for (t=0;t<x.n_quadros;t++)
          {
            aux = alfa[i*HMM->n_estados+k][t]*beta[i*HMM->n_estados+k][t];
            aux *= Nm[j][t][i*HMM->n_estados+k][l];
            aux /= c[t];

            // Coeficientes das gaussianas
            num_B[j][modelo[i]][k].c[l] += aux;

            // Medias das gaussianas
            for (m=0;m<HMM->ordem[j];m++)
            {
              aux1 = aux*x.par[j][t][m];
              num_B[j][modelo[i]][k].m[l][m] += aux1;
            }
            // Variancias das gaussianas
            for (m=0;m<HMM->ordem[j];m++)
            {
              aux1 = x.par[j][t][m] - b[j][i*HMM->n_estados+k].m[l][m];
              aux1 *= aux1*aux;
              num_B[j][modelo[i]][k].v[l][m] += aux1;
            }
          }

  // Atualizando den_B
  for (i=0;i<comp;i++)
    for (j=0;j<HMM->n_par;j++)
      for (k=0;k<HMM->n_estados;k++)
        for (l=0;l<HMM->n_gaussianas[j];l++)
        {
          for (t=0;t<x.n_quadros;t++)
          {
            aux = alfa[i*HMM->n_estados+k][t]*beta[i*HMM->n_estados+k][t]/c[t];

            // Coeficientes das gaussianas
            den_B_c[j][modelo[i]][k][l] += aux;
          }
        }
}

//------------------------------------------------------------------------------
// Funcao que atualiza as matrizes de transicao das subunidades depois de cada
// epoca de treinamento
void Atualiza_A
(
  double **den_A,    // var aux p/ atualizacao da matriz A
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  double ***num_A    // var aux p/ atualizacao da matriz A

)
{
  // Declaracao das variaveis locais
  int register i,j,k; // contadores

  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
      for (k=0;k<(HMM->salto_maximo+1);k++)
        HMM->A[i][j][k] = num_A[i][j][k]/den_A[i][j];

}

//------------------------------------------------------------------------------
// Funcao que atualiza as matrizes de emissao das subunidades depois de cada
// epoca de treinamento
void Atualiza_B
(
  double ****den_B_c,       // var aux p/ atualizacao de B
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  struct mistura ***num_B  // var aux p/ atualizacao de B
)
{
  // Declaracao das variaveis locais
  int register i,j,k,l,m; // contadores

  for (i=0;i<HMM->n_par;i++)
    for (j=0;j<HMM->n_fones;j++)
      for (k=0;k<HMM->n_estados;k++)
        for (l=0;l<HMM->n_gaussianas[i];l++)
        {
          // Coeficientes das gaussianas
          HMM->B[i][j][k].c[l] = num_B[i][j][k].c[l] / den_B_c[i][j][k][l];

          for (m=0;m<HMM->ordem[i];m++)
          {
            // Observe que as expressoes para num_B_c e den_B.m,den_B.v sao identicas
            // Medias das gaussianas
            HMM->B[i][j][k].m[l][m] = num_B[i][j][k].m[l][m] / num_B[i][j][k].c[l];

            // Variancias das gaussianas
            HMM->B[i][j][k].v[l][m] = num_B[i][j][k].v[l][m] / num_B[i][j][k].c[l];
          }
        }
}

//------------------------------------------------------------------------------

// Funcao que calcula a P(O/M) media para todas as locucoes
double CalculaProbMedia
(
	char **subunits,        // list of phonetic subunits
	struct config trainConfig, // training configuration settings  
	struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
	int nArquivosConv // number of files for convergence testing
)
{
  // Declaracao das variaveis locais
  double **a; // matriz de transicao p/ os modelos das locucoes
  double ***bj; // probabilidade de cada simbolo da locucao, dado o modelo
  struct mistura **b; // matriz de emissao para os modelos das locucoes
  int frase; // conta as locucoes de treinamento
  int register i,j,k,l,m; // contadores
  double ****Nm; // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
  double p_media; // P(O/M) media de todas as locucoes
  struct v_psi vpsi; // var aux p/ determinacao do caminho otimo (Viterbi)
  struct locucao x; // vetores de parametros de treinamento
  FILE *ptrUtterances;
  FILE *ptrTranscriptions;
  char wavFile[512]; // name of file with the utterance
  char transcFile[512]; // name of file with phonetic transcription of the utterance
  int *transcription; // phonetic transcription of the training utterance
  int length; // number of phonetic subunits in phonetic transcription

  // Inicializando verossimilhanca media
  p_media = 0.0;

  // Processando arquivos
  ptrUtterances = fopen(trainConfig.utterancesFiles,"rt");
  ptrTranscriptions = fopen(trainConfig.transcriptionsFiles,"rt");
  
  for (frase=0;frase<nArquivosConv;frase++)
  {
    // Mensagem ao usuario
    printf("%cVerificando convergência do treinamento. Arquivo %d de %d",13,frase+1,nArquivosConv);
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
      b[i] = malloc(sizeof(struct mistura)*HMM->n_estados*length);
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

    // Calculando verossimilhanca
    p_media += Viterbi(a,bj,length,*HMM,x.n_quadros,vpsi);

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

  } // fim loop for frase
  printf("\n");
  
  fclose(ptrUtterances);
  fclose(ptrTranscriptions);
  
  return (p_media/nArquivosConv); // P(O/M) media de todas as locucoes
}

//------------------------------------------------------------------------------
// Funcao que implementa o algoritmo de Viterbi
// Retorna P(O|M)
double Viterbi
(
  double **a,        // matriz de transicao
  double ***bj,      // probabilidade de cada simbolo da locucao, dado o modelo
  int comp,          // numero de subunidades na locucao
  struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  int tamanho,       // comprimento da sequencia de entrada
  struct v_psi vpsi         // estrutura para armazenar o caminho otimo
)
{
  double **a1;  // armazena as probabilidades de transicao em forma de logaritmo
  double ***b1; // armazena as probabilidades de emissao em forma de logaritmo
  int bInf; // verifica se a log-probabilidade de emissao e menos infinito (1) ou nao (0)
  int i,j,k,t;    // contadores
  double aux1,aux2,aux3; // variaveis auxiliares para verificar qual o caminho de maior probabilidade
  double maior;  // variavel que indica a maior probabilidade dentre os caminhos que chegam a um no
  int qual;     // caminho de maior probabilidade que chega ao no

  // Alocando memoria
  a1 = malloc(sizeof(double*)*HMM.n_estados*comp);
  for (i=0;i<(HMM.n_estados*comp);i++)
    a1[i] = malloc(sizeof(double)*(HMM.salto_maximo+1));
  b1 = malloc(sizeof(double**)*HMM.n_par);
  for (i=0;i<HMM.n_par;i++)
    b1[i] = malloc(sizeof(double*)*HMM.n_estados*comp);
  for (i=0;i<HMM.n_par;i++)
    for (j=0;j<(HMM.n_estados*comp);j++)
      b1[i][j] = malloc(sizeof(double)*tamanho);

  // Calculando log10 das probabilidades de transicao e de emissao
  for (i=0;i<HMM.n_estados*comp;i++)
    for (j=0;j<(HMM.salto_maximo+1);j++)
      if (a[i][j] != 0.0)
        a1[i][j] = log10(a[i][j]);
      else
        a1[i][j] = -Inf;

  for (i=0;i<HMM.n_par;i++)
    for (j=0;j<HMM.n_estados*comp;j++)
      for (k=0;k<tamanho;k++)
        if (bj[i][j][k] > Zero)
          b1[i][j][k] = log10(bj[i][j][k]);
        else
          b1[i][j][k] = -Inf;  

  // Inicializacao das matrizes v e psi
  for (i=0;i<HMM.n_estados*comp;i++)
    for (j=0;j<tamanho;j++)
    {
      vpsi.v[i][j] = -Inf;
      vpsi.psi[i][j] = -2;
    }

  // Calculo das probabilidades das transicoes para o estado 0.
  t = 0;
  vpsi.v[0][t] = a1[0][0];
  bInf = 0;
  i = 0;
  while ((i < HMM.n_par) && (bInf == 0))
  {
    if (b1[i][0][t] != -Inf)
      vpsi.v[0][t] += b1[i][0][t];
    else
      bInf = 1;
    i++;
  }
  if (bInf == 1)
    vpsi.v[0][t] = -Inf;
  vpsi.psi[0][t] = -1;

  for (t=1;t<tamanho;t++)
  {
    vpsi.v[0][t] = vpsi.v[0][t-1] + a1[0][0];

    bInf = 0;
    i = 0;
    while ((i < HMM.n_par) && (bInf == 0))
    {
      if (b1[i][0][t] != -Inf)
        vpsi.v[0][t] += b1[i][0][t];
      else
        bInf = 1;
      i++;
    }
    if (bInf == 1)
      vpsi.v[0][t] = -Inf;
    vpsi.psi[0][t] = 0;
  }

  // Calculo das probabilidades para o estado 1.
  for (t=1;t<tamanho;t++)
  {
    // Verificando se veio do estado anterior
    if ((vpsi.v[0][t-1] != -Inf) && (a1[0][1] != -Inf))
    {
      aux1 = vpsi.v[0][t-1] + a1[0][1];
      bInf = 0;
      i = 0;
      while ((i < HMM.n_par) && (bInf == 0))
      {
        if (b1[i][1][t] != -Inf)
          aux1 += b1[i][1][t];
        else
          bInf = 1;
        i++;
      }
      if (bInf == 1)
        aux1 = -Inf;
    }
    else
      aux1 = -Inf;

    // Verificando se veio do mesmo estado
    if ((vpsi.v[1][t-1] != -Inf) && (a1[1][0] != -Inf))
    {
      aux2 = vpsi.v[1][t-1] + a1[1][0];
      bInf = 0;
      i = 0;
      while ((i < HMM.n_par) && (bInf == 0))
      {
        if (b1[i][1][t] != -Inf)
          aux2 += b1[i][1][t];
        else
          bInf = 1;
        i++;
      }    
      if (bInf == 1)
        aux2 = -Inf;
    }
    else
      aux2 = -Inf;

    if (aux1 > aux2) // veio do estado anterior
    {
      maior = aux1;
      qual = 0;
    }
    else // veio do mesmo estado
    {
      maior = aux2;
      qual = 1;
    }
    if (maior != -Inf)
    {
      vpsi.v[1][t] = maior;
      vpsi.psi[1][t] = qual;
    }
  }

  // Calculo das probabilidades para o estado 2.
  for (t=1;t<tamanho;t++)
  {
    // Verificando se veio de 2 estados anteriores
    if ((vpsi.v[0][t-1] != -Inf) && (a1[0][2] != -Inf))
    {
      aux1 = vpsi.v[0][t-1] + a1[0][2];
      bInf = 0;
      i = 0;
      while ((i < HMM.n_par) && (bInf == 0))
      {
        if (b1[i][2][t] != -Inf)
          aux1 += b1[i][2][t];
        else
          bInf = 1;
        i++;
      }
      if (bInf == 1)
        aux1 = -Inf;
    }
    else
      aux1 = -Inf;

    // Verificando se veio do estado anterior
    if ((vpsi.v[1][t-1] != -Inf) && (a1[1][1] != -Inf))
    {
      aux2 = vpsi.v[1][t-1] + a1[1][1];
      bInf = 0;
      i = 0;
      while ((i < HMM.n_par) && (bInf == 0))
      {
        if (b1[i][2][t] != -Inf)
          aux2 += b1[i][2][t];
        else
          bInf = 1;
        i++;
      }
      if (bInf == 1)
        aux2 = -Inf;
    }
    else
      aux2 = -Inf;

    // Verificando se veio do mesmo estado
    if ((vpsi.v[2][t-1] != -Inf) && (a1[2][0] != -Inf))
    {
      aux3 = vpsi.v[2][t-1] + a1[2][0];
      bInf = 0;
      i = 0;
      while ((i < HMM.n_par) && (bInf == 0))
      {
        if (b1[i][2][t] != -Inf)
          aux3 += b1[i][2][t];
        else
          bInf = 1;
        i++;
      }
      if (bInf == 1)
        aux3 = -Inf;
    }
    else
      aux3 = -Inf;

    if (aux1 > aux2) // veio de 2 estados anteriores
    {
      maior = aux1;
      qual = 0;
    }
    else // veio do estado anterior
    {
      maior = aux2;
      qual = 1;
    }
    if (aux3 > maior) // veio do mesmo estado
    {
      maior = aux3;
      qual = 2;
    }
    if (maior != -Inf)
    {
      vpsi.v[2][t] = maior;
      vpsi.psi[2][t] = qual;
    }
  }

  // Calculo das probabilidades para os demais estados.
  if (HMM.n_estados*comp > 3)
  {
    for (i=3;i<HMM.n_estados*comp;i++)
    {
      for (t=2;t<tamanho;t++)
      {
        // Verificando se veio de 2 estados anteriores
        if ((vpsi.v[i-2][t-1] != -Inf) && (a1[i-2][2] != -Inf))
        {
          aux1 = vpsi.v[i-2][t-1] + a1[i-2][2];

          bInf = 0;
          j = 0;
          while ((j < HMM.n_par) && (bInf == 0))
          {
            if (b1[j][i][t] != -Inf)
              aux1 += b1[j][i][t];
            else
              bInf = 1;
            j++;
          }
          if (bInf == 1)
            aux1 = -Inf;
        }
	    else
	      aux1 = -Inf;

        // Verificando se veio do estado anterior
        if ((vpsi.v[i-1][t-1] != -Inf) && (a1[i-1][1] != -Inf))
        {
          aux2 = vpsi.v[i-1][t-1] + a1[i-1][1];
	      bInf = 0;
          j = 0;
          while ((j < HMM.n_par) && (bInf == 0))
          {
            if (b1[j][i][t] != -Inf)
              aux2 += b1[j][i][t];
            else
              bInf = 1;
            j++;
          }
          if (bInf == 1)
            aux2 = -Inf;
        }
        else
          aux2 = -Inf;

        // Verificando se veio do mesmo estado
        if ((vpsi.v[i][t-1] != -Inf) && (a1[i][0] != -Inf))
        {
          aux3 = vpsi.v[i][t-1] + a1[i][0];
          bInf = 0;
          j = 0;
          while ((j < HMM.n_par) && (bInf == 0))
          {
            if (b1[j][i][t] != -Inf)
              aux3 += b1[j][i][t];
            else
              bInf = 1;
            j++;
          }
          if (bInf == 1)
            aux3 = -Inf;
        }
        else
          aux3 = -Inf;

        if (aux1 > aux2) // veio de 2 estados anteriores
        {
          maior = aux1;
	      qual = i - 2;
        }
        else // veio do estado anterior
        {
          maior = aux2;
          qual = i - 1;
        }
        if (aux3 > maior) // veio do mesmo estado
        {
          maior = aux3;
	      qual = i;
        }
        if (maior != -Inf)
        {
          vpsi.v[i][t] = maior;
	      vpsi.psi[i][t] = qual;
        }
      }
    }
  }

  // Desalocando ponteiros
  for (i=0;i<(HMM.n_estados*comp);i++)
    free(a1[i]);
  free(a1);

  for (i=0;i<HMM.n_par;i++)
    for (j=0;j<(HMM.n_estados*comp);j++)
      free(b1[i][j]);
  for (i=0;i<HMM.n_par;i++)
    free(b1[i]);
  free(b1);

  // Retorno de P(O|M)
  return vpsi.v[HMM.n_estados*comp-1][tamanho-1];
}

// Atribui um valor minimo (definido pela variavel MINIMO em HMM.cpp) p/
// as variancias das gaussianas.
int ValoresMinimos
(
  struct modelo_HMM* HMM // modelos HMM das subunidades foneticas
)
{
  int register i,j,k,l,m; // contadores
  int nSubstituicoes = 0; // numero de substituicoes realizadas
  //double soma; // var aux

  for (i=0;i<HMM->n_par;i++)
    for (j=0;j<HMM->n_fones;j++)
      for (k=0;k<HMM->n_estados;k++)
      {
        /*
        // Substituindo coeficientes de ponderacao menores que o limiar
        soma = 0.0;
        for (l=0;l<HMM->n_gaussianas[i];l++)
        {
          if (HMM->B[i][j][k].c[l] < HMM->MINIMO)
          {
            HMM->B[i][j][k].c[l] = HMM->MINIMO;
            nSubstituicoes++;
          }
          soma += HMM->B[i][j][k].c[l];
        }
        for (l=0;l<HMM->n_gaussianas[i];l++)
          HMM->B[i][j][k].c[l] /= soma;
        */
        // Substituindo variancias menores que o limiar
        for (l=0;l<HMM->n_gaussianas[i];l++)
          for (m=0;m<HMM->ordem[i];m++)
            if (HMM->B[i][j][k].v[l][m] < MINIMO)
            {
              HMM->B[i][j][k].v[l][m] = MINIMO;
              nSubstituicoes++;
            }
      }
  return nSubstituicoes;
}

// Verifica se a soma das probabilidades de transicao e igual a 1 e
// se a soma dos coeficientes de ponderacao das gaussianas e igual a 1
// retorna 1 se a verificacao nao encontrou erro
int TestaHMM
(
  struct modelo_HMM* HMM
)
{
  int register i,j,k,l; // contadores
  double soma; // var aux p/ verificar se as somas analisadas dao 1.0
  int sucesso = 1; // verifica se houve ou nao erro

  // Verificando matriz de transicao
  for (i=0;i<HMM->n_fones;i++)
    for (j=0;j<HMM->n_estados;j++)
    {
      soma = 0.0;
      for (k=0;k<HMM->salto_maximo+1;k++)
        soma += HMM->A[i][j][k];
      if ((soma < 0.9999) || (soma > 1.0001))
      {
        printf("Erro na matriz de transição! Fone %d; Estado %d; Transição %d;  Soma: %f",i,j,k,soma);
        sucesso = 0;
      }
    }

  // Verificando coeficientes de ponderacao das gaussianas
  for (i=0;i<HMM->n_par;i++)
    for (j=0;j<HMM->n_fones;j++)
      for (k=0;k<HMM->n_estados;k++)
      {
        soma = 0.0;
        for (l=0;l<HMM->n_gaussianas[i];l++)
          soma += HMM->B[i][j][k].c[l];
        if ((soma < 0.9999) || (soma > 1.0001))
        {
          printf("Erro na matriz de emissão! Parâmetros %d ,Fone %d, Estado %d, Gaussiana %d, Soma: %f",i,j,k,l,soma);
          sucesso = 0;
        }
      }
  return sucesso;
}
