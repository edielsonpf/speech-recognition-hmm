//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "oneStep.h"
#include "constantes.h"
#include "gausFunc.h"
#include "grammar.h"
//---------------------------------------------------------------------------
void oneStep
(
	struct locucao x, // locucao a ser reconhecida
  	char *frase, // frase reconhecida
  	struct configBusca searchConfig, // search algorithm configurations
	struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  	struct bigram *gramatica, // armazena gramatica bigram de palavras
  	int n_termos, // numero de pares da gramatica bigram
  	struct vocab vocabulario // estrutura que armazena o vocabulario
)
{
  int *ativog; // indica se o no gramatico esta ativo ou nao
  int **bp; // tempo decorrido p/ os nos gramaticos
  int ***elapse; // duracao do melhor caminho para o no i  desde que ele entra o LG
  int g; // numero de palavras na locucao
  double **glike; // verossimilhanca dos nos gramaticos
  int *idur=NULL; // duracoes das palavras na frase reconhecida
  int *iword=NULL; // palavras da sentenca reconhecida
  int i,j,k;  // contadores
  double ***like; // verossimilhanca acumulada do melhor caminho para o no i
  double ***log_local; // contribuicao local de cada arco p/ a verossinilhanca total
  double maximo; // var aux p/ determinar o comprimento da locucao
  double max_like; // verossimilhanca maxima de todos os estados em um determinado
                  // instante de tempo
  int nivel; // conta os niveis
  int n_words = vocabulario.n_palavras; // var aux
  int palavra; // conta as palavras do vocabulario
  double penalizado; // verossimilhanca da palavra penalizada pelo modelo de duracao
  double **prob_ant; // probabilidades de transicao dos nos gramaticos
  double *prob_emissao; // probabilidade de emissao p/ cada um dos estados de um fone
  int t; // conta os frames processados
  int **word; // palavra com a maior verossimilhanca
  //int cont;
  //----------------------------------------------------------------------------

  // Alocando memoria para ponteiros
	log_local = malloc(sizeof(double **)*HMM.n_fones);
	for (i=0;i<HMM.n_fones;i++)
		log_local[i] = malloc(sizeof(double *)*HMM.n_estados);
	for (i=0;i<HMM.n_fones;i++)
		for (j=0;j<HMM.n_estados;j++)
			log_local[i][j] = malloc(sizeof(double)*(HMM.salto_maximo+1));

	prob_emissao = malloc(sizeof(double)*HMM.n_estados);

	like = malloc(sizeof(double**)*vocabulario.n_palavras);
	for (i=0;i<vocabulario.n_palavras;i++)
		like[i] = malloc(sizeof(double*)*searchConfig.numberOfLevels);
	for (i=0;i<vocabulario.n_palavras;i++)
		for (j=0;j<searchConfig.numberOfLevels;j++)
			like[i][j] = malloc(sizeof(double)*(HMM.n_estados*vocabulario.Mw[i][0]));

	elapse = malloc(sizeof(int**)*vocabulario.n_palavras);
	for (i=0;i<vocabulario.n_palavras;i++)
		elapse[i] = malloc(sizeof(int*)*searchConfig.numberOfLevels);
	for (i=0;i<vocabulario.n_palavras;i++)
		for (j=0;j<searchConfig.numberOfLevels;j++)
			elapse[i][j] = malloc(sizeof(int)*(HMM.n_estados*vocabulario.Mw[i][0]));

	ativog = malloc(sizeof(int)*(searchConfig.numberOfLevels+1));
	glike = malloc(sizeof(double *)*2);
	word = malloc(sizeof(int *)*(searchConfig.numberOfLevels+1));
	bp = malloc(sizeof(int *)*(searchConfig.numberOfLevels+1));
	prob_ant = malloc(sizeof(double *)*(searchConfig.numberOfLevels+1));
	for (i=0;i<(searchConfig.numberOfLevels+1);i++)
	{
		word[i] = malloc(sizeof(int)*((int)x.n_quadros+1));
		bp[i] = malloc(sizeof(int)*((int)x.n_quadros+1));
		prob_ant[i] = malloc(sizeof(double)*2);
	}
	for (i=0;i<2;i++)
		glike[i] = malloc(sizeof(double)*(searchConfig.numberOfLevels+1));


//------------------------------------------------------------------------------
  // Inicializando variaveis
  for (i=0;i<vocabulario.n_palavras;i++)
  	for (j=0;j<searchConfig.numberOfLevels;j++)
  	  for (k=0;k<(HMM.n_estados*vocabulario.Mw[i][0]);k++)
  	    like[i][j][k] = -Inf;

  for (i=0;i<vocabulario.n_palavras;i++)
  	for (j=0;j<searchConfig.numberOfLevels;j++)
  	  for (k=0;k<(HMM.n_estados*vocabulario.Mw[i][0]);k++)
  	    elapse[i][j][k] = 0;

  for (i=0;i<(searchConfig.numberOfLevels+1);i++)
    for (j=0;j<(x.n_quadros+1);j++)
    {
      word[i][j] = -1;
      bp[i][j] = 0;
    }
  for (i=0;i<(searchConfig.numberOfLevels+1);i++)
    for (j=0;j<2;j++)
    {
      glike[j][i] = -Inf;
      prob_ant[i][j] = -Inf;
    }
  glike[0][0] = 0.0;

  for (i=0;i<(searchConfig.numberOfLevels+1);i++)
    ativog[i] = 0;
  ativog[0] = 1;

  for (i=0;i<HMM.n_fones;i++)
    for (j=0;j<HMM.n_estados;j++)
      for (k=0;k<(HMM.salto_maximo+1);k++)
        log_local[i][j][k] = -Inf;
//******************************************************************************

  // Inicializando prob_ant para o primeiro nivel
  // Prob de autotransicao do primeiro estado do primeiro fonema da palavra
  prob_ant[0][0] = HMM.A[vocabulario.Mw[0][1]][0][0];

	//printf("n_quadros = %d\n",x.n_quadros);		
  // Algoritmo one step
  for (t=0;t<(int)x.n_quadros;t++)
  {
//***************************************************************************
    // Beam Search
    // Inicializando max_like
    max_like = -Inf;
//***************************************************************************
    //printf("t = %d\n",t);

    // Calculando contribuicoes locais de cada arco, p/ cada subunidade fonetica
    for (i=0;i<HMM.n_fones;i++)
    {
      for (j=0;j<HMM.n_estados;j++)
      {
        prob_emissao[j] = 0.0;
        for (k=0;k<HMM.n_par;k++) // LogProb e definido em Gausfunc.h
          prob_emissao[j] += LogProb(HMM.B[k][i][j],HMM.n_gaussianas[k],HMM.ordem[k],x.par[k][t]);
      }
      for (j=0;j<HMM.n_estados;j++)
        for (k=0;k<(HMM.salto_maximo+1);k++)
          if (((j+k) < HMM.n_estados) && (HMM.A[i][j][k] != -Inf) && (prob_emissao[j+k] != -Inf))
            log_local[i][j][k] = HMM.A[i][j][k] + prob_emissao[j+k];
          else
            log_local[i][j][k] = -Inf;
      log_local[i][0][0] = prob_emissao[0]; // probabilidade de emissao no primeiro estado do fone
      log_local[i][2][1] = HMM.A[i][2][1]; // probabilidade de transicao do ultimo estado do fone p/ proximo fone
    } // end for i

    for (nivel=0;nivel<searchConfig.numberOfLevels;nivel++)
    {
      if (nivel==0)
        vocabulario.n_palavras = 1;
      else
        vocabulario.n_palavras = n_words;

      if ((ativog[nivel] == 1) || (nivel==0))
      {

        for (palavra=0;palavra<vocabulario.n_palavras;palavra++)
        {
          // Calculando verossimilhancas
          viterbi(palavra,glike[0][nivel],word[nivel][t],nivel,prob_ant[nivel][0],like[palavra][nivel],HMM,searchConfig,elapse[palavra][nivel],t,x,&max_like,vocabulario,n_termos,gramatica,log_local);

          // Path node merging
          if (like[palavra][nivel][(HMM.n_estados*vocabulario.Mw[palavra][0]-1)] > glike[1][nivel+1]) // evita calculo do modelo de duracao p/ palavra que nao tem chance de ganhar
          {
            if ((searchConfig.useWordDurationModel == 1) && (nivel != 0) && (palavra != 0))
              penalizado = wordDurationModel(like[palavra][nivel][(HMM.n_estados*vocabulario.Mw[palavra][0]-1)],elapse[palavra][nivel][(HMM.n_estados*vocabulario.Mw[palavra][0]-1)],palavra,vocabulario.Ddur,vocabulario.Mdur);
            else
              penalizado = like[palavra][nivel][(HMM.n_estados*vocabulario.Mw[palavra][0]-1)];

            if (penalizado > glike[1][nivel+1])
	          {
              glike[1][nivel+1] = penalizado;
              word[nivel+1][t+1] = palavra;
            }
          } // end if like
        } // end for palavra

        // bp e prob_ant so precisam ser atualizados depois que se conhece a palavra vencedora
        if (glike[1][nivel+1] > -Inf)
        {
          bp[nivel+1][t+1] = elapse[word[nivel+1][t+1]][nivel][(HMM.n_estados*vocabulario.Mw[word[nivel+1][t+1]][0]-1)];
          prob_ant[nivel+1][1] = HMM.A[vocabulario.Mw[word[nivel+1][t+1]][vocabulario.Mw[word[nivel+1][t+1]][0]]][HMM.n_estados-1][HMM.salto_maximo-1];
        }
      } // end if ativog
    } // end for nivel

    // Atualizando vetores glike e prob_ant para o proximo instante de tempo
    for (nivel=0;nivel<(searchConfig.numberOfLevels+1);nivel++)
    {
      glike[0][nivel] = glike[1][nivel];
      prob_ant[nivel][0] = prob_ant[nivel][1];
      if ((glike[0][nivel] > -Inf) && (prob_ant[nivel][0] > -Inf))
        ativog[nivel] = 1;
      glike[1][nivel] = -Inf;
      prob_ant[nivel][1] = -Inf;
    }

    // Beam Search
    if (searchConfig.useBeamSearch == 1)
      beamSearch(&ativog,searchConfig,&like,max_like,HMM,vocabulario,glike[0]);
    //printf("t=%d\n",t);
    //for(cont=0;cont<searchConfig.numberOfLevels;cont++)
    //  printf("ativog[%d] = %d \n",cont,ativog[cont]);
    //printf("\n"); 
       
    /*
    // Salvando glike em arquivo
    if ((arquivo = FormArquivo->AbreArquivo("osgl.dat","a+b")) == NULL)
      exit(1);
    else
    {
      for (nivel=0;nivel<(searchConfig.numberOfLevels+1);nivel++)
        fwrite(&glike[0][nivel], sizeof(double), 1,arquivo);
      fclose(arquivo);
    }

    // Salvando prob_ant em arquivo
    if ((arquivo = FormArquivo->AbreArquivo("ospa.dat","a+b")) == NULL)
      exit(1);
    else
    {
      for (nivel=0;nivel<(searchConfig.numberOfLevels+1);nivel++)
        fwrite(&prob_ant[nivel][0], sizeof(double), 1,arquivo);
      fclose(arquivo);
    }

    // Salvando max_like em arquivo
    if ((arquivo = FormArquivo->AbreArquivo("osml.dat","a+b")) == NULL)
      exit(1);
    else
    {
      fwrite(&max_like, sizeof(double), 1,arquivo);
      fclose(arquivo);
    }

    // Salvando ativog em arquivo
    if ((arquivo = FormArquivo->AbreArquivo("osag.dat","a+b")) == NULL)
      exit(1);
    else
    {
      for (i=0;i<(searchConfig.numberOfLevels+1);i++)
      {
        if (ativog[i] == true)
          iog = 1;
        else
          iog = 0;
        fwrite(&iog, sizeof(int), 1,arquivo);
      }
      fclose(arquivo);
    }
    */
  } // enf for t
  /*
  // Salvando word em arquivo
  if ((arquivo = FormArquivo->AbreArquivo("osw.dat","a+b")) == NULL)
    exit(1);
  else
  {
    for (j=0;j<(x.n_quadros+1);j++)
      for (i=0;i<(searchConfig.numberOfLevels+1);i++)
        fwrite(&word[i][j], sizeof(int), 1,arquivo);
    fclose(arquivo);
  }

  // Salvando bp em arquivo
  if ((arquivo = FormArquivo->AbreArquivo("osbp.dat","a+b")) == NULL)
    exit(1);
  else
  {
    for (j=0;j<(x.n_quadros+1);j++)
      for (i=0;i<(searchConfig.numberOfLevels+1);i++)
        fwrite(&bp[i][j], sizeof(int), 1,arquivo);
    fclose(arquivo);
  }
  */
  // Traceback
  // Determinando numero de palavras na locucao
  maximo = -Inf;
  g = 1;
  for (i=1;i<(searchConfig.numberOfLevels+1);i++)
  {
    if (glike[0][i] > maximo)
    {
      maximo = glike[0][i];
      g = i;
    }
  }
  /*
  // Salvando maximo em arquivo
  if ((arquivo = FormArquivo->AbreArquivo("osmax.dat","a+b")) == NULL)
    exit(1);
  else
  {
    fwrite(&maximo, sizeof(double), 1,arquivo);
    fclose(arquivo);
  }
  */
  if (maximo <= -Inf)
    printf("Limiar de poda muito grande: %f\n",searchConfig.beamSearchThreshold);
  else
  {
    iword = malloc(sizeof(int)*g);
    idur = malloc(sizeof(int)*g);
  }

	i = g; // conta o numero de palavras a serem decodificadas
	t = (int)x.n_quadros;
	k = -1; // conta as palavras na frase
	while (t != 0)
	{
		k++;
		iword[k] = word[i][t];
		idur[k] = bp[i][t];
		i--;
		t -= idur[k];
	}

	
    //~ if (d_final.da_day != d_inic.da_day)
        //~ t_final.ti_hour = (unsigned char)(t_final.ti_hour + 24);

    //~ decorrido = (long)t_final.ti_hour*3600+(long)t_final.ti_min*60+(long)t_final.ti_sec
               //~ -(long)t_inic.ti_hour*3600-(long)t_inic.ti_min*60-(long)t_inic.ti_sec;

    //~ horas = (int)decorrido/3600;
    //~ decorrido -= horas*3600;
    //~ minutos = (int)decorrido/60;
    //~ segundos = (int)decorrido%60;

	// Frase reconhecida
	strcpy(frase,vocabulario.Mws[iword[g-1]]);
	strcat(frase," ");
	for(i=(g-2);i>=0;i--)
	{
		strcat(frase,vocabulario.Mws[iword[i]]);
		strcat(frase," ");
	}
	strcat(frase,"\n");
	
	//~ strcat(frase,";");
	//~ strcat(frase,FloatToStrF(glike[0][g],ffFixed,15,4).c_str());
	//~ strcat(frase,";");
	//~ if (horas < 10)
		//~ strcat(frase,"0");
	//~ strcat(frase,IntToStr(horas).c_str());
	//~ strcat(frase,":");
	//~ if (minutos < 10)
		//~ strcat(frase,"0");
	//~ strcat(frase,IntToStr(minutos).c_str());
	//~ strcat(frase,":");
	//~ if (segundos < 10)
		//~ strcat(frase,"0");
	//~ strcat(frase,IntToStr(segundos).c_str());

	free(iword);
	free(idur);


//******************************************************************************

  // Desalocando ponteiros
  for (i=0;i<HMM.n_fones;i++)
    for (j=0;j<HMM.n_estados;j++)
      free(log_local[i][j]);
  for (i=0;i<HMM.n_fones;i++)
    free(log_local[i]);
  free(log_local);

  free(prob_emissao);

  for (i=0;i<vocabulario.n_palavras;i++)
    for (j=0;j<searchConfig.numberOfLevels;j++)
    {
      free(like[i][j]);
      free(elapse[i][j]);
    }
  for (i=0;i<vocabulario.n_palavras;i++)
  {
    free(like[i]);
    free(elapse[i]);
  }
  free(like);
  free(elapse);

  free(ativog);

  for (i=0;i<(searchConfig.numberOfLevels+1);i++)
  {
 	  free(word[i]);
 	  free(bp[i]);
 	  free(prob_ant[i]);
  }
  free(word);
  free(bp);
  free(prob_ant);

  for (i=0;i<2;i++)
 	free(glike[i]);
  free(glike);

}

//------------------------------------------------------------------------------
void viterbi
(
  int palavra, // palavra sob analise
  double glike, // verossimilhanca do no gramatico
  int word, // armazena as palavras vencedoras a cada nivel e a cada instante de tempo
  int nivel, // nivel da busca
  double prob_ant, // prob de transicao do nivel anterior
  double *like, // verossimilhanca do no intra-palavra
  struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  struct configBusca searchConfig, // configuracoes de busca selecionadas pelo usuario
  int *elapse, // tempo decorrido pelo caminho desde que entrou no no gramatico
  int t, // instante atual
  struct locucao x, // locucao a ser reconhecida
  double *max_like, // verossimilhanca maxima no instante t
  struct vocab vocabulario, // dados do vocabulario
  long n_termos, // numero de pares da gramatica bigram de palavras
  struct bigram *gramatica, // pares da gramatica bigram
  double ***log_local // contribuicao local de cada arco p/ a verossinilhanca total
)
{
  int fone; // conta os fonemas da palavra
  int i,j,node; // contadores
  int *elapse_ant; // armazena os tempos dos caminhos no frame anterior
  double prob_anterior; // var aux
  double **scratch; // vetor inicio temporario p/ salvar like[i] no frame anterior

	scratch = malloc(sizeof(double *)*(int)(HMM.n_estados*vocabulario.Mw[palavra][0]));
	for (i=0;i<(int)(HMM.n_estados*vocabulario.Mw[palavra][0]);i++)
		scratch[i] = malloc(sizeof(double)*HMM.n_estados);
	elapse_ant = malloc(sizeof(int)*(int)(HMM.n_estados*vocabulario.Mw[palavra][0]));

  for (i=0;i<(HMM.n_estados*vocabulario.Mw[palavra][0]);i++)
  {
    elapse_ant[i] = elapse[i];
    for (j=0;j<HMM.n_estados;j++)
      scratch[i][j] = -Inf;
  }

  // Expandindo primeiro fone da palavra
  fone = 0;

  // Expandindo primeiro estado do fone
  node = HMM.n_estados*fone;

  //------------------------------------------------------------------------
  // Gramatica Bigram
  if ((searchConfig.useGrammar == 1) && (nivel > 1) && (palavra > 0) && (word > 0))
    if (Pode(word,palavra,&prob_anterior,n_termos,gramatica) == 0)
      prob_anterior = -Inf;
    else
      prob_anterior = prob_ant;
  else
    prob_anterior = prob_ant;
  //------------------------------------------------------------------------

  if (((glike > -Inf) && (prob_anterior > -Inf)) && (log_local[vocabulario.Mw[palavra][fone+1]][0][0] > -Inf))// veio do nivel anterior
    scratch[node][1] = glike + prob_anterior + log_local[vocabulario.Mw[palavra][fone+1]][0][0];

  if (like[node] > -Inf)
  {
    // ficou no mesmo estado
    if ((HMM.A[vocabulario.Mw[palavra][fone+1]][0][0] > -Inf) && (log_local[vocabulario.Mw[palavra][fone+1]][0][0] > -Inf))
      scratch[node][0] = like[node] + HMM.A[vocabulario.Mw[palavra][fone+1]][0][0] +
                         log_local[vocabulario.Mw[palavra][fone+1]][0][0];
    // transicionou p/ o estado seguinte
    if (log_local[vocabulario.Mw[palavra][fone+1]][0][1] > -Inf)
      scratch[node+1][0] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][0][1];
    //transicionou p/ 2 estados a frente
    if (log_local[vocabulario.Mw[palavra][fone+1]][0][2] > -Inf)
      scratch[node+2][0] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][0][2];
  }

  // Expandindo segundo estado do fone
  node++;
  if (like[node] > -Inf)
  {
    if (log_local[vocabulario.Mw[palavra][fone+1]][1][0] > -Inf)
      scratch[node][1] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][1][0];
    if (log_local[vocabulario.Mw[palavra][fone+1]][1][1] > -Inf)
      scratch[node+1][1] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][1][1];
  }

  // Expandindo terceiro estado do fone
  node++;

  if ((like[node] > -Inf) && (log_local[vocabulario.Mw[palavra][fone+1]][2][0] > -Inf))
    scratch[node][2] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][2][0];


  // Expandindo demais fones da palavra
  for (fone=1;fone<vocabulario.Mw[palavra][0];fone++)
  {
    // Expandindo primeiro estado do fone
    node = HMM.n_estados*fone;

    // veio do fone anterior
    if (like[node-1] > -Inf)
    scratch[node][1] = like[node-1] +
                       log_local[vocabulario.Mw[palavra][fone]][2][1] +
                       log_local[vocabulario.Mw[palavra][fone+1]][0][0];

    if (like[node] > -Inf)
    {
      // ficou no mesmo estado
      scratch[node][0] = like[node] +
                         HMM.A[vocabulario.Mw[palavra][fone+1]][0][0] +
                         log_local[vocabulario.Mw[palavra][fone+1]][0][0];
      // transicionou p/ o estado seguinte
      scratch[node+1][0] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][0][1];
      // transicionou p/ 2 estados posteriores
      scratch[node+2][0] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][0][2];
    }

    // Expandindo segundo estado do fone
    node++;
    if (like[node] > -Inf)
    {
      scratch[node][1] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][1][0];
      scratch[node+1][1] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][1][1];
    }

    // Expandindo terceiro estado do fone
    node++;

    if (like[node] > -Inf)
      scratch[node][2] = like[node] + log_local[vocabulario.Mw[palavra][fone+1]][2][0];
  }

  // Atualizando verossimilhancas e duracoes dos nos
  for (fone=0;fone<vocabulario.Mw[palavra][0];fone++)
  {
    node = HMM.n_estados*fone;

    // Atualizando like e elapse p/ o primeiro estado do fonema
    if (scratch[node][0] > scratch[node][1]) // ficou no mesmo estado
    {
      like[node] = scratch[node][0];
      if (like[node] > *max_like)
        *max_like = like[node];
      elapse[node] = elapse_ant[node] + 1;
    }
    else // veio do estado ou do nivel anterior
    {
      if (scratch[node][1] > -Inf)
      {
        like[node] = scratch[node][1];
        if (like[node] > *max_like)
    	  *max_like = like[node];
    	if (node == 0)
    	  elapse[node] = 1;
    	else
    	  elapse[node] = elapse_ant[node-1]+1;
      }
      else // nao veio de lugar nenhum
      {
        like[node] = -Inf;
        elapse[node] = 0;
      }
    }
    node++;

    // Atualizando like e elapse p/ o segundo estado do fonema
    if (scratch[node][0] > scratch[node][1]) // veio do estado anterior
    {
      like[node] = scratch[node][0];
      if (like[node] > *max_like)
        *max_like = like[node];
      elapse[node] = elapse_ant[node-1]+1;
    }
    else // ficou no mesmo estado
    {
      if (scratch[node][1] > -Inf)
      {
        like[node] = scratch[node][1];
    	if (like[node] > *max_like)
    	  *max_like = like[node];
    	elapse[node] = elapse_ant[node]+1;
      }
      else // nao veio de lugar nenhum
      {
        like[node] = -Inf;
        elapse[node] = 0;
      }
    }
    node++;

    // Atualizando like e elapse p/ o terceiro estado do fonema
    if (scratch[node][0] > scratch[node][1])
    {
      if (scratch[node][0] > scratch[node][2])
      {
        // veio de dois estados anteriores
        like[node] = scratch[node][0];
    	if (like[node] > *max_like)
    	  *max_like = like[node];
    	elapse[node] = elapse_ant[node-2]+1;
      }
      else
      {
        if (scratch[node][2] > -Inf)
    	{
    	  // veio do mesmo estado
    	  like[node] = scratch[node][2];
          if (like[node] > *max_like)
            *max_like = like[node];
          elapse[node] = elapse_ant[node]+1;
        }
      }
    }
    else
    {
      if (scratch[node][1] > scratch[node][2])
      {
        // veio do estado anterior
      	like[node]=scratch[node][1];
      	if (like[node] > *max_like)
      	  *max_like = like[node];
      	elapse[node] = elapse_ant[node-1]+1;
      }
      else
      {
        if (scratch[node][2] > -Inf)
        {
      	  // veio do mesmo estado
      	  like[node] = scratch[node][2];
          if (like[node] > *max_like)
            *max_like = like[node];
          elapse[node] = elapse_ant[node]+1;
        }
        else // nao veio de lugar nenhum
        {
          like[node] = -Inf;
          elapse[node] = 0;
        }
      }
    }
  }

  // Desalocando ponteiros
  for (i=0;i<(int)(HMM.n_estados*vocabulario.Mw[palavra][0]);i++)
    free(scratch[i]);
  free(scratch);
  free(elapse_ant);
}

// Elimina estados com verossimilhanca menor que o limiar
void beamSearch
(
  int **ativog1, // indica se o no gramatico esta ativo ou nao
  struct configBusca searchConfig, // opcoes de configuracao do modo de busca
  double ****like1, // verossimilhanca acumulada do melhor caminho para o no i
  double max_like, // verossimilhanca maxima de todos os estados em um determinado instante de tempo
  struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  struct vocab vocabulario, // estrutura que armazena o vocabulario
  double *glike // verossimilhanca dos nos gramaticos
)
{
  int achou; // var aux p/ desativar niveis gramaticais
  int *ativog; // indica se o no gramatico esta ativo ou nao
  int nivel,palavra,estado; // contadores
  double ***like; // verossimilhanca acumulada do melhor caminho para o no i
  int n_words; // numero de palavras ativas por nivel
  double limiar; // limiar de poda p/ o Beam Search

  like = *like1;
  ativog = *ativog1;

  limiar = max_like - searchConfig.beamSearchThreshold;

  // Desativando estados
  for (nivel=0;nivel<searchConfig.numberOfLevels;nivel++)
  {
    if (nivel == 0)
      n_words = 1;
    else
      n_words = vocabulario.n_palavras;
    for (palavra=0;palavra<n_words;palavra++)
    {
      for (estado=0;estado<HMM.n_estados*vocabulario.Mw[palavra][0];estado++)
        if (like[palavra][nivel][estado] < limiar)
          like[palavra][nivel][estado] = -Inf;
    } // end for (palavra)
  } // end for (nivel)

  // Desativando niveis gramaticais
  for (nivel=0;nivel<searchConfig.numberOfLevels;nivel++)
  {
    if (ativog[nivel] == 1)
    {
      if ((glike[nivel] == -Inf) && (glike[nivel+1] == -Inf))
      {
        if (nivel == 0)
          n_words = 1;
        else
          n_words = vocabulario.n_palavras;

        achou = 0;
        palavra = 0;
        while ((achou == 0) && (palavra < n_words))
        {
          estado = 0;
          while ((achou == 0) && (estado<HMM.n_estados*vocabulario.Mw[palavra][0]))
          {
            if (like[palavra][nivel][estado] > -Inf)
              achou = 1;
            estado++;
          }
          palavra++;
        } // end for (palavra)
        if (achou == 0)
          ativog[nivel] = 0;
      } // end if (glike)   
    } // end if(ativog)

  } // end for (nivel)

  *ativog1 = ativog;
  *like1 = like;
}
