//---------------------------------------------------------------------------
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "estruturas.h"
#include "constantes.h"
#include "arvore.h"
#include "wordDurationModel.h"
#include "gausFunc.h"
#include "herrmanNey.h"

// Implementa o algoritmo de busca p/ lexico organizado em arvore
void HerrmanNey
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
  // Variaveis para gerenciamento das arvores
  //TArvore* Arvore = new TArvore; // classe para manipulacao das arvores lexicas
  int MaxProfundidade; // profundidade maxima observada em todas as arvores lexicas
  int n_arvores; // numero de arvores geradas a partir do arquivo de vocabulario
  int* n_nos; // numero de nos em cada arvore
  struct no ***vocab; // armazena as palavras do vocabulario em arvores

  // Variaveis p/ busca em largura nas arvores
  int ***lista; // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int *profundidades; // indica a profundidade do no mais profundo em cada arvore
  int **quantos; // indica quantos elementos existem em cada nivel de cada uma das listas

  // Variaveis para controle da busca
  int arvores; //conta as arvores
  int *ativog; // indica se o no gramatico esta ativo (1) ou nao (0)
  int ***ativop; // indica quais niveis intra-arvore estao ativos (1) ou inativos (0)
  int **bp; // tempo decorrido p/ os nos gramaticos
  int *bCalculado; // verifica se as contribuicoes locais ja foram (1) ou nao (0) calculadas para um determinado fone
  int g; // numero de palavras na locucao
  double **glike; // verossimilhanca dos nos gramaticos
  int *idur; // duracoes das palavras na frase reconhecida
  int *iword; // palavras da sentenca reconhecida
  double limiar; // limiar de poda para o Beam Search
  double ***log_local; // contribuicao local de cada arco p/ a verossinilhanca total
  double maximo; // var aux p/ determinar o comprimento da locucao
  double max_like; // valor da maxima verossimilhanca em todo o espaco de busca, p/ um determinado instante de tempo
  int nivel; // conta os niveis dentro das arvores lexicas
  int NivelG; // conta os niveis gramaticais
  int nos; // conta os nos das arvores
  int n_trees; // numero de arvores a serem processadas a cada nivel
  double penalizado; // verossimilhanca da palavra penalizada pelo modelo de duracao
  double **prob_ant; // probabilidades de transicao dos nos gramaticos
  int **word; // palavra com a maior verossimilhanca a cada instante de tempo

  // Variaveis de uso geral
  int register i,j,k,t; // contadores

  GeraArvores(&MaxProfundidade,vocabulario.n_palavras,vocabulario.n_fones,vocabulario);
  CarregaArvores(searchConfig.numberOfLevels,&n_arvores,HMM.n_estados,&n_nos,&vocab);
  //ApresentaArvores(n_arvores,n_nos,vocab[1],vocabulario);
  GeraListas(MaxProfundidade,n_arvores,n_nos,vocab[1]);
  CarregaListas(&lista,n_arvores,&profundidades,n_nos,&quantos);
  //Arvore->ApresentaListas(lista,n_arvores,n_niveis,quantos,vocab[1]);

  // Alocando memoria para os ponteiros
  log_local = malloc(sizeof(double **)*HMM.n_fones);
  for (i=0;i<HMM.n_fones;i++)
    log_local[i] = malloc(sizeof(double *)*HMM.n_estados);
  for (i=0;i<HMM.n_fones;i++)
    for (j=0;j<HMM.n_estados;j++)
      log_local[i][j] = malloc(sizeof(double)*(HMM.salto_maximo+1));

  glike = malloc(sizeof(double *)*2);
  word = malloc(sizeof(int *)*(searchConfig.numberOfLevels+1));
  bp = malloc(sizeof(int *)*(searchConfig.numberOfLevels+1));
  prob_ant = malloc(sizeof(double *)*(searchConfig.numberOfLevels+1));

  for (i=0;i<(searchConfig.numberOfLevels+1);i++)
  {
    word[i] = malloc(sizeof(int)*(int)(x.n_quadros+1));
    bp[i] = malloc(sizeof(int)*(int)(x.n_quadros+1));
    prob_ant[i] = malloc(sizeof(double)*2);
  }
  for (i=0;i<2;i++)
    glike[i] = malloc(sizeof(double)*(searchConfig.numberOfLevels+1));
   bCalculado = malloc(sizeof(int)*vocabulario.n_fones);

  ativog = malloc(sizeof(int)*(searchConfig.numberOfLevels+1));

  ativop = malloc(sizeof(int **)*searchConfig.numberOfLevels);
  for (i=0;i<searchConfig.numberOfLevels;i++)
    ativop[i] = malloc(sizeof(int *)*n_arvores);
  for (i=0;i<searchConfig.numberOfLevels;i++)
    for (j=0;j<n_arvores;j++)
      ativop[i][j] = malloc(sizeof(int)*profundidades[j]);

  // Inicializando nos gramaticais
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
  
  // Inicializando log_local
  for (i=0;i<HMM.n_fones;i++)
    for (j=0;j<HMM.n_estados;j++)
      for (k=0;k<(HMM.salto_maximo+1);k++)
        log_local[i][j][k] = -Inf;

  // Inicializando arvores
  for (NivelG=0;NivelG<searchConfig.numberOfLevels;NivelG++)
  {
    if (NivelG == 0) // p/ o primeiro nivel tem-se apenas uma arvore (silencio)
      n_trees = 1;
    else
      n_trees = n_arvores;

    // Inicializando pt -> probabilidade de transicao do ultimo estado de cada no para
    // o no seguinte
    for (i=0;i<n_trees;i++)
      for (j=0;j<n_nos[i];j++)
        vocab[NivelG][i][j].pt = HMM.A[vocab[NivelG][i][j].fone][2][1];

    // Inicializando like e elapse
    for (i=0;i<n_trees;i++)
      for (j=0;j<n_nos[i];j++)
        for (k=0;k<HMM.n_estados;k++)
        {
          vocab[NivelG][i][j].like[k] = -Inf;
          vocab[NivelG][i][j].elapse[k] = 0;
        }

    // Inicializando like_ant e elapse_ant
    for (i=0;i<n_trees;i++)
      for (j=0;j<n_nos[i];j++)
      {
        vocab[NivelG][i][j].like_ant = -Inf;
        vocab[NivelG][i][j].elapse_ant = 0;
      }

    // Inicializando ativo
    for (i=0;i<n_trees;i++)
      for (j=0;j<n_nos[i];j++)
        vocab[NivelG][i][j].ativo = 1;

  } // end for NivelG

  // Inicializando ativop
  for (i=0;i<searchConfig.numberOfLevels;i++)
    for (j=0;j<n_arvores;j++)
    {
      ativop[i][j][0] = 1; // no raiz esta sempre ativo, pois se o no gramatical estiver desativado, ele nao vai ser calculado de qualquer forma
      for (k=1;k<profundidades[j];k++)
        ativop[i][j][k] = 0;
    }

  for (t=0;t<(int)x.n_quadros;t++)
  {
    // Inicializando max_like
    max_like = -Inf;

    // Inicializando bCalculado
    for (i=0;i<HMM.n_fones;i++)
      bCalculado[i] = 0;

    for (NivelG=0;NivelG<searchConfig.numberOfLevels;NivelG++)
    {
      if (ativog[NivelG] == 1)
      {
        // Numero de arvores a serem processadas a cada nivel
        if (NivelG == 0) // no primeiro nivel, apenas o silencio
          n_trees = 1;
        else // nos demais, todas as arvores
          n_trees = n_arvores;

        // Processando as arvores
        for (arvores=0;arvores<n_trees;arvores++)
        {
          // No raiz
          nos = 0;
          nivel = 0;

          // Inicializando prob_ant para o primeiro nivel (prob. de autotransicao
          // do primeiro estado do primeiro fonema da palavra)
          prob_ant[0][0] = HMM.A[vocab[NivelG][arvores][nos].fone][0][0];

          // Calculando contribuicoes locais
          if (bCalculado[vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone] == 0)
          {
            CalcContrLocal(vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone,&log_local,HMM,t,x);
            bCalculado[vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone] = 1;
          }

          // Processando no raiz
          Viterbi(arvores,
                  nos,
                  vocab[NivelG][arvores][nos],
                  0,
                  glike[0][NivelG],
                  HMM,
                  prob_ant[NivelG][0],
                  log_local[vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone],
                  &max_like);

          // Atualizando nos gramaticos
          if (vocab[NivelG][arvores][lista[arvores][nivel][nos]].palavra != -1) // verificando se este no corresponde a uma palavra do vocabulario
          {
            if (vocab[NivelG][arvores][lista[arvores][nivel][nos]].like[2] > glike[1][NivelG+1]) // evita calculo do modelo de duracao p/ palavra que nao tem chance de ganhar
            {
              // Aplicando modelo de duracao de palavras
              if ((searchConfig.useWordDurationModel == 1) && (NivelG != 0) && (arvores != 0))
                penalizado = wordDurationModel(vocab[NivelG][arvores][lista[arvores][nivel][nos]].like[2],
                                                    vocab[NivelG][arvores][lista[arvores][nivel][nos]].elapse[2],
                                                    vocab[NivelG][arvores][lista[arvores][nivel][nos]].palavra,
                                                    vocabulario.Ddur,
                                                    vocabulario.Mdur);
              else
                penalizado = vocab[NivelG][arvores][lista[arvores][nivel][nos]].like[2];

              if (penalizado > glike[1][NivelG+1])
              {
                glike[1][NivelG+1] = penalizado;
                word[NivelG+1][t+1] = vocab[NivelG][arvores][lista[arvores][nivel][nos]].palavra;
                bp[NivelG+1][t+1] = vocab[NivelG][arvores][lista[arvores][nivel][nos]].elapse[2];
                prob_ant[NivelG+1][1] = HMM.A[vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone][2][1];
              }
            } // end if (vocab ... > glike)
          } // end if vocab[nivel][arvores][nos].palavra

          // Atualizando ativop
          if ((vocab[NivelG][arvores][lista[arvores][nivel][nos]].like[2] > -Inf) && (profundidades[arvores] > (nivel+1)))
            ativop[NivelG][arvores][nivel+1] = 1;

          // Processando demais nos
          if (n_nos[arvores] > 1)
          {
            for (nivel=1;nivel<profundidades[arvores];nivel++)
            {
              if (ativop[NivelG][arvores][nivel] == 1)
              {

                for (nos=0;nos<quantos[arvores][nivel];nos++)
                {
                  if (vocab[NivelG][arvores][lista[arvores][nivel][nos]].ativo == 1)
                  {
                    // Calculando contribuicoes locais
                    if (bCalculado[vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone] == 0)
                    {
                      CalcContrLocal(vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone,&log_local,HMM,t,x);
                      bCalculado[vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone] = 1;
                    }

                    // Processando no
                    Viterbi(arvores,
                            lista[arvores][nivel][nos],
                            vocab[NivelG][arvores][lista[arvores][nivel][nos]],
                            vocab[NivelG][arvores][vocab[NivelG][arvores][lista[arvores][nivel][nos]].pai].elapse_ant,
                            vocab[NivelG][arvores][vocab[NivelG][arvores][lista[arvores][nivel][nos]].pai].like_ant,
                            HMM,vocab[NivelG][arvores][vocab[NivelG][arvores][lista[arvores][nivel][nos]].pai].pt,
                            log_local[vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone],
                            &max_like);

                    // Atualizando  gramaticos
                    if (vocab[NivelG][arvores][lista[arvores][nivel][nos]].palavra != -1) // verificando se este no corresponde a uma palavra do vocabulario
                    {
                      if (vocab[NivelG][arvores][lista[arvores][nivel][nos]].like[2] > glike[1][NivelG+1]) // evita calculo do modelo de duracao p/ palavra que nao tem chance de ganhar
                      {
                        // Aplicando modelo de duracao de palavras
                        if ((searchConfig.useWordDurationModel == 1) && (NivelG != 0))
                          penalizado = wordDurationModel(vocab[NivelG][arvores][lista[arvores][nivel][nos]].like[2],
                                                              vocab[NivelG][arvores][lista[arvores][nivel][nos]].elapse[2],
                                                              vocab[NivelG][arvores][lista[arvores][nivel][nos]].palavra,
                                                              vocabulario.Ddur,
                                                              vocabulario.Mdur);
                        else
                          penalizado = vocab[NivelG][arvores][lista[arvores][nivel][nos]].like[2];

                        if (penalizado > glike[1][NivelG+1])
                        {
                          glike[1][NivelG+1] = penalizado;
                          word[NivelG+1][t+1] = vocab[NivelG][arvores][lista[arvores][nivel][nos]].palavra;
                          bp[NivelG+1][t+1] = vocab[NivelG][arvores][lista[arvores][nivel][nos]].elapse[2];
                          prob_ant[NivelG+1][1] = HMM.A[vocab[NivelG][arvores][lista[arvores][nivel][nos]].fone][2][1];
                        }
                      } // end if (vocab ... > glike)
                    } // end if vocab[NivelG][arvores][nos].palavra

                    // Atualizando ativop
                    if ((vocab[NivelG][arvores][lista[arvores][nivel][nos]].like[2] > -Inf) && (profundidades[arvores] > (nivel+1)))
                      ativop[NivelG][arvores][nivel+1] = 1;
                  } // end if (vocab[].ativo == 1)
                } // end for nos
              } // end if ativop[NivelG][arvores][nivel] == 1
            } // end for nivel
          } // end if (n_nos[arvores] > 1)
        } // end for arvores

        // Atualizando like_ant e elapse_ant
        for (i=0;i<n_trees;i++)
          for (j=0;j<n_nos[i];j++)
          {
            vocab[NivelG][i][j].like_ant = vocab[NivelG][i][j].like[2];
            vocab[NivelG][i][j].elapse_ant = vocab[NivelG][i][j].elapse[2];
          }
      } // end if (ativog[NivelG]==1)
    } // end for NivelG

    // Atualizando vetores glike, prob_ant e ativog para o proximo instante de tempo
    for (i=0;i<(searchConfig.numberOfLevels+1);i++)
    {
      glike[0][i] = glike[1][i];
      prob_ant[i][0] = prob_ant[i][1];
      glike[1][i] = -Inf;
      prob_ant[i][1] = -Inf;

      if ((glike[0][i] > -Inf) && (prob_ant[i][0] > -Inf))
        ativog[i] = 1;
    }

    // Podando caminhos menos provaveis (Beam Search)
    if (searchConfig.useBeamSearch == 1)
    {
      limiar = max_like - searchConfig.beamSearchThreshold;
      BeamSearch(&ativop,lista,profundidades,quantos,n_arvores,&vocab,HMM,limiar,glike,&ativog,searchConfig.numberOfLevels,t);
    }
    /*
    // Salvando ativog em arquivo
    if ((arquivo = FormArquivo->AbreArquivo("hnag.dat","a+b")) == NULL)
      exit(1);
    else
    {
      for (i=0;i<(searchConfig.numberOfLevels+1);i++)
      {
        if (ativog[i] == 1)
          iog = 1;
        else
          iog = 0;
        fwrite(&iog, sizeof(int), 1,arquivo);
      }
      fclose(arquivo);
    }

    // Salvando glike em arquivo
    if ((arquivo = FormArquivo->AbreArquivo("hngl.dat","a+b")) == NULL)
      exit(1);
    else
    {
      for (i=0;i<(searchConfig.numberOfLevels+1);i++)
        fwrite(&glike[0][i], sizeof(double), 1,arquivo);
      fclose(arquivo);
    }

    // Salvando prob_ant em arquivo
    if ((arquivo = FormArquivo->AbreArquivo("hnpa.dat","a+b")) == NULL)
      exit(1);
    else
    {
      for (i=0;i<(searchConfig.numberOfLevels+1);i++)
        fwrite(&prob_ant[i][0], sizeof(double), 1,arquivo);
      fclose(arquivo);
    }

    // Salvando max_like em arquivo
    if ((arquivo = FormArquivo->AbreArquivo("hnml.dat","a+b")) == NULL)
      exit(1);
    else
    {
      fwrite(&max_like, sizeof(double), 1,arquivo);
      fclose(arquivo);
    }
    */
  } // end for t
  /*
  // Salvando word em arquivo
  if ((arquivo = FormArquivo->AbreArquivo("hnw.dat","a+b")) == NULL)
    exit(1);
  else
  {
    for (j=0;j<(x.n_quadros+1);j++)
      for (i=0;i<(searchConfig.numberOfLevels+1);i++)
        fwrite(&word[i][j], sizeof(int), 1,arquivo);
    fclose(arquivo);
  }

  // Salvando bp em arquivo
  if ((arquivo = FormArquivo->AbreArquivo("hnbp.dat","a+b")) == NULL)
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
  if ((arquivo = FormArquivo->AbreArquivo("hnmax.dat","a+b")) == NULL)
    exit(1);
  else
  {
    fwrite(&maximo, sizeof(double), 1,arquivo);
    fclose(arquivo);
  }
  */
  if (maximo <= -Inf)
    printf("Limiar: %f. Limiar de poda muito grande!\n",searchConfig.beamSearchThreshold);
  else
  {
    iword=malloc(sizeof(int)*g);
    idur=malloc(sizeof(int)*g);

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

    // Parando cronometragem
    //gettime(&t_final);
    //getdate(&d_final);

    // Calculando tempo decorrido
    //if (d_final.da_day != d_inic.da_day)
    //    t_final.ti_hour = (unsigned char)(t_final.ti_hour + 24);

    //decorrido = (long)t_final.ti_hour*3600+(long)t_final.ti_min*60+(long)t_final.ti_sec
    //           -(long)t_inic.ti_hour*3600-(long)t_inic.ti_min*60-(long)t_inic.ti_sec;

    //horas = (int)decorrido/3600;
    //decorrido -= horas*3600;
    //minutos = (int)decorrido/60;
    //segundos = (int)decorrido%60;

    // Frase reconhecida
    strcpy(frase,vocabulario.Mws[iword[g-1]]);
    strcat(frase," ");
    for(i=(g-2);i>=0;i--)
    {
      strcat(frase,vocabulario.Mws[iword[i]]);
      strcat(frase," ");
    }
    /*
    strcat(frase,";");
    strcat(frase,FloatToStrF(glike[0][g],ffFixed,15,4).c_str());
    strcat(frase,";");
    if (horas < 10)
      strcat(frase,"0");
    strcat(frase,IntToStr(horas).c_str());
    strcat(frase,":");
    if (minutos < 10)
      strcat(frase,"0");
    strcat(frase,IntToStr(minutos).c_str());
    strcat(frase,":");
    if (segundos < 10)
      strcat(frase,"0");
    strcat(frase,IntToStr(segundos).c_str());
    */
    free(iword);
    free(idur);
  }

  // Desalocando ponteiros
  for (i=0;i<searchConfig.numberOfLevels;i++)
    for (j=0;j<n_arvores;j++)
      free(ativop[i][j]);
  for (i=0;i<searchConfig.numberOfLevels;i++)
    free(ativop[i]);
  free(ativop);

  DesalocaArvores(searchConfig.numberOfLevels,n_arvores,&n_nos,&vocab);
  DesalocaListas(&lista,n_arvores,&profundidades,&quantos);
  
  for (i=0;i<HMM.n_fones;i++)
    for (j=0;j<HMM.n_estados;j++)
      free(log_local[i][j]);
  for (i=0;i<HMM.n_fones;i++)
    free(log_local[i]);
  free(log_local);

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

  free(bCalculado);
}

//------------------------------------------------------------------------------
// Implementa o algoritmo de Viterbi para uso com o algoritmo Herrman-Ney
void Viterbi
(
  int arvore, // arvore a ser processada
  int no, // no a ser processado
  struct no vocab, // dados do no a ser processado
  int dur_ant, // duracao do caminho otimo ate o ultimo estado do fone anterior
  double log_ant, // verossimilhanca acumulada ate o ultimo estado do fone anterior
  struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  double pt, // probabilidade de transicao do no anterior para este no
  double **log_local, // contribuicao local de cada arco p/ a verossimilhanca total
  double *max_like1 // max veross. neste instante de tempo
)
{
  double max_like; // max veross. neste instante de tempo
  int *elapse_ant; // armazena os tempos dos caminhos no frame anterior
  int register i,j; // contadores
  double **scratch; // armazena as verossimilhancas de cada arco (scratch[destino][origem])
                    // origem = 0 -> veio de dois estados anteriores
                    // origem = 1 -> veio do estado anterior
                    // origem = 2 -> veio do mesmo estado

  max_like = *max_like1;

  scratch = malloc(sizeof(double *)*HMM.n_estados);
  for (i=0;i<HMM.n_estados;i++)
    scratch[i] = malloc(sizeof(double)*HMM.n_estados);

  elapse_ant = malloc(sizeof(int)*HMM.n_estados);

  // Inicializando variaveis
  for (i=0;i<HMM.n_estados;i++)
  {
    elapse_ant[i] = vocab.elapse[i];
    for (j=0;j<HMM.n_estados;j++)
      scratch[i][j] = -Inf;
  }

  // Processando primeiro estado do fone
  if (log_ant != -Inf) // veio do fone anterior
    scratch[0][1] = log_ant + pt + log_local[0][0];
  if (vocab.like[0] != -Inf) // veio do mesmo estado
    scratch[0][0] = vocab.like[0] + HMM.A[vocab.fone][0][0] + log_local[0][0];

  // Processando segundo estado do fone
  if (vocab.like[0] != -Inf) // veio do estado anterior
    scratch[1][1] = vocab.like[0] + log_local[0][1];
  if (vocab.like[1] != -Inf) // veio do mesmo estado
    scratch[1][0] = vocab.like[1] + log_local[1][0];

  // Processando terceiro estado do fone
  if (vocab.like[0] != -Inf) // veio de dois estados anteriores
    scratch[2][2] = vocab.like[0] + log_local[0][2];
  if (vocab.like[1] != -Inf) // veio do estado anterior
    scratch[2][1] = vocab.like[1] + log_local[1][1];
  if (vocab.like[2] != -Inf) // veio do mesmo estado
    scratch[2][0] = vocab.like[2] + log_local[2][0];

  // Atualizando verossimilhancas e duracoes dos estados
  // Primeiro estado
  if (scratch[0][1] > scratch[0][0]) // veio do fone anterior
  {
    vocab.like[0] = scratch[0][1];
    if (vocab.like[0] != -Inf)
      vocab.elapse[0] = dur_ant + 1;
  }
  else // ficou no mesmo estado
  {
    vocab.like[0] = scratch[0][0];
    if (vocab.like[0] != -Inf)
      vocab.elapse[0] = elapse_ant[0] + 1;
  }
  // Segundo estado
  if (scratch[1][1] > scratch[1][0]) // veio do estado anterior
  {
    vocab.like[1] = scratch[1][1];
    vocab.elapse[1] = elapse_ant[0] + 1;
  }
  else // ficou no mesmo estado
  {
    if (scratch[1][0] > -Inf)
    {
      vocab.like[1] = scratch[1][0];
      vocab.elapse[1] = elapse_ant[1] + 1;
    }
  }
  // Terceiro estado
  if (scratch[2][2] > scratch[2][1])
  {
    if (scratch[2][2] > scratch[2][0]) // veio de dois estados anteriores
    {
      vocab.like[2] = scratch[2][2];
      vocab.elapse[2] = elapse_ant[0] + 1;
    }
    else // ficou no mesmo estado
    {
      if (scratch[2][0] > -Inf)
      {
        vocab.like[2] = scratch[2][0];
        vocab.elapse[2] = elapse_ant[2] + 1;
      }
    }
  }
  else
  {
    if (scratch[2][1] > scratch[2][0]) // veio do estado anterior
    {
      vocab.like[2] = scratch[2][1];
      vocab.elapse[2] = elapse_ant[1] + 1;
    }
    else // ficou no mesmo estado
    {
      if (scratch[2][0] > -Inf)
      {
        vocab.like[2] = scratch[2][0];
        vocab.elapse[2] = elapse_ant[2] + 1;
      }
    }
  }

  for (i=0;i<HMM.n_estados;i++)
    if (vocab.like[i] > max_like)
      max_like = vocab.like[i];

  // Desalocando ponteiros
  for (i=0;i<HMM.n_estados;i++)
    free(scratch[i]);
  free(scratch);

  free(elapse_ant);

  *max_like1 = max_like;
}

// Implementa a poda de caminhos via Beam Search
void BeamSearch
(
  int ****ativop1, // indica quais niveis intra-arvore estao ativos
  int ***lista, // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int *profundidades, // indica a profundidade do no mais profundo em cada arvore
  int **quantos, // indica quantos elementos existem em cada nivel de cada uma das listas
  int n_arvores, // numero de arvores
  struct no ****vocab1, // armazena as palavras do vocabulario em arvores
  struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  double limiar, // limiar de poda
  double **glike, // verossimilhanca dos nos gramaticos
  int **ativog1, // indica se o no gramatico esta ativo ou nao
  int numberOfLevels, // numero maximo de palavras da locucao
  int t
)
{
  int ***ativop; // indica quais niveis intra-arvore estao ativos
  struct no ***vocab; // armazena as palavras do vocabulario em arvores
  int *ativog; // indica se o no gramatico esta ativo ou nao
  int achou; // var aux p/ desligar niveis gramaticais
  int *desligar; // verifica se os 3 estados do no foram podados, e entao retira o no da busca
  int register nivelg,arvore,profundidade,noh,estado; // contadores
  int n_trees; // numero de arvores ativas
  FILE *arquivo;
  
  ativop = *ativop1;
  ativog = *ativog1;
  vocab = *vocab1;

  desligar = malloc(sizeof(int)*HMM.n_estados);

  // Eliminando estados com verossimilhanca abaixo do limiar e desativando nos
  for (nivelg=0;nivelg<numberOfLevels;nivelg++)
  {
    if (ativog[nivelg] == 1)
    {
      if (nivelg == 0)
        n_trees = 1;
      else
        n_trees = n_arvores;

      for (arvore=0;arvore<n_trees;arvore++)
        for (profundidade=0;profundidade<profundidades[arvore];profundidade++)
          if (ativop[nivelg][arvore][profundidade] == 1)
            for (noh=0;noh<quantos[arvore][profundidade];noh++)
            {
              for (estado=0;estado<HMM.n_estados;estado++)
                if (vocab[nivelg][arvore][lista[arvore][profundidade][noh]].like[estado] < limiar)
                {
                  vocab[nivelg][arvore][lista[arvore][profundidade][noh]].like[estado] = -Inf;
                  if (estado == (HMM.n_estados-1))
                    vocab[nivelg][arvore][lista[arvore][profundidade][noh]].like_ant = -Inf;
                  desligar[estado] = 1;
                }
                else
                  desligar[estado] = 0;

                // Desativando o no
                if ((desligar[0] == 1) && (desligar[1] == 1) && (desligar[2] == 1))
                {
                  if (profundidade != 0)
                  {
                    if (vocab[nivelg][arvore][vocab[nivelg][arvore][lista[arvore][profundidade][noh]].pai].like[2] == -Inf)
                      vocab[nivelg][arvore][lista[arvore][profundidade][noh]].ativo = 0;
                    else
                      vocab[nivelg][arvore][lista[arvore][profundidade][noh]].ativo = 1;
                  }
                  else
                    vocab[nivelg][arvore][lista[arvore][profundidade][noh]].ativo = 0;
                } // end if (desligar)
                else
                  vocab[nivelg][arvore][lista[arvore][profundidade][noh]].ativo = 1;
                  
            } // end for noh
    } // end if (ativog)
  } // end for (nivelg)

  // Desativando niveis dentro das arvores
  for (nivelg=0;nivelg<numberOfLevels;nivelg++)
  {
    if (ativog[nivelg] == 1)
    {
      if (nivelg == 0)
        n_trees = 1;
      else
        n_trees = n_arvores;

      for (arvore=0;arvore<n_trees;arvore++)
        for (profundidade=0;profundidade<profundidades[arvore];profundidade++)
          if (ativop[nivelg][arvore][profundidade] == 1)
          {
            noh=0;
            while ((noh<quantos[arvore][profundidade]) && (vocab[nivelg][arvore][lista[arvore][profundidade][noh]].ativo == 0))
              noh++;
            if (noh == quantos[arvore][profundidade])
              ativop[nivelg][arvore][profundidade] = 0;
          }
    } // end if(ativog)
  } // end for (nivelg)

  // Desativando niveis gramaticais
  for (nivelg=0;nivelg<numberOfLevels;nivelg++)
  {
    if (ativog[nivelg] == 1)
    {
      if ((glike[0][nivelg] == -Inf) && (glike[0][nivelg+1] == -Inf))
      {
        if (nivelg == 0)
          n_trees = 1;
        else
          n_trees = n_arvores;

        achou = 0;
        arvore = 0;
        while ((achou == 0) && (arvore < n_trees))
        {
          profundidade = 0;
          while ((achou == 0) && (profundidade < profundidades[arvore]))
          {
            if (ativop[nivelg][arvore][profundidade] == 1)
              achou = 1;
            profundidade++;
          }
          arvore++;
        } // end while (achou)
        if (achou == 0)
          ativog[nivelg] = 0;
      } // end if(glike)
    } // end if (ativog)
  } // end for (nivelg)

  // Reativando nos raizes dos niveis gramaticais que nao foram podados
  for (nivelg=0;nivelg<numberOfLevels;nivelg++)
    if (ativog[nivelg] == 1)
    {
      if (nivelg == 0)
        n_trees = 1;
      else
        n_trees = n_arvores;

      for (arvore=0;arvore<n_trees;arvore++)
      {
        ativop[nivelg][arvore][0] = 1;
        vocab[nivelg][arvore][0].ativo = 1;
      }
    }

  free(desligar);

  // Salvando numero de nos ativos em arquivo
  int nAtivos = 0;
  for (nivelg=0;nivelg<numberOfLevels;nivelg++)
    if (ativog[nivelg] == 1)
    {
      if (nivelg == 0)
        n_trees = 1;
      else
        n_trees = n_arvores;

      for (arvore=0;arvore<n_trees;arvore++)
        for (profundidade=0;profundidade<profundidades[arvore];profundidade++)
          if (ativop[nivelg][arvore][profundidade] == 1)
            for (noh=0;noh<quantos[arvore][profundidade];noh++)
              for (estado=0;estado<HMM.n_estados;estado++)
                if (vocab[nivelg][arvore][lista[arvore][profundidade][noh]].ativo == 1)
                  nAtivos++;
    }
  arquivo = fopen("hnnat.dat","a+b");  
  if (arquivo == NULL)
  {
  	puts("Error opening file hnnat.dat. Verify permissions!\n");
    exit(1);
  }
  else
  {
    fwrite(&nAtivos, sizeof(int), 1,arquivo);
    fclose(arquivo);
  }

  *ativop1 = ativop;
  *ativog1 = ativog;
  *vocab1 = vocab;
}


