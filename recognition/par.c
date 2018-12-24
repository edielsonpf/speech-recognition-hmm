//---------------------------------------------------------------------------
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <wav.h>

#include "par.h"
#include "preProc.h"
#include "fft.h"
#include "mel.h"
#include "energia.h"
#include "perfil.h"
#include "lpc.h"
#include "pca.h"

#define ORDEM_MEL 12 // ordem dos parametros mel-cepstrais (tb delta e delta-delta)
#define MAXPAR 20 // numero maximo de parametros permitidos
//---------------------------------------------------------------------------
// Funcao que calcula os parametros de um sinal de entrada 'x'.
// Calcula os coeficientes ate ordem 'ordem' e armazena os resultados em 'par1'.
void calcPar
(
  char* nome_arquivo, // nome do arquivo WAV a ser parametrizado
	struct modelo_HMM HMM, // acoustic models (contain information about which parameters to use, dimension of parameters, whether use or not PCA, etc.)
  int* Nwindow1,      // numero de quadros com o qual foi analisado o sinal
	double**** par1    // sinal parametrizado
)
{
  // Declaracao de variaveis locais
  int aux; // variavel auxiliar
  short bps; // bits por amostra
  double CP=0.95; // coeficiente de pre-enfase
  int nQuadrosAdjacentes[MAXPAR]; // numero de quadros adjacentes a serem considerados no calculo de parameros delta
                      // delta1[0]: mel-cepstrais
                      // delta1[1]: delta mel-cepstrais
                      // delta1[2]: delta-delta mel-cepstrais
                      // delta1[3]: log enregia normalizada + delta + delta-delta
                      // delta1[4]: mel + delta + delta-delta
                      // delta1[5]: perfil de energia
                      // delta1[6]: delta perfil de energia
                      // delta1[7]: delta-delta perfil de energia
                      // delta1[8]: lpc
  double e_max; // energia do quadro de maior energia do sinal
  int fs; // frequencia de amostragem do sinal
  int i,j,k; // contadores
  int Lframe; // tamanho de quadro de analise 10 ms
  int Lwindow; // tamanho da janela de analise 20 ms
  int M; // M=9 ==> FFT de 512 pontos; M=10 ==> FFT de 1024 pontos
  double *mod_fft=NULL; // modulo da FFT
  int N; // numero de pontos da FFT (N=2**M)
  //double* MediaMel; // vetor media dos parametros mel-cepstrais
  int Nam; // num de amostras do sinal
  int Nwindow; //numero de quadros com o qual foi analisado o sinal
  //int* ordem;  // ordem dos vetores de cada um dos parametros
  double ***par=NULL; // sinal parametrizado
  int posicao[MAXPAR]; // ordem em que os parametros devem ser organizados no vetor de saida
               // posicao[0]: mel-cepstrais
               // posicao[1]: delta mel-cepstrais
               // posicao[2]: delta-delta mel-cepstrais
               // posicao[3]: log energia normalizada + delta + delta-delta
               // posicao[4]: mel + delta + delta-delta
               // posicao[5]: perfil de energia
               // posicao[6]: delta perfil de energia
               // posicao[7]: delta-delta perfil de energia
               // posicao[8]: lpc
               // Exemplo: posicao[0] = 2 -> o parametro mel-cepstral e o terceiro parametro
               // no vetor par
  int quais[MAXPAR]; // verifica quais parametros devem ser calculados (0) nÃ£o, (1) sim
                 // quais[0]: mel-cepstrais
                 // quais[1]: delta mel-cepstrais
                 // quais[2]: delta-delta mel-cepstrais
                 // quais[3]: log enrgia normalizada + delta + delta-delta
                 // quais[4]: mel + delta + delta-delta
                 // quais[5]: perfil de energia
                 // quais[6]: delta perfil de energia
                 // quais[7]: delta-delta perfil de energia
                 // quais[8]: lpc
  int tamanho; // numero de bytes do sinal
  char *x=NULL; // ponteiro para o sinal a ser analisado
  short *x1=NULL; // sinal armazenado em formato short
  double *x2=NULL; // sinal armazenado em formato double
  double *xw=NULL; // ams da janela de trabalho (max 60 mseg)
  double **maux=NULL,**daux=NULL,**ddaux=NULL; // var aux p/ armazenar as componentes do vetor de parametros mdd (mel, delta e delta-delta)
  
	for (i=0;i<MAXPAR;i++)
		nQuadrosAdjacentes[i]=0;
	
	// Verificando quais parametros devem ser calculados
  for (i=0;i<MAXPAR;i++)
  {
    quais[i] = 0;
    posicao[i] = -1;
  }

  for (i=0;i<HMM.n_par;i++)
    if (strcmp(HMM.parametros[i],"mel") == 0)
    {
      quais[0] = 1;
      posicao[0] = i;
      nQuadrosAdjacentes[0] = HMM.delta[i];
    }
    else if (strcmp(HMM.parametros[i],"delta-mel") == 0)
    {
      quais[1] = 1;
      posicao[1] = i;
      nQuadrosAdjacentes[1] = HMM.delta[i];
    }
    else if (strcmp(HMM.parametros[i],"delta-delta-mel") == 0)
    {
      quais[2] = 1;
      posicao[2] = i;
      nQuadrosAdjacentes[2] = HMM.delta[i];
    }
    else if (strcmp(HMM.parametros[i],"log-energia") == 0)
    {
      quais[3] = 1;
      posicao[3] = i;
      nQuadrosAdjacentes[3] = HMM.delta[i];
    }
    else if (strcmp(HMM.parametros[i],"mdd") == 0)
    {
      HMM.ordem[0] = 36;
      quais[4] = 1;
      posicao[4] = i;
      nQuadrosAdjacentes[4] = HMM.delta[i];
    }
    else if (strcmp(HMM.parametros[i],"perfil-energia") == 0)
    {
      quais[5] = 1;
      posicao[5] = i;
      nQuadrosAdjacentes[5] = HMM.delta[i];
    }
    else if (strcmp(HMM.parametros[i],"delta-perfil-energia") == 0)
    {
      quais[6] = 1;
      posicao[6] = i;
      nQuadrosAdjacentes[6] = HMM.delta[i];
    }
    else if (strcmp(HMM.parametros[i],"delta-delta-perfil-energia") == 0)
    {
      quais[7] = 1;
      posicao[7] = i;
      nQuadrosAdjacentes[7] = HMM.delta[i];
    }
    else if (strcmp(HMM.parametros[i],"lpc") == 0)
    {
      quais[8] = 1;
      posicao[8] = i;
      nQuadrosAdjacentes[8] = HMM.delta[i];
    }

  // Abrindo e lendo sinal a ser parametrizado
  // esta funcao tambem aloca o ponteiro x
  leWav(&bps,&fs,nome_arquivo,&tamanho,&x);

  // Determinando tamanho da janela, quadro de analise
  Lframe = fs / 100;
  Lwindow = 2*Lframe;

  // Definindo numero de pontos da FFT
  switch(fs)
  {
    case 8000:
      M=9;
      break;
    case 11025:
      M=9;
      break;
    default:
      M=10;
  }
  N = 1.0;
  for (i=0;i<M;i++)
    N *= 2; // N = 2^M

  // Numero de amostras do sinal
  Nam = tamanho/sizeof(short int);

  // Determina o numero total de janelas.
  Nwindow = (int)((Nam-(Lframe)-1) / (Lframe));
  if((Nam%Lframe)!=0)
    Nwindow++;

  // Alocando memoria
	x2 = malloc(sizeof(double)*Nam);
	xw = malloc(sizeof(double)*Lwindow);

	par = malloc(sizeof(double **)*HMM.n_par);
	for (i=0;i<HMM.n_par;i++)
		par[i] = malloc(sizeof(double *)*Nwindow);
	for (i=0;i<HMM.n_par;i++)
		for (j=0;j<Nwindow;j++)
			par[i][j] = malloc(sizeof(double)*HMM.ordem[i]);

	if (quais[4] == 1)
	{
		maux = malloc(sizeof(double *)*Nwindow);
		daux = malloc(sizeof(double *)*Nwindow);
		ddaux = malloc(sizeof(double *)*Nwindow);
		for (i=0;i<Nwindow;i++)
		{
			maux[i] = malloc(sizeof(double)*ORDEM_MEL);
			daux[i] = malloc(sizeof(double)*ORDEM_MEL);
			ddaux[i] = malloc(sizeof(double)*ORDEM_MEL);
		}
	}

	mod_fft = malloc(sizeof(double)*N);

	//MediaMel = new double[ORDEM_MEL];

  // Inicializando MediaMel
  //for (i=0;i<ORDEM_MEL;i++)
  //  MediaMel[i] = 0.0;

  // Armazenando sinal em formato short
  x1 = (short *)x;

  // Armazenando sinal em formato double
  for (i=0;i<Nam;i++)
    x2[i] = (double)x1[i];

  // Retirando nivel DC do sinal
  removeDC(Nam,x2);

  // Filtragem de pre-enfase
  preEnfase(CP,Nam,x2);

  //----------------------------------------------------------------------------
  // Calculando parametros
  for(k=0;k<(Nwindow-1);k++)
  {
    // Selecionando trecho a ser analisado
    for(i=0;i<Lwindow;i++)
      xw[i] = x2[k*Lframe+i];

    // Janelamento
    hamming(Lwindow,xw);

    // Calculando FFT
    fft(xw,Lwindow,fs,M,N,mod_fft);

    // Calculando parametros mel cepstrais
    if (quais[0] == 1)
    {
      calcMel(fs,mod_fft,N,ORDEM_MEL,par[posicao[0]][k]);
      //for (i=0;i<ORDEM_MEL;i++)
      //  MediaMel[i] += par[posicao[0]][k][i];
    }
    // Calculando parte mel dos parametros mdd
    else if (quais[4] == 1)
    {
      calcMel(fs,mod_fft,N,ORDEM_MEL,maux[k]);
      for (i=0;i<ORDEM_MEL;i++)
      {
        par[0][k][i] = maux[k][i];
        //MediaMel[i] += par[0][k][i];
      }
    }
    // Calculando parametros log-energia normalizada
    if (quais[3] == 1)
      par[posicao[3]][k][0] = calcEnergia(xw,Lwindow);

    // Calculando parametros perfil-energia
    if (quais[5] == 1)
      calcPerfilEnergia(fs,mod_fft,N,HMM.ordem[posicao[5]],par[posicao[5]][k]);

    // Calculando parametros LPC (Alteracao em 06/10 por Alexander Coelho)
    if (quais[8] == 1)
      calcLPC(xw,Lwindow,HMM.ordem[posicao[8]],par[posicao[8]][k]);
  }

  //----------------------------------------------------------------------------
  // Ultima janela
  // Selecionando trecho a ser analisado
  aux = Nam - Lwindow;
  for(i=0;i<Lwindow;i++)
    xw[i] = x2[aux+i-2];

  // Janelamento
  hamming(Lwindow,xw);

  // Calculando FFT
  fft(xw,Lwindow,fs,M,N,mod_fft);

  // Calculo dos coeficientes Mel-cepstrais
  if (quais[0] == 1)
  {
    calcMel(fs,mod_fft,N,ORDEM_MEL,par[posicao[0]][k]);
    //for (i=0;i<ORDEM_MEL;i++)
    //  MediaMel[i] += par[posicao[0]][k][i];
  }

  // Calculando parte mel dos parametros mdd
  else if (quais[4] == 1)
  {
    calcMel(fs,mod_fft,N,ORDEM_MEL,maux[k]);
    for (i=0;i<ORDEM_MEL;i++)
    {
      par[0][k][i] = maux[k][i];
      //MediaMel[i] += par[0][k][i];
    }
  }

  // Calculando parametros log-energia normalizada
  if (quais[3] == 1)
    par[posicao[3]][k][0] = calcEnergia(xw,Lwindow);

  // Calculando parametros perfil-energia
  if (quais[5] == 1)
    calcPerfilEnergia(fs,mod_fft,N,HMM.ordem[posicao[5]],par[posicao[5]][k]);

  //----------------------------------------------------------------------------
  /*
  // Retirando media dos parametros mel-cepstrais
  for (i=0;i<ORDEM_MEL;i++)
    MediaMel[i] /= Nwindow;

  if (quais[0] == true)
  {
    for (i=0;i<Nwindow;i++)
      for (j=0;j<ORDEM_MEL;j++)
        par[posicao[0]][i][j] -= MediaMel[j];
  }
  else if (quais[4] == true)
  {
    for (i=0;i<Nwindow;i++)
      for (j=0;j<ORDEM_MEL;j++)
        par[0][i][j] -= MediaMel[j];
  }
  */
  //----------------------------------------------------------------------------
  // Processamentos finais para o calculo dos parametros log-energia normalizada

  // Normalizando energia
  if (quais[3] == 1)
  {
    e_max = 0.0;
    for (i=0;i<Nwindow;i++)
      if (par[posicao[3]][i][0] > e_max)
        e_max = par[posicao[3]][i][0];
    for (i=0;i<Nwindow;i++)
      par[posicao[3]][i][0] /= e_max;

    // Calculando log10 da energia
    for (i=0;i<Nwindow;i++)
      par[posicao[3]][i][0] = log10(par[posicao[3]][i][0]);
  }

  // Calculando parametros LPC (Alteracao em 06/10 por Alexander Coelho)
  if (quais[8] == 1)
    calcLPC(xw,Lwindow,HMM.ordem[posicao[8]],par[posicao[8]][k]);
  //----------------------------------------------------------------------------
  // Calculando parametros delta mel
  if (quais[1] == 1)
    calcDelta(nQuadrosAdjacentes[1],Nwindow,HMM.ordem[0],par[posicao[0]],par[posicao[1]]);

  // Calculando parte delta mel dos parametros mdd
  else if (quais[4] == 1)
  {
    calcDelta(nQuadrosAdjacentes[4],Nwindow,ORDEM_MEL,maux,daux);
    for (k=0;k<Nwindow;k++)
      for (i=0;i<ORDEM_MEL;i++)
        par[0][k][i+ORDEM_MEL] = daux[k][i];
  }

  // Calculando parametros delta-perfil-energia
  if (quais[6] == 1)
    calcDelta(nQuadrosAdjacentes[6],Nwindow,HMM.ordem[posicao[6]],par[posicao[5]],par[posicao[6]]);

  //----------------------------------------------------------------------------
  // Calculando parametros delta-delta mel
  if (quais[2] == 1)
    calcDelta(nQuadrosAdjacentes[2],Nwindow,HMM.ordem[1],par[posicao[1]],par[posicao[2]]);

  // Calculando parte delta mel dos parametros mdd
  else if (quais[4] == 1)
  {
    calcDelta(nQuadrosAdjacentes[4],Nwindow,ORDEM_MEL,daux,ddaux);
    for (k=0;k<Nwindow;k++)
      for (i=0;i<ORDEM_MEL;i++)
        par[0][k][i+24] = ddaux[k][i];
  }
  // Calculando parametros delta e delta-delta para os parametros log-energia
  if (quais[3] == 1)
    CalculaDeltaEnergia(par[posicao[3]],Nwindow,nQuadrosAdjacentes[3]);

  // Calculando parametros delta-perfil-energia
  if (quais[7] == 1)
    calcDelta(nQuadrosAdjacentes[7],Nwindow,HMM.ordem[posicao[7]],par[posicao[6]],par[posicao[7]]);

//----------------------------------------------------------------------------
  // Reducao da ordem do vetor para parametros mdd (Analise PCA)
  if ((quais[4] == 1) && (HMM.usaPCA == 1))
    pca(HMM.nomePCA,HMM.ordemRed,&HMM.ordem,HMM.n_par,Nwindow,&par);
  //----------------------------------------------------------------------------

  // Desalocando ponteiros
  free(x1);
  free(x2);
  free(xw);

  if (quais[4] == 1)
  {
    for (i=0;i<Nwindow;i++)
    {
      free(maux[i]);
      free(daux[i]);
      free(ddaux[i]);
    }
    free(maux);
    free(daux);
    free(ddaux);
  }

  free(mod_fft);

  //free(MediaMel);
  /*
  // Salvando parametros em arquivo
  arquivo = fopen("c:\\users\\iog\\continuo\\comum\\perfil.dat","wb");
  for (i=0;i<n_par;i++)
    for (j=0;j<Nwindow;j++)
      for (k=0;k<ordem[i];k++)
        fwrite(&par[i][j][k],sizeof(double),1,arquivo);
  fclose(arquivo);
  */

  // Valores de retorno da funcao
  *par1 = par;
  *Nwindow1 = Nwindow;
}
//------------------------------------------------------------------------------
// Funcao que calcula parametros delta
void calcDelta
(
  int delta, // numero de janelas adjacentes a serem utilizadas no calculo dos parametros delta
  int n_quadros, // numero de quadros do sinal
  int ordem, // ordem os vetores de parametros
  double **param, // parametros dos quais se deseja calcular os parametros delta
  double **par_delta // armazena os parametros delta
)
{
  int register i,j,k; // contadores

  for (i=0;i<n_quadros;i++)
    for (j=0;j<ordem;j++)
      par_delta[i][j]=0.0;

  // Calculando parametros delta
  for (i=0;i<n_quadros;i++)
  {
    for (j=-delta;j<=delta;j++)
      if (((i+j)>=0) && ((i+j)<n_quadros))
        for (k=0;k<ordem;k++)
          par_delta[i][k] += -(double)j*param[i+j][k];
  }
  for (i=0;i<n_quadros;i++)
    for (j=0;j<ordem;j++)
      par_delta[i][j] /= (2.0*(double)delta+1.0);

}
//------------------------------------------------------------------------------
