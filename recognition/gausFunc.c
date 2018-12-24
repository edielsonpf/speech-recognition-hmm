//---------------------------------------------------------------------------
#include <stdlib.h>
#include <math.h>

#include "constantes.h"
#include "estruturas.h"
#include "gausFunc.h"

// Funcao que implementa uma fdp gaussiana multidimensional
double Gauss
(
  double *media, // vetor com as medias
  int ordem,     // ordem do vetor de parametros
  double *var,   // vetor com as variancias
  double *x      // vetor de parametros
)
{
  // Declaracao das variaveis locais
  double aux=0.0,aux1,aux2 = 0.0; // variaveis auxiliares
  double det; // determinante da matriz de covariancia
  double dif; // var aux p/ o calculo da gaussiana
  int i; // contador
  double prob=0.0; // probabilidade do vetor dada a distribuicao

  // (2*pi)^(ordem/2)
  aux1 = 2.0*PI;
  aux2 = ordem / 2.0;
  aux1 = pow(aux1,aux2);

  // Calculando determinante da matriz de covariancia
  det = 1.0;
  for (i=0;i<ordem;i++)
    det *= var[i];

  if (det != 0)
  {
    aux2 = fabs(det);
    aux2 = pow(aux2,0.5);

    for (i=0;i<ordem;i++)
    {
      dif = x[i] - media[i];
      aux += (dif*dif)/var[i];
    }
    aux *= (-0.5);
    aux = exp(aux);

    prob = aux/(aux1*aux2);
  }

  return(prob);

}

//------------------------------------------------------------------------------
// Funcao que calcula a probabilidade de um simbolo para uma mistura de gaussianas
void ProbSimbolo
(
  struct mistura **b, // matriz de emissao para os modelos das locucoes
  double ****bj1,     // probabilidade de cada simbolo da locucao, dado o modelo
  int *comprimentos,  // numero de subunidades foneticas das transcricoes de cada locucao
  int frase,          // conta as locucoes de treinamento
  int n_estados,      // numero de estados para cada modelo HMM de subunidade
  int *n_gaussianas,  // numero de gaussianas por mistura
  double *****Nm1,    // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
  int n_par,          // numero de parametros de treinamento
  int *ordem,         // dimensao do vetor para cada um dos parâmetros
  int tamanho,        // numero de quadros da locucao de treinamento
  double ***x         // vetores de parametros de treinamento
)
{
  // Declaracao das variaveis locais
  double ***bj; // probabilidade de cada simbolo da locucao, dado o modelo
  int i,j,k,l; // contadores
  double ****Nm; // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo

  // Alocando memoria
	bj = malloc(sizeof(double **)*n_par);
	for (i=0;i<n_par;i++)
		bj[i] = malloc(sizeof(double *)*(n_estados*comprimentos[frase]));
	for (i=0;i<n_par;i++)
		for (j=0;j<n_estados*comprimentos[frase];j++)
			bj[i][j] = malloc(sizeof(double)*tamanho);

	Nm = malloc(sizeof(double ***)*n_par);
	for (i=0;i<n_par;i++)
		Nm[i] = malloc(sizeof(double **)*tamanho);
	for (i=0;i<n_par;i++)
		for (j=0;j<tamanho;j++)
			Nm[i][j] = malloc(sizeof(double *)*(n_estados*comprimentos[frase]));
	for (i=0;i<n_par;i++)
		for (j=0;j<tamanho;j++)
			for (k=0;k<n_estados*comprimentos[frase];k++)
				Nm[i][j][k] = malloc(sizeof(double)*n_gaussianas[i]);

  // Zerando as variaveis alocadas
  for (i=0;i<n_par;i++)
    for (j=0;j<n_estados*comprimentos[frase];j++)
      for (k=0;k<tamanho;k++)
      {
        bj[i][j][k] = 0.0;
        for (l=0;l<n_gaussianas[i];l++)
          Nm[i][k][j][l] = 0.0;
      }

  // Calculando probabilidades
  for (i=0;i<n_par;i++)
    for (j=0;j<n_estados*comprimentos[frase];j++)
      for (k=0;k<tamanho;k++)
        for (l=0;l<n_gaussianas[i];l++)
        {
          Nm[i][k][j][l] = b[i][j].c[l]*Gauss(b[i][j].m[l],ordem[i],b[i][j].v[l],x[i][k]);
          bj[i][j][k] += Nm[i][k][j][l];
        }

  for (i=0;i<n_par;i++)
    for (j=0;j<n_estados*comprimentos[frase];j++)
      for (k=0;k<tamanho;k++)
        if (bj[i][j][k] != 0.0)
          for (l=0;l<n_gaussianas[i];l++)
            Nm[i][k][j][l] /= bj[i][j][k];
        else
          for (l=0;l<n_gaussianas[i];l++)
            Nm[i][k][j][l] = 0.0;

  *bj1 = bj;
  *Nm1 = Nm;
}

//------------------------------------------------------------------------------
// Funcao que calcula a probabilidade de um simbolo para uma mistura de gaussianas
// alternativa para o aplicativo de Segmentação (Ynoguti, 19/09/2003)
void ProbSimbolo1
(
  struct mistura **b, // matriz de emissao para os modelos das locucoes
  double ****bj1,     // probabilidade de cada simbolo da locucao, dado o modelo
  int comprimento,  // numero de subunidades foneticas das transcricoes de cada locucao
  int n_estados,      // numero de estados para cada modelo HMM de subunidade
  int *n_gaussianas,  // numero de gaussianas por mistura
  double *****Nm1,    // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
  int n_par,          // numero de parametros de treinamento
  int *ordem,         // dimensao do vetor para cada um dos parâmetros
  int tamanho,        // numero de quadros da locucao de treinamento
  double ***x         // vetores de parametros de treinamento
)
{
  // Declaracao das variaveis locais
  double ***bj; // probabilidade de cada simbolo da locucao, dado o modelo
  int i,j,k,l; // contadores
  double ****Nm; // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo

  // Alocando memoria
	bj = malloc(sizeof(double **)*n_par);
	for (i=0;i<n_par;i++)
		bj[i] = malloc(sizeof(double *)*n_estados*comprimento);
	for (i=0;i<n_par;i++)
		for (j=0;j<n_estados*comprimento;j++)
			bj[i][j] = malloc(sizeof(double)*tamanho);

	Nm = malloc(sizeof(double ***)*n_par);
	for (i=0;i<n_par;i++)
		Nm[i] = malloc(sizeof(double **)*tamanho);
	for (i=0;i<n_par;i++)
		for (j=0;j<tamanho;j++)
			Nm[i][j] = malloc(sizeof(double *)*n_estados*comprimento);
	for (i=0;i<n_par;i++)
		for (j=0;j<tamanho;j++)
			for (k=0;k<n_estados*comprimento;k++)
				Nm[i][j][k] = malloc(sizeof(double)*n_gaussianas[i]);

  // Zerando as variaveis alocadas
  for (i=0;i<n_par;i++)
    for (j=0;j<n_estados*comprimento;j++)
      for (k=0;k<tamanho;k++)
      {
        bj[i][j][k] = 0.0;
        for (l=0;l<n_gaussianas[i];l++)
          Nm[i][k][j][l] = 0.0;
      }

  // Calculando probabilidades
  for (i=0;i<n_par;i++)
    for (j=0;j<n_estados*comprimento;j++)
      for (k=0;k<tamanho;k++)
        for (l=0;l<n_gaussianas[i];l++)
        {
          Nm[i][k][j][l] = b[i][j].c[l]*Gauss(b[i][j].m[l],ordem[i],b[i][j].v[l],x[i][k]);
          bj[i][j][k] += Nm[i][k][j][l];
        }

  for (i=0;i<n_par;i++)
    for (j=0;j<n_estados*comprimento;j++)
      for (k=0;k<tamanho;k++)
        if (bj[i][j][k] != 0.0)
          for (l=0;l<n_gaussianas[i];l++)
            Nm[i][k][j][l] /= bj[i][j][k];
        else
          for (l=0;l<n_gaussianas[i];l++)
            Nm[i][k][j][l] = 0.0;

  *bj1 = bj;
  *Nm1 = Nm;
}



/*
// Funcao que calcula a probabilidade de uma mistura de gaussianas
double TGauss::Prob
(
  struct mistura b, // coef. de ponderacao, media e variancia de cada gaussiana
  int n_gaussianas, // numero de gaussianas na mistura
  int ordem, // ordem do vetor x
  double *x // vetor de dados
)
{
  double p = 0.0;
  int register i;

  for (i=0;i<n_gaussianas;i++)
    p += b.c[i]*Gauss(b.m[i],ordem,b.v[i],x);

  return p;
}
*/


// Funcao que calcula a log-probabilidade de uma mistura de gaussianas
double LogProb
(
  struct mistura b, // coef. de ponderacao, media e variancia de cada gaussiana
  int n_gaussianas, // numero de gaussianas na mistura
  int ordem, // ordem do vetor x
  double *x // vetor de dados
)
{
  double aux; // var aux p/ calculo de p
  double p = 0.0; // probabilidade da mistura de gaussianas
  int i; // contador

  for (i=0;i<n_gaussianas;i++)
  {
    aux = Gauss(b.m[i],ordem,b.v[i],x);
    if (aux > Zero)
      p += b.c[i]*aux;
  }
  if (p > Zero)
    p = log10(p);
  else
    p = -Inf;

  return p;
}

// Funcao que calcula contribuicoes locais de cada arco, p/ cada subunidade fonetica
void CalcContrLocal
(
  int fone, // fone sob analise
  double ****log_local1, // contribuicao local de cada arco p/ a verossinilhanca total
  struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  int t, // frame sob analise
  struct locucao x // locucao a ser reconhecida
)
{
  int j,k; // contadores
  double ***log_local; // contribuicao local de cada arco p/ a verossinilhanca total
  double *prob_emissao; // probabilidade de emissao p/ cada um dos estados de um fone

  prob_emissao = malloc(sizeof(double)*HMM.n_estados);

  log_local = *log_local1;

  for (j=0;j<HMM.n_estados;j++)
  {
    prob_emissao[j] = 0.0;
    for (k=0;k<HMM.n_par;k++) // LogProb e definido em Gausfunc.h
      prob_emissao[j] += LogProb(HMM.B[k][fone][j],HMM.n_gaussianas[k],HMM.ordem[k],x.par[k][t]);
  }
  for (j=0;j<HMM.n_estados;j++)
    for (k=0;k<(HMM.salto_maximo+1);k++)
      if ((j+k) < HMM.n_estados)
        log_local[fone][j][k] = HMM.A[fone][j][k] + prob_emissao[j+k];
  log_local[fone][0][0] = prob_emissao[0]; // probabilidade de emissao no primeiro estado do fone
  log_local[fone][2][1] = HMM.A[fone][2][1]; // probabilidade de transicao do ultimo estado do fone p/ proximo fone

  free(prob_emissao);

  *log_local1 = log_local;
}
