#ifndef PAR_H_
#define PAR_H_

#include "estruturas.h"

// funcao que calcula os parametros mel-cepstrais de um sinal de entrada 'x'.
// Calcula os coeficientes ate ordem 'ordem' e armazena os resultados em 'y'.
void calcPar
(
  char* nome_arquivo, // nome do arquivo WAV a ser parametrizado
	struct modelo_HMM HMM, // acoustic models (contain information about which parameters to use, dimension of parameters, whether use or not PCA, etc.)
  int* Nwindow1,      // numero de quadros com o qual foi analisado o sinal
	double**** par1    // sinal parametrizado
);
	
// Funcao que calcula parametros delta
void calcDelta
(
	int delta, // numero de janelas adjacentes a serem utilizadas no calculo dos parametros delta
	int n_quadros, // numero de quadros do sinal
	int ordem, // ordem os vetores de parametros
	double **param, // parametros dos quais se deseja calcular os parametros delta
	double **par_delta // armazena os parametros delta
);

#endif /*PAR_H_*/
