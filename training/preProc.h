#ifndef PREPROC_H_
#define PREPROC_H_

// Retira componente DC de sinal de voz
void removeDC
(
	int Nam,       // numero de amostras do sinal
	double *x // sinal a ser analisado
);

// Implementa filtro de pre-enfase de primeira ordem
void preEnfase
(
	double CP, // coeficiente de pre-enfase
	int tamanho, // numero de amostras de x
	double *x // sinal a ser pre-enfatizado
);

// Aplica uma janela de Hamming ao sinal
void hamming
(
	int tamanho, // numero de amostras do sinal
	double *x // sinal a ser janelado
);

#endif /*PREPROC_H_*/
