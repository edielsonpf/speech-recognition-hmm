#ifndef LPC_H_
#define LPC_H_

// Calculo dos paramentros LPC, através do algoritmo de Durbin.
int calcLPC
(
  double *xw, // amostras do sinal
  int Lwindow, // numero de amostras de xw
  int ordem,   // ordem dos parametros
  double *lpc // ponteiro onde serao armazenados os parametros calculados
);

// Função para cálculo da função de autocorrelação
void autoCorrelate
(
  double *xw, // amostras do sinal
  double *r, // vetor de autocorrelacao
  int order, // ordem da analise lpc
  int Lwindow // numero de amostras de xw
);

int normalize_corr
(
  double *r, // vetor de autocorrelacao
  int order // ordem da analise lpc
);

#endif /*LPC_H_*/
