#ifndef ENERGIA_H_
#define ENERGIA_H_

// Calcula a energia um frame
double calcEnergia
(
  double *xw,  // sinal a ser analisado
  int Lwindow // numero de amostras do sinal
);

// Funcao que calcula parametros delta e delta-delta p/ os parametros log energia normalizada
void CalculaDeltaEnergia
(
  double **energia, // vetor com os parametros log-energia, delta e delta-delta
  int n_quadros, // numero de quadros na locucao
  int delta // janelas a esquerda e a direita a serem consideradas para o calculo dos parametros delta
);

#endif /*ENERGIA_H_*/
