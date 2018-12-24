#ifndef PERFIL_H_
#define PERFIL_H_

// Funcao para calculo de perfis de energia
void calcPerfilEnergia
(
  int fs,          // frequencia de amostragem do sinal
  double *mod_fft, // modulo da FFT do sinal
  int N,           // numero de pontos da FFT
  int nPerfis,     // numero de perfis desejado (ordem do parametro)
  double *Perfis   // retorna os perfis de energia
);

// Dados dois pontos (x1,y1) e (x2,y2), calcula o valor de x, que passa pela reta
// dada por estes pontos, na posicao intermediaria y.
double InterpLinear
(
  double x1, // componente x do ponto 1
  double y1, // componente y do ponto 1
  double x2, // componente x do ponto 2
  double y2, // componente y do ponto 2
  double y // componente y do ponto intermediario
);

#endif /*PERFIL_H_*/
