#include <math.h>

#include "constantes.h"

// Funcao que aplica penalizacao segundo modelo de duracao de palavras
double wordDurationModel
(
  double like, // verossimilhanca
  int elapse, // tempo de duracao (em quadros) da palavra
  int palavra, // palavra sob analise
  int *Ddur, // desvio padrao da duracao das palavras
  int *Mdur // duracao media das palavras
)
{

  int aux; // var intermediaria p/ calculo da fdp gaussiana
  int desvp; // desvio padrao da duracao da palavra
  int duracao; // duracao da palavra
  double prob; // prob da duracao segundo o modelo de duracao

  duracao = elapse*10;

  if (Ddur[palavra] < (int)(Mdur[palavra]/3))
    desvp = (int)(Mdur[palavra]/3);
  else
	desvp = Ddur[palavra];
  if (desvp == 0)
	desvp = 10;

  aux = (int)((duracao-Mdur[palavra])/desvp);
  aux *= aux;

  if (aux < 690)
	prob = c1+log10(desvp)+0.5*aux*c2; // c1 e c2 sao definidos em Gausfunc.cpp
  else
	prob = -Inf;

  if (prob > -Inf)
	prob = like-prob;

  return prob;

}
