#include "estruturas.h"
#include "constantes.h"
#include "grammar.h"

// Aplica gramatica bigram de palavras
int Pode
(
  int palavra1, // primeira palavra
  int palavra2, // segunda palavra
  double *probabilidade, // p(segunda|primeira)
  long n_termos, // numero de pares na gramatica
  struct bigram *gram // pares da gramatica
)
{
  // Declaracao das variaveis locais
  int achou_primeira; // verifica se palavra1 consta (1) ou nao (0) do vocabulario
  int achou_segunda; // verifica se palavra2 consta (1) ou nao (0) do vocabulÃ¡rio
  int L, H, I, C;
  int pos; // posicao da busca no vetor gramatica

  achou_primeira = 0;
  L = 0;
  H = n_termos - 1;
  while (L <= H)
  {
    I = (L + H) >> 1;
    C = gram[I].primeira - palavra1;
    if (C < 0)
      L = I + 1;
    else
    {
      H = I - 1;
      if (C == 0)
      {
        achou_primeira = 1;
        L = I;
      }
    }
  }
  pos = L;

  // Procurando ocorrencia da segunda palavra (palavra2)
  if (achou_primeira == 1)
  {
    achou_segunda = 0;
    while ((achou_segunda == 0) && (gram[pos].primeira == palavra1) && (pos < n_termos))
    {
      if (gram[pos].segunda == palavra2)
        achou_segunda = 1;
      else
      {
        if (gram[pos].segunda > palavra2)
        {
          pos--;
          while ((achou_segunda == 0) && (gram[pos].primeira == palavra1) && (pos >= 0))
          {
            if ((gram[pos].segunda == palavra2) && (gram[pos].primeira == palavra1))
              achou_segunda = 1;
            pos--;
          }
        }
        else
        {
          pos++;
          while ((achou_segunda == 0) && (gram[pos].primeira == palavra1) && (pos < n_termos))
          {
            if ((gram[pos].segunda == palavra2) && (gram[pos].primeira == palavra1))
              achou_segunda = 1;
            pos++;
          }
        }
      }
    }
    // Probabilidade da sequencia de palavras
    if (achou_segunda == 1)
      *probabilidade = (double)gram[pos].prob;
    else
      *probabilidade = -Inf;
  }

	if (achou_primeira==1 && achou_segunda==1)
    return 1;
	else
		return 0;
}

//------------------------------------------------------------------------------
