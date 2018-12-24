// Nesta unidade estao definidas funcoes para pre-processamentos do sinal de
// voz a ser parametrizado.
//---------------------------------------------------------------------------
#include <math.h>
#include <stdlib.h>

//#include "preProc.h"

//------------------------------------------------------------------------------
// Funcao que remove a componente DC de um sinal de voz
void removeDC
(
  int Nam,
  double *x
)
{
  double DC = 0.0;
  int i;

  for (i=0;i<Nam;i++)
    DC += x[i];
  DC /= (double)Nam;
  for (i=0;i<Nam;i++)
    x[i] -= DC;

}

//------------------------------------------------------------------------------
// Implementa filtro de pre-enfase de primeira ordem
void preEnfase
(
  double CP, // coeficiente de pre-enfase
  int tamanho, // numero de amostras de x
  double *x // sinal a ser pre-enfatizado
)
{
  int register i; // contador
  double *y; // sinal depois da pre-enfase

  y = malloc(sizeof(double)*tamanho);

  for (i=0;i<(tamanho-1);i++)
    y[i] = x[i + 1] - CP * x[i];

  for (i=0;i<tamanho;i++)
    x[i] = y[i];

  free(y);
}

//------------------------------------------------------------------------------
// Aplica uma janela de Hamming ao sinal
void hamming
(
  int tamanho, // numero de amostras do sinal
  double *x // sinal a ser janelado
)
{
  double a; // var aux para o calculo da janela
  int register i; // contador

  // Calculando parametros
  for(i=0;i<tamanho;i++)
  {
    a = 2 * M_PI * i / (tamanho - 1);
    x[i] *= (0.54 - 0.46 * cos(a));
  }
}
