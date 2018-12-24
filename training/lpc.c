#include <stdlib.h>

#include "lpc.h"
//---------------------------------------------------------------------------
// Calculo dos paramentros LPC, através do algoritmo de Durbin.
int calcLPC
(
  double *xw, // amostras do sinal
  int Lwindow, // numero de amostras de xw
  int ordem,   // ordem dos parametros
  double *lpc // ponteiro onde serao armazenados os parametros calculados
)
{
  int 	i,j;
  double sum,error;
  double *b;
  double *r; // vetor de autocorrelacao
  double *acf; // coeficientes LPC
  double *kcf; // parametros K

	// Allocating memory
  b = malloc(sizeof(double)*(ordem+1));
  r = malloc(sizeof(double)*(ordem+1));
	acf = malloc(sizeof(double)*(ordem+1));
	kcf = malloc(sizeof(double)*(ordem+1));
	
  // Calculo do vetor de autocorrelacao
  autoCorrelate(xw,r,ordem,Lwindow);

  /* initialize predictor coeffs. */
  acf[0] = 1.0;
  for(i=1; i<=ordem; i++)
    acf[i]=0.0;

  /* normalize autocorrs */
  if(!normalize_corr(r,ordem))
    return 1;

  /* durbin's recursion */
  for(error = r[0], i = 1; i <= ordem; i++)
  {
    for(sum = r[i], j=1; j < i; j++)
      sum -= acf[j] * r[i - j];
    acf[i] = kcf[i] = sum / error;
    for(j = 1; j < i; j++)
      b[j] = acf[i - j];
    for(j = 1; j < i; j++)
      acf[j] -= kcf[i] * b[j];
    error = (1 - kcf[i] * kcf[i]) * error;
  }

  free(b);
  free(r);
  free(acf);
	free(kcf);
	
  return 1;
}
//------------------------------------------------------------------------------
// Função para cálculo da função de autocorrelação
void autoCorrelate
(
  double *xw, // amostras do sinal
  double *r, // vetor de autocorrelacao
  int order, // ordem da analise lpc
  int Lwindow // numero de amostras de xw
)
{
	int i,j; // counters
  //double energy;

  for (i=0;i<=order;i++)
  {
    double sum = 0.0;
    for (j=1;j<=Lwindow-i;j++)
      sum += xw[j]*xw[j+i];
    r[i] = sum;

    /*
    if (i==0)
      energy = sum;
    else
      r[i] = sum;
    */
  }
  //return energy;
}

/*------------------------------------------------------------------------*
* 	normalize r's to prevent round off error
*------------------------------------------------------------------------*/
int normalize_corr
(
  double *r, // vetor de autocorrelacao
  int order // ordem da analise lpc
)
{
  int	i;
  if(r[0] == (double)0)
    return 0;
  for(i=1; i<=order; i++)
    r[i] /= r[0];
  r[0] = 1;
  return 1;
}
//..............................................................................
