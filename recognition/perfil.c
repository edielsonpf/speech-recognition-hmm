#include "perfil.h"

//---------------------------------------------------------------------------
// Funcao para calculo de perfis de energia
void calcPerfilEnergia
(
  int fs,          // frequencia de amostragem do sinal
  double *mod_fft, // modulo da FFT do sinal
  int N,           // numero de pontos da FFT
  int nPerfis,     // numero de perfis desejado (ordem do parametro)
  double *Perfis   // retorna os perfis de energia
)
{
  double ETotal; // energia total do frame
  int register i,j; // contadores
  double Limiar; // limiar para cada um dos perfis
  double Parcial; // energia parcial ate a raia atual
  double ParcialAnterior; // energia parcial ate a raia anterior
  double Passo; // intervalo entre um perfil e outro
  double x1, x2; // variaveis auxiliares para interpolacao linear dos valores dos perfis

  // Calculando energia total do frame
  ETotal = mod_fft[0] + mod_fft[N/2];
  for (i=1;i<(N/2);i++)
    ETotal += 2*mod_fft[i];

  // Calculando perfis de energia
  Passo = ETotal / (double)(nPerfis + 1);
  Limiar = 0.0;
  ParcialAnterior = 0.0;
  Parcial = mod_fft[0];
  j = 1;
  for (i=0;i<nPerfis;i++)
  {
    Limiar += Passo;
    while (Parcial <= Limiar)
    {
      ParcialAnterior = Parcial;
      Parcial += 2*mod_fft[j];
      j++;
    }
    if (Parcial > Limiar) // passou do ponto
    {
      x1 = ((double)j-2.0)*(double)fs/(double)N; // em Hertz
      x2 = ((double)j-1.0)*(double)fs/(double)N; // em Hertz
      Perfis[i] = InterpLinear(x1,ParcialAnterior,x2,Parcial,Limiar);
    }
    else
      Perfis[i] = ((double)j-1.0)*(double)fs/(double)N; // em Hertz
  }
}

// Dados dois pontos (x1,y1) e (x2,y2), calcula o valor de x, que passa pela reta
// dada por estes pontos, na posicao intermediaria y.
double InterpLinear
(
  double x1, // componente x do ponto 1
  double y1, // componente y do ponto 1
  double x2, // componente x do ponto 2
  double y2, // componente y do ponto 2
  double y // componente y do ponto intermediario
)
{
  double a; // coeficiente angular da reta que passa pelos pontos 1 e 2
  double b; // termo independente da reta que passa pelos pontos 1 e 2
  double aux_x, aux_y; // variaveis auxiliares
  double x; // valor de retorno da funcao

  aux_x = x1 - x2;
  aux_y = y1 - y2;

  a = aux_y / aux_x;
  b = y1-a*x1;

  x = (y - b)/a;

  return x;
}
