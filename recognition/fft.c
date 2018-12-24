//---------------------------------------------------------------------------
#include <math.h>
#include <stdlib.h>
//---------------------------------------------------------------------------
// Calculo da FFT de N pontos do vetor xw[]. Usa-se o algoritmo proposto no
// Prob 6.4, pagina 333, do livro "Digital Signal Processing" de Oppenheim &
// Schaffer. Nomenclatura usada no livro e aqui: LE= n1, LE1= n2, NV2= n3.
// Retorna o modulo da FFT em mod_fft.
void fft
(
  double *xw, // janela de analise
  int Lwindow, // tamanho da janela de analise 20 ms
  int fs, // frequencia de amostragem do sinal
  int M, // M=9 ==> FFT de 512 pontos; M=10 ==> FFT de 1024 pontos
  int N, // numero de pontos da FFT (N=2**M)
  double *mod_fft // modulo da FFT
)
{
  // Declaracao das variaveis locais
  int i,j,k; // contadores
  double Tr, Ti, Ur, Ui, Wr, Wi, b; // var auxiliares para o calculo da FFT
  int ip, n1, n2, n3; // var aux para o calculo da FFT
  double *Xr,*Xi;     // partes real e imaginaria da FFT

  Tr = Ti = Ur = Ui = Wr = Wi = b = 0.0;
  ip = n1 = n2 = n3 = 0;



  // Alocando memoria
  Xr = malloc(sizeof(double)*N);
  Xi = malloc(sizeof(double)*N);

  for(k=0; k<N; k++) // inicializacao
  {
	  Xr[k] = 0.0;
    Xi[k] = 0.0;
  }
  for(k=0; k<Lwindow; k++)
    Xr[k] = xw[k];

  // Embaralhamento da sequencia de entrada
  n3 = N / 2;
  j = 0;
  for(i=0; i<N-1; i++) // comeco do "for"
  {
    if(i < j)
    {
	  Tr = Xr[j];
      Xr[j] = Xr[i];
      Xr[i] = Tr;
   /* Ti = Xi[j];      // para seq.complexas
      Xi[j] = Xi[i];
      Xi[i] = Ti;  */
	}
	k = n3;
	while(k <= j)
    {
      j -= k;
      k /= 2;
    }
	j += k;
  }     // fim do "for" (e do embaralhamento)

  //-------------------------------
  // Calculo da FFT
  //-------------------------------
  for(k=1; k<=M; k++) // "for 1"
  {
    n1 = 1;
	for(i=1; i<=k; i++) // equivale a fazer 2**k
      n1 *= 2;
	n2 = n1 / 2;
	Ur = 1.0;
    Ui = 0.0;
	b = M_PI / n2;                  // M_PI=3.1415... (definido em <math.h>)
	Wr = cos(b);
    Wi = -sin(b);

	for(j=0; j<n2; j++) // "for 2"
    {
	  for(i=j; i<N; i+=n1) // "for 3"
      {
	    ip = i + n2;
		Tr = Xr[ip] * Ur - Xi[ip] * Ui;
        Ti = Xr[ip] * Ui + Xi[ip] * Ur;
		Xr[ip] = Xr[i] - Tr;
        Xi[ip] = Xi[i] - Ti;
		Xr[i] =  Xr[i] + Tr;
        Xi[i]  = Xi[i] + Ti;
      }     // fim "for 3"
	  Tr = Ur;
	  Ur = Ur * Wr - Ui * Wi;
      Ui = Tr * Wi + Ui * Wr;
	}        // fim "for 2"
  }           // fim "for 1"

  // Calculo do quadrado do modulo da FFT (equivale a energia)
  for(k=0; k<N; k++)
	 mod_fft[k] = Xr[k] * Xr[k] + Xi[k] * Xi[k];

  // Desalocando memoria
  free(Xr);
  free(Xi);
     
  // Fim do calculo da FFT. A seq. obtida esta na ordem sequencial
}
