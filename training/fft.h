#ifndef FFT_H_
#define FFT_H_

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
);

#endif /*FFT_H_*/
