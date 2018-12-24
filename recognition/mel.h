#ifndef MEL_H_
#define MEL_H_

//  Funcao para calculo dos coeficientes mel-cepstrais
//  Os coef mel-cepstrais sao obtidos segundo o algoritmo de Mermelstein;
//  os passos a seguir sao:
//
//  - Calculo da FFT das am's pertencentes a janela de analise.
//  - Calculo do quadrado do modulo da FFT (equivale a energia).
//  - Filtragem do modulo acima pelo banco de filtros na escala mel,
//    obtendo-se a energia por filtro.
//  - Calculo do log da energia por filtro.
//  - Calculo da transformada inversa, obtendo-se os coef's mel-cepstrais
//    (esto e' implementado calculando-se a DCT: discrete cosine transform)
//
//  OBS: Lwindow: num de ams/janela; xw: am's janeladas
// +------------------------------------------------------------------------+
void calcMel
(
  int fs,          // frequencia de amostragem do sinal
  double *mod_fft, // modulo ao quadrado da FFT (energia)
  int N,           // numero de pontos da FFT
  int ordem,       // numero de parametros mel a serem calculados
  double *Mel_ceps // ponteiro onde serao armazenados os parametros calculados
);

#endif /*MEL_H_*/
