#ifndef WAV_H_
#define WAV_H_

// Funcao que abre um arquivo wav de 16 bits e devolve um ponteiro contendo
// os dados, bem como o numero de amostras
void leWav
(
  short *BitsPorAmostra1, // numero de bits por amostra
  int *fs1, // frequencia de amostragem
  char *nome_arquivo, // nome do arquivo a ser lido
  int *nSamples, // number of samples
  short **x1 // samples of the audio file
);

#endif /*WAV_H_*/
