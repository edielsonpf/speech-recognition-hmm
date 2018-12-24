#ifndef PCA_H_
#define PCA_H_

// Realiza reducao da dimensao do vetor de parametros atraves da Analise de
// Componente Principal (por enquanto usado somente para parametros mdd)
void pca
(
  char* NomePCA,  // nome da matriz de transformacao gerada pelo programa KL.exe
  int OrdemRed,   // ordem desejada para o vetor de parametros apos a reducao pela PCA
  int** ordem1,   // ordem dos vetores de cada um dos parametros
  int n_par,      // numero de parametros a serem calculados
  int Nwindow,    //numero de quadros com o qual foi analisado o sinal
  double ****par1 // sinal parametrizado
);

#endif /*PCA_H_*/
