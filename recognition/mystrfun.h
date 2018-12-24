#ifndef MYSTRFUN_H_
#define MYSTRFUN_H_

// Funcao que extrai palavras de uma frase
// As palavras sao devolvidas em "palav" e o numero de palavras em "numpal".
// O valor de retorno deve ser 1 para palavras corretas
int codpal
(
  char *frase, // frase a ser analisada
  char **palav, // estrutura onde serao devolvidas as palavras
  int *numpal, // numero de palavras na frase
  char sep1[], // separador
  char sep2[], // separador
  int MaxCrtPlv, // numero maximo de caracteres permitido em uma palavra
  int NPF // numero maximo de palavras na frase
);

// Funcao que extrai frases de um arquivo texto.
// A frase obtida eh retornada em "frase".
// O valor de retorno deve ser 1 para uma frase correta.
// Valor zero indica excesso de caracteres ou falta do ponto final.
int codfrs
(
  FILE *nomarq,
  char *frase,
  int MaxCrtFrs,
  char fim1[],
  char coment[]
);



#endif /*MYSTRFUN_H_*/
