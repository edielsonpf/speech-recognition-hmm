#ifndef LOADVOCAB_H_
#define LOADVOCAB_H_

//------------------------------------------------------------------------------
// Programa que compila descricao do vocabulario de palavras.
// A entrada deve ser um arquivo texto seguindo o padrao:
//
// *fonemas
// k
// ...
// u
// *fim
// *classes
// sub
// ...
// prep
// *fim
// *vocab
// casa/k a z a/<media>/<desvio padrao>/classes
// ...
// *fim
//
// O arquivo de sai­da tera o seguinte formato
//
//  #palavras(long)
//  #caracteres_na_palavra(char),"c","a","s","a",
//  #fones(char),<fone1>(char),<fone2>(char),...,
//  <media>(double),<desvio>(double),
//  #classes(char),<classe1>(int),<classe2>(int),...
//  ...
int CompilaVocabulario
(
  char *nome // nome do arquivo com o vocabulario
);

//---------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// FunÃ§Ã£o que extrai frases de um arquivo texto.
// A frase obtida e retornada em "frase".
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

//------------------------------------------------------------------------------
// Funcao que carrega as subunidades, os modelos foneticos, modelos de duracao e
// classes das palavras do vocabulario codificado pela funcao 'compila_vocabulario()'
void CarregaVocabulario
(
  char *nome, // nome do arquivo de vocabulario compilado
  struct vocab *vocabulario1 // estrutura que armazena o vocabulario
);
	
//------------------------------------------------------------------------------
// Funcao que desaloca o ponteiro vocabulario
void DesalocaVocabulario
(
  struct vocab *vocabulario1 // estrutura que armazena o vocabulario
);

#endif /*LOADVOCAB_H_*/
