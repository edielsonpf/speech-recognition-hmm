#ifndef LOADGRAMMAR_H_
#define LOADGRAMMAR_H_

// Funcao que carrega modelo de gramatica bigram de palavras, gerado pelo programa Gramatica.exe

void loadGrammar
(
  char* nome, // nome do arquivo com os dados da gramatica
  long* n_termos, // numero de pares na gramatica
  struct bigram **grammar1 // pares da gramatica
);

#endif /*LOADGRAMMAR_H_*/
