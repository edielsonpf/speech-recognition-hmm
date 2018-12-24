#ifndef ARVORE_H_
#define ARVORE_H_

// Gera arvores que representam as palavras do vocabulario
void GeraArvores
(
  int *MaxProfundidade1, // numero maximo de niveis observado
  long n_palavras, // numero de palavras no vocabulario
  int n_fones, // numero de fones utilizados p/ transcricao fonetica das palavras
  struct vocab vocabulario // vocabulario
);

// Gera arvore com fone inicial diferente
void NovaArvore
(
  int *n_arvores1, // numero de arvores geradas
  int **n_nos1, // numero de nos de cada arvore
  struct no ***vocab1, // arvores geradas
  int *MW, // modelos foneticos de cada palavra
  int palavra, // palavra sob analise
  int n_estados // numero de estados dos modelos HMM dos fones
);

// Adiciona palavra a arvore ja existente
void AtualizaArvore
(
  int *n_nos1, // numero de nos desta arvore
  int *n_arvores1, // numero de arvores
  struct no **vocab1, // arvores geradas
  int qual, // arvore a ser atualizada
  int *MW, // modelos foneticos de cada palavra
  int palavra, // palavra sob analise
  int n_estados // numero de estados dos modelos HMM dos fones
);

// Carrega arvores geradas pela funcao GeraArvores
void CarregaArvores
(
  int n_niveis, // numero maximo de palavras na locucao
  int *n_arvores1, // numero de arvores criadas no vocabulario
  int n_estados, // numero de estados dos modelos HMM dos fones
  int **n_nos1, // numero de nos em cada arvore
  struct no ****vocab1 // armazena as palavras do vocabulario em arvore
);


// Apresenta as arvores geradas na tela
void ApresentaArvores
(
  int n_arvores, // numero de arvores geradas
  int* n_nos, // numero de nos em cada arvore
  struct no **vocab, // armazena as palavras do vocabulario em arvore
  struct vocab vocabulario // vocabulario
);


// Desaloca arvores alocadas na funcao 'CarregaArvores
void DesalocaArvores
(
  int n_niveis, // numero maximo de palavras na locucao
  int n_arvores, // numero de arvores geradas
  int **n_nos1, // numero de nos em cada arvore
  struct no ****vocab1 // armazena as palavras do vocabulario em arvore
);

// Gera listagem dos nos pertencentes a uma dada profundidade dentro da arvore
void GeraListas
(
  int MaxProfundidade, // profundidade maxima observada em todas as arvores
  int n_arvores, // numero de arvores geradas
  int *n_nos, // numero de nos em cada arvore
  struct no **vocab // armazena as palavras do vocabulario em arvore
);

// Funcao que que apresenta na tela as listas geradas pela funcao 'GeraListas'
void ApresentaListas
(
  int ***lista, // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int n_arvores, // numero de arvores geradas
  int *profundidades, // profundidade maxima observada em cada uma das arvores
  int **quantos, // indica quantos elementos existem em cada nivel de cada uma das listas
  struct no **vocab // armazena as palavras do vocabulario em arvore
);
  
// Carrega as listas geradas pela funcao GeraLista
void CarregaListas
(
  int ****lista1, // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int n_arvores, // numero de arvores geradas
  int **profundidades1, // profundidade maxima observada em cada uma das arvores
  int *n_nos, // numero de nos em cada arvore
  int ***quantos1 // indica quantos elementos existem em cada nivel de cada uma das listas
);

// Desaloca ponteiros 'lista', 'profundidades' e 'quantos', alocados na funcao 'CarregaListas'
void DesalocaListas
(
  int ****lista1, // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int n_arvores, // numero de arvores geradas
  int **profundidades1, // profundidade maxima observada em cada uma das arvores
  int ***quantos1 // indica quantos elementos existem em cada nivel de cada uma das listas
);

#endif /*ARVORE_H_*/
