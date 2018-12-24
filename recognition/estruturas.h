#ifndef ESTRUTURAS_H_
#define ESTRUTURAS_H_

// Estrutura para armazenar o arquivo de vocabulario
struct vocab
{
  int**  classes;    // Classes gramaticais das palavras.
  int*   Ddur;       // Desvio padrao das duracoes das palavras do vocabulario
  int*   Mdur;       // Media das duracoes das palavras do vocabulario
  int**  Mw;         // Transcricao das palavras
  char** Mws;        // Modelo fonetico das palavras
  int    n_fones;    // numero de subunidades foneticas utilizadas p/ a transcricao das palavras
  int    n_palavras; // numero de palavras no vocabulario
};

// Estrutura que armazena os dados das misturas referentes as fdp's de emissao de simbolos
struct mistura
{
  double* c; // coeficientes de ponderacao para cada gaussiana
  double** m; // medias para cada uma das gaussianas
  double** v; // variancias para cada uma das gaussianas
};


// Estrutura para armazenar os modelos HMM das subunidades acusticas
struct modelo_HMM
{
  int* delta;          // numero de janelas adjacentes utilizadas para calculo de parametros delta
  int n_estados;       // numero de estados para cada subunidade fonetica
  int n_fones;         // numero de subunidades acusticas no modelo
  int* n_gaussianas;   // numero de gaussianas na mistura, para cada um dos parametros
  int n_par;           // numero de parametros utilizados no reconhecimento
  char** parametros;   // descricao dos parametros utilizados no reconhecimento (p. ex. 'mel')
  int* ordem;          // ordem dos parametros
  double*** A;         // matrizes de transicao
  struct mistura*** B; // matrizes de emissao
  int salto_maximo;    // salto maximo permitido entre estados no modelo HMM das subunidades
  int usaPCA;          // usa (1) ou nao (0) a PCA para reducao da dimensao dos vetores acusticos
  char nomePCA[512];   // nome do arquivo com a matriz de transformacao PCA
  int ordemRed;        // nova dimensao dos vetores acusticos, depois da transformacao PCA
};

// Estrutura p/ armazenar gramatica bigram de palavras
struct bigram
{
  int primeira; // primeira palavra da seq. (p1)
  int segunda;  // segunda palavra da seq. (p2)
  float prob;   // p(p1/p2)
};

//Estrutura para armazenar as locucoes a serem reconhecidas
struct locucao
{
  int n_quadros; // numero de quadros com o qual foi parametrizada a locucao
  double*** par;  // locucao parametrizada (n_par x n_janelas x ordem)
};

//---------------------------------------------------------------------------
// Estrutura para armazenar as opcoes de configuracao do modo de busca
struct configBusca
{
  int searchAlgorithm;         // 0 - LB, 1 - OS, 2 - HN
  int numberOfLevels;           // numero de niveis com que e feita a busca (# palavras na locucao)
  //int parada_automatica; // verifica se o usuario quer que o sistema determine automaticamente
                          // o numero de palavras na locucao (somente p/ level building)
  int useBeamSearch;   // verifica se o usuario deseja usar Viterbi Beam Search (1) ou nao (0)
  float beamSearchThreshold;          // limiar de poda p/ o Beam Search
  int useWordDurationModel;       // verifica se o usuario deseja usar modelo de duracao de palavras (1) ou nao (0)
  int useGrammar;     // verifica se o usuario deseja usar gramatica bigram de palavras (1) ou nao
};

// Estrutura de um no da arvore p/ utilizacao com o algoritmo Herrmann-Ney
struct no
{
  // Variaveis p/ gerenciamento das arvores
  int profundidade; // profundidade do no
  int fone; // fone a que corresponde o no
  int palavra; // se este no for o ultimo da palavra, indica a qual palavra corresponde, cc. o seu valor e -1
  int pai; // no pai referente a este no
  int n_filhos; // numero de descendentes deste no
  int *filhos; // identificadores dos descendentes deste no

  // Variaveis p/ gerenciamento da busca
  double* like; // verossimilhanca acumulada do melhor caminho ate este no
  int* elapse; // duracao do melhor caminho para o no
  double pt; // probabilidade de transicao do ultimo estado deste fone para o seguinte
  double like_ant; // prob. acumulada do ultimo estado no instante anterior
  double elapse_ant; // duracao do caminho otimo ate o ultimo estado no instante anterior
  int ativo; // indica que este no esta ativo (1) ou nao (0) na hora da busca
};


#endif /*ESTRUTURAS_H_*/
