#ifndef HERRMANNEY_H_
#define HERRMANNEY_H_

//---------------------------------------------------------------------------
// Implementa o algoritmo de busca p/ lexico organizado em arvore
void HerrmanNey
(
	struct locucao x, // locucao a ser reconhecida
  	char *frase, // frase reconhecida
  	struct configBusca searchConfig, // search algorithm configurations
	struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  	struct bigram *gramatica, // armazena gramatica bigram de palavras
  	int n_termos, // numero de pares da gramatica bigram
  	struct vocab vocabulario // estrutura que armazena o vocabulario
);

// Implementa o algoritmo de Viterbi para uso com Herrman-Ney
void Viterbi
(
  int arvore, // arvore a ser processada
  int no, // no a ser processado
  struct no vocab, // dados do no a ser processado
  int dur_ant, // duracao do caminho otimo ate o ultimo estado do fone anterior
  double log_ant, // verossimilhanca acumulada ate o ultimo estado do fone anterior
  struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  double pt, // probabilidade de transicao do no anterior para este no
  double **log_local, // contribuicao local de cada arco p/ a verossimilhanca total
  double *max_like1 // max veross. neste instante de tempo
);

// Implementa a poda de caminhos via Beam Search
void BeamSearch
(
  int ****ativop1, // indica quais niveis intra-arvore estao ativos
  int ***lista, // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int *profundidades, // indica a profundidade do no mais profundo em cada arvore
  int **quantos, // indica quantos elementos existem em cada nivel de cada uma das listas
  int n_arvores, // numero de arvores
  struct no ****vocab1, // armazena as palavras do vocabulario em arvores
  struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  double limiar, // limiar de poda
  double **glike, // verossimilhanca dos nos gramaticos
  int **ativog1, // indica se o no gramatico esta ativo ou nao
  int numberOfLevels, // numero maximo de palavras da locucao
  int t
);

#endif /*HERRMANNEY_H_*/
