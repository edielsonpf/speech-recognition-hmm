#ifndef ONESTEP_H_
#define ONESTEP_H_

#include "estruturas.h"
#include "wordDurationModel.h"

void oneStep
(
	struct locucao x, // locucao a ser reconhecida
  char *frase, // frase reconhecida
  struct configBusca searchConfig, // search algorithm configurations
	struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
  struct bigram *gramatica, // armazena gramatica bigram de palavras
  int n_termos, // numero de pares da gramatica bigram
  struct vocab vocabulario // estrutura que armazena o vocabulario
);
	
void viterbi
(
	int palavra, // palavra sob analise
	double glike, // verossimilhanca do no gramatico
	int word, // armazena as palavras vencedoras a cada nivel e a cada instante de tempo
	int nivel, // nivel da busca
	double prob_ant, // prob de transicao do nivel anterior
	double *like, // verossimilhanca do no intra-palavra
	struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
	struct configBusca Configuracao, // configuracoes de busca selecionadas pelo usuario
	int *elapse, // tempo decorrido pelo caminho desde que entrou no no gramatico
	int t, // instante atual
	struct locucao x, // locucao a ser reconhecida
	double *max_like, // verossimilhanca maxima no instante t
	struct vocab vocabulario, // dados do vocabulario
	long n_termos, // numero de pares da gramatica bigram de palavras
	struct bigram *gramatica, // pares da gramatica bigram
	double ***log_local // contribuicao local de cada arco p/ a verossinilhanca total
);

// Elimina estados com verossimilhanca menor que o limiar
void beamSearch
(
	int **ativog1, // indica se o no gramatico esta ativo ou nao
	struct configBusca searchConfig, // opcoes de configuracao do modo de busca
	double ****like1, // verossimilhanca acumulada do melhor caminho para o no i
	double max_like, // verossimilhanca maxima de todos os estados em um determinado instante de tempo
	struct modelo_HMM HMM, // modelos HMM das subunidades foneticas
	struct vocab vocabulario, // estrutura que armazena o vocabulario
	double *glike // verossimilhanca dos nos gramaticos
);


#endif /*ONESTEP_H_*/
