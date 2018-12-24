#ifndef WORDDURATIONMODEL_H_
#define WORDDURATIONMODEL_H_

// Funcao que aplica penalizacao segundo modelo de duracao de palavras
double wordDurationModel
(
	double like, // verossimilhanca
	int elapse, // tempo de duracao (em quadros) da palavra
	int palavra, // palavra sob analise
	int *Ddur, // desvio padrao da duracao das palavras
	int *Mdur // duracao media das palavras
);

#endif /*WORDDURATIONMODEL_H_*/
