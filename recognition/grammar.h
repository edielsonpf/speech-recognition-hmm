#ifndef GRAMMAR_H_
#define GRAMMAR_H_

// Aplica gramatica bigram de palavras
int Pode
(
	int palavra1, // primeira palavra
	int palavra2, // segunda palavra
	double *probabilidade, // p(segunda|primeira)
	long n_termos, // numero de pares na gramatica
	struct bigram *gram // pares da gramatica
);

#endif /*GRAMMAR_H_*/
