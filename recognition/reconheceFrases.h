#ifndef RECONHECEFRASES_H_
#define RECONHECEFRASES_H_

void reconheceFrases
(
	struct configBusca searchConfig, // search algorithm configurations	
	char *utterancesList,
	struct modelo_HMM HMM, // acoustic models
	struct bigram *grammar, // stores the bigram information
	long nTerms,	// number of entries in the grammar
	struct vocab vocabulary, // stores the vocabulary
	char *name_out
);

#endif /*RECONHECEFRASES_H_*/
