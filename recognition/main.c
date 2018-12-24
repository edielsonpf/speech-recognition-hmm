#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "estruturas.h"
#include "loadConfig.h"
#include "loadVocab.h"
#include "loadHMM.h"
#include "loadGrammar.h"
#include "reconheceFrases.h"

int main(int argc, char *argv[])
{
	//int searchAlgorithm; // LevelBuilding(0), OneStep(1), HerrmannNey(2)
	//int useBeamSearch; // use (0) or not (1) the beam search pruning strategy
	//float beamSearchThreshold; // threshold for node pruning
	//int numberOfLevels; // max number of words in the utterance
	//int useWordDurationModel; // use (0) or not (1) the word duration model
	//int useGrammar; // use (0) or not (1) the language model
	struct configBusca searchConfig; // search algorithm configurations
	char grammarFile[512]; // file with the language model
	char vocabFile[512]; // file that stores the vocabulary
	char HMMFile[512]; // file with the HMM parameters
	int nUtterances; // number of utterances to be recognized
	char utterancesList[512]; // file with the list of utterances to be recognized
	struct vocab vocabulary; // stores the vocabulary
	struct modelo_HMM HMM; // stores the acoustic models
	struct bigram *grammar=NULL; // stores the bigram information
	long nTerms;	// number of entries in the grammar
	char name_out[512]; // text filename with the recognized utterances
	// Testing arguments
	if (argc < 2)
	{
		puts("Usage: reco config.txt [log.txt]\n");
		puts("where:\n");
		puts("- config.txt is a textfile with configurations settings.\n");
		puts("- log.txt is an optional textfile where the recognized sentences are stored.\n");
		return (1);
	}
	
	if (argc==3)
	  strcpy(name_out,argv[2]);
	  
	// Loading configurations
	loadConfig(argv[1],&searchConfig,grammarFile,vocabFile,HMMFile,&nUtterances,utterancesList);

  //~ // Presenting configurations	
	//~ printf("Algoritmo de busca: ");
	//~ switch(searchConfig.searchAlgorithm)
	//~ {
		//~ case 0: 
			//~ printf("Level Building\n");
			//~ break;
		//~ case 1: 
			//~ printf("One Step\n");
			//~ break;
		//~ default: 
		  //~ printf("Herrmann Ney\n");	
	//~ }
	//~ printf("Número de níveis: %d\n",searchConfig.numberOfLevels);
	//~ printf("Usa Beam Search: ");
	//~ if(searchConfig.useBeamSearch == 0)
	//~ {
    //~ printf("sim\n");
		//~ printf("Limiar de poda: %f\n",searchConfig.beamSearchThreshold);
	//~ }
	//~ else	
		//~ printf("não\n");
	//~ printf("Usa modelo de duração de palavras: ");
	//~ if(searchConfig.useWordDurationModel == 0)
    //~ printf("sim\n");
	//~ else	
		//~ printf("não\n");
	//~ printf("Usa gramática: ");
	//~ if(searchConfig.useGrammar == 0)
	//~ {
    //~ printf("sim\n");
		//~ printf("Gramática: %s",grammarFile);
	//~ }
	//~ else	
		//~ printf("não\n");

	// Loading Vocabulary
	CarregaVocabulario(vocabFile,&vocabulary);
	
	// Loading HMMs
	loadHMM(HMMFile,&HMM);
  
	// Loading Grammar
	if (searchConfig.useGrammar == 1)
		loadGrammar(grammarFile,&nTerms,&grammar);

	// Recognizing utterances
    reconheceFrases(searchConfig,utterancesList,HMM,grammar,nTerms,vocabulary,name_out);
	
    // Freeing memory
	DesalocaVocabulario(&vocabulary);
    freeHMM(&HMM);
	free(grammar);
	
	return (0);
}
