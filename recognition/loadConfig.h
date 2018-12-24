#ifndef LOADCONFIG_H_
#define LOADCONFIG_H_

int loadConfig
(
  char *configFile,
  //~ int *searchAlgorithm, // 0: Level Building, 1: One Step, 2: Herrmann Ney
  //~ int *useBeamSearch, // 0 yes, 1 no
	//~ float *beamSearchThreshold, // threshold for node pruning
	//~ int *numberOfLevels, // maximum number of words per utterance
	//~ int *useWordDurationModel, // use (0) or not (1) the word duration model
	//~ int *useGrammar, // use (0) or not (1) the language model
	struct configBusca *searchConfig,
	char *grammarFile, // file with the language model
	char *vocabFile, // file with the vocabulary
	char *HMMFile, // file with the HMM parameters
	int *nUtterances, // number of utterances to be recognized
	char *utterancesList // file with the list of utterances to be recognized
);

#endif /*LOADCONFIG_H_*/
