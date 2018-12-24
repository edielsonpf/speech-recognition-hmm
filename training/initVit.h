#ifndef INITVIT_H_
#define INITVIT_H_

//---------------------------------------------------------------------------
// Funcao para inicializacao dos modelos das subunidades foneticas via
// Algoritmo de Viterbi

void Inicializa_Viterbi
(
	char **subunits,        // list of phonetic subunits
	struct config trainConfig, // training configuration settings  
	struct modelo_HMM *HMM // modelos HMM das subunidades foneticas
);

#endif /*INITVIT_H_*/
