#ifndef LOADHMM_H_
#define LOADHMM_H_

// Le informacoes do arquivo .inf referente aos modelos HMM
int loadHMMInfo
(
  char *nome_HMM, // nome do arquivo com os modelos HMM das subunidades
  struct modelo_HMM *ModFones// modelos HMM das subunidades foneticas
);
	
//-------------------------------------------------------------------------
// Carrega modelos HMM das subunidades foneticas
void loadHMM
(
  char *nome_HMM, // nome do arquivo com os modelos HMM das subunidades
  struct modelo_HMM *ModFones// modelos HMM das subunidades foneticas
);

//-------------------------------------------------------------------------
// Funcao que desaloca os ponteiros de ModFones
void freeHMM
(
  struct modelo_HMM *ModFones// modelos HMM das subunidades foneticas
);

#endif /*LOADHMM_H_*/
