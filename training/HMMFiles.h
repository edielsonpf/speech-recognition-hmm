#ifndef HMMFILES_H_
#define HMMFILES_H_

// Le informacoes do arquivo .inf referente aos modelos HMM
int loadHMMInfo
(
  char *nome_HMM, // nome do arquivo com os modelos HMM das subunidades
  struct modelo_HMM *HMM// modelos HMM das subunidades foneticas
);
	
//-------------------------------------------------------------------------
// Carrega modelos HMM das subunidades foneticas
void loadHMM
(
  char *nome_HMM, // nome do arquivo com os modelos HMM das subunidades
  struct modelo_HMM *HMM// modelos HMM das subunidades foneticas
);

//-------------------------------------------------------------------------
// Funcao que desaloca os ponteiros de HMM
void freeHMM
(
  struct modelo_HMM *HMM// modelos HMM das subunidades foneticas
);

//-------------------------------------------------------------------------
// Funcao que guarda em arquivo modelos HMM treinados
//
// - numero de subunidades foneticas (int)
// - numero de parametros (int)
// - Para cada parametro: - numero de gaussianas na mistura (int)
//                        - ordem (int)
//                        - extensao dos arquivos de parametros (char [MAXEXT])
void SalvaHMM
(
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  char nome_saida[256], // nome do arquivo com os parâmetros HMM das sub-unidades
  int epoca // numero da epoca de treinamento
);

// Gera arquivo de informacao sobre o arquivo HMM sendo treinado
void GeraArquivoInformacao
(
  struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
  char* nome_saida // nome do arquivo de saida
);

#endif /*HMMFILES_H_*/
