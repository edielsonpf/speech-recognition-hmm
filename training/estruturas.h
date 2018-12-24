#ifndef ESTRUTURAS_H_
#define ESTRUTURAS_H_

//---------------------------------------------------------------------------
// Estrutura para armazenar as opcoes de configuracao do modo de busca
struct config
{
	int initializationMode; // 0 - m = 1.0; var = 1.0; c = 1.0
                          // 1 - Segmental K-Means
                          // 2 - Modelos pré-treinados
  	char HMMFile[512]; // name of file that will store the trained HMMs
  	char subunitsFile[512]; // file with phonetic subunits listing
  	char transcriptionsFiles[512]; // path of transcriptions files 
	char utterancesFiles[512]; // path of transcriptions files
  	int nUtterances; // number of training utterances
  	char preTrained[512]; // file with pre-trained HMMs
};

// Estrutura que armazena os dados das misturas referentes as fdp's de emissao de simbolos
struct mistura
{
	double* c; // coeficientes de ponderacao para cada gaussiana
	double** m; // medias para cada uma das gaussianas
  	double** v; // variancias para cada uma das gaussianas
};


// Estrutura para armazenar os modelos HMM das subunidades acusticas
struct modelo_HMM
{
  	int* delta;          // numero de janelas adjacentes utilizadas para calculo de parametros delta
  	int n_estados;       // numero de estados para cada subunidade fonetica
  	int n_fones;         // numero de subunidades acusticas no modelo
  	int* n_gaussianas;   // numero de gaussianas na mistura, para cada um dos parametros
  	int n_par;           // numero de parametros utilizados no reconhecimento
  	char** parametros;   // descricao dos parametros utilizados no reconhecimento (p. ex. 'mel')
  	int* ordem;          // ordem dos parametros
  	double*** A;         // matrizes de transicao
  	struct mistura*** B; // matrizes de emissao
  	int salto_maximo;    // salto maximo permitido entre estados no modelo HMM das subunidades
	int usaPCA;          // usa (1) ou nao (0) a PCA para reducao da dimensao dos vetores acusticos
	char nomePCA[512];   // nome do arquivo com a matriz de transformacao PCA
	int ordemRed;        // nova dimensao dos vetores acusticos, depois da transformacao PCA
};

//Estrutura para armazenar as locucoes a serem reconhecidas
struct locucao
{
  int n_quadros; // numero de quadros com o qual foi parametrizada a locucao
  double*** par;  // locucao parametrizada (n_par x n_janelas x ordem)
};

// Estrutura que armazena as verossimilhancas e os backpointers durante o processamento
// pelo algoritmo de Viterbi (versao p/ o Level Building)
struct v_psi
{
  double** v;
  int**   psi;
};

#endif /*ESTRUTURAS_H_*/
