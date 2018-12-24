#ifndef FORBACK_H_
#define FORBACK_H_


// Funcao que implementa o algoritmo Forward - Backward para o treinamento dos HMM's
void forwardBackward
(
	char **subunits,        // list of phonetic subunits
	struct config trainConfig, // training configuration settings  
	struct modelo_HMM *HMM // modelos HMM das subunidades foneticas
)
;

//---------------------------------------------------------------------------
// Funcao que implementa o passo FORWARD do algoritmo Baum-Welch
double Foward
(
  double **a,         // matriz de transicao do modelo HMM da locucao
  double **alfa,      // variavel FORWARD
  struct mistura **b, // matriz de emissao do modelo HMM da locucao
  double ***bj,       // probabilidade de cada simbolo da locucao, dado o modelo
  double *c,          // fator de escala
  int comp,           // numero de subunidades na locucao
  struct modelo_HMM *ModFones, // modelos HMM das subunidades foneticas
  long tamanho        // comprimento da sequencia de entrada
);

//---------------------------------------------------------------------------
// Funcao que implementa o passo BACKWARD do algoritmo Baum-Welch
void Backward
(
  double **a,         // matriz de transicao do modelo HMM da locucao
  struct mistura **b, // matriz de emissao do modelo HMM da locucao
  double **beta,      // variavel BACKWARD
  double ***bj,       // probabilidade de cada simbolo da locucao, dado o modelo
  double *c,          // fator de escala
  int comp,           // numero de subunidades na locucao
  struct modelo_HMM *ModFones, // modelos HMM das subunidades foneticas
  long tamanho        // numero de quadros na locucao
);

//------------------------------------------------------------------------------
// Funcao que soma as contagens de cada locucao para atualizacao das matrizes de
// transicao e emissao dos modelos HMM das subunidades

void SomaContagens
(
  double **a,                // matriz de transicao do modelo HMM da locucao
  double **alfa,             // variavel FORWARD
  struct mistura **b,        // matriz de emissao do modelo HMM da locucao
  double **beta,             // variavel BACKWARD
  double ***bj,              // probabilidade de cada simbolo da locucao, dado o modelo
  double *c,                 // fator de escala
  int comp,                  // numero de subunidades na locucao
  double **den_A,            // var aux p/ atualizacao da matriz A
  double ****den_B_c,        // var aux para atualizacao das probs de emissao dos modelos das subunidades
  int *modelo,               // modelo fonetico da locucao de treinamento
  struct modelo_HMM *ModFones, // modelos HMM das subunidades foneticas
  double ****Nm,             // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
  double ***num_A,           // var aux p/ atualizacao da matriz A
  struct mistura ***num_B,   // var aux para atualizacao das probs de emissao dos modelos das subunidades
  struct locucao x           // locucao parametrizada
);

//------------------------------------------------------------------------------
// Funcao que atualiza as matrizes de transicao das subunidades depois de cada
// epoca de treinamento
void Atualiza_A
(
  double **den_A,    // var aux p/ atualizacao da matriz A
  struct modelo_HMM *ModFones, // modelos HMM das subunidades foneticas
  double ***num_A    // var aux p/ atualizacao da matriz A

);

//---------------------------------------------------------------------------
// Funcao que atualiza as matrizes de emissao das subunidades depois de cada
// epoca de treinamento
void Atualiza_B
(
  double ****den_B_c,       // var aux p/ atualizacao de B
  struct modelo_HMM *ModFones, // modelos HMM das subunidades foneticas
  struct mistura ***num_B  // var aux p/ atualizacao de B
);

//---------------------------------------------------------------------------
// Funcao que calcula a P(O/M) media para todas as locucoes

double CalculaProbMedia
(
	char **subunits,        // list of phonetic subunits
	struct config trainConfig, // training configuration settings  
	struct modelo_HMM *HMM, // modelos HMM das subunidades foneticas
	int nArquivosConv // number of files for convergence testing
);

//------------------------------------------------------------------------------
// Funcao que implementa o algoritmo de Viterbi
double Viterbi
(
  double **a,        // matriz de transicao
  double ***bj,      // probabilidade de cada simbolo da locucao, dado o modelo
  int comp,          // numero de subunidades na locucao
  struct modelo_HMM ModFones, // modelos HMM das subunidades foneticas
  int tamanho,       // comprimento da sequencia de entrada
  struct v_psi vpsi         // estrutura para armazenar o caminho otimo
);

//------------------------------------------------------------------------------
// Atribui um valor minimo (definido pela variavel MINIMO em HMM.cpp) p/ os
// coeficientes de ponderacao e variancias das gaussianas.
int ValoresMinimos
(
  struct modelo_HMM* ModFones // modelos HMM das subunidades foneticas
);

//------------------------------------------------------------------------------
// Verifica se a soma das probabilidades de transicao e igual a 1 e
// se a soma dos coeficientes de ponderacao das gaussianas e igual a 1
// retorna true se a verificacao nao encontrou erro
int TestaHMM
(
  struct modelo_HMM* ModFones
);

#endif /*FORBACK_H_*/
