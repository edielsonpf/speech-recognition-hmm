#ifndef GAUSFUNC_H_
#define GAUSFUNC_H_

// Funcao que implementa uma fdp gaussiana multidimensional
double Gauss
(
	double *media, // vetor com as medias
	int ordem,     // ordem do vetor de parametros
	double *var,   // vetor com as variancias
	double *x      // vetor de parametros
);

// Funcao que calcula a probabilidade de um simbolo para uma mistura de gaussianas
void ProbSimbolo
(
	struct mistura **b, // matriz de emissao para os modelos das locucoes
	double ****bj1,     // probabilidade de cada simbolo da locucao, dado o modelo
	int *comprimentos,  // numero de subunidades foneticas das transcricoes de cada locucao
	int frase,          // conta as locucoes de treinamento
	int n_estados,      // numero de estados para cada modelo HMM de subunidade
	int *n_gaussianas,  // numero de gaussianas por mistura
	double *****Nm1,    // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
	int n_par,          // numero de parametros de treinamento
	int *ordem,         // dimensao do vetor para cada um dos parâmetros
	int tamanho,        // numero de quadros da locucao de treinamento
	double ***x         // vetores de parametros de treinamento
);

	
// Funcao que calcula a probabilidade de um simbolo para uma mistura de gaussianas
// alternativa para o aplicativo de Segmentação (Ynoguti, 19/09/2003)
void ProbSimbolo1
(
	struct mistura **b, // matriz de emissao para os modelos das locucoes
	double ****bj1,     // probabilidade de cada simbolo da locucao, dado o modelo
	int comprimento,  // numero de subunidades foneticas das transcricoes de cada locucao
	int n_estados,      // numero de estados para cada modelo HMM de subunidade
	int *n_gaussianas,  // numero de gaussianas por mistura
	double *****Nm1,    // contribuicao de cada gaussiana da mistura na formacao da probabilidade de simbolo
	int n_par,          // numero de parametros de treinamento
	int *ordem,         // dimensao do vetor para cada um dos parâmetros
	int tamanho,        // numero de quadros da locucao de treinamento
	double ***x         // vetores de parametros de treinamento
);
/*
// Funcao que calcula a probabilidade de uma mistura de gaussianas
double TGauss::Prob
(
	struct mistura b, // coef. de ponderacao, media e variancia de cada gaussiana
	int n_gaussianas, // numero de gaussianas na mistura
	int ordem, // ordem do vetor x
	double *x // vetor de dados
);
*/
// Funcao que calcula a log-probabilidade de uma mistura de gaussianas
double LogProb
(
	struct mistura b, // coef. de ponderacao, media e variancia de cada gaussiana
	int n_gaussianas, // numero de gaussianas na mistura
	int ordem, // ordem do vetor x
	double *x // vetor de dados
);

// Funcao que calcula contribuicoes locais de cada arco, p/ cada subunidade fonetica
void CalcContrLocal
(
	int fone, // fone sob analise
	double ****log_local1, // contribuicao local de cada arco p/ a verossinilhanca total
	struct modelo_HMM ModFones, // modelos HMM das subunidades foneticas
	int t, // frame sob analise
	struct locucao x // locucao a ser reconhecida
);

#endif /*GAUSFUNC_H_*/
