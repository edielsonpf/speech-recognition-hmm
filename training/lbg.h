#ifndef LBG_H_
#define LBG_H_

//---------------------------------------------------------------------------
// Funcao que calcula 'n_codebooks' codebooks a partir de 'n_vetores' vetores
// de treinamento. Os vetores sao de ordem 'ordem', e sao armazenados em 'x'.
// A funcao retorna os 'n_codebooks' codebooks calculados
void lbg
(
  double ***codebook1, // armazena os codebooks
  int n_codebooks, // numero de codebooks desejado para o quantizador
  int n_vetores,   // numero de vetores exemplo para o calculo dos codebooks
  int ordem,       // ordem dos vetores exemplo
  double **x       // vetores exemplo para o calculo dos codebooks
);

//---------------------------------------------------------------------------
// Funcao que calcula a distancia euclidiana entre dois vetores
double deucl
(
  double *p1, // vetor 1
  double *p2, // vetor 2
  int ordem        // ordem dos vetores
);

//---------------------------------------------------------------------------
// Funcao que retorna o centroide de um conjunto de vetores
void calcula_centroide
(
  double **centroide1, // centroide
  int n_vetores,            // numero de vetores exemplo para o calculo dos codebooks
  int ordem,                // ordem dos vetores exemplo
  double **x                // vetores exemplo para o calculo dos codebooks
);

//---------------------------------------------------------------------------
// Funcao que ordena os vetores do codebook 'codebook_in' em ordem decrescente,
// de acordo com o numero de vetores de treinamento associado a cada um dos
// vetores codigo.
// Posteriormente, retorna em 'codebook_out' apenas os 'n_codevectors' com maior
// numero de vetores de treinamento associados
void ordena_elimina
(
  double **codebook_in, // codebook a ser reordenado
  double ***codebook_out1, // codebook reordenado com apenas 'n_codevectors' vetores codigo
  int n_codevectors, // numero de vetores codigo desejado
  int *numero, // numero de vetores de treinamento associados a cada vetor codigo
  int num_vet, // numero de vetores codigo no codebook atual
  int ordem // ordem dos vetores exemplo
);

//---------------------------------------------------------------------------
// Realiza operacoes para o calculo de um novo codebook, depois do splitting
void novo_codebook
(
  double ***codebook1, // armazena o codebook a ser processado
  int **numero1, // armazena o numero de vetores em cada particao
  //int num_vet, // menor numero maior ou igual a n_codevectors que e tambem potencia de 2.
  int n_vet, // tamanho atual do codebook
  int n_vetores,   // numero de vetores exemplo para o calculo dos codebooks
  int ordem,       // ordem dos vetores exemplo
  double ***somas1, // variavel auxiliar
  double **x       // vetores exemplo para o calculo dos codebooks
);

//---------------------------------------------------------------------------
// Funcao que associa vetores de treinamento ao vetor codigo mais proximo
void associa
(
  double **codebook, // codebook
  double *d_media1, // distancia media com este codebook
  int **numero1, // armazena o numero de vetores em cada particao  
  int n_vet,         // tamanho atual do codebook
  int n_vetores,   // numero de vetores exemplo para o calculo dos codebooks
  int ordem,       // ordem dos vetores exemplo
  double ***somas1, // variavel auxiliar
  double **x       // vetores exemplo para o calculo dos codebooks
);

//-----------------------------------------------------------------
// Funcao que, dado um vetor exemplo e um codebook, verifica qual vetor
// codigo esta mais proximo do vetor exemplo.
// O valor de retorno indica o indice do vetor codigo mais proximo do
// vetor exemplo
int mais_proximo
(
  double **codebook, // codebook
  double *minimo1, // distancia euclidiana entre o vetor exemplo e o vetor codigo mais proximo
  int n_codevectors, // numero de vetores no codebook
  int ordem, // ordem dos vetores
  double *x // vetor exemplo
);

#endif /*LBG_H_*/
