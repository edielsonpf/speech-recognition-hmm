//---------------------------------------------------------------------------
// Programa para gerar codebook para quantizacao vetorial
// Utiliza o algoritmo LBG em sua versao splitting
//---------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "lbg.h"
//------------------------------------------------------------------------------
// Funcao que calcula 'n_codevectors' codebooks a partir de 'n_vetores' vetores
// de treinamento. Os vetores sao de ordem 'ordem', e sao armazenados em 'x'.
// A funcao retorna os 'n_codevectors' codebooks calculados
void lbg
(
  double ***codebook1, // armazena o codebook
  int n_codevectors,   // numero de codebooks desejado para o quantizador
  int n_vetores,       // numero de vetores exemplo para o calculo dos codebooks
  int ordem,           // ordem dos vetores exemplo
  double **x           // vetores exemplo para o calculo dos codebooks
)
{
  // Declaracao de variaveis locais
  double **codebook; // armazena o codebook
  double **codebook2; // codebook gerado com 2^n vetores
  double epsilon=0.005; // perturbacao (ver demonstracao do algoritmo)
  int expoente; // var aux para calcular num_vet
  int register i,j; // contadores
  int *numero; // armazena o numero de vetores em cada particao
  int num_vet; // menor numero maior ou igual a n_codevectors que e tambem potencia de 2.
  int n_vet; // tamanho atual do codebook
  double *p1,*p2; // variaveis auxiliares no calculo de distancias euclidianas
  double **somas; // variavel auxiliar

//******************************************************************************

  codebook = *codebook1;

  // Determinando menor numero maior ou igual a n_codevectors, e que seja potencia
  // de 2 (pois na versao splitting, o LBG calcula apenas codebooks com numero
  // de vetores que seja potencia de 2)
  expoente = 0;
  num_vet = pow(2,expoente);
  while (num_vet < n_codevectors)
  {
    expoente++;
    num_vet = pow(2,expoente);
  }

  // Alocando memoria p/ os ponteiros

  // Codebook
  codebook2 = malloc(sizeof(double *)*num_vet);
  for (i=0;i<num_vet;i++)
      codebook2[i] = malloc(sizeof(double)*ordem);

  // Matriz auxiliar no calculo do codebook
  somas = malloc(sizeof(double *)*num_vet);
  for (i=0;i<num_vet;i++)
    somas[i] = malloc(sizeof(double)*ordem);

  // Var's auxiliares no calculo da dist. Euclidiana
  p1 = malloc(sizeof(double)*ordem);
  p2 = malloc(sizeof(double)*ordem);

  // Numero de vetores em cada particao
  numero = malloc(sizeof(int)*num_vet);

  // Calculando centroide de todos os vetores (codebook inicial)
  calcula_centroide(&codebook2[0],n_vetores,ordem,x);
  n_vet = 1;

  // Calculando demais codebooks
  while (n_vet < num_vet)
  {

    // Splitting matriz codebook
	for (i=0; i < n_vet; i++)
	  for (j=0; j < ordem; j++)
	    somas[i][j] = codebook2[i][j];

	for (i=1; i<=n_vet; i++)
	  for (j=0; j < ordem; j++)
      {
	    codebook2[2*(n_vet-i)][j] = somas[n_vet-i][j] - epsilon;
		codebook2[2*(n_vet-i)+1][j] = somas[n_vet-i][j] + epsilon;
	  }

	// Proximo codebook
	n_vet *= 2;

    // Calculando novo codebook
    novo_codebook(&codebook2,&numero,n_vet,n_vetores,ordem,&somas,x);
  }

//------------------------------------------------------------------------------
  // Temos agora um codebook otimizado com num_vet vetores codigo
  // Se (n_codevectors < num_vet), devemos retirar alguns vetores codigo do codebook
  // O criterio de retirada e o numero de vetores de treinamento associados a cada
  // um dos vetores codigo. Aqueles com menor numero de vetores de treinamento
  // associados sao eliminados primeiro.

  if (n_codevectors < num_vet)
  {
    // Eliminando vetores codigo de 'codebook2' para gerar 'codebook' com 'n_codevectors'
    // vetores codigo
    ordena_elimina(codebook2,&codebook,n_codevectors,numero,num_vet,ordem);

    // Atualizando codebook
    novo_codebook(&codebook,&numero,n_codevectors,n_vetores,ordem,&somas,x);
  }
  else
    for (i=0;i<n_codevectors;i++)
      for (j=0;j<ordem;j++)
        codebook[i][j] = codebook2[i][j];
  
//******************************************************************************

  // Desalocando ponteiros
  free(p1);
  free(p2);

  for (i=0;i<num_vet;i++)
  {
    free(somas[i]);
    free(codebook2[i]);
  }
  free(somas);
  free(codebook2);

  free(numero);

  // Retornando codebook calculado
  *codebook1 = codebook;
}

//******************************************************************************
// Funcao que calcula a distancia euclidiana entre dois vetores
double deucl
(
  double *p1, // vetor 1
  double *p2, // vetor 2
  int ordem   // ordem dos vetores
)
{
  int i;
  double soma;

  soma = 0.0;
  for (i=0; i<ordem; i++)
	soma += (p1[i]-p2[i])*(p1[i]-p2[i]);

  return sqrt(soma);
}

//******************************************************************************
// Funcao que retorna o centroide de um conjunto de vetores
void calcula_centroide
(
  double **centroide1, // centroide
  int n_vetores,            // numero de vetores exemplo para o calculo dos codebooks
  int ordem,                // ordem dos vetores exemplo
  double **x                // vetores exemplo para o calculo dos codebooks
)
{
  // Declaracao das variaveis locais
  double *centroide; // centroide dos vetores exemplo
  int i,j; // contadores

  centroide = *centroide1;

  // Inicializando variavel que ira conter o centroide
  for (i=0;i<ordem;i++)
	centroide[i]=0.0;

  // Somando os quadros
  for (i=0;i<n_vetores;i++)
    for (j=0;j<ordem;j++)
      centroide[j] += x[i][j];

  for (i=0;i<ordem;i++)
    centroide[i] /= (float)n_vetores;

  *centroide1 = centroide;

}
//******************************************************************************
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
)
{
  double **codebook_out; // codebook reordenado com apenas 'n_codevectors' vetores codigo
  int done = 0; // flag auxilliar para fim de processamento
  register int i,j,k; // contadores
  int *index; // indices dos vetores codigo do codebook (var aux para ordenacao)

  codebook_out = *codebook_out1;

  // Alocando memoria
  index = malloc(sizeof(int)*num_vet);

  for (i=0;i<num_vet;i++)
    index[i] = i;

  while (!done)
  {
    done = 1;
    for (i=0;i<(num_vet-1);i++)
    {
      j = index[i];
      k = index[i+1];
      if (numero[j] < numero[k])
      {
        index[i] = k;
        index[i+1] = j;
        done = 0;
      }
    }
  }

  // Salvando 'n_codevectors' em 'codebook_out'
  for (i=0;i<n_codevectors;i++)
    for (j=0;j<ordem;j++)
      codebook_out[i][j] = codebook_in[index[i]][j];

  // Desalocando ponteiro
  free(index);

  // Codebook reduzido
  *codebook_out1 = codebook_out;
}

//------------------------------------------------------------------------------
// Realiza operacoes para o calculo de um novo codebook, depois do splitting
void novo_codebook
(
  double ***codebook1, // armazena o codebook a ser processado
  int **numero1, // armazena o numero de vetores em cada particao
  int n_vet, // tamanho atual do codebook
  int n_vetores,   // numero de vetores exemplo para o calculo dos codebooks
  int ordem,       // ordem dos vetores exemplo
  double ***somas1, // variavel auxiliar
  double **x       // vetores exemplo para o calculo dos codebooks
)
{
  double **codebook; // // armazena o codebook a ser processado
  double distorcao; // verifica se o algoritmo convergiu
  double d_media[2]; // distancia media ao realizar a quantizacao com o codebook
  int register i,j,n; // contadores
  int *numero; // armazena o numero de vetores em cada particao
  double **somas; // variavel auxiliar

  numero = *numero1;
  somas = *somas1;
  codebook = *codebook1;

  // Inicializando vetor que armazena a distancia media
  d_media[0] = d_media[1] = 1e20;

  // Inicializando a variavel numero (conta quantos vetores em cada particao)
  for (i=0;i<n_vet;i++)
    numero[i] = 0;

  // Inicializando a variavel somas (var aux p/ calcular o centroide de cada particao)
  for (i=0;i<ordem;i++)
    for (j=0;j<n_vet;j++)
      somas[j][i] = 0.0;

  // Associando vetores a particoes
  associa(codebook,&d_media[1],&numero,n_vet,n_vetores,ordem,&somas,x);

  distorcao = 1;
  n=0;
  while (distorcao > 0.005)
  {
    // Calculando o centroide de cada particao
    for (i=0;i<n_vet;i++)
      for (j=0;j<ordem;j++)
        codebook[i][j] = somas[i][j]/numero[i];

    // Inicializando a variavel numero (conta quantos vetores em cada particao)
    for (i=0;i<n_vet;i++)
      numero[i] = 0;

    // Inicializando a variavel somas (var aux p/ calcular o centroide de cada particao)
    for (i=0;i<ordem;i++)
      for (j=0;j<n_vet;j++)
        somas[j][i] = 0.0;

    // Calculando a distancia media com este codebook
    d_media[0] = d_media[1];
    d_media[1] = 0.0;

    // Associando vetores a particoes
    associa(codebook,&d_media[1],&numero,n_vet,n_vetores,ordem,&somas,x);

    d_media[1] /= (float)n_vetores;

    // Calculando a distorcao
    if (d_media[1] != 0)
      distorcao = (d_media[0]-d_media[1])/d_media[1];
    else
      distorcao = 0; // se a distancia media e zero, temos um quantizador otimo
    if (distorcao < 0.0)
      distorcao=1;
    n++;
  }

  // Valores de retorno da funcao
  *codebook1 = codebook;
  *numero1 = numero;
  *somas1 = somas;
}

//------------------------------------------------------------------------------
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
)
{
  //double distancia; // distancia euclidiana entre dois vetores
  double d_media; // distancia media com este codebook
  int register i,l; // contadores
  double minimo; // menor distancia
  int *numero; // armazena o numero de vetores em cada particao
  int qual; // vetor codigo associado ao vetor de treinamento atual
  double **somas; // variavel auxiliar

  somas = *somas1;
  numero = *numero1;
  d_media = *d_media1;
  
  // Calculando distancias dos parametros aos centroides
  for (i=0;i<n_vetores;i++)
  {
    // Verificando vetor codigo mais proximo de x[i]
    qual = mais_proximo(codebook,&minimo,n_vet,ordem,x[i]);

    // Atualizando contagens dos clusters
    for (l=0;l<ordem;l++)
      somas[qual][l] += x[i][l];
    numero[qual]++;

    // Atualizando d_media
    d_media += minimo;
  }

  // Valores de retorno da funcao
  *d_media1 = d_media;
  *somas1 = somas;
  *numero1 = numero;
}

//-----------------------------------------------------------------
// Funcao que, dado um vetor 'x' e uma lista de vetores 'lista', verifica qual vetor
// de 'lista' esta mais proximo de 'x'.
// O valor de retorno indica o indice do vetor codigo mais proximo do
// vetor exemplo
int mais_proximo
(
  double **lista, // lista de vetores com a qual o vetor 'x' ira ser analisado
  double *minimo1, // distancia euclidiana entre o vetor exemplo e o vetor codigo mais proximo
  int n_codevectors, // numero de vetores no codebook
  int ordem, // ordem dos vetores
  double *x // vetor sob analise
)
{
  double distancia; // distancia euclidiana entre dois vetores
  int register l; //contadores
  double minimo; // distancia euclidiana entre o vetor exemplo e o vetor codigo mais proximo
  int qual; // vetor codigo mais proximo do vetor exemplo

  minimo = *minimo1;

  minimo = 1e20;

  for (l=0;l<n_codevectors;l++)
  {
    distancia = deucl(x,lista[l],ordem);
    if (distancia < minimo)
    {
      minimo = distancia;
      qual = l;
    }
  }

  *minimo1 = minimo;
  return qual;
}


