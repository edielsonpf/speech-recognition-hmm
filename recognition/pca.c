#include <stdio.h>
#include <stdlib.h>
//---------------------------------------------------------------------------
// Realiza reducao da dimensao do vetor de parametros atraves da Analise de
// Componente Principal (por enquanto usado somente para parametros mdd)
void pca
(
  char* NomePCA,  // nome da matriz de transformacao gerada pelo programa KL.exe
  int OrdemRed,   // ordem desejada para o vetor de parametros apos a reducao pela PCA
  int** ordem1,   // ordem dos vetores de cada um dos parametros
  int n_par,      // numero de parametros a serem calculados
  int Nwindow,    //numero de quadros com o qual foi analisado o sinal
  double ****par1 // sinal parametrizado
)
{
  FILE* arquivo; // handle p/ leitura de arquivo
  int register i,j,k; // contadores
  double **KL; // matriz de transformacao
  double *media; // vetor media dos parametros
  int* ordem;    // ordem dos vetores de cada um dos parametros
  double*** par; // sinal parametrizado
  double** ParRed; // sinal parametrizado reduzido pela PCA
  double* S; // raiz quadrada da diagonal principal da matriz de covariancia
  par = *par1;
  ordem = *ordem1;

  // Alocando memoria
	media = malloc(sizeof(double)*ordem[0]);

	S = malloc(sizeof(double)*ordem[0]);

	KL = malloc(sizeof(double *)*ordem[0]);
	for (i=0;i<ordem[0];i++)
		KL[i] = malloc(sizeof(double)*ordem[0]);

	ParRed = malloc(sizeof(double *)*Nwindow);
	for (i=0;i<Nwindow;i++)
		ParRed[i] = malloc(sizeof(double)*OrdemRed);

  // Inicializando matriz que ira armazenar os parametros reduzidos
  for (i=0;i<Nwindow;i++)
    for (j=0;j<OrdemRed;j++)
      ParRed[i][j] = 0.0;

  // Lendo dados
	arquivo = fopen(NomePCA,"rb");	
  if (arquivo == NULL)
	{
		printf("Error opening file %s.\n",NomePCA);
    exit(1);
	}
  // Media dos parametros
  for (i=0;i<ordem[0];i++)
    fread(&media[i],sizeof(double),1,arquivo);

  // Raiz quadrada da diagonal principal da matriz de covariancia
  for (i=0;i<ordem[0];i++)
    fread(&S[i],sizeof(double),1,arquivo);

  // Matriz de transformacao
  for (i=0;i<ordem[0];i++) // autovetor
    for (j=0;j<ordem[0];j++) // elem. do autovetor
      fread(&KL[i][j],sizeof(double),1,arquivo);
  fclose(arquivo);

  // Gerando parametros reduzidos
  for (i=0;i<Nwindow;i++)
    for (j=0;j<ordem[0];j++)
      for (k=0;k<OrdemRed;k++)
        ParRed[i][k] += ((par[0][i][j]-media[j])/S[j])*KL[k][j];

  // Redimensionando matriz de parametros para armazenar os parametros reduzidos
  for (i=0;i<n_par;i++)
    for (j=0;j<Nwindow;j++)
      free(par[i][j]);

  for (i=0;i<n_par;i++)
    for (j=0;j<Nwindow;j++)
      par[i][j] = malloc(sizeof(double)*OrdemRed);

  // Armazenando parametros reduzidos na matriz de parametros
  for (i=0;i<Nwindow;i++)
    for (j=0;j<OrdemRed;j++)
      par[0][i][j] = ParRed[i][j];

  // Desalocando ponteiros
  free(media);
  free(S);
  for (i=0;i<ordem[0];i++)
    free(KL[i]);
  free(KL);
  for (i=0;i<Nwindow;i++)
    free(ParRed[i]);
  free(ParRed);

  // Atualizando ordem do vetor de parametros
  ordem[0] = OrdemRed;

  // Valores de retorno da funcao
  *par1 = par;
  *ordem1 = ordem;
}
