#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"estruturas.h"
#include"mystrfun.h"
#include"loadVocab.h"
//-------------------------------------------------------------------------
// Programa que compila descricao do vocabulario de palavras.
// A entrada deve ser um arquivo texto seguindo o padrao:
//
// *fonemas
// k
// ...
// u
// *fim
// *classes
// sub
// ...
// prep
// *fim
// *vocab
// casa/k a z a/<media>/<desvio padrao>/classes
// ...
// *fim
//
// O arquivo de sai­da tera o seguinte formato
//
//  #palavras(int)
//  #caracteres_na_palavra(char),"c","a","s","a",
//  #fones(char),<fone1>(char),<fone2>(char),...,
//  <media>(double),<desvio>(double),
//  #classes(char),<classe1>(int),<classe2>(int),...
//  ...
int CompilaVocabulario
(
  char *nome // nome do arquivo com o vocabulario
)
{
  int  i,numero,j;   // contador
  char *frase; // Armazena frase.
  long int TamArq;
  char **palavras;
  int tipo1; // verifica se a palavra e fonema (1), classe(2), ou vocabulario (3)
  int tipo2;
  char palavra[20]; // armazena a palavra do vocabulario
  int modelo[25]; // armazena a transcricao fonetica da palavra
  int  clspal[5],nclspal=0,fone=0;
  double media,desvp;
  char fonemas[5000][10];
  char classes[50][15];
  int  nfon,ncla;
  int nplv;
  char tam;
  int MaxCrtPlv=50,	  // Numero mi¡ximo de caracteres na palavra (si­mbolo).
      MaxCrtFrs=500,  // Numero mi¡ximo de caracteres por frase.
      NPF=50;		  // Numero mi¡ximo de palavras por frase.

  char sep1[]=" ",   // Separador
       sep2[]="/",   // Separador
       fim1[]="\n",  // Marcador de final de frase
       coment[]="%", // Comenti¡rio
       pchv1[]="*fonemas", // Palavras-chave.
       pchv2[]="*classes",
       pchv3[]="*vocab",
       pchv4[]="*fim";
  FILE *stream, *stream2;
  int vocab_OK = 1; // verifica se o vocabulario esta OK (1) ou nao (0)
  char NomeArquivo[256]; // var aux p/ formacao do nome do arquivo de saida
	int t; // aux var - length of strings
//---------------------------------------------------------------------------


  palavras = malloc(NPF*sizeof(char*));
  for (i=0;i<NPF;i++)
    palavras[i] = malloc(MaxCrtPlv*sizeof(char));
  frase = malloc(MaxCrtFrs*sizeof(char));

  // Acessando o arquivo de entrada
	stream = fopen(nome,"rt");
  if (stream == NULL)
	{
		printf("Error opening vocabulary file %s.\n",nome);
    exit(1);
	}
  else
  {
    fseek(stream,0,2);
    TamArq = ftell(stream);
    fseek(stream,0,0);
  }

  // Criando arquivo de sai­da
  strcpy(NomeArquivo,"vocab.vcb");
	stream2 = fopen(NomeArquivo,"wb");
  if (stream2 == NULL)
	{
		printf("Error opening file vocab.vcb\n");
    exit(1);
	}
  nfon = 0;
  fwrite(&nfon,sizeof(nfon),1,stream2);
  nplv = 0;
  fwrite(&nplv,sizeof(nplv),1,stream2);

  // Inicializa variaveis de controle e contadores
  tipo1 = 0;
  nfon = 0;
  ncla = 0;

  // Rotina principal
  while (ftell(stream) < TamArq)
  {
    // Busca frase
    if (codfrs(stream,frase,MaxCrtFrs,fim1,coment)!=1)
    {
      printf("Erro! Mi¡ximo de caracteres na frase excedido!");
      exit(1);
    }
		
		// removing control characters (\r \n) from the end of string
		t = strlen(frase);
		while ((t>0)&&(frase[t-1]=='\r' || frase[t-1]=='\n'))
		{
			frase[t-1] = '\0';
			t--;
		}
				
    for (i=0;i<20;i++)
      palavra[i]='\0';

    // Busca frase
    if (strlen(frase)==0)
      continue;

    // Busca palavras na frase
    if (codpal(frase,palavras,&numero,sep1,sep2,MaxCrtPlv,NPF)!=1)
    {
      printf("Erro! Mi¡ximo de caracteres na palavra excedido!");
      exit(1);
    }

    // Identifica palavras chave: 1(fonemas), 2(classes), 3(vocabuli¡rio).
    // Varii¡vel de controle "tipo1".
    if (strcmp(palavras[0],pchv1)==0)
    {
      tipo1=1;
      continue;
    }
    if (strcmp(palavras[0],pchv2)==0)
    {
      tipo1=2;
      continue;
    }
    if (strcmp(palavras[0],pchv3)==0)
    {
      tipo1=3;
      continue;
    }
    if (strcmp(palavras[0],pchv4)==0)
    {
      tipo1=0;
      continue;
    }

    tipo2=0;
    for (i=0;i<numero;i++)
    {

      // Verifica se palavra atual i© o separador de tipo "/".
      // A varii¡vel de controle i© "tipo2".
      if (palavras[i][0]==*sep2)
      {
        tipo2++;
        continue;
      }

      // Efetua acao correspondente a construcao da tabela dos fonemas(1),
      // das classes(2) e do vocabulario(3), controlada por "tipo1".
      if (tipo1==1)
        strcpy(fonemas[nfon++],palavras[0]);
      if (tipo1==2)
        strcpy(classes[ncla++],palavras[0]);
      if (tipo1==3)
      {

        // Armazena as palavras do vocabulario
        if (tipo2==0)
        {
          fone=0;
          nclspal=0;
          media=0;
          desvp=0;
          if (strlen(palavra)>0)
            strcat(palavra," ");
          strcat(palavra,palavras[i]);
        }

        // Armazena fonemas do modelo de cada palavra
        if (tipo2==1)
        {
          for (j=0;j<nfon;j++)
            if (strcmp(palavras[i],fonemas[j])==0)
            {
              modelo[fone++] = j;
              break;
            }
          if (strcmp(palavras[i],fonemas[j])!=0)
          {
						printf("Fone %s da palavra %s desconhecido!\n",palavras[i],palavra);
            vocab_OK = 0;
          }
        }

        // Armazena a mi©dia e o desvio padrao
        if (tipo2==2)
          media=atof(palavras[i]);
        if (tipo2==3)
          desvp=atof(palavras[i]);

        // Codifica as classes de cada palavra
        if (tipo2==4)
        {
          for (j=0;j<ncla;j++)
            if (strcmp(palavras[i],classes[j])==0)
            {
              clspal[nclspal++]=j;
	          break;
            }
          if (strcmp(palavras[i],classes[j])!=0)
          {
						printf("Nao encontrei classe da palavra %d: %s %s",(int)(nplv-1),palavra,palavras[i]);
	        exit(1);
          }
        }
      }
    }

    // Grava os dados no arquivo de sai­da
    if (tipo1==3)
    {
      nplv++;
      tam=(char)strlen(palavra);
      fwrite(&tam,sizeof(char),1,stream2);
      fwrite(&palavra,strlen(palavra),1,stream2);
      if (strlen(palavra)==0)
        printf("Palavra %d i© vazia!\n",nplv);

      tam=(char)fone;
      fwrite(&tam,sizeof(char),1,stream2);
      for (i=0;i<fone;i++)
        fwrite(&modelo[i],sizeof(int),1,stream2);
      fwrite(&media,sizeof(double),1,stream2);
      fwrite(&desvp,sizeof(double),1,stream2);
      if (nclspal==0)
        printf("Palavra %s (%d) nao tem classe definida!\n",palavra,nplv);
      tam=(char)nclspal;
      fwrite(&tam,sizeof(char),1,stream2);
      for (i=0;i<nclspal;i++)
        fwrite(&clspal[i],sizeof(int),1,stream2);
      strcpy(palavra,"");
    }
  }

  // Retorna i  posicao inicial do arquivo de sai­da e grava numero de palavras.
  fseek(stream2,0,SEEK_SET);
  fwrite(&nfon,sizeof(nfon),1,stream2);
  fwrite(&nplv,sizeof(nplv),1,stream2);

  // Fechando arquivos
  fclose(stream);
  fclose(stream2);

  // Desalocando ponteiros
  for (i=0;i<NPF;i++)
    free(palavras[i]);
  free(palavras);
  free(frase);
  
	// Valor de retorno da funcao
  return vocab_OK;
}

//------------------------------------------------------------------------------
// Funcao que carrega as subunidades, os modelos foneticos, modelos de duracao e
// classes das palavras do vocabulario codificado pela funcao 'compila_vocabulario()'
void CarregaVocabulario
(
  char *vocabFile, // nome do arquivo de vocabulario compilado
  struct vocab *vocabulario1 // estrutura que armazena o vocabulario
)
{
  FILE *stream;
  int  i; // contadores
	char TamPal,Nfon,NumCla;
  double media,desvp;
  struct vocab vocabu;

	// Compilando vocabulario. O resultado i© armazenado no arquivo vocab.vcb
	CompilaVocabulario(vocabFile);
	
  // Abrindo arquivo
	stream = fopen("vocab.vcb","rb");
  if (stream == NULL)
	{
		printf("Error opening file vocab.vcb.\n");
    exit(1);
	}
  fread(&vocabu.n_fones,sizeof(int),1,stream);
  fread(&vocabu.n_palavras,sizeof(int),1,stream);

	vocabu.Mws  = malloc(vocabu.n_palavras*sizeof(int *));
	vocabu.Mw   = malloc(vocabu.n_palavras*sizeof(int *));
	vocabu.Mdur = malloc(vocabu.n_palavras*sizeof(int));
	vocabu.Ddur = malloc(vocabu.n_palavras*sizeof(int));
	vocabu.classes = malloc(vocabu.n_palavras*sizeof(int *));

  for (i=0;i<vocabu.n_palavras;i++)
  {
		fread(&TamPal,sizeof(char),1,stream);
		vocabu.Mws[i] = malloc((TamPal+1)*sizeof(char));
		fread(vocabu.Mws[i],sizeof(char),TamPal,stream);
		vocabu.Mws[i][(int)TamPal]=0;
		fread(&Nfon,sizeof(char),1,stream);
		vocabu.Mw[i] = malloc((Nfon+1)*sizeof(int));
		fread(&vocabu.Mw[i][1],sizeof(int),Nfon,stream);
		vocabu.Mw[i][0]=Nfon;
		fread(&media,sizeof(double),1,stream);
		fread(&desvp,sizeof(double),1,stream);
		vocabu.Mdur[i]=(int)media;

		vocabu.Ddur[i]=(int)desvp;
		fread(&NumCla,sizeof(char),1,stream);
		vocabu.classes[i] = malloc((NumCla+1)*sizeof(int));
		fread(&vocabu.classes[i][1],sizeof(int),NumCla,stream);
		vocabu.classes[i][0]=NumCla;
  }
  fclose(stream);

  *vocabulario1 = vocabu;
}
//------------------------------------------------------------------------------
// Funcao que desaloca o ponteiro vocabulario
void DesalocaVocabulario
(
  struct vocab *vocabulario1 // estrutura que armazena o vocabulario
)
{
  int register i; // contador
  struct vocab vocabu; // estrutura que armazena o vocabulario

  vocabu = *vocabulario1;

  for (i=0;i<vocabu.n_palavras;i++)
  {
    free(vocabu.Mw[i]);
    free(vocabu.Mws[i]);
    free(vocabu.classes[i]);
  }
  free(vocabu.Mws);
  free(vocabu.Mw);
  free(vocabu.Mdur);
  free(vocabu.Ddur);
  free(vocabu.classes);

  vocabu.classes = NULL;
  vocabu.Ddur = NULL;
  vocabu.Mdur = NULL;
  vocabu.Mw = NULL;
  vocabu.Mws = NULL;

  *vocabulario1 = vocabu;
}
