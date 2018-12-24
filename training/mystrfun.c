#include <stdio.h>
#include <string.h>
//---------------------------------------------------------------------------
// Funcao que extrai palavras de uma frase
// As palavras sao devolvidas em "palav" e o numero de palavras em "numpal".
// O valor de retorno deve ser 1 para palavras corretas
int codpal
(
  char *frase, // frase a ser analisada
  char **palav, // estrutura onde serao devolvidas as palavras
  int *numpal, // numero de palavras na frase
  char sep1[], // separador
  char sep2[], // separador
  int MaxCrtPlv, // numero maximo de caracteres permitido em uma palavra
  int NPF // numero maximo de palavras na frase
)
{
  int i,pos,npl=0;

  pos=0;
  for (i=0;i<=(int)strlen(frase);i++)
  {
    if ((frase[i]!=*sep1)&&(frase[i]!=0)&&(frase[i]!=*sep2))
    {
      palav[npl][pos]=frase[i];
      if (pos<MaxCrtPlv)
        pos++;
      else
      {
        palav[npl][pos]=0;
        *numpal=npl+1;
        return(0);
      }
    }
    else
    if (frase[i]==*sep2)
    {
      if (pos!=0)
      {
        palav[npl][pos]=0;
        if (npl<NPF)
          npl++;
        else
        {
          *numpal=npl+1;
          return(0);
        }
        pos=0;
      }
      palav[npl][pos++]=frase[i];
      palav[npl][pos]=0;
      if (npl<NPF)
        npl++;
      else
      {
        *numpal=npl+1;
        return(0);
      }
      pos=0;
    }
    else
      if (pos!=0)
      {
        palav[npl][pos]=0;
        if (npl<NPF)
          npl++;
        else
        {
          *numpal=npl+1;
          return(0);
        }
        pos=0;
      }
  }
  *numpal=npl;
  return(1);
}
//------------------------------------------------------------------------------
// Funcao que extrai frases de um arquivo texto.
// A frase obtida eh retornada em "frase".
// O valor de retorno deve ser 1 para uma frase correta.
// Valor zero indica excesso de caracteres ou falta do ponto final.
int codfrs
(
  FILE *nomarq,
  char *frase,
  int MaxCrtFrs,
  char fim1[],
  char coment[]
)
{
  char buf[]="";
  int  k=0;
  char frs[MaxCrtFrs*sizeof(char)];

  while (fread(buf,1,1,nomarq)==1)
  {
    if (*buf!=*fim1)		// Separadores de frase
    {
      frs[k]=*buf;
      if (k<MaxCrtFrs)
        k++;
      else
      {
        frs[k]=0;
        return(0);
      }
    }
    else if (k!=0)
    {
      frs[k]=0;
      k=0;
      if (*frs==*coment)
        continue;
      strcpy(frase,frs);
      return(1);
    }
  }
  frs[k]=0;
  strcpy(frase,frs);

  return(1);
}
