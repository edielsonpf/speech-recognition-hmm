#include <math.h>
#include <stdlib.h>

//--------------------------------------------------------------------------+
//  Os coef mel-cepstrais sao obtidos segundo o algoritmo de Mermelstein;
//  os passos a seguir sao:
//
//  - Calculo da FFT das am's pertencentes a janela de analise.
//  - Calculo do quadrado do modulo da FFT (equivale a energia).
//  - Filtragem do modulo acima pelo banco de filtros na escala mel,
//    obtendo-se a energia por filtro.
//  - Calculo do log da energia por filtro.
//  - Calculo da transformada inversa, obtendo-se os coef's mel-cepstrais
//    (esto e' implementado calculando-se a DCT: discrete cosine transform)
//
//  OBS: Lwindow: num de ams/janela; xw: am's janeladas
// +------------------------------------------------------------------------+
void calcMel
(
  int fs,           // frequencia de amostragem do sinal
  double *mod_fft, // modulo ao quadrado da FFT (energia)
  int N,           // numero de pontos da FFT
  int ordem,       // dimensao dos vetores acusticos (12 neste caso)
  double *Mel_ceps // retorna os parametros mel-cepstrais

)
{
  // Declaracao das variaveis locais
  int *fc;            // freq chaves dos filtros mel para um banco de filtros variável
  int *f;             // freq chaves dos filtros mel para um banco de filtros fixo
  double *epb;         // energia/banda e log da energia/banda
  register int i, j, k;               // var auxiliares para loops
  int Nj;             // numero de bandas do banco de filtros
  double b;
  float indice;       // var que armazena os novos valores das frequencias (freq. adaptadas) em float


  // Definindo valores de variaveis que mudam com a frequencia de amostragem do sinal
  switch(fs)
  {
    case 8000:
    {
      Nj=18;
      break;
    }
    case 11025:
    {
      Nj=21;
      break;
    }
    case 16000:
    {
      Nj=24;
      break;
    }
  }

  // Alocando memoria
  fc = malloc(sizeof(int)*(Nj+2));
  f = malloc(sizeof(int)*26);
  epb = malloc(sizeof(double)*(N/2));

  f[0]  = 0 ;    f[1]  = 100;   f[2] = 200;     f[3] = 300;    f[4] = 400;
  f[5] = 500;    f[6]  = 600;   f[7] = 700;     f[8] = 800;    f[9] = 900;
  f[10] = 1000;  f[11] = 1149;  f[12] = 1320;   f[13] = 1516;  f[14] = 1741;
  f[15] = 2000;  f[16] = 2297;  f[17] = 2639;   f[18] = 3031;  f[19] = 3482;
  f[20] = 4000;  f[21] = 4595;  f[22] = 5278;   f[23] = 6063;  f[24] = 6964;
  f[25] = 8000;


  // Definicao das freq chaves que definem as freq centrais e as larguras de banda
  // dos filtros na escala mel.
  switch(fs)
  {
    case 8000:
     {
       for (i=0;i<Nj+2;i++)
       {
        indice=(f[i]*512)/fs;
        fc[i] = floor(indice + 0.5 );
        if (i == (Nj+1))
           fc[i] = 223;
       }
      break;
     }
    case 11025:
     {
       for (i=0;i<Nj+2;i++)
       {
        indice=(f[i]*512)/fs;
        fc[i] = floor(indice + 0.5 );
        if (i==(Nj+1))
          fc[i] = 245;
       }
      break;
     }
    case 16000:
     {
       for (i=0;i<Nj+2;i++)
       {
         indice=(f[i]*1024)/fs;
         fc[i] = floor(indice + 0.5 );
         if (i==(Nj+1))
           fc[i] = 512;
       }
      break;
     }
  }


  // Filtragem da FFT pelo banco de filtros na escala mel e calculo do log da
  // saida dos filtros
  for(k=0; k<N/2; k++)  // inicializacao
    epb[k] = 0.0;

  for(k=1; k<=Nj; k++) // k: filtro mel
  {
    for(j=1; j<=(fc[k]-fc[k-1]); j++)      // energia da semi-banda inf + freq central
      epb[fc[k]] += mod_fft[fc[k-1]+j] * (j * 1.0/(fc[k]-fc[k-1]));

	for(j=1; j<(fc[k+1]-fc[k]); j++)       // energia da semi-banda sup
	  epb[fc[k]] += mod_fft[fc[k]+j] * (1.0 - j * 1.0/(fc[k+1]-fc[k]));

	if(epb[fc[k]] > 0)
      epb[fc[k]] = log10(epb[fc[k]]);  // log10 da epb (<math.h>)
	else
      epb[fc[k]] = -100.0;
  }

  // Calculo dos coef mel-cepstrais 1 a 12 (Mel_ceps[0] = ganho ==> inutil)
  for(j=1; j<=ordem; j++) // discrete cosine transform (DCT)
  {
	b = 0.0;               // var aux para calculo dos coef com precissao dupla
	for(k=1; k<=Nj; k++)
      b += epb[fc[k]] * cos(j*(k-0.5)*M_PI/Nj);
	Mel_ceps[j-1] = (double) b;   // armazenamento com precissao simples
  }

  // Desalocando ponteiros
  free(fc);
  free(f);
  free(epb);

}       // fim da funcao "Calculo_Mel_Cepstrais()"
