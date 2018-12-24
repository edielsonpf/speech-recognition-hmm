//---------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>

//---------------------------------------------------------------------------
// Estrutura que armazena o cabecalho de arquivo .wav
struct soundhdr {
char riff[4]; /* "RIFF" */
long flength; /* file length in bytes */
char wave[4]; /* "WAVE" */
char fmt[4]; /* "fmt " */
long block_size; /* in bytes (generally 16) */
short format_tag; /* 1=PCM, 257=Mu-Law, 258=A-Law, 259=ADPCM */
short num_chans; /* 1=mono, 2=stereo */
long srate; /* Sampling rate in samples per second */
long bytes_per_sec; /* bytes per second */
short bytes_per_samp; /* 2=16-bit mono, 4=16-bit stereo */
short bits_per_samp; /* Number of bits per sample */
char data[4]; /* "data" */
long dlength; /* data length in bytes (filelength - 44) */
};

//---------------------------------------------------------------------------
// Funcao que abre um arquivo wav de 16 bits e devolve um ponteiro contendo
// os dados, bem como o numero de amostras
void leWav
(
  short *BitsPorAmostra1, // numero de bits por amostra
  int *fs1, // frequencia de amostragem
  char *nome_arquivo, // nome do arquivo a ser lido
  int *nSamples, // number of samples
  short **x1 // samples of the audio file
)
{
  FILE *arquivo; // handle para leitura do arquivo wav
  short BitsPorAmostra; // numero de bits por amostra
  struct soundhdr cabecalho; // cabecalho do arquivo wav
  int fs; // frequencia de amostragem
  int tamanho; // numero de amostras do sinal
  short *x; // buffer onde sera armazenado o sinal
  int nLidos;
  int tam;
	
  // Abrindo arquivo
	arquivo = fopen(nome_arquivo,"rb");
	
	if (arquivo == NULL)
	{
		printf("Error opening file %s.\n",nome_arquivo);
		exit(1);
	}

  // Lendo tamanho do arquivo (sem o cabecalho) e alocando vetor de dados
  fseek (arquivo,0,SEEK_END);
  tamanho = ftell(arquivo);
  tamanho -= 44;
  tamanho /= sizeof(short);
    
  x = malloc(sizeof(short)*tamanho);

  // Posicionando ponteiro no inicio o cabecalho
  rewind(arquivo);

  // Lendo cabecalho
	nLidos = fread(&cabecalho,sizeof(struct soundhdr),1,arquivo);
	tam = sizeof(struct soundhdr);
	
  // Informacoes sobre frequencia de amostragem e no. de bits por amostra
  BitsPorAmostra = cabecalho.bytes_per_samp;
  fs = cabecalho.srate;
	
  // Posicionando ponteiro no inicio dos dados
	fseek(arquivo,44,SEEK_SET);

  // Lendo dados
	nLidos = fread(x,sizeof(short),tamanho,arquivo);
  
  // Fechando handle do arquivo
  fclose(arquivo);

  // Valores de retorno da funcao
  *BitsPorAmostra1 = BitsPorAmostra;
  *fs1 = fs;
  *nSamples = tamanho;
  *x1 = x;

}
