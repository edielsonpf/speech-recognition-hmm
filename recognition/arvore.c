//---------------------------------------------------------------------------
#include <stdlib.h>
#include <stdio.h>

#include "estruturas.h"
#include "arvore.h"
//---------------------------------------------------------------------------

// Gera arvores que representam as palavras do vocabulario
// Estas arvores sao armazenadas no arquivo 'arvores.dat'
void GeraArvores
(
  int *MaxProfundidade1, // numero maximo de niveis observado
  long n_palavras, // numero de palavras no vocabulario
  int n_fones, // numero de fones utilizados p/ transcricao fonetica das palavras
  struct vocab vocabulario // vocabulario

)
{
  FILE *arquivo; // handle para o arquivo de saida
  int arvore; // conta as arvores
  int register i,j,k; // contadores
  int ja_tem; // verifica se ja existe (1) ou nao (0) arvore com fone sob analise
  const int max_filhos = 100; // numero maximo de filhos que um no pode ter
  int MaxProfundidade = 0; // profundidade do no mais profundo observado em todo o vocabulario
  const int max_nos = 1000; // numero maximo de nos por arvore
  int n_arvores = 0; // numero de arvores criadas no vocabulario
  int *n_nos; // numero de nos em cada arvore
  int palavra; // conta as palavras
  int qual; // verifica a qual arvore pertence a palavra sob analise
  struct no **vocab; // armazena as palavras do vocabulario em arvore

  //AnsiString mensagem; // var aux p/ apresentacao de menssagens

  // Alocando memoria
  vocab = malloc(sizeof(struct no *)*(n_fones-1));
  for (i=0;i<(n_fones-1);i++)
    vocab[i] = malloc(sizeof(struct no)*max_nos);
  for (i=0;i<(n_fones-1);i++)
    for (j=0;j<max_nos;j++)
      vocab[i][j].filhos = malloc(sizeof(int)*max_filhos);

  n_nos = malloc(sizeof(int)*(n_fones-1));

  // Inicializando variavel vocab
  for (i=0;i<(n_fones-1);i++)
    for (j=0;j<max_nos;j++)
    {
      vocab[i][j].n_filhos = 0;
      vocab[i][j].palavra = -1;
    }

  // Inicializando contagens dos nos nas arvores
  for (i=0;i<(n_fones-1);i++)
    n_nos[i] = 0;

  for (palavra=0;palavra<n_palavras;palavra++)
  {
    // Verifica se ja existe arvore com o fone inicial da palavra atual
    ja_tem = 0;
    arvore = 0;

    while ((arvore < n_arvores) && (ja_tem == 0))
    {
      if (vocabulario.Mw[palavra][1] == vocab[arvore][0].fone)
      {
        ja_tem = 1;
        qual = arvore;
      }
      arvore++;
    }

    if (ja_tem == 0)

      // Cria arvore se nao existir uma
      NovaArvore(&n_arvores,&n_nos,&vocab,vocabulario.Mw[palavra],palavra,3);
    else

      // Atualiza arvore existente
      AtualizaArvore(&n_nos[qual],&n_arvores,&vocab[qual],qual,vocabulario.Mw[palavra],palavra,3);

      // Atualizando max_niveis
      if (vocabulario.Mw[palavra][0] > MaxProfundidade)
        MaxProfundidade = vocabulario.Mw[palavra][0];
  }

  // O numero de niveis e na verdade (max_niveis+1) pois o primeiro nivel e o nivel 0
  MaxProfundidade++;

  // Salvando arvores geradas em arquivo
  arquivo = fopen("arvores.dat","wb");
  if (arquivo == NULL)
  {
    printf("Error creating file arvores.dat. Verify your permissions!\n");
    exit(1);
  }

  // Numero de arvores geradas
  fwrite(&n_arvores,sizeof(n_arvores),1,arquivo);

  // Dados das arvores
  for (i=0;i<n_arvores;i++)
  {
    // Numero de nos
    fwrite(&n_nos[i],sizeof(n_nos[i]),1,arquivo);

    // Dados de cada no
    for (j=0;j<n_nos[i];j++)
    {
      // Nivel do no
      fwrite(&vocab[i][j].profundidade,sizeof(vocab[i][j].profundidade),1,arquivo);
      // Fone
      fwrite(&vocab[i][j].fone,sizeof(vocab[i][j].fone),1,arquivo);
      // Palavra
      fwrite(&vocab[i][j].palavra,sizeof(vocab[i][j].palavra),1,arquivo);
      // No pai
      fwrite(&vocab[i][j].pai,sizeof(vocab[i][j].pai),1,arquivo);
      // Numero de filhos
      fwrite(&vocab[i][j].n_filhos,sizeof(vocab[i][j].n_filhos),1,arquivo);
      // Filhos
      for (k=0;k<vocab[i][j].n_filhos;k++)
        fwrite(&vocab[i][j].filhos[k],sizeof(vocab[i][j].filhos[k]),1,arquivo);
    }
  }
  // Fechando arquivo
  fclose(arquivo);

  // Desalocando memoria
  free(n_nos);

  for (i=0;i<(n_fones-1);i++)
    for (j=0;j<max_nos;j++)
      free(vocab[i][j].filhos);
  for (i=0;i<(n_fones-1);i++)
    free(vocab[i]);
  free(vocab);

  // Valor de retorno da funcao
  *MaxProfundidade1 = MaxProfundidade;
}

//------------------------------------------------------------------------------
// Gera arvore com fone inicial diferente
void NovaArvore
(
  int *n_arvores1, // numero de arvores geradas
  int **n_nos1, // numero de nos de cada arvore
  struct no ***vocab1, // arvores geradas
  int *MW, // modelos foneticos de cada palavra
  int palavra, // palavra sob analise
  int n_estados // numero de estados dos modelos HMM dos fones
)
{
  int i; // contador
  int n_arvores; // numero de arvores
  int *n_nos; // numero de nos em cada arvore
  struct no **vocab;

  n_nos = *n_nos1;
  n_arvores = *n_arvores1;
  vocab = *vocab1;

  // Criando no raiz
  vocab[n_arvores][0].profundidade = 0;
  vocab[n_arvores][0].fone = MW[1];
  vocab[n_arvores][0].pai = -1;

  // Alocando memoria p/ estruturas like e elapse do no
  vocab[n_arvores][0].like = malloc(sizeof(double)*n_estados);
  vocab[n_arvores][0].elapse = malloc(sizeof(int)*n_estados);

  if (MW[0] != 1) // a palavra tem mais de 1 fone
  {
    vocab[n_arvores][0].filhos[0] = 1;
    vocab[n_arvores][0].n_filhos++;
  }
  else // a palavra tem apenas 1 fone
    vocab[n_arvores][0].palavra = palavra;

  // Atualizando contagens de nos na arvore
  n_nos[n_arvores]++;

  // Criando nos filhos p/ o caso de a palavra ter mais de 1 fone
  for (i=2;i<=MW[0];i++)
  {
    vocab[n_arvores][i-1].profundidade = i-1;
    vocab[n_arvores][i-1].fone = MW[i];
    vocab[n_arvores][i-1].pai = n_nos[n_arvores]-1;

    // Alocando memoria p/ estruturas like e elapse do no
    vocab[n_arvores][i-1].like = malloc(sizeof(double)*n_estados);
    vocab[n_arvores][i-1].elapse = malloc(sizeof(int)*n_estados);

    if (i == MW[0]) // ultimo fone da palavra
      vocab[n_arvores][i-1].palavra = palavra;
    else // a palavra tem mais fones
    {
      vocab[n_arvores][i-1].filhos[0] = i;
      vocab[n_arvores][i-1].n_filhos++;
    }

    // Atualizando contagens de nos na arvore
    n_nos[n_arvores]++;
  }

  n_arvores++;

  *n_nos1 = n_nos;
  *n_arvores1 = n_arvores;
  *vocab1 = vocab;
}

//------------------------------------------------------------------------------
// Adiciona palavra a arvore ja existente
void AtualizaArvore
(
  int *n_nos1, // numero de nos desta arvore
  int *n_arvores1, // numero de arvores
  struct no **vocab1, // arvores geradas
  int qual, // arvore a ser atualizada
  int *MW, // modelos foneticos de cada palavra
  int palavra, // palavra sob analise
  int n_estados // numero de estados dos modelos HMM dos fones
)
{
  int i,j; // contadores
  int ident; // var aux p/ identificar os nos
  int ja_tem; // verifica se um determinado no ja consta da arvore
  int n_nos; // numero de nos na arvore
  int pai; // pai do no atual
  struct no *vocab; // armazena a arvore correspondente ao fone selecionado

  n_nos = *n_nos1;
  vocab = *vocab1;

  if (MW[0] == 1) // a palavra tem apenas 1 fone
    vocab[0].palavra = palavra;
  else // a palavra tem mais de 1 fone
  {
    pai = 0;
    for (i=2;i<=MW[0];i++)
    {
      // Verificando se o fone ja consta da arvore
      ja_tem = 0;
      j = 0; // identifica o no filho
      while ((j<vocab[pai].n_filhos) && (ja_tem == 0))
      {
        if (vocab[vocab[pai].filhos[j]].fone == MW[i])
        {
          ja_tem = 1;
          ident = vocab[pai].filhos[j]; // identificando no reconhecido
        }
        j++;
      }

      if (ja_tem == 1) // o fone ja consta da arvore
      {
        pai = ident; // o no pai passa a ser o no identificado acima
        if (i == MW[0]) // testa se eh fone final da palavra
          vocab[ident].palavra = palavra;
      }
      else // adiciona novo no se nao houver
      {
        // Atualizando contagem dos nos da arvore
        n_nos++;

        // Atualizando no pai
        vocab[pai].n_filhos++;
        vocab[pai].filhos[vocab[pai].n_filhos-1] = n_nos-1;

        // Criando no filho
        vocab[n_nos-1].profundidade = i-1;
        vocab[n_nos-1].fone = MW[i];
        vocab[n_nos-1].pai = pai;

        // Alocando memoria p/ estruturas like e elapse do no
        vocab[n_nos-1].like = malloc(sizeof(double)*n_estados);
        vocab[n_nos-1].elapse = malloc(sizeof(int)*n_estados);

        if (i == MW[0]) // se eh fone final da palavra, entao este no corresponde a uma palavra
          vocab[n_nos-1].palavra = palavra;

        // Atualizando o no pai
        pai = n_nos - 1;
      }
    }
  }

  *n_nos1 = n_nos;
  *vocab1 = vocab;
}

// Carrega arvores geradas pela funcao GeraArvores
void CarregaArvores
(
  int n_niveis, // numero maximo de palavras na locucao
  int *n_arvores1, // numero de arvores criadas no vocabulario
  int n_estados, // numero de estados dos modelos HMM dos fones
  int **n_nos1, // numero de nos em cada arvore
  struct no ****vocab1 // armazena as palavras do vocabulario em arvore
)
{
  FILE* arquivo; // handle para leitura do arquivo com as arvores geradas
  int register arvore; // conta as arvores
  int register filho; // conta os nos filhos
  int register NivelG; // conta os niveis gramaticos
  int register noh; // conta os nos
  int n_arvores; // numero de arvores criadas no vocabulario
  int *n_nos; // numero de nos em cada arvore
  struct no ***vocab; // armazena as palavras do vocabulario em arvore
  //char nome[256]; // nome do arquivo onde foram armazenados os dados das arvores
  int num_nos; // numero de nos da primeira arvore (silencio)

  n_arvores = *n_arvores1;
  n_nos = *n_nos1;
  vocab = *vocab1;

  // Alocando memoria
  vocab = malloc(sizeof(struct no**)*n_niveis);

  //----------------------------------------------------------------------------
  // Para o primeiro nivel gramatico so e carregada a primeira arvore (silencio)
  NivelG = 0;

  // Carregando arvores do arquivo
  arquivo = fopen("arvores.dat","rb");
  if (arquivo == NULL)
  {
    puts("Error opening file arvores.dat. Verify whether if it exists.\n");
    exit(1);
  }

  // Numero de arvores geradas
  fread(&n_arvores,sizeof(n_arvores),1,arquivo);
  n_arvores = 1;

  // Alocando memoria
  vocab[NivelG] = malloc(sizeof(struct no*)*n_arvores);

  // Dados da arvore
  for (arvore=0;arvore<n_arvores;arvore++)
  {
    // Numero de nos
    fread(&num_nos,sizeof(num_nos),1,arquivo);

    // Alocando memoria
    vocab[NivelG][arvore] = malloc(sizeof(struct no)*num_nos);

    // Dados de cada no
    for (noh=0;noh<num_nos;noh++)
    {
      // Nivel do no
      fread(&vocab[NivelG][arvore][noh].profundidade,sizeof(vocab[NivelG][arvore][noh].profundidade),1,arquivo);
      // Fone
      fread(&vocab[NivelG][arvore][noh].fone,sizeof(vocab[NivelG][arvore][noh].fone),1,arquivo);
      // Palavra
      fread(&vocab[NivelG][arvore][noh].palavra,sizeof(vocab[NivelG][arvore][noh].palavra),1,arquivo);
      // No pai
      fread(&vocab[NivelG][arvore][noh].pai,sizeof(vocab[NivelG][arvore][noh].pai),1,arquivo);
      // Numero de filhos
      fread(&vocab[NivelG][arvore][noh].n_filhos,sizeof(vocab[NivelG][arvore][noh].n_filhos),1,arquivo);

      // Alocando memoria
      vocab[NivelG][arvore][noh].filhos = malloc(sizeof(int)*vocab[NivelG][arvore][noh].n_filhos);
      vocab[NivelG][arvore][noh].like = malloc(sizeof(double)*n_estados);
      vocab[NivelG][arvore][noh].elapse = malloc(sizeof(int)*n_estados);

      // Filhos
      for (filho=0;filho<vocab[NivelG][arvore][noh].n_filhos;filho++)
        fread(&vocab[NivelG][arvore][noh].filhos[filho],sizeof(vocab[NivelG][arvore][noh].filhos[filho]),1,arquivo);
    } // end for noh
  } // end for arvore

  // Fechando arquivo
  fclose(arquivo);

  //----------------------------------------------------------------------------

  // Carregando dados para os demais niveis
  for (NivelG=1;NivelG<n_niveis;NivelG++)
  {
    // Carregando arvores do arquivo
    arquivo = fopen("arvores.dat","rb");
    if (arquivo == NULL)
    {
      puts("Error opening file arvores.dat. Verify your permissions or if it exists.\n");
      exit(1);
    }

    // Numero de arvores geradas
    fread(&n_arvores,sizeof(n_arvores),1,arquivo);

    // Alocando memoria
    vocab[NivelG] = malloc(sizeof(struct no*)*n_arvores);
    if (NivelG == 1)
      n_nos = malloc(sizeof(int)*n_arvores);

    // Dados das arvores
    for (arvore=0;arvore<n_arvores;arvore++)
    {
      // Numero de nos
      fread(&n_nos[arvore],sizeof(n_nos[arvore]),1,arquivo);

      // Alocando memoria
      vocab[NivelG][arvore] = malloc(sizeof(struct no)*n_nos[arvore]);

      // Dados de cada no
      for (noh=0;noh<n_nos[arvore];noh++)
      {
        // Nivel do no
        fread(&vocab[NivelG][arvore][noh].profundidade,sizeof(vocab[NivelG][arvore][noh].profundidade),1,arquivo);
        // Fone
        fread(&vocab[NivelG][arvore][noh].fone,sizeof(vocab[NivelG][arvore][noh].fone),1,arquivo);
        // Palavra
        fread(&vocab[NivelG][arvore][noh].palavra,sizeof(vocab[NivelG][arvore][noh].palavra),1,arquivo);
        // No pai
        fread(&vocab[NivelG][arvore][noh].pai,sizeof(vocab[NivelG][arvore][noh].pai),1,arquivo);
        // Numero de filhos
        fread(&vocab[NivelG][arvore][noh].n_filhos,sizeof(vocab[NivelG][arvore][noh].n_filhos),1,arquivo);

        // Alocando memoria
        vocab[NivelG][arvore][noh].filhos = malloc(sizeof(int)*vocab[NivelG][arvore][noh].n_filhos);
        vocab[NivelG][arvore][noh].like = malloc(sizeof(double)*n_estados);
        vocab[NivelG][arvore][noh].elapse = malloc(sizeof(int)*n_estados);

        // Filhos
        for (filho=0;filho<vocab[NivelG][arvore][noh].n_filhos;filho++)
          fread(&vocab[NivelG][arvore][noh].filhos[filho],sizeof(vocab[NivelG][arvore][noh].filhos[filho]),1,arquivo);
      } // end for noh
    } // end for arvore

    // Fechando arquivo
    fclose(arquivo);
  } // end for nivel

  // Valores de retorno da funcao
  *n_arvores1 = n_arvores;
  *n_nos1 = n_nos;
  *vocab1 = vocab;
}


// Apresenta as arvores geradas na tela
void ApresentaArvores
(
  int n_arvores, // numero de arvores geradas
  int *n_nos, // numero de nos em cada arvore
  struct no **vocab, // armazena as palavras do vocabulario em arvore
  struct vocab vocabulario // vocabulario
)
{
  int register i,j,k; // contadores

  // Apresentando resultados na tela
  for (i=0;i<n_arvores;i++)
  {
    printf("ÁRVORE %d\n",i);
    for (j=0;j<n_nos[i];j++)
    {
      printf("Nó %d, Nível %d, Fone %d, Palavra ",j,vocab[i][j].profundidade,vocab[i][j].fone);
      if (vocab[i][j].palavra != -1)
        printf("%s",vocabulario.Mws[vocab[i][j].palavra]);
      else 
        printf(" ");
      printf(", Pai %d, # Filhos %d, Filhos: ",vocab[i][j].pai,vocab[i][j].n_filhos);
      for (k=0;k<vocab[i][j].n_filhos;k++)
        printf("%d ",vocab[i][j].filhos[k]);
       puts("\n"); 
    }
    puts("\n");
  }
}


// Desaloca arvores alocadas na funcao 'CarregaArvores
void DesalocaArvores
(
  int n_niveis, // numero maximo de palavras na locucao
  int n_arvores, // numero de arvores geradas
  int **n_nos1, // numero de nos em cada arvore
  struct no ****vocab1 // armazena as palavras do vocabulario em arvore
)
{
  int register i,j,k; // contadores
  int* n_nos; // numero de nos em cada arvore
  int n_trees; // var aux p/ desalocar ponteiros
  struct no*** vocab; // armazena as palavras do vocabulario em arvore

  n_nos = *n_nos1;
  vocab = *vocab1;

  // Desalocando memoria
  for (k=0;k<n_niveis;k++)
  {
    if (k == 0)
      n_trees = 1;
    else
      n_trees = n_arvores;
    for (i=0;i<n_trees;i++)
      for (j=0;j<n_nos[i];j++)
      {
        free(vocab[k][i][j].filhos);
        free(vocab[k][i][j].like);
        free(vocab[k][i][j].elapse);
      }
  }
  for (k=0;k<n_niveis;k++)
  {
    if (k == 0)
      n_trees = 1;
    else
      n_trees = n_arvores;
    for (i=0;i<n_trees;i++)
      free(vocab[k][i]);
  }
  for (k=0;k<n_niveis;k++)
    free(vocab[k]);
  free(vocab);


  free(n_nos);

  *n_nos1 = n_nos;
  *vocab1 = vocab;
}

// Gera listagem dos nos pertencentes a uma dada profundidade dentro da arvore
void GeraListas
(
  int MaxProfundidade, // profundidade maxima observada em todas as arvores
  int n_arvores, // numero de arvores geradas
  int *n_nos, // numero de nos em cada arvore
  struct no **vocab // armazena as palavras do vocabulario em arvore
)
{
  FILE *arquivo; // handle para manipulacao de arquivos
  int arvore; // conta as arvores
  int register i,j,k; // contadores
  int ***lista; // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int max_nos = 100; // numero maximo de nos permitido a cada nivel
  int profundidade; // profundidade de cada no dentro da arvore
  int *n_niveis; // indica o numero de niveis em cada arvore
  int noh; // conta os nos
  int **quantos; // indica quantos elementos existem em cada nivel de cada uma das listas

  // Alocando memoria e inicializando estrutura que ira armazenar a listagem
  // dos nos
  lista = malloc(sizeof(int **)*n_arvores);
  for (i=0;i<n_arvores;i++)
    lista[i] = malloc(sizeof(int *)*MaxProfundidade);
  for (i=0;i<n_arvores;i++)
    for (j=0;j<MaxProfundidade;j++)
      lista[i][j] = malloc(sizeof(int)*max_nos);

  quantos = malloc(sizeof(int *)*n_arvores);
  for (i=0;i<n_arvores;i++)
    quantos[i] = malloc(sizeof(int)*MaxProfundidade);

  n_niveis = malloc(sizeof(int)*n_arvores);
  

  for (i=0;i<n_arvores;i++)
    for (j=0;j<MaxProfundidade;j++)
      for (k=0;k<max_nos;k++)
        lista[i][j][k] = -1;

  for (i=0;i<n_arvores;i++)
    for (j=0;j<MaxProfundidade;j++)
      quantos[i][j] = 0;

  for (i=0;i<n_arvores;i++)
    n_niveis[i] = 0;

  for (arvore=0;arvore<n_arvores;arvore++)
    for (noh=0;noh<n_nos[arvore];noh++)
    {
      // Armazenando no na lista
      lista[arvore][vocab[arvore][noh].profundidade][quantos[arvore][vocab[arvore][noh].profundidade]] = noh;

      // Atualizando numero de nos no nivel
      quantos[arvore][vocab[arvore][noh].profundidade]++;

      // Atualizando numero maximo de niveis para esta arvore
      if ((vocab[arvore][noh].profundidade+1) > n_niveis[arvore])
        n_niveis[arvore] =  vocab[arvore][noh].profundidade+1;
    }

  // Salvando listas geradas em arquivo
  arquivo = fopen("listas.dat","wb");
  if (arquivo == NULL)
  {
    puts("Error creating file listas.dat. Verify permissions!\n");
    exit(1);
  }

  for (arvore=0;arvore<n_arvores;arvore++)
  {
    // Numero de niveis da arvore
    fwrite(&n_niveis[arvore],sizeof(n_niveis[arvore]),1,arquivo);

    for (profundidade=0;profundidade<n_niveis[arvore];profundidade++)
    {
      fwrite(&quantos[arvore][profundidade],sizeof(quantos[arvore][profundidade]),1,arquivo);

      for (noh=0;noh<quantos[arvore][profundidade];noh++)
        fwrite(&lista[arvore][profundidade][noh],sizeof(lista[arvore][profundidade][noh]),1,arquivo);
    }
  }
  fclose(arquivo);

  // Desalocando ponteiros
  for (i=0;i<n_arvores;i++)
    for (j=0;j<MaxProfundidade;j++)
      free(lista[i][j]);
  for (i=0;i<n_arvores;i++)
    free(lista[i]);
  free(lista);

  for (i=0;i<n_arvores;i++)
    free(quantos[i]);
  free(quantos);

  free(n_niveis);
}

// Carrega as listas geradas pela funcao 'GeraListas'
// Aloca os ponteiros 'lista', 'profundidades' e 'quantos'
void CarregaListas
(
  int ****lista1, // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int n_arvores, // numero de arvores geradas
  int **profundidades1, // profundidade maxima observada em cada uma das arvores
  int *n_nos, // numero de nos em cada arvore
  int ***quantos1 // indica quantos elementos existem em cada nivel de cada uma das listas
)
{
  FILE *arquivo; // handle para manipulacao de arquivos
  int register arvore; // conta as arvores
  int ***lista; // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int register profundidade; // profundidade do no
  int *profundidades; // indica o numero de niveis em cada arvore
  int register noh; // conta os nos
  int **quantos; // indica quantos elementos existem em cada nivel de cada uma das listas

  arquivo = fopen("listas.dat","rb");
  if (arquivo == NULL)
  {
    puts("Error opening file listas.dat!\n");
    exit(1);
  }
  profundidades = malloc(sizeof(int)*n_arvores);
  quantos = malloc(sizeof(int *)*n_arvores);
  lista = malloc(sizeof(int **)*n_arvores);

  for (arvore=0;arvore<n_arvores;arvore++)
  {
    // Profundidade do no mais profundo da arvore
    fread(&profundidades[arvore],sizeof(profundidades[arvore]),1,arquivo);

    // Alocando memoria
    quantos[arvore] = malloc(sizeof(int)*profundidades[arvore]);
    lista[arvore] = malloc(sizeof(int *)*profundidades[arvore]);

    for (profundidade=0;profundidade<profundidades[arvore];profundidade++)
    {
      // Numero de nos em cada nivel de profundidade
      fread(&quantos[arvore][profundidade],sizeof(quantos[arvore][profundidade]),1,arquivo);

      // Alocando memoria
      lista[arvore][profundidade] = malloc(sizeof(int)*quantos[arvore][profundidade]);

      // Listagem dos nos de cada nivel
      for (noh=0;noh<quantos[arvore][profundidade];noh++)
        fread(&lista[arvore][profundidade][noh],sizeof(lista[arvore][profundidade][noh]),1,arquivo);
    }
  }
  fclose(arquivo);

  // Atualizando dados de saida
  *lista1 = lista;
  *profundidades1 = profundidades;
  *quantos1 = quantos;
}

/*
// Funcao que que apresenta na tela as listas geradas pela funcao 'GeraListas'
void ApresentaListas
(
  int ***lista, // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int n_arvores, // numero de arvores geradas
  int *profundidades, // profundidade maxima observada em cada uma das arvores
  int **quantos, // indica quantos elementos existem em cada nivel de cada uma das listas
  struct no **vocab // armazena as palavras do vocabulario em arvore
)
{
  int register arvore; // conta as arvores
  //AnsiString mensagem; // var aux p/ apresentar dados na tela
  int register profundidade; // profundidade do no dentro da arvore
  int register noh; // conta os nos

  for (arvore=0;arvore<n_arvores;arvore++)
  {
    mensagem = "Arvore " + IntToStr(arvore);
    FormMain->RichEdit1->Lines->Add(mensagem);
    for (profundidade=0;profundidade<profundidades[arvore];profundidade++)
    {
      mensagem = "Profundidade " + IntToStr(profundidade);
      mensagem += " " + IntToStr(quantos[arvore][profundidade]);
      mensagem += " nos";
      FormMain->RichEdit1->Lines->Add(mensagem);

      mensagem = "";
      for (noh=0;noh<quantos[arvore][profundidade];noh++)
        mensagem += IntToStr(vocab[arvore][lista[arvore][profundidade][noh]].fone) + " ";
      FormMain->RichEdit1->Lines->Add(mensagem);
      MessageBox(NULL," "," ",MB_OK);
    }
  }
}
*/
// Desaloca ponteiros 'lista', 'profundidades' e 'quantos', alocados na funcao 'CarregaListas'
void DesalocaListas
(
  int ****lista1, // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int n_arvores, // numero de arvores geradas
  int **profundidades1, // profundidade maxima observada em cada uma das arvores
  int ***quantos1 // indica quantos elementos existem em cada nivel de cada uma das listas
)
{
  int register arvore; // conta as arvores
  int ***lista; // estrutura que armazena os nos pertencentes a cada nivel, para cada arvore
  int register profundidade; // profundidade do no na arvore
  int *profundidades; // profundidade do no mais profundo em cada arvore
  int **quantos; // indica quantos elementos existem em cada nivel de cada uma das listas

  lista = *lista1;
  profundidades = *profundidades1;
  quantos = *quantos1;

  for (arvore=0;arvore<n_arvores;arvore++)
    for (profundidade=0;profundidade<profundidades[arvore];profundidade++)
      free(lista[arvore][profundidade]);
  for (arvore=0;arvore<n_arvores;arvore++)
    free(lista[arvore]);
  free(lista);

  for (arvore=0;arvore<n_arvores;arvore++)
    free(quantos[arvore]);
  free(quantos);

  free(profundidades);

  *lista1 = NULL;
  *profundidades1 = NULL;
  *quantos1 = NULL;

}
