#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void loadSubunits
(
  char nome_subunidades[512], // nome do arquivo com as subunidades foneticas
  int *nSubunits1, // number of phonetic subunits
  char ***subunits1
)
{
	FILE *arquivo;
	char data[512];
	int i;
	char **subunits;
	int nSubunits = 0;
	
	arquivo = fopen(nome_subunidades,"rt");
	if (arquivo == NULL)
  	{
		printf("Error opening file %s.\n",nome_subunidades);
		exit(1);
  	}
	
  	// Counting phonetic subunits
  	fgets(data, 500, arquivo);
  	while (!feof(arquivo))
  	{
  	    nSubunits += 1;
  		fgets(data, 500, arquivo);
  	}
  	
  	// Allocating memory
  	subunits = malloc(sizeof(char *)*nSubunits);
  	for (i=0;i<nSubunits;i++)
  		subunits[i] = malloc(sizeof(char)*10);

    // Reading subunits
    rewind(arquivo);
    fgets(data, 500, arquivo);
  
  	// removing control characters (\r \n) frm the end of string
  	while (data[strlen(data)-1]=='\r' || data[strlen(data)-1]=='\n')
		data[strlen(data)-1] = '\0';
	
	i=0;
  	while (!feof(arquivo))
  	{
  		strcpy(subunits[i++],data);
  		
		// Reading new entry from file
		if (!feof(arquivo))
			fgets(data, 500, arquivo);
		
  		// removing control characters (\r \n) frm the end of string
		while (data[strlen(data)-1]=='\r' || data[strlen(data)-1]=='\n')
			data[strlen(data)-1] = '\0';
  	} // end of configuration file reading
   
	fclose(arquivo);

	*subunits1 = subunits;
	*nSubunits1 = nSubunits;
}

void freeSubunits
(
	char ***subunits1,
	int nSubunits
)
{
	char **subunits;
	int i;
	
	subunits = *subunits1;
	
	for(i=0;i<nSubunits;i++)
		free(subunits[i]);
	free(subunits);
}
