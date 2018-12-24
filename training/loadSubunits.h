#ifndef LOADSUBUNITS_H_
#define LOADSUBUNITS_H_

void loadSubunits
(
  char nome_subunidades[512], // nome do arquivo com as subunidades foneticas
  int *nSubunits, // number of phonetic subunits
  char ***subunits1
);

void freeSubunits
(
	char ***subunits1,
	int nSubunits
);

#endif /*LOADSUBUNITS_H_*/
