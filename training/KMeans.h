#ifndef KMEANS_H_
#define KMEANS_H_

void K_Means
(
	char **subunits,        // list of phonetic subunits
	struct config trainConfig, // training configuration settings  
	struct modelo_HMM *HMM // modelos HMM das subunidades foneticas
);

#endif /*KMEANS_H_*/
