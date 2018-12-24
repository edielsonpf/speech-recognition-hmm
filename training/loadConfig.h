#ifndef LOADCONFIG_H_
#define LOADCONFIG_H_

// Function that load the training configuration
int loadConfig
(
  	char *configFile,
  	struct config *trainConfig,
  	struct modelo_HMM *HMM
);

// Function that shows the training configuration in the screen (just for debug)
void showConfig
(
  	struct config trainConfig,
  	struct modelo_HMM HMM
);

#endif /*LOADCONFIG_H_*/
