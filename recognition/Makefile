#Makefile projeto Reconhecimento de Fala
#Edielson Prevato Frigieri
#Atualização: 30/07/07

CC=gcc 
#CC=/media/sda4/projetos/arm/atmel/buildroot-v23434/build_arm/staging_dir/usr/bin/arm-linux-gcc-4.2.4 
CFLAGS= -c -Wall
#LDFLAGS= -lm -ldl -lasound
LDFLAGS=-lm

SOURCES=main.c arvore.c energia.c fft.c gausFunc.c grammar.c herrmanNey.c loadConfig.c loadGrammar.c loadHMM.c loadVocab.c lpc.c mel.c mystrfun.c oneStep.c par.c pca.c perfil.c preProc.c reconheceFrases.c wav.c wordDurationModel.c

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=reco

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *o
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(EXECUTABLE)
         
