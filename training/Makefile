#Makefile projeto Reconhecimento de Fala
#Edielson Prevato Frigieri
#Atualização: 30/07/07

CC=gcc 
#CC=/media/sda4/projetos/arm/atmel/buildroot-v23434/build_arm/staging_dir/usr/bin/arm-linux-gcc-4.2.4 
CFLAGS= -c -Wall
#LDFLAGS= -lm -ldl -lasound
LDFLAGS=-lm

SOURCES=main.c energia.c fft.c forback.c gausFunc.c HMMFiles.c initVit.c KMeans.c lbg.c loadConfig.c loadSubunits.c loadTranscription.c lpc.c mel.c mystrfun.c par.c pca.c perfil.c preProc.c wav.c

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=treina

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *o
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(EXECUTABLE)
         
