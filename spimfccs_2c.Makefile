SCRIPTS= spimfccs_2c_0.ini\
         spimfccs_2c_1.ini\
         spimfccs_2c_2.ini\
         spimfccs_2c_3.ini\
         spimfccs_2c_4.ini\
         spimfccs_2c_5.ini\
         spimfccs_2c_6.ini\
         spimfccs_2c_7.ini\
         spimfccs_2c_8.ini\
         spimfccs_2c_9.ini\
		 
SHELL = sh

SCRIPTS_TARGET = $(subst .ini,.target,$(SCRIPTS))

TERMINAL_COMMAND=		

ifeq ($(findstring Msys,$(OS)),Msys)
EXE_SUFFIX=.exe
TERMINAL_COMMAND=
else
EXE_SUFFIX=
#TERMINAL_COMMAND=konsole -e 
#TERMINAL_COMMAND=x-terminal-emulator -e 
endif		
		 
all: ${SCRIPTS_TARGET}

%.target: %.ini
	@echo -e "starting on $< ..."
	${TERMINAL_COMMAND} ./diffusion4${EXE_SUFFIX} $< > $<.log
	@echo -e "work on $< DONE!"
