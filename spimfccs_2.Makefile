SCRIPTS= spimfccs2_3.ini\
         spimfccs2_4.ini\
         spimfccs2_2.ini\
         spimfccs2_5.ini\
         spimfccs2_1.ini\
         spimfccs2_0.ini
		 
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
