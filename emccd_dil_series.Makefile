SCRIPTS= emccd_dil_series1.ini\
         emccd_dil_series2.ini
		 
SHELL = sh

SCRIPTS_TARGET = $(subst .ini,.target,$(SCRIPTS))

		

ifeq ($(findstring Msys,$(OS)),Msys)
EXE_SUFFIX=.exe
TERMINAL_COMMAND=
else
EXE_SUFFIX=
TERMINAL_COMMAND=konsole -e 
TERMINAL_COMMAND=x-terminal-emulator -e 
endif		
		 
all: ${SCRIPTS_TARGET}

%.target: %.ini
	@echo -e "starting on $< ..."
	${TERMINAL_COMMAND} ./diffusion4${EXE_SUFFIX} $<
	@echo -e "work on $< DONE!"
