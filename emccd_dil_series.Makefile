SCRIPTS= emccd_dil_series14.ini\
         emccd_dil_series13.ini \
		 emccd_dil_series12.ini \
		 emccd_dil_series11.ini\
         emccd_dil_series5.ini\
         emccd_dil_series1.ini\
         emccd_dil_series7.ini\
         emccd_dil_series2.ini\
         emccd_dil_series9.ini\
         emccd_dil_series10.ini\
         emccd_dil_series4.ini\
         emccd_dil_series3.ini\
         emccd_dil_series8.ini\
         emccd_dil_series6.ini
		 
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
