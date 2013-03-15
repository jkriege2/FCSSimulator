SCRIPTS= emccd_testgain1.ini\
         emccd_testgain2.ini\
         emccd_testgain3.ini\
         emccd_testgain4.ini\
         emccd_testgain5.ini\
         emccd_testgain6.ini\
         emccd_testgain7.ini\
         emccd_testgain8.ini\
         emccd_testgain9.ini
		 
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
