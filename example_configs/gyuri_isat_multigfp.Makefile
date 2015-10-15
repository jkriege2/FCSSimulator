SCRIPTS= gyuri_isat_multigfp0.ini\
	     gyuri_isat_multigfp4.ini\
         gyuri_isat_multigfp8.ini\
         gyuri_isat_multigfp1.ini\
         gyuri_isat_multigfp5.ini\
         gyuri_isat_multigfp7.ini\
         gyuri_isat_multigfp2.ini\
         gyuri_isat_multigfp3.ini\
         gyuri_isat_multigfp6.ini\
         gyuri_isat_multigfp9.ini\
         gyuri_isat_multigfp10.ini\
         gyuri_isat_multigfp11.ini\
         gyuri_isat_multigfp12.ini\
         gyuri_isat_multigfp13.ini\
         gyuri_isat_multigfp14.ini\
         gyuri_isat_multigfp15.ini\
		 gyuri_isat_multigfp16.ini\
         gyuri_isat_multigfp17.ini\
         gyuri_isat_multigfp18.ini\
         gyuri_isat_multigfp19.ini


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
	@echo -e "random wait before starting $< ..."
	@bash random_sleep.sh
	@bash -c "sleep $$[ ( $$RANDOM % 10 )  + 1 ]s"
	@echo -e "starting on $< ..."
	${TERMINAL_COMMAND} ./diffusion4${EXE_SUFFIX} $< > $<.log
	@echo -e "work on $< DONE!"
