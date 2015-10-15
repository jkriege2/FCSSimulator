
SCRIPTS = $(filter-out %btsplot.plt %btsplotwith0.plt %btsplotwith1.plt %btsplotwith2.plt %btsplotwith3.plt %btsplotwith4.plt %btsplotwith5.plt,$(wildcard ${DIR}/*.plt))
SCRIPTS_BTS = $(wildcard ${DIR}/*btsplot.plt) $(wildcard ${DIR}/*btsplotwith0.plt ${DIR}/*btsplotwith1.plt ${DIR}/*btsplotwith2.plt ${DIR}/*btsplotwith3.plt ${DIR}/*btsplotwith4.plt ${DIR}/*btsplotwith5.plt)

GNUPLOT =gnuplot

SHELL = sh -c 

SCRIPTS_TARGET = $(subst .plt,.target,$(SCRIPTS))
SCRIPTSBTS_TARGET = $(subst .plt,.target,$(SCRIPTS_BTS))

TERMINAL_COMMAND=

all: show_nobts ${SCRIPTS_TARGET}

show_withbts:
	@echo -e "INCLUDING BTS PLOTS"

show_nobts:
	@echo -e "NO BTS PLOTS"


withbts: show_withbts ${SCRIPTS_TARGET} ${SCRIPTSBTS_TARGET}

THISDIR=$(shell pwd)

%.target: %.plt
	@echo -e "starting on $< ..."
	-${SHELL} 'cd $(DIR) && pwd && ${TERMINAL_COMMAND} ${GNUPLOT}${EXE_SUFFIX} $(notdir $<) < /bin/true '
	@echo -e "work on $< DONE!"
