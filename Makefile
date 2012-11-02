all: Release

LIB=../../../LIB/trunk
CC=g++
CFLAGS =  -Wall #--enable-auto-import
LDFLAGS = -lgsl -lgslcblas -lm

Release: CFLAGS += -O2 -march=native -ffast-math

Debug: CC += -DDEBUG -g

EXECUTABLE=diffusion4
SRC_FILE= browniandynamics.cpp \
          diffusiontools.cpp \
          dynamicsfromfiles2.cpp \
          fcsmeasurement.cpp \
          fluorescenceimaging.cpp \
          fluorescencemeasurement.cpp \
          fluorophordynamics.cpp \
          main.cpp \
          gridrandomwalkdynamics.cpp \
          $(LIB)/datatable.cpp \
          $(LIB)/highrestimer.cpp \
          $(LIB)/jkiniparser2.cpp \
          $(LIB)/jkmathparser.cpp \
          $(LIB)/statistics_tools.cpp \
          $(LIB)/tools.cpp \


SRC_FILE_O = $(subst .cpp,.o,$(SRC_FILE))


ifeq ($(findstring Msys,$(OS)),Msys)
PREFIX=/mingw
EXE_SUFFIX=.exe
SO_SUFFIX=.dll
SO_PATH=$(PREFIX)/bin
else
PREFIX=/usr/local
EXE_SUFFIX=
SO_SUFFIX=.so
SO_PATH=$(PREFIX)/lib
endif


Debug: ${EXECUTABLE}$(EXE_SUFFIX)

Release: ${EXECUTABLE}$(EXE_SUFFIX)

${EXECUTABLE}: ${SRC_FILE_O}
	$(CC) $(CFLAGS) -o $(EXECUTABLE)$(EXE_SUFFIX) ${SRC_FILE_O} $(LDFLAGS)

$(SRC_FILE_O): %.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f *.exe
	rm -f ${EXECUTABLE}$(EXE_SUFFIX)

	rm -f *.o



