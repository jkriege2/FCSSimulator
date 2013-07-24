all: Release

LIB=../../../LIB/trunk
CC=g++
CFLAGS =  -Wall -DHAVE_INLINE -DNO_LIBTIFF -DJKIMAGE_USES_TINYTIFF -m64 -I../../../LIB/trunk/ #--enable-auto-import
LDFLAGS = -lgsl -lgslcblas -lm

Release: CFLAGS += -O3 -mtune=native -march=native -ffast-math -msse -msse2 -mfpmath=both -malign-double -mmmx -m3dnow  -ftree-vectorize -ftree-vectorizer-verbose=1

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
          msdmeasurement.cpp \
          childdynamics.cpp \
          trajectoryplot.cpp \
          fretdynamics.cpp \
          $(LIB)/datatable.cpp \
          $(LIB)/highrestimer.cpp \
          $(LIB)/jkiniparser2.cpp \
          $(LIB)/jkmathparser.cpp \
          $(LIB)/statistics_tools.cpp \
          $(LIB)/tools.cpp \
          $(LIB)/image_tools.cpp \
          $(LIB)/tinytiffreader.cpp \
          $(LIB)/tinytiffwriter.cpp


SRC_FILE_O = $(subst .cpp,.o,$(SRC_FILE))


ifeq ($(findstring Msys,$(OS)),Msys)
PREFIX=/mingw
EXE_SUFFIX=.exe
SO_SUFFIX=.dll
SO_PATH=$(PREFIX)/bin
else
ifeq ($(findstring Windows,$(OS)),Windows)
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
endif


Debug:  ${EXECUTABLE}$(EXE_SUFFIX)

Release:  ${EXECUTABLE}$(EXE_SUFFIX)

${EXECUTABLE}$(EXE_SUFFIX): ${SRC_FILE_O}
	$(CC) $(CFLAGS) -o $(EXECUTABLE)$(EXE_SUFFIX) ${SRC_FILE_O} $(LDFLAGS)

$(SRC_FILE_O): %.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	echo $(OS)
	rm -f *.exe
	rm -f ${EXECUTABLE}$(EXE_SUFFIX)
	rm -f ${SRC_FILE_O}
	rm -f *.o


