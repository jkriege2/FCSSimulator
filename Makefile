all: Release
#debug: Debug
#release: Release

LIB=../../../LIB/trunk
CC=g++

CFLAGS =  -Wall -DQF_DONT_USE_ALIGNED_MALLOC -DNO_LIBTIFF -DJKIMAGE_USES_TINYTIFF -I. -I./extlibs/TinyTIFF/ -I./extlibs/StatisticsTools/  -fexceptions
#--enable-auto-import
LDFLAGS = -lm
CFLAGS += -I../../extlibs/gsl/include/ -I./extlibs/gsl/include/
LDFLAGS += -L../../extlibs/gsl/lib/ -L./extlibs/gsl/lib/ -lgsl -lgslcblas


# Release: CFLAGS += -O2 -mtune=native -march=native -ffast-math -msse -msse2 -mfpmath=sse -malign-double -ftree-vectorize -ftree-vectorizer-verbose=0

Release: CFLAGS += -O2 -mtune=native -march=native -ffast-math -msse -msse2 -mfpmath=sse -ftree-vectorize 

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
          nulldynamics.cpp \
          trajectoryplot.cpp \
          datatable.cpp \
          highrestimer.cpp \
          jkiniparser2.cpp \
          jkmathparser.cpp \
          extlibs/StatisticsTools/statistics_tools.cpp \
          tools.cpp \
          image_tools.cpp \
          extlibs/TinyTIFF/tinytiffreader.cpp \
          extlibs/TinyTIFF/tinytiffwriter.cpp \
          gnuplot_tools.cpp \
          alvtools.cpp


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


