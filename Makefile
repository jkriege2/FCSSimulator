LIB=../../../../LIB/trunk
CC=g++
OPTIMIZATION_OPTS = -O2
CPPFLAGS= -Wall $(OPTIMIZATION_OPTS) -I${LIB}
LDFLAGS = -lgsl -lcblas -lm

EXECUTEABLE=diffusion4
SRC_FILE= browniandynamics.cpp \
          diffusiontools.cpp \
          dynamicsfromfiles2.cpp \
          fcsmeasurement.cpp \
          fluorescenceimaging.cpp \
          fluorescencemeasurement.cpp \
          fluorophordynamics.cpp \
          main.cpp \
          $(LIB)/datatable.cpp \
          $(LIB)/highrestimer.cpp \
          $(LIB)/jkiniparser2.cpp \
          $(LIB)/jkmathparser.cpp \
          $(LIB)/statistics_tools.cpp \
          $(LIB)/tools.cpp \


SRC_FILE_O = $(subst .cpp,.o,$(SRC_FILE))


${EXECUTEABLE}: ${SRC_FILE_O}
	$(CC) $(CPPFLAGS) $(LDFLAGS) -o $(EXECUTABLE) ${SRC_FILE_O}

$(SRC_FILE_O): %.o: %.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	rm -f ${EXECUTEABLE} *.o
	rm -f $(SRC_FILE_O)



