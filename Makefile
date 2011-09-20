LIB=../../../LIB/trunk
CC=g++
CFLAGS =   -Wall -O2 #-march=core2 -msse4 -mcx16 -mpopcnt -msahf -O2 -pipe -ggdb -Wall
LDFLAGS = -lgsl -lgslcblas -lm

EXECUTABLE=diffusion4
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


${EXECUTABLE}: ${SRC_FILE_O}
	$(CC) $(CFLAGS) -o $(EXECUTABLE) ${SRC_FILE_O} $(LDFLAGS)

$(SRC_FILE_O): %.o: %.cpp
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f ${EXECUTABLE} *.o
	rm -f $(SRC_FILE_O)



