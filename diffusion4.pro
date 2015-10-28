TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    alvtools.cpp \
    browniandynamics.cpp \
    childdynamics.cpp \
    datatable.cpp \
    diffusiontools.cpp \
    dynamicsfromfiles2.cpp \
    fcsmeasurement.cpp \
    fluorescenceimaging.cpp \
    fluorescencemeasurement.cpp \
    fluorophordynamics.cpp \
    gnuplot_tools.cpp \
    gridrandomwalkdynamics.cpp \
    highrestimer.cpp \
    image_tools.cpp \
    jkiniparser2.cpp \
    jkmathparser.cpp \
    msdmeasurement.cpp \
    tools.cpp \
    trajectoryplot.cpp \
    nulldynamics.cpp



HEADERS += \
    alvtools.h \
    browniandynamics.h \
    childdynamics.h \
    datatable.h \
    diffusiontools.h \
    dynamicsfromfiles2.h \
    fcsmeasurement.h \
    fluorescenceimaging.h \
    fluorescencemeasurement.h \
    fluorophordynamics.h \
    gnuplot_tools.h \
    gridrandomwalkdynamics.h \
    highrestimer.h \
    histogram.h \
    image_tools.h \
    jkimage.h \
    jkiniparser2.h \
    jkmathparser.h \
    lib_imexport.h \
    msdmeasurement.h \
    teebuf.h \
    textcolor.h \
    ticktock.h \
    tools.h \
    trajectoryplot.h \
    nulldynamics.h

DISTFILES += \
    Makefile

