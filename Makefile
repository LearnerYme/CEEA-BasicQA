# Compiler: g++
CXX = $(shell root-config --cxx)

# Root Compiler flags
ROOTFLAGS = $(shell root-config --cflags)

VMCWORKDIR ?= $(error VMCWORKDIR is not set!)
FAIRROOTPATH ?= $(error FAIRROOTPATH is not set!)

CEEROOT_INCLUDES = -I${VMCWORKDIR} \
                   -I${VMCWORKDIR}/CeeData \
                   -I${VMCWORKDIR}/CEETPC_Tracker \
                   -I${VMCWORKDIR}/tpc \
                   -I${VMCWORKDIR}/tpc/edm \
                   -I${VMCWORKDIR}/itof \
                   -I${VMCWORKDIR}/itof/edm \
                   -I${VMCWORKDIR}/etof \
                   -I${VMCWORKDIR}/etof/edm \
                   -I${VMCWORKDIR}/zdc \
                   -I${VMCWORKDIR}/zdc/edm \
                   -I${VMCWORKDIR}/mwdc \
                   -I${VMCWORKDIR}/mwdc/edm 

CFLAGS = ${ROOTFLAGS} -g -I${FAIRROOTPATH}/include ${CEEROOT_INCLUDES}

# Root Library
ROOTLIBS = $(shell root-config --libs)

FAIRLIBDIR = ${FAIRROOTPATH}/lib
CEELIBDIR = ${VMCWORKDIR}/build/lib

LDFLAGS = ${ROOTLIBS} \
          -L${FAIRLIBDIR} -lBase \
					-L${CEELIBDIR} -lCeeData -lCeeTpc -lCeeeTOF -lCeeiTOF -lCeeZDC -lCeeMWDC -lCEETPC_Tracker

# Source files
SOURCES = analysis.cxx AnaMaker.cxx

# Object files
OBJDIR = obj
OBJECTS = $(OBJDIR)/analysis.o $(OBJDIR)/AnaMaker.o

# Executable
EXECUTABLE = run

$(shell mkdir -p $(OBJDIR))

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	${CXX} ${CFLAGS} $^ -o$@ ${LDFLAGS}

$(OBJDIR)/analysis.o: analysis.cxx
	${CXX} -fPIC -c ${CFLAGS} $< -o $@

$(OBJDIR)/AnaMaker.o: AnaMaker/AnaMaker.cxx AnaMaker/AnaMaker.h
	${CXX} -fPIC -c ${CFLAGS} $< -o $@

clean:
	rm -vf $(OBJECTS) $(EXECUTABLE)
	rm -rf $(OBJDIR)

.PHONY: clean all
