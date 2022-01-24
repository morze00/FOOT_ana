CXX=clang++
LINKER=clang++

ROOTCONFIG := root-config

CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I${FAIRROOTPATH}/include
CFLAGS += -I$(SIMPATH)/include
CFLAGS += -I$(VMCWORKDIR)/r3bdata
CFLAGS += -I$(VMCWORKDIR)/tracking
CFLAGS += -I$(VMCWORKDIR)/r3bdata/fibData
CFLAGS += -I$(VMCWORKDIR)/field
CFLAGS += -I$(SIMPATH)/include
CFLAGS += -I$(ROOT_INCLUDE_PATH)
CFLAGS += -I$(ROOT_INCLUDE_DIR)
CFLAGS += --std=c++11 -g -O0 -fexceptions
#CFLAGS += -g -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings 

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += -lEG $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(ROOT_LIBRARY_DIR) -L$(FAIRROOTPATH)/lib
LDFLAGS += -L/Users/vpanin/r3broot/build_r3broot/lib -lR3BTracking
LDFLAGS += -L/Users/vpanin/r3broot/build_r3broot/lib -lField
LDFLAGS += -g

INCLUDEDIR=include
DIR_INC=-I$(INCLUDEDIR)
EXEC=Pedestals_Slice_Fit
#OBJ= main.o Propagation.o
OBJ= Pedestals_Slice_Fit.o

HEADERS= libs.hh

DEPS = $(patsubst %,$(INCLUDEDIR)/%,$(HEADERS))

%.o : %.cpp $(DEPS)
	$(MAKEDEPEND)
	${CXX} ${CFLAGS} $(DIR_INC) -c $< -o $@
	echo "	CXX $@"


default: $(OBJ)
	${LINKER} ${LDFLAGS} ${CFLAGS} -o $(EXEC) $(DIR_INC) $(OBJ)
	echo " COMP $(EXEC)"


clean:
	rm -f *.o
	rm -f Pedestals_Slice_Fit
	rm -rf *.dSYM
