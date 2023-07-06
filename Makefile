CCC = g++       # compiler
DICT = rootcint # root dictionary-generator
AR = ar         # archiver
LD = g++        # linker

COAST_DIR=/net/minus/home/shchegolev/CORSIKA/corsika-75600

ARFLAGS = rcs

LDFLAGS += -fPIC -ggdb3 -Wl,--no-as-needed
LDFLAGS += $(shell root-config --libs)

LDFLAGS += -L$(COAST_DIR)/lib -L$(COAST_DIR)/lib/unknown
LDFLAGS += -lCorsikaFileIO
LDFLAGS += -lCorsikaIntern
LDFLAGS += -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6/

CPPFLAGS += -c -fPIC  -ggdb3 -std=gnu++0x -lpthread
CPPFLAGS += -I$(COAST_DIR)/include
CPPFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

#EXE = CorsikaPlotter
EXE = Cors_P64nmuh.e


all: ${EXE}

curvout: ${CURVOUT}

.cc.o:
	${CCC} ${CPPFLAGS} $^

#${EXE}: CorsikaPlotter.o
${EXE}: CorsikaReader.o
	${CCC} $^ -o $@ ${LDFLAGS} generator.o -std=gnu++0x -lpthread
	#${CCC} $^ -o $@ ${LDFLAGS} P64_muh.o -lg2c


clean: 
	@rm -f *.o *.a *.so *~ core *.root $(EXE)
