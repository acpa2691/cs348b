PBRT_TOOLS = exrassemble exravg exrtotiff tifftoexr ply2pbrt samplepat

ARCH = $(shell uname)

PBRTDIR=..
PBRTLIB=$(PBRTDIR)/core/libpbrt.a
EXRDIR=/afs/ir/class/cs348b/exr
TIFFDIR=/afs/ir/class/cs348b/tiff
LIBS = -L/usr/lib -L/usr/local/lib -L$(EXRDIR)/lib -L$(TIFFDIR)/lib
INCLUDES = -I/usr/include -I/usr/local/include -I$(EXRDIR)/include/OpenEXR -I$(TIFFDIR)/include
OPT=-g
CXXFLAGS = $(OPT) $(INCLUDES) -I$(PBRTDIR)/core

ifeq ($(ARCH),Darwin)
LDFLAGS = $(LIBS)
LDFLAGS += -Wl,-all_load $(PBRTDIR)/core/libpbrt.a 
CXXFLAGS += -Wno-long-double
else
LDFLAGS = $(LIBS)
LDFLAGS += -Wl,--export-dynamic -Wl,-whole-archive $(PBRTLIB) -Wl,-no-whole-archive 
endif

LDFLAGS += -lIlmImf -lImath -lIex -lHalf -lz -ltiff -ldl

default: $(PBRT_TOOLS)

%: %.cpp $(PBRTLIB)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

ply2pbrt: ply.h ply.c ply2pbrt.c
	$(CC) $(OPT) -I. -o ply2pbrt ply2pbrt.c ply.c -lm
clean:
	rm -f $(PBRT_TOOLS) *.o *~
