CXX							=	g++-12

INCDIRS						=	. cnpy fluid input geodesics metric
LIBDIR						=	/usr/local/lib
BINDIR						=	build

DEPFLAGS					=	-MP -MD
LDFLAGS						= 	-L$(LIBDIR)
LDLIBS						= 	-lz
OMPFLAGS					= 	-fopenmp
OPTFLAGS					= 	-O3
CXXFLAGS					=	-fno-common -Wall -Wno-unused-variable -Wno-unused-but-set-variable -std=c++17	\
								$(foreach D,$(INCDIRS),-I$(D)) $(DEPFLAGS) $(OPTFLAGS) #$(OMPFLAGS)

BINARY						=	$(BINDIR)/CurrentSheetGeodesics
SOURCES						=	$(foreach D,$(INCDIRS),$(wildcard $(D)/*.cpp))
OBJECTS						=	$(patsubst %.cpp,%.o,$(SOURCES))
DEPFILES					=	$(patsubst %.cpp,%.d,$(SOURCES))


all:	$(BINARY)

$(BINARY):	$(OBJECTS)
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS) $(CXXFLAGS)

%.o: %.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm $(BINARY) $(OBJECTS) $(DEPFILES)