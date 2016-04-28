

CXX     = g++ # STANDARD COMPILER; OSX + LINUX
CC      = gcc
#CXX     = x86_64-w64-mingw32-g++ -static        # CROSS PLATFORM COMPILER WINDOWS
#CC      = x86_64-w64-mingw32-gcc -static        # CROSS PLATFORM COMPILER WINDOWS
#CXX      = icpc
#CC       = icc
CFLAGS  = -DISUNIX -O3 -w
#CFLAGS = -DISUNIX -g -w
LDFLAGS = 
SRCS=ding.cxx modelomatic.cxx bionj.cxx       data.cxx        impsamp.cxx     interface.cxx   model.cxx       new_process.cxx optimise.cxx    parsimony.cxx   process.cxx     tools.cxx       tree.cxx        treelist.cxx ini/cpp/INIReader.cxx codon_model.cxx
OBJS=$(SRCS:%.cxx=%.o) 
OBJS+=ini/ini.o

all:	modelomatic 

modelomatic: $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.cxx: %.cpp
	cp $< $@

%.o: %.cxx
	$(CXX) -c $(CFLAGS) $< -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CFLAGS) -MM $^ >> ./.depend;

include .depend

.PHONY: clean
clean:
		rm -f *.o
		rm -f ini/*.o
		rm -f ini/cpp/*.o
		rm -f modelomatic
		rm -f .depend
