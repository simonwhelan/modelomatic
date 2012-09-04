CXX     = g++
#CXX     = icpc
CFLAGS  = -DISUNIX -O3
LDFLAGS = 
SRCS=ModelAssess.cxx bionj.cxx       data.cxx        impsamp.cxx     interface.cxx   model.cxx       new_process.cxx optimise.cxx    parsimony.cxx   process.cxx     tools.cxx       tree.cxx        treelist.cxx ini/cpp/INIReader.cxx
OBJS=$(SRCS:%.cxx=%.o) 
OBJS+=ini/ini.o

all: modelassess

modelassess: $(OBJS)
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
		rm -f modelassess
		rm -f .depend
