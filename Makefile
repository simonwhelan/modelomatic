CXX     = g++
CFLAGS  = -DISUNIX -O3
LDFLAGS = 
SRCS=ModelAssess.cxx bionj.cxx       data.cxx        impsamp.cxx     interface.cxx   model.cxx       new_process.cxx optimise.cxx    parsimony.cxx   process.cxx     tools.cxx       tree.cxx        treelist.cxx
OBJS=$(SRCS:%.cxx=%.o) 

all: modelassess

modelassess: $(OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

%.o: %.cxx
	$(CXX) -c $(CFLAGS) $<
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
