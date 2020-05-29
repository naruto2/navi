PROG = navi
CC   = /opt/cuda/bin/nvcc
SRCS =  op.cpp \
	solver.cpp \
	navi.cpp

OBJS = $(SRCS:.cpp=.o)

.SUFFIXES: .cpp .o

LDFLAGS = -Llib -lforeach -lxmesh -lcusolver

all: $(PROG)
$(PROG) : $(OBJS) main.cpp
	(cd lib/foreach; make)
	(cd lib/xmesh; make)
	$(CC) -g main.cpp $(OBJS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) -g -c -I. $<

clean: ;
	(cd lib/foreach; make clean)
	(cd lib/xmesh; make clean)
	rm -f $(PROG) $(OBJS) 0000?.gz

test: $(PROG)
	./navi -test
	../test/test.sh
