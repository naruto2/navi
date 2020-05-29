PROG = navi
CC   = /opt/cuda/bin/nvcc
SRCS = 	solver.cpp \
	navierstokes.cpp

OBJS = $(SRCS:.cpp=.o)

.SUFFIXES: .cpp .o

LDFLAGS = -Llib -lforeach -lop -lxmesh -lcusolver

all: $(PROG)
$(PROG) : $(OBJS) main.cpp
	(cd lib/op; make)
	(cd lib/foreach; make)
	(cd lib/xmesh; make)
	$(CC) -g main.cpp $(OBJS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) -g -c -I. $<

clean: ;
	(cd lib/op; make clean)
	(cd lib/foreach; make clean)
	(cd lib/xmesh; make clean)
	rm -f $(PROG) $(OBJS)

test: $(PROG)
	../test/test.sh
