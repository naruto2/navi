PROG = navi
CC   = /opt/cuda/bin/nvcc
SRCS =  op.cpp \
	solver.cpp \
	foreach.cpp \
	f2mesh.cpp \
	navi.cpp

OBJS = $(SRCS:.cpp=.o)

.SUFFIXES: .cpp .o

LDFLAGS = -lcusolver

all: $(PROG)
$(PROG) : $(OBJS) main.cpp
	$(CC) -g main.cpp $(OBJS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) -g -c -I. $<

clean: ;
	rm -f $(PROG) $(OBJS) 0000?.gz

test: $(PROG)
	./navi -test
	../test/test.sh
