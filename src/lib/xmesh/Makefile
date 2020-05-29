all:  obj
SRCS =	f2mesh.c++ \
	dimp2.c++ \
	delta.c++


OBJS = $(SRCS:.c++=.o)
.SUFFIXES: .c++ .o

CPPFLAGS =
LDFLAGS = -Llib -lop

.PHONY: all clean

$(PROG) : $(OBJS)
	c++ xmesh.c++ $(OBJS) -o $@ $(LDFLAGS)

.c++.o:
	c++ $(CPPFLAGS) -c $< 

obj: $(OBJS)
	ar cr ../libxmesh.a $(OBJS)
	rm *.o

clean :
	rm -f $(OBJS) $(PROG) ../libxmesh.a 


