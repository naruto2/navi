all: obj
a.out: op.hpp op.c++ main.c++
	c++ -Wall *.c++

obj: op.hpp op.c++
	c++ -Wall -c op.c++
	ar cr ../libop.a *.o
	rm *.o


test: a.out
	./a.out -v aaa -f bbb

clean:
	rm -f a.out ../libop.a *.o
