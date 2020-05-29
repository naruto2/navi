obj: foreach.c++ Rf.c++ ary.c++ forgammap1.c++
	c++ -c foreach.c++
	c++ -c Rf.c++
	c++ -c ary.c++
	c++ -c forgammap1.c++
	ar cr ../libforeach.a *.o
	rm *.o

clean: ;
	rm -f ../libforeach.a *.o
