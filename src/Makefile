navierstokes: main.cpp solver.cpp
	/opt/cuda/bin/nvcc -g solver.cpp navierstokes.cpp main.cpp -Llib  -lforeach -lop  -l cusolver -lxmesh -o navi

clean: ;
	rm -f navi

