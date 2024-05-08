objs = quadpack.o prequel.o rkns.o setup.o parameters.o nsq.o

flags = -Wall -Wno-conversion -Wno-unused-variable -Wno-unused-dummy-argument -g -fcheck=all # LAPACK/BLAS/QUADPACK calls may produce INF/NAN and appear to handle such situations appropriately. Do not compile to catch IEE FPEs

link_flags = -llapack -lblas

comp = gfortran

all : $(objs)
	$(comp) $(flags) wfn_plot.f90 -o wfn_plot nsq.o setup.o rkns.o prequel.o parameters.o quadpack.o $(link_flags)
	$(comp) $(flags) nsq_main.f90 -o nsq_main nsq.o setup.o rkns.o prequel.o parameters.o quadpack.o $(link_flags)

quadpack.o: quadpack.f90
	$(comp) $(flags) -c quadpack.f90 -o quadpack.o 

parameters.o: params.f90
	$(comp) $(flags) -c params.f90 -o parameters.o

prequel.o: prequel.f90 parameters.o
	$(comp) $(flags) -c prequel.f90 -o prequel.o 

setup.o: setup_H.f90 parameters.o
	$(comp) $(flags) -c setup_H.f90 -o setup.o

rkns.o: prequel.o quadpack.o rkns.f90 parameters.o
	$(comp) $(flags) -c rkns.f90 -o rkns.o

nsq.o: nsq.f90 parameters.o
	$(comp) $(flags) -c nsq.f90 -o nsq.o $(link_flags)
