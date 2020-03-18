CC := nvcc
CFLAGS := -llapack -O3 -lcublas --gpu-architecture=sm_61

objects 	  = atom.o transport_e-ph_sub.o transport_e-ph.o matmul_cublas.o
executable = myprogram

$(executable): $(objects)
			$(CC) -o $@ $^ ${CFLAGS}

atom.o: atom.cpp atom.h
			$(CC) -o $@ -c $< ${CFLAGS}

matmul_cublas.o: matmul_cublas.cu matmul_cublas.h
			$(CC) -o $@ -c $< ${CFLAGS}

transport_e-ph_sub.o: transport_e-ph_sub.cpp funct.h
			$(CC) -o $@ -c $< ${CFLAGS}

transport_e-ph.o: transport_e-ph.cpp
			$(CC) -o $@ -c $< ${CFLAGS}

.PHONY: clean
clean:
			rm $(objects) $(executable)
