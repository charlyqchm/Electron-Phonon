CC := icc
CFLAGS := -llapack -O3

objects 	  = atom.o transport_e-ph_sub.o transport_e-ph.o
executable = myprogram

$(executable): $(objects)
			$(CC) -o $@ $^ ${CFLAGS}

atom.o: atom.cpp atom.h
			$(CC) -o $@ -c $< ${CFLAGS}

transport_e-ph_sub.o: transport_e-ph_sub.cpp funct.h
			$(CC) -o $@ -c $< ${CFLAGS}

transport_e-ph.o: transport_e-ph.cpp
			$(CC) -o $@ -c $< ${CFLAGS}

.PHONY: clean
clean:
			rm $(objects) $(executable)
