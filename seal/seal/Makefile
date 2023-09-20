HEADERS = bootstrapping.h contextchain.h karatsuba.h polyarith.h slots.h transform.h
SOURCES = bootstrapping.cpp contextchain.cpp polyarith.cpp slots.cpp transform.cpp seal.cpp
SEAL_DIR = /home/feanor/seal_lib

galois_bootstrapping:
	g++ -I$(SEAL_DIR)/include/SEAL-4.1 -L$(SEAL_DIR)/lib $(SOURCES) -lseal-4.1 -o galois_bootstrapping

clean:
	rm -f galois_bootstrapping