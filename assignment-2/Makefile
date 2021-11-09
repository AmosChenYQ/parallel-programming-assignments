MPICXX= mpic++
MPICXXFLAGS=-O2 -std=c++11 -Wall

all: determinant_mpi

determinant_mpi: determinant_mpi.cpp
	$(MPICXX) $(MPICXXFLAGS) determinant_mpi.cpp -o determinant_mpi

clean:
	rm -rf determinant_mpi
