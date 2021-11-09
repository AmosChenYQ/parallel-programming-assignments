#! /bin/bash -x
mpirun -np 1 ./determinant_mpi 4
mpirun -np 1 ./determinant_mpi 6
mpirun -np 2 ./determinant_mpi 8
mpirun -np 4 ./determinant_mpi 8
mpirun -np 6 ./determinant_mpi 8
mpirun -np 8 ./determinant_mpi 8
mpirun -np 4 ./determinant_mpi 10
mpirun -np 6 ./determinant_mpi 10
mpirun -np 8 ./determinant_mpi 10
mpirun -np 4 ./determinant_thread 12
mpirun -np 6 ./determinant_thread 12
mpirun -np 8 ./determinant_thread 12
# mpirun -np 4 ./determinant_thread 14
# mpirun -np 6 ./determinant_thread 14
# mpirun -np 8 ./determinant_thread 14
# mpirun -np 4 ./determinant_thread 16
# mpirun -np 6 ./determinant_thread 16
# mpirun -np 8 ./determinant_thread 16