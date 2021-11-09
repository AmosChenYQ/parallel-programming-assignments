# parallel-programming-assignments

## Assignment 1

Make executable
```bash
make
```

Clean logs

```bash
make clean
```

Submit job under TorQue

```bash
qsub pthreads
```

Or you can use `nohub` to run bash

```bash
nohub ./run_batch.sh > threads.log &
```

## Assignment 2

Make executable
```bash
make
```

Clean logs

```bash
make clean
```

Submit job under TorQue

```bash
qsub mpi
```

Or you can use `nohub` to run mpirun directly

```bash
 mpirun -np 4 ./determinant_mpi 12 >& mpi.log
```

