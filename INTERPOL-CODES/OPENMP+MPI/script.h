rm a.out b fort.*
mpif77 -openmp -O0 -mcmodel medium -shared-intel splin-2d-int1p.f -limf -lmpi
mpirun -machinefile hosts -np 3 -x KMP_STACKSIZE=20480000000 a.out>b&
