mpif77 -openmp -O0 -mcmodel medium -shared-intel trifluorocyndnew.f -limf -lmpi
mpirun -machinefile hosts -np 8 -x KMP_STACKSIZE=2048000000 a.out>b&
