mpif77 -openmp -O3 -o abc.exe -I/usr/include/fftw3.h abc.f -L/usr/lib64/libfftw3.so.3.2.4 -lfftw3 -lm -limf -lmpi
mpirun -machinefile hosts -np 4 -x KMP_STACKSIZE=2048000000 abc.exe>b&
