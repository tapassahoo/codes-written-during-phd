rm test.exe MATRIX* fort.* b
#ifort -openmp -O3 trifluorocynd.f -o test.exe 
ifort -openmp -O0 -mcmodel medium -shared-intel trifluorocynd.f -o test.exe
./test.exe>b&
