FC = gfortran
FFLAGS = -march=native -mtune=native -O3
FDEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal
# IDIR = -I/usr/local/include/fgsl/
# LIBS = -lgsl -lfgsl -lopenblas -lm

objects: solver_rk4.exe

all: $(objects)

solver_rk4.exe: rk4.f90 parameters.o
	$(FC) $(FFLAGS) -o $@ $(IDIR) $^ $(LIBS)

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $(IDIR) $<

debug: rk4.f90 parameters.f90
	$(FC) -g $(FDEBUGFLAGS) -o dbg_rk4.exe $(IDIR) $^ $(LIBS)

clean:
	rm -f *.o *.mod *.exe
