FC = gfortran
FFLAGS = -march=native -mtune=native -O3
FDEBUGFLAGS = -g -fbacktrace -ffpe-trap=zero,overflow,underflow,denormal
# IDIR = -I/usr/local/include/fgsl/
# LIBS = -lgsl -lfgsl -lopenblas -lm

objects: blasius.exe

all: $(objects)

blasius.exe: blasius.f90 parameters.o
	$(FC) $(FFLAGS) -o $@ $(IDIR) $^ $(LIBS)

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $@ $(IDIR) $<

debug: blasius.f90 parameters.f90
	$(FC) -g $(FDEBUGFLAGS) -o dbg_blasius.exe $(IDIR) $^ $(LIBS)

clean:
	rm -f *.o *.mod *.exe
