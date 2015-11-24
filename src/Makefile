.SUFFIXES: .f90 .f

.f90.o:
	$(F90) -fpp $(FOPT) $(INCLUDE) -c $<

F90=ifort
FOPT=-openmp -msse3 -axSSE4.2
INCLUDE=-I$(MKLROOT)/include -I.

LOPT=-mkl=parallel

fit_para.x: model.o fitting.o banddata.o constants.o cg.o
	$(F90) -o fit_para.x main.f90 $(LOPT) model.o fitting.o banddata.o constants.o cg.o -lmkl_lapack95_lp64

test.x: test.f90 model.o fitting.o banddata.o constants.o
	$(F90) -o test.x test.f90 $(LOPT) model.o fitting.o banddata.o constants.o -lmkl_lapack95_lp64

clean:
	rm -f *.x *.mod *.o

constants.o: constants.f90
banddata.o : banddata.f90 constants.o
model.o    : model.f90 constants.o
fitting.o  : fitting.f90 constants.o model.o banddata.o
cg.o       : cg.f
