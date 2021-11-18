.PHONY: clean all

FC = gfortran

OBJS = function_NMCT.o NMCT.o complex_bessel.o

target = ./../NMCT.exe

$(target): $(OBJS)
	$(FC) $(OBJS) ./../lapack-3.10.0/liblapack.a  ./../lapack-3.10.0/librefblas.a -o $@

%.o: %.f90
	$(FC) -c $< -o $@
	
complex_bessel.o: complex_bessel.f
	$(FC) -c $< -o $@
    
clean:
	-rm -rf *.o *.mod *.exe
