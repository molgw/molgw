OPENMP= -openmp -DOPENMP

FC=ifort -fpp $(OPENMP) -DHAVE_LIBXC -DMOLECULES #-DCHI0  #-DSYMMETRIZED #-DCHI0 # -DAUXIL_BASIS
#FC=ifort -fpp  -DHAVE_LIBXC -DCHI0 # -DAUXIL_BASIS
#FC=ifort -fpp -DHAVE_LIBXC  #-DCHI0 # -DAUXIL_BASIS
#FC=gfortran -cpp  -ffree-line-length-none
CC=g++
FCOPTS= -O0 -g #-check all
#FCOPTS= -O3

#LAPACK= -L/home/bruneval/src/lapack-3.3.1/ -llapack_LINUX -lblas_LINUX 
#LAPACK=  -Wl,--start-group  $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a $(MKLROOT)/lib/em64t/libmkl_intel_thread.a $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -openmp -lpthread
#LAPACK=  -Wl,--start-group  $(MKLROOT)/lib/em64t/libmkl_intel_lp64.a $(MKLROOT)/lib/em64t/libmkl_sequential.a $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread
LAPACK= -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack95_lp64 -lmkl_core

LIBXC_INCS= -I/home/bruneval/src/libxc-1.0_ifort/include
LIBXC=  -L/home/bruneval/src/libxc-1.0_ifort/lib -Bstatic -lxc -Bdynamic

LINKS= -lstdc++ -lm
LIBINT_ROOT=/home/bruneval/src/libint-2.0.0_install/
LIBINT=-L$(LIBINT_ROOT)/lib -lint2
LIBINT_INCS=-I$(LIBINT_ROOT)/include/libint2

../atom: m_definitions.o m_timing.o m_warning.o m_calculation_type.o m_tools.o \
       m_gaussian.o m_basis_set.o m_dft.o m_eri.o m_gw.o wrapper_libint.o lebedev_quadrature.o mp2_selfenergy.o mp2.o io.o atom.o
	$(FC) $(FCOPTS) m_definitions.o m_timing.o m_warning.o m_calculation_type.o m_tools.o m_gaussian.o m_basis_set.o m_gw.o m_dft.o m_eri.o \
                        wrapper_libint.o lebedev_quadrature.o mp2_selfenergy.o mp2.o io.o atom.o -o ../atom \
                        $(LAPACK) \
                        $(LIBXC) $(LIBXC _INCS) \
                        $(LIBINT) $(LIBINT_INCS) \
                        $(LINKS)

atom.o: m_definitions.f90 m_timing.f90 m_warning.f90 m_calculation_type.f90 m_tools.f90 m_gaussian.f90 m_basis_set.f90 m_gw.f90 m_dft.f90 m_eri.f90 atom.f90
	$(FC) $(FCOPTS) -c atom.f90 

io.o: m_definitions.f90 m_warning.f90 io.f90
	$(FC) $(FCOPTS) -c io.f90 

m_eri.o: m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_eri.f90
	$(FC) $(FCOPTS) -c m_eri.f90

mp2.o: m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_eri.f90 mp2.f90 
	$(FC) $(FCOPTS) -c mp2.f90

mp2_selfenergy.o: m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_eri.f90 mp2_selfenergy.f90
	$(FC) $(FCOPTS) -c mp2_selfenergy.f90

lebedev_quadrature.o: lebedev_quadrature.f
	$(FC) $(FCOPTS) -c lebedev_quadrature.f

wrapper_libint.o: wrapper_libint.cc
	$(CC) -c wrapper_libint.cc $(LIBINT_INCS)

m_gw.o: m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_calculation_type.f90 m_gw.f90 m_eri.f90
	$(FC) $(FCOPTS) -c m_gw.f90

m_dft.o: m_definitions.f90 m_dft.f90
	$(FC) $(FCOPTS) -c m_dft.f90 $(LIBXC_INCS)

m_basis_set.o: m_definitions.f90 m_gaussian.f90 m_basis_set.f90
	$(FC) $(FCOPTS) -c m_basis_set.f90

m_gaussian.o: m_definitions.f90 m_tools.f90 m_gaussian.f90
	$(FC) $(FCOPTS) -c m_gaussian.f90

m_tools.o: m_definitions.f90 m_tools.f90
	$(FC) $(FCOPTS) -c m_tools.f90

m_calculation_type.o: m_definitions.f90 m_calculation_type.f90
	$(FC) $(FCOPTS) -c m_calculation_type.f90

m_warning.o: m_definitions.f90 m_warning.f90
	$(FC) $(FCOPTS) -c m_warning.f90

m_timing.o: m_definitions.f90 m_timing.f90
	$(FC) $(FCOPTS) -c m_timing.f90

m_definitions.o: m_definitions.f90
	$(FC) $(FCOPTS) -c m_definitions.f90

archive:
	tar czf atom.tgz ../../atom/README ../../atom/src/Makefile ../../atom/src/*.f90 ../../atom/src/*.f ../../atom/src/*.cc ../../atom/basis/* ../../atom/util/* ../../atom/tests/
#	tar czf atom.tgz ../src/README ../src/Makefile ../src/*.f90 ../src/*.f ../src/*.cc ../basis/* ../doc/notes.tex ../doc/figs/* ../util/* ../tests/
	scp atom.tgz 132.166.68.166:

doc:
	pdflatex ../doc/notes.tex	
	pdflatex ../doc/notes.tex	
#	rm -f notesNotes.bib notes.aux notes.log
	mv notes.pdf ../doc/notes.pdf
clean:
	rm -f *.o *.mod *.log *.bib *.aux
