OPENMP=

FC=ifort $(OPENMP)
CC=g++
FCOPTS= -O0
CCOPTS= 

LAPACK=

SCALAPACK=

LIBXC_ROOT= 

LIBINT_ROOT=

include ./my_machine.arch

LINKS= -lstdc++ -lm
LIBXC_INCS=-I$(LIBXC_ROOT)/include
LIBXC=-L$(LIBXC_ROOT)/lib -Bstatic -lxcf90 -lxc -Bdynamic
LIBINT=-L$(LIBINT_ROOT)/lib -lint2
LIBINT_INCS=-I$(LIBINT_ROOT)/include/libint2




../molgw: macros.h m_definitions.o m_mpi.o m_timing.o m_warning.o m_calculation_type.o m_tools.o m_scf.o m_atoms.o \
         m_gaussian.o m_basis_set.o m_eri.o m_spectral_function.o m_dft_grid.o gw.o dft_xc.o wrapper_libint.o lebedev_quadrature.o mp2_selfenergy.o mp2.o io.o hamiltonian.o molgw.o
	$(FC) $(FCOPTS) m_definitions.o m_mpi.o m_timing.o m_warning.o m_calculation_type.o m_tools.o m_scf.o m_atoms.o \
                        m_gaussian.o m_basis_set.o m_eri.o m_dft_grid.o m_spectral_function.o \
                        wrapper_libint.o lebedev_quadrature.o mp2_selfenergy.o mp2.o gw.o dft_xc.o io.o hamiltonian.o molgw.o -o ../molgw \
                        $(LAPACK) $(SCALAPACK) \
                        $(LIBXC) $(LIBXC _INCS) \
                        $(LIBINT) $(LIBINT_INCS) \
                        $(LINKS)

molgw.o: macros.h m_definitions.f90 m_timing.f90 m_warning.f90 m_calculation_type.f90 m_tools.f90 m_scf.f90 \
  m_gaussian.f90 m_basis_set.f90  m_eri.f90 m_dft_grid.f90 m_spectral_function.f90 molgw.f90
	$(FC) $(FCOPTS) -c molgw.f90 $(LIBXC_INCS)

hamiltonian.o: macros.h m_definitions.f90 m_timing.f90 m_eri.f90 hamiltonian.f90
	$(FC) $(FCOPTS) -c hamiltonian.f90 

gw.o: macros.h m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_calculation_type.f90 m_spectral_function.f90 m_eri.f90 gw.f90
	$(FC) $(FCOPTS) -c gw.f90 $(LIBXC_INCS)

dft_xc.o: macros.h m_definitions.f90 m_dft_grid.f90 m_tools.f90 m_atoms.f90 dft_xc.f90
	$(FC) $(FCOPTS) -c dft_xc.f90 $(LIBXC_INCS)

io.o: macros.h m_definitions.f90 m_warning.f90 m_calculation_type.f90 m_atoms.f90 m_scf.f90 m_basis_set.f90 io.f90
	$(FC) $(FCOPTS) -c io.f90 $(LIBXC_INCS)

mp2.o: macros.h m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_eri.f90 mp2.f90 
	$(FC) $(FCOPTS) -c mp2.f90

mp2_selfenergy.o: macros.h m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_eri.f90 m_calculation_type.f90 mp2_selfenergy.f90
	$(FC) $(FCOPTS) -c mp2_selfenergy.f90 $(LIBXC_INCS)

lebedev_quadrature.o: macros.h lebedev_quadrature.f
	$(FC) $(FCOPTS) -c lebedev_quadrature.f

wrapper_libint.o: macros.h wrapper_libint.cc
	$(CC) $(CCOPTS) -c wrapper_libint.cc $(LIBINT_INCS)

m_dft_grid.o: macros.h m_definitions.f90 m_dft_grid.f90 m_tools.f90 m_atoms.f90
	$(FC) $(FCOPTS) -c m_dft_grid.f90 $(LIBXC_INCS)

m_spectral_function.o: macros.h m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_calculation_type.f90 m_spectral_function.f90 m_eri.f90
	$(FC) $(FCOPTS) -c m_spectral_function.f90 $(LIBXC_INCS)

m_eri.o: macros.h m_definitions.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_eri.f90
	$(FC) $(FCOPTS) -c m_eri.f90

m_basis_set.o: macros.h m_definitions.f90 m_gaussian.f90 m_basis_set.f90
	$(FC) $(FCOPTS) -c m_basis_set.f90

m_gaussian.o: macros.h m_definitions.f90 m_tools.f90 m_gaussian.f90
	$(FC) $(FCOPTS) -c m_gaussian.f90

m_atoms.o: macros.h m_definitions.f90 m_atoms.f90
	$(FC) $(FCOPTS) -c m_atoms.f90

m_scf.o: macros.h m_definitions.f90 m_scf.f90
	$(FC) $(FCOPTS) -c m_scf.f90

m_tools.o: macros.h m_definitions.f90 m_tools.f90
	$(FC) $(FCOPTS) -c m_tools.f90

m_calculation_type.o: macros.h m_definitions.f90 m_calculation_type.f90
	$(FC) $(FCOPTS) -c m_calculation_type.f90 $(LIBXC_INCS)

m_warning.o: macros.h m_definitions.f90 m_warning.f90
	$(FC) $(FCOPTS) -c m_warning.f90

m_timing.o: macros.h m_definitions.f90 m_timing.f90
	$(FC) $(FCOPTS) -c m_timing.f90

m_mpi.o: macros.h m_definitions.f90 m_mpi.f90
	$(FC) $(FCOPTS) -c m_mpi.f90

m_definitions.o: macros.h m_definitions.f90
	$(FC) $(FCOPTS) -c m_definitions.f90

archive:
	tar czf molgw.tgz ../../molgw/README ../../molgw/src/Makefile ../../molgw/src/*.f90 ../../molgw/src/*.f ../../molgw/src/*.cc \
                          ../../molgw/src/macros.h \
                          ../../molgw/basis/* ../../molgw/util/* ../../molgw/tests/ ../../molgw/src/my_machine.arch
#	tar czf molgw.tgz ../src/README ../src/Makefile ../src/*.f90 ../src/*.f ../src/*.cc ../basis/* ../doc/notes.tex ../doc/figs/* ../util/* ../tests/
	scp molgw.tgz 132.166.68.166:

doc:
	pdflatex ../doc/notes.tex	
	pdflatex ../doc/notes.tex	
#	rm -f notesNotes.bib notes.aux notes.log
	mv notes.pdf ../doc/notes.pdf
clean:
	rm -f *.o *.mod *.log *.bib *.aux
