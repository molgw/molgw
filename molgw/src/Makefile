OPENMP=

FC=ifort 
CC=g++
FCOPTS= -O0
CCOPTS= 

LAPACK=

SCALAPACK=

LIBXC_ROOT= 

LIBINT_ROOT=

include ./my_machine.arch

FC_W_OPTS= $(FC) $(FCOPTS) $(OPENMP)
CC_W_OPTS= $(CC) $(CCOPTS)

LINKS= -lstdc++ -lm
LIBXC_INCS=-I$(LIBXC_ROOT)/include
LIBXC=-L$(LIBXC_ROOT)/lib -Bstatic -lxcf90 -lxc -Bdynamic
LIBINT=-L$(LIBINT_ROOT)/lib -lint2
LIBINT_INCS=-I$(LIBINT_ROOT)/include/libint2


.PHONY: all clean archive doc input prepare_sourcecode


####################################################
# Parallel make execution requires a bit of ordering
all: 
	+make  prepare_sourcecode
	+make  level000
	+make  level005
	+make  level010
	+make  level012
	+make  level015
	+make  level020
	+make  level025
	+make  level030
	+make  level032
	+make  level035
	+make  level040
	+make  ../molgw

level000: m_definitions.o libint_wrapper.o libint_wrapper_nounderscore.o lebedev_quadrature.o pdbssolver1.o

level005: m_warning.o m_gos.o

level010: m_elements.o m_timing.o m_memory.o m_tools.o 

level012: m_mpi.o

level015: m_atoms.o m_gaussian.o m_block_diago.o

level020: m_basis_set.o

level025: m_inputparam.o

level030: m_eri.o m_scf.o m_dft_grid.o

level032: m_eri_calculate.o m_eri_lr_calculate.o m_eri_ao_mo.o

level035: m_spectral_function.o m_hamiltonian.o m_hamiltonian_sca.o m_hamiltonian_dist.o mp2_selfenergy.o

level040: linear_response.o static_polarizability.o gw.o bse.o dft_xc.o io.o mp2.o scf_loop.o molgw.o



####################################################



../molgw: m_definitions.o libint_wrapper.o libint_wrapper_nounderscore.o lebedev_quadrature.o pdbssolver1.o \
        m_warning.o m_gos.o \
        m_elements.o m_timing.o m_memory.o m_tools.o m_mpi.o \
        m_atoms.o m_gaussian.o m_block_diago.o \
        m_basis_set.o \
        m_inputparam.o \
        m_eri.o m_scf.o m_dft_grid.o \
        m_eri_calculate.o m_eri_lr_calculate.o m_eri_ao_mo.o \
        m_spectral_function.o m_hamiltonian.o m_hamiltonian_sca.o m_hamiltonian_dist.o mp2_selfenergy.o \
        linear_response.o static_polarizability.o gw.o bse.o dft_xc.o io.o mp2.o  \
        scf_loop.o molgw.o
	$(FC_W_OPTS) \
                        m_definitions.o m_warning.o m_elements.o m_mpi.o m_timing.o m_memory.o m_inputparam.o m_tools.o m_scf.o m_atoms.o \
                        m_gos.o m_gaussian.o m_block_diago.o m_basis_set.o m_eri.o m_dft_grid.o m_eri_calculate.o m_eri_lr_calculate.o m_eri_ao_mo.o \
                        m_spectral_function.o \
                        libint_wrapper.o libint_wrapper_nounderscore.o lebedev_quadrature.o mp2_selfenergy.o \
                        mp2.o gw.o bse.o linear_response.o static_polarizability.o dft_xc.o io.o \
                        m_hamiltonian.o m_hamiltonian_sca.o m_hamiltonian_dist.o scf_loop.o pdbssolver1.o  \
                        molgw.o \
                        -o ../molgw \
                        $(SCALAPACK) $(LAPACK)   \
                        $(LIBXC) $(LIBXC _INCS)  \
                        $(LIBINT) $(LIBINT_INCS) \
                        $(LINKS)

molgw.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_memory.f90 m_inputparam.f90 m_tools.f90 m_scf.f90 \
  m_gos.f90 m_gaussian.f90 m_basis_set.f90  m_eri.f90 m_dft_grid.f90 m_spectral_function.f90 m_eri_calculate.f90 m_eri_lr_calculate.f90 m_eri_ao_mo.o \
  m_hamiltonian.f90 m_hamiltonian_sca.f90 m_hamiltonian_dist.o \
  molgw.f90
	$(FC_W_OPTS) -c molgw.f90 $(LIBXC_INCS)

scf_loop.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_inputparam.f90 m_tools.f90 m_scf.f90 \
            m_gaussian.f90 m_basis_set.f90  m_eri.f90 m_dft_grid.f90 m_spectral_function.f90 m_eri_calculate.f90 m_eri_lr_calculate.f90 m_eri_ao_mo.o \
            m_hamiltonian.f90 m_hamiltonian_sca.f90 m_hamiltonian_dist.o scf_loop.f90
	$(FC_W_OPTS) -c scf_loop.f90 $(LIBXC_INCS)

gw.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_spectral_function.f90 m_eri.f90 m_eri_calculate.f90 m_eri_ao_mo.f90 \
      m_inputparam.f90 gw.f90
	$(FC_W_OPTS) -c gw.f90 $(LIBXC_INCS)

bse.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_tools.f90 m_basis_set.f90 m_spectral_function.f90 m_eri.f90 m_inputparam.f90 \
        bse.f90
	$(FC_W_OPTS) -c bse.f90 $(LIBXC_INCS)

static_polarizability.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_tools.f90 m_basis_set.f90 m_spectral_function.f90 m_eri.f90 m_inputparam.f90 \
        static_polarizability.f90
	$(FC_W_OPTS) -c static_polarizability.f90 $(LIBXC_INCS)

linear_response.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_block_diago.f90 m_basis_set.f90 m_spectral_function.f90 m_eri.f90 \
                   m_hamiltonian.f90 linear_response.f90
	$(FC_W_OPTS) -c linear_response.f90 $(LIBXC_INCS)

m_hamiltonian_dist.o: m_definitions.f90 m_mpi.f90 m_warning.f90 m_timing.f90 m_eri.f90 m_basis_set.f90 m_hamiltonian_dist.f90
	$(FC_W_OPTS) -c m_hamiltonian_dist.f90 $(LIBXC_INCS)

m_hamiltonian_sca.o: m_definitions.f90 m_mpi.f90 m_warning.f90 m_timing.f90 m_eri.f90 m_basis_set.f90 m_hamiltonian_sca.f90
	$(FC_W_OPTS) -c m_hamiltonian_sca.f90 $(LIBXC_INCS)

m_hamiltonian.o: m_definitions.f90 m_mpi.f90 m_warning.f90 m_timing.f90 m_eri.f90 m_basis_set.f90 m_hamiltonian.f90
	$(FC_W_OPTS) -c m_hamiltonian.f90 $(LIBXC_INCS)

dft_xc.o: m_definitions.f90 m_dft_grid.f90 m_tools.f90 m_atoms.f90 m_mpi.f90 dft_xc.f90
	$(FC_W_OPTS) -c dft_xc.f90 $(LIBXC_INCS)

io.o: m_definitions.f90 m_mpi.f90 m_warning.f90 m_atoms.f90 m_scf.f90 m_basis_set.f90 io.f90 revision.f90 m_hamiltonian.f90
	$(FC_W_OPTS) -c io.f90 $(LIBXC_INCS)

mp2.o: m_definitions.f90 m_timing.f90 m_warning.f90 m_mpi.f90 m_basis_set.f90 m_eri.f90 mp2.f90 
	$(FC_W_OPTS) -c mp2.f90 $(LIBXC_INCS)

mp2_selfenergy.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_eri.f90 mp2_selfenergy.f90
	$(FC_W_OPTS) -c mp2_selfenergy.f90 $(LIBXC_INCS)

pdbssolver1.o: pdbssolver1.f
	$(FC_W_OPTS) -c pdbssolver1.f $(LIBXC_INCS)

lebedev_quadrature.o: lebedev_quadrature.f
	$(FC_W_OPTS) -c lebedev_quadrature.f $(LIBXC_INCS)

libint_wrapper.o: libint_wrapper.c
	$(CC_W_OPTS) -c libint_wrapper.c $(LIBINT_INCS)

libint_wrapper_nounderscore.o: libint_wrapper_nounderscore.c
	$(CC_W_OPTS) -c libint_wrapper_nounderscore.c $(LIBINT_INCS)

m_dft_grid.o: m_definitions.f90 m_memory.f90 m_inputparam.f90 m_tools.f90 m_atoms.f90 m_mpi.f90 m_basis_set.f90 m_dft_grid.f90
	$(FC_W_OPTS) -c m_dft_grid.f90 $(LIBXC_INCS)

m_spectral_function.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_basis_set.f90 m_spectral_function.f90 m_eri.f90
	$(FC_W_OPTS) -c m_spectral_function.f90 $(LIBXC_INCS)

m_eri_ao_mo.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_basis_set.f90 m_inputparam.f90 \
                   m_eri.f90 m_eri_ao_mo.f90
	$(FC_W_OPTS) -c m_eri_ao_mo.f90 $(LIBXC_INCS)

m_eri_lr_calculate.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_basis_set.f90 m_inputparam.f90 \
                   m_eri.f90 m_eri_lr_calculate.f90
	$(FC_W_OPTS) -c m_eri_lr_calculate.f90 $(LIBXC_INCS)

m_eri_calculate.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_basis_set.f90 m_inputparam.f90 \
                   m_eri.f90 m_eri_calculate.f90
	$(FC_W_OPTS) -c m_eri_calculate.f90 $(LIBXC_INCS)

m_eri.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_basis_set.f90 m_inputparam.f90 m_eri.f90 
	$(FC_W_OPTS) -c m_eri.f90 $(LIBXC_INCS)

m_basis_set.o: m_definitions.f90 m_elements.f90 m_gaussian.f90 m_mpi.f90 m_basis_set.f90 
	$(FC_W_OPTS) -c m_basis_set.f90 $(LIBXC_INCS)

m_gaussian.o: m_definitions.f90 m_tools.f90 m_gos.f90 m_mpi.f90 m_gaussian.f90
	$(FC_W_OPTS) -c m_gaussian.f90 $(LIBXC_INCS)

m_scf.o: m_definitions.f90 m_inputparam.f90 m_scf.f90
	$(FC_W_OPTS) -c m_scf.f90 $(LIBXC_INCS)

m_inputparam.o: m_definitions.f90 m_elements.f90 m_warning.f90 m_atoms.f90 m_basis_set.f90 m_inputparam.f90 input_variables.f90 basis_path.f90
	$(FC_W_OPTS) -c m_inputparam.f90 $(LIBXC_INCS)

m_block_diago.o: m_definitions.f90 m_warning.f90 m_memory.f90 m_mpi.f90 m_tools.f90 m_block_diago.f90
	$(FC_W_OPTS) -c m_block_diago.f90 $(LIBXC_INCS)

m_atoms.o: m_definitions.f90 m_warning.f90 m_elements.f90 m_atoms.f90
	$(FC_W_OPTS) -c m_atoms.f90 $(LIBXC_INCS)

m_mpi.o: m_definitions.f90 m_warning.f90 m_mpi.f90
	$(FC_W_OPTS) -c m_mpi.f90 $(LIBXC_INCS)

m_tools.o: m_definitions.f90 m_warning.f90 m_tools.f90
	$(FC_W_OPTS) -c m_tools.f90 $(LIBXC_INCS)

m_memory.o: m_definitions.f90 m_warning.f90 m_memory.f90
	$(FC_W_OPTS) -c m_memory.f90 $(LIBXC_INCS)

m_timing.o: m_definitions.f90 m_warning.f90 m_timing.f90
	$(FC_W_OPTS) -c m_timing.f90 $(LIBXC_INCS)

m_elements.o: m_definitions.f90 m_warning.f90 m_elements.f90
	$(FC_W_OPTS) -c m_elements.f90 $(LIBXC_INCS)

m_gos.o: m_definitions.f90 m_gos.f90
	$(FC_W_OPTS) -c m_gos.f90 $(LIBXC_INCS)

m_warning.o: m_definitions.f90 m_warning.f90
	$(FC_W_OPTS) -c m_warning.f90 $(LIBXC_INCS)

m_definitions.o: m_definitions.f90
	$(FC_W_OPTS) -c m_definitions.f90 $(LIBXC_INCS)

archive:
	tar czf molgw.tgz ../../molgw/README.md ../../molgw/src/Makefile \
                          ../../molgw/src/*.py ../../molgw/src/*.f90 ../../molgw/src/*.f ../../molgw/src/*.c \
                          ../../molgw/basis/* ../../molgw/util/* ../../molgw/tests/ ../../molgw/src/my_machine.arch

prepare_sourcecode:
	python prepare_sourcecode.py

input:
	python ../util/input_variables.py

doc:
	pdflatex ../doc/notes.tex	
	pdflatex ../doc/notes.tex	
	rm -f notes.aux notes.log notesNotes.bib
	mv notes.pdf ../doc/notes.pdf

clean:
	rm -f *.o *.mod *.log *.bib *.aux
