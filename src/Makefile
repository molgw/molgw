# This file is part of MOLGW
# Author: Fabien Bruneval

PYTHON=python
ifeq ($(FC),f77)
  FC=gfortran
endif
ifeq ($(FCFLAGS),)
  FCFLAGS=-O2
endif
ifeq ($(CXXFLAGS),)
  CXXFLAGS=-O2
endif
ifeq ($(CPPFLAGS),)
  CPPFLAGS=-cpp -DHAVE_LIBXC -DHAVE_LIBINT_ONEBODY -DFORTRAN2008
endif
ifeq ($(LAPACK),)
  LAPACK=-lblas -llapack
endif
#SCALAPACK=
#FFTW_ROOT=
#ELPA=
#ELPA_INCS=
#LIBXC_ROOT=
#LIBINT_ROOT=


-include ./my_machine.arch

FC_W_OPTS= $(FC) $(FCFLAGS) $(OPENMP) $(CPPFLAGS)
CXX_W_OPTS= $(CXX) $(CXXFLAGS) $(CPPFLAGS)
LINKER= $(LIBTOOL) $(FC_W_OPTS)

LINKS= -lstdc++ -lm
ifeq ($(LIBXC_ROOT),)
  LIBXC_INCS=
  ifneq (,$(findstring -DHAVE_LIBXC,$(CPPFLAGS)))
    LIBXC= -lxc
  else
    LIBXC=
  endif
else
  LIBXC_INCS=-I$(LIBXC_ROOT)/include
  LIBXC=-L$(LIBXC_ROOT)/lib -Bstatic -lxc -Bdynamic
endif
ifeq ($(LIBINT_ROOT),)
  LIBINT=-lint2
  LIBINT_INCS=
else
  LIBINT=-L$(LIBINT_ROOT)/lib -lint2
  LIBINT_INCS=-I$(LIBINT_ROOT)/include/libint2 -I$(LIBINT_ROOT)/include
endif
ifeq ($(FFTW_ROOT),)
  FFTW=
  FFTW_INCS=
else
  FFTW=-L$(FFTW_ROOT)/lib/ -lfftw3
  FFTW_INCS=-I$(FFTW_ROOT)/include/
endif

.PHONY: all clean archive tarball docs input prepare_sourcecode install uninstall


####################################################
# Parallel make execution (make -j) requires a bit of ordering
molgw:
	+make  prepare_sourcecode
	+make  level000
	+make  level005
	+make  level008
	+make  level010
	+make  level015
	+make  level020
	+make  level022
	+make  level025
	+make  level030
	+make  level040
	+make  level050
	+make  level060
	+make  level065
	+make  level070
	+make  level100
	+make  ../molgw

spectrum:
	+make  prepare_sourcecode
	+make  level000
	+make  ../spectrum

all: molgw spectrum

level000: m_definitions.o libint_onebody.o libint_onebody_grad.o libint_twobody.o libint_twobody_grad.o \
          lebedev_quadrature.o pdbssolver1.o boys_function.o

level005: m_warning.o real_spherical_harmonics.o

level008: m_mpi_tools.o

level010: m_tddft_variables.o m_elements.o m_timing.o m_memory.o \
          m_string_tools.o m_numerical_tools.o m_linear_algebra.o  \
	  m_mpi.o m_lbfgs.o m_cart_to_pure.o m_libxc_tools.o

level015: m_scalapack.o m_atoms.o m_ecp.o m_gaussian.o

level020: m_inputparam.o m_basis_set.o

level022: m_libint_tools.o m_block_diago.o m_hamiltonian_tools.o

level025: m_eri.o m_scf.o m_dft_grid.o

level030: m_eri_calculate.o m_density_tools.o m_io.o

level040: m_spectral_function.o m_hamiltonian_twobodies.o m_hamiltonian_onebody.o

level050: m_hamiltonian_wrapper.o m_eri_ao_mo.o m_restart.o m_multipole.o

level060: m_selfenergy_tools.o m_virtual_orbital_space.o m_tddft_fxc.o

level065: m_pt_density_matrix.o m_fourier_quadrature.o m_gw_selfenergy_grid.o 

level070: m_scf_loop.o m_ci.o m_dm_analysis.o m_selfenergy_evaluation.o m_dm_mbpt.o m_tddft_propagator.o 

level100: pt_selfenergy.o linear_response.o static_polarizability.o gw_selfenergy.o \
          gwgamma_selfenergy.o bse.o \
          mp2.o force.o molgw.o



####################################################


../spectrum: m_definitions.o m_warning.o spectrum.o
	$(LINKER) m_definitions.o m_warning.o spectrum.o -o ../spectrum $(FFTW) $(FFTW_INCS) $(LINKS) $(LAPACK)

../molgw: m_definitions.o libint_onebody.o libint_onebody_grad.o libint_twobody.o libint_twobody_grad.o \
        lebedev_quadrature.o pdbssolver1.o boys_function.o \
        m_warning.o real_spherical_harmonics.o m_lbfgs.o \
        m_mpi_tools.o m_libxc_tools.o \
        m_elements.o m_tddft_variables.o m_timing.o m_memory.o m_string_tools.o m_numerical_tools.o m_linear_algebra.o m_mpi.o m_cart_to_pure.o \
        m_scalapack.o m_atoms.o m_ecp.o m_gaussian.o m_block_diago.o \
        m_basis_set.o \
        m_inputparam.o \
        m_eri.o m_scf.o m_dft_grid.o m_libint_tools.o m_selfenergy_tools.o \
        m_eri_calculate.o m_eri_ao_mo.o m_density_tools.o m_tddft_fxc.o m_dm_mbpt.o \
        m_spectral_function.o m_hamiltonian_tools.o m_hamiltonian_twobodies.o m_hamiltonian_onebody.o m_ci.o  \
        m_hamiltonian_wrapper.o \
        m_virtual_orbital_space.o pt_selfenergy.o m_pt_density_matrix.o m_fourier_quadrature.o \
        linear_response.o static_polarizability.o gw_selfenergy.o m_gw_selfenergy_grid.o gwgamma_selfenergy.o \
        bse.o m_multipole.o m_io.o m_restart.o mp2.o m_scf_loop.o m_dm_analysis.o m_selfenergy_evaluation.o force.o m_tddft_propagator.o\
        molgw.o
	$(LINKER) \
                        m_definitions.o m_warning.o m_lbfgs.o m_elements.o m_cart_to_pure.o \
                        m_mpi.o m_mpi_tools.o \
                        m_libxc_tools.o m_scalapack.o m_tddft_variables.o m_timing.o m_memory.o m_inputparam.o \
                        m_string_tools.o m_numerical_tools.o m_linear_algebra.o m_scf.o m_selfenergy_tools.o m_atoms.o m_ecp.o \
                        m_gaussian.o m_block_diago.o m_basis_set.o m_eri.o m_dft_grid.o m_libint_tools.o \
                        m_eri_calculate.o m_eri_ao_mo.o m_density_tools.o m_tddft_fxc.o m_dm_mbpt.o \
                        m_spectral_function.o m_ci.o \
                        libint_onebody.o libint_onebody_grad.o libint_twobody.o libint_twobody_grad.o \
                        boys_function.o lebedev_quadrature.o m_virtual_orbital_space.o pt_selfenergy.o m_pt_density_matrix.o \
                        mp2.o gw_selfenergy.o m_gw_selfenergy_grid.o gwgamma_selfenergy.o bse.o \
                        linear_response.o static_polarizability.o m_multipole.o m_io.o m_restart.o \
                        m_hamiltonian_tools.o m_hamiltonian_twobodies.o m_hamiltonian_onebody.o m_scf_loop.o m_dm_analysis.o \
                        m_hamiltonian_wrapper.o m_fourier_quadrature.o \
                        pdbssolver1.o m_selfenergy_evaluation.o real_spherical_harmonics.o force.o m_tddft_propagator.o \
                        molgw.o \
                        -o ../molgw \
                        $(SCALAPACK) $(LAPACK)   \
                        $(LIBXC) $(LIBXC_INCS)  \
                        $(LIBINT) $(LIBINT_INCS) \
                        $(ELPA) $(ELPA_INCS) \
                        $(LINKS)

molgw.o: m_definitions.f90 m_warning.f90 m_tddft_variables.f90 m_timing.f90 m_memory.f90 m_inputparam.f90 m_lbfgs.f90 \
         m_atoms.f90 m_scf.f90 m_scalapack.f90  m_ecp.f90 m_io.f90 \
         m_gaussian.f90 m_basis_set.f90  m_eri.f90 m_dft_grid.f90 m_spectral_function.f90 m_ci.f90 m_eri_calculate.f90 m_eri_ao_mo.f90 \
         m_hamiltonian_onebody.f90 m_selfenergy_tools.f90 m_scf_loop.f90 \
         m_dm_analysis.f90 m_tddft_propagator.f90 m_restart.f90 m_dm_mbpt.f90 m_multipole.f90 \
         m_virtual_orbital_space.f90 m_selfenergy_evaluation.f90 \
         molgw.f90
	$(FC_W_OPTS) -c molgw.f90 $(LIBXC_INCS)

m_fourier_quadrature.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_atoms.f90 m_basis_set.f90 m_inputparam.f90 m_hamiltonian_onebody.f90 \
        m_fourier_quadrature.f90
	$(FC_W_OPTS) -c m_fourier_quadrature.f90  $(LIBXC_INCS)

force.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_atoms.f90 m_basis_set.f90 m_inputparam.f90 m_hamiltonian_onebody.f90 m_eri.f90 m_eri_calculate.f90 \
         force.f90
	$(FC_W_OPTS) -c force.f90  $(LIBXC_INCS)

m_selfenergy_evaluation.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_inputparam.f90 m_memory.f90 m_scf.f90 \
                         m_eri.f90 m_dft_grid.f90 m_spectral_function.f90 m_eri_calculate.f90 m_eri_ao_mo.f90 m_io.f90 \
                         m_selfenergy_tools.f90 m_virtual_orbital_space.f90 m_gw_selfenergy_grid.f90 m_selfenergy_evaluation.f90
	$(FC_W_OPTS) -c m_selfenergy_evaluation.f90 $(LIBXC_INCS)

m_scf_loop.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_memory.f90 m_atoms.f90 m_inputparam.f90 m_scf.f90 m_scalapack.f90 m_mpi.f90 \
              m_basis_set.f90 m_dft_grid.f90 m_spectral_function.f90 m_eri_ao_mo.f90 \
              m_hamiltonian_tools.f90 m_hamiltonian_wrapper.f90 m_hamiltonian_twobodies.f90 \
              m_selfenergy_tools.f90 m_linear_algebra.f90 m_restart.f90 \
              m_scf_loop.f90
	$(FC_W_OPTS) -c m_scf_loop.f90 $(LIBXC_INCS)

m_dm_mbpt.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_memory.f90 m_inputparam.f90 m_scf.f90 \
             m_spectral_function.f90 m_atoms.f90 m_basis_set.f90 m_multipole.f90 m_hamiltonian_tools.f90 \
             m_hamiltonian_wrapper.f90 m_selfenergy_tools.f90 m_gw_selfenergy_grid.f90 m_dm_mbpt.f90
	$(FC_W_OPTS) -c m_dm_mbpt.f90 $(LIBXC_INCS)

m_dm_analysis.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_memory.f90 m_atoms.f90 m_inputparam.f90 m_scf.f90 m_scalapack.f90 \
              m_gaussian.f90 m_lbfgs.f90 m_ecp.f90 m_density_tools.f90 m_multipole.f90 m_io.f90 \
              m_basis_set.f90  m_eri.f90 m_dft_grid.f90 m_spectral_function.f90 m_eri_calculate.f90 m_eri_ao_mo.f90 m_hamiltonian_tools.f90 \
              m_hamiltonian_onebody.f90 m_selfenergy_tools.f90 m_virtual_orbital_space.f90 \
              m_dm_analysis.f90
	$(FC_W_OPTS) -c m_dm_analysis.f90 $(LIBXC_INCS)

gwgamma_selfenergy.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 m_spectral_function.f90 m_eri_ao_mo.f90 \
      m_inputparam.f90 m_tddft_fxc.f90 m_selfenergy_tools.f90 gwgamma_selfenergy.f90
	$(FC_W_OPTS) -c gwgamma_selfenergy.f90 $(LIBXC_INCS)

gw_selfenergy.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_scalapack.f90 \
                 m_basis_set.f90 m_spectral_function.f90 m_eri_ao_mo.f90 m_linear_algebra.f90 \
                 m_inputparam.f90 m_selfenergy_tools.f90 gw_selfenergy.f90
	$(FC_W_OPTS) -c gw_selfenergy.f90 $(LIBXC_INCS)

m_gw_selfenergy_grid.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_scalapack.f90 \
      m_basis_set.f90 m_spectral_function.f90 m_eri_ao_mo.f90 m_linear_algebra.f90 \
      m_inputparam.f90 m_selfenergy_tools.f90 m_gw_selfenergy_grid.f90
	$(FC_W_OPTS) -c m_gw_selfenergy_grid.f90 $(LIBXC_INCS)

bse.o: m_definitions.f90 m_mpi.f90 m_scalapack.f90 m_timing.f90 m_warning.f90 m_memory.f90 \
       m_basis_set.f90 m_spectral_function.f90 m_inputparam.f90 m_tddft_fxc.f90 m_eri_ao_mo.f90 m_dft_grid.f90 \
       bse.f90
	$(FC_W_OPTS) -c bse.f90 $(LIBXC_INCS)

static_polarizability.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 \
	m_basis_set.f90 m_spectral_function.f90 m_eri_ao_mo.f90 m_inputparam.f90 \
        static_polarizability.f90
	$(FC_W_OPTS) -c static_polarizability.f90 $(LIBXC_INCS)

linear_response.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_cart_to_pure.f90 m_scalapack.f90 m_io.f90 \
                   m_inputparam.f90 m_block_diago.f90 m_basis_set.f90 m_spectral_function.f90 m_eri_ao_mo.f90 m_dft_grid.f90 m_hamiltonian_onebody.f90 \
                   linear_response.f90
	$(FC_W_OPTS) -c linear_response.f90 $(LIBXC_INCS)

m_hamiltonian_onebody.o: m_definitions.f90 m_mpi.f90 m_scalapack.f90 m_warning.f90 m_memory.f90 m_timing.f90 m_atoms.f90 m_cart_to_pure.f90 m_ecp.f90 m_dft_grid.f90 \
                        m_basis_set.f90 m_libint_tools.f90 m_inputparam.f90 m_tddft_variables.f90 m_io.f90 m_hamiltonian_onebody.f90
	$(FC_W_OPTS) -c m_hamiltonian_onebody.f90 $(LIBXC_INCS)

m_hamiltonian_tools.o: m_definitions.f90 m_mpi.f90 m_scalapack.f90 m_warning.f90 m_timing.f90 m_memory.f90 m_cart_to_pure.f90 \
	               m_linear_algebra.f90 m_inputparam.f90 m_basis_set.f90 m_hamiltonian_tools.f90
	$(FC_W_OPTS) -c m_hamiltonian_tools.f90 $(LIBXC_INCS)

m_hamiltonian_twobodies.o: m_definitions.f90 m_mpi.f90 m_scalapack.f90 m_warning.f90 m_timing.f90 m_memory.f90 m_cart_to_pure.f90 m_libxc_tools.f90 \
	         m_inputparam.f90 m_eri_calculate.f90 m_basis_set.f90 m_dft_grid.f90 m_density_tools.f90 m_io.f90 m_hamiltonian_twobodies.f90
	$(FC_W_OPTS) -c m_hamiltonian_twobodies.f90 $(LIBXC_INCS)

m_hamiltonian_wrapper.o: m_definitions.f90 m_mpi.f90 m_scalapack.f90 m_warning.f90 m_timing.f90 m_basis_set.f90 m_scf.f90 \
	         m_inputparam.f90 m_hamiltonian_twobodies.f90 m_hamiltonian_tools.f90 m_hamiltonian_wrapper.f90
	$(FC_W_OPTS) -c m_hamiltonian_wrapper.f90 $(LIBXC_INCS)

m_density_tools.o: m_definitions.f90 m_atoms.f90 m_mpi.f90 m_gaussian.f90 m_inputparam.f90 m_basis_set.f90 m_hamiltonian_tools.f90 m_dft_grid.f90 \
                   m_density_tools.f90
	$(FC_W_OPTS) -c m_density_tools.f90 $(LIBXC_INCS)

m_multipole.o: m_definitions.f90 m_atoms.f90 m_inputparam.f90 m_basis_set.f90 m_hamiltonian_tools.f90 m_hamiltonian_onebody.f90 m_multipole.f90
	$(FC_W_OPTS) -c m_multipole.f90 $(LIBXC_INCS)

m_io.o: m_definitions.f90 m_mpi.f90 m_warning.f90 m_atoms.f90 m_memory.f90 m_cart_to_pure.f90 \
	      m_basis_set.f90 m_libint_tools.f90 m_inputparam.f90 m_eri.f90 \
        m_dft_grid.f90 git_sha.f90 m_tddft_variables.f90 m_timing.f90 m_linear_algebra.f90 m_string_tools.f90 \
			  m_libxc_tools.f90 m_hamiltonian_tools.f90 m_io.f90
	$(FC_W_OPTS) -c m_io.f90 $(LIBXC_INCS)

m_restart.o: m_definitions.f90 m_timing.f90 m_mpi.f90 m_atoms.f90 m_basis_set.f90 m_hamiltonian_onebody.f90 m_inputparam.f90 m_linear_algebra.f90 m_hamiltonian_tools.f90 \
             m_restart.f90
	$(FC_W_OPTS) -c m_restart.f90 $(LIBXC_INCS)

mp2.o: m_definitions.f90 m_mpi.f90 m_cart_to_pure.f90 \
       m_basis_set.f90 m_inputparam.f90 m_eri_ao_mo.f90 \
       mp2.f90
	$(FC_W_OPTS) -c mp2.f90 $(LIBXC_INCS)

m_virtual_orbital_space.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_scalapack.f90 \
	                  m_basis_set.f90 m_eri_ao_mo.f90 m_hamiltonian_tools.f90 m_linear_algebra.f90 m_io.f90 \
                          m_inputparam.f90 m_hamiltonian_onebody.f90 m_virtual_orbital_space.f90
	$(FC_W_OPTS) -c m_virtual_orbital_space.f90 $(LIBXC_INCS)

pt_selfenergy.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_basis_set.f90 \
                  m_spectral_function.f90 m_inputparam.f90 m_selfenergy_tools.f90 m_eri_ao_mo.f90 pt_selfenergy.f90
	$(FC_W_OPTS) -c pt_selfenergy.f90 $(LIBXC_INCS)

m_pt_density_matrix.o: m_definitions.f90 m_mpi.f90 m_memory.f90 m_timing.f90 m_warning.f90 m_scalapack.f90 m_basis_set.f90 \
                  m_spectral_function.f90 m_inputparam.f90 m_selfenergy_tools.f90 m_eri_ao_mo.f90 m_hamiltonian_tools.f90 m_hamiltonian_onebody.f90 m_io.f90 m_pt_density_matrix.f90
	$(FC_W_OPTS) -c m_pt_density_matrix.f90 $(LIBXC_INCS)

pdbssolver1.o: pdbssolver1.f
	$(FC_W_OPTS) -c pdbssolver1.f $(LIBXC_INCS) $(ELPA_INCS)

boys_function.o: boys_function.f90
	$(FC_W_OPTS) -c boys_function.f90 $(LIBXC_INCS)

real_spherical_harmonics.o: m_definitions.f90 real_spherical_harmonics.f90
	$(FC_W_OPTS) -c real_spherical_harmonics.f90 $(LIBXC_INCS)

lebedev_quadrature.o: lebedev_quadrature.f
	$(FC_W_OPTS) -c lebedev_quadrature.f $(LIBXC_INCS)

libint_onebody.o: libint_molgw.h libint_onebody.cc
	$(CXX_W_OPTS) -c libint_onebody.cc $(LIBINT_INCS)

libint_onebody_grad.o: libint_molgw.h libint_onebody_grad.cc
	$(CXX_W_OPTS) -c libint_onebody_grad.cc $(LIBINT_INCS)

libint_twobody.o: libint_molgw.h libint_twobody.cc
	$(CXX_W_OPTS) -c libint_twobody.cc $(LIBINT_INCS)

libint_twobody_grad.o: libint_molgw.h libint_twobody_grad.cc
	$(CXX_W_OPTS) -c libint_twobody_grad.cc $(LIBINT_INCS)

m_libint_tools.o: m_definitions.f90 m_cart_to_pure.f90 m_basis_set.f90 m_libint_tools.f90
	$(FC_W_OPTS) -c m_libint_tools.f90 $(LIBXC_INCS)

m_dft_grid.o: m_definitions.f90 m_memory.f90 m_inputparam.f90 m_warning.f90 m_cart_to_pure.f90 m_timing.f90 m_elements.f90 m_atoms.f90 m_mpi.f90 m_basis_set.f90 m_numerical_tools.f90 \
              m_dft_grid.f90
	$(FC_W_OPTS) -c m_dft_grid.f90 $(LIBXC_INCS)

m_ci.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_basis_set.f90 m_eri_ao_mo.f90 m_inputparam.f90 m_scalapack.f90 m_selfenergy_tools.f90 \
        m_linear_algebra.f90 m_hamiltonian_onebody.f90 m_ci.f90
	$(FC_W_OPTS) -c m_ci.f90 $(LIBXC_INCS)

m_spectral_function.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_scalapack.f90 m_eri.f90 m_eri_calculate.f90 m_inputparam.f90 m_numerical_tools.f90 \
                       m_spectral_function.f90
	$(FC_W_OPTS) -c m_spectral_function.f90 $(LIBXC_INCS)

m_eri_ao_mo.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_basis_set.f90 m_inputparam.f90 \
               m_scalapack.f90 m_eri.f90 m_eri_ao_mo.f90
	$(FC_W_OPTS) -c m_eri_ao_mo.f90 $(LIBXC_INCS)

m_tddft_fxc.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_basis_set.f90 m_inputparam.f90 m_dft_grid.f90 m_density_tools.f90 m_libxc_tools.f90 \
                m_tddft_fxc.f90
	$(FC_W_OPTS) -c m_tddft_fxc.f90 $(LIBXC_INCS)

m_eri_calculate.o: m_definitions.f90 m_mpi.f90 m_scalapack.f90 m_timing.f90 m_warning.f90 m_memory.f90 m_cart_to_pure.f90 m_basis_set.f90 m_inputparam.f90 \
                   m_libint_tools.f90 m_eri.f90 m_eri_calculate.f90
	$(FC_W_OPTS) -c m_eri_calculate.f90 $(LIBXC_INCS)

m_eri.o: m_definitions.f90 m_mpi.f90 m_timing.f90 m_memory.f90 m_warning.f90 m_scalapack.f90 m_cart_to_pure.f90 m_basis_set.f90 m_inputparam.f90 m_libint_tools.f90 m_eri.f90
	$(FC_W_OPTS) -c m_eri.f90 $(LIBXC_INCS)

m_basis_set.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_elements.f90 m_gaussian.f90 m_mpi.f90 m_ecp.f90 m_cart_to_pure.f90 m_atoms.f90 m_string_tools.f90 \
               m_basis_set.f90
	$(FC_W_OPTS) -c m_basis_set.f90 $(LIBXC_INCS)

m_gaussian.o: m_definitions.f90 m_mpi.f90 m_cart_to_pure.f90 m_string_tools.f90 m_gaussian.f90
	$(FC_W_OPTS) -c m_gaussian.f90 $(LIBXC_INCS)

m_scf.o: m_definitions.f90 m_mpi.f90 m_warning.f90 m_timing.f90 m_memory.f90 m_inputparam.f90 m_scalapack.f90 m_lbfgs.f90 m_linear_algebra.f90 \
	m_scf.f90
	$(FC_W_OPTS) -c m_scf.f90 $(LIBXC_INCS)

m_selfenergy_tools.o: m_definitions.f90 m_warning.f90 m_inputparam.f90 m_mpi.f90 m_atoms.f90 m_lbfgs.f90 m_basis_set.f90 m_dft_grid.f90 m_numerical_tools.f90 \
	              m_hamiltonian_tools.f90 m_hamiltonian_wrapper.f90 m_selfenergy_tools.f90
	$(FC_W_OPTS) -c m_selfenergy_tools.f90 $(LIBXC_INCS)

m_inputparam.o: m_definitions.f90 m_mpi.f90  m_elements.f90 m_warning.f90 m_atoms.f90 m_ecp.f90 \
	        input_variables.f90 input_variable_declaration.f90 basis_path.f90 echo_input_variables.f90 echo_input_variables_yaml.f90  m_string_tools.f90 m_libxc_tools.f90 \
                m_inputparam.f90
	$(FC_W_OPTS) -c m_inputparam.f90 $(LIBXC_INCS)

m_block_diago.o: m_definitions.f90 m_warning.f90 m_timing.f90 m_memory.f90 m_mpi.f90 m_scalapack.f90 m_inputparam.f90 m_linear_algebra.f90 m_block_diago.f90
	$(FC_W_OPTS) -c m_block_diago.f90 $(LIBXC_INCS) $(ELPA_INCS)

m_atoms.o: m_definitions.f90 m_warning.f90 m_elements.f90 m_linear_algebra.f90 m_lbfgs.f90 m_atoms.f90
	$(FC_W_OPTS) -c m_atoms.f90 $(LIBXC_INCS)

m_ecp.o: m_definitions.f90 m_warning.f90 m_elements.f90 m_string_tools.f90 m_ecp.f90
	$(FC_W_OPTS) -c m_ecp.f90 $(LIBXC_INCS)

m_scalapack.o: m_definitions.f90 m_warning.f90 m_linear_algebra.f90 m_mpi.f90 m_scalapack.f90
	$(FC_W_OPTS) -c m_scalapack.f90 $(LIBXC_INCS)

m_mpi.o: m_definitions.f90 m_warning.f90 m_mpi_tools.f90  \
         m_mpi.f90
	$(FC_W_OPTS) -c m_mpi.f90 $(LIBXC_INCS)

m_mpi_tools.o: m_definitions.f90 m_warning.f90  m_mpi_tools.f90
	$(FC_W_OPTS) -c m_mpi_tools.f90 $(LIBXC_INCS)

m_string_tools.o: m_definitions.f90 m_warning.f90 m_string_tools.f90
	$(FC_W_OPTS) -c m_string_tools.f90 $(LIBXC_INCS)

m_numerical_tools.o: m_definitions.f90 m_warning.f90 m_numerical_tools.f90
	$(FC_W_OPTS) -c m_numerical_tools.f90 $(LIBXC_INCS)

m_linear_algebra.o: m_definitions.f90 m_warning.f90 m_linear_algebra.f90
	$(FC_W_OPTS) -c m_linear_algebra.f90 $(LIBXC_INCS)

m_memory.o: m_definitions.f90 m_warning.f90 m_memory.f90
	$(FC_W_OPTS) -c m_memory.f90 $(LIBXC_INCS)

m_tddft_variables.o: m_tddft_variables.f90
	$(FC_W_OPTS) -c m_tddft_variables.f90 $(LIBXC_INCS)

m_timing.o: m_definitions.f90 m_warning.f90 m_timing.f90
	$(FC_W_OPTS) -c m_timing.f90 $(LIBXC_INCS)

m_elements.o: m_definitions.f90 m_warning.f90 m_elements.f90
	$(FC_W_OPTS) -c m_elements.f90 $(LIBXC_INCS)

m_cart_to_pure.o: m_definitions.f90 m_warning.f90 m_cart_to_pure.f90
	$(FC_W_OPTS) -c m_cart_to_pure.f90 $(LIBXC_INCS)

m_warning.o: m_definitions.f90 m_warning.f90
	$(FC_W_OPTS) -c m_warning.f90 $(LIBXC_INCS)

m_libxc_tools.o: m_definitions.f90 m_warning.f90 m_libxc_tools.f90
	$(FC_W_OPTS) -c m_libxc_tools.f90 $(LIBXC_INCS)

m_lbfgs.o: m_definitions.f90 m_lbfgs.f90
	$(FC_W_OPTS) -c m_lbfgs.f90 $(LIBXC_INCS)

m_definitions.o: m_definitions.f90
	$(FC_W_OPTS) -c m_definitions.f90 $(LIBXC_INCS)

m_tddft_propagator.o: m_definitions.f90 m_basis_set.f90 m_hamiltonian_wrapper.f90 m_memory.f90 m_tddft_variables.f90 m_linear_algebra.f90 \
                      m_hamiltonian_tools.f90 m_hamiltonian_onebody.f90 m_string_tools.f90 m_io.f90 m_multipole.f90 \
                      m_inputparam.f90 m_dft_grid.f90 m_scf.f90 m_warning.f90 m_timing.f90 m_density_tools.f90 m_tddft_propagator.f90
	 $(FC_W_OPTS) -c m_tddft_propagator.f90 $(LIBXC_INCS)

spectrum.o: m_definitions.f90 m_warning.f90 spectrum.f90
	$(FC_W_OPTS) -c spectrum.f90  $(FFTW_INCS)

tarball:
	+make archive

archive:
	+make prepare_sourcecode
	mkdir -p ../tmp/molgw/src
	cp -rp  ../README.md ../COPYING ../RELEASE.md  \
            ../basis/ ../utils/ ../tests/            \
            ../docs/ ../config/   \
            ../tmp/molgw
	cp -rp  ../src/*.py ../src/*.f90 ../src/*.f  \
            ../src/*.cc ../src/*.h  \
            ../src/my_machine.arch ../src/Makefile \
            ../tmp/molgw/src
	cd ../tmp && tar czf ../molgw.tgz molgw/
	rm -rf ../tmp/

prepare_sourcecode:
	$(PYTHON) prepare_sourcecode.py

input:
	$(PYTHON) ../utils/input_variables.py

docs:
	pdflatex ../docs/notes.tex
	pdflatex ../docs/notes.tex
	rm -f notes.aux notes.log notesNotes.bib
	mv notes.pdf ../docs/notes.pdf

clean:
	$(RM) *.o *.mod ../molgw.tgz ../molgw ../spectrum

install:
	install -d $(DESTDIR)$(bindir)
	install  ../molgw $(DESTDIR)$(bindir)
	install  ../utils/run_molgw.py $(DESTDIR)$(bindir)
	install  ../utils/scf_status.py $(DESTDIR)$(bindir)
	install -d $(DESTDIR)$(bindir)/basis/
	install  ../basis/* $(DESTDIR)$(bindir)/basis/

uninstall:
	$(RM) $(DESTDIR)$(bindir)/molgw
	$(RM) $(DESTDIR)$(bindir)/run_molgw.py
	$(RM) $(DESTDIR)$(bindir)/scf_status.py
	$(RM) -r $(DESTDIR)$(bindir)/basis/
