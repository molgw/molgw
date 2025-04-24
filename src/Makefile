# This file is part of MOLGW
# Author: Fabien Bruneval

ifeq ($(FC),f77)
  FC=gfortran
endif
-include ../my_machine.arch
-include ./my_machine.arch

PREFIX ?=
LDFLAGS ?=
PYTHON ?= python3
FCFLAGS ?= -O2 -cpp
CXXFLAGS ?=-O2
CPPFLAGS ?=
#ELPA=
#ELPA_INCS=



FC_W_OPTS= $(FC) $(FCFLAGS) $(OPENMP) $(CPPFLAGS)
CXX_W_OPTS= $(CXX) $(CXXFLAGS) $(CPPFLAGS)
LINKER= $(LIBTOOL) $(FC_W_OPTS) $(LDFLAGS)

LINKS= -lstdc++ -lm

##################################
# LAPACK /SCALAPACK
# if both LAPACK and SCALAPACK are empty, then set the default
#
##################################
ifeq ($(SCALAPACK),)
  ifeq ($(LAPACK),)
    LAPACK ?=-lblas -llapack
  endif
endif

##################################
#
# LIBXC variables
# 
#   if LIBXC is already set, don't do anything
#   if LIBXC is not set, check out LIBXC_ROOT and set LIBXC accordingly
#   if LIBXC_INCS is already set, don't do anything
#   if LIBXC_INCS is not set, check out LIBXC_ROOT and set LIBXC_INCS accordingly
#
##################################
ifeq ($(LIBXC),)
  ifeq ($(LIBXC_ROOT),)
    ifeq (,$(findstring -DNO_LIBXC,$(CPPFLAGS)))
      LIBXC=-lxc
    else
      LIBXC=
    endif
  else
    LIBXC=-L$(LIBXC_ROOT)/lib/ -L$(LIBXC_ROOT)/lib64/ -Wl,-rpath,$(LIBXC_ROOT)/lib/ -Wl,-rpath,$(LIBXC_ROOT)/lib64/ -lxc
  endif
endif
ifeq ($(LIBXC_INCS),)
  ifneq ($(LIBXC_ROOT),)
    LIBXC_INCS=-I$(LIBXC_ROOT)/include
  endif
endif

##################################
#
# LIBCINT variables
#
#   if LIBCINT is already set, don't do anything
#   if LIBCINT is not set, check out LIBCINT_ROOT and set LIBCINT accordingly
#
##################################
ifeq ($(LIBCINT),)
  ifeq ($(LIBCINT_ROOT),)
    LIBCINT= -lcint
  else
    LIBCINT= -L$(LIBCINT_ROOT)/lib64/ -Wl,-rpath,$(LIBCINT_ROOT)/lib64/ -lcint
  endif
endif

##################################
#
# LIBINT variables
#
#   if LIBINT is already set, don't do anything
#   if LIBINT is not set, check out LIBINT_ROOT and set LIBINT accordingly
#   if LIBINT_INCS is already set, don't do anything
#   if LIBINT_INCS is not set, check out LIBINT_ROOT and set LIBINT_INCS accordingly
#
##################################
ifeq ($(LIBINT),)
  ifeq ($(LIBINT_ROOT),)
    LIBINT=
  else
    LIBINT=-L$(LIBINT_ROOT)/lib -lint2
  endif
endif
ifeq ($(LIBINT_INCS),)
  ifeq ($(LIBINT_ROOT),)
    LIBINT_INCS=
  else
    LIBINT_INCS=-I$(LIBINT_ROOT)/include/libint2 -I$(LIBINT_ROOT)/include
  endif
endif

##################################
#
# HDF5 variables
#   if HDF5 is already set, don't do anything
#   if HDF5 is not set, check out HDF5_ROOT and set HDF5 accordingly
#   if HDF5_INCS is already set, don't do anything
#   if HDF5_INCS is not set, check out HDF5_ROOT and set HDF5_INCS accordingly
#
##################################
ifeq ($(HDF5),)
  ifeq ($(HDF5_ROOT),)
    HDF5=
  else
    HDF5=-L$(HDF5_ROOT)/lib/ -Wl,-rpath,$(HDF5_ROOT)/lib/ -lhdf5 -lhdf5_fortran 
  endif
endif
ifeq ($(HDF5_INCS),)
  ifeq ($(HDF5_ROOT),)
    HDF5_INCS=
  else
    HDF5_INCS=-I$(HDF5_ROOT)/include/
  endif
endif


####################################################
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
	+make  level035
	+make  level040
	+make  level045
	+make  level050
	+make  level060
	+make  level065
	+make  level068
	+make  level070
	+make  level100
	+make  ../molgw

all: molgw

level000: m_definitions.o libint_onebody.o libint_onebody_grad.o libint_twobody.o libint_twobody_grad.o \
          lebedev_quadrature.o pdbssolver1.o boys_function.o m_nofoutput.o

level005: m_warning.o real_spherical_harmonics.o m_lbfgs_intern.o m_integd.o m_rdmd.o m_anti2unit.o m_scalapack_interface.o

level008: m_mpi_tools.o m_elag.o

level010: m_elements.o m_timing.o m_hdf5_tools.o m_memory.o \
          m_string_tools.o m_numerical_tools.o m_linear_algebra.o  \
	  m_mpi.o m_lbfgs.o m_cart_to_pure.o m_libxc_tools.o m_hessian.o m_gammatodm2.o

level015: m_scalapack.o m_atoms.o m_ecp.o m_gaussian.o m_diagf.o m_e_grad_occ.o m_e_grad_occ_cpx.o  

level020: m_inputparam.o m_basis_set.o

level022: m_libint_tools.o m_libcint_tools.o m_block_diago.o m_hamiltonian_tools.o

level025: m_eri.o m_scf.o m_dft_grid.o m_optocc.o m_optorb.o m_tz_pccd_amplitudes.o 

level030: m_eri_calculate.o m_density_tools.o m_io.o

level035: m_eri_ao_mo.o

level040: m_hamiltonian_twobodies.o m_hamiltonian_onebody.o m_tddft_fxc.o

level045: m_spectral_function.o

level050: m_hamiltonian_wrapper.o m_restart.o m_multipole.o m_spectra.o m_build_bse.o

level060: m_selfenergy_tools.o m_virtual_orbital_space.o m_linear_response.o

level065: m_pt_density_matrix.o m_fourier_quadrature.o m_gw_selfenergy_grid.o m_noft_driver.o m_g3w2_selfenergy.o m_tdhf_selfenergy.o m_relativistic.o \
          m_gw_selfenergy_analytic.o

level068: m_pt_selfenergy.o m_force.o

level070: m_scf_loop.o m_ci.o m_dm_analysis.o m_selfenergy_evaluation.o m_acfd.o m_dm_mbpt.o m_realtime_tddft.o m_noft.o m_mp2_energy.o

level100: molgw.o



####################################################


../molgw: m_definitions.o m_hdf5_tools.o libint_onebody.o libint_onebody_grad.o libint_twobody.o libint_twobody_grad.o \
        lebedev_quadrature.o pdbssolver1.o boys_function.o \
        m_warning.o real_spherical_harmonics.o m_lbfgs.o \
        m_mpi_tools.o m_libxc_tools.o \
        m_elements.o m_timing.o m_memory.o m_string_tools.o m_numerical_tools.o m_linear_algebra.o m_mpi.o m_cart_to_pure.o \
        m_scalapack.o m_atoms.o m_ecp.o m_gaussian.o m_block_diago.o \
        m_basis_set.o \
        m_inputparam.o \
        m_eri.o m_scf.o m_dft_grid.o m_libint_tools.o m_libcint_tools.o m_selfenergy_tools.o \
        m_eri_calculate.o m_eri_ao_mo.o m_density_tools.o m_tddft_fxc.o m_dm_mbpt.o m_spectra.o \
        m_spectral_function.o m_hamiltonian_tools.o m_hamiltonian_twobodies.o m_hamiltonian_onebody.o m_ci.o  \
        m_hamiltonian_wrapper.o \
        m_virtual_orbital_space.o m_pt_selfenergy.o m_pt_density_matrix.o m_fourier_quadrature.o \
        m_linear_response.o m_acfd.o m_gw_selfenergy_analytic.o m_gw_selfenergy_grid.o m_g3w2_selfenergy.o m_tdhf_selfenergy.o \
        m_nofoutput.o m_lbfgs_intern.o m_integd.o m_rdmd.o m_anti2unit.o m_elag.o m_hessian.o m_gammatodm2.o m_diagf.o \
        m_e_grad_occ.o m_e_grad_occ_cpx.o m_optocc.o m_optorb.o m_tz_pccd_amplitudes.o m_noft_driver.o gitver.o \
        m_build_bse.o m_multipole.o m_io.o m_restart.o m_mp2_energy.o m_relativistic.o m_noft.o m_scf_loop.o m_dm_analysis.o m_selfenergy_evaluation.o \
        m_force.o m_realtime_tddft.o\
        molgw.o
	$(LINKER) \
                        m_definitions.o m_warning.o m_lbfgs.o m_elements.o m_cart_to_pure.o \
                        m_mpi.o m_mpi_tools.o m_hdf5_tools.o m_scalapack_interface.o \
                        m_libxc_tools.o m_scalapack.o m_timing.o m_memory.o m_inputparam.o \
                        m_string_tools.o m_numerical_tools.o m_linear_algebra.o m_scf.o m_selfenergy_tools.o m_atoms.o m_ecp.o \
                        m_gaussian.o m_block_diago.o m_basis_set.o m_eri.o m_dft_grid.o m_libint_tools.o m_libcint_tools.o \
                        m_eri_calculate.o m_eri_ao_mo.o m_density_tools.o m_tddft_fxc.o m_dm_mbpt.o m_spectra.o \
                        m_spectral_function.o m_ci.o \
                        libint_onebody.o libint_onebody_grad.o libint_twobody.o libint_twobody_grad.o \
                        boys_function.o lebedev_quadrature.o m_virtual_orbital_space.o m_pt_selfenergy.o m_pt_density_matrix.o \
                        m_nofoutput.o m_lbfgs_intern.o m_integd.o m_rdmd.o m_anti2unit.o m_elag.o m_hessian.o m_gammatodm2.o m_diagf.o \
                        m_e_grad_occ.o m_e_grad_occ_cpx.o m_optocc.o m_optorb.o m_tz_pccd_amplitudes.o m_noft_driver.o gitver.o \
                        m_mp2_energy.o m_relativistic.o m_noft.o m_gw_selfenergy_analytic.o m_gw_selfenergy_grid.o m_g3w2_selfenergy.o m_tdhf_selfenergy.o m_build_bse.o \
                        m_linear_response.o m_acfd.o m_multipole.o m_io.o m_restart.o \
                        m_hamiltonian_tools.o m_hamiltonian_twobodies.o m_hamiltonian_onebody.o m_scf_loop.o m_dm_analysis.o \
                        m_hamiltonian_wrapper.o m_fourier_quadrature.o \
                        pdbssolver1.o m_selfenergy_evaluation.o real_spherical_harmonics.o m_force.o m_realtime_tddft.o \
                        molgw.o \
                        -o ../molgw \
                        $(SCALAPACK) $(LAPACK)   \
                        $(LIBXC) $(LIBXC_INCS)  \
                        $(LIBINT) $(LIBINT_INCS) \
                        $(ELPA) $(ELPA_INCS) \
                        $(LIBCINT) \
			$(HDF5) $(HDF5_INCS) \
                        $(LINKS)

molgw.o: m_definitions.o m_warning.o m_hdf5_tools.o m_timing.o m_memory.o m_inputparam.o m_lbfgs.o \
         m_atoms.o m_scf.o m_scalapack.o  m_ecp.o m_io.o \
         m_gaussian.o m_basis_set.o  m_eri.o m_dft_grid.o m_spectral_function.o m_ci.o m_eri_calculate.o m_eri_ao_mo.o \
         m_hamiltonian_onebody.o m_selfenergy_tools.o m_scf_loop.o \
         m_dm_analysis.o m_realtime_tddft.o m_restart.o m_dm_mbpt.o m_multipole.o \
         m_virtual_orbital_space.o m_selfenergy_evaluation.o m_noft.o m_relativistic.o \
         m_fourier_quadrature.o m_libcint_tools.o m_linear_response.o m_acfd.o m_hamiltonian_twobodies.o m_mp2_energy.o m_force.o \
         molgw.f90
	$(FC_W_OPTS) -c molgw.f90 $(LIBXC_INCS)

m_fourier_quadrature.o: m_definitions.o m_warning.o m_timing.o m_atoms.o m_basis_set.o m_inputparam.o m_hamiltonian_onebody.o \
        m_fourier_quadrature.f90
	$(FC_W_OPTS) -c m_fourier_quadrature.f90  $(LIBXC_INCS)

m_force.o: m_definitions.o m_warning.o m_timing.o m_atoms.o m_basis_set.o m_inputparam.o m_hamiltonian_tools.o \
         m_hamiltonian_onebody.o m_eri.o m_eri_calculate.o \
         m_force.f90
	$(FC_W_OPTS) -c m_force.f90  $(LIBXC_INCS)

m_selfenergy_evaluation.o: m_definitions.o m_warning.o m_timing.o m_inputparam.o m_memory.o m_scf.o \
                           m_eri.o m_spectral_function.o m_eri_calculate.o m_eri_ao_mo.o m_io.o \
                           m_selfenergy_tools.o m_virtual_orbital_space.o m_gw_selfenergy_grid.o m_linear_response.o \
			   m_g3w2_selfenergy.o m_tdhf_selfenergy.o m_gw_selfenergy_analytic.o m_pt_selfenergy.o \
                           m_selfenergy_evaluation.f90
	$(FC_W_OPTS) -c m_selfenergy_evaluation.f90 $(LIBXC_INCS)

m_scf_loop.o: m_definitions.o m_warning.o m_timing.o m_memory.o m_atoms.o m_inputparam.o m_scf.o m_scalapack.o m_mpi.o \
              m_basis_set.o m_dft_grid.o m_spectral_function.o m_eri_ao_mo.o \
              m_hamiltonian_tools.o m_hamiltonian_wrapper.o m_hamiltonian_twobodies.o \
              m_selfenergy_tools.o m_linear_algebra.o m_restart.o m_linear_response.o m_gw_selfenergy_analytic.o m_pt_selfenergy.o \
              m_scf_loop.f90
	$(FC_W_OPTS) -c m_scf_loop.f90 $(LIBXC_INCS)

m_dm_mbpt.o: m_definitions.o m_warning.o m_timing.o m_memory.o m_inputparam.o m_scf.o \
             m_spectral_function.o m_atoms.o m_basis_set.o m_multipole.o m_hamiltonian_tools.o \
             m_hamiltonian_wrapper.o m_selfenergy_tools.o m_gw_selfenergy_grid.o m_pt_density_matrix.o \
             m_linear_response.o m_dm_mbpt.f90
	$(FC_W_OPTS) -c m_dm_mbpt.f90 $(LIBXC_INCS)

m_dm_analysis.o: m_definitions.o m_warning.o m_timing.o m_memory.o m_atoms.o m_inputparam.o m_scf.o m_scalapack.o \
              m_gaussian.o m_lbfgs.o m_ecp.o m_density_tools.o m_multipole.o m_io.o \
              m_basis_set.o  m_eri.o m_dft_grid.o m_spectral_function.o m_eri_calculate.o m_eri_ao_mo.o m_hamiltonian_tools.o \
              m_hamiltonian_onebody.o m_selfenergy_tools.o m_virtual_orbital_space.o \
              m_dm_analysis.f90
	$(FC_W_OPTS) -c m_dm_analysis.f90 $(LIBXC_INCS)

m_g3w2_selfenergy.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_basis_set.o m_spectral_function.o m_eri_ao_mo.o \
      m_inputparam.o m_tddft_fxc.o m_selfenergy_tools.o m_linear_response.o m_g3w2_selfenergy.f90
	$(FC_W_OPTS) -c m_g3w2_selfenergy.f90 $(LIBXC_INCS)

m_tdhf_selfenergy.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_basis_set.o m_spectral_function.o m_eri_ao_mo.o \
      m_inputparam.o m_tddft_fxc.o m_selfenergy_tools.o m_linear_response.o m_tdhf_selfenergy.f90
	$(FC_W_OPTS) -c m_tdhf_selfenergy.f90 $(LIBXC_INCS)

m_gw_selfenergy_analytic.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_scalapack.o \
                 m_spectral_function.o m_eri_ao_mo.o m_linear_algebra.o \
                 m_inputparam.o m_selfenergy_tools.o m_gw_selfenergy_analytic.f90
	$(FC_W_OPTS) -c m_gw_selfenergy_analytic.f90 $(LIBXC_INCS)

m_gw_selfenergy_grid.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_scalapack.o \
      m_basis_set.o m_spectral_function.o m_eri_ao_mo.o m_linear_algebra.o \
      m_inputparam.o m_selfenergy_tools.o m_linear_response.o m_gw_selfenergy_grid.f90
	$(FC_W_OPTS) -c m_gw_selfenergy_grid.f90 $(LIBXC_INCS)

m_build_bse.o: m_definitions.o m_mpi.o m_scalapack.o m_timing.o m_warning.o m_memory.o \
       m_basis_set.o m_spectral_function.o m_inputparam.o m_tddft_fxc.o m_eri_ao_mo.o m_dft_grid.o \
       m_build_bse.f90
	$(FC_W_OPTS) -c m_build_bse.f90 $(LIBXC_INCS)

m_linear_response.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_cart_to_pure.o m_scalapack.o \
                   m_inputparam.o m_block_diago.o m_basis_set.o m_spectral_function.o m_eri_ao_mo.o \
                   m_spectra.o m_build_bse.o m_linear_response.f90
	$(FC_W_OPTS) -c m_linear_response.f90 $(LIBXC_INCS)

m_acfd.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_cart_to_pure.o m_scalapack.o m_numerical_tools.o \
          m_inputparam.o m_basis_set.o m_spectral_function.o m_gw_selfenergy_grid.o m_selfenergy_tools.o \
          m_spectra.o m_scf.o m_linear_response.o m_dft_grid.o m_hamiltonian_twobodies.o \
	  m_acfd.f90
	$(FC_W_OPTS) -c m_acfd.f90 $(LIBXC_INCS)

m_hamiltonian_onebody.o: m_definitions.o m_mpi.o m_scalapack.o m_warning.o m_memory.o m_timing.o m_atoms.o m_cart_to_pure.o m_ecp.o m_dft_grid.o \
                        m_basis_set.o m_libint_tools.o m_libcint_tools.o m_inputparam.o m_io.o m_hamiltonian_onebody.f90
	$(FC_W_OPTS) -c m_hamiltonian_onebody.f90 $(LIBXC_INCS) $(LIBINT_INCS)

m_hamiltonian_tools.o: m_definitions.o m_mpi.o m_scalapack.o m_warning.o m_timing.o m_memory.o m_cart_to_pure.o \
	               m_linear_algebra.o m_inputparam.o m_basis_set.o m_hamiltonian_tools.f90
	$(FC_W_OPTS) -c m_hamiltonian_tools.f90 $(LIBXC_INCS)

m_hamiltonian_twobodies.o: m_definitions.o m_mpi.o m_scalapack.o m_warning.o m_timing.o m_memory.o m_cart_to_pure.o m_libxc_tools.o \
	         m_inputparam.o m_eri_calculate.o m_basis_set.o m_dft_grid.o m_density_tools.o m_io.o m_hamiltonian_tools.o \
                 m_hamiltonian_twobodies.f90
	$(FC_W_OPTS) -c m_hamiltonian_twobodies.f90 $(LIBXC_INCS)

m_hamiltonian_wrapper.o: m_definitions.o m_mpi.o m_scalapack.o m_warning.o m_timing.o m_basis_set.o m_scf.o \
	         m_inputparam.o m_hamiltonian_twobodies.o m_hamiltonian_tools.o m_hamiltonian_wrapper.f90
	$(FC_W_OPTS) -c m_hamiltonian_wrapper.f90 $(LIBXC_INCS)

m_density_tools.o: m_definitions.o m_atoms.o m_mpi.o m_gaussian.o m_inputparam.o m_basis_set.o m_hamiltonian_tools.o m_dft_grid.o \
                   m_density_tools.f90
	$(FC_W_OPTS) -c m_density_tools.f90 $(LIBXC_INCS)

m_multipole.o: m_definitions.o m_atoms.o m_inputparam.o m_basis_set.o m_hamiltonian_tools.o m_hamiltonian_onebody.o m_multipole.f90
	$(FC_W_OPTS) -c m_multipole.f90 $(LIBXC_INCS)

m_io.o: m_definitions.o m_mpi.o m_warning.o m_atoms.o m_memory.o m_cart_to_pure.o m_hdf5_tools.o \
	m_basis_set.o m_libint_tools.o m_libcint_tools.o m_inputparam.o m_eri.o \
        m_dft_grid.o git_sha.f90 m_timing.o m_linear_algebra.o m_string_tools.o \
        m_elements.o m_libxc_tools.o m_hamiltonian_tools.o m_io.f90
	$(FC_W_OPTS) -c m_io.f90 $(LIBXC_INCS) $(LIBINT_INCS)

m_restart.o: m_definitions.o m_timing.o m_mpi.o m_atoms.o m_basis_set.o m_hamiltonian_onebody.o m_inputparam.o m_linear_algebra.o m_hamiltonian_tools.o \
             m_restart.f90
	$(FC_W_OPTS) -c m_restart.f90 $(LIBXC_INCS)

m_nofoutput.o: m_definitions.o ./noft/m_nofoutput.F90
	$(FC_W_OPTS) -c ./noft/m_nofoutput.F90

m_lbfgs_intern.o: m_definitions.o ./noft/m_lbfgs_intern.F
	$(FC_W_OPTS) -c ./noft/m_lbfgs_intern.F

m_integd.o: m_definitions.o ./noft/m_integd.F90
	$(FC_W_OPTS) -c ./noft/m_integd.F90

m_rdmd.o: m_definitions.o ./noft/m_rdmd.F90
	$(FC_W_OPTS) -c ./noft/m_rdmd.F90

m_anti2unit.o: m_definitions.o ./noft/m_anti2unit.F90
	$(FC_W_OPTS) -c ./noft/m_anti2unit.F90

m_elag.o: m_definitions.o ./noft/m_elag.F90
	$(FC_W_OPTS) -c ./noft/m_elag.F90

m_hessian.o: m_definitions.o ./noft/m_hessian.F90
	$(FC_W_OPTS) -c ./noft/m_hessian.F90

m_gammatodm2.o: m_definitions.o ./noft/m_gammatodm2.F90
	$(FC_W_OPTS) -c ./noft/m_gammatodm2.F90

m_diagf.o: m_definitions.o ./noft/m_diagf.F90
	$(FC_W_OPTS) -c ./noft/m_diagf.F90

m_e_grad_occ_cpx.o: m_definitions.o ./noft/m_e_grad_occ_cpx.F90
	$(FC_W_OPTS) -c ./noft/m_e_grad_occ_cpx.F90

m_e_grad_occ.o: m_definitions.o ./noft/m_e_grad_occ.F90
	$(FC_W_OPTS) -c ./noft/m_e_grad_occ.F90

m_optocc.o: m_definitions.o ./noft/m_optocc.F90
	$(FC_W_OPTS) -c ./noft/m_optocc.F90

m_optorb.o: m_definitions.o ./noft/m_optorb.F90
	$(FC_W_OPTS) -c ./noft/m_optorb.F90

m_tz_pccd_amplitudes.o: m_definitions.o ./noft/m_tz_pccd_amplitudes.F90
	$(FC_W_OPTS) -c ./noft/m_tz_pccd_amplitudes.F90

m_noft_driver.o: m_definitions.o ./noft/m_noft_driver.F90
	$(FC_W_OPTS) -c ./noft/gitver.F90
	$(FC_W_OPTS) -c ./noft/m_noft_driver.F90
	
m_mp2_energy.o: m_definitions.o m_mpi.o m_cart_to_pure.o \
       m_basis_set.o m_inputparam.o m_eri_ao_mo.o \
       m_mp2_energy.f90
	$(FC_W_OPTS) -c m_mp2_energy.f90 $(LIBXC_INCS)

m_relativistic.o: m_definitions.o m_mpi.o m_cart_to_pure.o \
       m_basis_set.o m_inputparam.o m_eri_ao_mo.o \
       m_warning.o m_timing.o m_elements.o m_string_tools.o m_atoms.o m_ecp.o  m_gaussian.o m_hamiltonian_onebody.o m_hamiltonian_tools.o m_libcint_tools.o \
       m_relativistic.f90
	$(FC_W_OPTS) -c m_relativistic.f90 $(LIBXC_INCS)

m_noft.o: m_definitions.o m_mpi.o m_cart_to_pure.o \
       m_basis_set.o m_inputparam.o m_eri_ao_mo.o \
       m_hamiltonian_onebody.o m_hamiltonian_twobodies.o m_hamiltonian_tools.o m_hamiltonian_wrapper.o m_noft_driver.o \
       m_noft.f90
	$(FC_W_OPTS) -c m_noft.f90 $(LIBXC_INCS)

m_virtual_orbital_space.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_scalapack.o \
	                  m_basis_set.o m_eri_ao_mo.o m_hamiltonian_tools.o m_linear_algebra.o m_io.o \
                          m_inputparam.o m_hamiltonian_onebody.o m_virtual_orbital_space.f90
	$(FC_W_OPTS) -c m_virtual_orbital_space.f90 $(LIBXC_INCS)

m_pt_selfenergy.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_basis_set.o \
                 m_spectral_function.o m_inputparam.o m_selfenergy_tools.o m_eri_ao_mo.o \
                 m_linear_response.o m_gw_selfenergy_analytic.o m_pt_selfenergy.f90
	$(FC_W_OPTS) -c m_pt_selfenergy.f90 $(LIBXC_INCS)

m_pt_density_matrix.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_basis_set.o \
                  m_spectral_function.o m_inputparam.o m_selfenergy_tools.o m_eri_ao_mo.o m_hamiltonian_tools.o m_hamiltonian_onebody.o m_io.o m_pt_density_matrix.f90
	$(FC_W_OPTS) -c m_pt_density_matrix.f90 $(LIBXC_INCS)

pdbssolver1.o: pdbssolver1.f
	$(FC_W_OPTS) -c pdbssolver1.f $(LIBXC_INCS) $(ELPA_INCS)

boys_function.o: boys_function.f90
	$(FC_W_OPTS) -c boys_function.f90 $(LIBXC_INCS)

real_spherical_harmonics.o: m_definitions.o real_spherical_harmonics.f90
	$(FC_W_OPTS) -c real_spherical_harmonics.f90 $(LIBXC_INCS)

lebedev_quadrature.o: lebedev_quadrature.f
	$(FC_W_OPTS) -c lebedev_quadrature.f $(LIBXC_INCS)

libint_onebody.o: libint_molgw.h libint_onebody.cpp
	$(CXX_W_OPTS) -c libint_onebody.cpp $(LIBINT_INCS)

libint_onebody_grad.o: libint_molgw.h libint_onebody_grad.cpp
	$(CXX_W_OPTS) -c libint_onebody_grad.cpp $(LIBINT_INCS)

libint_twobody.o: libint_molgw.h libint_twobody.cpp
	$(CXX_W_OPTS) -c libint_twobody.cpp $(LIBINT_INCS)

libint_twobody_grad.o: libint_molgw.h libint_twobody_grad.cpp
	$(CXX_W_OPTS) -c libint_twobody_grad.cpp $(LIBINT_INCS)

m_libint_tools.o: m_definitions.o m_cart_to_pure.o m_basis_set.o m_libint_tools.f90
	$(FC_W_OPTS) -c m_libint_tools.f90 $(LIBXC_INCS) $(LIBINT_INCS)

m_libcint_tools.o: m_definitions.o m_cart_to_pure.o m_basis_set.o m_atoms.o m_libcint_tools.f90
	$(FC_W_OPTS) -c m_libcint_tools.f90 $(LIBXC_INCS) $(LIBINT_INCS)

m_dft_grid.o: m_definitions.o m_memory.o m_inputparam.o m_warning.o m_cart_to_pure.o m_timing.o m_elements.o m_atoms.o m_mpi.o m_basis_set.o m_numerical_tools.o \
              m_ecp.o m_dft_grid.f90
	$(FC_W_OPTS) -c m_dft_grid.f90 $(LIBXC_INCS)

m_ci.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_basis_set.o m_eri_ao_mo.o m_inputparam.o m_scalapack.o m_selfenergy_tools.o \
        m_linear_algebra.o m_hamiltonian_onebody.o m_ci.f90
	$(FC_W_OPTS) -c m_ci.f90 $(LIBXC_INCS)

m_spectral_function.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_scalapack.o m_eri.o m_eri_calculate.o m_inputparam.o m_numerical_tools.o \
                       m_basis_set.o m_eri_ao_mo.o m_spectral_function.f90
	$(FC_W_OPTS) -c m_spectral_function.f90 $(LIBXC_INCS)

m_eri_ao_mo.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_basis_set.o m_inputparam.o \
               m_scalapack.o m_eri.o m_eri_ao_mo.f90
	$(FC_W_OPTS) -c m_eri_ao_mo.f90 $(LIBXC_INCS)

m_spectra.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_basis_set.o m_inputparam.o m_mpi.o m_cart_to_pure.o \
               m_scalapack.o m_spectral_function.o m_hamiltonian_onebody.o m_spectra.f90
	$(FC_W_OPTS) -c m_spectra.f90 $(LIBXC_INCS)

m_tddft_fxc.o: m_definitions.o m_mpi.o m_timing.o m_warning.o m_memory.o m_basis_set.o m_inputparam.o m_dft_grid.o m_density_tools.o m_libxc_tools.o \
                m_tddft_fxc.f90
	$(FC_W_OPTS) -c m_tddft_fxc.f90 $(LIBXC_INCS)

m_eri_calculate.o: m_definitions.o m_mpi.o m_scalapack.o m_timing.o m_warning.o m_memory.o m_cart_to_pure.o m_basis_set.o m_inputparam.o \
                   m_libint_tools.o m_eri.o m_eri_calculate.f90
	$(FC_W_OPTS) -c m_eri_calculate.f90 $(LIBXC_INCS)

m_eri.o: m_definitions.o m_mpi.o m_timing.o m_memory.o m_warning.o m_scalapack.o m_cart_to_pure.o m_basis_set.o m_inputparam.o m_libint_tools.o m_libcint_tools.o \
         m_eri.f90
	$(FC_W_OPTS) -c m_eri.f90 $(LIBXC_INCS)

m_basis_set.o: m_definitions.o m_warning.o m_timing.o m_elements.o m_gaussian.o m_mpi.o m_ecp.o m_cart_to_pure.o m_atoms.o m_string_tools.o \
               m_basis_set.f90
	$(FC_W_OPTS) -c m_basis_set.f90 $(LIBXC_INCS)

m_gaussian.o: m_definitions.o m_mpi.o m_cart_to_pure.o m_string_tools.o m_gaussian.f90
	$(FC_W_OPTS) -c m_gaussian.f90 $(LIBXC_INCS)

m_scf.o: m_definitions.o m_mpi.o m_warning.o m_timing.o m_memory.o m_inputparam.o m_scalapack.o m_lbfgs.o m_linear_algebra.o \
	m_scf.f90
	$(FC_W_OPTS) -c m_scf.f90 $(LIBXC_INCS)

m_selfenergy_tools.o: m_definitions.o m_warning.o m_inputparam.o m_mpi.o m_atoms.o m_lbfgs.o m_basis_set.o m_dft_grid.o m_numerical_tools.o \
	              m_hamiltonian_tools.o m_hamiltonian_onebody.o m_hamiltonian_wrapper.o m_selfenergy_tools.f90
	$(FC_W_OPTS) -c m_selfenergy_tools.f90 $(LIBXC_INCS)

m_inputparam.o: m_definitions.o m_mpi.o  m_elements.o m_warning.o m_atoms.o m_ecp.o \
	        input_variables.f90 input_variable_declaration.f90 basis_path.f90 echo_input_variables.f90 echo_input_variables_yaml.f90 m_string_tools.o m_libxc_tools.o \
                m_inputparam.f90
	$(FC_W_OPTS) -c m_inputparam.f90 $(LIBXC_INCS) $(LIBINT_INCS)

m_block_diago.o: m_definitions.o m_warning.o m_timing.o m_memory.o m_mpi.o m_scalapack.o m_inputparam.o m_linear_algebra.o m_block_diago.f90
	$(FC_W_OPTS) -c m_block_diago.f90 $(LIBXC_INCS) $(ELPA_INCS)

m_atoms.o: m_definitions.o m_warning.o m_elements.o m_linear_algebra.o m_lbfgs.o m_atoms.f90
	$(FC_W_OPTS) -c m_atoms.f90 $(LIBXC_INCS)

m_ecp.o: m_definitions.o m_warning.o m_elements.o m_string_tools.o m_ecp.f90
	$(FC_W_OPTS) -c m_ecp.f90 $(LIBXC_INCS)

m_scalapack.o: m_definitions.o m_warning.o m_linear_algebra.o m_mpi.o m_scalapack_interface.o m_scalapack.f90
	$(FC_W_OPTS) -c m_scalapack.f90 $(LIBXC_INCS)

m_mpi.o: m_definitions.o m_warning.o m_mpi_tools.o m_mpi.f90
	$(FC_W_OPTS) -c m_mpi.f90 $(LIBXC_INCS)

m_mpi_tools.o: m_definitions.o m_warning.o  m_mpi_tools.f90
	$(FC_W_OPTS) -c m_mpi_tools.f90 $(LIBXC_INCS)

m_string_tools.o: m_definitions.o m_warning.o m_string_tools.f90
	$(FC_W_OPTS) -c m_string_tools.f90 $(LIBXC_INCS)

m_numerical_tools.o: m_definitions.o m_warning.o m_numerical_tools.f90
	$(FC_W_OPTS) -c m_numerical_tools.f90 $(LIBXC_INCS)

m_linear_algebra.o: m_definitions.o m_warning.o m_linear_algebra.f90
	$(FC_W_OPTS) -c m_linear_algebra.f90 $(LIBXC_INCS)

m_memory.o: m_definitions.o m_warning.o m_memory.f90
	$(FC_W_OPTS) -c m_memory.f90 $(LIBXC_INCS)

m_timing.o: m_definitions.o m_warning.o m_timing.f90
	$(FC_W_OPTS) -c m_timing.f90 $(LIBXC_INCS)

m_hdf5_tools.o: m_definitions.o m_hdf5_tools.f90
	$(FC_W_OPTS) -c m_hdf5_tools.f90 $(LIBXC_INCS) $(HDF5_INCS)

m_elements.o: m_definitions.o m_warning.o m_elements.f90
	$(FC_W_OPTS) -c m_elements.f90 $(LIBXC_INCS)

m_cart_to_pure.o: m_definitions.o m_warning.o m_cart_to_pure.f90
	$(FC_W_OPTS) -c m_cart_to_pure.f90 $(LIBXC_INCS)

m_warning.o: m_definitions.o m_warning.f90
	$(FC_W_OPTS) -c m_warning.f90 $(LIBXC_INCS)

m_libxc_tools.o: m_definitions.o m_warning.o m_libxc_tools.f90
	$(FC_W_OPTS) -c m_libxc_tools.f90 $(LIBXC_INCS)

m_scalapack_interface.o: m_scalapack_interface.f90
	$(FC_W_OPTS) -c m_scalapack_interface.f90 $(LIBXC_INCS)

m_lbfgs.o: m_definitions.o m_lbfgs.f90
	$(FC_W_OPTS) -c m_lbfgs.f90 $(LIBXC_INCS)

m_definitions.o: m_definitions.f90
	$(FC_W_OPTS) -c m_definitions.f90 $(LIBXC_INCS)

m_realtime_tddft.o: m_definitions.o m_basis_set.o m_hamiltonian_wrapper.o m_memory.o m_linear_algebra.o \
                      m_hamiltonian_tools.o m_hamiltonian_onebody.o m_string_tools.o m_io.o m_multipole.o \
                      m_inputparam.o m_dft_grid.o m_scf.o m_warning.o m_timing.o m_density_tools.o m_atoms.o m_hdf5_tools.o \
                      m_realtime_tddft.f90
	 $(FC_W_OPTS) -c m_realtime_tddft.f90 $(LIBXC_INCS)

tarball:
	+make archive

archive:
	+make prepare_sourcecode
	mkdir -p ../tmp/molgw/src
	mkdir -p ../tmp/molgw/src/noft
	cp -rp  ../README.md ../COPYING ../RELEASE.md  \
            ../basis/ ../utils/ ../tests/            \
            ../docs/ ../config/ ../Makefile  \
            ../tmp/molgw
	cp -rp  ../src/*.py ../src/*.f90 ../src/*.f  \
            ../src/*.cpp ../src/*.h  \
            ../src/Makefile \
            ../tmp/molgw/src
	cp -rp  ../src/noft/*.F90 ../src/noft/*.F  \
            ../src/noft/*.sh  \
            ../tmp/molgw/src/noft/
	cd ../tmp && tar czf ../molgw.tgz molgw/
	rm -rf ../tmp/

prepare_sourcecode:
	chmod 777 ./noft/gitversion.sh
	./noft/gitversion.sh > /dev/null
	mv gitver.F90 ./noft 
	$(PYTHON) prepare_sourcecode.py $(PREFIX)

input:
	$(PYTHON) ../utils/input_variables.py

docs:
	pdflatex ../docs/notes.tex
	pdflatex ../docs/notes.tex
	rm -f notes.aux notes.log notesNotes.bib
	mv notes.pdf ../docs/notes.pdf

clean:
	$(RM) *.o *.mod ../molgw.tgz ../molgw ./noft/gitver.F90

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
