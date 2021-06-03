all:
	gfortran -c m_vars.F90
	gfortran -c m_cg.F90
	gfortran -w -ffixed-form -fd-lines-as-comments -c m_lbfgs.F
	gfortran -c m_integd.F90
	gfortran -c m_rdmd.F90
	gfortran -c m_elag.F90
	gfortran -c m_diagf.F90
	gfortran -c m_gammatodm2.F90
	gfortran -c m_e_grad_occ.F90
	gfortran -c m_optocc.F90
	gfortran -c m_optorb.F90
	gfortran -c m_noft_driver.F90
clean:
	/bin/rm -rf *.o *.mod
tar:
	tar -cvf module_noft.tar *F90 *F
