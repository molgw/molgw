#!/usr/bin/python

# namelist /molgw/ scf,postscf,                                                                          &
#                  alpha_hybrid,beta_hybrid,gamma_hybrid,                                                &
#                  basis,auxil_basis,basis_path,gaussian_type,no_4center,                                &
#                  nspin,charge,magnetization,                                                           &
#                  grid_quality,integral_quality,                                                        &
#                  nscf,alpha_mixing,mixing_scheme,tolscf,                                               &
#                  tda,triplet,eta,frozencore,ncoreg,ncorew,nvirtualg,nvirtualw,nomega_sigma,step_sigma, &
#                  ignore_restart,ignore_bigrestart,print_matrix,print_eri,print_wfn,print_w,print_sigma,&
#                  length_unit,natom

class variable:
        keyword  =''
	family   =''
        datatype =''
	mandatory='no'
        default  =''
        comment  =''
	def printhtml(self,f):
		f.write('<hr>')
		f.write('<a name='+self.keyword+'>')
		f.write('<li><b>'+self.keyword+'</b><br><br>')
		if self.mandatory == 'yes':
			f.write('<i>Mandatory</i><br>')
        	else:
			f.write('<i>Optional</i><br>')
		if self.default == '':
			f.write('Default: None<br><br>')
	 	else:
			f.write('Default: '+str(self.default)+'<br><br>')
		f.write(self.comment+'</li><br>')
		


vl = []


#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='scf'
vl[i].family   ='general'
vl[i].datatype ='characters'
vl[i].mandatory='yes'
vl[i].comment  ='Contains the self-consistent scheme name. \n\
Try LDA, PBE, HSE06, or HF for instance'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='postscf'
vl[i].family   ='post'
vl[i].datatype ='characters'
vl[i].comment  ='Contains the post-processing scheme name. \n\
TD stands for TD-DFT\n\
BSE stands for Bethe-Salpeter\n\
GW stands for perturbative G0W0\n\
GnW0 stands for GW with eigenvalue self-consistentcy on G\n\
GnWn stands for GW with eigenvalue self-consistentcy on both G and W\n\
MP2 stands for guess what.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='alpha_mixing'
vl[i].family   ='scf'
vl[i].default  =0.25
vl[i].datatype ='real'
vl[i].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the amount of range-independent exact-exchange'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='beta_mixing'
vl[i].family   ='scf'
vl[i].default  =0.
vl[i].datatype ='real'
vl[i].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the amount of long-range exact-exchange'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='gamma_mixing'
vl[i].family   ='scf'
vl[i].default  =1000000.
vl[i].datatype ='real'
vl[i].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the separation between long-range and short-range. It is input in bohr^-1.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='basis'
vl[i].family   ='general'
vl[i].datatype ='character'
vl[i].mandatory='yes'
vl[i].comment  ='Sets the basis set \
For Pople sets, use 6-31G for instance or 6-31pGs, where p stands for + and s for *. \
For Dunning sets, use aug-cc-pVTZ for instance. \
Note that Pople sets are to be used with gaussian_type=\'cart\' \
One may use ones own basis sets provided that the files are labeled X_mybasisset where X is the element.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='auxil_basis'
vl[i].family   ='general'
vl[i].datatype ='character'
vl[i].comment  ='Sets the auxiliary basis set. \
For instance, cc-pVDZ-RI for a Weigend basis set. \
If present, the auxiliary basis will be used for postscf calculations (TD-DFT, BSE, or GW) \
If specifically requested with no_4center, the auxiliary basis can be used for scf cycles too.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='basis_path'
vl[i].family   ='io'
vl[i].default  ='./'
vl[i].datatype ='character'
vl[i].comment  ='Sets the path pointing to the basis functions files.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='gaussian_type'
vl[i].family   ='general'
vl[i].default  ='pure'
vl[i].datatype ='character'
vl[i].comment  ='Asks for pure or spherical Gaussian type orbitals with \'pure\' \
or for Cartesian Gaussian orbital with \'cart\'.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='no_4center'
vl[i].family   ='scf'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='If switched on, the auxiliary basis set is used in both SCF cycles and in post-scf methods.\
This avoids the calculation and the storage of the 4-center Coulomb integrals.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='nspin'
vl[i].family   ='system'
vl[i].default  =1
vl[i].datatype ='integer'
vl[i].comment  ='Sets the number of spin channels. 1 enforces spin-restricted calculations. \
2 means spin-unrestricted.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='charge'
vl[i].family   ='system'
vl[i].default  =0.0
vl[i].datatype ='real'
vl[i].comment  ='Sets the total charge of the system. 0 is a neutral system. -2 is a doubly charged anion etc.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='magnetization'
vl[i].family   ='system'
vl[i].default  ='0.0'
vl[i].datatype ='real'
vl[i].comment  ='Sets the number of unpaired electrons. In other words, this is the difference between \
the spin up and spin down occupation. For instance, a spin-doublet calculation is obtained with magnetization=1.0. \
Only meaningful when nspin=2.' 

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='grid_quality'
vl[i].family   ='scf'
vl[i].default  ='high'
vl[i].datatype ='quality'
vl[i].comment  ='Sets the number of grid points use to evaluate the exchange-correlation integrals in real space. \
Possible values are \'low\', \'medium\', \'high\', \'very high\', \'insane\'. \
It could be abbreviated in \'l\', \'m\', \'h\', \'vh\', \'i\'. \
\'high\' is usually fine. \'insane\' is only meant for debugging since it is overdoing a lot.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='integral_quality'
vl[i].family   ='scf'
vl[i].default  ='high'
vl[i].datatype ='quality'
vl[i].comment  ='Sets the tolerance value for the screening of the negligible integrals. \
Possible values are \'low\', \'medium\', \'high\', \'very high\', \'insane\'. \
It could be abbreviated in \'l\', \'m\', \'h\', \'vh\', \'i\'. \
\'high\' is usually fine. \'insane\' is only meant for debugging since it is overdoing a lot.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='nscf'
vl[i].family   ='scf'
vl[i].default  =30
vl[i].datatype ='integer'
vl[i].comment  ='Sets the maximum number of SCF cycles'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='alpha_mixing'
vl[i].family   ='scf'
vl[i].default  =0.7
vl[i].datatype ='real'
vl[i].comment  ='Sets the amount of output density-matrix for the next iteration. \
When the SCF cycles have difficulties to converge, one may try to lower this value.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='mixing_scheme'
vl[i].family   ='scf'
vl[i].default  ='pulay'
vl[i].datatype ='characters'
vl[i].comment  ='Sets the density-matrix update method for SCF cycles. \
Possible choices are \'pulay\' for Pulay DIIS method or \'simple\' for a simple linear mixing between input and output density-matrices.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='tolscf'
vl[i].family   ='scf'
vl[i].default  =1.0E-7
vl[i].datatype ='real'
vl[i].comment  ='Sets the residual norm target for the SCF cycles.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='tda'
vl[i].family   ='post'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Triggers the use of Tamm-Dancoff approximation in TD-DFT or BSE.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='triplet'
vl[i].family   ='post'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Triggers the calculation of the triplet final state in TD-DFT or BSE.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='frozencore'
vl[i].family   ='post'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Triggers the neglect of core states in GW. \
H, He, Li, Be have no core states. B-Na have the 1s. \
Al-Ca have the 1s2s2p. Manual tuning could be achieved with ncoreg, ncorew.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='ncoreg'
vl[i].family   ='post'
vl[i].default  =0
vl[i].datatype ='integer'
vl[i].comment  ='Sets the number of frozen core states in the Green\'s function G.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='ncorew'
vl[i].family   ='post'
vl[i].default  =0
vl[i].datatype ='integer'
vl[i].comment  ='Sets the number of frozen core states in the screened Coulomb interaction W, in TD-DFT, and in BSE.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='nvirtualg'
vl[i].family   ='post'
vl[i].default  =100000
vl[i].datatype ='integer'
vl[i].comment  ='Sets the state beyond which they are excluded from the sum in the Green\'s function G.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='nvirtualw'
vl[i].family   ='post'
vl[i].default  =100000
vl[i].datatype ='integer'
vl[i].comment  ='Sets the state beyond which they are excluded from the sum in the screened Coulomb interaction W, in TD-DFT, and in BSE.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='nomega_sigma'
vl[i].family   ='post'
vl[i].default  =51
vl[i].datatype ='integer'
vl[i].comment  ='Sets the number of frequencies used to solve the quasiparticle equation in the GW self-energy.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='step_sigma'
vl[i].family   ='post'
vl[i].default  =0.01
vl[i].datatype ='real'
vl[i].comment  ='Sets the spacing between frequencies in the GW self-energy evaluation.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='ignore_restart'
vl[i].family   ='io'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Ignore the RESTART file and restart from scratch.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='ignore_bigrestart'
vl[i].family   ='io'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Considers a big RESTART as if it was a small RESTART.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='print_matrix'
vl[i].family   ='io'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Prints some matrices for debugging purposes.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='print_eri'
vl[i].family   ='io'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Dumps the Electron Repulsion Integral on a file.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='print_wfn'
vl[i].family   ='io'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Prints some wavefunctions along some selected lines.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='print_w'
vl[i].family   ='io'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Dumps the spectral function of the screened Coulomb W. This is necessary for a subsequent BSE run.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='print_sigma'
vl[i].family   ='io'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='Prints the value of the GW self-energy on the sampling frequencies in files.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='length_unit'
vl[i].family   ='system'
vl[i].default  ='angstrom'
vl[i].datatype ='characters'
vl[i].comment  ='Chooses the units of the atomic coordinates. Can be \'angstrom\' or \'bohr\'. \
Could be abbreviated in \'A\' or \'au\'.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='natom'
vl[i].family   ='system'
vl[i].mandatory='yes'
vl[i].default  =''
vl[i].datatype ='integer'
vl[i].comment  ='Sets the number of atoms in the molecule. This is the number of lines to be read in the following section of the input file.'

















#================================
fhtml = open('../doc/input_variables.html','w')

fhtml.write('<html>\n')
fhtml.write('<head>\n')
fhtml.write('<link rel="stylesheet" type="text/css" href="molgw.css">\n')
fhtml.write('</head>\n')

fhtml.write('<a name=top>')
fhtml.write('<h1>Input variable list</h1>')
fhtml.write('<br><br>')

# Mandatory
fhtml.write('<h3>Mandatory input variables</h3>')
for i in range(0,len(vl)):
	if vl[i].mandatory =='yes':
		fhtml.write('<a href=#'+vl[i].keyword+'>'+vl[i].keyword+'</a> ')

# System
fhtml.write('<h3>Set up of the physical system input variables</h3>')
for i in range(0,len(vl)):
	if vl[i].family =='system':
		fhtml.write('<a href=#'+vl[i].keyword+'>'+vl[i].keyword+'</a> ')

# General
fhtml.write('<h3>General input variables</h3>')
for i in range(0,len(vl)):
	if vl[i].family =='general':
		fhtml.write('<a href=#'+vl[i].keyword+'>'+vl[i].keyword+'</a> ')

# SCF
fhtml.write('<h3>Self-consistency input variables</h3>')
for i in range(0,len(vl)):
	if vl[i].family =='scf':
		fhtml.write('<a href=#'+vl[i].keyword+'>'+vl[i].keyword+'</a> ')

# Post 
fhtml.write('<h3>Correlation and excited states post-treatment input variables</h3>')
for i in range(0,len(vl)):
	if vl[i].family =='post':
		fhtml.write('<a href=#'+vl[i].keyword+'>'+vl[i].keyword+'</a> ')



# IO family
fhtml.write('<h3>IO input variables</h3>')
for i in range(0,len(vl)):
	if vl[i].family =='io':
		fhtml.write('<a href=#'+vl[i].keyword+'>'+vl[i].keyword+'</a> ')


# Start the complete list
fhtml.write('<br><br><br><br>')
fhtml.write('<h2>Complete list of input variables</h2>')
fhtml.write('<br><br><ul>')
for i in range(0,len(vl)):
	vl[i].printhtml(fhtml)
fhtml.write('</ul>')
fhtml.write('<br><br><br><br><br><br>')
fhtml.write('<a href=#top>Back to the top of the page</a> ')
fhtml.write('</html>\n')




