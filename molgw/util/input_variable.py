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
        datatype =''
        default  =''
        comment  =''
	def printhtml(self,f):
		f.write('<hr>')
		f.write('<a name='+self.keyword+'>')
		f.write('<li><b>'+self.keyword+'</b><br><br>')
		f.write('Default: '+str(self.default)+'<br><br>')
		f.write(self.comment+'</li><br>')
		


vl = []
#vl = [variable(),variable()]


#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='scf'
vl[i].datatype ='characters'
vl[i].comment  ='Contains the self-consistent scheme name. \n\
Try LDA, PBE, HSE06, or HF for instance'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='postscf'
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
vl[i].default  =0.25
vl[i].datatype ='real'
vl[i].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the amount of range-independent exact-exchange'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='beta_mixing'
vl[i].default  =0.
vl[i].datatype ='real'
vl[i].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the amount of long-range exact-exchange'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='gamma_mixing'
vl[i].default  =1000000.
vl[i].datatype ='real'
vl[i].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the separation between long-range and short-range. It is input in bohr^-1.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='basis'
vl[i].datatype ='character'
vl[i].comment  ='Sets the basis set \
For Pople sets, use 6-31G for instance or 6-31pGs, where p stands for + and s for *. \
For Dunning sets, use aug-cc-pVTZ for instance. \
Note that Pople sets are to be used with gaussian_type=\'cart\' \
One may use ones own basis sets provided that the files are labeled X_mybasisset where X is the element.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='auxil_basis'
vl[i].datatype ='character'
vl[i].comment  ='Sets the auxiliary basis set. \
For instance, cc-pVDZ-RI for a Weigend basis set. \
If present, the auxiliary basis will be used for postscf calculations (TD-DFT, BSE, or GW) \
If specifically requested with no_4center, the auxiliary basis can be used for scf cycles too.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='basis_path'
vl[i].default  ='./'
vl[i].datatype ='character'
vl[i].comment  ='Sets the path pointing to the basis functions files.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='gaussian_type'
vl[i].default  ='pure'
vl[i].datatype ='character'
vl[i].comment  ='Asks for pure or spherical Gaussian type orbitals with \'pure\' \
or for Cartesian Gaussian orbital with \'cart\'.'

#================================
vl.append(variable())
i = len(vl) - 1
vl[i].keyword  ='no_4center'
vl[i].default  ='no'
vl[i].datatype ='yes/no'
vl[i].comment  ='If switched on, the auxiliary basis set is used in both SCF cycles and in post-scf methods.\
This avoids the calculation and the storage of the 4-center Coulomb integrals.'







#================================
fhtml = open('../doc/input_variables.html','w')

fhtml.write('<html>\n')
fhtml.write('<head>\n')
fhtml.write('<link rel="stylesheet" type="text/css" href="molgw.css">\n')
fhtml.write('</head>\n')

fhtml.write('<a name=top>')
fhtml.write('<h1>Input variable list</h1>')
fhtml.write('<br><br>')
for i in range(0,len(vl)):
	fhtml.write('<a href=#'+vl[i].keyword+'>'+vl[i].keyword+'</a> ')
fhtml.write('<br><br><ul>')
for i in range(0,len(vl)):
	vl[i].printhtml(fhtml)
fhtml.write('</ul>')
fhtml.write('<br><br><br><br><br><br>')
fhtml.write('<a href=#top>Back to the top of the page</a> ')
fhtml.write('</html>\n')




