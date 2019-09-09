#!/usr/bin/python3
# This file is part of MOLGW
# Author: Fabien Bruneval

import time,os
from yaml        import load,dump
from yaml        import CLoader as Loader, CDumper as Dumper

today=time.strftime("%d")+' '+time.strftime("%B")+' '+time.strftime("%Y")


class variable:
  keyword  =''
  family   =''
  datatype =''
  mandatory='no'
  experimental='no'
  default  =''
  comment  =''
  def printhtml(self,f):
    f.write('<hr>\n')
    f.write('<a name='+self.keyword+'>')
    f.write('<li>    \
             <span style="display:inline-block;background:#EEEEEE;width:400px">  \
             <b>'+self.keyword+'</b>  </span>  \n')

    if self.experimental == 'yes':
      f.write('<b><font color="red">EXPERIMENTAL</font> </b> \n')
    f.write('<br><br>\n')

    if self.mandatory == 'yes':
      f.write('<i>Mandatory</i><br>\n')
    else:
      f.write('<i>Optional</i><br>\n')
    if self.default == '':
      f.write('Default: None<br><br>\n')
    else:
      f.write('Default: '+str(self.default)+'<br><br>\n')
    f.write(self.comment+'</li><br>\n')


vl = []

#================================
vl.append(variable())
vl[-1].keyword  ='comment'
vl[-1].family   ='general'
vl[-1].datatype ='characters'
vl[-1].default  =''
vl[-1].comment  ='This is a free expression place. Use it as you wish for commenting, naming, labeling etc. (140 character max just as twitter)'

#================================
vl.append(variable())
vl[-1].keyword  ='scf'
vl[-1].family   ='general'
vl[-1].datatype ='characters'
vl[-1].mandatory='yes'
vl[-1].comment  ='Contains the self-consistent scheme name. \n\
Try LDA, PBE, HSE06, or HF for instance'

#================================
vl.append(variable())
vl[-1].keyword  ='postscf'
vl[-1].family   ='post'
vl[-1].datatype ='characters'
vl[-1].comment  ='Contains the post-processing scheme name. \n\
TD stands for TD-DFT or TD-HF.\n\
BSE stands for Bethe-Salpeter.\n\
GW stands for perturbative G0W0.\n\
GnW0 stands for GW with eigenvalue self-consistentcy on G.\n\
GnWn stands for GW with eigenvalue self-consistentcy on both G and W.\n\
MP2 stands for guess what.\n\
GWGAMMA (EXPERIMENTAL) stands for vertex corrections.'

#================================
vl.append(variable())
vl[-1].keyword  ='memory_evaluation'
vl[-1].family   ='general'
vl[-1].datatype ='yes/no'
vl[-1].default  ='no'
vl[-1].comment  ='Requests a memory evaluation. MOLGW will start normaly, evaluate the size of the arrays, and exit without performing an actual calculation.'

#================================
vl.append(variable())
vl[-1].keyword  ='move_nuclei'
vl[-1].family   ='general'
vl[-1].default  ='no'
vl[-1].datatype ='characters'
vl[-1].comment  ='Tells the code to move or not the position of the nuclei. \
Available options are \'no\' or \'relax\'.'

#================================
vl.append(variable())
vl[-1].keyword  ='nstep'
vl[-1].family   ='general'
vl[-1].default  = 50
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of steps when moving the nuclei.'

#================================
vl.append(variable())
vl[-1].keyword  ='tolforce'
vl[-1].family   ='general'
vl[-1].default  = 1.0e-5
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the target threshold for the maximum force component after nuclei relaxation.'

#================================
vl.append(variable())
vl[-1].keyword  ='alpha_hybrid'
vl[-1].family   ='scf'
vl[-1].default  =0.
vl[-1].datatype ='real'
vl[-1].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the amount of range-independent exact-exchange'

#================================
vl.append(variable())
vl[-1].keyword  ='beta_hybrid'
vl[-1].family   ='scf'
vl[-1].default  =0.
vl[-1].datatype ='real'
vl[-1].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the amount of long-range exact-exchange'

#================================
vl.append(variable())
vl[-1].keyword  ='gamma_hybrid'
vl[-1].family   ='scf'
vl[-1].default  =1000000.
vl[-1].datatype ='real'
vl[-1].comment  ='Only works for Range-Separated hybrid functionals scf=\'rsh\' \
Sets the separation between long-range and short-range. It is input in bohr^-1.'

#================================
vl.append(variable())
vl[-1].keyword  ='basis'
vl[-1].family   ='general'
vl[-1].datatype ='characters'
vl[-1].mandatory='yes'
vl[-1].comment  ='Sets the basis set \
For Pople sets, use 6-31G for instance or 6-31pGs, where p stands for + and s for *. \
For Dunning sets, use aug-cc-pVTZ for instance. \
Note that Pople sets are to be used with gaussian_type=\'cart\' \
One may use ones own basis sets provided that the files are labeled X_mybasisset where X is the element.'

#================================
vl.append(variable())
vl[-1].keyword  ='auxil_basis'
vl[-1].family   ='general'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the auxiliary basis set. \
For instance, cc-pVDZ-RI for a Weigend basis set. \
If present, the auxiliary basis will be used for both the scf cycles and the postscf calculations (TD-DFT, BSE, or GW).'

#================================
vl.append(variable())
vl[-1].keyword  ='auto_auxil_fsam'
vl[-1].family   ='general'
vl[-1].default  = 1.5
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the F_SAM parameter in the automatic generation of the auxiliary basis set. \
The closer to 1.0 the more auxiliary basis functions it will generate. \
See Yang-Rendell-Frisch, JChemPhys 127, 074102 (2007) for more details.'

#================================
vl.append(variable())
vl[-1].keyword  ='auto_auxil_lmaxinc'
vl[-1].family   ='general'
vl[-1].default  = 1
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the l_MAXINC parameter in the automatic generation of the auxiliary basis set. The larger the more auxiliary basis functions it will generate. \
See Yang-Rendell-Frisch, JChemPhys 127, 074102 (2007) for more details.'

#================================
vl.append(variable())
vl[-1].keyword  ='incore'
vl[-1].family   ='general'
vl[-1].datatype ='yes/no'
vl[-1].default  ='yes'
vl[-1].comment  ='Specify if the 4-center integrals are all calculated at once and stored or if they are calculated on-the-fly.'

#================================
vl.append(variable())
vl[-1].keyword  ='basis_path'
vl[-1].family   ='io'
vl[-1].default  =''
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the path pointing to the basis functions files. \
                 If not specified, then the basis set files will be searched in folder ~molgw/basis/.'

#================================
vl.append(variable())
vl[-1].keyword  ='small_basis'
vl[-1].family   ='postscf'
vl[-1].default  =''
vl[-1].experimental='yes'
vl[-1].datatype ='characters'
vl[-1].comment  ='Calls for a smaller basis set used to represent the virtual orbital space with fewer functions. \
Only meaningful for GW.'

#================================
vl.append(variable())
vl[-1].keyword  ='ecp_small_basis'
vl[-1].family   ='postscf'
vl[-1].default  =''
vl[-1].experimental='yes'
vl[-1].datatype ='characters'
vl[-1].comment  ='Calls for a smaller basis set used to represent the virtual orbital space with fewer functions. \
This is the small basis set used for elements with an effective core potential. Only meaningful for GW.'


#================================
vl.append(variable())
vl[-1].keyword  ='gaussian_type'
vl[-1].family   ='general'
vl[-1].default  ='pure'
vl[-1].datatype ='characters'
vl[-1].comment  ='Asks for pure or spherical Gaussian type orbitals with \'pure\' \
or for Cartesian Gaussian orbital with \'cart\'.'

#================================
vl.append(variable())
vl[-1].keyword  ='nspin'
vl[-1].family   ='system'
vl[-1].default  =1
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of spin channels. 1 enforces spin-restricted calculations. \
2 means spin-unrestricted.'

#================================
vl.append(variable())
vl[-1].keyword  ='charge'
vl[-1].family   ='system'
vl[-1].default  =0.0
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the total charge of the system. 0 is a neutral system. -2 is a doubly charged anion etc.'

#================================
vl.append(variable())
vl[-1].keyword  ='magnetization'
vl[-1].family   ='system'
vl[-1].default  ='0.0'
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the number of unpaired electrons. In other words, this is the difference between \
the spin up and spin down occupation. For instance, a spin-doublet calculation is obtained with magnetization=1.0. \
Only meaningful when nspin=2.' 

#================================
vl.append(variable())
vl[-1].keyword  ='temperature'
vl[-1].family   ='system'
vl[-1].default  ='0.0'
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the electronic temperature in the Fermi-Dirac functions. Helps the convergence for some systems. \
The value is input in Hartree atomic units.'

#================================
vl.append(variable())
vl[-1].keyword  ='grid_quality'
vl[-1].family   ='scf'
vl[-1].default  ='high'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the number of grid points use to evaluate the exchange-correlation integrals in real space for the DFT potential and energy. \
Possible values are \'low\', \'medium\', \'high\', \'very high\', \'insane\'. \
It could be abbreviated in \'l\', \'m\', \'h\', \'vh\', \'i\'. \
\'high\' is usually fine. \'insane\' is only meant for debugging since it is overdoing a lot.'

#================================
vl.append(variable())
vl[-1].keyword  ='tddft_grid_quality'
vl[-1].family   ='post'
vl[-1].default  ='high'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the number of grid points use to evaluate the exchange-correlation integrals in real space for the TDDFT kernel. \
Possible values are \'low\', \'medium\', \'high\', \'very high\', \'insane\'. \
It could be abbreviated in \'l\', \'m\', \'h\', \'vh\', \'i\'. \
\'high\' is usually fine. \'insane\' is only meant for debugging since it is overdoing a lot.'

#================================
vl.append(variable())
vl[-1].keyword  ='integral_quality'
vl[-1].family   ='scf'
vl[-1].default  ='high'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the tolerance value for the screening of the negligible integrals. \
Possible values are \'low\', \'medium\', \'high\', \'very high\', \'insane\'. \
It could be abbreviated in \'l\', \'m\', \'h\', \'vh\', \'i\'. \
\'high\' is usually fine. \'insane\' is only meant for debugging since it is overdoing a lot.'

#================================
vl.append(variable())
vl[-1].keyword  ='partition_scheme'
vl[-1].family   ='scf'
vl[-1].default  ='ssf'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the partition scheme for the xc quadrature. \
Possible choices are \'becke\' or \'ssf\' (Stratmann-Scuseria-Frisch).'

#================================
vl.append(variable())
vl[-1].keyword  ='nscf'
vl[-1].family   ='scf'
vl[-1].default  =50
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the maximum number of SCF cycles'

#================================
vl.append(variable())
vl[-1].keyword  ='alpha_mixing'
vl[-1].family   ='scf'
vl[-1].default  =0.7
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the amount of output density-matrix for the next iteration. \
When the SCF cycles have difficulties to converge, one may try to lower this value.'

#================================
vl.append(variable())
vl[-1].keyword  ='mixing_scheme'
vl[-1].family   ='scf'
vl[-1].default  ='pulay'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the density-matrix update method for SCF cycles. \
Possible choices are \'pulay\' for Pulay DIIS method, \'adiis\' for Hu-Yang method, or \'simple\' for a simple linear mixing between input and output density-matrices.'

#================================
vl.append(variable())
vl[-1].keyword  ='diis_switch'
vl[-1].family   ='scf'
vl[-1].default  =0.05
vl[-1].datatype ='real'
vl[-1].comment  ='When running ADIIS, sets the residue value below which the DIIS method is used to finalize the convergence.'

#================================
vl.append(variable())
vl[-1].keyword  ='level_shifting_energy'
vl[-1].family   ='scf'
vl[-1].default  =0.
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the energy shift up of the unoccupied states. Should help the convergence in the case of small HOMO-LUMO gaps.'

#================================
vl.append(variable())
vl[-1].keyword  ='init_hamiltonian'
vl[-1].family   ='scf'
vl[-1].default  ='guess'
vl[-1].datatype ='characters'
vl[-1].comment  ='Selects how to initiate the first hamiltonian for SCF cycles. Today, two options are available: \'guess\' for an educated guess based on approximate atomic densities \
or \'core\' for the core hamiltonian.'

#================================
vl.append(variable())
vl[-1].keyword  ='tolscf'
vl[-1].family   ='scf'
vl[-1].default  =1.0e-7
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the residual norm target for the density matrix for the SCF cycles.'

#================================
vl.append(variable())
vl[-1].keyword  ='density_matrix_damping'
vl[-1].family   ='scf'
vl[-1].default  =0.0
vl[-1].datatype ='real'
vl[-1].comment  ='Adds an additional linear mixing on the density matrix in combination with the Hamiltonian mixing in order to damp out the charge oscillations. Especially useful for metallic systems.'

#================================
vl.append(variable())
vl[-1].keyword  ='kerker_k0'
vl[-1].family   ='scf'
vl[-1].default  =0.0
vl[-1].datatype ='real'
vl[-1].experimental ='yes'
vl[-1].comment  ='Analog to k0 in Kerker preconditioning for metallic systems. Helps to damp charge oscillations to ensure better SCF convergence.'

#================================
vl.append(variable())
vl[-1].keyword  ='min_overlap'
vl[-1].family   ='scf'
vl[-1].default  =1.0e-5
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the minimal eigenvalue of the overlap matrix S. Small eigenvalues imply overcompleteness of the basis set.'

#================================
vl.append(variable())
vl[-1].keyword  ='npulay_hist'
vl[-1].family   ='scf'
vl[-1].default  =6
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the history record length for Pulay DIIS.'

#================================
vl.append(variable())
vl[-1].keyword  ='scf_diago_flavor'
vl[-1].family   ='scf'
vl[-1].default  =' '
vl[-1].datatype ='characters'
vl[-1].comment  ='Selects the LAPACK/ScaLAPACK diagonalization routines in the SCF cycles. Available choices are \' \', \'R\', \'D\', and \'X\'.'

#================================
vl.append(variable())
vl[-1].keyword  ='postscf_diago_flavor'
vl[-1].family   ='postscf'
vl[-1].default  ='R'
vl[-1].datatype ='characters'
vl[-1].comment  ='Selects the LAPACK/ScaLAPACK diagonalization routines in the post SCF calculations. Available choices are \' \', \'R\', \'D\', and \'X\'.'

#================================
vl.append(variable())
vl[-1].keyword  ='nstep_gw'
vl[-1].family   ='post'
vl[-1].default  ='1'
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of GW iterations for eigenvalue self-consistent GW calculations (GnWn or GnW0).'

#================================
vl.append(variable())
vl[-1].keyword  ='tda'
vl[-1].family   ='post'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Triggers the use of Tamm-Dancoff approximation in TD-DFT or BSE.'

#================================
vl.append(variable())
vl[-1].keyword  ='triplet'
vl[-1].family   ='post'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Triggers the calculation of the triplet final state in TD-DFT or BSE.'

#================================
vl.append(variable())
vl[-1].keyword  ='nexcitation'
vl[-1].family   ='post'
vl[-1].default  =0
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of neutral excitations to be calculated in TD-DFT or BSE. \
                 0 stands for all the states and triggers the full diagonalization.'

#================================
vl.append(variable())
vl[-1].keyword  ='toldav'
vl[-1].family   ='post'
vl[-1].default  =1.0e-4
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the tolerance criterium for the maximum norm of the residual in the Davidson diagonalization of TD-DFT, BSE, and full CI.'

#================================
vl.append(variable())
vl[-1].keyword  ='nstep_dav'
vl[-1].family   ='post'
vl[-1].default  =15
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the maximum number of Davidson partial diagonalization steps. Used for TD-DFT, BSE, and full CI.'

#================================
vl.append(variable())
vl[-1].keyword  ='frozencore'
vl[-1].family   ='post'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Triggers the neglect of core states in GW. \
H, He, Li, Be have no core states. B-Na have the 1s. \
Al-Ca have the 1s2s2p. Manual tuning could be achieved with ncoreg, ncorew.'

#================================
vl.append(variable())
vl[-1].keyword  ='ncoreg'
vl[-1].family   ='post'
vl[-1].default  =0
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of frozen core states in the Green\'s function G.'

#================================
vl.append(variable())
vl[-1].keyword  ='ncorew'
vl[-1].family   ='post'
vl[-1].default  =0
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of frozen core states in the screened Coulomb interaction W, in TD-DFT, and in BSE.'

#================================
vl.append(variable())
vl[-1].keyword  ='nvirtualg'
vl[-1].family   ='post'
vl[-1].default  =100000
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the starting state beyond which states are excluded from the sum in the Green\'s function G.'

#================================
vl.append(variable())
vl[-1].keyword  ='nvirtualw'
vl[-1].family   ='post'
vl[-1].default  =100000
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the starting state beyond which states are excluded from the sum in the screened Coulomb interaction W, in TD-DFT, and in BSE.'

#================================
vl.append(variable())
vl[-1].keyword  ='nomega_imag'
vl[-1].family   ='post'
vl[-1].default  = 0
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of frequencies used to perform the integral on the imaginary axis'

#================================
vl.append(variable())
vl[-1].keyword  ='selfenergy_state_min'
vl[-1].family   ='post'
vl[-1].default  =1
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the starting states for the range of the self-energy evaluation'

#================================
vl.append(variable())
vl[-1].keyword  ='selfenergy_state_max'
vl[-1].family   ='post'
vl[-1].default  =100000
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the final states for the range of the self-energy evaluation'

#================================
vl.append(variable())
vl[-1].keyword  ='selfenergy_state_range'
vl[-1].family   ='post'
vl[-1].default  =100000
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the range of states around the HOMO level for the self-energy evaluation. For instance, selfenergy_state_range=0 will trigger the calculation of the HOMO only. \
            selfenergy_state_range=1 will trigger the evaluation of the HOMO-1, HOMO, HOMO+1. etc.'

#================================
vl.append(variable())
vl[-1].keyword  ='nomega_sigma'
vl[-1].family   ='post'
vl[-1].default  =51
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of frequencies used to solve the quasiparticle equation in the GW self-energy.'

#================================
vl.append(variable())
vl[-1].keyword  ='step_sigma'
vl[-1].family   ='post'
vl[-1].default  =0.01
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the spacing between frequencies in the GW self-energy evaluation.'

#================================
vl.append(variable())
vl[-1].keyword  ='pt3_a_diagrams'
vl[-1].family   ='post'
vl[-1].default  ='yes'
vl[-1].datatype ='characters'
vl[-1].comment  ='Switch whether to calculate the A diagrams family in PT3. A diagrams are the self-consistent diagrams (PT2 inclusions in the Green\'s function). \
                Valid choices include: \'yes\', \'no\', or \'only\'.'

#================================
vl.append(variable())
vl[-1].keyword  ='pt_density_matrix'
vl[-1].family   ='post'
vl[-1].default  ='no'
vl[-1].datatype ='characters'
vl[-1].comment  ='Triggers the calculation of a correlated density matrix within MBPT. \
                Valid choices include: \'no\', \'PT2\', \'ONE-RING\', or \'GW\'.'

#================================
vl.append(variable())
vl[-1].keyword  ='use_correlated_density_matrix'
vl[-1].family   ='post'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Chooses to use another density matrix for the Fock hamiltonian to be employed in self-energy calculations. \
                 Used in conjonction with \'pt_density_matrix\' or with \'read_fchk\' or read an existing DENSITY_MATRIX file.'

#================================
vl.append(variable())
vl[-1].keyword  ='read_restart'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Read the RESTART file and restart from it.'

#================================
vl.append(variable())
vl[-1].keyword  ='ignore_bigrestart'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Considers a big RESTART as if it was a small RESTART.'

#================================
vl.append(variable())
vl[-1].keyword  ='force_energy_qp'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Force the reading of the ENERGY_QP file whatever the postscf choice.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_hartree'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the Hartree potential and exchange expectation value on eigenstates.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_density_matrix'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the density matrix in the DENSITY_MATRIX file'

#================================
vl.append(variable())
vl[-1].keyword  ='print_rho_grid'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the electronic density discretized on the DFT grid into a file \'rho_grid.dat\'. The density is calculated from the DENSITY_MATRIX file or from a Gaussian file using \'read_fchk\'.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_eri'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Dumps the Electron Repulsion Integral on a file.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_wfn'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints some wavefunctions along some selected lines.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_cube'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints some wavefunctions in a 3D volumetric file with cube format'

#================================
vl.append(variable())
vl[-1].keyword  ='read_fchk'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='characters'
vl[-1].comment  ='Triggers the reading of an external Gaussian formatted checkpoint file (named gaussian.fchk) that contains density matrices. \
Basis sets have to be precisely the same in MOLGW and in Gaussian, which requires a manual input of the basis set in both codes. \
Options are \'no\' (no reading), \'SCF\' (for self-consistent field), \'CC\' (for coupled-cluster), or \'MP2\' (for MP2). \
Today, only works for Cartesian Gaussian and for spin restricted calculations.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_w'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Dumps the spectral function of the screened Coulomb W. This is necessary for a subsequent BSE run.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_restart'
vl[-1].family   ='io'
vl[-1].default  ='yes'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints a small RESTART file at each SCF cycle. \
There are two kinds of RESTART files: the small RESTART and the big RESTART. \
The former contains only the information about the occupied wavefunctions. \
This is a very small file and the writing should not hit too much on performance.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_bigrestart'
vl[-1].family   ='io'
vl[-1].default  ='yes'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the big RESTART file at the end of the SCF loop. \
There are two kinds of RESTART files: the small RESTART and the big RESTART. \
The latter is written only when self-consistency has been reached. \
It contains all the states and the Hamiltonian and allows one to completely skip the scf loop \
or to start over with another basis set.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_sigma'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the value of the GW self-energy on the sampling frequencies in files.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_pdos'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the Mulliken weight of each eigenvector on a given atom or a given series of atoms.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_spatial_extension'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the wavefunction extension calculated as &lt;r<sup>2</sup>&gt; - &lt;r&gt;<sup>2</sup>'

#================================
vl.append(variable())
vl[-1].keyword  ='print_multipole'
vl[-1].family   ='io'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the electric multipole expansion for the electronic density and the nuclei.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_yaml'
vl[-1].family   ='io'
vl[-1].default  ='yes'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Creates an output file in YAML format. Easier to read for python post-processing.'

#================================
vl.append(variable())
vl[-1].keyword  ='length_unit'
vl[-1].family   ='system'
vl[-1].default  ='angstrom'
vl[-1].datatype ='characters'
vl[-1].comment  ='Chooses the units of the atomic coordinates. Can be \'angstrom\' or \'bohr\'. \
Could be abbreviated in \'A\' or \'au\'.'

#================================
vl.append(variable())
vl[-1].keyword  ='natom'
vl[-1].family   ='system'
vl[-1].default  =0
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of atoms in the molecule. This is the number of lines to be read in the following section of the input file if no xyz file is provided.'

#================================
vl.append(variable())
vl[-1].keyword  ='xyz_file'
vl[-1].family   ='system'
vl[-1].default  =''
vl[-1].datatype ='characters'
vl[-1].comment  ='Specifies the location of the xyz file that contains the atomic positions. It can be used as an alternate route to set atomic coordinate.'

#================================
vl.append(variable())
vl[-1].keyword  ='nghost'
vl[-1].family   ='system'
vl[-1].default  =0
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of ghost atoms in the molecule. Used to place basis function where there is no atom. Useful for Basis Set Superposition Error'

#================================
vl.append(variable())
vl[-1].keyword  ='eta'
vl[-1].family   ='post'
vl[-1].default  = 0.001
vl[-1].datatype ='real'
vl[-1].comment  ='Is a the tiny imaginary part used in the denominator of the Green\'s function to shift the pole off the axis, so to avoid divergences.\
This is an energy in Hartree. \
It should be set to the lowest value possible in theory. However, in practice, a too low value of eta would induce huge and unstable GW corrections. \
The default value is usually very accurate and there is no need to use a lower value. But for states apart from the band gap, a large value of eta may be beneficial \
for stability. eta=0.01 is already much more stable. Note that for QSGW increasing eta is most often unavoidable.'

#================================
vl.append(variable())
vl[-1].keyword  ='scissor'
vl[-1].family   ='post'
vl[-1].default  = 0.
vl[-1].datatype ='real'
vl[-1].comment  ='Sets a rigid energy shift of the unoccupied states, so to mimick a GW calculation without actually doing it.'

#================================
vl.append(variable())
vl[-1].keyword  ='grid_memory'
vl[-1].family   ='hardware'
vl[-1].default  = 400.0
vl[-1].datatype ='real'
vl[-1].comment  ='Sets the maximum memory usage in Mb allowed to store the wavefunctions on the quadrature points for XC integrals.'

#================================
vl.append(variable())
vl[-1].keyword  ='scalapack_block_min'
vl[-1].family   ='hardware'
vl[-1].default  = 1000
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the minimum block size to distribute a non-distributed matrix with SCALAPACK. \
If scalapack_block_min=400, then a 900x900 matrix will be distributed on a 2x2 processor grid. \
If scalapack_block_min=500, then a 900x900 matrix will no be distributed.'

#================================
vl.append(variable())
vl[-1].keyword  ='eri3_nprow'
vl[-1].family   ='hardware'
vl[-1].default  = 1
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets number of row processors for the distribution of the 3-center integrals.  \
eri3_nprow X eri3_npcol must be equal to the number of MPI threads else MOLGW decides on its own.'

#===============================
vl.append(variable())
vl[-1].keyword  ='eri3_npcol'
vl[-1].family   ='hardware'
vl[-1].default  = 1
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets number of column processors for the distribution of the 3-center integrals.  \
eri3_nprow X eri3_npcol must be equal to the number of MPI threads else MOLGW decides on its own.'

#===============================
vl.append(variable())
vl[-1].keyword  ='eri3_nbatch'
vl[-1].family   ='hardware'
vl[-1].default  = 1
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of batches when calculating the 3-center integrals. \
Having a large eri3_nbatch reduces the memory foot print, however it may lower the performance.'

#================================
vl.append(variable())
vl[-1].keyword  ='mpi_nproc_ortho'
vl[-1].family   ='hardware'
vl[-1].default  =1
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of processors left to parallelize on other directions. The main direction (auxiliary basis or DFT grid points) is obtained by \
<b>mpi_nproc</b> / <b>mpi_nproc_ortho</b>, which must be an integer.'

#================================
vl.append(variable())
vl[-1].keyword  ='dft_core'
vl[-1].family   ='post'
vl[-1].default  =0
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of states considered as core in &lt;&Sigma;<sub>x</sub>-<i>v</i><sub>xc</sub>&gt. This options is meant to mimic the pseudopotential approximation.'

#================================
vl.append(variable())
vl[-1].keyword  ='virtual_fno'
vl[-1].family   ='post'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Activates the Frozen Natural Orbitals technique to span the virtual orbitals subspace with fewer orbitals. \
The dimension of the space is set up with the input keyword nvirtualg or nvirtualw. \
Actually the virtual orbital space is determined by the minimum MIN(nvirtualg,nvirtualw).'

#================================
vl.append(variable())
vl[-1].keyword  ='rcut_mbpt'
vl[-1].family   ='post'
vl[-1].default  ='1.0'
vl[-1].datatype ='real'
vl[-1].experimental  ='yes'
vl[-1].comment  ='EXPERIMENTAL'

#================================
vl.append(variable())
vl[-1].keyword  ='gwgamma_tddft'
vl[-1].family   ='post'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].experimental  ='yes'
vl[-1].comment  ='EXPERIMENTAL. Calculates the vertex using the DFT flavor specified in the ground-state calculation.'

#================================
vl.append(variable())
vl[-1].keyword  ='ecp_type'
vl[-1].family   ='general'
vl[-1].default  =''
vl[-1].datatype ='characters'
vl[-1].comment  ='Name of the Effective Core Potential. For instance, Gold using the cc-pVDZ-PP basis set should have ecp_type=\'PP\', \
so that MOLGW looks for the file Au_PP in the basis_path folder.'

#================================
vl.append(variable())
vl[-1].keyword  ='ecp_elements'
vl[-1].family   ='general'
vl[-1].default  =''
vl[-1].datatype ='characters'
vl[-1].comment  ='Contains the list of elements (separated by spaces) that should be treated with an Effective Core Potential.'

#================================
vl.append(variable())
vl[-1].keyword  ='ecp_quality'
vl[-1].family   ='general'
vl[-1].default  ='high'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the number of grid points use to evaluate the Effective Core Potential integrals in real space. \
Possible values are \'low\', \'medium\', \'high\', \'very high\', \'insane\'. \
It could be abbreviated in \'l\', \'m\', \'h\', \'vh\', \'i\'. \
\'high\' is usually fine. \'insane\' is only meant for debugging since it is overdoing a lot.'

#================================
vl.append(variable())
vl[-1].keyword  ='ecp_basis'
vl[-1].family   ='general'
vl[-1].default  =''
vl[-1].datatype ='characters'
vl[-1].comment  ='Name of the basis set to be used for elements specified in list ecp_elements.'

#================================
vl.append(variable())
vl[-1].keyword  ='ecp_auxil_basis'
vl[-1].family   ='general'
vl[-1].default  =''
vl[-1].datatype ='characters'
vl[-1].comment  ='Name of the auxiliary basis set to be used for elements specified in list ecp_elements.'

#================================
vl.append(variable())
vl[-1].keyword  ='time_step'
vl[-1].family   ='rt_tddft'
vl[-1].default  =1.
vl[-1].datatype ='real'
vl[-1].comment  ='Time step for real-time dynamics in atomic units.'

#================================
vl.append(variable())
vl[-1].keyword  ='time_sim'
vl[-1].family   ='rt_tddft'
vl[-1].default  =10.
vl[-1].datatype ='real'
vl[-1].comment  ='Duration of a real-time dynamics in atomic units.'

#================================
vl.append(variable())
vl[-1].keyword  ='prop_type'
vl[-1].family   ='rt_tddft'
vl[-1].default  ='MAG2'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the type of propagation algorithm in the real-time dynamics. \
                 \'CN\' stands for Crank-Nickolson \
                 \'MAG2\' stands for Magnus 2nd order.'
#================================
vl.append(variable())
vl[-1].keyword  ='ci_greens_function'
vl[-1].family   ='post'
vl[-1].default  ='holes'
vl[-1].datatype ='characters'
vl[-1].experimental  ='yes'
vl[-1].comment  ='EXPERIMENTAL. Selects which part of the Green\'s function is to be calculated: holes, electrons, or both.'

#================================
vl.append(variable())
vl[-1].keyword  ='excit_name'
vl[-1].family   ='rt_tddft'
vl[-1].default  ='NO'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the type of excitation of a system in the real-time dynamics. \
                         \'GAU stands for a linearly polarized uniform Gaussian electric field\''

#================================
vl.append(variable())
vl[-1].keyword  ='ci_type'
vl[-1].family   ='post'
vl[-1].default  ='all'
vl[-1].datatype ='characters'
vl[-1].comment  ='Selects which excitations will be included in the CI expansion. Valid choices are \'all\', \'CISD\', \'CISDT\', \'CISDTQ\'.'

# For excitation field parameters follow notations from article Lopata et al. Modeling Fast Electron ... J. Chem Theory Comput, 2011
#================================
vl.append(variable())
vl[-1].keyword  ='excit_kappa'
vl[-1].family   ='rt_tddft'
vl[-1].default  =2.e-5
vl[-1].datatype ='real'
vl[-1].comment  ='Maximum Gaussian excitation field strength in atomic units.'

#================================
vl.append(variable())
vl[-1].keyword  ='ci_nstate'
vl[-1].family   ='post'
vl[-1].default  = 1
vl[-1].datatype ='integer'
vl[-1].comment  ='Selects how many CI states should be calculated in the diagonalization. If ci_nstate is lower than the number of configuration, \
 a Davidson partial diagonalization is performed, else a full (SCA)LAPACK diagonalization is triggered.'

#================================
vl.append(variable())
vl[-1].keyword  ='excit_omega'
vl[-1].family   ='rt_tddft'
vl[-1].default  =0.2
vl[-1].datatype ='real'
vl[-1].comment  ='The excitation pulse width in atomic units for the real-time dynamics.'

#================================
vl.append(variable())
vl[-1].keyword  ='ci_nstate_self'
vl[-1].family   ='post'
vl[-1].default  = 1
vl[-1].datatype ='integer'
vl[-1].comment  ='Selects how many CI states in the N+1 or N-1 electron calculations. If ci_nstate_self is lower than the number of configuration, \
 a Davidson partial diagonalization is performed, else a full (SCA)LAPACK diagonalization is triggered.'

#================================
vl.append(variable())
vl[-1].keyword  ='excit_time0'
vl[-1].family   ='rt_tddft'
vl[-1].default  =3.
vl[-1].datatype ='real'
vl[-1].comment  ='Center of the excitation pulse in atomic units for the real-time dynamics.'

#================================
vl.append(variable())
vl[-1].keyword  ='ci_spin_multiplicity'
vl[-1].family   ='post'
vl[-1].default  = 1
vl[-1].datatype ='integer'
vl[-1].comment  ='Spin multiplicity in CI calculations.'

#================================
vl.append(variable())
vl[-1].keyword  ='excit_dir'
vl[-1].family   ='rt_tddft'
vl[-1].default  =( 1.0 , 0.0, 0.0 )
vl[-1].datatype ='vector_1d_3'
vl[-1].comment  ='Excitation direction for the real-time dynamics.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_tddft_matrices'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints some matrices of the real-time dynamics into the file check_matrix.dat.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_cube_rho_tddft'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints electronic density in a 3D volumetric file with cube format for each simulation step in the real-time dynamics'

#================================
vl.append(variable())
vl[-1].keyword  ='print_line_rho_tddft'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints electronic density along a line, which parameters must be provided in manual_plot_rho_tddft file.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_line_rho_diff_tddft'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints electronic density difference along a line, which parameters must be provided in manual_plot_rho_tddft file.'

#================================
vl.append(variable())
vl[-1].keyword  ='write_step'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='1'
vl[-1].datatype ='real'
vl[-1].comment  ='Determines the time step for data recording in the real-time dynamics'

#================================
vl.append(variable())
vl[-1].keyword  ='pred_corr'
vl[-1].family   ='rt_tddft'
vl[-1].default  ='PC1'
vl[-1].datatype ='characters'
vl[-1].comment  ='Sets the predictor-corrector scheme in the real-time dynamics.'

#================================
vl.append(variable())
vl[-1].keyword  ='n_hist'
vl[-1].family   ='rt_tddft'
vl[-1].default  =2
vl[-1].datatype ='integer'
vl[-1].comment  ='Number of memorised previous hamiltonian values for its extrapolation in the real-time dynamics. n_hist=1 means that H(t_i+1)=H(t_i); n_hist=2 : H(t_i+1)=a*H(t_i)+b*(t_i-1); etc.'

#================================
vl.append(variable())
vl[-1].keyword  ='n_iter'
vl[-1].family   ='rt_tddft'
vl[-1].default  ='2'
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of iterations for the PC7 in the real-time dynamics'

#================================
vl.append(variable())
vl[-1].keyword  ='calc_spectrum'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Calculates absorption spectrum in the real-time dynamics'

#================================
vl.append(variable())
vl[-1].keyword  ='read_tddft_restart'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Ignore the RESTART_TDDFT file.'

#================================
vl.append(variable())
vl[-1].keyword  ='print_tddft_restart'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='yes'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints a RESTART_TDDFT file which contains wavefunction coefficients for the last time moment of a simulation.'

#================================
vl.append(variable())
vl[-1].keyword  ='vel_projectile'
vl[-1].family   ='rt_tddft'
vl[-1].default  =( 0.0 , 0.0, 1.0 )
vl[-1].datatype ='vector_1d_3'
vl[-1].comment  ='Projectile initial velocity.'

#================================
vl.append(variable())
vl[-1].keyword  ='n_restart_tddft'
vl[-1].family   ='rt_tddft'
vl[-1].default  = 50
vl[-1].datatype ='integer'
vl[-1].comment  ='RESTART_TDDFT file will be written during simulation each n_retart_tddft iteration (provided that print_tddft_restart is yes)'

#================================
vl.append(variable())
vl[-1].keyword  ='print_cube_diff_tddft'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the difference of electronic density with respect to initial density in a 3D volumetric file with cube format for each simulation step in the real-time dynamics'

#================================
vl.append(variable())
vl[-1].keyword  ='print_dens_traj_tddft'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the electronic density along the projectile trajectory for several impact parameters in the real-time dynamics'

#================================
vl.append(variable())
vl[-1].keyword  ='calc_q_matrix'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Calculate and print q_matrix which is the projection of a propagated state psi(t) onto the initial state psi(0) in the real-time dynamics'

#================================
vl.append(variable())
vl[-1].keyword  ='print_dens_traj'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the electronic density along the projectile trajectory for several impact parameters using real wave function'

#================================
vl.append(variable())
vl[-1].keyword  ='print_dens_traj_points_set'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Prints the electronic density between pairs of points given in manual_dens_points_set file.'

#================================
vl.append(variable())
vl[-1].keyword  ='calc_dens_disc'
vl[-1].family   ='io_rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Calculate electronic density in the discs during the real-time dynamics'

#================================
vl.append(variable())
vl[-1].keyword  ='r_disc'
vl[-1].family   ='rt_tddft'
vl[-1].default  = 200.0
vl[-1].datatype ='real'
vl[-1].comment  ='Radius of the disc for denisty calculations (option calc_dens_disc) for the real-time dynamics.'

#================================
vl.append(variable())
vl[-1].keyword  ='tddft_frozencore'
vl[-1].family   ='rt_tddft'
vl[-1].default  ='no'
vl[-1].datatype ='yes/no'
vl[-1].comment  ='Do not "propagate" states mentioned in the manual_tddft_frozencore file in the real-time dynamics.'

#================================
vl.append(variable())
vl[-1].keyword  ='ncore_tddft'
vl[-1].family   ='rt_tddft'
vl[-1].default  = 0
vl[-1].datatype ='integer'
vl[-1].comment  ='Sets the number of frozen core states in the real-time dynamics.'




#============================================================================
# Sort by alphabetical order
#============================================================================
vl.sort(key=lambda x: x.keyword)


#============================================================================
#            Fortran output: input variable declaration
#============================================================================

object_dict = dict()

for v in vl:
    object_dict[v.keyword]=dict()
    object_dict[v.keyword]['family']   = v.family
    object_dict[v.keyword]['default']  = v.default
    object_dict[v.keyword]['datatype'] = v.datatype
    object_dict[v.keyword]['comment']  = v.comment

with open('input_variables.yaml', 'w') as outfile:
    dump(object_dict, outfile, Dumper=Dumper)

#============================================================================
#            Fortran output: input variable declaration
#============================================================================

script_file_name = os.path.basename(__file__)
header = '!======================================================================\n'            + \
         '! The following lines have been generated by a python script: '+script_file_name+'\n' + \
         '! Do not alter them directly: they will be overriden sooner or later by the script\n' + \
         '! To add a new input variable, modify the script directly\n'                          + \
         '! Generated by '+script_file_name+' on '+today+'\n'                                   + \
         '!======================================================================\n\n'

print("Set up file: ../src/input_variable_declaration.f90")
ffor = open('../src/input_variable_declaration.f90','w')

ffor.write(header)

for v in vl:
  # Exclude a few input variable due to name clash
  if v.keyword in ['basis','auxil_basis','natom','nghost','read_restart']:
    continue

  if v.datatype =='integer':
    ffor.write(' integer,protected :: ' + v.keyword +'\n')
  elif v.datatype =='real':
    ffor.write(' real(dp),protected :: ' + v.keyword +'\n')
  elif v.datatype =='vector_1d_3':
    ffor.write(' real(dp),protected :: ' + v.keyword +'(3)\n')
  elif v.datatype =='yes/no':
    ffor.write(' character(len=3),protected :: ' + v.keyword +'\n')
  elif v.datatype =='characters':
    ffor.write(' character(len=140),protected :: ' + v.keyword +'\n')


ffor.write('\n\n!======================================================================\n')
ffor.close()

#============================================================================
#            Fortran output: input variable declaration, input variable namelist, and their default value
#============================================================================

print("Set up file: ../src/input_variables.f90")
ffor = open('../src/input_variables.f90','w')

ffor.write(header)

ffor.write(' namelist /molgw/   &\n')
for i in range(len(vl)-1):
  ffor.write('    '+vl[i].keyword+',       &\n')
ffor.write('    '+vl[-1].keyword+'\n\n')

ffor.write('!=====\n\n')


for v in vl:
  if v.datatype =='integer':
    ffor.write(' '+v.keyword+'='+str(v.default)+'\n')
  elif v.datatype =='real':
    ffor.write(' '+v.keyword+'='+str(v.default)+'_dp \n')
  elif v.datatype =='vector_1d_3':
     x,y,z=str(v.default).strip("()").split(',')
     ffor.write(' '+v.keyword+'='+'(/ '+x+'_dp ,'+y+'_dp ,'+z+'_dp'+' /)'+'\n')
  elif v.datatype =='yes/no' or v.datatype =='characters':
    ffor.write(' '+v.keyword+'=\''+str(v.default)+'\'\n')


ffor.write('\n\n!======================================================================\n')
ffor.close()

#============================================================================
#            Fortran output: Echoing of all the input variable values
#============================================================================


print("Set up file: ../src/echo_input_variables.f90")
ffor = open('../src/echo_input_variables.f90','w')

ffor.write(header)

for v in vl:
  if 'real' in v.datatype:
    fortran_format = '\'(1x,a24,2x,es16.8)\''
  elif 'integer' in v.datatype:
    fortran_format = '\'(1x,a24,2x,i8)\''
  elif 'characters' in v.datatype:
    fortran_format = '\'(1x,a24,6x,a)\''
  elif 'yes' in v.datatype:
    fortran_format = '\'(1x,a24,6x,a)\''
  elif 'vector_1d_3' in v.datatype:
    fortran_format = '\'(1x,a24,2x,"(",3(es16.8,2x),")")\''
  else:
    fortran_format = 'ERROR'
  ffor.write(' write(stdout,'+fortran_format+') \''+v.keyword+'\','+v.keyword+' \n')



ffor.write('\n\n!======================================================================\n')
ffor.close()

#============================================================================
#            Fortran output: Echoing of all the input variable values in YAML format
#============================================================================

print("Set up file: ../src/echo_input_variables_yaml.f90")
ffor = open('../src/echo_input_variables_yaml.f90','w')

ffor.write(header)

for v in vl:
  keyvalue = v.keyword
  right_spaces = str( 30 - len(v.keyword) )
  if 'real' in v.datatype:
    fortran_format = '\'(4x,a,' + right_spaces + 'x,es16.8)\''
  elif 'integer' in v.datatype:
    fortran_format = '\'(4x,a,' + right_spaces + 'x,i8)\''
  elif 'characters' in v.datatype:
    keyvalue = 'TRIM(' + keyvalue + ')'
    fortran_format = '\'(4x,a,' + right_spaces + 'x,a)\''
  elif 'yes' in v.datatype:
    if 'y' in v.keyword.lower():
        keyvalue = "'"+str(True)+"'"
    else:
        keyvalue = "'"+str(False)+"'"
    fortran_format = '\'(4x,a,' + right_spaces + 'x,a)\''
  elif 'vector_1d_3' in v.datatype:
    fortran_format = '\'(4x,a,' + right_spaces + 'x,"[",es16.8,", ",es16.8,", ",es16.8,"]")\''
  else:
    fortran_format = 'ERROR'
  ffor.write(' write(unit_yaml,'+fortran_format+') \''+v.keyword+':\','+keyvalue+' \n')



ffor.write('\n\n!======================================================================\n')
ffor.close()


#============================================================================
#            HTML output
#============================================================================
print("Set up file: ../docs/input_variables.html")
fhtml = open('../docs/input_variables.html','w')

fhtml.write('<html>\n')
fhtml.write('<head>\n')
fhtml.write('<link rel="stylesheet" type="text/css" href="molgw.css">\n')
fhtml.write('</head>\n')

fhtml.write('<body>\n')
fhtml.write('<a name=top>\n')
fhtml.write('<h1>Input variable list</h1>\n')
fhtml.write('<hr>\n<br>\n')

# Mandatory
fhtml.write('<h3>Mandatory input variables</h3>\n<p>\n')
for v in vl:
  if v.mandatory =='yes':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')

# System
fhtml.write('<h3>Physical system setup input variables</h3>\n<p>\n')
for v in vl:
  if v.family =='system':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')

# General
fhtml.write('<h3>General input variables</h3>\n<p>\n')
for v in vl:
  if v.family =='general':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')

# SCF
fhtml.write('<h3>Self-consistency input variables</h3>\n<p>\n')
for v in vl:
  if v.family =='scf':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')

# Post 
fhtml.write('<h3>Correlation and excited states post-treatment input variables</h3>\n<p>\n')
for v in vl:
  if v.family =='post':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')

# IO family
fhtml.write('<h3>IO input variables</h3>\n<p>\n')
for v in vl:
  if v.family =='io':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')

# Parallelization family
fhtml.write('<h3>Hardware input variables</h3>\n<p>\n')
for v in vl:
  if v.family =='hardware':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')

# Real-time TDDFT
fhtml.write('<h3>Real time TDDFT input variables</h3>\n<p>\n')
for v in vl:
  if v.family =='rt_tddft':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')

# IO Real-time TDDFT
fhtml.write('<h3>IO Real time TDDFT input variables</h3>\n<p>\n')
for v in vl:
  if v.family =='io_rt_tddft':
    fhtml.write('<a href=#'+v.keyword+'>'+v.keyword+'</a> ')



fhtml.write('<br><br><br><hr>\n')

# Start the complete list
fhtml.write('<br><br><br>\n')
fhtml.write('<h2>Complete list of input variables</h2>\n')
fhtml.write('<br><br>\n<ul>\n')
for v in vl:
  v.printhtml(fhtml)
fhtml.write('</ul>\n')
fhtml.write('<br><br><br><br><br><br><br><br>\n')
fhtml.write('<a href=#top>Back to the top of the page</a> ')
fhtml.write('<div style="float: right"><a href=molgw_manual.html>Back to the manual</a></div>')
fhtml.write('<br><br>')
fhtml.write('<i>Generated by '+script_file_name+' on '+today+'</i>')
fhtml.write('<br><br>')
fhtml.write('</body>\n')
fhtml.write('</html>\n')

fhtml.close()

print("Done!")


