#!/usr/bin/python3
##################################################
#
# This file is part of MOLGW
# Author: Fabien Bruneval
#
# This python submodule provides useful functions to prepare/analyze ml23 runs
#
#
##################################################

import numpy as np
import json, yaml
import pathlib
import matplotlib.pyplot as plt
from molgw import __version__, Molecule

default_input_parameters = {
      'scf': 'hf',
      'postscf': 'g0w0',
      'basis': 'aug-cc-pVQZ',
      'auxil_basis': 'pauto',
      'selfenergy_state_range': 3,
      'selfenergy_state_range': 3,
              }


# Number of frozen MO: none
frozencore = {
    'H': 0,
    'Li': 0,
    'Be': 0,
    'B': 0,
    'C': 0,
    'N': 0,
    'O': 0,
    'F': 0,
    'Si': 0,
    'S': 0,
    'P': 0,
    'Cl': 0,
    'Ne': 0,
    'Ar': 0,
}

def get_frozencore(cas):
    structure_string = structures[cas]
    lines = structure_string.split('\n')
    natom = int(lines[0])
    frozen_orbitals = 0
    for iatom in range(natom):
        element = lines[iatom+2].split()[0]
        frozen_orbitals += frozencore[element]
    return frozen_orbitals

structures = {
 'co': '2\n\nC  0.0000  0.0000  0.0000\nO  0.0000  0.0000  1.1335 \n', 'cs': '2\n\nC  0.0000  0.0000  0.0000\nS  0.0000  0.0000  1.5485 \n', 'licl': '2\n\nLi 0.0000  0.0000  0.0000\nCl  0.0000  0.0000  2.0265\n', 'c2': '2\nCarbon dimer,^1\\Sigma_g^+,CC3,aug-cc-pVTZ\nC    0.00000000            0.00000000            0.62402126\nC    0.00000000            0.00000000           -0.62402126\n', 'sih4': '5\n\nSi 0.0000  0.0000  0.0000\nH  1.4807  0.0000  0.0000\nH -0.4936  1.3960  0.0000\nH -0.4936 -0.6980 -1.2090\nH -0.4936 -0.6980  1.2090\n', 'ph3': '4\n\nP  0.0000  0.0000  0.0000\nH  0.7746  1.1826  0.0000\nH  0.7746 -0.5913  1.0242\nH  0.7746 -0.5913 -1.0242\n', 'h2o': '3\n\nO  0.0000  0.0000  0.0000\nH  0.9591  0.0000  0.0000\nH -0.2373  0.9293  0.0000\n', 'n2': '2\n\nN  0.0000  0.0000  0.0000\nN  0.0000  0.0000  1.1007\n', 'h2s': '3\n\nS    0.00000000            0.00000000           -0.26652056\nH    0.00000000            0.96219289            0.66259489\nH    0.00000000           -0.96219289            0.66259489\n', 'bn': '2\n\nB  0.0000  0.0000  0.0000\nN  0.0000  0.0000  1.2765\n', 'ch4': '5\n\nC  0.0000   0.0000  0.0000\nH  1.0879   0.0000  0.0000\nH -0.3626   1.0257  0.0000\nH -0.3626  -0.5128 -0.8883\nH -0.3626  -0.5128  0.8883\n', 'ch2o': '4\n\nC  0.0000  0.0000 -0.6030\nO  0.0000  0.0000  0.6054\nH  0.0000  0.9347 -1.1822\nH  0.0000 -0.9347 -1.1822\n', 'co2': '3\n\nO -1.1652  0.0000  0.0000\nC  0.0000  0.0000  0.0000\nO  1.1652  0.0000  0.0000\n', 'hf': '2\n\nH  0.0000  0.0000  0.0000\nF  0.0000  0.0000  0.9196\n', 'f2': '2\n\nF  0.0000  0.0000  0.0000\nF  0.0000  0.0000  1.4137\n', 'ar': '1\n\nAr  0.0000  0.0000  0.0000\n', 'bh3': '4\n\nH  0.0000  0.0000   0.0000\nB  1.1848  0.0000   0.0000\nH  1.7772  1.0260   0.0000\nH  1.7772 -1.0260   0.0000\n', 'hcl': '2\n\nH  0.0000  0.0000  0.0000\nCl 0.0000  0.0000  1.2751\n', 'nh3': '4\n\nN  0.0000  0.0000  0.0000\nH  0.3816  0.9375  0.0000\nH  0.3816 -0.4687  0.8119\nH  0.3816 -0.4687 -0.8119\n', 'beo': '2\n\nBe 0.0000  0.0000  0.0000\nO  0.0000  0.0000  1.3621\n', 'lif': '2\n\nLi 0.0000  0.0000  0.0000\nF  0.0000  0.0000  1.5783\n', 'bf': '2\n\nB  0.0000  0.0000  0.0000\nF  0.0000  0.0000  1.2685\n', 'ne': '1\n\nNe  0.0000  0.0000  0.0000\n'}

names = [ k for k in structures ]

#
# ml23.molecule: dictionary of molgw.Molecule objects
#
molecule = dict()
for k, v in structures.items():
    molecule[k] = Molecule(v)


############################################
# Results from Marie-Loos
# https://doi.org/10.1021/acs.jctc.4c00216
############################################
#  aug-cc-pVTZ

g0w0_avt_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "QUACK",
    "basis_size": "3",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVQZ",
    "calc_type": "G0W0@HF",
    "data": {'ch4': -14.753, 'h2s': -10.5, 'bn': -11.752, 'licl': -9.9844, 'c2': -12.928, 'sih4': -13.214, 'co2': -14.221, 'ch2o': -11.38, 'beo': -9.7273, 'lif': -11.384, 'f2': -16.334, 'hf': -16.237, 'co': -14.777, 'n2': -16.35, 'bf': -11.325, 'ne': -21.432, 'ph3': -10.787, 'h2o': -12.884, 'hcl': -12.778, 'cs': -12.378, 'ar': -15.711, 'bh3': -13.678, 'nh3': -11.201}
}

eomccsd_atq_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "CFOUR",
    "basis_size": "3",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVTZ",
    "calc_type": "eom-CCSD",
    "data": {'ch4': -14.387, 'h2s': -10.421, 'bn': -11.971, 'licl': -9.9564, 'c2': -12.978, 'sih4': -12.844, 'co2': -13.787, 'ch2o': -10.848, 'beo': -9.875, 'lif': -11.398, 'f2': -15.616, 'hf': -16.021, 'co': -14.19, 'n2': -15.641, 'bf': -11.25, 'ne': -21.326, 'ph3': -10.623, 'h2o': -12.594, 'hcl': -12.712, 'cs': -11.553, 'ar': -15.672, 'bh3': -13.342, 'nh3': -10.862}
}

eomccsdt_atq_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "CFOUR",
    "basis_size": "3",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVTZ",
    "calc_type": "eom-CCSDT",
    "data": {'ch4': -14.365, 'h2s': -10.39, 'bn': -11.98, 'licl': -9.8825, 'c2': -12.54, 'sih4': -12.794, 'co2': -13.716, 'ch2o': -10.84, 'beo': -9.8642, 'lif': -11.379, 'f2': -15.688, 'hf': -16.077, 'co': -13.952, 'n2': -15.517, 'bf': -11.157, 'ne': -21.398, 'ph3': -10.599, 'h2o': -12.629, 'hcl': -12.667, 'cs': -11.346, 'ar': -15.606, 'bh3': -13.304, 'nh3': -10.876}
}

eomccsdtq_atq_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "CFOUR",
    "basis_size": "3",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVTZ",
    "calc_type": "eom-CCSDTQ",
    "data": {'ch4': -14.376, 'h2s': -10.393, 'bn': -11.966, 'licl': -9.8964, 'c2': -12.471, 'sih4': -12.795, 'co2': -13.734, 'ch2o': -10.887, 'beo': -9.9385, 'lif': -11.453, 'f2': -15.725, 'hf': -16.14, 'co': -13.935, 'n2': -15.487, 'bf': -11.15, 'ne': -21.455, 'ph3': -10.6, 'h2o': -12.673, 'hcl': -12.676, 'cs': -11.316, 'ar': -15.614, 'bh3': -13.307, 'nh3': -10.899}
}

selectedci_atq_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "QUANTUM PACKAGE",
    "basis_size": "3",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVTZ",
    "calc_type": "selected CI",
    "data": {'ch4': -14.377, 'h2s': -10.393, 'bn': -11.98, 'licl': -9.897, 'c2': -12.463, 'sih4': -12.793, 'co2': -13.733, 'ch2o': -10.899, 'beo': -9.972, 'lif': -11.468, 'f2': -15.729, 'hf': -16.149, 'co': -13.925, 'n2': -15.486, 'bf': -11.149, 'ne': -21.461, 'ph3': -10.596, 'h2o': -12.679, 'hcl': -12.676, 'cs': -11.3, 'ar': -15.613, 'bh3': -13.307, 'nh3': -10.901}
}

#  aug-cc-pVQZ
g0w0_avq_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "QUACK",
    "basis_size": "4",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVQZ",
    "calc_type": "G0W0@HF",
    "data": {'ch4': -14.872, 'h2s': -10.66, 'bn': -11.907, 'licl': -10.18, 'c2': -13.082, 'sih4': -13.312, 'co2': -14.425, 'ch2o': -11.556, 'beo': -9.9066, 'lif': -11.594, 'f2': -16.559, 'hf': -16.453, 'co': -14.915, 'n2': -16.519, 'bf': -11.42, 'ne': -21.655, 'ph3': -10.911, 'h2o': -13.08, 'hcl': -12.967, 'cs': -12.523, 'ar': -15.926, 'bh3': -13.78, 'nh3': -11.362}
}

eomccsd_avq_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "CFOUR",
    "basis_size": "4",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVQZ",
    "calc_type": "eom-CCSD",
    "data": {'ch4': -14.428, 'h2s': -10.49, 'bn': -12.018, 'licl': -10.067, 'c2': -13.03, 'sih4': -12.877, 'co2': -13.879, 'ch2o': -10.924, 'beo': -9.941, 'lif': -11.493, 'f2': -15.722, 'hf': -16.117, 'co': -14.246, 'n2': -15.709, 'bf': -11.279, 'ne': -21.432, 'ph3': -10.662, 'h2o': -12.675, 'hcl': -12.812, 'cs': -11.609, 'ar': -15.802, 'bh3': -13.375, 'nh3': -10.923}
}

eomccsdt_avq_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "CFOUR",
    "basis_size": "4",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVQZ",
    "calc_type": "selected CI",
    "data": {'ch4': -14.395, 'h2s': -10.455, 'bn': -12.019, 'licl': -9.993, 'c2': -12.585, 'sih4': -12.82, 'co2': -13.794, 'ch2o': -10.898, 'beo': -9.9158, 'lif': -11.452, 'f2': -15.767, 'hf': -16.145, 'co': -14.005, 'n2': -15.574, 'bf': -11.185, 'ne': -21.473, 'ph3': -10.634, 'h2o': -12.689, 'hcl': -12.764, 'cs': -11.404, 'ar': -15.735, 'bh3': -13.33, 'nh3': -10.922}
}

selectedci_avq_homo = {
    "orbital": "HOMO",
    "remark": "",
    "code": "QUANTUM PACKAGE",
    "basis_size": "4",
    "basis": "gaussian",
    "code_version": "",
    "qpe": "",
    "DOI": "https://doi.org/10.1021/acs.jctc.4c00216",
    "basis_name": "aug-cc-pVQZ",
    "calc_type": "selected CI",
    "data": {'ch4': -14.407, 'h2s': -10.456, 'bn': -12.019, 'licl': -10.007, 'c2': -12.497, 'sih4': -12.818, 'co2': -13.823, 'ch2o': -10.954, 'beo': -10.018, 'lif': -11.538, 'f2': -15.808, 'hf': -16.214, 'co': -13.975, 'n2': -15.541, 'bf': -11.175, 'ne': -21.533, 'ph3': -10.628, 'h2o': -12.737, 'hcl': -12.77, 'cs': -11.355, 'ar': -15.739, 'bh3': -13.332, 'nh3': -10.945}
}


########################################################################
def diff(data1, data2):
    """
        Returns a dictionary containing data2 - data1
    """
    errors = dict()
    shared_keys = set(data1.keys()) & set(data2.keys())
    for cas in shared_keys:
        errors[cas] = data2[cas] - data1[cas]
    return errors

########################################################################
def mae_mse_max(data1, data2):
    """
        Returns MAE, MSE, Max errors with data1 being the reference
    """
    errors = list( diff(data1, data2).values() )

    ndata = len(errors)
    mse = np.sum(errors) / float(ndata)
    mae = np.sum(np.abs(errors)) / float(ndata)
    mxe = np.max(np.abs(errors))
    return mae, mse, mxe
