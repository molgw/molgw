#!/usr/bin/python3
##################################################
#
# This file is part of MOLGW
# Author: Fabien Bruneval
#
# This python module provides useful functions to read molgw.yaml files
# 
#
##################################################

import os, sys
import json
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

bohr_ang =  0.52917721092
Ha_eV    = 27.21138505
periodic_table = [ 'H',                                                                                                  'He', \
                   'Li', 'Be',                                                              'B',  'C',  'N',  'O',  'F', 'Ne', \
                   'Na', 'Mg',                                                             'Al', 'Si',  'P',  'S', 'Cl', 'Ar', \
                    'K', 'Ca', 'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', \
                   'Rb', 'Sr',  'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',  'I', 'Xe', \
                   'Cs', 'Ba', 'La', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', \
                   'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', \
                   'Th', 'Pa',  'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', \
                   'Fr', 'Ra', 'Ac', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' ]


########################################################################
def get_chemical_formula(calc):
    # Try to have a somewhat conventional ordering
    element_list  = [ 'He', 'Ne', 'Ar', 'Kr', 'Xe' ]
    element_list += [ 'Li', 'Na', 'K', 'Rb' ]
    element_list += [ 'Be', 'Mg', 'Ca', 'Sr', 'Ba' ]
    element_list += [ 'B', 'C', 'N', 'P', 'S', 'Si', 'As', 'Se', 'Cu', 'Zn', 'Al', 'Ag', 'Ge', 'Ti', 'Ga']
    element_list += [ 'O', 'F', 'Cl', 'Br', 'I', 'H']

    numbers = [0 for x in element_list]
    
    for atom in calc["physical system"]["atom list"]:
        i = element_list.index(atom[0])
        numbers[i] += 1
    formula = ''
    for e,n in zip(element_list,numbers):
        if n > 1:
            formula += e + str(n)
        elif n > 0:
            formula += e
    return formula

 
########################################################################
def print_xyz_file(calc,filename):
    atom_list = calc["physical system"]["atom list"]
    with open(filename,'w') as f:
        f.write('{}\n\n'.format(len(atom_list)))
        for atom in atom_list:
            f.write('{:2}   {:14.8f} {:14.8f} {:14.8f}\n'.format(atom[0],float(atom[1])*bohr_ang,float(atom[2])*bohr_ang,float(atom[3])*bohr_ang))


########################################################################
def get_homo_energy(approx,calc):
    homo = int(calc["physical system"]["electrons"] * 0.50)
    key = approx + " energy" 
    try:
        ehomo = -9999.9
        for state in calc[key]["spin channel 1"].keys():
            if int(state) <= homo:
                ehomo = max(ehomo,float(calc[key]["spin channel 1"][state]))
        return ehomo
    except:
        return None


########################################################################
def get_lumo_energy(approx,calc):
    homo = int(calc["physical system"]["electrons"] * 0.50)
    key = approx + " energy" 
    try:
        elumo = 9999.9
        for state in calc[key]["spin channel 1"].keys():
            if int(state) > homo:
                elumo = min(elumo,float(calc[key]["spin channel 1"][state]))
        return elumo
    except:
        return None


########################################################################
# returns a list of dictionaries: one dictionary per yaml file in the directory
def parse_yaml_files(directory):
    # List all the yaml files in the directory
    yaml_files = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if '.yaml' in filename:
                yaml_files.append(os.path.join(root, filename))
    print('{} molgw.yaml files identified in directory '.format(len(yaml_files)) + directory)

    # Read all the yaml files -> dictionary
    calc = []
    for yaml_file in yaml_files:
        with open(yaml_file, 'r') as stream:
            try:
                calc.append(load(stream,Loader=Loader))
            except:
                print(yaml_file + ' is corrupted')
                pass
    return calc


########################################################################
# check if a calculation dictionary is valid
# - "scf is converged" exists
# - "scf is converged" is True
# - "run" exists ("run" key is written in YAML file at the very end of a MOLGW calculation)
def check_calc(calc):
    valid = True
    try:
        calc["scf is converged"]
        calc["run"]
    except KeyError:
        valid = False
    except:
        sys.exit(1)
    return valid and calc["scf is converged"]


########################################################################
def create_gw100_json(filename,data,**kwargs):
    dict_gw100 = dict()
    dict_gw100["code"]= "MOLGW"
    dict_gw100["code_version"]= "2.E"
    dict_gw100["basis"]= "gaussian"
    dict_gw100["qpe"]= "solved"
    dict_gw100["DOI"]= "unpublished"

    dict_gw100.update(kwargs)
    dict_gw100["data"] = data

    with open(filename, 'w') as json_file:
        json.dump(dict_gw100,json_file,indent=2,separators=(',', ': '))


########################################################################
def kev_to_au(mass,e_kev):
    return (1000.*e_kev*2.0/Ha_eV/(1836.1253*mass))**0.5

def au_to_kev(mass,v_au):
    return 0.5*mass*1836.1253*v_au**2*Ha_eV/1000.
    

########################################################################
class gaussian_cube:
    atoms_element = []    # list of elements
    atoms_position = []   # list of positions
    nx = 0                # number of grid points along 1st vector
    ny = 0                # number of grid points along 2nd vector
    nz = 0                # number of grid points along 3rd vector
    dx = []               # 1st vector
    dy = []               # 2nd vector
    dz = []               # 3rd vector
    rr = []               # grid points list
    data = []             # volumetric data
    dv = 0.0              # grid point associated volume

    # Initialize class with a file "filename"
    def __init__(self,filename):

        cf=open(filename,'r')
        cf.readline()
        cf.readline()

        line = cf.readline().split()
        natom = int(line[0])
        r0 = [float(x) for x in line[1:5]]
        line = cf.readline().split()
        self.nx = int(line[0])
        self.dx = [float(x) for x in line[1:5]]
        line = cf.readline().split()
        self.ny = int(line[0])
        self.dy = [float(x) for x in line[1:5]]
        line = cf.readline().split()
        self.nz = int(line[0])
        self.dz = [float(x) for x in line[1:5]]
        # atom list
        self.atoms_element = []
        self.atoms_position = []
        for i in range(natom):
            line = cf.readline().split()
            self.atoms_element.append(int(line[0]))
            self.atoms_position.append([float(line[2]), float(line[3]), float(line[4])])

        # volumetric data
        self.data = []
        for line in cf:
             self.data.extend( float(x) for x in line.split() )
        cf.close()

        self.rr = []
        for ix in range(self.nx):
            for iy in range(self.ny):
                for iz in range(self.nz):
                    x = r0[0] + ix * self.dx[0] + iy * self.dy[0] + iz * self.dz[0]
                    y = r0[1] + ix * self.dx[1] + iy * self.dy[1] + iz * self.dz[1]
                    z = r0[2] + ix * self.dx[2] + iy * self.dy[2] + iz * self.dz[2]
                    self.rr.append([x,y,z])
        self.dv = self.dx[0] * self.dy[1] * self.dz[2] \
                 +self.dx[1] * self.dy[2] * self.dz[0] \
                 +self.dx[2] * self.dy[0] * self.dz[1] \
                 -self.dx[2] * self.dy[1] * self.dz[0] \
                 -self.dx[0] * self.dy[2] * self.dz[1] \
                 -self.dx[1] * self.dy[0] * self.dz[2]
        return

########################################################################
