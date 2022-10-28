#!/usr/bin/python3
##################################################
#
# This file is part of MOLGW
# Author: Fabien Bruneval
#
# This python module provides useful functions to automate MOLGW
#
#
##################################################

import os, sys, subprocess
import difflib
import json
import copy
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

path = __file__.split("/utils")[0]
exe  = path + "/molgw"

#class mgwi:
#
#    def __init__(self,inp):
#        self.input_parameters = inp

########################################################################
def check_input(pyinput):
    sanity = True
    yml = path + '/src/input_variables.yaml'
    with open(yml, 'r') as stream:
        try:
            input_vars = load(stream,Loader=Loader)
        except:
            print('input_variables.yaml file is corrupted')
            pass
    # Check keywords exist
    keywords = [k for k in input_vars.keys() ]
    for k in pyinput:
        if not k.lower() in keywords:
            print('Wrong input variable:    ' + k)
            similar = ''
            for kk in keywords:
                if difflib.SequenceMatcher(None,k,kk).ratio() > 0.6:
                    #print(kk,difflib.SequenceMatcher(None,k,kk).ratio())
                    similar += ' ' + kk
            if len(similar) > 0:
                print(' -> did you mean:   ' + similar)
            else:
                print(' -> no similar keyword found')
            sanity = False

    # Check all mandatory keywords are there
    mandatory = [k for k in input_vars.keys() if input_vars[k]["mandatory"]=="yes" ]
    for k in mandatory:
        if not k in [key.lower() for key in pyinput]:
            print('Mandatory keyword not present:   ' + k)
            sanity = False

    return sanity


########################################################################
def run(inputfile="molgw.in",**kwargs):
    if "pyinput" in kwargs:
        print_input_file(kwargs["pyinput"],inputfile)       
    process = subprocess.Popen([exe,inputfile],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = process.communicate()
    if len(error) > 100:
        print(error.decode("utf-8"))
    with open('molgw.yaml', 'r') as stream:
        try:
            results = load(stream,Loader=Loader)
        except:
            print('molgw.yaml file is corrupted')
            results = {}
            pass
    return results


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
def read_xyz_file(filename):
    with open(filename,'r') as f:
        line = f.readline()
        natom = int(line.strip())
        f.readline()
        structure = []
        for i in range(natom):
            l = f.readline().split()[0:4]
            structure.append([l[0],float(l[1]),float(l[2]),float(l[3])])
        return structure

########################################################################
# structure class is mostly a list of atoms in angstrom
class structure:
    def __init__(self,strucin):
        if type(strucin) == str:
            self.list = read_xyz_file(strucin)
        else:
            self.list = copy.deepcopy(strucin)
    def print_xyz_file(self,filename,comment=""):
        with open(filename,'w') as f:
            f.write('{}\n'.format(len(self.list)))
            f.write(comment.strip()+'\n')
            f.write(self.string())
            #for atom in self.list:
            #    f.write('{:2}   {:14.8f} {:14.8f} {:14.8f}\n'.format(atom[0],float(atom[1]),float(atom[2]),float(atom[3])))
    def __repr__(self):
        return "MOLGW structure (angstrom units)"
    def string(self):
        s = ''
        for atom in self.list:
            s += "{:<2} {:.6f} {:.6f} {:.6f} \n".format(*atom[0:4])
        return s
    def __str__(self):
        return self.string()


########################################################################
def get_homo_energy(approx,calc):
    key = approx + " energy"
    energies = [ ei for ei in calc[key]["spin channel 1"].values()]
    if calc["input parameters"]["nspin"] == 1:
        energies += [ ei for ei in calc[key]["spin channel 1"].values()]
    else:
        energies += [ ei for ei in calc[key]["spin channel 2"].values()]
    energies.sort()
    return energies[int(calc["physical system"]["electrons"])-1 - 2*(min(list(calc[key]["spin channel 1"].keys()))-1) ]


########################################################################
def get_lumo_energy(approx,calc):
    key = approx + " energy"
    energies = [ ei for ei in calc[key]["spin channel 1"].values()]
    if calc["input parameters"]["nspin"] == 1:
        energies += [ ei for ei in calc[key]["spin channel 1"].values()]
    else:
        energies += [ ei for ei in calc[key]["spin channel 2"].values()]
    energies.sort()
    return energies[int(calc["physical system"]["electrons"]) - 2*(min(list(calc[key]["spin channel 1"].keys()))-1) ]


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
def print_input_file(calc,filename="molgw.in"):
    with open(filename,'w') as f:
        f.write('&molgw\n')
        for key, value in calc.items():
            if key == "xyz":
                f.write('  {:30} = {}\n'.format("natom",value.count("\n")) )
            elif type(value) in [type(int()),type(float())]:
                f.write('  {:30} = {}\n'.format(key,value) )
            else:
                f.write('  {:30} = \'{}\'\n'.format(key,value) )
        f.write('/\n')
        if "xyz" in calc:
            f.write(calc["xyz"])


########################################################################
# Conversions for stopping power
def kev_to_au(e_kev,mass=1.0):
    return (1000.*e_kev*2.0/Ha_eV/(1836.1253*mass))**0.5

def au_to_kev(v_au,mass=1.0):
    return 0.5*mass*1836.1253*v_au**2*Ha_eV/1000.


########################################################################
# Load a gaussian cube file into a class
class gaussian_cube:
    atoms_element = []    # list of elements
    atoms_position = []   # list of positions
    nx = 0                # number of grid points along 1st vector
    ny = 0                # number of grid points along 2nd vector
    nz = 0                # number of grid points along 3rd vector
    dx = []               # 1st vector in bohr
    dy = []               # 2nd vector in bohr
    dz = []               # 3rd vector in bohr
    rr = []               # grid points list
    data = []             # volumetric data
    dv = 0.0              # grid point associated volume in bohr^3

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
