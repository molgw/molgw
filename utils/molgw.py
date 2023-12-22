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

"""
molgw module contains classes and modules to automate running and reading of MOLGW.
"""

__author__  = "Fabien Bruneval"
__version__ = "3.2"

import math
import os, sys, shutil, subprocess
import difflib
import json
import copy
import pathlib, glob
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

bohr_ang =  0.52917721092
Ha_eV    = 27.21138505
c_speedlight = 137.035999074   # speed of light in atomic units (fine-structure constant inverse)
periodic_table = [ 'H',                                                                                                  'He',
                   'Li', 'Be',                                                              'B',  'C',  'N',  'O',  'F', 'Ne',
                   'Na', 'Mg',                                                             'Al', 'Si',  'P',  'S', 'Cl', 'Ar',
                    'K', 'Ca', 'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                   'Rb', 'Sr',  'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',  'I', 'Xe',
                   'Cs', 'Ba', 
                   'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 
                   'Lu', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
                   'Fr', 'Ra',
                   'Ac', 'Th', 'Pa',  'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 
                   'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt' 
                ]
z_element = {element: index+1 for index, element in enumerate(periodic_table)}

molgw_rootfolder = str(pathlib.Path(__file__).resolve().parent.parent)
exe  = molgw_rootfolder + "/molgw"


########################################################################
def check_input(pyinput):
    sanity = True
    yml = molgw_rootfolder + '/src/input_variables.yaml'
    with open(yml, 'r') as stream:
        try:
            input_vars = load(stream,Loader=Loader)
        except:
            sys.exit('input_variables.yaml file is corrupted')
            pass

    valid_keywords = [k for k in input_vars.keys() ]
    additional_keywords = ["xyz", "rawxyz"]
    valid_keywords += additional_keywords

    # Check keywords exist
    pyinput_lower = [ k.lower() for k in pyinput ]
    for k in pyinput_lower:
        if not k in valid_keywords:
            print('Wrong input variable:    ' + k)
            similar = ''
            for kk in valid_keywords:
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
        if not k in [key for key in pyinput_lower]:
            print('Mandatory keyword not present:   ' + k)
            sanity = False
    # Check that some sort of structure is there
    structure_kw = [ "natom", "xyz_file", "xyz", "rawxyz"]
    if not any(kw in pyinput_lower for kw in structure_kw):
        print("No structural data given")
        sanity = False

    return sanity


########################################################################
def run(inputfile="molgw.in",outputfile="molgw.out",pyinput={},mpirun="",executable_path="",openmp=1,tmp="",**kwargs):
    if len(tmp) > 0:
        os.makedirs(tmp,exist_ok=True)
        current_directory = os.getcwd()
        new_working_directory = current_directory + '/' + tmp
        #os.chdir(new_working_directory)
        os.chdir(tmp)
    if len(executable_path) > 0:
        exe_local = executable_path
    else:
        exe_local = exe
    if len(pyinput) > 0:
        print_input_file(pyinput,inputfile)
    os.environ['OMP_NUM_THREADS'] = str(openmp)
    if len(mpirun) == 0:
        process = subprocess.Popen([exe_local,inputfile],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    else:
        process = subprocess.Popen(mpirun.split()+[exe_local,inputfile],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    output, error = process.communicate()
    if len(outputfile) >0:
        with open(outputfile,'w') as f:
            f.write(output.decode("utf-8"))
    if len(error) > 100:
        print(error.decode("utf-8"))
    with open('molgw.yaml', 'r') as stream:
        try:
            results = load(stream,Loader=Loader)
        except:
            print('molgw.yaml file is corrupted')
            results = {}
            pass
    if len(tmp) > 0:
        os.chdir(current_directory)
        shutil.rmtree(tmp)
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
# with a few method to read, print, transform
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
    key = approx + " energies"
    if key not in calc.keys():
        print(f"Problem reading calculation: {calc['input parameters']['comment']}")
        print(f"{key} not found")
        sys.exit("Problem")
    energies = [ float(ei) for ei in calc[key]["spin channel 1"].values()]
    if calc["input parameters"]["nspin"] == 1:
        energies += [ float(ei) for ei in calc[key]["spin channel 1"].values()]
    else:
        energies += [ float(ei) for ei in calc[key]["spin channel 2"].values()]
    energies.sort()
    return energies[int(calc["physical system"]["electrons"])-1 - 2*(min(list(calc[key]["spin channel 1"].keys()))-1) ]


########################################################################
def get_lumo_energy(approx,calc):
    key = approx + " energies"
    if key not in calc.keys():
        print(f"Problem reading calculation: {calc['input parameters']['comment']}")
        print(f"{key} not found")
        sys.exit("Problem")
    energies = [ float(ei) for ei in calc[key]["spin channel 1"].values()]
    if calc["input parameters"]["nspin"] == 1:
        energies += [ float(ei) for ei in calc[key]["spin channel 1"].values()]
    else:
        energies += [ float(ei) for ei in calc[key]["spin channel 2"].values()]
    energies.sort()
    index = int(calc["physical system"]["electrons"]) - 2*(min(list(calc[key]["spin channel 1"].keys()))-1)
    if index >= len(energies):
        sys.exit(f"not enough states to get LUMO in {approx} energies")
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
    dict_gw100["code_version"]= __version__
    dict_gw100["basis"]= "gaussian"
    dict_gw100["qpe"]= "solved"
    dict_gw100["DOI"]= "unpublished"

    dict_gw100.update(kwargs)
    dict_gw100["data"] = data

    with open(filename, 'w') as json_file:
        json.dump(dict_gw100,json_file,indent=2,separators=(',', ': '))


########################################################################
def print_input_file(pyinput,filename="molgw.in"):
    with open(filename,'w') as f:
        f.write('&molgw\n')
        for key, value in pyinput.items():
            if key == "xyz":
                f.write('  {:30} = {}\n'.format("natom",value.count("\n")) )
            elif key == "rawxyz":
                continue
            elif key == "vel_projectile":
                f.write('  {:30} = {}\n'.format(key,value) )
            elif type(value) in [type(int()),type(float())]:
                f.write('  {:30} = {}\n'.format(key,value) )
            else:
                f.write('  {:30} = \'{}\'\n'.format(key,value) )
        f.write('/\n')
        if "xyz" in pyinput:
            f.write(pyinput["xyz"])
        if "rawxyz" in pyinput:
            f.write(pyinput["rawxyz"])


########################################################################
# Conversions for stopping power

# velocity conversion
def vel_kev_to_au(e_kev,mass=1.0):
    return (1000.*e_kev*2.0/Ha_eV/(1836.1253*mass))**0.5

def velrel_kev_to_au(e_kev,mass=1.0):
    mc2 = 1836.1253*mass*c_speedlight**2
    e_au = 1000.*e_kev/Ha_eV
    return c_speedlight*(1.0-(mc2/(e_au+mc2))**2)**0.5

def vel_au_to_kev(v_au,mass=1.0):
    return 0.5*mass*1836.1253*v_au**2*Ha_eV/1000.

# stopping cross section (S/rho) conversion
# Srim is in 1e-15 eV cm**2 /atom
def scs_srim_to_au(scs_srim):
    return scs_srim / Ha_eV * (bohr_ang * 1.0e-8) / ( 1.0e15 * (bohr_ang * 1.0e-8)**3 )

def scs_au_to_srim(scs_au):
    return scs_au * Ha_eV / (bohr_ang * 1.0e-8) * ( 1.0e15 * (bohr_ang * 1.0e-8)**3 )

def se_au_to_kevpernm(se_au):
    return se_au * (Ha_eV * 1.0e-3) / (bohr_ang * 0.10)

def se_kevpernm_to_au(se_kevpernm):
    return se_kevpernm / (Ha_eV * 1.0e-3) * (bohr_ang * 0.10)

# Classic Bethe-Formula in atomic units
def scs_bethe_au(v_au,nelec,ionization,charge=1.0):
     return 4.0*math.pi*nelec*charge**2/v_au**2 * math.log( 2.0*v_au**2/ionization)

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
class Molgw_input:
    """MOLGW input"""
    def __init__(self,dict_in):
        self.d = dict_in
    def __str__(self):
        return str(self.d)
    def __getitem__(self, key):
        return self.d[key]
    def set(self,key,value):
        self.d[key] = value
    def check(self):
        return check_input(self.d)
    def to_dict(self):
        return self.d
    def to_file(self,filename):
        return print_input_file(self.d,filename)
    def run(self,**kwargs):
        return Molgw_output(run(pyinput=self.d,**kwargs))


########################################################################
class Molgw_output:
    """MOLGW output"""
    def __init__(self,dict_in):
        self.d = dict_in
    def __str__(self):
        return str(self.d)
    def __getitem__(self, key):
        return self.d[key]
    def keys(self):
        return [ k for k in self.d.keys() ]
    def homo_energy(self,approx):
        return get_homo_energy(approx,self.d)
    def lumo_energy(self,approx):
        return get_lumo_energy(approx,self.d)
    def check(self):
        return check_calc(self.d)
    def to_dict(self):
        return self.d
    def chemical_formula(self):
        return get_chemical_formula(self.d)


########################################################################
class Molgw_outputs:
    """MOLGW collection of outputs"""
    def __init__(self, origin=''):
        self.files = []
        self.data = []
        if len(origin) == 0:
            return
        if isinstance(origin, list):
            origins = origin
        else:
            origins = [origin]
        for orig in origins:
            if not os.path.isdir(orig):
                sys.exit(orig + "is not a valid folder")
            self.files += glob.glob(orig+"/**/*.yaml",recursive=True)
        self.files = list(set(self.files))
        for file in self.files:
            with open(file,'r') as f:
                self.data.append(Molgw_output(load(f,Loader=Loader)))
    def __len__(self):
        return len(self.files)
    def __str__(self):
        s = f'MOLGW results from {len(self.files)} files'
        for i,file in enumerate(self.files):
            s+= f"\n - {file:<30}: {self.data[i].get('input parameters')['comment']}"
        return s
    def __iter__(self):
        self.current = 0
        return self
    def __next__(self):
        if self.current < len(self):
            result = self.data[self.current]
            self.current += 1
            return result
        else:
            raise StopIteration
        return self.data
    def __getitem__(self, index):
        return self.data[index]

    # Returns a copy of self containing all those calculations
    # that match the input parameters mentioned in "filters" dictionary
    def filtering(self, filters, verbose=False):
       mlgo_filtered = Molgw_outputs()
       if verbose:
           print("Selection rules:")
           for key, value in filters.items():
               print(f"  {key} == {value}?")
       for f, mlgo in zip(self.files,self.data):
           corresponds = True
           for key, value in filters.items():
               if mlgo.get("input parameters")[key] != value:
                   corresponds = False
           if corresponds:
               mlgo_filtered.files.append(f)
               mlgo_filtered.data.append(mlgo)
       if verbose:
           print(f"Found {len(mlgo_filtered)} corresponding calculations")
       return mlgo_filtered


########################################################################
