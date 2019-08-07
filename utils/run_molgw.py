#!/usr/bin/python3
##################################################
#
# This file is part of MOLGW
# Author: Fabien Bruneval
#
# This python script prepares a series of MOLGW input files
# for all the XYZ files found in directory named 'structures'
# 
#
##################################################
  

import os, sys, shutil, stat, subprocess

##################################################
#
# Definitions
#
def longname(basis):
    if basis == 'DZ':
        return 'cc-pVDZ'
    elif basis == 'TZ':
        return 'cc-pVTZ'
    elif basis == 'QZ':
        return 'cc-pVQZ'
    elif basis == '5Z':
        return 'cc-pV5Z'
    elif basis == 'ADZ':
        return 'aug-cc-pVDZ'
    elif basis == 'ATZ':
        return 'aug-cc-pVTZ'
    elif basis == 'AQZ':
        return 'aug-cc-pVQZ'
    elif basis == 'A5Z':
        return 'aug-cc-pV5Z'

class input_parameters:
    def __init__(self, scf='bhlyp', postscf='gw', selfenergy_state_range=20, frozencore='yes', basis='ADZ'):
        self.scf = scf
        self.postscf = postscf
        self.selfenergy_state_range = selfenergy_state_range
        self.frozencore = frozencore
        self.basis = basis

    def print_input_file(self,f):
        for attr, value in self.__dict__.items():
            if not attr == 'basis':
                if type(value) in [type(int()),type(float())] :
                    f.write('  {:30} = {}\n'.format(attr,value) )
                else:
                    f.write('  {:30} = \'{}\'\n'.format(attr,value) )
            else:
                f.write('  {:30} = \'{}\'\n'.format(attr,longname(value)) )
                f.write('  {:30} = \'{}-RI\'\n'.format('auxil_basis',longname(value)) )



##################################################
#
# Hard-coded information
#
directory       = 'run_001'
executable      = '/home/bruneval/devel/molgw-devel/molgw'
run_it          = False   # run it from python script or wait 
atom_number_max = 14      # limit the size of the calculated molecules


# Create the calculation list here
ip = []
for basis in ['ADZ','ATZ','AQZ','A5Z']:
    ip.append(input_parameters(basis=basis))


#########################################
# Implement a size limit
#
#molecule_list = [ filexyz.replace('.xyz','') for filexyz in os.listdir("structures") ]

molecule_list = []
for filexyz in os.listdir("structures"):
    with open('structures/'+filexyz) as f:
        if int(f.readline()) <= atom_number_max:
            molecule_list.append(filexyz.replace('.xyz',''))



# Molecule list
print('=========== Molecule list')
print(molecule_list)
print('==========================')


#########################################
#
#
os.makedirs(directory,exist_ok=True)
 
script = open('run.sh','w')

for molecule in molecule_list:
    for calc in ip:

        folder_name = molecule + '_' + calc.basis + '_' + calc.scf
        folder = directory + '/' + folder_name
        os.makedirs(folder,exist_ok=True)
  
        #
        # Check if the molgw.yaml is already there and finalized
        #
        try:
            with open(folder + '/molgw.yaml', 'r') as f:
                last_line = f.readlines()[-1]
                if '...' in last_line:
                    print('{:24} {:5} is already calculated. Skip it'.format(molecule,calc.basis))
                    continue
        except:
            pass
            
        print('{:24} {:5} is being calculated'.format(molecule,calc.basis))
  
        os.chdir(folder)
  
        script.write('cd ' + folder + '\n')
        
        with open('molgw.in','w') as fin:
            fin.write('&molgw\n')
            fin.write('  comment                 = \'' + molecule + '\'\n\n')
            calc.print_input_file(fin)
            fin.write('  print_yaml              = \'yes\'\n')
            fin.write('  print_spatial_extension = \'yes\'\n')
            fin.write('  print_sigma             = \'yes\'\n')
            fin.write('  xyz_file                = \'../../structures/' + molecule + '.xyz\'\n')
            fin.write('/\n')
  
        script.write(executable + ' molgw.in > molgw.out\n')
        if run_it:
            process = subprocess.Popen([executable,'molgw.in'],stdout=subprocess.PIPE)
            output, error = process.communicate()
  
        script.write('cd ../.. \n')
  
        os.chdir('../../')

script.close()
os.chmod('run.sh',stat.S_IRUSR+stat.S_IWUSR+stat.S_IXUSR)


