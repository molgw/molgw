#!/usr/bin/python3

import os, sys, numpy

list_am       = 'S P D F G H I       SP'
list_elements = [ 'H',                                                                                                 'He', \
                 'Li', 'Be',                                                              'B',  'C',  'N',  'O',  'F', 'Ne', \
                 'Na', 'Mg',                                                             'Al', 'Si',  'P',  'S', 'Cl', 'Ar', \
                  'K', 'Ca', 'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', \
                 'Rb', 'Sr',  'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',  'I', 'Xe', \
                 'Cs', 'Ba', 'La', 'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn']
                 


def fortran_float( string ):
  return float(string.replace('D','E').replace('d','e'))


###################################################################
# First check command line input
#
if len(sys.argv) > 1:
  print('Reading file: '+sys.argv[1])
else:
  print('Give a file name as an argument')
  sys.exit(0)

filename=sys.argv[1]

if not os.path.exists(filename):
  print('File '+filename+' does not exist')
  sys.exit(0)




###################################################################
# Second read input file
#

nwfile = open(filename,'r')

alpha  = [] 
angmom = []
coeff  = []
basis_set = []
new_basisfunction = False

for line in nwfile:

  #
  # Should we skip the line?
  # 

  # Skip empty lines
  if len(line) == 1:
    continue

  # If the line is not empty, split it
  parsing = line.split()

  # If the line starts with #, it is ignored
  if parsing[0][0] == '#':
    continue
  # If the line starts with BASIS, it is ignored
  if parsing[0] == 'BASIS':
    continue
  #
  # If the line starts with END, the file is finished
  if parsing[0] == 'END':
    break


  #
  # Check if the first string is an element
  if any(parsing[0] == element for element in list_elements):

    element = parsing[0]

    # Find if the found element is already in the list
    already = False
    for index_basis,basis in enumerate(basis_set):
      if element == basis:  
        already = True
        break
    # If not found, then create a new basis in basis_set
    if not already:
      basis_set.append(element)
      index_basis = len(basis_set)-1
      alpha.append([])
      coeff.append([])
      angmom.append([])

    angmom_tmp = int(list_am.index(parsing[1])/2)

    new_basisfunction = True

 
  #
  # It is not an element, then it has to be a coefficient list: alpha, c1, c2, c3...
  else:

    if new_basisfunction == True:
      new_basisfunction = False

      functions         = len(alpha[index_basis])
      new_functions     = len(parsing)-1

      for i in range(new_functions):
        alpha[index_basis].append([])
        coeff[index_basis].append([])
        angmom[index_basis].append(angmom_tmp)

      alpha_tmp = fortran_float(parsing[0])
      for i in range(functions,functions+new_functions):
        j = i - functions + 1
        coeff_tmp = fortran_float(parsing[j])
        if abs(coeff_tmp) > 1.e-8 :
          alpha[index_basis][i].append(alpha_tmp)
          coeff[index_basis][i].append(coeff_tmp)

    else:

      alpha_tmp = fortran_float(parsing[0])
      for i in range(functions,functions+new_functions):
        j = i - functions + 1
        coeff_tmp = fortran_float(parsing[j])
        if abs(coeff_tmp) > 1.e-8 :
          alpha[index_basis][i].append(alpha_tmp)
          coeff[index_basis][i].append(coeff_tmp)


nwfile.close()

print('File read')
print()
print(str(len(basis_set)) + ' elements found')
print()




###################################################################
# Thrid write MOLGW basis files
#

newfiles = 0 
for index_basis,basis in enumerate(basis_set):

  outfilename = '../basis/' + basis + '_' + filename.partition('.')[0].partition('/')[2]
  # If file already exists, then skip writing
  if os.path.exists(outfilename):
    print("{:36s}".format(outfilename) + ' skipped since it already exits')
    continue

  print("{:36s}".format(outfilename) + ' generated')
  newfiles = newfiles + 1
  outfile=open(outfilename,'w')
  outfile.write("{:d}".format(len(alpha[index_basis])) + '\n' )

  for func in range(len(alpha[index_basis])):
    nprim = len(alpha[index_basis][func])
    outfile.write("{:d}  {:d}".format(nprim, angmom[index_basis][func]) + '\n' )

    for iprim in range(nprim):
      outfile.write("  {:16.8e}  {:16.8e}".format(alpha[index_basis][func][iprim],coeff[index_basis][func][iprim]) + '\n' )


  outfile.close()


print()
print(str(newfiles) + ' files have been generated')
print()


###################################################################
