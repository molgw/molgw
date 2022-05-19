#!/usr/bin/python3
# This file is part of MOLGW
# Author: Fabien Bruneval

import os, sys
import numpy as np
import molgw



def readdata(f,n,data_type='float'):
    data = []
    if data_type == 'float':
        for i in range((n-1)//5+1):
            line = f.readline()
            data.extend([float(x) for x in line.split()])
    else:
        for i in range((n-1)//6+1):
            line = f.readline()
            data.extend([int(x) for x in line.split()])
    return data


if len(sys.argv) > 1:
    print('Reading file: '+sys.argv[1])
else:
    print('Please provide a Gaussian fchk file')
    sys.exit(1)


filename=sys.argv[1]

f=open(filename,'r')


for line in f:

    if 'Number of atoms' in line:
        natom = int(line.split()[-1])
        if natom != 1:
            print('just coded for single atoms')
            sys.exit(1)

    if 'Atomic numbers' in line:
        element = molgw.periodic_table[int(f.readline())-1]
        print('Element:              ' + element)

    if 'Number of contracted shells' in line:
        nshell = int(line.split()[-1])
        print('Number of shells:     {}'.format(nshell))
    if 'Number of primitive shells' in line:
        nprim  = int(line.split()[-1])

    if 'Shell types' in line:
        angmom = readdata(f,nshell,data_type='int')
        print('Max angular momentum: {:2d}'.format(max(np.abs(angmom))))
    if 'Number of primitives per shell' in line:
        prim = readdata(f,nshell,data_type='int')
    if 'Primitive exponents' in line:
        expo = readdata(f,nprim)
    if 'Contraction coefficients' in line:
        coef = readdata(f,nprim)

f.close()

print('Output basis file: ' + element + '_gaussian')
fb=open(element + '_gaussian','w')

fb.write(str(nshell) + '\n')
iprim=0
for l,p in zip(angmom,prim):
    fb.write('{} {}\n'.format(p,abs(l)))
    for j in range(iprim,iprim+p):
        fb.write('{:16.8e} {:16.8e}\n'.format(expo[j],coef[j]))
    iprim+=p


fb.close()



