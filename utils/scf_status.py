#!/usr/bin/python3
# This file is part of MOLGW
# Author: Fabien Bruneval

import os, sys, numpy


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



outfile = open(filename,'r')

scf_cycle = []
etot = []
homo = []
homo_lumo = []
residual = []
for line in outfile:
  if 'SCF cycle No:' in line:
    parsing = line.split(':')
    scf_cycle.append(int(parsing[1]))
    
  if 'Total Energy    (Ha):' in line:
    parsing = line.split(':')
    etot.append(float(parsing[1]))

  if 'gKS HOMO energy    (eV):' in line:
    parsing = line.split(':')
    homo.append(float(parsing[1]))

  if 'gKS HOMO-LUMO gap  (eV):' in line:
    parsing = line.split(':')
    homo_lumo.append(float(parsing[1]))

  if 'Convergence criterium on the density matrix' in line:
    parsing = line.split(' ')
    ctmp = parsing[len(parsing)-1]
    ctmp = ctmp.replace('\n','')
    residual.append(float(ctmp))

ncycle = min( len(scf_cycle),len(etot),len(homo),len(residual),len(homo_lumo) )

ediff = []
ediff.append(0.0)
for iscf in range(1,ncycle):
  ediff.append(etot[iscf] - etot[iscf-1])

eerror = []
for iscf in range(0,ncycle):
  eerror.append(etot[iscf] - etot[ncycle-1])


print('{0:d} completed cycles found \n'.format(ncycle))

print('SCF cycle          Etotal               Diff                 Deviation         HOMO          H-L Gap        Residual')
for iscf in range(ncycle):
  print('  {0:4d}     {1:18.10f}     {2:+16.10f}      {3:+16.10f}     {4:8.4f}       {5:8.4f}      {6:12.6f}  '  \
              .format(scf_cycle[iscf],etot[iscf],ediff[iscf],eerror[iscf],homo[iscf],homo_lumo[iscf],residual[iscf]))


sys.exit(0)



