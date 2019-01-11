#!/usr/bin/python3
# This file is part of MOLGW
# Author: Fabien Bruneval

import os, sys, numpy


###################################################################
# First check command line input
#
if len(sys.argv) > 2:
  print('Reading file: '+sys.argv[1])
  print('Reading file: '+sys.argv[2])
else:
  print('Give two file names as arguments')
  sys.exit(0)

filename1 = sys.argv[1]
filename2 = sys.argv[2]

if not os.path.exists(filename1):
  print('File '+filename1+' does not exist')
  sys.exit(0)
if not os.path.exists(filename2):
  print('File '+filename2+' does not exist')
  sys.exit(0)


rho1 = []
w1   = []
f1 = open(filename1,'r')
for line in f1:
  rho1.append( line.split()[0] )
  w1.append( line.split()[1] )
f1.close()

rho2 = []
f2 = open(filename2,'r')
for line in f2:
  rho2.append( line.split()[0] )
f2.close()

if len(rho1) != len(rho2):
  print('Files '+filename1+' and '+filename2+' have different lengths') 
  sys.exit(0)

norm1 = 0.0 
norm2 = 0.0 
chi1  = 0.0
chi2  = 0.0
for i in range(len(w1)):
  norm1 = norm1 + float(w1[i]) * float(rho1[i])
  norm2 = norm2 + float(w1[i]) * float(rho2[i])
  chi1  = chi1  + float(w1[i]) * abs( float(rho1[i]) - float(rho2[i]) )
  chi2  = chi2  + float(w1[i]) *    ( float(rho1[i]) - float(rho2[i]) )**2
chi2 = numpy.sqrt(chi2)

print('Normalization 1: {:.6f}'.format(norm1))
print('Normalization 2: {:.6f}'.format(norm2))
print()
print('           Chi1: {:.6f}'.format(chi1))
print('           Chi2: {:.6f}'.format(chi2))
print()



