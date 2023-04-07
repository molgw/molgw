#!/usr/bin/python3
###################################
# This file is part of MOLGW.
# Author: Fabien Bruneval
# This python script generates the Fortran files: git_sha.f90 and basis_path.f90
###################################


import sys, os, time, shutil, subprocess

today=time.strftime("%d")+' '+time.strftime("%B")+' '+time.strftime("%Y")


###################################
#
# Create the default basis set folder
#
###################################

# Get the location of the ~molgw root
basis_path = __file__.split("/src")[0]

# If an installation path is given in environment variable PREFIX, then override basis_path
basis_path = os.getenv('PREFIX',basis_path)
basis_path = basis_path + '/basis/'

# Write first in a tmp file that will be moved later
# only if it differs from the existing file
# The purpose is to preserve the date of the existing file
# to avoid trigerring re-compilation all the time
fbp = open('basis_path.f90_tmp','w')

fbp.write('!======================================================================\n')
fbp.write('! This file is part of MOLGW.\n')
fbp.write('! The following lines have been generated by python script: prepare_sourcecode.py \n')
fbp.write('! Do not alter them directly: they will be overriden sooner or later by the script\n')
fbp.write('! To add a new input variable, modify the script directly\n')
fbp.write('! Generated by prepare_sourcecode.py\n')
fbp.write('!======================================================================\n\n')
fbp.write(' default_basis_path = &\n  \'' + basis_path + '\' \n\n')

fbp.close()


do_copy = False

if os.path.isfile('basis_path.f90'):
  fout = open('tmp','w')
  subprocess.call(['diff','basis_path.f90_tmp','basis_path.f90'],stdout=fout)
  fout.close()

  fout = open('tmp','r')
  do_copy = ( len(fout.readlines()) != 0 )
  fout.close()
  os.remove('tmp')

else:
  do_copy = True


if do_copy:
  print('Update file: basis_path.f90')
  shutil.copyfile('basis_path.f90_tmp','basis_path.f90')

os.remove('basis_path.f90_tmp')



###################################
#
# Create the version numbering
#
###################################

sha = []
fout = open('tmp_prepare_source','w')
try:
  subprocess.call(['git','rev-parse','HEAD'],stdout=fout,stderr=fout)
  fout.close()
  fout = open('tmp_prepare_source','r')
  for line in fout:
    sha = line.strip()
  if 'fatal' in sha:
    sha = []
except:
  pass

fout.close()

os.remove('tmp_prepare_source')


# Write first in a tmp file that will be moved later
# only if it differs from the existing file
# The purpose is to preserve the date of the existing file
# to avoid trigerring re-compilation all the time
frev = open('git_sha.f90_tmp','w')

frev.write('!======================================================================\n')
frev.write('! This file is part of MOLGW.\n')
frev.write('! The following lines have been generated by python script: prepare_sourcecode.py \n')
frev.write('! Do not alter them directly: they will be overriden sooner or later by the script\n')
frev.write('! To add a new input variable, modify the script directly\n')
frev.write('! Generated by prepare_sourcecode.py\n')
frev.write('!======================================================================\n\n')
if len(sha) > 0:
  frev.write(' git_sha = \'' + sha + '\'\n\n')
else:
  frev.write(' git_sha = \'' + ' -1 ' + '\'\n\n')

frev.close()


do_copy = False

if os.path.isfile('git_sha.f90'):
  fout = open('tmp','w')
  subprocess.call(['diff','git_sha.f90_tmp','git_sha.f90'],stdout=fout)
  fout.close()

  fout = open('tmp','r')
  do_copy = ( len(fout.readlines()) != 0 ) and len(sha) > 0
  fout.close()
  os.remove('tmp')

else:
  do_copy = True


if do_copy:
  print('Update file: git_sha.f90')
  shutil.copyfile('git_sha.f90_tmp','git_sha.f90')

os.remove('git_sha.f90_tmp')
