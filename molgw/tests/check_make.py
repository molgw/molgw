#!/usr/bin/python3

###################################
# This file is part of MOLGW


###################################
# Running the MOLGW test suite
###################################


import sys, os, time, shutil, subprocess

today=time.strftime("%Y")+'_'+time.strftime("%m")+'_'+time.strftime("%d")
start_time = time.time()
keeptmp = False


     
###################################
# Parse the command line


if len(sys.argv) > 1:
  if '--help' in sys.argv:
    print('Run several compilations of MOLGW to check the compilation options')
    print('  --keep             Keep the temporary folder')
    sys.exit(0)
  if '--keep' in sys.argv:
    keeptmp = True

print('\n===============================')
print('Starting MOLGW make check\n')


###################################
tmpfolder='tmp_'+today

try:
  os.mkdir(tmpfolder)
except OSError:
  pass



###################################
# Loop over the CPPFLAGS
###################################

success = 0
tested = 0
cppflags_list = []

cppflags_list.append('-DHAVE_LIBXC')
cppflags_list.append('-DHAVE_LIBXC -DHAVE_MPI')
cppflags_list.append('-DHAVE_LIBXC -DHAVE_MPI -DSCALAPACK')
cppflags_list.append('-DHAVE_LIBXC -DHAVE_MPI -DSCALAPACK -DSELECT_PDSYEVX')

os.chdir('../src')
shutil.move('my_machine.arch','my_machine.arch.bak')

for imake in range(len(cppflags_list)):

  cppflags = cppflags_list[imake]
  print('Compilation: ',imake+1)
  print('Compiling with options: '+cppflags)

  fref = open('my_machine.arch.bak','r')
  fnew = open('my_machine.arch','w')
  fout = open('../tests/'+tmpfolder+'/make.log'+str(imake),'w')

  fnew.write('CPPFLAGS= '+cppflags+'\n')

  for line in fref:
    if not 'CPPFLAGS' in line:
      fnew.write(line)
  fref.close()
  fnew.close()
  
  print('  make clean')
  subprocess.call(['make','clean'],stdout=fout,stderr=subprocess.STDOUT)
  print('  make -j 2')
  subprocess.call(['make','-j','2'],stdout=fout,stderr=subprocess.STDOUT)
  tested+=1

  fout.close()

  if os.path.isfile('../molgw'):
    print('Compiltation:  [ \033[92m\033[1mOK\033[0m ] \n')
    success+=1
  else:
    print('Compiltation:  [\033[91m\033[1mFAIL\033[0m] \n')


fout = open('../tests/'+tmpfolder+'/make.log_final','w')
print('make clean')
subprocess.call(['make','clean'],stdout=fout,stderr=subprocess.STDOUT)
fout.close()

shutil.move('my_machine.arch.bak','my_machine.arch')
os.chdir('../tests')




###################################

print('\n\n===============================')
print('      Test Summary \n')
print('  Succesful tests:   ',success,' / ',tested,'\n')
print(' Elapsed time (s):   ','{:.2f}'.format(time.time() - start_time) )
print('===============================\n')

###################################
# Erase the tmp folder by default
if not keeptmp:
  shutil.rmtree(tmpfolder)

###################################
