#!/usr/bin/python3

###################################
# This file is part of MOLGW
# Author: Fabien Bruneval


###################################
# Running the MOLGW test suite for several compilation options
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
tmpfolder='tmp_manymakes_'+today
tmptestsuite='tmp_'+today

try:
  shutil.rmtree(tmpfolder)
except OSError:
  pass

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
cppflags_list.append('-DHAVE_LIBXC -DHAVE_LIBINT_ONEBODY')
cppflags_list.append('-DHAVE_LIBXC -DHAVE_MPI -DHAVE_SCALAPACK')
cppflags_list.append('-DHAVE_LIBXC -DHAVE_MPI -DHAVE_SCALAPACK -DSELECT_PDSYEVX')
cppflags_list.append('-DHAVE_LIBXC -DHAVE_MPI -DHAVE_SCALAPACK -DHAVE_LIBINT_ONEBODY')

try:
  os.remove(tmpfolder+'/tests.log')
except:
  pass

os.chdir('../src')
shutil.move('my_machine.arch','my_machine.arch.bak')

for imake in range(len(cppflags_list)):

  cppflags = cppflags_list[imake]
  print('\n\nCompilation: {0:} \n'.format(imake+1) )
  print('Compiling with options: '+cppflags)

  fref = open('my_machine.arch.bak','r')
  fnew = open('my_machine.arch','w')
  fout = open('../tests/'+tmpfolder+'/make.log'+str(imake+1),'w')

  fnew.write('CPPFLAGS= '+cppflags+'\n')

  for line in fref:
    if not 'CPPFLAGS' in line:
      fnew.write(line)
  fref.close()
  fnew.close()
  
  print('  make clean')
  subprocess.call(['make','clean'],stdout=fout,stderr=subprocess.STDOUT)
  print('  make -j 4')
  subprocess.call(['make','-j','4'],stdout=fout,stderr=subprocess.STDOUT)
  tested+=1

  fout.close()

  if 'HAVE_MPI' in cppflags :
    nproc = 2
  else :
    nproc = 1

  os.chdir('../tests')
  if os.path.isfile('../molgw'):
    print('  Compilation:  [ \033[92m\033[1mOK\033[0m ] \n')
    success+=1
    print('\nRunning test suite')
    if keeptmp :
      print('  ./run_testsuite.py --np '+str(nproc)+' --keep')
      subprocess.call('./run_testsuite.py --np '+str(nproc)+' --keep >> '+tmpfolder+'/tests.log',shell=True)
    else :
      print('  ./run_testsuite.py --np '+str(nproc))
      subprocess.call('./run_testsuite.py --np '+str(nproc)+' >> '+tmpfolder+'/tests.log',shell=True)


    for line in reversed(open(tmpfolder+'/tests.log','r').readlines()) :
      if 'Succesful tests' in line:
        tests_ok    = int( line.split(':')[1].split('/')[0] )
        tests_total = int( line.split(':')[1].split('/')[1] )
        print( '  Succesful tests: {0:} / {1:}'.format(tests_ok,tests_total) )
        if tests_ok == tests_total :
          print('  Test Suite:   [ \033[92m\033[1mOK\033[0m ] \n')
        else :
          print('  Test Suite:   [\033[91m\033[1mFAIL\033[0m] \n')

        break

    

  else:
    print('  Compilation:  [\033[91m\033[1mFAIL\033[0m] \n')

  if keeptmp :
    shutil.move(tmptestsuite,tmpfolder+'/'+tmptestsuite+'_'+str(imake+1))

  os.chdir('../src')

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
