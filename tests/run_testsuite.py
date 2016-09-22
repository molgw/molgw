#!/usr/bin/python3

###################################
# This file is part of MOLGW
# Author: Fabien Bruneval


###################################
# Running the MOLGW test suite
###################################


import sys, os, time, shutil, subprocess

today=time.strftime("%Y")+'_'+time.strftime("%m")+'_'+time.strftime("%d")
start_time = time.time()
keeptmp = False

selected_input_file=''
mpirun=''
nprocs=1
debug=False



###################################
def clean_run(inp,out,restart):
  shutil.copy('inputs/'+inp,tmpfolder+'/'+inp)
  os.chdir(tmpfolder)
  if not restart:
    try:
      os.remove('RESTART')
    except FileNotFoundError:
      pass
    try:
      os.remove('ENERGY_QP')
    except FileNotFoundError:
      pass
    try:
      os.remove('SCREENED_COULOMB')
    except FileNotFoundError:
      pass
  fout = open(out, 'w')
  if len(mpirun) < 1:
    subprocess.call(['../../molgw',inp],stdout=fout,stderr=subprocess.STDOUT)
  else:
    subprocess.call([mpirun,'-n',str(nprocs),'../../molgw',inp],stdout=fout,stderr=subprocess.STDOUT)
  fout.close()
  os.chdir('..')


###################################
def check_output(out,testinfo):
  global success,tested

  #
  # First check if the test was aborted because of parallelization not available
  #
  for line in open(tmpfolder+'/'+out,'r').readlines():
    if 'one CPU only' in line:
      print('Test not functional in parallel => skip it')
      return
  #
  # Second check if there is a memory leak
  #
  key = '    Memory ('
  ref = 0.000
  tol = 0.001
  key_found = False
  tested += 1
  for line in reversed(open(tmpfolder+'/'+out,'r').readlines()):
    if key in line:
      key_found = True
      parsing  = line.split(':')
      parsing2 = parsing[1].split()
      if abs( float(parsing2[0]) - ref ) < tol:
        print('No memory leak'.rjust(30)+'[ \033[92m\033[1mOK\033[0m ]'.rjust(30))
        success += 1
        fdiff.write(str(tested).rjust(6) + parsing2[0].rjust(30) \
              + str(ref).rjust(30)+str(float(parsing2[0]) - ref).rjust(30)+'  OK  \n')
        break
      else:
        print('No memory leak'.rjust(30)+'[\033[91m\033[1mFAIL\033[0m]'.rjust(30))
        fdiff.write(str(tested).rjust(6) + parsing2[0].rjust(30) \
              + str(ref).rjust(30)+str(float(parsing2[0]) - ref).rjust(30)+' FAIL \n')
        break

  if not key_found:
    print('No memory leak'.rjust(30)+'[\033[91m\033[1mNOT FOUND\033[0m]'.rjust(30))

  #
  # Then, parse the output and perform the checks
  # 
  for itest in range(len(testinfo)):
    tested += 1
    key = testinfo[itest][0].strip()
    ref = float(testinfo[itest][1])
    pos = int(  testinfo[itest][2])
    tol = float(testinfo[itest][3])

    if debug:
      print('===debug:')
      print('open file: '+tmpfolder+'/'+out)
      print('===end debug')

    key_found = False

    for line in reversed(open(tmpfolder+'/'+out,'r').readlines()):
      if key in line:
        key_found = True
        parsing  = line.split(':')
        if debug:
          print('===debug:')
          print(parsing)
          print('===end debug')
        parsing2 = parsing[1].split()
        if debug:
          print('===debug:')
          print(parsing2)
          print('===end debug')


        if abs( float(parsing2[pos]) - ref ) < tol:
          print(key.rjust(30)+'[ \033[92m\033[1mOK\033[0m ]'.rjust(30))
          success += 1
          fdiff.write(str(tested).rjust(6) + parsing2[pos].rjust(30) \
                + str(ref).rjust(30)+str(float(parsing2[pos]) - ref).rjust(30)+'  OK  \n')
          break
        else:
          print(key.rjust(30)+'[\033[91m\033[1mFAIL\033[0m]'.rjust(30))
          fdiff.write(str(tested).rjust(6) + parsing2[pos].rjust(30) \
                + str(ref).rjust(30)+str(float(parsing2[pos]) - ref).rjust(30)+' FAIL \n')
          break
    if not key_found:
      print(key.rjust(30)+'[\033[91m\033[1mNOT FOUND\033[0m]'.rjust(30))
     
###################################
# Parse the command line


if len(sys.argv) > 1:
  if '--help' in sys.argv:
    print('Run the complete test suite of MOLGW')
    print('  --keep             Keep the temporary folder')
    print('  --np     n         Set the number of cores to n')
    print('  --mpirun launcher  Set the MPI launcher name')
    print('  --input  file      Only run this input file')
    print('  --debug            Outputs debug information for this script')
    sys.exit(0)
  if '--keep' in sys.argv:
    keeptmp = True
  if '--np' in sys.argv:
    i = sys.argv.index('--np') + 1
    nprocs = int( sys.argv[i] )
    mpirun='mpirun'
  if '--mpirun' in sys.argv:
    i = sys.argv.index('--mpirun') + 1
    mpirun = sys.argv[i]
  if '--input' in sys.argv:
    i = sys.argv.index('--input') + 1
    selected_input_file = sys.argv[i]
  if '--debug' in sys.argv:
    debug = True

print('\n===============================')
print('Starting MOLGW test suite\n')

if len(mpirun) < 1:
  if( nprocs > 1 ):
    print('No MPI launcher has been provided. Set the number of cores back to 1')
  nprocs = 1


print('Running with ',nprocs,'cores')

###################################
# Parse the file testsuite
###################################
ninput = 0
input_files = []
restarting  = []
test_names  = []
testinfo    = []

ftestsuite = open('inputs/testsuite','r')
for line in ftestsuite:
  # Removing the comments
  parsing0 = line.split('#')
  # Find the number of comas
  parsing  = parsing0[0].split(',')

  if len(parsing) == 2:
    ninput+=1
    input_files.append(parsing[0].strip())
    test_names.append(parsing[1].strip())
    testinfo.append([])
    restarting.append(False)

  if len(parsing) == 3:
    ninput+=1
    input_files.append(parsing[0].strip())
    test_names.append(parsing[1].strip())
    testinfo.append([])
    restarting.append(True)

  elif len(parsing) == 4:
    testinfo[ninput-1].append(parsing)

ftestsuite.close()

if len(selected_input_file) == 0:
  ninput2 = ninput
else:
  ninput2 = 1

  if not '.in' in selected_input_file:
    selected_input_file = selected_input_file + '.in'

  if '/' in selected_input_file:
    parsing  = selected_input_file.split('/')
    selected_input_file = parsing[len(parsing)-1]

  if not selected_input_file in input_files:
    print('Input file name:',selected_input_file,'not present in the test suite')
    sys.exit(1)


print('Input files to be run:',ninput2)


###################################
tmpfolder='tmp_'+today

try:
  os.mkdir(tmpfolder)
except OSError:
  pass
#  print('Temporary folder already exists: '+tmpfolder)
###################################

success = 0
tested = 0

fdiff = open(tmpfolder+'/diff', 'w')
fdiff.write('#  test index          calculated                   reference                   difference        test status \n')

for iinput in range(ninput):

  if len(selected_input_file) != 0 and selected_input_file != input_files[iinput]:
    continue

  inp     = input_files[iinput]
  out     = input_files[iinput]+'_out'
  restart = restarting[iinput]

  print('\nRunning test file: '+inp)
  print(test_names[iinput])
  
  clean_run(inp,out,restart)

  check_output(out,testinfo[iinput])


fdiff.close()

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
