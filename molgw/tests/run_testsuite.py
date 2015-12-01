#!/usr/bin/python3

###################################
# Running the MOLGW test suite
###################################


import sys, os, time, shutil, subprocess

today=time.strftime("%Y")+'_'+time.strftime("%m")+'_'+time.strftime("%d")
start_time = time.time()
keeptmp = False

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
  subprocess.call(['../../molgw',inp],stdout=fout)
  fout.close()
  os.chdir('..')


###################################
def check_output(out,testinfo):
  global success,tested
  for itest in range(0,len(testinfo)):
    tested += 1
    key = testinfo[itest][0].strip()
    ref = float(testinfo[itest][1])
    pos = int(  testinfo[itest][2])
    tol = float(testinfo[itest][3])

    for line in reversed(open(tmpfolder+'/'+out,'r').readlines()):
      if key in line:
        parsing  = line.split(':')
        parsing2 = parsing[1].split()

        fdiff.write(str(tested).rjust(6) + parsing2[pos].rjust(30)+str(ref).rjust(30)+str(float(parsing2[pos]) - ref).rjust(30)+'\n')

        if abs( float(parsing2[pos]) - ref ) < tol:
          print(key.rjust(30)+'[ \033[92m\033[1mOK\033[0m ]'.rjust(30))
          success += 1
          break
        else:
          print(key.rjust(30)+'[\033[91m\033[1mFAIL\033[0m]'.rjust(30))
          break
     
###################################
# Parse the command line

if len(sys.argv) > 1:
  if '--help' in sys.argv[1]:
    print('Run the complete test suite of MOLGW')
    print('--keep  Keep the temporary folder')
    sys.exit(0)
  if '--keep' in sys.argv[1]:
    keeptmp = True




print('\n===============================')
print('Starting MOLGW test suite\n')

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

print('Input files to be run:',ninput)


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

for iinput in range(0,ninput):

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
