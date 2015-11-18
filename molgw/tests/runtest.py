#!/usr/bin/python3

###################################
# Running the MOLGW test suite
###################################


import sys, os, time, shutil, subprocess

today=time.strftime("%Y")+'_'+time.strftime("%m")+'_'+time.strftime("%d")
start_time = time.time()
keeptmp = False

###################################
def clean_run(inp,out):
  shutil.copy('inputs/'+inp,tmpfolder+'/'+inp)
  os.chdir(tmpfolder)
  fout = open(out, 'w')
  try:
    os.remove('RESTART')
  except FileNotFoundError:
    pass
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
    tol = float(testinfo[itest][2])

    for line in reversed(open(tmpfolder+'/'+out,'r').readlines()):
      if key in line:
        parsing=line.split(':')

        if abs( float(parsing[1]) - ref ) < tol:
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
test_names = []
testinfo = []

ftestsuite = open('inputs/testsuite','r')
for line in ftestsuite:
  parsing=line.split(',')

  if len(parsing) == 2:
    ninput+=1
    input_files.append(parsing[0].strip())
    test_names.append(parsing[1].strip())
    testinfo.append([])

  elif len(parsing) == 3:
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

for iinput in range(0,ninput):

  inp = input_files[iinput]
  out = input_files[iinput]+'_out'
  print('\nRunning test file: '+inp)
  print(test_names[iinput])
  
  clean_run(inp,out)

  check_output(out,testinfo[iinput])

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
