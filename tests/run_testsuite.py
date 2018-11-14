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
run_tddft = False

selected_input_files= []
excluded_input_files= []
mpirun=''
nprocs=1
debug=False

in_timing_section=True
sections_separator="--- Timings in (s) and # of calls ---"

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
    try:
      os.remove('RESTART_TDDFT')
    except FileNotFoundError:
      pass
    try:
      os.remove('EIGVEC_CI_0')
    except FileNotFoundError:
      pass
    try:
      os.remove('EIGVEC_CI_P')
    except FileNotFoundError:
      pass
    try:
      os.remove('EIGVEC_CI_M')
    except FileNotFoundError:
      pass
  fout = open(out, 'w')
  if len(mpirun) < 1:
    subprocess.call(['../../molgw',inp],stdout=fout,stderr=subprocess.STDOUT)
  else:
    subprocess.call([mpirun,'-n',str(nprocs),'-oversubscribe','../../molgw',inp],stdout=fout,stderr=subprocess.STDOUT)
  fout.close()
  os.chdir('..')


###################################
def check_output(out,testinfo):
  global success,tested,test_files_skipped

  #
  # First check if the test was aborted because of some limitation at compilation
  #
  for line in open(tmpfolder+'/'+out,'r').readlines():
    if 'one CPU only' in line:
      print('Test not functional in parallel => skip test')
      test_files_skipped += 1
      return
    if 'Need to compile MOLGW with HAVE_LIBINT_ONEBODY' in line:
      print('Test not functional without gradients => skip test')
      test_files_skipped += 1
      return
    if  'Angular momentum is too high' in line:
      print('LIBINT installation does not have the needed high angular momenta => skip test')
      test_files_skipped += 1
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

    in_timing_section=True
    for line in reversed(open(tmpfolder+'/'+out,'r').readlines()):
      if sections_separator in line:
        in_timing_section=False
      if key in line and not in_timing_section:
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

option_list = ['--keep','--np','--mpirun','--input','--exclude','--debug','--tddft']

if len(sys.argv) > 1:
  if '--help' in sys.argv:
    print('Run the complete test suite of MOLGW')
    print('  --keep             Keep the temporary folder')
    print('  --tddft            Run tddft tests also')
    print('  --np     n         Set the number of MPI threads to n')
    print('  --mpirun launcher  Set the MPI launcher name')
    print('  --input files      Only run these input files')
    print('  --exclude files    Run all input files but these ones')
    print('  --debug            Output debug information for this script')
    sys.exit(0)

  for argument in sys.argv:
    if '--' in argument and  argument not in option_list:
      print('Unknown option: ' + argument)
      sys.exit(1)

  if '--keep' in sys.argv or '-keep' in sys.argv:
    keeptmp = True

  if '--tddft' in sys.argv or '-tddft' in sys.argv:
    run_tddft = True

  if '--np' in sys.argv:
    i = sys.argv.index('--np') + 1
    nprocs = int( sys.argv[i] )
    mpirun='mpirun'

  if '--mpirun' in sys.argv:
    i = sys.argv.index('--mpirun') + 1
    mpirun = sys.argv[i]

  if '--input' in sys.argv:
    i = sys.argv.index('--input') + 1
    for j in range(i,len(sys.argv)):
      if '--' in sys.argv[j]:
        break
      selected_input_files.append(sys.argv[j])

  if '--exclude' in sys.argv:
    i = sys.argv.index('--exclude') + 1
    for j in range(i,len(sys.argv)):
      if '--' in sys.argv[j]:
        break
      excluded_input_files.append(sys.argv[j])

  if '--debug' in sys.argv:
    debug = True

  if len(selected_input_files) * len(excluded_input_files) > 0:
    print('--input and --exclude options are mutually exclusive. Select the one you really want.')
    sys.exit(1)



print('\n===============================')
print('Starting MOLGW test suite\n')

if len(mpirun) < 1:
  if( nprocs > 1 ):
    print('No MPI launcher has been provided. Set the number of MPI threads back to 1')
  nprocs = 1

if not os.path.isfile('../molgw') :
  print('molgw executable not found!\nMay be you should compile it first? May be you moved it around?')
  sys.exit(1)

try:
  ncores = int(os.environ['OMP_NUM_THREADS'])
except:
  ncores = 1
  pass
print('Running with \033[91m\033[1m{:3d}\033[0m MPI    threads'.format(nprocs))
print('Running with \033[91m\033[1m{:3d}\033[0m OPENMP threads'.format(ncores))
print()


###################################
# Create the temporary folder
###################################
tmpfolder='tmp_'+today

try:
  os.mkdir(tmpfolder)
except OSError:
  pass


###################################
# Run the fake.in input to get MOLGW compilation options
###################################
clean_run('fake.in','fake.out',False)
#ffake = open(tmpfolder+'/fake.out','r')
#for line in ffake:
#ffake.close()
have_openmp           = 'Running with OPENMP' in open(tmpfolder+'/fake.out').read()
have_libxc            = 'Running with LIBXC' in open(tmpfolder+'/fake.out').read()
have_mpi              = 'Running with MPI' in open(tmpfolder+'/fake.out').read()
have_scalapack        = 'Running with SCALAPACK' in open(tmpfolder+'/fake.out').read()
have_libint_onebody   = 'Running with external LIBINT calculation of the one-body operators' in open(tmpfolder+'/fake.out').read()
have_libint_gradients = 'Running with external LIBINT calculation of the gradients of the integrals' in open(tmpfolder+'/fake.out').read()
with open(tmpfolder+'/fake.out','r') as ffake:
  for line in ffake:
    if 'Perform diagonalizations with LAPACK routines' in line:
      lapack_diago_flavor = line.split(':')[1].strip()
print('MOLGW compilation details:')
print('                   OPENMP: {}'.format(have_openmp) )
print('                      MPI: {}'.format(have_mpi) )
print('                SCALAPACK: {}'.format(have_scalapack) )
print('                    LIBXC: {}'.format(have_libxc) )
print('            1-body LIBINT: {}'.format(have_libint_onebody) )
print('         gradients LIBINT: {}'.format(have_libint_gradients) )
print('        (SCA)LAPACK diago: {}'.format(lapack_diago_flavor) )
print()

os.remove(tmpfolder+'/fake.out')



###################################
# Parse the file testsuite
###################################
ninput = 0
input_files    = []
restarting     = []
parallel       = []
need_scalapack = []
tddft          = []
test_names     = []
testinfo       = []

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
    parallel.append(True)
    need_scalapack.append(False)
    tddft.append(False)

  if len(parsing) == 3:
    ninput+=1
    input_files.append(parsing[0].strip())
    test_names.append(parsing[1].strip())
    testinfo.append([])
    if 'restart' in parsing[2].lower():
      restarting.append(True)
    else:
      restarting.append(False)
    if 'noparallel' in parsing[2].lower():
      parallel.append(False)
    else:
      parallel.append(True)
    need_scalapack.append( 'need_scalapack' in parsing[2].lower() )
    if 'tddft' in parsing[2].lower():
      tddft.append(True)
    else:
      tddft.append(False)

  elif len(parsing) == 4:
    testinfo[ninput-1].append(parsing)

ftestsuite.close()


###################################
# Check the selection and the exclusion lists
###################################
if len(selected_input_files) == 0 and len(excluded_input_files) == 0:
  ninput2 = ninput

elif len(selected_input_files) > 0:
  ninput2 = len(selected_input_files)

  for i in range(len(selected_input_files)):
    if not '.in' in selected_input_files[i]:
      selected_input_files[i] = selected_input_files[i] + '.in'

    if '/' in selected_input_files[i]:
      parsing  = selected_input_files[i].split('/')
      selected_input_files[i] = parsing[len(parsing)-1]

    if not selected_input_files[i] in input_files:
      print('Input file name:',selected_input_files[i],'not present in the test suite')
      sys.exit(1)

else:
  ninput2 = ninput - len(excluded_input_files)

  for i in range(len(excluded_input_files)):
    if not '.in' in excluded_input_files[i]:
      excluded_input_files[i] = excluded_input_files[i] + '.in'

    if '/' in excluded_input_files[i]:
      parsing  = excluded_input_files[i].split('/')
      excluded_input_files[i] = parsing[len(parsing)-1]

    if not excluded_input_files[i] in input_files:
      print('Input file name:',excluded_input_files[i],'not present in the test suite')
      sys.exit(1)



print('Input files found in the test suite: {}'.format(ninput2))


###################################

success            = 0
tested             = 0
test_files_skipped = 0

fdiff = open(tmpfolder+'/diff', 'w')
fdiff.write('#  test index          calculated                   reference                   difference        test status \n')

for iinput in range(ninput):

  if len(selected_input_files) != 0 and not input_files[iinput] in selected_input_files:
    continue
  if len(excluded_input_files) != 0 and input_files[iinput] in excluded_input_files:
    continue
  if need_scalapack[iinput]:
    test_files_skipped += 1
    print('\nSkipping test file: '+inp)
    print('  because this compilation of MOLGW does not have SCALAPACK')
    continue
  if not parallel[iinput] and nprocs > 1:
    print('\nSkipping test file: '+inp)
    print('  because this test is only serial')
    test_files_skipped += 1
    continue
  if not run_tddft and tddft[iinput]:
    print('\nSkipping test file: '+inp)
    print('  because the RT-TDDFT needs to be specifically activated with --tddft')
    test_files_skipped += 1
    continue


  inp     = input_files[iinput]
  out     = input_files[iinput].split('.in')[0]+'.out'
  restart = restarting[iinput]

  print('\nRunning test file: '+inp)
  print(test_names[iinput])

  clean_run(inp,out,restart)

  check_output(out,testinfo[iinput])


fdiff.close()

print('\n\n===============================')
print('      Test Summary \n')
print('        Succesful tests:   {0:} / {1:}\n'.format(success,tested))
if test_files_skipped > 0 :
  print('  Test files not tested:    {:}\n'.format(test_files_skipped))
print('       Elapsed time (s):   ','{:.2f}'.format(time.time() - start_time) )
print('===============================\n')


###################################
# Erase the tmp folder by default
if not keeptmp:
  shutil.rmtree(tmpfolder)


###################################
