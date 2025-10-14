#!/usr/bin/python3

###################################
# This file is part of MOLGW
# Author: Fabien Bruneval


###################################
# Running the MOLGW test suite
###################################


import sys, os, time, shutil, subprocess
import re

today = time.strftime("%Y") + '_' + time.strftime("%m") + '_' + time.strftime("%d")
start_time = time.time()
keeptmp = False

selected_input_files= []
excluded_input_files= []
input_param_selection= []
mpirun = ''
nprocs = 1
ncores = 1
debug = False
listing = False

in_timing_section = True
sections_separator = "--- Timings in (s) and # of calls ---"

def extract_between_quotes(s):
    match = re.search(r'"(.*?)"', s)
    return match.group(1) if match else ''

###################################
def clean_run(inp, out, restart, command=""):
  shutil.copy('inputs/' + inp, tmpfolder + '/' + inp)
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

  if len(command) > 0:
    #result = subprocess.run(command.split(), capture_output=True, text=True)
    result = subprocess.run(command.split())

  fout = open(out, 'w')
  if len(mpirun) < 1:
    subprocess.call(['../../molgw', inp], stdout=fout, stderr=subprocess.STDOUT)
  else:
    # mpirun from openmpi may need '-oversubscribe'
    subprocess.call(mpirun.split() + ['-n', str(nprocs), '../../molgw', inp], stdout=fout, stderr=subprocess.STDOUT)
  fout.close()

  with open(out, 'r') as fout:
      functional = "Welcome to the fascinating world of MOLGW" in fout.read()

  os.chdir('..')
  return functional
  


###################################
def check_output(out, testinfo):
  global success, tested, test_files_skipped, test_files_success, skipping_reason

  #
  # First check if the test was aborted because of some limitation at compilation
  #
  for line in open(tmpfolder + '/' + out, 'r').readlines():
    if  'Angular momentum is too high' in line:
      test_files_skipped += 1
      print('LIBINT or LIBCINT installation does not have the high enough angular momenta => skip test')
      skipping_reason.append('LIBINT or LIBCINT installation does not have high enough angular momenta')
      return
  #
  # Second check if there is a memory leak
  #
  key = '    Memory ('
  ref = 0.000
  tol = 0.001
  key_found = False
  tested += 1
  success_in_this_file = 0
  for line in reversed(open(tmpfolder + '/' + out, 'r').readlines()):
    if key in line:
      key_found = True
      parsing  = line.split(':')
      parsing2 = parsing[1].split()
      if abs( float(parsing2[0]) - ref ) < tol:
        print('No memory leak'.rjust(30) + '[ \033[92m\033[1mOK\033[0m ]'.rjust(30))
        success += 1
        success_in_this_file += 1
        fdiff.write(str(tested).rjust(6) + "Memory leak".rjust(30) + parsing2[0].rjust(30) \
              + str(ref).rjust(30) + str(float(parsing2[0]) - ref).rjust(30) + '  OK  \n')
        break
      else:
        print('No memory leak'.rjust(30) + '[\033[91m\033[1mFAIL\033[0m]'.rjust(30))
        fdiff.write(str(tested).rjust(6) + "Memory leak".rjust(30) + parsing2[0].rjust(30) \
              + str(ref).rjust(30) + str(float(parsing2[0]) - ref).rjust(30) + ' FAIL \n')
        break

  if not key_found:
    print('No memory leak'.rjust(30) + '[\033[91m\033[1mNOT FOUND\033[0m]'.rjust(30))

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
      print('open file: ' + tmpfolder + '/' + out)
      print('===end debug')

    key_found = False

    in_timing_section = True
    for line in reversed(open(tmpfolder + '/' + out, 'r').readlines()):
      if sections_separator in line:
        in_timing_section = False
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
          print(key.rjust(30) + '[ \033[92m\033[1mOK\033[0m ]'.rjust(30))
          success += 1
          success_in_this_file += 1
          fdiff.write(str(tested).rjust(6) + key.rjust(30) + parsing2[pos].rjust(30) \
                + str(ref).rjust(30) + str(float(parsing2[pos]) - ref).rjust(30) + '  OK  \n')
          break
        else:
          print(key.rjust(30) + '[\033[91m\033[1mFAIL\033[0m]'.rjust(30))
          fdiff.write(str(tested).rjust(6) + key.rjust(30) + parsing2[pos].rjust(30) \
                + str(ref).rjust(30) + str(float(parsing2[pos]) - ref).rjust(30) + ' FAIL \n')
          break
    if not key_found:
      print(key.rjust(30) + '[\033[91m\033[1mNOT FOUND\033[0m]'.rjust(30))

  failures_in_this_file = len(testinfo) + 1 - success_in_this_file

  if failures_in_this_file == 0:
     test_files_success += 1

  return failures_in_this_file

###################################
# Parse the command line

option_list = ['--keep', '--np', '--nc', '--mpirun', '--input', '--exclude', '--input-parameter', '--debug', '--list']

if len(sys.argv) > 1:
  if '--help' in sys.argv:
    print('Run the complete test suite of MOLGW')
    print('  --keep             Keep the temporary folder')
    print('  --np     n         Set the number of MPI processes to n')
    print('  --nc     n         Set the number of OPENMP threads to n')
    print('  --mpirun launcher  Set the MPI launcher name')
    print('  --input files      Only run these input files')
    print('  --list             list all the input files')
    print('  --exclude files    Run all input files but these ones')
    print('  --input-parameter  Only run input files that contain this input parameter. Example:')
    print('                     --input-parameter scf = \'LDA\' ')
    print('  --debug            Output debug information for this script')
    sys.exit(0)

  for argument in sys.argv:
    if '--' in argument and  argument not in option_list:
      print('Unknown option: ' + argument)
      sys.exit(1)

  if '--keep' in sys.argv or '-keep' in sys.argv:
    keeptmp = True

  if '--np' in sys.argv:
    i = sys.argv.index('--np') + 1
    nprocs = int( sys.argv[i] )
    mpirun = 'mpirun'

  if '--mpirun' in sys.argv:
    i = sys.argv.index('--mpirun') + 1
    mpirun = sys.argv[i]

  if '--nc' in sys.argv:
    i = sys.argv.index('--nc') + 1
    ncores = int( sys.argv[i] )

  if '--input' in sys.argv:
    i = sys.argv.index('--input') + 1
    for j in range(i, len(sys.argv)):
      if '--' in sys.argv[j]:
        break
      selected_input_files.append(sys.argv[j])

  if '--input-parameter' in sys.argv:
    i = sys.argv.index('--input-parameter') + 1
    for j in range(i, len(sys.argv)):
      if '--' in sys.argv[j]:
        break
      input_param_selection.append(sys.argv[j])

  if '--exclude' in sys.argv:
    i = sys.argv.index('--exclude') + 1
    for j in range(i, len(sys.argv)):
      if '--' in sys.argv[j]:
        break
      excluded_input_files.append(sys.argv[j])

  if '--debug' in sys.argv:
    debug = True

  if '--list' in sys.argv:
    listing = True

  if len(selected_input_files) * len(excluded_input_files) > 0:
    print('--input and --exclude options are mutually exclusive. Select the one you really want.')
    sys.exit(1)



print('\n===============================')
print('Starting MOLGW test suite\n')

if len(mpirun) < 1:
  if( nprocs > 1 ):
    print('No MPI launcher has been provided. Set the number of MPI processes back to 1')
  nprocs = 1

if not os.path.isfile('../molgw') :
  print('molgw executable not found!\nMay be you should compile it first? May be you moved it around?')
  sys.exit(1)

if ncores > 1:
  os.environ["OMP_NUM_THREADS"]      = str(ncores)
  os.environ["MKL_NUM_THREADS"]      = str(ncores)
  os.environ["OPENBLAS_NUM_THREADS"] = str(ncores)
  os.environ["OMP_STACKSIZE"]   = "128M"
else:
  os.environ["OMP_NUM_THREADS"] = '1'
  os.environ["MKL_NUM_THREADS"] = '1'

try:
  ncores = int(os.environ['OMP_NUM_THREADS'])
except:
  ncores = 1
  pass
print('Running with \033[91m\033[1m{:3d}\033[0m MPI  processes'.format(nprocs))
print('Running with \033[91m\033[1m{:3d}\033[0m OPENMP threads'.format(ncores))
print()


###################################
# Create the temporary folder
###################################
tmpfolder = 'tmp'

try:
  os.mkdir(tmpfolder)
except OSError:
  pass

###################################
# Parse molgw.h to obtain MOLGW version
###################################
with open('../src/molgw.h', 'r') as stream:
    version = ''
    for line in stream:
        words = line.split()
        if len(words) > 1:
            if words[0] == '#define' and words[1] == 'MOLGW_VERSION':
                version = words[2].replace('\"', '').replace('\'', '')

###################################
# Run the fake.in input to get MOLGW compilation options
###################################
molgw_executable_functional = clean_run('fake.in', 'fake.out', False)

if not molgw_executable_functional:
    print("MOLGW executable is not functional")
    print("Dump last output:")
    with open(tmpfolder + '/fake.out', 'r') as f:
        print(f.read())
    sys.exit("MOLGW executable is not functional")

have_openmp           = 'Running with OPENMP' in open(tmpfolder + '/fake.out').read()
have_libxc            = 'Running with LIBXC' in open(tmpfolder + '/fake.out').read()
have_mpi              = 'Running with MPI' in open(tmpfolder + '/fake.out').read()
have_scalapack        = 'Running with SCALAPACK' in open(tmpfolder + '/fake.out').read()
have_onebody          = 'Running with external LIBINT or LIBCINT calculation of the one-body operators' in open(tmpfolder + '/fake.out').read()
have_gradients        = 'Running with external LIBINT calculation of the gradients of the one-body integrals' in open(tmpfolder + '/fake.out').read() \
                        or 'Code compiled with LIBCINT' in open(tmpfolder + '/fake.out').read()
#have_libint_forces    = 'Running with external LIBINT calculation of the gradients of the Coulomb integrals' in open(tmpfolder + '/fake.out').read()
have_libint_forces    = False  # FIXME force calculation is broken as of today
is_libcint            = 'Code compiled with LIBCINT support' in open(tmpfolder + '/fake.out').read()
have_hdf5             = 'Running with HDF5' in open(tmpfolder + '/fake.out').read()

#with open(tmpfolder + '/fake.out', 'r') as ffake:
#  for line in ffake:
#    if 'Perform diagonalizations with (Sca)LAPACK routines' in line:
#      lapack_diago_flavor = line.split(':')[1].strip()

print('MOLGW compilation details:')
print('              MOLGW version: ' + version)
print('                     OPENMP: {}'.format(have_openmp) )
print('                        MPI: {}'.format(have_mpi) )
print('                  SCALAPACK: {}'.format(have_scalapack) )
print('                      LIBXC: {}'.format(have_libxc) )
print('                       HDF5: {}'.format(have_hdf5) )
if is_libcint:
    print('                  integrals: LIBCINT' )
else:
    print('                  integrals: LIBINT' )
print('           1-body integrals: {}'.format(have_onebody) )
print('        gradients integrals: {}'.format(have_gradients) )
print(' LIBINT gradients integrals: {}'.format(have_libint_forces) )
#print('        (Sca)LAPACK diago: {}'.format(lapack_diago_flavor) )
print()

#os.remove(tmpfolder + '/fake.out')


# The default OMP stacksize is generally too small:
if have_openmp:
  os.environ["OMP_STACKSIZE"]   = '32M'

###################################
# Parse the file testsuite
###################################
ninput = 0
input_files    = []
restarting     = []
parallel       = []
need_scalapack = []
need_gradients = []
need_forces    = []
need_libcint   = []
test_names     = []
testinfo       = []
command        = []

ftestsuite = open('inputs/testsuite', 'r')
for line in ftestsuite:
  # Removing the comments
  parsing0 = line.split('#')
  # Find the number of comas
  parsing  = parsing0[0].split(', ')

  if len(parsing) == 2:
    ninput +=1
    input_files.append(parsing[0].strip())
    test_names.append(parsing[1].strip())
    testinfo.append([])
    restarting.append(False)
    parallel.append(True)
    need_scalapack.append(False)
    need_gradients.append(False)
    need_forces.append(False)
    need_libcint.append(False)
    command.append("")

  if len(parsing) == 3:
    ninput +=1
    input_files.append(parsing[0].strip())
    test_names.append(parsing[1].strip())
    testinfo.append([])
    command.append(extract_between_quotes(parsing[2]))

    if 'restart' in parsing[2].lower():
      restarting.append(True)
    else:
      restarting.append(False)
    if 'noparallel' in parsing[2].lower():
      parallel.append(False)
    else:
      parallel.append(True)
    need_scalapack.append( 'need_scalapack' in parsing[2].lower() )
    need_gradients.append( 'need_gradients' in parsing[2].lower() )
    need_forces.append( 'need_forces' in parsing[2].lower() )
    need_libcint.append( 'need_libcint' in parsing[2].lower() )

  elif len(parsing) == 4:
    testinfo[ninput-1].append(parsing)

ftestsuite.close()

###################################
# List all the input files and exit
###################################
if listing:
    print('=== List of input files in the Suite ===')
    for i, inpfile in enumerate(input_files):
        print('{:04}: {}'.format(i + 1, inpfile))
    print('========================================')
    sys.exit(0)


###################################
# Check the selection by input variable value
###################################
if len(input_param_selection) > 0:
  print('\nUser asked for a specific subset of input files containing:')
  key1 = input_param_selection[0].split('=')[0]
  key2 = input_param_selection[-1].split('=')[-1]
  print('    ' + key1 + ' = ' + key2 + '\n')

  for iinput in range(ninput):
    inp = input_files[iinput]
  
    present = False

    fin = open('inputs/' + inp, 'r')
    for line in fin:
      if key1.lower() in line.lower() and key2.lower() in line.lower():
        present = True
    fin.close()

    if present:
      selected_input_files.append(inp)
  if len(selected_input_files) == 0:
    print('User selected an input parameter or a value that is not present in any input file')
    sys.exit(1)


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
      print('Input file name:', selected_input_files[i], 'not present in the test suite')
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
      print('Input file name:', excluded_input_files[i], 'not present in the test suite')
      sys.exit(1)


  


print('Input files to be executed: {}'.format(ninput2))


###################################

success            = 0
tested             = 0
test_files_skipped = 0
test_files_success = 0
skipping_reason    = []

fdiff = open(tmpfolder + '/diff', 'w')
fdiff.write('#   test index       property tested                     calculated                 reference                 difference        test status \n')
ffailed = open(tmpfolder + '/failed_tests', 'w')

for iinput in range(ninput):

  if len(selected_input_files) != 0 and not input_files[iinput] in selected_input_files:
    continue
  if len(excluded_input_files) != 0 and input_files[iinput] in excluded_input_files:
    continue

  inp     = input_files[iinput]
  out     = input_files[iinput].split('.in')[0] + '.out'
  restart = restarting[iinput]

  if need_libcint[iinput] and not is_libcint:
    test_files_skipped += 1
    print('\nSkipping test file: ' + inp)
    print('  because this compilation of MOLGW does not have LIBCINT')
    skipping_reason.append('this compilation of MOLGW does not have LIBCINT')
    continue
  if need_scalapack[iinput] and not have_scalapack:
    test_files_skipped += 1
    print('\nSkipping test file: ' + inp)
    print('  because this compilation of MOLGW does not have SCALAPACK')
    skipping_reason.append('this compilation of MOLGW does not have SCALAPACK')
    continue
  if need_gradients[iinput] and not have_gradients:
    test_files_skipped += 1
    print('\nSkipping test file: ' + inp)
    print('  because this compilation of MOLGW does not have the integral gradients')
    skipping_reason.append('this compilation of MOLGW does not have the integral gradients')
    continue
  if need_forces[iinput] and not have_libint_forces:
    test_files_skipped += 1
    print('\nSkipping test file: ' + inp)
    print('  because this compilation of MOLGW does not have the force integrals')
    skipping_reason.append('this compilation of MOLGW does not have the force integrals')
    continue
  if not parallel[iinput] and nprocs > 1:
    test_files_skipped += 1
    print('\nSkipping test file: ' + inp)
    print('  because this test is only serial')
    skipping_reason.append('this test is only serial')
    continue


  print('\nRunning test file: ' + inp)
  print(test_names[iinput])
  fdiff.write("# " + inp + "\n")

  molgw_executable_functional = clean_run(inp, out, restart, command=command[iinput])
  
  failures = check_output(out, testinfo[iinput])
  if failures != 0:
    ffailed.write(out + "\n")

    # displays the content of the .out file associated with the failed test
    output_path = os.path.join(tmpfolder, out)
    print(f"\n===== TEST failure : {inp} =====")
    print(f"----- Displaying the content of {out} associated with the failed test -----\n")
    with open(output_path, "r", encoding="utf-8", errors="replace") as f:
        print(f.read())
    print("=================================\n")


fdiff.close()
ffailed.close()

print('\n\n===============================')
print('      Test Suite Summary \n')
print('      Test files tested:   {:4d} / {:4d}\n'.format(ninput2-test_files_skipped, ninput2))
if success == tested:
  print('     Test files success:   \033[92m\033[1m{:4d} / {:4d}\033[0m  '.format(test_files_success, ninput2-test_files_skipped))
  print('       Successful tests:   \033[92m\033[1m{:4d} / {:4d}\033[0m\n'.format(success, tested))
else:
  print('     Test files success:   \033[91m\033[1m{:4d} / {:4d}\033[0m  '.format(test_files_success, ninput2-test_files_skipped))
  print('       Successful tests:   \033[91m\033[1m{:4d} / {:4d}\033[0m\n'.format(success, tested))
print('       Elapsed time (s):   ', '{:.2f}'.format(time.time() - start_time) )
print('===============================\n')
if test_files_skipped > 0 :
  print(' Some tests have been skipped for the following reasons:')
  for reason in list(set(skipping_reason)):
    ireason = skipping_reason.count(reason)
    print('   * {:<80}  ({:=4d} tests)'.format(reason, ireason))
  print('===============================\n')


###################################
# Erase the tmp folder by default
if not keeptmp:
  shutil.rmtree(tmpfolder)

sys.exit(abs(success-tested))
###################################
