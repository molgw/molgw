#!/usr/bin/env python3

###################################
# This file is part of MOLGW
# Author: Fabien Bruneval
###################################
# Running the MOLGW test suite
###################################

import sys
import os
import time
import shutil
import subprocess
import re
import argparse
from dataclasses import dataclass, field
from pathlib import Path

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class TestCase:
    input_file: str
    name: str
    checks: list = field(default_factory=list)
    restart: bool = False
    parallel: bool = True
    need_scalapack: bool = False
    need_gradients: bool = False
    need_forces: bool = False
    need_libcint: bool = False
    need_fftw3: bool = False
    command: str = ""

    @property
    def output_file(self):
        return self.input_file.replace('.in', '.out')


@dataclass
class RunStats:
    success: int = 0
    tested: int = 0
    files_skipped: int = 0
    files_success: int = 0
    skipping_reasons: list = field(default_factory=list)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TODAY = time.strftime("%Y_%m_%d")
SECTIONS_SEPARATOR = "--- Timings in (s) and # of calls ---"
FILES_TO_CLEAN = [
    'RESTART', 'ENERGY_QP', 'SCREENED_COULOMB',
    'RESTART_TDDFT', 'EIGVEC_CI_0', 'EIGVEC_CI_P', 'EIGVEC_CI_M',
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def extract_between_quotes(s: str) -> str:
    match = re.search(r'"(.*?)"', s)
    return match.group(1) if match else ''


def read_fake_out(path: Path) -> str:
    return path.read_text(encoding='utf-8', errors='replace')


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def clean_run(test: TestCase, tmpfolder: Path, mpirun: str, nprocs: int, debug: bool) -> bool:
    """Copy input file to tmpfolder, optionally remove restart files, then run molgw."""
    shutil.copy(Path('inputs') / test.input_file, tmpfolder / test.input_file)

    if not test.restart:
        for fname in FILES_TO_CLEAN:
            (tmpfolder / fname).unlink(missing_ok=True)

    if test.command:
        subprocess.run(test.command.split(), check=False)

    out_path = tmpfolder / test.output_file
    molgw_bin = str(Path('../../molgw'))
    cmd = (mpirun.split() + ['-n', str(nprocs), molgw_bin, test.input_file]
           if mpirun else [molgw_bin, test.input_file])

    with out_path.open('w') as fout:
        subprocess.run(cmd, stdout=fout, stderr=subprocess.STDOUT, cwd=tmpfolder)

    content = out_path.read_text(encoding='utf-8', errors='replace')
    return "Welcome to the fascinating world of MOLGW" in content


def check_output(test: TestCase, tmpfolder: Path, stats: RunStats,
                 fdiff, ffailed, verbose: bool, debug: bool) -> int:
    """Parse the output file and compare against reference values."""
    out_path = tmpfolder / test.output_file
    lines = out_path.read_text(encoding='utf-8', errors='replace').splitlines()
    reversed_lines = list(reversed(lines))

    # --- Skip check: angular momentum ---
    if any('Angular momentum is too high' in ln for ln in lines):
        stats.files_skipped += 1
        print('LIBINT or LIBCINT installation does not have the high enough angular momenta => skip test')
        stats.skipping_reasons.append(
            'LIBINT or LIBCINT installation does not have high enough angular momenta')
        return 0

    stats.tested += 1
    success_in_file = 0

    # --- Memory-leak check ---
    mem_key = '    Memory ('
    mem_ref = 0.0
    mem_tol = 0.001
    mem_found = False

    for ln in reversed_lines:
        if mem_key in ln:
            mem_found = True
            val = float(ln.split(':')[1].split()[0])
            diff = val - mem_ref
            if abs(diff) < mem_tol:
                print('No memory leak'.rjust(30) + '[ \033[92m\033[1mOK\033[0m ]'.rjust(30))
                stats.success += 1
                success_in_file += 1
                fdiff.write(f"{stats.tested:6d}{'Memory leak':>30}{val:>30}{mem_ref:>30}{diff:>30}  OK  \n")
            else:
                print('No memory leak'.rjust(30) + '[\033[91m\033[1mFAIL\033[0m]'.rjust(30))
                fdiff.write(f"{stats.tested:6d}{'Memory leak':>30}{val:>30}{mem_ref:>30}{diff:>30} FAIL \n")
            break

    if not mem_found:
        print('No memory leak'.rjust(30) + '[\033[91m\033[1mNOT FOUND\033[0m]'.rjust(30))

    # --- Per-quantity checks ---
    for check in test.checks:
        stats.tested += 1
        key  = check[0].strip()
        ref  = float(check[1])
        pos  = int(check[2])
        tol  = float(check[3])

        if debug:
            print(f'===debug: looking for key="{key}" in {out_path}')

        key_found = False
        in_timing  = True  # start True; becomes False once separator seen

        for ln in reversed_lines:
            if SECTIONS_SEPARATOR in ln:
                in_timing = False
            if key in ln and not in_timing:
                key_found = True
                parts = ln.split(':')[1].split()
                if debug:
                    print(f'===debug: parts={parts}')
                val  = float(parts[pos])
                diff = val - ref

                if abs(diff) < tol:
                    print(key.rjust(30) + '[ \033[92m\033[1mOK\033[0m ]'.rjust(30))
                    stats.success += 1
                    success_in_file += 1
                    fdiff.write(f"{stats.tested:6d}{key:>30}{val:>30}{ref:>30}{diff:>30}  OK  \n")
                else:
                    print(key.rjust(30) + '[\033[91m\033[1mFAIL\033[0m]'.rjust(30))
                    fdiff.write(f"{stats.tested:6d}{key:>30}{val:>30}{ref:>30}{diff:>30} FAIL \n")
                break

        if not key_found:
            print(key.rjust(30) + '[\033[91m\033[1mNOT FOUND\033[0m]'.rjust(30))

    failures_in_file = len(test.checks) + 1 - success_in_file

    if failures_in_file == 0:
        stats.files_success += 1
    else:
        ffailed.write(test.output_file + '\n')
        if verbose:
            print(f"\n===== TEST failure : {test.input_file} =====")
            print(f"----- Displaying the content of {test.output_file} -----\n")
            print(out_path.read_text(encoding='utf-8', errors='replace'))
            print("=================================\n")

    return failures_in_file


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Run the complete test suite of MOLGW',
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument('-k', '--keep', action='store_true',
                        help='Keep the temporary folder')
    parser.add_argument('-np', '--np', type=int, default=1, metavar='N',
                        help='Set the number of MPI processes to N')
    parser.add_argument('-nc', '--nc', type=int, default=1, metavar='N',
                        help='Set the number of OpenMP threads to N')
    parser.add_argument('-mpirun', '--mpirun', default='', metavar='LAUNCHER',
                        help='Set the MPI launcher name')
    parser.add_argument('-i', '--input', nargs='+', default=[], metavar='FILE',
                        help='Only run these input files')
    parser.add_argument('-e', '--exclude', nargs='+', default=[], metavar='FILE',
                        help='Run all input files but these ones')
    parser.add_argument('-p', '--input-parameter', nargs='+', default=[], metavar='PARAM',
                        help="Only run inputs that contain this parameter.\n"
                             "Example: --input-parameter scf = 'LDA'")
    parser.add_argument('-d', '--debug', action='store_true',
                        help='Output debug information for this script')
    parser.add_argument('-l', '--list', action='store_true',
                        help='List all input files and exit')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Display the output file for any failed test')
    args = parser.parse_args()

    if args.input and args.exclude:
        parser.error('--input and --exclude are mutually exclusive.')

    return args


# ---------------------------------------------------------------------------
# Testsuite file parser
# ---------------------------------------------------------------------------

def parse_testsuite(path: Path) -> list[TestCase]:
    """Return a list of TestCase objects parsed from the testsuite file."""
    tests: list[TestCase] = []

    with path.open('r') as fts:
        for line in fts:
            # Strip comments
            core = line.split('#')[0]
            parts = [p.strip() for p in core.split(', ')]

            if len(parts) == 2:
                tests.append(TestCase(
                    input_file=parts[0],
                    name=parts[1],
                ))

            elif len(parts) == 3:
                flags = parts[2].lower()
                tests.append(TestCase(
                    input_file=parts[0],
                    name=parts[1],
                    command=extract_between_quotes(parts[2]),
                    restart='restart' in flags,
                    parallel='noparallel' not in flags,
                    need_scalapack='need_scalapack' in flags,
                    need_gradients='need_gradients' in flags,
                    need_forces='need_forces' in flags,
                    need_libcint='need_libcint' in flags,
                    need_fftw3='need_fftw3' in flags,
                ))

            elif len(parts) == 4 and tests:
                tests[-1].checks.append(parts)

    return tests


# ---------------------------------------------------------------------------
# Selection / exclusion helpers
# ---------------------------------------------------------------------------

def normalize_input_name(name: str) -> str:
    """Ensure the name ends with .in and has no directory prefix."""
    p = Path(name)
    stem = p.name
    return stem if stem.endswith('.in') else stem + '.in'


def apply_selection(tests: list[TestCase], selected: list[str],
                    excluded: list[str]) -> list[TestCase]:
    all_names = {t.input_file for t in tests}

    if selected:
        normalized = [normalize_input_name(s) for s in selected]
        for name in normalized:
            if name not in all_names:
                sys.exit(f'Input file name: {name} not present in the test suite')
        return [t for t in tests if t.input_file in normalized]

    if excluded:
        normalized = {normalize_input_name(e) for e in excluded}
        for name in normalized:
            if name not in all_names:
                sys.exit(f'Input file name: {name} not present in the test suite')
        return [t for t in tests if t.input_file not in normalized]

    return tests


def resolve_restart_dependencies(selected: list[TestCase], all_tests: list[TestCase]) -> list[TestCase]:
    """For every selected test tagged as restart=True, ensure the closest preceding
    non-restart test in the full suite is also included (it produces the files the
    restart test depends on).  Original testsuite order is preserved."""
    selected_files = {t.input_file for t in selected}
    to_add = []

    for test in selected:
        if not test.restart:
            continue
        idx = all_tests.index(test)
        # Walk backwards to find the nearest preceding non-restart test
        for predecessor in reversed(all_tests[:idx]):
            if not predecessor.restart:
                if predecessor.input_file not in selected_files:
                    to_add.append(predecessor)
                    selected_files.add(predecessor.input_file)
                    print(f'  Adding prerequisite test "{predecessor.input_file}" '
                          f'required by restart test "{test.input_file}"')
                break

    # Rebuild in original suite order
    return [t for t in all_tests if t.input_file in selected_files]


def apply_param_filter(tests: list[TestCase], param_tokens: list[str],
                       all_tests: list[TestCase]) -> list[TestCase]:
    """Keep only tests whose input file contains the key=value pair,
    then pull in any prerequisite tests needed by restart-tagged entries."""
    if not param_tokens:
        return tests

    raw = ' '.join(param_tokens)
    key1, _, key2 = raw.partition('=')
    key1 = key1.strip()
    key2 = key2.strip().strip("'\"")
    print(f'\nUser asked for a specific subset of input files containing:\n    {key1} = {key2}\n')

    selected = []
    for test in tests:
        inp_path = Path('inputs') / test.input_file
        content  = inp_path.read_text(encoding='utf-8', errors='replace').lower()
        if key1.lower() in content and key2.lower() in content:
            selected.append(test)

    if not selected:
        sys.exit('User selected an input parameter or a value that is not present in any input file')

    return resolve_restart_dependencies(selected, all_tests)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    start_time = time.time()
    args = parse_args()

    mpirun = args.mpirun
    nprocs = args.np
    ncores = args.nc

    # --np implies mpirun
    if nprocs > 1 and not mpirun:
        mpirun = 'mpirun'

    print('\n===============================')
    print('Starting MOLGW test suite\n')

    if not mpirun and nprocs > 1:
        print('No MPI launcher has been provided. Setting the number of MPI processes back to 1')
        nprocs = 1

    molgw_bin = Path('../molgw')
    if not molgw_bin.is_file():
        sys.exit('molgw executable not found!\nMay be you should compile it first? May be you moved it around?')

    # OpenMP / MKL environment
    if ncores > 1:
        for var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
            os.environ[var] = str(ncores)
        os.environ['OMP_STACKSIZE'] = '128M'
    else:
        os.environ['OMP_NUM_THREADS'] = '1'
        os.environ['MKL_NUM_THREADS'] = '1'

    try:
        ncores = int(os.environ['OMP_NUM_THREADS'])
    except (KeyError, ValueError):
        ncores = 1

    print(f'Running with \033[91m\033[1m{nprocs:3d}\033[0m MPI  processes')
    print(f'Running with \033[91m\033[1m{ncores:3d}\033[0m OPENMP threads\n')

    # Temporary folder
    tmpfolder = Path('tmp')
    tmpfolder.mkdir(exist_ok=True)

    # Read MOLGW version
    version = ''
    molgw_h = Path('../src/molgw.h')
    with molgw_h.open('r') as fh:
        for ln in fh:
            words = ln.split()
            if len(words) > 2 and words[0] == '#define' and words[1] == 'MOLGW_VERSION':
                version = words[2].strip('"\'')

    # Run fake.in to detect compilation options
    fake_test = TestCase(input_file='fake.in', name='fake')
    molgw_ok = clean_run(fake_test, tmpfolder, mpirun, nprocs, args.debug)
    if not molgw_ok:
        print('MOLGW executable is not functional\nDump last output:')
        print((tmpfolder / 'fake.out').read_text(encoding='utf-8', errors='replace'))
        sys.exit('MOLGW executable is not functional')

    fake_content = read_fake_out(tmpfolder / 'fake.out')

    have_openmp        = 'Running with OPENMP' in fake_content
    have_libxc         = 'Running with LIBXC' in fake_content
    have_mpi           = 'Running with MPI' in fake_content
    have_scalapack     = 'Running with SCALAPACK' in fake_content
    have_onebody       = 'Running with external LIBINT or LIBCINT calculation of the one-body operators' in fake_content
    have_gradients     = ('Running with external LIBINT calculation of the gradients of the one-body integrals' in fake_content
                          or 'Code compiled with LIBCINT' in fake_content)
    have_libint_forces = False  # FIXME: force calculation is broken as of today
    is_libcint         = 'Code compiled with LIBCINT support' in fake_content
    have_hdf5          = 'Running with HDF5' in fake_content
    have_fftw3         = 'Running with FFTW3' in fake_content

    print('MOLGW compilation details:')
    print(f'              MOLGW version: {version}')
    print(f'                     OPENMP: {have_openmp}')
    print(f'                        MPI: {have_mpi}')
    print(f'                  SCALAPACK: {have_scalapack}')
    print(f'                      LIBXC: {have_libxc}')
    print(f'                       HDF5: {have_hdf5}')
    print(f'                      FFTW3: {have_fftw3}')
    print(f'                  integrals: {"LIBCINT" if is_libcint else "LIBINT"}')
    print(f'           1-body integrals: {have_onebody}')
    print(f'        gradients integrals: {have_gradients}')
    print(f' LIBINT gradients integrals: {have_libint_forces}')
    print()

    if have_openmp:
        os.environ['OMP_STACKSIZE'] = '32M'

    # Parse testsuite file
    all_tests = parse_testsuite(Path('inputs/testsuite'))

    # --list
    if args.list:
        print('=== List of input files in the Suite ===')
        for i, test in enumerate(all_tests, start=1):
            print(f'{i:04d}: {test.input_file}')
        print('========================================')
        sys.exit(0)

    # Apply filters
    tests = apply_param_filter(all_tests, args.input_parameter, all_tests)
    tests = apply_selection(tests, args.input, args.exclude)

    print(f'Input files to be executed: {len(tests)}')

    # Run tests
    stats = RunStats()

    with (tmpfolder / 'diff').open('w') as fdiff, \
         (tmpfolder / 'failed_tests').open('w') as ffailed:

        fdiff.write(
            '#   test index       property tested'
            '                     calculated'
            '                 reference'
            '                 difference        test status \n'
        )

        for test in tests:
            # Skip checks
            skip_msg = None
            if test.need_libcint and not is_libcint:
                skip_msg = 'this compilation of MOLGW does not have LIBCINT'
            elif test.need_fftw3 and not have_fftw3:
                skip_msg = 'this compilation of MOLGW does not have FFTW3'
            elif test.need_scalapack and not have_scalapack:
                skip_msg = 'this compilation of MOLGW does not have SCALAPACK'
            elif test.need_gradients and not have_gradients:
                skip_msg = 'this compilation of MOLGW does not have the integral gradients'
            elif test.need_forces and not have_libint_forces:
                skip_msg = 'this compilation of MOLGW does not have the force integrals'
            elif not test.parallel and nprocs > 1:
                skip_msg = 'this test is only serial'

            if skip_msg:
                stats.files_skipped += 1
                print(f'\nSkipping test file: {test.input_file}')
                print(f'  because {skip_msg}')
                stats.skipping_reasons.append(skip_msg)
                continue

            print(f'\nRunning test file: {test.input_file}')
            print(test.name)
            fdiff.write(f'# {test.input_file}\n')

            clean_run(test, tmpfolder, mpirun, nprocs, args.debug)
            check_output(test, tmpfolder, stats, fdiff, ffailed, args.verbose, args.debug)

    # Summary
    n_run = len(tests) - stats.files_skipped
    color = '\033[92m' if stats.success == stats.tested else '\033[91m'
    reset = '\033[0m'
    bold  = '\033[1m'

    print('\n\n===============================')
    print('      Test Suite Summary \n')
    print(f'      Test files tested:   {n_run:4d} / {len(tests):4d}\n')
    print(f'     Test files success:   {color}{bold}{stats.files_success:4d} / {n_run:4d}{reset}  ')
    print(f'       Successful tests:   {color}{bold}{stats.success:4d} / {stats.tested:4d}{reset}\n')
    print(f'       Elapsed time (s):   {time.time() - start_time:.2f}')
    print('===============================\n')

    if stats.files_skipped > 0:
        print(' Some tests have been skipped for the following reasons:')
        for reason in set(stats.skipping_reasons):
            count = stats.skipping_reasons.count(reason)
            print(f'   * {reason:<80}  ({count:=4d} tests)')
        print('===============================\n')

    if not args.keep:
        shutil.rmtree(tmpfolder)

    sys.exit(abs(stats.success - stats.tested))


if __name__ == '__main__':
    main()
