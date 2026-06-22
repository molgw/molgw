#!/usr/bin/env python3
"""
CP2K -> MOLGW basis set converter
Handles multiple CP2K blocks in one basis definition.
"""

import sys
from collections import OrderedDict


def parse_single_block(lines, start):
    """
    Parse one CP2K basis block starting at index 'start' in lines.
    Returns:
       next_index, exponents, coeff_columns, shell_counts
    """

    parts = [int(x) for x in lines[start].split()]
    n, lmin, lmax, nexp = parts[:4]
    nshells = parts[4:]

    shell_counts = OrderedDict()
    for ell, count in zip(range(lmin, lmax + 1), nshells):
        if count > 0:
            shell_counts[ell] = count

    exponents = []
    coeff_matrix = []

    idx = start + 1
    for _ in range(nexp):
        nums = [float(x) for x in lines[idx].split()]
        exponents.append(nums[0])
        coeff_matrix.append(nums[1:])
        idx += 1

    # Transpose: columns = contracted functions
    coeff_columns = list(map(list, zip(*coeff_matrix)))

    return idx, exponents, coeff_columns, shell_counts


def parse_cp2k_file(lines):
    """
    Parse an entire CP2K basis entry, possibly containing several blocks.
    Returns combined list of (exponents, coefficients, l)
    """

    # Strip comments and blanks
    data = [
        l.strip() for l in lines
        if l.strip() and not l.strip().startswith("#")
    ]

    # First line: element and basis name -> skip
    data = data[1:]

    # Number of blocks
    nblocks = int(data[0])
    data = data[1:]

    idx = 0

    all_functions = []

    for _ in range(nblocks):
        idx, exponents, coeff_cols, shell_counts = parse_single_block(data, idx)

        col = 0
        for l, count in shell_counts.items():
            for _ in range(count):
                all_functions.append((l, exponents, coeff_cols[col]))
                col += 1

    return all_functions


def write_molgw(functions, out=sys.stdout):
    """
    Write MOLGW format from list of (l, exponents, coeffs)
    """

    print(len(functions), file=out)

    for l, exps, coeffs in functions:
        print(f"{len(exps)} {l}", file=out)
        for e, c in zip(exps, coeffs):
            print(f"  {e:18.12f} {c:18.12f}", file=out)


def convert_file(filename):
    with open(filename) as f:
        functions = parse_cp2k_file(f.readlines())

    write_molgw(functions)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: basisset_cp2k2molgw.py <cp2k_basis_file>")
        sys.exit(1)

    convert_file(sys.argv[1])

