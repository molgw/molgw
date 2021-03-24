#!/usr/bin/python3
##################################################
#
# This file is part of MOLGW
# Author: Fabien Bruneval
#
# This python module provides useful functions to read molgw.yaml files
# 
#
##################################################

import os
from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper



def get_chemical_formula(calc):
    element_list  = [ 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'Si', 'As', 'Se', 'I', 'Li', 'Br', 'H', 'Mg', 'Rb', 'Cu', 'Na', 'Al', 'Ag', 'Kr', 'Be', 'Ge', 'K', 'He', 'Xe', 'Ne', 'He', 'Ti', 'Ga', 'Ar']

    numbers = [0 for x in element_list]
    
    for atom in calc["physical system"]["atom list"]:
        i = element_list.index(atom[0])
        numbers[i] += 1
    formula = ''
    for e,n in zip(element_list,numbers):
        if n > 1:
            formula += e + str(n)
        elif n > 0:
            formula += e
    return formula


def get_homo_energy(approx,calc):
    homo = int(calc["physical system"]["electrons"] * 0.50)
    key = approx + " energy" 
    try:
        ehomo = -9999.9
        for state in calc[key]["spin channel 1"].keys():
            if int(state) <= homo:
                ehomo = max(ehomo,float(calc[key]["spin channel 1"][state]))
        return ehomo
    except:
        return None

def parse_yaml_files(directory):
    ########################################################################
    # List all the yaml files in the directory
    #
    yaml_files = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if '.yaml' in filename:
                yaml_files.append(os.path.join(root, filename))
    print('{} molgw.yaml files identified in directory '.format(len(yaml_files)) + directory)

    ########################################################################
    # Read all the yaml files -> dictionary
    #
    calc = []
    for yaml_file in yaml_files:
        with open(yaml_file, 'r') as stream:
            try:
                calc.append(load(stream,Loader=Loader))
            except:
                print(yaml_file + ' is corrupted')
                pass
    return calc



