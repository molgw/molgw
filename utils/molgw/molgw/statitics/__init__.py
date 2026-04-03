#!/usr/bin/python3
##################################################
#
# This file is part of MOLGW
# Author: Fabien Bruneval
#
# This python submodule provides useful functions to analyze data
#
#
##################################################

import numpy as np
import json
import matplotlib.pyplot as plt
from molgw import __version__, Molecule


########################################################################
def diff(data1, data2):
    """
        Returns a dictionary containing data2 - data1
    """
    errors = dict()
    shared_keys = set(data1.keys()) & set(data2.keys())
    for cas in shared_keys:
        errors[cas] = data2[cas] - data1[cas]
    return errors

########################################################################
def mae_mse_max(data1, data2):
    """
        Returns MAE, MSE, Max errors with data1 being the reference
    """
    errors = list( diff(data1, data2).values() )

    ndata = len(errors)
    mse = np.sum(errors) / float(ndata)
    mae = np.sum(np.abs(errors)) / float(ndata)
    mxe = np.max(np.abs(errors))
    return mae, mse, mxe

########################################################################
def correlation_plot(data1, data2, labels=False):
    fig = plt.figure(figsize=(6, 6))
    xxx = []
    yyy = []
    cass = []
    shared_keys = set(data1.keys()) & set(data2.keys())
    for cas in shared_keys:
        xxx.append(data1[cas])
        yyy.append(data2[cas])
        cass.append(cas)
    xmin = min(xxx + yyy) - 0.1
    xmax = max(xxx + yyy) + 0.1
    plt.plot([xmin, xmax], [xmin, xmax], '-', color='black')
    plt.scatter(xxx, yyy, marker='o')
    if labels:
        delta = abs(xmax - xmin) * 0.02
        for x, y, cas in zip(xxx, yyy, cass):
            plt.text(x + delta, y - delta, cas)
    plt.xlabel("HOMO (eV)")
    plt.ylabel("HOMO (eV)")
    return fig

