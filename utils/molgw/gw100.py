#!/usr/bin/python3
##################################################
#
# This file is part of MOLGW
# Author: Fabien Bruneval
#
# This python submodule provides useful functions to prepare/analyze gw100 runs
#
#
##################################################

import numpy as np
import json
import matplotlib.pyplot as plt
from . import __version__

chemical_formulas = {'7446-09-5': 'SO2', '71-30-7': 'C4H5N3O', '66-22-8': 'C4H4N2O2', '74-84-0': 'C2H6', '1603-84-5': 'OCSe', '74-90-8': 'NCH', '7784-18-1': 'AlF3', '7784-23-8': 'AlI3', '7727-37-9': 'N2', '75-73-0': 'CF4', '507-25-5': 'CI4', '7722-84-1': 'HOOH', '7783-63-3': 'TiF4', '56-23-5': 'CCl4', '558-13-4': 'CBr4', '13768-60-0': 'BF', '302-01-2': 'N2H4', '7553-56-2': 'I2', '12184-80-4': 'C4', '74-98-6': 'C3H8', '19287-45-7': 'B2H6', '7726-95-6': 'Br2', '7786-30-3': 'MgCl2', '7782-50-5': 'Cl2', '13283-31-3': 'BH3', '1304-56-9': 'BeO', '7664-39-3': 'HF', '7784-42-1': 'AsH3', '74-86-2': 'C2H2', '60-29-7': 'C4H10O', '65-71-4': 'C5H6N2O2', '73-40-5': 'C5H5N5O', '7439-90-9': 'Kr', '463-58-1': 'OCS', '67-56-1': 'CH4O', '64-17-5': 'C2H6O', '64-18-6': 'CH2O2', '1333-74-0': 'H2', '12190-70-4': 'Cu2', '7440-37-1': 'Ar', '75-15-0': 'CS2', '7664-41-7': 'NH3', '7782-79-8': 'HN3', '7758-02-3': 'KBr', '50-00-0': 'H2CO', '7440-59-7': 'He', '7647-01-0': 'HCl', '17108-85-9': 'GaCl', '7693-26-7': 'KH', '7803-51-2': 'PH3', '73-24-5': 'C5H5N5', '630-08-0': 'CO', '1590-87-0': 'Si2H6', '7732-18-5': 'H2O', '75-07-0': 'C2H4O', '124-38-9': 'CO2', '75-02-5': 'C2H3F', '7783-06-4': 'SH2', '75-01-4': 'C2H3Cl', '7783-60-0': 'SF4', '7647-14-5': 'NaCl', '74-85-1': 'C2H4', '12187-06-3': 'Ag2', '75-19-4': 'C3H6', '7789-24-4': 'LiF', '7440-01-9': 'Ne', '14868-53-2': 'Si5H12', '74-82-8': 'CH4', '25681-79-2': 'Na2', '7782-65-2': 'GeH4', '7803-62-5': 'SiH4', '25681-80-5': 'K2', '25681-81-6': 'Rb2', '39297-86-4': 'Na4', '7580-67-8': 'LiH', '542-92-7': 'C5H6', '7783-40-6': 'MgF2', '39297-88-6': 'Na6', '629-20-9': 'C8H8', '62-53-3': 'C6H5NH2', '12185-09-0': 'P2', '100-41-4': 'C8H10', '71-43-2': 'C6H6', '23878-46-8': 'As2', '108-88-3': 'C7H8', '110-86-1': 'C5H5N', '14452-59-6': 'Li2', '17739-47-8': 'PN', '392-56-3': 'C6F6', '7782-41-4': 'F2', '10028-15-6': 'O3', '106-97-8': 'C4H10', '544-92-3': 'NCCu', '10043-11-5': 'BN', '1309-48-4': 'MgO', '57-13-6': 'CH4N2O', '593-66-8': 'C2H3I', '108-95-2': 'C6H5OH', '593-60-2': 'C2H3Br', '7440-63-3': 'Xe'}

molecules = list(chemical_formulas.keys())

default_input_parameters = { 
      'scf': 'pbeh',
      'alpha_hybrid': 0.75,
      'postscf': 'g0w0',
      'basis': 'Def2-TZVPP',
      'ecp_basis': 'Def2-TZVPP',
      'auxil_basis': 'pauto',
      'ecp_auxil_basis': 'pauto',
      'ecp_elements': 'Xe Ag Rb I',
      'ecp_type': 'Def2-ECP',
      'selfenergy_state_range': 2,
      'frozencore': 'yes'
              }


structures = {
    "7440-37-1": "1\nArgon; atom; s\nAr 0.0 0.0 0.0\n",
    "7440-01-9": "1\nNeon; atom; s\nNe 0.0 0.0 0.0\n",
    "7440-63-3": "1\nXenon; atom; s\nXe 0.0 0.0 0.0\n",
    "7439-90-9": "1\nKrypton; atom; s\nKr 0.0 0.0 0.0\n",
    "7440-59-7": "1\nHelium; atom; s\nHe 0.0 0.0 0.0\n",
    "12187-06-3": "2\nSilver dimer; experimental structure form simard01; s\nAg 0 0 0\nAg 0 0 2.5335\n",
    "12190-70-4": "2\nCopper dimer; experimental structure form HCP92; s\nCu 0.0 0.0 0.0\nCu 0.0 0.0 2.2197\n",
    "7553-56-2": "2\nIodine; experimental structure from HCP92; s\nI 0.0000 0.0000 0.0000\nI 0.0000 0.0000 2.6663\n",
    "7726-95-6": "2\nBromine; experimental structure from HCP92; s\nBr 0.0000 0.0000 0.0000\nBr 0.0000 0.0000 2.2811\n",
    "7782-41-4": "2\nFluorine; experimental structure from HCP92; s\nF 0.0000 0.0000 0.0000\nF 0.0000 0.0000 1.4119\n",
    "7727-37-9": "2\nNitrogen; experimental structure from HCP92; s\nN 0.0000 0.0000 0.0000\nN 0.0000 0.0000 1.0977\n",
    "1333-74-0": "2\nHydrogen; experimental structure from HCP92; s\nH 0.0000 0.0000 0.0000\nH 0.0000 0.0000 0.74144\n",
    "7782-50-5": "2\nChlorine; experimental structure from HCP92; s\nCl 0.0000 0.0000 0.0000\nCl 0.0000 0.0000 1.9878\n",
    "13768-60-0": "2\nFluoroborane; experimental structure from HCP92; s\nB 0.0000 0.0000 0.0000\nF 0.0000 0.0000 1.2626\n",
    "25681-80-5": "2\nDipotassium; experimental structure form HCP92; s\nK 0.0000 0.0000 0.0000\nK 0.0000 0.0000 3.9051\n",
    "10043-11-5": "2\nBoron nitride; experimental structure from HCP92; s\nB 0.0000 0.0000 0.0000\nN 0.0000 0.0000 1.281\n",
    "25681-81-6": "2\nDirubidium; experimental structure from JANAF; s\nRb 0.0000 0.0000 0.0000\nRb 0.0000 0.0000 4.12256\n",
    "14452-59-6": "2\nLithium dimer; experimental structure from HCP92; s\nLi 0.0000 0.0000 0.0000\nLi 0.0000 0.0000 2.6729\n",
    "7693-26-7": "2\nPotassium hydride; experimental structure from HCP92; s\nK 0.0000 0.0000 0.0000\nH 0.0000 0.0000 2.244\n",
    "630-08-0": "2\nCarbon monoxide; experimental structure from HCP92; s\nC 0.0000 0.0000 0.0000\nO 0.0000 0.0000 1.283\n",
    "25681-79-2": "2\nSodium dimer; experimental structure from HCP92; s\nNa 0.0000 0.0000 0.0000\nNa 0.0000 0.0000 3.0789\n",
    "23878-46-8": "2\nArsenic dimer; experimental structure from HCP92; s\nAs 0.0000 0.0000 0.0000\nAs 0.0000 0.0000 2.1026\n",
    "7758-02-3": "2\nPotassium bromide; experimental structure from HCP92; s\nBr 0.0000 0.0000 0.0000\nK  0.0000 0.0000 2.8208\n",
    "7789-24-4": "2\nLithium fluoride; experimental structure from HCP92; s\nLi 0.0000 0.0000 0.0000\nF  0.0000 0.0000 1.5639\n",
    "7580-67-8": "2\nLithium hydride; experimental structure from HCP92; s\nLi 0.0000 0.0000 0.0000\nH  0.0000 0.0000 1.5949\n",
    "7664-39-3": "2\nHydrogen fluoride; experimental structure from HCP92; s\nH 0.0000 0.0000 0.0000\nF 0.0000 0.0000 0.9169\n",
    "1309-48-4": "2\nMagnesium monoxide; experimental structure from HCP92; s\nMg 0.0000 0.0000 0.0000\nO  0.0000 0.0000 1.749\n",
    "7647-14-5": "2\nSodium chloride; experimental structure from HCP92; s\nNa 0.0000 0.0000 0.0000\nCl 0.0000 0.0000 2.3609\n",
    "1304-56-9": "2\nBeryllium monoxide; experimental structure from HCP92; s\nBe 0.0000 0.0000 0.0000\nO  0.0000 0.0000 1.3308\n",
    "12185-09-0": "2\nPhosphorus dimer; experimental structure from HCP92; s\nP 0.0000 0.0000 0.0000\nP 0.0000 0.0000 1.8931 \n",
    "10028-15-6": "3\nOzon; experimental structure from HCP92; s\nO  0.0000 0.0000 0.0000\nO  1.0869 0.0000 0.6600\nO -1.0869 0.0000 0.6600",
    "7647-01-0": "2\nHydrogen chloride; experimental structure from HCP92; s\nH  0.0000 0.0000 0.0000\nCl 0.0000 0.0000 1.2746\n",
    "7732-18-5": "3\nWater; experimental structure from HCP92; s\nO  0.0000 0.0000 0.0000\nH  0.7571 0.0000 0.5861\nH -0.7571 0.0000 0.5861",
    "17739-47-8": "2\nPhosphorus mononitride; experimental structure from HCP92; s\nP 0.0000 0.0000 0.0000\nN 0.0000 0.0000 1.49087\n",
    "17108-85-9": "2\nGallium monochloride; experimental structure from HCP92; s\nGa 0.0000 0.0000 0.0000\nCl 0.0000 0.0000 2.2017\n",
    "74-90-8": "3\nHydrogen cyanide; Experimental structure form HCP92; s\nC 0.0000 0.0000 0.0000\nH 0.0000 0.0000 1.0655\nN 0.0000 0.0000 -1.1532\n",
    "544-92-3": "3\nCopper cyanide; Experimental structure from HCP92; s\nC 0.0000 0.0000 0.0000\nN 0.0000 0.0000 1.158\nCu 0.0000 0.0000 -1.832\n",
    "124-38-9": "3\nCarbon dioxide; experimental structure from HCP92; s\nO 0.0000 0.0000  1.16\nC 0.0000 0.0000  0.0000\nO 0.0000 0.0000 -1.16\n",
    "463-58-1": "3\nCarbon oxysulfide; Experimental structure from HCP92; s\nO 0.0000 0.0000 1.1578\nC 0.0000 0.0000 0.0000\nS 0.0000 0.0000 -1.5601\n",
    "1603-84-5": "3\nCarbon oxyselenide; Experimental structure from HCP92; s\nO 0.0000 0.0000 1.159\nC 0.0000 0.0000 0.0000\nSe 0.0000 0.0000 -1.709\n",
    "7446-09-5": "3\nSulfer dioxide; experimental structure from HCP92; m\nS  0.0000 0.0000 0.0000\nO  1.2349 0.0000 0.7226\nO -1.2349 0.0000 0.7226",
    "75-15-0": "3\nCarbon disulfide; experimental structure from HCP95; s\nC  0.0000 0.0000  0.0000\nS  0.0000 0.0000  1.5526\nS  0.0000 0.0000 -1.5526\n",
    "7783-06-4": "3\nHydrogen sulfide; experimental structure from HCP92; s\nS  0.0000 0.0000 0.0000\nH  0.9617 0.0000 0.9268\nH -0.9617 0.0000 0.9268\n",
    "7783-40-6": "3\nMagnesium fluoride; experimental structure from HCP92; s\nF  0.0000 0.0000  1.771\nMg 0.0000 0.0000  0.0000\nF  0.0000 0.0000 -1.771\n",
    "7786-30-3": "3\nMagnesium chloride; experimental structure from HCP92; s\nMg  0.0000  0.0000  0.0000\nCl  0.0000  0.0000  2.179 \nCl  0.0000  0.0000 -2.179\n",
    "74-86-2": "4\nAcetylene; experimental structure from HCP92; s\nC 0.0000 0.0000  0.6015\nC 0.0000 0.0000 -0.6015\nH 0.0000 0.0000  1.6615\nH 0.0000 0.0000 -1.6615",
    "7664-41-7": "4\nAmonia; experimental structure from HCP92; s\nN  0.0000  0.0000  0.0000\nH  0.0000 -0.9377 -0.3816\nH  0.8121  0.4689 -0.3816\nH -0.8121  0.4689 -0.3816",
    "7803-51-2": "4\nPhosphine; experimental structure from HCP92; s\nP  0.0000  0.0000  0.0000\nH  0.0000 -1.1932 -0.7717\nH  1.0333  0.5966 -0.7717\nH -1.0333  0.5966 -0.7717",
    "13283-31-3": "4\nBorane; experimental structure from HCP92; s\nB  0.0000  0.0000  0.0000\nH   0.0000  0.0000  1.19 \nH   0.0000  1.0306 -0.595\nH   0.0000 -1.0306 -0.595\n",
    "50-00-0": "4\nFormaldehyde; experimental structure from HCP92; s\nC  0.0000 0.0000  0.0000\nO  0.0000 0.0000  1.208\nH  0.9490 0.0000 -0.5873\nH -0.9490 0.0000 -0.5873",
    "7782-79-8": "4\nHydrogen azide; experimental structure from HCP92; s\nH -0.9585 0.0000 -0.3338\nN  0.0000 0.0000  0.0000\nN  0.0000 0.0000  1.2450\nN  0.1617 0.0000  2.3674\n",
    "7722-84-1": "4\nHydrogen peroxide; experimental structure from HCP92; s\nO  0.0000  0.7375 -0.0528\nO  0.0000 -0.7375 -0.0528\nH  0.8190  0.8170  0.4220\nH -0.8190 -0.8170  0.4220\n",
    "7784-42-1": "4\nArsine; experimental structure from HCP92; s\nAs  0.0000  0.0000  0.0000\nH   0.0000  1.2561  0.8398 \nH   1.0878 -0.6281  0.8398\nH  -1.0878 -0.6281  0.8398",
    "7784-23-8": "4\nAluminum triiodide; experimental structure from HCP92; s\nAl  0.0000  0.0000  0.0000\nI   0.0000  0.0000  2.461 \nI   0.0000  2.1313 -1.2305\nI   0.0000 -2.1313 -1.2305\n",
    "7784-18-1": "4\nAluminum trifluoride; experimental structure from HCP92; s\nAl  0.0000  0.0000  0.0000\nF   0.0000  0.0000  1.633\nF   0.0000  1.4142 -0.8165\nF   0.0000 -1.4142 -0.8165\n",
    "12184-80-4": "4\nTetracarbon; experimental structure from weltner89 (calculated); m\nC  1.2247  0.0000 0.0000\nC -1.2247  0.0000 0.0000\nC  0.0000 -0.7286 0.0000 \nC  0.0000  0.7286 0.0000",
    "64-18-6": "5\nFormic Acid; experimental structure from HCP92; m\nO  0.9858 0.0000  2.0307 \nH -1.0241 0.0000  1.7361\nC  0.0000 0.0000  1.3430\nO  0.0000 0.0000  0.0000\nH  0.9329 0.0000 -0.2728",
    "74-82-8": "5\nMethane; experimental structure from HCP92; m\nC  0.0000  0.0000  0.0000\nH  0.6276 -0.6275  0.6276\nH -0.6276  0.6276  0.6276     \nH -0.6276 -0.6276 -0.6276 \nH  0.6276  0.6276 -0.6276",
    "7803-62-5": "5\nSilane; experimental structure from HCP92; m\nSi  0.0000  0.0000  0.0000\nH   0.8544 -0.8544  0.8544  \nH  -0.8544  0.8544  0.8544    \nH  -0.8544 -0.8544 -0.8544 \nH   0.8544  0.8544 -0.8544 ",
    "7783-63-3": "5\nTitanium fluoride; experimental structure from HCP92; m\nTi  0.0000  0.0000  0.0000\nF   1.0127 -1.0127  1.0127  \nF  -1.0127  1.0127  1.0127    \nF  -1.0127 -1.0127 -1.0127 \nF   1.0127  1.0127 -1.0127 ",
    "7782-65-2": "5\nGermane; experimental structure from HCP92; m\nGe  0.0000  0.0000  0.0000\nH   0.8805 -0.8805  0.8805  \nH  -0.8805  0.8805  0.8805    \nH  -0.8805 -0.8805 -0.8805 \nH   0.8805  0.8805 -0.8805 ",
    "7783-60-0": "5\nSulfer tetrafluoride; experimental structure from tolles61; m\nS  0.0000  0.0000  0.3726\nF  0.0000  1.6430  0.2731\nF  0.0000 -1.6430  0.2731\nF  1.1969  0.0000 -0.6044\nF -1.1969  0.0000 -0.6044 ",
    "74-85-1": "6\nEthylene; experimental structure from HCP92; s\nC  0.0000 0.0000  0.0000\nC  0.0000 0.0000  1.3290\nH  0.9235 0.0000 -0.5637\nH -0.9235 0.0000 -0.5637\nH  0.9235 0.0000  1.8927\nH -0.9235 0.0000  1.8927",
    "75-73-0": "5\nCarbon tetrafluoride; experimental structure from HCP92; m\nC  0.0000  0.0000  0.0000\nF  0.7638 -0.7638  0.7638  \nF -0.7638  0.7638  0.7638    \nF -0.7638 -0.7638 -0.7638 \nF  0.7638  0.7638 -0.7638 ",
    "507-25-5": "5\nCarbon tetraiodide; experimental structure from HCP92; m\nC  0.0000  0.0000  0.0000\nI  1.2411 -1.2411  1.2411  \nI -1.2411  1.2411  1.2411    \nI -1.2411 -1.2411 -1.2411 \nI  1.2411  1.2411 -1.2411 ",
    "302-01-2": "6\nHydrazene; experimental structure from HCP92; m\nN  0.0000  0.7230 -0.1123\nN  0.0000 -0.7230 -0.1123\nH -0.4470  1.0031  0.7562\nH  0.4470 -1.0031  0.7562\nH  0.9663  1.0031  0.0301\nH -0.9663 -1.0031  0.0301",
    "56-23-5": "5\nCarbon tetrachloride; experimental structure from HCP92; m\nC   0.0000  0.0000  0.0000\nCl  1.0202 -1.0202  1.0202  \nCl -1.0202  1.0202  1.0202    \nCl -1.0202 -1.0202 -1.0202 \nCl  1.0202  1.0202 -1.0202 ",
    "558-13-4": "5\nCarbon tetrabromide; experimental structure from HCP92; m\nC   0.0000  0.0000  0.0000\nBr  1.1172 -1.1172  1.1172  \nBr -1.1172  1.1172  1.1172    \nBr -1.1172 -1.1172 -1.1172 \nBr  1.1172  1.1172 -1.1172 ",
    "39297-86-4": "4\nSodium tetramer; TM QZVP pbe optimized; s\nNa    0.0002445   -0.0998053    1.5471126 \nNa   -0.0002444    3.1776586    0.0486374 \nNa    0.0002444    0.0997722   -1.5472150 \nNa   -0.0002444   -3.1776254   -0.0485350 \n",
    "67-56-1": "6\nMethanol; experimental structure from HCP92; m\nC -0.722791 -0.007039 0.000000\nO 0.701687 0.011691 0.000000\nH -1.022488 -1.059189 0.000000\nH -1.162308 0.455180 -0.888335\nH -1.148266 0.468270 0.888264\nH 0.990233 0.911667 0.000000\n",
    "593-66-8": "6\nVynil iodide; experimental structure from HCP92; m\nC 0.000000 0.000000 0.000000\nC 0.000000 0.000000 1.328000\nH -0.903971 0.000000 -0.595154\nH -0.899347 0.000000 1.924168\nH 0.912883 0.000000 -0.577101\nI 1.747544 0.000000 2.461568\n",
    "75-02-5": "6\nVynil fluoride; experimental structure from HCP92; m\nC 0.000000 0.000000 0.000000\nC 0.000000 0.000000 1.321000\nH -0.942589 0.000000 -0.521841\nH -0.874292 0.000000 1.955045\nH 0.922424 0.000000 -0.561725\nF 1.142469 0.000000 2.026603\n",
    "75-01-4": "6\nVynil chloride; experimental structure from HCP92; m\nC -0.554265 -0.445361 0.111076\nC 0.372254 0.438035 -0.234540\nH -1.322093 -0.210763 0.831940\nH -0.543697 -1.425246 -0.340536\nH 1.153254 0.241370 -0.951543\nCl 0.440430 2.028766 0.431795\n",
    "593-60-2_old": "6\nVynil bromide; experimental structure from HCP92; m\nC 0.000000 0.000000 0.000000\nC 0.000000 0.000000 1.325600\nH -0.895976 0.000000 -0.602298\nH -0.894897 0.000000 1.927173\nH 0.908386 0.000000 -0.581003\nBr 1.357668 0.000000 2.194533\n",
    "593-60-2": "6\nVynil bromide; experimental structure from HCP92 revised structure; m\nBr 6.004480 4.499660 6.133020\nC 4.420800 4.500000 3.789300\nC 4.420800 4.500000 5.114900\nH 3.477470 4.501520 3.258560\nH 3.525900 4.500000 5.716470\nH 5.329190 4.500000 3.208300\n",
    "75-07-0": "7\nAcetaldehyde; structure from HCP 92;m\nC  0.000000  0.000000  0.000000\nC  0.000000  0.000000  1.515000\nO  1.001953  0.000000  2.193373\nH -1.019805  0.000000  1.997060\nH -0.905700 -0.522900 -0.363000\nH  0.000000  1.045800 -0.363000\nH  0.905700 -0.522900 -0.363000\n",
    "57-13-6": "8\nUrea; experimental structure from godfrey97; m\nO  0.0000  1.3049  0.0000\nC  0.0000  0.0838  0.0000\nN  1.1603 -0.6595  0.0000\nN -1.1603 -0.6595 -0.0000\nH  1.1383 -1.5964  0.3424\nH  1.9922 -0.0940  0.1760\nH -1.1383 -1.5964 -0.3424\nH -1.9922 -0.0940 -0.1760",
    "64-17-5": "9\nEthanol; Experimental structure from cou98; m\nC  1.1879 -0.3829  0.0000\nC  0.0000  0.5526  0.0000\nO -1.1867 -0.2472  0.0000\nH -1.9237  0.3850  0.0000\nH  2.0985  0.2306  0.0000\nH  1.1184 -1.0093  0.8869\nH  1.1184 -1.0093 -0.8869\nH -0.0227  1.1812  0.8852\nH -0.0227  1.1812 -0.8852\n",
    "19287-45-7": "8\nDiborane6; experimental structure from HCP92; m\nB  0.0000  0.0000  0.8870\nB  0.0000  0.0000 -0.8870\nH  0.9960  0.0000  0.0000  \nH -0.9960  0.0000  0.0000  \nH  0.0000  1.0408  1.4639\nH  0.0000 -1.0408  1.4639     \nH  0.0000  1.0408 -1.4639     \nH  0.0000 -1.0408 -1.4639  \n",
    "1590-87-0": "8\nDisilane; structure form HCP92; m\nSi 0.000000  0.000000 -1.165500\nSi 0.000000  0.000000  1.165500\nH  1.399330  0.000000  1.683128\nH -1.399330  0.000000 -1.683128\nH  0.699500  1.211600 -1.683128\nH  0.699500 -1.211600 -1.683128\nH -0.699500 -1.211600  1.683128\nH -0.699500  1.211600  1.683128\n",
    "74-84-0": "8\nEthane; experimental structure from HCP92; m\nC 0.000100 0.765700 0.000000\nC -0.000100 -0.769300 0.000000\nH -0.134892 1.160801 -1.011192\nH -0.808178 1.160889 0.622373\nH 0.943359 1.160660 0.388723\nH 0.808348 -1.164489 -0.622151\nH 0.134598 -1.164401 1.011231\nH -0.943239 -1.164260 -0.389014\n",
    "39297-88-6": "6\nSodium hexamer; TM def2-QZVP PBE optimized; m\nNa   -2.4732949   -1.7969539   -0.2367313 \nNa   -2.4732949    1.7969539   -0.2367313 \nNa    0.9447146   -2.9075325   -0.2367313 \nNa    0.9447146    2.9075325   -0.2367313 \nNa    3.0571606    0.0000000   -0.2367313 \nNa    0.0000000    0.0000000    1.1836565 \n",
    "75-19-4": "9\nCyclopropane; experimental structure from HCP92; m\nC 0.036473 0.859901 -0.182257\nC -0.275674 -0.564282 -0.600227\nC 0.268426 -0.283892 0.786852\nH -0.412409 -0.380348 1.623541\nH 1.285827 -0.580215 1.010404\nH -0.796804 1.523474 0.013213\nH 0.896783 1.350428 -0.620593\nH -1.331415 -0.803627 -0.632093\nH 0.315621 -1.079884 -1.346828\n",
    "74-98-6": "11\nPropane; experimental structure form hellwege76; m\nC  0.0000  0.5863 -0.0000\nC -1.2681 -0.2626  0.0000\nC  1.2681 -0.2626 -0.0000\nH  0.0000  1.2449  0.8760\nH -0.0003  1.2453 -0.8758\nH -2.1576  0.3742  0.0000\nH  2.1576  0.3743  0.0000\nH -1.3271 -0.9014  0.8800\nH -1.3271 -0.9014 -0.8800\nH  1.3271 -0.9014 -0.8800\nH  1.3272 -0.9014  0.8800",
    "71-43-2": "12\nBenzene; experimental structure from HCP92; l\nC  0.0000  1.3990 0.0000\nC  1.2115  0.6995 0.0000\nC  1.2115 -0.6995 0.0000\nC  0.0000 -1.3990 0.0000\nC -1.2115 -0.6995 0.0000\nC -1.2115  0.6995 0.0000\nH  0.0000  2.5000 0.0000\nH  2.1651  1.2500 0.0000\nH  2.1651 -1.2500 0.0000\nH  0.0000 -2.5000 0.0000\nH -2.1651 -1.2500 0.0000\nH -2.1651  1.2500 0.0000\n",
    "110-86-1": "11\nPyridine; expoerimental from HCP92; l\nN 0.000000 0.000000 0.000000\nC -0.476428 -1.252444 0.000000\nC -0.903103 0.989952 0.000000\nC -2.282876 0.784403 0.000000\nC -1.835282 -1.567988 0.000000\nC -2.760265 -0.525306 0.000000\nH -0.532213 2.008528 0.000000\nH 0.266630 -2.041697 0.000000\nH -2.958369 1.628364 0.000000\nH -2.153556 -2.601071 0.000000\nH -3.818275 -0.726658 0.000000\n",
    "542-92-7": "11\nCyclopentadiene; experimental structure from HCP92; l\nC 0.735000 0.000000 0.000000\nC -0.735000 0.000000 0.000000\nC 1.180760 0.000000 1.265805\nC -1.180760 0.000000 1.265805\nC -0.003091 0.000000 2.209296\nH 2.228506 0.000000 1.566354\nH 1.364644 0.000000 -0.889746\nH -1.364960 0.000000 -0.889523\nH -2.228838 0.000000 1.565193\nH -0.001086 0.885447 2.844969\nH -0.005702 -0.882815 2.848618\n",
    "392-56-3": "12\nHexafluorobenzene; esperimental structure form hellwege76; l\nC -1.092692 -0.875775 0.000000\nC 0.212112 -1.384183 0.000000\nC 1.304703 -0.508408 0.000000\nC 1.092692 0.875775 0.000000\nC -0.212112 1.384083 0.000000\nC -1.304804 0.508308 0.000000\nF -2.126549 -1.704487 0.000000\nF 0.412876 -2.693885 0.000000\nF 2.539354 -0.989305 0.000000\nF 2.126549 1.704487 0.000000\nF -0.412858 2.693787 0.000000\nF -2.539356 0.989457 0.000000\n",
    "108-95-2_old": "13\nPhenol; structure form HCP92; l\nC -1.046085 -0.892147 -0.000000\nC 0.257414 -1.400049 0.000000\nC 1.331097 -0.503372 -0.000000\nC 1.091601 0.874901 -0.000000\nC -0.121202 1.347367 0.000000\nC -1.230130 0.494536 0.000000\nH -1.874522 -1.578781 0.001032\nH 0.431625 -2.467932 0.001287\nH 2.336451 -0.886828 0.000832\nH 1.937051 1.553332 -0.000803\nO -0.312498 2.697885 -0.001272\nH -2.244385 0.877082 -0.000405\nH -1.251130 2.879281 -0.002039\n",
    "108-95-2": "13\nPhenol; structure form HCP92 revised structure; l\nC 4.555420 5.661760 4.489060\nC 4.584960 2.843420 4.498910\nC 3.378760 3.548270 4.498890\nC 3.359200 4.943960 4.500070\nC 5.758650 4.948670 4.502800\nC 5.789840 3.550210 4.500330\nH 4.622340 1.760920 4.495600\nH 3.631910 7.321310 4.523190\nH 2.474610 2.963630 4.498170\nH 2.384290 5.419790 4.496670\nH 6.695040 5.492670 4.498160\nH 6.727200 3.021210 4.497160\nO 4.537700 7.024010 4.500450\n",
    "60-29-7": "15\nEthoxy ethane; experimental structure from kuc98; l\nO  0.0000  0.0000  0.2696\nC  0.0000  1.1705 -0.5184\nC  0.0000 -1.1705 -0.5184\nC  0.0000  2.3716  0.4082\nC  0.0000 -2.3716  0.4082\nH -0.8879  1.1870 -1.1676\nH  0.8879  1.1870 -1.1676\nH  0.8879 -1.1870 -1.1676\nH -0.8879 -1.1870 -1.1676\nH  0.0000  3.2961 -0.1729\nH  0.0000 -3.2961 -0.1729\nH  0.8840  2.3552  1.0456\nH -0.8840  2.3552  1.0456\nH -0.8840 -2.3552  1.0456\nH  0.8840 -2.3552  1.0456\n",
    "62-53-3": "14\nAniline; Structure from HCP92; l\nC -1.086143 -0.870526 0.000000\nC 0.210840 -1.375886 0.000000\nC 1.296883 -0.505360 0.000000\nC 1.086142 0.870526 0.000000\nC -0.210841 1.375786 0.000000\nC -1.296983 0.505261 0.000000\nH -1.932704 -1.548451 0.000000\nH 0.374980 -2.447943 0.000000\nH 2.307474 -0.899002 0.000000\nH 1.932383 1.548850 0.000000\nH -2.307832 0.898242 0.000000\nN -0.428501 2.790157 0.000000\nH 0.323567 3.340664 0.346569\nH -1.327164 3.070170 0.333980\n",
    "629-20-9": "16\nCyclooctatetraene; experimental structure from krummli08, l\nC -0.2627 -1.6663  0.3833\nC  1.0331 -1.3409  0.3827\nC -1.0297  1.3407  0.3845\nC  0.2630  1.6666  0.3835\nC -1.3424 -1.0272 -0.3796\nC  1.6770 -0.2635 -0.3837\nC -1.6841  0.2615 -0.3911\nC  1.3455  1.0312 -0.3823\nH -0.5690 -2.5492  0.9659\nH  1.7214 -1.9720  0.9656\nH -1.7184  1.9710  0.9702\nH  0.5706  2.5481  0.9698\nH -1.9688 -1.7275 -0.9575\nH  2.5603 -0.5678 -0.9650\nH -2.5689  0.5611 -0.9710\nH  1.9751  1.7236 -0.9613 \n",
    "108-88-3": "15\nTuloene; structure from HCP92; l\nC -1.091600 -0.874900 0.000000\nC 0.211900 -1.382800 0.000000\nC 1.303400 -0.507900 0.000000\nC 1.091600 0.874900 0.000000\nC -0.211900 1.382700 0.000000\nC -1.303500 0.507800 0.000000\nH -1.957699 -1.569142 0.000000\nH 0.380087 -2.479984 0.000000\nH 1.957699 1.569142 0.000000\nH -0.380072 2.479886 0.000000\nH -2.337729 0.910876 0.000000\nC 2.723481 -1.061023 0.000000\nH 3.397184 -0.350813 0.523286\nH 2.739366 -2.039786 0.523326\nH 3.068240 -1.195306 -1.046523\n",
    "66-22-8": "12\nUracil; TM def2-QZVP optimized pbe structure; l\nH     1.7181113   -0.0000002   -2.1244327 \nC     1.1438667   -0.0000001   -1.2038069 \nC     1.7446080   -0.0000000    0.0097511 \nH     2.8271212   -0.0000001    0.1298240 \nC    -0.3090145    0.0000000   -1.2970091 \nO    -0.9725128   -0.0000001   -2.3249912 \nN    -0.9523747    0.0000001   -0.0343227 \nH    -1.9701793    0.0000002   -0.0511160 \nC    -0.3751368    0.0000001    1.2248796 \nO    -0.9991833   -0.0000001    2.2733353 \nN     1.0223980    0.0000001    1.1766034 \nH     1.4813173    0.0000001    2.0806653 \n",
    "71-30-7": "13\nCytosine; TM def2-QZVP pbe optimized; l\nH    -2.0638946    1.7581987   -0.0048606 \nC    -1.1537688    1.1649096    0.0014411 \nC     0.0751951    1.7542640    0.0005025 \nH     0.2153507    2.8348293    0.0000624 \nN     1.1910674    0.9892600   -0.0011780 \nH     2.1178150    1.4035814   -0.0026938 \nC     1.1663885   -0.4437838   -0.0000952 \nO     2.2347026   -1.0410121   -0.0014923 \nN    -0.0756682   -1.0234974    0.0071744 \nC    -1.1667470   -0.2723068    0.0031814 \nN    -2.3615150   -0.9269214   -0.0250037 \nH    -2.3405765   -1.9329874    0.0897975 \nH    -3.2271802   -0.4358623    0.1455697 \n",
    "106-97-8": "14\nButhane; TM def2-QZVP pbe optimized structure; m\nC    -0.5698992    0.0010721   -0.5106280 \nC    -1.9574388   -0.0010272    0.1310139 \nC     0.5699130    0.0010684    0.5106802 \nH    -2.1019869   -0.8898549    0.7620598 \nH    -2.1011958    0.8826120    0.7695094 \nH    -2.7542204    0.0025732   -0.6252247 \nH    -0.4643309   -0.8776931   -1.1680830 \nH    -0.4658369    0.8818734   -1.1655696 \nH     0.4658491    0.8818113    1.1656828 \nH     0.4644008   -0.8777422    1.1680681 \nC     1.9574251   -0.0010272   -0.1310574 \nH     2.7542844    0.0069401    0.6250702 \nH     2.1035257   -0.8919972   -0.7587064 \nH     2.0995094    0.8804516   -0.7729107 \n",
    "65-71-4": "15\nThymine; TM optimized def2-QZVP pbe; l\nC     0.5676469    0.0000748   -1.1919386 \nC     1.5304659    0.0002388   -0.2353148 \nH     2.5926489    0.0003404   -0.4816273 \nC    -0.8332110    0.0000968   -0.7623841 \nO    -1.8038589    0.0001316   -1.5105895 \nN    -1.0160281    0.0001744    0.6372389 \nH    -1.9813565    0.0002431    0.9610129 \nC    -0.0490742    0.0002212    1.6309425 \nO    -0.2873837   -0.0005911    2.8289393 \nN     1.2432728    0.0002954    1.1110067 \nH     1.9834347    0.0004963    1.8040615 \nC     0.8557798   -0.0004840   -2.6602502 \nH     1.9357502   -0.0001999   -2.8533074 \nH     0.4097225    0.8785617   -3.1465795 \nH     0.4104462   -0.8804337   -3.1455944 \n",
    "73-24-5": "15\nAdenine; TM optimized def2-QZVP pbe; l\nC    -0.7843450    0.0024660    0.6854944 \nC     0.5301823    0.0009168    0.1951947 \nN    -1.9138024    0.0022584   -0.0334730 \nC    -1.6479109   -0.0019330   -1.3472026 \nH    -2.5169904   -0.0044917   -2.0095241 \nN    -0.4540095   -0.0051210   -1.9657824 \nC     0.6617763   -0.0009654   -1.2107858 \nN     1.8686255    0.0181736   -1.8266445 \nH     2.7131398   -0.0801668   -1.2801876 \nH     1.9011412   -0.0813319   -2.8324875 \nN     1.4556976   -0.0028403    1.2254693 \nH     1.1016503   -0.0046108    3.3314279 \nC     0.7173645   -0.0023279    2.3155892 \nN    -0.6380338    0.0012265    2.0572467 \nH    -1.3931273    0.0022322    2.7328517 \n",
    "73-40-5": "16\nGuanine; TM optimized def2-QZVP pbe functional; l\nC    -0.8909136    0.0022495    0.4879726 \nC     0.4648657    0.0066008    0.8428894 \nN    -1.4589692    0.0149180   -0.7432767 \nC    -0.5646417    0.0081792   -1.7097505 \nN     0.6167192    0.0012409    2.2159144 \nH    -0.8858269   -0.0099074    3.7340662 \nC    -0.6098052   -0.0051706    2.6839376 \nN    -1.5660059   -0.0044593    1.6825133 \nH    -2.5739287   -0.0104122    1.7887665 \nC     1.4519352    0.0039931   -0.2048012 \nO     2.6727207   -0.0060475   -0.1670731 \nN     0.7873252    0.0064149   -1.4864710 \nH     1.4327355   -0.0565403   -2.2706380 \nN    -0.9871915   -0.0557168   -3.0197562 \nH    -1.9795941    0.1251659   -3.1279347 \nH    -0.4045996    0.3813093   -3.7247680 \n",
    "14868-53-2": "17\nPentasilane; def2-QZVP pbe optimized; m\nSi   -0.0048335   -3.8969717   -0.5238439 \nH     1.1887767   -3.9134075   -1.4252499 \nH     0.0223993   -5.1288186    0.3252902 \nH    -1.2385400   -3.9282984   -1.3688912 \nSi    0.0104464   -1.9536989    0.7969029 \nSi   -0.0004053    0.0000130   -0.5100885 \nH    -1.1904082   -1.9457891    1.6943990 \nH     1.2296666   -1.9529052    1.6693294 \nSi   -0.0104987    1.9536727    0.7969268 \nH    -1.2088685   -0.0034247   -1.3978287 \nH     1.2088478    0.0037073   -1.3964714 \nSi    0.0051940    3.8969653   -0.5238481 \nH    -0.0191436    5.1288161    0.3253672 \nH     1.2377333    3.9262510   -1.3707348 \nH    -1.1896928    3.9156075   -1.4234610 \nH     1.1910193    1.9496877    1.6935745 \nH    -1.2290887    1.9491200    1.6702050 \n",
    "100-41-4": "18\nEthybenzene; TM optimized def2-QZVP pbe structure; l\nC    -2.2693535   -0.0000389   -0.2398724 \nC    -1.5797056   -1.2063053   -0.1005171 \nC    -0.2110721   -1.2030110    0.1749198 \nC     0.4952604    0.0000401    0.3167968 \nC    -0.2111286    1.2030501    0.1748811 \nC    -1.5797652    1.2062634   -0.1005510 \nH    -3.3389696   -0.0000690   -0.4503195 \nH    -2.1102702   -2.1538484   -0.2012097 \nH     0.3195630   -2.1509209    0.2885225 \nH     0.3194582    2.1509892    0.2884519 \nH    -2.1103902    2.1537683   -0.2012753 \nC     2.8147071   -0.0000376   -0.7147854 \nC     1.9829690    0.0000498    0.5777693 \nH     2.5926203   -0.8860058   -1.3248301 \nH     2.5930639    0.8861475   -1.3246747 \nH     3.8904081   -0.0003286   -0.4914695 \nH     2.2473809   -0.8820350    1.1803103 \nH     2.2474159    0.8821780    1.1802303 \n"
}


ccsdt_homo = {
    "orbital": "HOMO", 
    "remark": "", 
    "code": "CFOUR and Gaussian16", 
    "basis_size": "3", 
    "basis": "gaussian", 
    "code_version": "", 
    "qpe": "",
    "DOI": "https://doi.org/10.3389/fchem.2021.749779",
    "basis_name": "def2-TZVPP", 
    "calc_type": "CCSD(T)", 
    "data": {
        "106-97-8": -11.567, 
        "558-13-4": -10.41199, 
        "7722-84-1": -11.51763, 
        "65-71-4": -9.081, 
        "7782-41-4": -15.58481, 
        "108-88-3": -8.96840, 
        "100-41-4": -8.91875, 
        "23878-46-8": -9.84959, 
        "74-84-0": -12.71117, 
        "71-43-2": -9.36006, 
        "593-66-8": -9.327, 
        "7783-60-0": -12.58832, 
        "12184-80-4": -11.24071, 
        "13768-60-0": -11.08635, 
        "7758-02-3": -8.127, 
        "71-30-7": -8.77093,
        "7439-90-9": -13.94, 
        "19287-45-7": -12.25030, 
        "7732-18-5": -12.565, 
        "64-18-6": -11.421, 
        "74-98-6": -12.03385, 
        "7440-63-3": -12.260, 
        "7803-51-2": -10.523, 
        "62-53-3": -8.04342, 
        "630-08-0": -14.209, 
        "7580-67-8": -7.961, 
        "39297-88-6": -4.39181, 
        "7784-23-8": -9.71683, 
        "7783-40-6": -13.710,
        "7446-09-5": -12.29831, 
        "1590-87-0": -10.645, 
        "12190-70-4": -7.566, 
        "10028-15-6": -12.70550, 
        "7727-37-9": -15.4767, 
        "7664-41-7": -10.807, 
        "75-07-0": -10.207, 
        "1333-74-0": -16.403, 
        "7783-06-4": -10.31, 
        "7647-14-5": -9.02732, 
        "12185-09-0": -10.52987, 
        "74-85-1": -10.666, 
        "74-82-8": -14.373, 
        "73-24-5": -8.33, 
        "542-92-7": -8.71208, 
        "302-01-2": -9.68764, 
        "75-19-4": -10.865, 
        "593-60-2_old": -9.269, 
        "593-60-2": -9.84109,
        "17108-85-9": -9.771, 
        "25681-81-6": -3.94089, 
        "7693-26-7": -6.128, 
        "7803-62-5": -12.796, 
        "12187-06-3": -7.494, 
        "73-40-5": -8.034, 
        "7440-01-9": -21.32107, 
        "74-90-8": -13.72282,
        "10043-11-5": -11.98447,
        "57-13-6": -10.053, 
        "7789-24-4": -11.32066,
        "7440-59-7": -24.512, 
        "507-25-5": -9.19399, 
        "74-86-2": -11.424, 
        "60-29-7": -9.816, 
        "7664-39-3": -16.026, 
        "629-20-9": -8.39799, 
        "110-86-1": -9.72955, 
        "56-23-5": -11.49885, 
        "14868-53-2": -9.27474, 
        "66-22-8": -9.48267,
        "67-56-1": -11.042, 
        "13283-31-3": -13.27449, 
        "75-73-0": -16.22752, 
        "75-01-4": -10.093, 
        "108-95-2_old": -8.702, 
        "108-95-2": -8.64963,
        "7784-18-1": -15.32499,
        "75-15-0": -9.981, 
        "7553-56-2": -9.48436, 
        "14452-59-6": -5.266, 
        "124-38-9": -13.7112, 
        "17739-47-8": -11.81239, 
        "7440-37-1": -15.544, 
        "7782-79-8": -10.676, 
        "392-56-3": -9.93, 
        "75-02-5": -10.554, 
        "1304-56-9": -9.94321, 
        "463-58-1": -11.173, 
        "50-00-0": -10.84, 
        "64-17-5": -10.685, 
        "7726-95-6": -10.52389, 
        "7782-65-2": -12.497, 
        "7786-30-3": -11.65515, 
        "1603-84-5": -10.47061, 
        "7783-63-3": -15.41022, 
        "544-92-3": -10.854,
        "39297-86-4": -4.24246, 
        "7647-01-0": -12.593, 
        "7784-42-1": -10.398, 
        "1309-48-4": -7.90916,
        "25681-80-5": -4.07463, 
        "25681-79-2": -4.95498, 
        "7782-50-5": -11.40512
    }
}

g0w0_hf_homo = {
  "code": "MOLGW",
  "code_version": "2.E",
  "basis": "gaussian",
  "qpe": "solved",
  "DOI": "https://doi.org/10.3389/fchem.2021.749779",
  "orbital": "HOMO",
  "remark": "RIJK with AUTO",
  "basis_size": "3",
  "basis_name": "Def2-TZVPP",
  "calc_type": "gw@hf",
  "data": {
    "75-19-4": -11.2305368,
    "75-02-5": -10.7636477,
    "74-98-6": -12.5701914,
    "7803-51-2": -10.7601539,
    "60-29-7": -10.4288651,
    "73-40-5": -8.35658192,
    "7783-63-3": -16.0316473,
    "50-00-0": -11.3133176,
    "23878-46-8": -9.75926305,
    "7784-18-1": -15.5858456,
    "1309-48-4": -7.82498746,
    "25681-79-2": -4.92991808,
    "7446-09-5": -12.8883403,
    "74-85-1": -10.7060223,
    "19287-45-7": -12.7751074,
    "7647-14-5": -9.19474855,
    "39297-86-4": -4.24252074,
    "108-88-3": -9.07936335,
    "64-17-5": -11.2289359,
    "1590-87-0": -11.0667642,
    "7439-90-9": -14.0486091,
    "558-13-4": -10.8061096,
    "7664-41-7": -11.1385369,
    "7803-62-5": -13.2239795,
    "64-18-6": -11.8897917,
    "56-23-5": -11.9611963,
    "7782-65-2": -12.8820956,
    "7553-56-2": -9.6713214,
    "7784-23-8": -9.97956835,
    "302-01-2": -10.1090177,
    "7784-42-1": -10.5982697,
    "630-08-0": -14.999212,
    "75-15-0": -10.2752377,
    "507-25-5": -9.49823318,
    "14452-59-6": -5.28761143,
    "100-41-4": -9.04536646,
    "62-53-3": -8.26916184,
    "12185-09-0": -10.5030503,
    "542-92-7": -8.78808362,
    "7783-40-6": -13.7883045,
    "7693-26-7": -6.27565455,
    "67-56-1": -11.5106517,
    "12187-06-3": -7.12181345,
    "17739-47-8": -12.3172275,
    "74-86-2": -11.5367367,
    "7789-24-4": -11.3046419,
    "7440-37-1": -15.7287019,
    "66-22-8": -10.0171033,
    "1603-84-5": -10.7052901,
    "7440-01-9": -21.3474972,
    "74-90-8": -13.8188036,
    "7783-60-0": -13.2470766,
    "7647-01-0": -12.7700318,
    "7664-39-3": -16.1675236,
    "17108-85-9": -9.91300774,
    "124-38-9": -14.1607081,
    "25681-80-5": -4.02647921,
    "106-97-8": -12.1316673,
    "13283-31-3": -13.6345296,
    "7783-06-4": -10.4803831,
    "463-58-1": -11.5022976,
    "7782-79-8": -11.0336329,
    "7580-67-8": -8.15233216,
    "10043-11-5": -11.6850439,
    "39297-88-6": -4.41335012,
    "75-07-0": -10.7507378,
    "7786-30-3": -11.8904667,
    "57-13-6": -10.6121237,
    "14868-53-2": -9.76794696,
    "7782-50-5": -11.733006,
    "1333-74-0": -16.4746358,
    "71-43-2": -9.44701351,
    "75-01-4": -10.3095695,
    "629-20-9": -8.60231915,
    "7782-41-4": -16.2634732,
    "12184-80-4": -11.5537833,
    "7440-59-7": -24.6044038,
    "108-95-2_old": -8.82944531,
    "108-95-2": -8.95045374,
    "25681-81-6": -3.84937685,
    "12190-70-4": -7.18166439,
    "7722-84-1": -12.0026476,
    "74-84-0": -13.1385927,
    "7732-18-5": -12.8149513,
    "7727-37-9": -16.2961489,
    "10028-15-6": -13.4928835,
    "75-73-0": -16.7908246,
    "7440-63-3": -12.3234026,
    "593-60-2_old": -9.43712268,
    "593-60-2": -10.0392676,
    "13768-60-0": -11.260267,
    "71-30-7": -9.20527097,
    "392-56-3": -10.5493927,
    "73-24-5": -8.62673442,
    "544-92-3": -10.8314476,
    "7726-95-6": -10.7907363,
    "593-66-8": -9.46048413,
    "110-86-1": -9.83795375,
    "7758-02-3": -8.24074958,
    "1304-56-9": -9.75811412,
    "74-82-8": -14.7316617,
    "65-71-4": -9.60739852
  }
}

########################################################################
def create_official_json(filename = "", data = {}, **kwargs):
    """
        dump a json file named filename in the official format of the GW100 website
    """
    dict_gw100 = dict()
    dict_gw100["code"]= "MOLGW"
    dict_gw100["code_version"]= __version__
    dict_gw100["basis"]= "gaussian"
    dict_gw100["qpe"]= "solved"
    dict_gw100["DOI"]= "unpublished"

    dict_gw100.update(kwargs)
    dict_gw100["data"] = data

    if filename != "":
        with open(filename, 'w') as json_file:
            json.dump(dict_gw100,json_file,indent=2,separators=(',', ': '))
    return dict_gw100


########################################################################
def correlation_plot(data1,data2,labels=False):
    fig = plt.figure(figsize=(6, 6))
    xxx = []
    yyy = []
    cass = []
    shared_keys = set(data1.keys()) & set(data2.keys())
    for cas in shared_keys:
        xxx.append(data1[cas])
        yyy.append(data2[cas])
        cass.append(cas)
    xmin = min(xxx+yyy) - 0.1
    xmax = max(xxx+yyy) + 0.1
    plt.plot([xmin,xmax], [xmin,xmax], '-', color='black')
    plt.scatter(xxx, yyy, marker='o')
    if labels:
        delta = abs(xmax-xmin) * 0.02
        for x, y, cas in zip(xxx,yyy,cass):
            plt.text(x+delta, y-delta, chemical_formulas[cas])
    plt.xlabel("HOMO (eV)")
    plt.ylabel("HOMO (eV)")
    return fig

########################################################################
def diff(data1,data2):
    """
        Returns a dictionary containing data2 - data1
    """
    errors = dict()
    shared_keys = set(data1.keys()) & set(data2.keys())
    for cas in shared_keys:
        errors[cas] = data2[cas] - data1[cas]
    return errors

########################################################################
def mae_mse_max(data1,data2):
    """
        Returns MAE, MSE, Max errors with data1 being the reference
    """
    errors = list( diff(data1,data2).values() )

    ndata = len(errors)
    mse = sum(errors) / float(ndata)
    mae = sum(np.abs(errors)) / float(ndata)
    mxe = np.max(np.abs(errors))
    return mae, mse, mxe
