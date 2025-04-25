#!/usr/bin/python3
###################################
# Check whether the Fortran module dependencies
# are correctly stated in the Makefile
###################################

import sys, os, time, shutil, subprocess

fout = open('tmp','w')
subprocess.call(["grep 'use m_'  ../src/*.f90"],shell=True,stdout=fout)
fout.close()

fortran_file_list = []
module_file_list = []
tmp_list = []
fout = open('tmp','r')
for line in fout:
    parse0 = line.split("src/")
    parse1 = parse0[1].split(":")
    parse2 = parse1[1].split("use ")
    if not parse1[0] in fortran_file_list  :
        if len(tmp_list) > 0 :
            module_file_list.append(tmp_list)
        fortran_file_list.append(parse1[0])
        tmp_list = []
    module_name = parse2[1].split(",")[0].rstrip()
    if not module_name in tmp_list :
        tmp_list.append(module_name)

module_file_list.append(tmp_list)
fout.close()
os.remove('tmp')

#print(len(fortran_file_list))
#print(len(module_file_list))
#print(fortran_file_list)
#print(module_file_list)

for i in range(len(fortran_file_list)) :
    object_file = fortran_file_list[i].split(".")[0]+'.o:'
    fout = open('../src/Makefile','r')
    for line in fout :
        if object_file in line :
            string = line.split(':')[1]
            current_line = line
            while '\\' in current_line :
                current_line = next(fout)
                string = string + current_line

            string = string.replace('\\',' ')
            parse = string.split()
            #Remove the last one (itself)
            parse.pop()
            #Remove the file that are not modules
            parse = [ x for x in parse if "m_" in x]

            print(object_file)
            for j in range(len(module_file_list[i])) :
                if not module_file_list[i][j]+'.o' in parse :
                    print('{0:s} is missing'.format((module_file_list[i][j]+'.o').ljust(24)))
            for j in range(len(parse)) :
                if not parse[j].split('.o')[0] in module_file_list[i] :
                    print('{0:s} is not needed'.format((parse[j]).ljust(24)))
            print()
    fout.close()

exit(0)

