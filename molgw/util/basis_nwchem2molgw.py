#!/usr/bin/python

import os,sys

if len(sys.argv) > 1:
        print 'Opening file: ',sys.argv[1]
else:
        print 'Give a file name as an argument'
        sys.exit(0)

filename=sys.argv[1]

if not os.path.exists(filename):
        print 'File ',filename,' does not exist'
        sys.exit(0)


nwfile = open(filename,'r')

list_am='SPDFGHI'
am=0
#alpha = []
coeff = []
firstline=0
total_function=0
final_string=''

for line in nwfile:
	parsing=line.split()
	# test whether the line contains coefficents or characters
	try:
		alpha_tmp=float(parsing[0])
		alpha.append(alpha_tmp)
		nfunctions=len(parsing)
		coeff_list = []
		for function in range(1,nfunctions):
			coeff_list.append(float(parsing[function]))
		coeff.append(coeff_list)
	except ValueError:
		# Dump out the values only if it is not the first line
		if firstline==0:
			firstline=1
		else:
			for function in range(0,len(coeff[0])):
				total_function+=1
				#
				# retain only those primitives with non zero coefficient
				primitive_local=0
				for primitive in range(0,len(alpha)):
					if abs(coeff[primitive][function]) > 1.0e-8:
						primitive_local+=1
#				print primitive_local,am
				final_string+='%2d   %2d \n' % (primitive_local,am)
		 		for primitive in range(0,len(alpha)):
					if abs(coeff[primitive][function]) > 1.0e-8:
						string='%14.7f    ' %  alpha[primitive] 
						string+='%14.7E    ' % coeff[primitive][function]
#						print string
						final_string+=string+'\n'

		# if the file is finished, exit loop
		if parsing[0]=='END':
			print 'Found %2d functions' % total_function
			string='%2d \n' % total_function
			final_string=string+final_string
			nwfile.close()
			print 'Closing file: ',filename
			outfilename=filename.partition('.')
			print 'Opening file: ',outfilename[0]
			outfile=open(outfilename[0],'w')
			outfile.write(final_string)
			print 'Closing file: ',outfilename[0]
			outfile.close()
			sys.exit(1)
		# read the angular momentum and transform it to the numerical version S->0, P->1, D->2
 		am=list_am.find(parsing[1])
		alpha = []
		coeff = []


