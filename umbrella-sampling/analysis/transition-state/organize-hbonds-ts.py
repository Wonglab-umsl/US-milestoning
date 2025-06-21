#!/usr/bin/env python

import sys
import numpy as np
import math

aadict={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','GLH':'E','PHE':'F','GLY':'G','HID':'H','HIE':'H','HIP':'H','ILE':'I','LYS':'K',
        'LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}


def splitname(name): #,cpd):
	words=name.split('-')
	resname=words[0]
	if (resname=="CPD272"):
		resname="ligand" # \\textbf\{{}"
	else:
		resaa=aadict[resname[0:3]]
		ires=int(resname[3:])+419
		resname=resaa+str(ires)
	atomname=words[2]
	return (resname,atomname)


#def getres(key):
	
#def createkey(donor,acceptor):


#name=sys.argv[1]
names=['pyk2-cpd1-amber','pyk2-cpd2-amber','pyk2-cpd8-amber']
states=['bound','transition']
min_occ=1.0
if (len(sys.argv)>1):
	min_occ=float(sys.argv[1])
data={}
hbondlist=[]
#data=[]
for name in names:
	for state in states:
		fname='an/hbonds-counts-{0}-{1}-state'.format(name,state)
		input=open(fname,'r')
		for line in input:
			words=line.split()
			if ((words[0]=='Found') or (words[0]=='donor')):
				continue
			#donkey=splitname(words[0])
			#acckey=splitname(words[1])
			hbondkey=splitname(words[0])+splitname(words[1])
			#exclude non-conventional hydrogen-bonds  (carbon, sulfur, etc.)
			if ((hbondkey[1][0] in 'NOF') and (hbondkey[3][0] in 'NOF')):
				#if (hbondkey not in hbondlist):
				#	hbondlist.append(hbondkey)
				key=(name,state,)+hbondkey
				occ=float(words[2].replace('%',''))
				if (occ>=min_occ):
					data[key]=occ
					if (hbondkey not in hbondlist):
						hbondlist.append(hbondkey)

#print hbondlist
#try printing out in latex format.
cpd8_missing_atoms=['O32','O33','N30']
print('\\hline')
for hbondkey in hbondlist:
	line='{0} & {1} & {2} & {3} &'.format(hbondkey[0],hbondkey[1],hbondkey[2],hbondkey[3])
	for name in names:
		for state in states:
			key=(name,state,)+hbondkey
			if ((name=='pyk2-cpd8-amber') and ((hbondkey[1] in cpd8_missing_atoms) or (hbondkey[3] in cpd8_missing_atoms))):
				occ='n/a'
			else:
				if (key in data):
					occ='{0:.2f}'.format(data[key])
				else:
					occ=''
			line=line+' '+str(occ)
			if ((name==names[-1]) and (state==states[-1])):
				line=line+' \\\\'
			else:
				line=line+' &'
	print(line)
