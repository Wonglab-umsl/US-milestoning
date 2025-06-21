#!/usr/bin/env python
import sys
import os
import numpy as np
import math
import mdtraj as md
import pdb

name=sys.argv[1]
milestonefile=sys.argv[2]
refpdb=sys.argv[3]
distance_cutoff=float(sys.argv[4])

initial_milestones=[]
input=open(milestonefile,'r')
for line in input:
	words=line.split()
	ianchor=int(words[1])
	janchor=int(words[2])
	initial_milestones.append((ianchor,janchor))
input.close()
print(initial_milestones)
nmilestone=len(initial_milestones)
#autodetect product milestones based on COM cutoff

print("Detecting product milestones with COM distance cutoff ",distance_cutoff," A ",flush=True)
psf='prmtop/{0}-solvated.prmtop'.format(name)
#load the trajectory and compute the aligned ligand COM
reftraj=md.load(refpdb)
top=md.load_prmtop(psf)
#backbone_atoms=top.select('backbone')
backbone_atoms=top.select('name N or name CA or name C')
ligand_atoms=top.select('resname CPD')
ref_com=md.compute_center_of_mass(reftraj.atom_slice(atom_indices=ligand_atoms))
ref_com=ref_com[0]*10 #convert from nm to A
#autodetect product milestones based on COM cutoff
outfname='{0}-product-milestones'.format(name)
output=open(outfname,'w')
for milestone in initial_milestones:
	ianchor=milestone[0]
	janchor=milestone[1]
	print("examining milestone {0:d} {1:d}".format(ianchor,janchor),flush=True)
	pdbfname='../ref/{0}-milestone-{1:d}-{2:d}-solvated-ref.pdb'.format(name,ianchor,janchor)
	traj=md.load(pdbfname)
	traj=traj.superpose(reftraj,atom_indices=backbone_atoms)
	com_data=md.compute_center_of_mass(traj.atom_slice(atom_indices=ligand_atoms))
	com_data=com_data[0]*10 #convert from nm to A
	#pdb.set_trace()
	diff=com_data[:]-ref_com[:]
	dist=np.sqrt(np.dot(diff,diff))
	if (dist>distance_cutoff):
		line='{0:d} {1:d}\n'.format(ianchor,janchor)
		output.write(line)

output.close()

