#!/usr/bin/env python
import sys
import os
import numpy as np
import math
import mdtraj as md
#import pdb

#load the list of anchors
nanchor=int(sys.argv[1])
anchorsfname=sys.argv[2]
ianchor=int(sys.argv[3])
janchor=int(sys.argv[4])
psf=sys.argv[5]
refpdb=sys.argv[6]
dcdname=sys.argv[7]

#load the list of anchors
#format: ianchor x y z energy, assumed in order
input=open(anchorsfname,'r')
anchors=np.full((nanchor,3),np.nan,dtype=np.float64)
for line in input:
	words=line.split()
	kanchor=int(words[0])-1
	for l in range(0,3):
		anchors[kanchor,l]=float(words[l+1])
input.close()
anchors=np.array(anchors)
#print(anchors)
nanchor=anchors.shape[0]
#print(nanchor)

#load the trajectory and compute the aligned ligand COM
reftraj=md.load(refpdb)
top=md.load_prmtop(psf)
#backbone_atoms=top.select('backbone')
backbone_atoms=top.select('name N or name CA or name C')
traj=md.load(dcdname,top=top)
traj=traj.superpose(reftraj,atom_indices=backbone_atoms)

ligand_atoms=top.select('resname CPD')
com_data=md.compute_center_of_mass(traj.atom_slice(atom_indices=ligand_atoms))
com_data=com_data*10 #convert from nm to A
#print(com_data.shape)
#print(com_data)

nframes=com_data.shape[0]

#determine the distances from the COM in each frame to all the anchors

dist=np.full((nframes,nanchor),np.nan,dtype=np.float64)
for kanchor in range(0,nanchor):
	diff=com_data[:,:]-np.expand_dims(anchors[kanchor,:],axis=0)
	dist[:,kanchor]=np.sqrt(np.einsum('ij,ij->i',diff,diff)) #intended shape (nframes)	

#choose the two anchors that are closest for each frame -- +1 for conversion to 1-based 
closest_anchor=np.argsort(dist,axis=1)+1
#closest_anchor[iframe,0] is the closest anchor (zero based), closest_anchor[iframe,1] is the second closest, etc.
#pdb.set_trace()
#search for the first frame for which the closest anchor is not one of the two starting anchors for this milestone.
#if we find such a frame, then the trajectory has left the two voronoi cells for these anchors, 
#has struck another milestone, and should stop.  The entries in closest_anchor for this frame correspond to the 
#new milestone that has been struck, (which may not be on the list).
#report the frame number and new milestone back to the calling script.

#mask is "false" if frame is still on milestone, "true" otherwise
mask=np.all([(closest_anchor[:,0]!=ianchor),(closest_anchor[:,0]!=janchor)],axis=0)
#are there any frames off milestone?
if (np.any(mask)):
	#we moved off this milestone! -- frame is 0-based
	first_frame=np.where(mask)[0][0]
	#these entries correspond to the new milestone.  Sort to ensure that the lower anchor is first
	new_milestone=np.sort(closest_anchor[first_frame,0:2])
	#write out the stop message with a 1-based frame (the calling script assumes this)
	result='STOP {0:d} {1:d} {2:d}'.format(first_frame+1,new_milestone[0],new_milestone[1])
else:
	#we stayed on milestone.  Continue the simulation
	result='CONTINUE'
	
print(result)

