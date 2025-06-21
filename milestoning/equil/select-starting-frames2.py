#!/usr/bin/env python
import sys
import os
import numpy as np
import math
import mdtraj as md

#load the list of anchors
name=sys.argv[1]
nanchor=int(sys.argv[2])
anchorsfname=sys.argv[3]
ianchor=int(sys.argv[4])
janchor=int(sys.argv[5])
#psf=sys.argv[6]
#refpdb=sys.argv[7]
#dcdname=sys.argv[8]
#scratch='/scratch/jmsc87/pyk2/milestoning/'
#storage='/storage/hpc/group/wong/scratch/Justin/pyk2/milestoning'
scratch=os.environ['scratch']
storage=os.environ['storage']
if (name.find('charmm')>=0):
	psf='../solvate-charmm/psf/{1}-milestone-{2:d}-{3:d}-solvated.psf'.format(storage,name,ianchor,janchor)
	top=md.load_psf(psf)
else:
	psf='../solvate/prmtop/{1}-milestone-{2:d}-{3:d}-solvated.prmtop'.format(storage,name,ianchor,janchor)
	top=md.load_prmtop(psf)
if ((name=='pyk2-cpd1-amber') or (name=='pyk2-cpd10-amber')):
	refpdb=os.environ['HOME']+'/pyk2/umbrella-amber/ref/{0}-tramd-window-1-ref.pdb'.format(name)
elif (name.find('charmm')>=0):
	refpdb='../ref/{0}-window-1-ref.pdb'.format(name)
else:
	refpdb=os.environ['HOME']+'/pyk2/umbrella-amber/ref/{0}-window-1-ref.pdb'.format(name)

dcdname='{0}/equil-dcd/{1}-milestone-{2:d}-{3:d}-equil.dcd'.format(storage,name,ianchor,janchor)
nsample=int(sys.argv[6])
#start_frame_dir=sys.argv[7]
start_frame_dir='{0}/start-frames/'.format(scratch)


#load the list of anchors
#format: ianchor x y z energy, assumed in order
input=open(anchorsfname,'r')
anchors=np.full((nanchor,3),np.nan,dtype=np.float64)
for line in input:
	words=line.split()
	i=int(words[0])-1
	for k in range(0,3):
		anchors[i,k]=float(words[k+1])
input.close()
anchors=np.array(anchors)
#print(anchors)
nanchor=anchors.shape[0]
#print(nanchor)

#determine if any previous frames have been selected
prevsample=[]
try:
	infofname='data/{0}-milestone-{1:d}-{2:d}-selected-frames'.format(name,ianchor,janchor)
	input=open(infofname,'r')
except OSError as err:
	#not found, skip it
	pass
else:
	print("reading previous sample from ",infofname,flush=True)
	for line in input:
		words=line.split()
		prevsample.append(int(words[1])-1)
prevsample=np.sort(np.array(prevsample))
nprevsample=len(prevsample)
if (nprevsample>0):
	print("previous sample: ",prevsample,flush=True)



#load the trajectory and compute the aligned ligand COM
reftraj=md.load(refpdb)
print(len(list(top.atoms)))
#backbone_atoms=top.select('backbone')
backbone_atoms=top.select('name N or name CA or name C')
print("loading trajectory ",dcdname,flush=True)
traj=md.load(dcdname,top=top)
traj=traj.superpose(reftraj,atom_indices=backbone_atoms)
print("finished loading trajectory",flush=True)
ligand_atoms=top.select('resname CPD')
com_data=md.compute_center_of_mass(traj.atom_slice(atom_indices=ligand_atoms))
com_data=com_data*10 #convert from nm to A
print(com_data.shape)
print(com_data)

nframes=com_data.shape[0]

#determine the distances from the COM in each frame to all the anchors

dist=np.full((nframes,nanchor),np.nan,dtype=np.float64)
for kanchor in range(0,nanchor):
	diff=com_data[:,:]-np.expand_dims(anchors[kanchor,:],axis=0)
	dist[:,kanchor]=np.sqrt(np.einsum('ij,ij->i',diff,diff)) #intended shape (nframes)	

#choose the two anchors that are closest for each frame (+1 for 1-based anchor numbering)
info=np.argsort(dist,axis=1)+1
milestones=np.sort(info[:,0:2],axis=1) #so that the lower-numbered milestone is first
#print(milestones)
#these are the frames from the second half of the trajectory, which are closer to the two anchor
#than any others
start=int(np.ceil(nframes/2))
mask=np.all([(milestones[start:,0]==ianchor),(milestones[start:,1]==janchor)],axis=0)
#print(milestones[:,0]==1)
print(len(mask))
#select only frames from the second half of the trajectory
indices=np.where(mask)[0]+start
#remove any frames that were previously sampled
indices=np.setdiff1d(indices,prevsample)
if (len(indices)<0.1*nsample):
	print('warning: found too few frames, trying last 75% of the trajectory')
	start=int(np.ceil(nframes/4))
	mask=np.all([(milestones[start:,0]==ianchor),(milestones[start:,1]==janchor)],axis=0)
	indices=np.where(mask)[0]+start
if (len(indices)<0.1*nsample):
	print('warning: found too few frames, trying last 90% of the trajectory')
	start=int(np.ceil(0.1*nframes))
	mask=np.all([(milestones[start:,0]==ianchor),(milestones[start:,1]==janchor)],axis=0)
	indices=np.where(mask)[0]+start


print(indices)
with_replacement=(len(indices)<nsample)
if (with_replacement):
	print('warning: sampling frames with replacement because of insufficient frames meeting criteria')
	print('number of usable frames: ',len(indices),' samples ',nsample)
sample=np.random.choice(indices,size=nsample,replace=with_replacement)
sample=np.sort(sample)
print(sample)

#sys.exit(-1)
#reload the trajectory without superposition
#otherwise changing the orientation of the periodic unit cell leads to bad overlaps
#when trying to start the simulation
del traj
traj=md.load(dcdname,top=top)

#write the sample to a data file
infofname='data/{0}-milestone-{1:d}-{2:d}-selected-frames'.format(name,ianchor,janchor)
if (nprevsample>0):
	output=open(infofname,'w')
else:
	output=open(infofname,'a')
#for id in range(0,nprevsample):
#	iframe=prevsample[id]
#	output.write('{0:d} {1:d}\n'.format(id+1,iframe+1))
for id in range(0,nsample):
	iframe=sample[id]
	print('writing frame {0:d} id {1:d}'.format(iframe+1,id+nprevsample+1),flush=True)
	output.write('{0:d} {1:d}\n'.format(id+nprevsample+1,iframe+1))
	#print('unit cell: ',traj.unitcell_lengths[iframe,:])
	outfname='{0}/{1}-milestone-{2:d}-{3:d}-start-id-{4:d}.pdb'.format(start_frame_dir,name,ianchor,janchor,id+nprevsample+1)
	newtraj=md.formats.PDBTrajectoryFile(outfname,mode='w')
	newtraj.write(traj.xyz[iframe,:,:]*10,top)
	newtraj.close()

output.close()


##for frame 
#print(dist[0,info[0,0]],dist[0,info[0,1]])
