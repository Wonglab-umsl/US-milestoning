#!/usr/bin/env python3.6

#must use python 2.7.16
import sys
import os
import pyemma
import numpy as np
import math
#import matplotlib
#matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 22})
#import matplotlib.pyplot as plt
import pickle
import mdtraj as md
from operator import itemgetter

name=sys.argv[1]
iter=int(sys.argv[2])
nseg=int(sys.argv[3])
ncontinue=int(sys.argv[4])
#alpha=float(sys.argv[4])
logfname=sys.argv[5]
#the environment variables must be exported from the running script
scratch=os.environ.get('scratch')
#logfile=os.environ.get('logfile')
#todo: check they are not None and bail with error message
if (scratch is None):
	print('need to specify scratch directory for trajectories')
	sys.exit(-1)


#print('setting up featurizer',flush=True)
psf='prmtop/{0}-solvated.prmtop'.format(name)
refpdb='pdb/{0}-solvated-ref.pdb'.format(name)
reftraj=md.load(refpdb)
top=md.load_prmtop(psf)
dcdlist=['{0}/dcd/{1}-{2:d}-{3:d}.dcd'.format(scratch,name,iter,seg) for seg in range(1,nseg+1)]
#print(dcdlist)
traj_list=[]
backbone_atoms=top.select('backbone')
#print(atom_indices)
for dcdname in dcdlist:
	traj=md.load(dcdname,top=top)
	traj=traj.superpose(reftraj,atom_indices=backbone_atoms)
	traj_list.append(traj)

ligand_atoms=top.select('resname CPD')
final_com_distances=np.full(nseg,np.nan)
orig_com=md.compute_center_of_mass(reftraj.atom_slice(atom_indices=ligand_atoms))[0]
#print(orig_com)
#sys.exit(-1)
for (iseg,traj) in enumerate(traj_list):
	#print(iseg,traj.n_frames)
	com_data=md.compute_center_of_mass(traj.atom_slice(atom_indices=ligand_atoms))
	final_com=com_data[-1,:]
	#print(final_com.shape)
	diff=final_com[:]-orig_com[:]
	#10.0 converts from nm to A
	dist=np.sqrt(np.dot(diff,diff))*10.0
	#print(iseg,final_com,dist)
	final_com_distances[iseg]=dist


#the distances are in A
#print(final_com_distances)
mindist=np.min(final_com_distances)
maxdist=np.max(final_com_distances)
print('iteration {0:d} ligand COM distance range: {1:.2f} - {2:.2f} A'.format(iter,mindist,maxdist),flush=True)

order=np.flip(np.argsort(final_com_distances))
print(final_com_distances)
print(order)

#this is the number of simulations to start from each parent
nsim_to_start=nseg//ncontinue
if (nseg%ncontinue!=0):
	nsim_to_start+=1
numseg=np.zeros(nseg,dtype=np.int32)
count=nseg
i=0
while (count>0):
	actual=min(count,nsim_to_start)
	numseg[order[i]]=actual
	count-=actual
	i+=1


#print(numseg)
print('parent segment, final distance (A), number of segments to start',flush=True)
for iseg in range(0,nseg):
	print('{0:3d} {1:.2f} {2:d}'.format(iseg+1,final_com_distances[iseg],numseg[iseg]),flush=True)
assert(np.sum(numseg)==nseg)
#numseg now contains the number of segments to be started from each segment.  
#sys.exit(-1)
#time to append to the logfile which contains the information on which segments started from which 
#format: name, next iteration, next segment, iteration, segment to start from
nextiter=iter+1
logfile=open(logfname,'a')
nextseg=1
for iseg in range(0,nseg):
	for i in range(0,numseg[iseg]):
		line='{0} {1:d} {2:d} {3:d} {4:d}\n'.format(name,nextiter,nextseg,iter,iseg+1)
		logfile.write(line)
		nextseg+=1
logfile.close()
sys.exit(0)
