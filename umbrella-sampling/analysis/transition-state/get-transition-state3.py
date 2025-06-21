#!/usr/bin/env python
import sys
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 17})
import matplotlib.pyplot as plt
from datetime import datetime
import pdb
import itertools as it
import mdtraj as md

#aadict={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','GLH':'E','PHE':'F','GLY':'G','HID':'H','HIE':'H','HIP':'H','ILE':'I','LYS':'K',
#        'LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}


windowcounts={'pyk2-cpd1-amber':421,'pyk2-cpd2-amber':353,'pyk2-cpd8-amber':187,'pyk2-cpd10-amber':452,'pyk2-cpd6-amber':368}
#windowcounts={'pyk2-cpd1-amber':300,'pyk2-cpd2-amber':50}
dir='/mnt/stor/ceph/scratch/jmsc87/pyk2/umbrella-amber3/dcd-joined/'
scratch='/lustre/scratch/jmsc87/pyk2/umbrella-amber/'
name=sys.argv[1]
nanchor=int(sys.argv[2])
nwindow=windowcounts[name]

if (name=='pyk2-cpd1-amber'):
	anchorsfname='{0}-pathway-4'.format(name)
	committorsfname='{0}-committors-product-milestones'.format(name)
elif (name=='pyk2-cpd2-amber'):
	anchorsfname='{0}-pathway-15'.format(name)
	committorsfname='{0}-committors-product-milestones2'.format(name)
elif (name=='pyk2-cpd8-amber'):
	anchorsfname='{0}-pathway-2-inserted'.format(name)
	committorsfname='{0}-committors-product-milestones2'.format(name)
else:
	assert False, 'unrecognized name'



input=open(anchorsfname,'r')
anchors=np.full((nanchor,3),np.nan,dtype=np.float64)
for line in input:
	words=line.split()
	kanchor=int(words[0])-1
	if (kanchor<nanchor):
		for l in range(0,3):
			anchors[kanchor,l]=float(words[l+1])
input.close()
anchors=np.array(anchors)

bound_center=anchors[0,:]
input=open(committorsfname,'r')
#committors_by_anchor=np.full(nanchor,np.nan,dtype=np.float64)
ts_milestones=[]
for line in input:
	words=line.split()
	#format: imilestone ianchor janchor committor
	ianchor=int(words[1])
	janchor=int(words[2])
	if (int(words[0])==1):
		bound_milestone=(ianchor,janchor)
	c=float(words[3])
	if ((c>0.4) and (c<0.6)):
		#pdb.set_trace()
		#midpoint=0.5*(anchors[ianchor-1,:]+anchors[janchor-1,:])
		print('transition state milestone: ianchor, janchor, committor = {0:d} {1:d} {2:.3f}'.format(
			ianchor,janchor,c),flush=True)
		ts_milestones.append((ianchor,janchor))

        
	#if ((ianchor<nanchor) and ((janchor-ianchor)==1)):
		#committors_by_anchor[ianchor-1]=c


#pdb.set_trace()
nts_milestones=len(ts_milestones)
psffname='prmtop/{0}.prmtop'.format(name)
top=md.load_prmtop(psffname)
ligand_com_data=[]
#need to read the COM data, and slice by interval
null_xyz=np.zeros((0,top.n_atoms,3))
bound_traj_all=md.Trajectory(null_xyz,top,unitcell_lengths=np.zeros((0,3)),unitcell_angles=np.zeros((0,3)))
ts_traj_all=md.Trajectory(null_xyz,top,unitcell_lengths=np.zeros((0,3)),unitcell_angles=np.zeros((0,3)))
for iwindow in range(1,nwindow+1):
	this_lig_com_data=[]
	if ((name=='pyk2-cpd1-amber') or (name=='pyk2-cpd10-amber')):
		fname='{0}/data-rel-to-first/{1}-tramd-com-{2:d}'.format(scratch,name,iwindow)
	else:
		fname='{0}/data-rel-to-first/{1}-com-{2:d}'.format(scratch,name,iwindow)
	print('reading ligand COM data from ',fname,flush=True)
	input=open(fname,'r')
	for line in input:
		words=line.split()
		#format: time x y z rmsd
		this_lig_com_data.append(np.array([float(words[1]),float(words[2]),float(words[3])]))
	input.close()
	this_lig_com_data=np.stack(this_lig_com_data)
	print(this_lig_com_data.shape)
	nframes=this_lig_com_data.shape[0]
	dist=np.full((nframes,nanchor),np.nan,dtype=np.float64)
	for kanchor in range(0,nanchor):
		diff=this_lig_com_data[:,:]-np.expand_dims(anchors[kanchor,:],axis=0)
		dist[:,kanchor]=np.sqrt(np.einsum('ij,ij->i',diff,diff)) #intended shape (nframes)

	#choose the two anchors that are closest for each frame -- +1 for conversion to 1-based 
	closest_anchor=np.argsort(dist,axis=1)+1
	closest_anchor=np.sort(closest_anchor[:,0:2],axis=1) #find the two closest anchors for each frame, ensure ianchor<janchor
	#identify those that are on the bound state milestone
	this_lig_com_bound=np.all([(closest_anchor[:,0]==bound_milestone[0]),(closest_anchor[:,1]==bound_milestone[1])],axis=0)
	this_lig_com_bound_count=np.count_nonzero(this_lig_com_bound)
	#search for frames that are on any ts milestone
	this_lig_com_ts=[]
	for m in ts_milestones:
		this_lig_com_ts.append(np.all([(closest_anchor[:,0]==m[0]),(closest_anchor[:,1]==m[1])],axis=0))
	this_lig_com_ts=np.any(np.stack(this_lig_com_ts),axis=0)
	this_lig_com_ts_count=np.count_nonzero(this_lig_com_ts)
	print('window {0:d} has {1:d} bound state and {2:d} transition state frames'.format(
		iwindow,this_lig_com_bound_count,this_lig_com_ts_count),flush=True)
	if ((this_lig_com_bound_count>0) or (this_lig_com_ts_count>0)):
		if ((name=='pyk2-cpd1-amber') or (name=='pyk2-cpd10-amber')):
			dcdfname='{0}/{1}-tramd-window-{2:d}-aligned.dcd'.format(dir,name,iwindow)
		else:
			dcdfname='{0}/{1}-window-{2:d}-aligned.dcd'.format(dir,name,iwindow)
		print('reading trajectory ',dcdfname,flush=True)
		traj=md.load(dcdfname,top=top)
		#pdb.set_trace()
		if (this_lig_com_bound_count>0):
			bound_traj=traj.slice(this_lig_com_bound,copy=True)
			bound_traj_all=bound_traj_all.join(bound_traj)
		if (this_lig_com_ts_count>0):
			ts_traj=traj.slice(this_lig_com_ts,copy=True)
			ts_traj_all=ts_traj_all.join(ts_traj)


print('total {0:d} bound state and {1:d} transition state frames'.format(bound_traj_all.n_frames,ts_traj_all.n_frames),flush=True)
outfname='{0}/transition-state/{1}-bound-state.dcd'.format(scratch,name)
print('saving bound state to ',outfname,flush=True)
bound_traj_all.save_dcd(outfname)
outfname='{0}/transition-state/{1}-transition-state.dcd'.format(scratch,name)
print('saving transition state to ',outfname,flush=True)
ts_traj_all.save_dcd(outfname)

