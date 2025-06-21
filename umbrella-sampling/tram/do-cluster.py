#!/usr/bin/env python3.6

#must use python 2.7.16
import sys
import pyemma
import pyemma.coordinates.transform
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import pickle
import mdtraj as md
from operator import itemgetter
from datetime import datetime

def load_numpy_array_showing_progress(fname,blocksize,desc):
	try:
		mmap = np.load(fname, mmap_mode='r')
		output = np.empty_like(mmap)
		#print(mmap.shape)
		n_blocks = int(np.ceil(mmap.shape[0] / blocksize))
		#n_blocks = int(np.ceil(mmap.size / blocksize))
		#print(n_blocks)
		pg = ProgressReporter()
		assert pg.show_progress
		pg.register(n_blocks,description=desc)
		for b in range(n_blocks):
			#print('progress: {}/{}'.format(b, n_blocks))  # use any progress indicator
			pg.update(1)
			output[b*blocksize : (b+1) * blocksize] = np.copy(mmap[b*blocksize : (b+1) * blocksize])
		pg.finish()
	finally:
		del mmap  # make sure file is closed again
		del pg
	return output



#to be passed to MDFeaturizer.add_custom_func.  First parameter is an MDTrajectory object.  Must return
#a numpy array of shape (nframes,3)
#check alignment of original trajectories
#def get_com(traj,selstr):
#	xyz=traj.xyz
#	top=traj.topology
#	indices=top.select(selstr)
#	totalmass=0.0
#	result=np.zeros((traj.n_frames,3),dtype=np.float32)
#	for iatom in indices:
#		mass=top.atom(iatom).element.mass
#		totalmass+=mass
#		result+=mass*xyz[:,iatom,:]
#	result=result/totalmass
#	#print(result)
	#sys.exit(-1)
#	return result

name=sys.argv[1]
nunbiased=int(sys.argv[2])
nwindow=int(sys.argv[3])
nsim=nunbiased+nwindow
tica_lag_time=float(sys.argv[4])
ncluster=int(sys.argv[5])
itrial=int(sys.argv[6])
#nev_plot=int(sys.argv[3])
#short=((len(sys.argv)>2) and (sys.argv[3]=="short"))
short=False
scratch='/mnt/stor/ceph/scratch/jmsc87/pyk2/'
#sys.exit(-1)
#dcd='/mnt/stor/scratch/jmsc87/fak/anal/fak-cpd32-cgenff-explicit.dcd'
#tica_lag=int(0.05/dt)

#print('setting up featurizer',flush=True)
#psf='../anal/{0}.prmtop'.format(name)
#featurizer = pyemma.coordinates.featurizer(psf)
#the below adds too many coordinates (a total of 9485)
#drug_indices=featurizer.select("resname CPD and (mass>=2.0)")
#print(drug_indices)
#featurizer.add_distances_ca()
#featurizer.add_distances(indices=featurizer.select_Ca(),indices2=drug_indices)
#startres=1
#endres=270
#drugres=271
#pairlist=[[ires,drugres] for ires in range(startres,endres+1)]
#pairlist=np.array(pairlist)
#print(pairlist)
#print(pairlist.shape)
#sys.exit(-1)
#featurizer.add_residue_mindist(residue_pairs=pairlist,ignore_nonprotein=False,periodic=False)
#featurizer.add_selection(featurizer.select_Ca())
#featurizer.add_selection(drug_indices)
#featurizer.add_backbone_torsions(cossin=True, periodic=False)
#featurizer.add_sidechain_torsions(cossin=True,periodic=False)
#featurizer.add_custom_func(get_com,3,"resname CPD")
#drug_atoms=featurizer.select("resname CPD")
#drug_res=271
#featurizer.add_residue_COM([drug_res],scheme='all',mass_weighted=True)

#dih_input=open('{0}-dihedrals'.format(name),'r')
#all_dih_indices=[]
#for line in dih_input:
#	dih_names=line.split()
#	if (len(dih_names)>=4):
#		dih_indices=[]
#		for atname in dih_names[0:4]:
#			idx=featurizer.select("resname CPD and name {0}".format(atname))
#			#print(atname,idx)
#			dih_indices.append(idx)
#		dih_indices=np.array(dih_indices)
#		#print(dih_indices)
#	all_dih_indices.append(dih_indices)
#all_dih_indices=np.array(all_dih_indices,dtype=np.int)
#all_dih_indices=np.hstack(all_dih_indices).transpose()
#print(all_dih_indices.shape)
#featurizer.add_dihedrals(indexes=all_dih_indices,cossin=True,periodic=False)

#atoms=featurizer.select("residue 420")
#print(atoms)
#sys.exit(-1)
#featurizer.add_all()
#res_start=420
#res_end=690
#nres=res_end-res_start+1
#print(res_start,res_end,nres)




print('loading TICA object for TICA lag time {0:.3f} ns'.format(tica_lag_time),flush=True)
infname='{0}/umbrella-amber3/tram-anal2b/{1}-tica-lag-time-{2:.3f}.dat'.format(scratch,name,tica_lag_time)
tica=pyemma.coordinates.transform.TICA.load(infname,model_name='tica')

print('loading TICA projection',flush=True)
infname2=infname[0:-3]+'npy'
#tica_output=load_numpy_array_showing_progress(infname2,1,'loading TICA projection')
tica_output=np.load(infname2,allow_pickle=True)
n_tica_dim=np.argmax(tica.cumvar > 0.95)+1
print("Dimension for 95% of variance:",n_tica_dim)
print(tica_output.shape)
print(tica_output[0].shape)
#we need to exclude the unneeded dims
tica_output_reduced=[tica_output[i][:,0:n_tica_dim] for i in range(0,len(tica_output))]
#print(tica_output_reduced.shape)
print(tica_output_reduced[0].shape)
tica_output_restacked=np.concatenate(tica_output_reduced)
print(tica_output_restacked.shape)
#sys.exit(-1)

n_plot_dim=5

print('beginning clustering with {0:d} clusters trial {1:d} at {2}'.format(ncluster,itrial,datetime.now()),flush=True)
cluster = pyemma.coordinates.cluster_kmeans(tica_output_reduced, k=ncluster, max_iter=500, stride=50)
cluster.overwrite_dtrajs=True
#cluster.save_dtrajs(prefix='dtrajs')
dtrajs=cluster.dtrajs

objfname='{0}/umbrella-amber3/tram-anal2b/{1}-tram-objects-{2:d}-clusters.dat'.format(scratch,name,ncluster)
objname='clusters-trial-{0:d}'.format(itrial)
cluster.save(objfname,objname,overwrite=True)

print(len(dtrajs))
print(dtrajs[0].shape)
#for (i,dtraj) in enumerate(dtrajs):
#	print(i,dtraj.shape)
unique, counts=np.unique(np.concatenate(dtrajs[0:nunbiased]),return_counts=True)
d=dict(zip(unique,counts))
l=sorted(d.items(), key = lambda kv:(kv[1], kv[0]),reverse=True)
print('unbiased',l)
#sys.exit(0)
for iwindow in range(1,nwindow+1):
	unique, counts=np.unique(dtrajs[nunbiased+iwindow-1],return_counts=True)
	d=dict(zip(unique,counts))
	l=sorted(d.items(), key = lambda kv:(kv[1], kv[0]),reverse=True)
	print(iwindow,l)

fname='{0}/umbrella-amber3/tram-anal2b/{1}-ligand-com-data.npy'.format(scratch,name)

print('loading ligand COM data',flush=True)
ligand_com=np.load(fname,allow_pickle=True)
#convert all COM's from nm to A
for i in range(0,len(ligand_com)):
        ligand_com[i]=ligand_com[i]*10
print(ligand_com.shape,ligand_com[0],ligand_com[nunbiased])
ligand_com_restacked=np.concatenate(ligand_com)
print(ligand_com_restacked.shape)

print(ligand_com_restacked[0,0],ligand_com_restacked[0,1],ligand_com_restacked[0,2])
#sys.exit(-1)

outfname='data/{0}-lcom-{1:d}-clusters-trial-{2:d}'.format(name,ncluster,itrial)

dtrajs_restacked=np.concatenate(dtrajs)
output=open(outfname,'w')
for i in range(0,len(dtrajs_restacked)):
	line='{0:d} {1:.6f} {2:.6f} {3:.6f} {4:d}\n'.format(i,
		ligand_com_restacked[i,0],ligand_com_restacked[i,1],ligand_com_restacked[i,2],dtrajs_restacked[i])
	output.write(line)
output.close()

sys.exit(0)
