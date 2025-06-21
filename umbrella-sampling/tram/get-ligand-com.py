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
#nev_plot=int(sys.argv[3])
#short=((len(sys.argv)>2) and (sys.argv[3]=="short"))
short=False
scratch='/mnt/stor/ceph/scratch/jmsc87/pyk2/'

#sys.exit(-1)
#dcd='/mnt/stor/scratch/jmsc87/fak/anal/fak-cpd32-cgenff-explicit.dcd'
#tica_lag=int(0.05/dt)

print('setting up featurizer',flush=True)
psf='../anal3/{0}.prmtop'.format(name)
featurizer = pyemma.coordinates.featurizer(psf)
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
drug_res=271
featurizer.add_residue_COM([drug_res],scheme='all',mass_weighted=True)


#atoms=featurizer.select("residue 420")
#print(atoms)
#sys.exit(-1)
#featurizer.add_all()
#res_start=420
##res_end=690
#nres=res_end-res_start+1
#print(res_start,res_end,nres)


desc=featurizer.describe()
#features_res_list=[None]*len(desc)
#descfname='{0}-features'.format(name)
#descfile=open(descfname,'w')
ifeat_com=[]
for i, s in enumerate(desc):
	#descfile.write('{0:d} {1}\n'.format(i,s))
	if (s.find("COM")>=0):
		ifeat_com.append(i)
#descfile.close()
print("Total number of features",featurizer.dimension())
print(ifeat_com)
#sys.exit(-1)


fname='{0}/umbrella-amber3/tram-anal2b/{1}-ligand-com-data.npy'.format(scratch,name)
if (short):
	dir_unbiased=scratch+'unbiased-amber/dcd-short/'
	dir_umbrella=scratch+'umbrella-amber3/dcd-short/'
	start=0
	dt=0.01
else:
	dir_unbiased=scratch+'unbiased-amber/dcd-joined/'
	dir_umbrella=scratch+'umbrella-amber3/dcd-joined/'
	start=0
	dt=0.001
#idlist=list(range(1,6))+[11,22,46,64,73,87,98,121,138,143,164,347,409]
if (name=='pyk2-cpd1-amber'):
        idlist=list(range(1,6))+[11,22,46,64,73,87,98,121,138,143,164,347,409]
elif (name=='pyk2-cpd10-amber'):
        idlist=list(range(1,6))+[10,24,52,72,82,95,136,164,166,190,239,270,310,393,411,418]
else:
        print('need to specify unbiased trajectories for name {0}'.format(name))
        sys.exit(-1)
dcdlist_unbiased=['{0}/{1}-tramd-id-{2:d}-aligned.dcd'.format(dir_unbiased,name,id) for id in idlist]
dcdlist_umbrella=['{0}/{1}-tramd-window-{2:d}-aligned.dcd'.format(dir_umbrella,name,iwindow) for iwindow in range(1,nwindow+1)]
dcdlist=dcdlist_unbiased+dcdlist_umbrella
print(dcdlist)
features_data = pyemma.coordinates.source(dcdlist, features=featurizer, dtype=np.float64)
#the unbiased trajectories are first in the list
nframes_unbiased=0
nframes_umbrella=0
for itraj in range(0,nunbiased):
	nframes_unbiased+=features_data.trajectory_length(itraj,skip=start)
for itraj in range(nunbiased,nsim):
	nframes_umbrella+=features_data.trajectory_length(itraj,skip=start)
nframes_total=nframes_unbiased+nframes_umbrella
frac_unbiased=float(nframes_unbiased)/float(nframes_total)
print("Unbiased, umbrella frames: {0:d} {1:d}".format(nframes_unbiased,nframes_umbrella),flush=True)
print("Total frames: {0:d}".format(nframes_total),flush=True)
print("Fraction unbiased: {0:.3f}".format(frac_unbiased),flush=True)
#sys.exit(-1)
ligand_com=features_data.get_output(dimensions=ifeat_com,skip=start)
print(len(ligand_com),ligand_com[0].shape)
print(ligand_com[0])
print(ligand_com[nunbiased]) #first umbrella window
np.save(fname,ligand_com)
sys.exit(0)
