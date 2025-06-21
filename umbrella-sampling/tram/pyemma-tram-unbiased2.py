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
msm_lag_time=float(sys.argv[7])
connectivity_factor=float(sys.argv[8])
#nev_plot=int(sys.argv[3])
#short=((len(sys.argv)>2) and (sys.argv[3]=="short"))
short=False
scratch='/mnt/stor/ceph/scratch/jmsc87/pyk2/'
if (short):
	start=0
	dt=0.01
        #tica_lag_times=[0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0]
else:
	start=0
	dt=0.001
        #tica_lag_times=[0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0]

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


#desc=featurizer.describe()
##features_res_list=[None]*len(desc)
#descfname='{0}-features'.format(name)
#descfile=open(descfname,'w')
#ifeat_com=[]
#for i, s in enumerate(desc):
#	descfile.write('{0:d} {1}\n'.format(i,s))
#	if (s.find("COM")>=0):
#		ifeat_com.append(i)
#descfile.close()
#print("Total number of features",featurizer.dimension())
#print(ifeat_com)
#sys.exit(-1)


#fname='{0}-ligand-com-data.npy'.format(name)
fname='{0}/umbrella-amber3/tram-anal2b/{1}-ligand-com-data.npy'.format(scratch,name)
print('loading ligand COM data at ',datetime.now(),flush=True)
#this has been converted to float32 -- need to convert it back (better would be to write it as float64,
#but not sure how to control the type of get_output above.
#ligand_com=np.load(fname).astype(np.float64)
#ligand_com=load_numpy_array_showing_progress(fname,1,'loading ligand COM data')
ligand_com=np.load(fname,allow_pickle=True)
ligand_com=[ligand_com[i].astype(np.float64) for i in range(0,len(ligand_com))]
#ligand_com=ligand_com.astype(np.float64)
#ligand_com=list(ligand_com)
print(len(ligand_com),ligand_com[0].shape,ligand_com[0].dtype)
        #print(ligand_com.shape)
ligand_com_unbiased=ligand_com[0:nunbiased]
ligand_com_umbrella=ligand_com[nunbiased:nsim]
#sys.exit(-1)

umbrella_centers=np.full((nwindow,3),np.nan,dtype=np.float64)
fname='{0}-tramd-windows-rel-to-first'.format(name)
input=open(fname,'r')
for line in input:
	words=line.split()
	i=int(words[0])-1
	for k in range(0,3):
		#pyemma works in nm, my centers are in A
		umbrella_centers[i,k]=0.1*float(words[k+1])
input.close()
umbrella_centers=list(umbrella_centers)



print('loading TICA object for TICA lag time {0:.3f} ns at {1}'.format(tica_lag_time,datetime.now()),flush=True)
infname='{0}/umbrella-amber3/tram-anal2b/{1}-tica-lag-time-{2:.3f}.dat'.format(scratch,name,tica_lag_time)
tica=pyemma.coordinates.transform.TICA.load(infname,model_name='tica')

print('loading TICA projection at',datetime.now(),flush=True)
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

print('loading clustering with {0:d} clusters trial {1:d} at {2}'.format(ncluster,itrial,datetime.now()),flush=True)
#print('loading msm model {0:d} clusters trial {1:d}'.format(ncluster,itrial),flush=True)
if (tica_lag_time>=1.0):
        objfname='{0}/umbrella-amber3/tram-anal2b//{1}-tram-objects-{2:d}-clusters.dat'.format(scratch,name,ncluster)
else:
        objfname='{0}/umbrella-amber3/tram-anal2b//{1}-tram-objects-short-tica-lag-{2:d}-clusters.dat'.format(scratch,name,ncluster)
objname='clusters-trial-{0:d}'.format(itrial)
cluster=pyemma.load(objfname,objname)
dtrajs=cluster.dtrajs
print(len(dtrajs))
print(dtrajs[0].shape)

dtrajs_unbiased=dtrajs[0:nunbiased]
dtrajs_umbrella=dtrajs[nunbiased:nsim]
#the weight is exp[-beta*k*(r-r0)/2 - ln(z)]
#output will be in kt units
kt=(1.987e-3*300)
beta=1.0/kt
#factor of 100 is to convert from (kcal/mol/A^2) to (kcal/mol/nm^2) !
k=2.0*100*beta
force_constants=[k]*nwindow
print(force_constants)
#There is an undocumented parameter "width" that specifies periodicity of umbrella coordinates.  Zero=no periodicity.
#width=[0.0]*3
#width=[None]*3
#width=np.zeros(3,dtype=np.float64)
#width=np.require(width, requirements='C')
#print(width)
msm_lag=int(msm_lag_time/dt)
print('beginning tram estimation for lag time {0:.3f} ns {1:d} frames at {2}'.format(msm_lag_time,msm_lag,datetime.now()),flush=True)
tram_model=pyemma.thermo.estimate_umbrella_sampling(us_trajs=ligand_com_umbrella,us_dtrajs=dtrajs_umbrella,
        us_centers=umbrella_centers,us_force_constants=force_constants,
	md_trajs=ligand_com_unbiased,md_dtrajs=dtrajs_unbiased,
	estimator='tram',lag=msm_lag,connectivity_factor=connectivity_factor,
	init='mbar',init_maxerr=1.5e-2,maxerr=1.0e-2,save_convergence_info=1,
	direct_space=True)
#if (tica_lag_time>=1.0):
#       objfname='{0}/umbrella-amber3/tram-anal2/{1}-tram-objects-{2:d}-clusters-msm-lag-{3:.3f}.dat'.format(
#		scratch,name,ncluster,msm_lag_time)
#else:
#        objfname='{0}/umbrella-amber3/tram-anal2/{1}-tram-objects-short-tica-lag-{2:d}-clusters-msm-lag-{3:.3f}.dat'.format(
#		scratch,name,ncluster,msm_lag_time)
objname='tram-{0:d}-clusters-trial-{1:d}-msm-lag-{2:.3f}-cf-{3:.1f}'.format(ncluster,itrial,msm_lag_time,connectivity_factor)
print('tram estimation complete at',datetime.now(),flush=True)
tram_model.save(objfname,objname,overwrite=True)
print(tram_model.f)
sys.exit(-1)

