#!/usr/bin/env python3.6

#must use python 2.7.16
import sys
import pyemma
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import pickle
import mdtraj as md
from operator import itemgetter

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
#short=((len(sys.argv)>2) and (sys.argv[3]=="short"))
short=False
scratch='/mnt/stor/ceph/scratch/jmsc87/pyk2/'
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
psf='../anal3/prmtop/{0}.prmtop'.format(name)
if (name=='pyk2-cpd1-amber'):
	idlist=list(range(1,6))+[11,22,46,64,73,87,98,121,138,143,164,347,409]
elif (name=='pyk2-cpd10-amber'):
	idlist=list(range(1,6))+[10,24,52,72,82,95,136,164,166,190,239,270,310,393,411,418]
elif (name=='pyk2-cpd2-amber'):
	idlist=list(range(1,6))+[66,72,79,84,93,100,102,115,118,128,137,142,143,149]
elif (name=='pyk2-cpd6-amber'):
	idlist=list(range(1,6))+[68,78,83,86,99,101,108,114,119,141,143,147,148,152,158,162]
elif (name=='pyk2-cpd8-amber'):
	idlist=list(range(1,6))+[63,68,74,77,85,96,98,104,110,113,121,131,143,147,156,158,168,176,178,184,187]
else:
	print('need to specify unbiased trajectories for name {0}'.format(name))
	sys.exit(-1)

if ((name=='pyk2-cpd1-amber') or (name=='pyk2-cpd10-amber')):
	dcdlist_unbiased=['{0}/{1}-tramd-id-{2:d}-aligned.dcd'.format(dir_unbiased,name,id) for id in idlist]
	dcdlist_umbrella=['{0}/{1}-tramd-window-{2:d}-aligned.dcd'.format(dir_umbrella,name,iwindow) for iwindow in range(1,nwindow+1)]
else:
	dcdlist_unbiased=['{0}/{1}-id-{2:d}-aligned.dcd'.format(dir_unbiased,name,id) for id in idlist]
	dcdlist_umbrella=['{0}/{1}-window-{2:d}-aligned.dcd'.format(dir_umbrella,name,iwindow) for iwindow in range(1,nwindow+1)]

dcdlist=dcdlist_unbiased+dcdlist_umbrella
print(dcdlist)
#sys.exit(-1)
#dcd='/mnt/stor/scratch/jmsc87/fak/anal/fak-cpd32-cgenff-explicit.dcd'
#tica_lag=int(0.05/dt)

print('setting up featurizer',flush=True)
featurizer = pyemma.coordinates.featurizer(psf)
#the below adds too many coordinates (a total of 9485)
#drug_indices=featurizer.select("resname CPD and (mass>=2.0)")
#print(drug_indices)
#featurizer.add_distances_ca()
#featurizer.add_distances(indices=featurizer.select_Ca(),indices2=drug_indices)
startres=1
endres=270
drug_res=271
pairlist=[[ires,drug_res] for ires in range(startres,endres+1)]
pairlist=np.array(pairlist)
print(pairlist)
print(pairlist.shape)
#sys.exit(-1)
featurizer.add_residue_mindist(residue_pairs=pairlist,ignore_nonprotein=False,periodic=False)
#featurizer.add_selection(featurizer.select_Ca())
#featurizer.add_selection(drug_indices)
#featurizer.add_backbone_torsions(cossin=True, periodic=False)
#featurizer.add_sidechain_torsions(cossin=True,periodic=False)
#featurizer.add_custom_func(get_com,3,"resname CPD")
#drug_atoms=featurizer.select("resname CPD")
#drug_res=271
#featurizer.add_residue_COM([drug_res],scheme='all',mass_weighted=True)

dih_input=open('{0}-dihedrals'.format(name),'r')
all_dih_indices=[]
for line in dih_input:
	dih_names=line.split()
	if (len(dih_names)>=4):
		dih_indices=[]
		for atname in dih_names[0:4]:
			idx=featurizer.select("resname CPD and name {0}".format(atname))
			#print(atname,idx)
			dih_indices.append(idx)
		dih_indices=np.array(dih_indices)
		#print(dih_indices)
	all_dih_indices.append(dih_indices)
#all_dih_indices=np.array(all_dih_indices,dtype=np.int)
all_dih_indices=np.hstack(all_dih_indices).transpose()
print(all_dih_indices.shape)
featurizer.add_dihedrals(indexes=all_dih_indices,cossin=True,periodic=False)

#atoms=featurizer.select("residue 420")
#print(atoms)
#sys.exit(-1)
#featurizer.add_all()
res_start=420
res_end=690
nres=res_end-res_start+1
print(res_start,res_end,nres)


features_res_list=[]
desc=featurizer.describe()
features_res_list=[None]*len(desc)
descfname='{0}-features'.format(name)
descfile=open(descfname,'w')
for i, s in enumerate(desc):
	descfile.write('{0:d} {1}\n'.format(i,s))
	if ((s.find("PHI")>=0) or (s.find("PSI")>=0) or (s.find("CHI")>=0)):
		#this feature came from an amino acid residue
		words=s.replace("("," ").replace(")"," ").split()
		res=int(words[-1])
		#amber numbers residues starting from 1
		features_res_list[i]=(res-1)+res_start
	#TODO: deal with the drug com and dihedral cases 
	elif (s.find("COM")>=0):
		#it's a ligand COM
		features_res_list[i]=res_end+1
	elif ((s.find("CPD")>=0) and (s.find("COM")<0)):
		features_res_list[i]=res_end+2
	else:
		pass
descfile.close()
print(featurizer.dimension())
#print(features_res_list)
#sys.exit(-1)
#res_start=min(filter(None,features_res_list))
#res_end=max(filter(None,features_res_list))
#sys.exit(-1)
print('loading trajectory',flush=True)
#sys.exit(-1)
#check to make sure the get_com function worked properly
features_data = pyemma.coordinates.source(dcdlist, features=featurizer, dtype=np.float64)
#the unbiased trajectories are first in the list
nframes_unbiased=0
nframes_umbrella=0
for itraj in range(0,nunbiased):
	nframes_unbiased+=features_data.trajectory_length(itraj,skip=start)
for itraj in range(nunbiased,nsim):
	nframes_umbrella+=features_data.trajectory_length(itraj,skip=start)
#for some reason they came out as tuples
nframes_total=nframes_unbiased+nframes_umbrella
frac_unbiased=float(nframes_unbiased)/float(nframes_total)
print("Unbiased, umbrella frames: {0:d} {1:d}".format(nframes_unbiased,nframes_umbrella),flush=True)
print("Total frames: {0:d}".format(nframes_total),flush=True)
print("Fraction unbiased: {0:.3f}".format(frac_unbiased),flush=True)
#sys.exit(-1)
#output=features_data.get_output()
#print(len(output),output[0].shape)
#print(output)
#can change to "pyemma.coordinates.source" to stream the data instead
#labels = ['backbone\ntorsions']
#sys.exit(-1)
#fig=plt.figure()
#ax=fig.add_subplot(1,1,1)
ndim_plot=20
tica_lag=int(tica_lag_time/dt)
print('performing TICA for lag time {0:.3f} ns {1:d} frames\n'.format(tica_lag_time,tica_lag),flush=True)
tica = pyemma.coordinates.tica(features_data, lag=tica_lag, dtype=np.float64,skip=start,var_cutoff=1.0)
print("Timescales: ",tica.timescales[0:ndim_plot])
#tica_output = np.array(tica.get_output())
#print("Dimension: ",tica_output.shape)
#nframes=tica.n_frames_total()//nwindow
#x=range(1,len(tica.cumvar)+1)
#ax.plot(x,tica.cumvar,label='lag time {0:.3f} ns'.format(tica_lag_time))
ndim=np.argmax(tica.cumvar > 0.95)+1
print("Dimension for 95% of variance:",ndim)
#dump the tica to a pickle file
if (short):
	outfname='{0}/umbrella-amber3/tram-anal2b/{1}-short-tica-lag-time-{2:.3f}.dat'.format(scratch,name,tica_lag_time)
else:
	outfname='{0}/umbrella-amber3/tram-anal2b/{1}-tica-lag-time-{2:.3f}.dat'.format(scratch,name,tica_lag_time)
#output=open(outfname,'wb')
tica.save(outfname,model_name='tica',overwrite=True)
#pickle.dump(tica,output)
#output.close()
#print(nframes)
outfname2=outfname[0:-3]+'npy'
np.save(outfname2,tica.get_output())
#ax.set_xlabel('TICA dimension')
#ax.set_ylabel('cumulative variance')
#ax.set_ylim(bottom=0.8,top=1)
#ax.grid()
#ax.legend(loc='lower right',fontsize=11)
#outfname='{0}-tica-cumvar.png'.format(name)
#fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
#plt.close(fig)
sys.exit(0)
