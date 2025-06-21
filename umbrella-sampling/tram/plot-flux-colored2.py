#!/usr/bin/env python3.6

#must use python 2.7.16
import sys
import pyemma
import pyemma.coordinates.transform
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
#matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import pickle
import mdtraj as md
from pyemma._base.progress import ProgressReporterMixin, ProgressReporter
from functools import reduce
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
#	#sys.exit(-1)
#	return result

name=sys.argv[1]
#nwindow=int(sys.argv[2])
tica_lag_time=float(sys.argv[2])
dcutoff=10.0
if (len(sys.argv)>3):
	dcutoff=float(sys.argv[3])
#ncluster=int(sys.argv[3])
#itrial=int(sys.argv[4])
#msm_lag_time=float(sys.argv[5])
#connectivity_factor=float(sys.argv[6])
#neigvec=int(sys.argv[7]) #number of eigenvectors
#n_plot_dim=int(sys.argv[8])
short=False
scratch='/mnt/stor/ceph/scratch/jmsc87/pyk2/'
if (short):
        dir=scratch+'dcd-short/'
        start=0
        dt=0.01
        #tica_lag_times=[0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0]
else:
        dir=scratch+'dcd-joined/'
        start=0
        dt=0.001
        #tica_lag_times=[0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0,2.0,5.0]
psf='../anal3/prmtop/{0}.prmtop'.format(name)



print('setting up featurizer',flush=True)
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
#featurizer.add_selection(featurizer.select_Ca())#featurizer.add_selection(drug_indices)
featurizer.add_backbone_torsions(cossin=True, periodic=False)
featurizer.add_sidechain_torsions(cossin=True,periodic=False)
#featurizer.add_custom_func(get_com,3,"resname CPD")
drug_res=271
featurizer.add_residue_COM([drug_res],scheme='all',mass_weighted=True)

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
#features_res_list=[]
desc=featurizer.describe()
#features_res_list=[None]*len(desc)
descfname='{0}-features'.format(name)
descfile=open(descfname,'w')
for i, s in enumerate(desc):
	descfile.write('{0:d} {1}\n'.format(i,s))
	#if ((s.find("PHI")>0) or (s.find("PSI")>0) or (s.find("CHI")>0)):
	#	#this feature came from an amino acid residue
	#	words=s.replace("("," ").replace(")"," ").split()
	#	res=int(words[-1])
	#	features_res_list[i]=res
descfile.close()
print(featurizer.dimension())
#print(features_res_list)
#res_start=min(mask(None,features_res_list))
#res_end=max(mask(None,features_res_list))
#nres=res_end-res_start+1
#print(res_start,res_end,nres)
#sys.exit(-1)
#print('loading trajectory',flush=True)
#dcdlist=['{0}/{1}-tramd-window-{2:d}-aligned.dcd'.format(dir,name,iwindow) for iwindow in range(1,nwindow+1)]
#print(dcdlist)
#sys.exit(-1)
#check to make sure the get_com function worked properly
#features_data = pyemma.coordinates.source(dcdlist, features=featurizer, dtype=np.float64)
#output=features_data.get_output()
#print(len(output),output[0].shape)
#print(output)
#can change to "pyemma.coordinates.source" to stream the data instead
#labels = ['backbone\ntorsions']
#sys.exit(-1)
#print('performing TICA',flush=True)
#print('loading TICA object with lag time {:f} ns'.format(tica_lag_time),flush=True)
#tica_lag=int(0.1/dt)
#tica = pyemma.coordinates.tica(features_data, lag=tica_lag, dtype=np.float64,skip=start)
#infname='{0}/umbrella-amber3/tram-anal2/{1}-tica-lag-time-{2:.3f}.dat'.format(scratch,name,tica_lag_time)
#input=open(infname,'rb')
#tica=pickle.load(input)
#input.close()
#tica=pyemma.coordinates.transform.TICA.load(infname,model_name='tica')
#tica.data_producer = features_data

#print('loading TICA projection with lag time {:f} ns'.format(tica_lag_time),flush=True)
#infname2=infname[0:-3]+'npy'
#tica_output=load_numpy_array_showing_progress(infname2,1,'loading TICA projection')
#tica_output=np.load(infname2,allow_pickle=True)
#n_tica_dim=np.argmax(tica.cumvar > 0.95)+1
#print("Dimension for 95% of variance:",n_tica_dim)
#print(tica_output.shape)
#print(tica_output[0].shape)
#we need to exclude the unneeded dims
#tica_output_reduced=[tica_output[i][:,0:n_tica_dim] for i in range(0,len(tica_output))]
#print(tica_output_reduced.shape)
#print(tica_output_reduced[0].shape)
#tica_output_restacked=np.concatenate(tica_output_reduced)
#print(tica_output_restacked.shape)
#assigned_free_energies=multi_therm_model.f[dtrajs_stacked]
#print(assigned_free_energies.shape)


print('loading ligand COM data',flush=True)
fname='{0}/umbrella-amber3/tram-anal2b/{1}-ligand-com-data.npy'.format(scratch,name)

ligand_com=np.load(fname,allow_pickle=True)
#convert all COM's from nm to A
for i in range(0,len(ligand_com)):
        ligand_com[i]=ligand_com[i]*10
print(ligand_com.shape,ligand_com[0]) #,ligand_com[nunbiased])
#sys.exit(-1)
N_k=[ligand_com[i].shape[0] for i in range(0,len(ligand_com))]
N_k=np.array(N_k)
print(N_k)
frames_by_state=[]
for i in range(0,len(N_k)):
        frames_by_state.append(np.full(N_k[i],i,dtype=np.int32))
frames_by_state=np.concatenate(frames_by_state)

print(frames_by_state)


ligand_com_restacked=np.concatenate(ligand_com)
print(ligand_com_restacked.shape)
nframes=ligand_com_restacked.shape[0]


mfpt_fig=plt.figure()
mfpt_ax=mfpt_fig.add_subplot(1,1,1)
if (tica_lag_time>=1.0):
        plotdir='plots'
else:
        plotdir='plots-short-tica-lag'
connectivity_factor=1.0
itrial=1
ncluster_list=[50,100,200]
msm_lag_time_list=[1.0,2.0,5.0,10.0]
mfpt_fname='{0}-mfpt-all'.format(name)
mfpt_file=open(mfpt_fname,'w')
for ncluster in ncluster_list:
	for msm_lag_time in msm_lag_time_list:
		print('loading msm model {0:d} clusters trial {1:d} MSM lag time {2:.3f} ns'.format(ncluster,itrial,msm_lag_time),flush=True)
		if (tica_lag_time>=1.0):
			objfname='{0}/umbrella-amber3/tram-anal2b/{1}-tram-objects-{2:d}-clusters.dat'.format(
				scratch,name,ncluster)
		else:
			objfname='{0}/umbrella-amber3/tram-anal2b/{1}-tram-objects-short-tica-lag-{2:d}-clusters.dat'.format(
				scratch,name,ncluster,msm_lag_time)

		try:
			objname='clusters-trial-{0:d}'.format(itrial)
			cluster=pyemma.load(objfname,objname)

			objname='tram-{0:d}-clusters-trial-{1:d}-msm-lag-{2:.3f}-cf-{3:.1f}'.format(
				ncluster,itrial,msm_lag_time,connectivity_factor)
			tram_model=pyemma.load(objfname,objname)
		except:
			print('error loading model file {0} name {1}'.format(objfname,objname))
			continue


		#It appears that the LAST model in the list represents the unbiased state,
		#see thermo/util/util.py line 171
		unbiased_model=tram_model.models[-1]
		#print(tram_model.models)
		#sys.exit(-1)
		#print('fraction of states used = {:f}'.format(unbiased_model.active_state_fraction))
		#print('fraction of counts used = {:f}'.format(unbiased_model.active_count_fraction))
		#unbiased_model.pcca(nstates)
		#print(len(cluster.clustercenters))
		#print(len(unbiased_model.metastable_assignments))
		#print(unbiased_model.metastable_assignments)
		print(unbiased_model.active_set)
		nactive=len(unbiased_model.active_set)
		active_state_fraction=float(len(unbiased_model.active_set))/float(ncluster)
		#print(unbiased_model.connectivity)
		eigvec=unbiased_model.eigenvectors_right()
		print(eigvec.shape)
		print('first eigenvector is one: {} (min={}, max={})'.format(
			np.allclose(eigvec[:, 0], 1, atol=1e-15), eigvec[:, 0].min(), eigvec[:, 0].max()))
		print(eigvec[:,1])
		#sys.exit(-1)

		#The "mapped_dtrajs" array contains microstate indices for each frame  mapped to the active set from the MSM.
		#It contains a -1 for each frame that does not belong to a state from the active set.
		#stable_dtrajs contains PCCA state assignments only for those  frames that have been mapped to the active set
		dtrajs_restacked=np.concatenate(cluster.dtrajs)
		#this boolean array shows for each frame whether it was in the active set
		mask=np.isin(dtrajs_restacked,unbiased_model.active_set)
		active_count_fraction=float(np.count_nonzero(mask))/float(len(dtrajs_restacked))
		#this is a list of the dtrajs for the frames that are in the active set
		active_dtrajs=dtrajs_restacked[mask]
		print(active_dtrajs)
		#this is the indices into the active set for them
		mapped_dtrajs=np.searchsorted(unbiased_model.active_set,active_dtrajs)
		print(mapped_dtrajs)
		#mapped_dtrajs=[unbiased_model.active_set.index(dtrajs_restacked[i]) for i in range(0,len(dtrajs_restacked))]
		print('fraction of states used = {:f}'.format(active_state_fraction))
		print('fraction of counts used = {:f}'.format(active_count_fraction))
		#sys.exit(-1)
		#mapped_dtrajs=np.where(mask,dtrajs_restacked,-1)
		#stable_dtrajs=unbiased_model.metastable_assignments[mapped_dtrajs[(mapped_dtrajs>=0)]]
		#print(stable_dtrajs.shape)
		#print(mapped_dtrajs.shape)
		#sys.exit(-1)
		#print(len(unbiased_model.discrete_trajectories_active),unbiased_model.discrete_trajectories_active[0].shape)
		#mask=np.where(mapped_dtrajs>=0)[0]
		#print(len(mask))
		#print(unbiased_model.discrete_trajectories_active)
		#sys.exit(-1)
		#stable_trajs=unbiased_model.metastable_assignments[cluster.dtrajs]

		#calculate the distance of each frame from the initial position 
		initial_center=ligand_com_restacked[0,:]
		print(initial_center)

		ligand_com_masked=ligand_com_restacked[mask,:]
		diff=ligand_com_masked[:,:]-np.expand_dims(initial_center[:],axis=0) #shape (nframes,3) t
		dist=np.sqrt(np.einsum('ij,ij->i',diff,diff)) #intended shape (nframes)
		print(dist.shape,dist)
		maxdist=np.max(dist)
		print(maxdist)
		#for each cluster in the active set, calculate the mean COM distance over all frames in that cluster
		dist_list=[]
		avgdist=np.full(nactive,np.nan)
		for iclus in range(0,nactive):
			indexlist=np.where(mapped_dtrajs==iclus)
			dist_list.append(dist[indexlist])
			avgdist[iclus]=np.mean(dist[indexlist])
		print(avgdist)
		distfrac=avgdist[:]/maxdist
		fig=plt.figure()
		ax=fig.add_subplot(1,1,1)
		ax.violinplot(dataset=dist_list,positions=range(0,nactive),widths=0.5,vert=False,showmeans=True,showextrema=False)
		ax.set_xlabel('ligand COM distance (A)')
		ax.set_ylabel('cluster number')
		plotdir='plots'
		outfname='{0}/{1}-lcom-dist-{2:d}-clusters-trial-{3:d}-msm-lag-{4:.3f}-cf-{5:.1f}.png'.format(
			plotdir,name,ncluster,itrial,msm_lag_time,connectivity_factor)
		fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
		plt.close(fig)



		#outfname='{0}-lcom-{1:d}-clusters-trial-{2:d}-msm-lag-{3:.3f}-cf-{4:.1f}'.format(name,ncluster,itrial,msm_lag_time,connectivity_factor)

		#output=open(outfname,'w')
		#for i in range(0,len(mapped_dtrajs)):
		#	line='{0:d} {1:f} {2:f} {3:f} {4:d}\n'.format(i,ligand_com_masked[i,0],ligand_com_masked[i,1],ligand_com_masked[i,2],
		#		mapped_dtrajs[i])
		#	output.write(line)
		#output.close()




		#Try to calculate MFPT !
		#cutoff=5.0
		cutoff_list=[]
		mfpt_list=[]
		dt_sec=1.0e-12 #one frmae is 1 ps
		cutoff=0.0
		while (cutoff<maxdist):
			A=np.where(avgdist<=cutoff)[0]
			B=np.where(avgdist>cutoff)[0]
			#print(cutoff,A,B)
			if ((len(A)>0) and (len(B)>0)):
				mfpt=unbiased_model.mfpt(A,B)*dt_sec
				print('MFPT: ',cutoff,mfpt)
				cutoff_list.append(cutoff)
				mfpt_list.append(mfpt)
				#line='{0:d} {1:.3f} {2:.1f} {3:f}\n'.format(ncluster,msm_lag_time,cutoff,mfpt)
				#no. clusters, msm lag time, off-rate = 1/mfpt
				if (cutoff==dcutoff):
					line='{0:d} & {1:.1f} & {2:.2f}\\\\\n'.format(ncluster,msm_lag_time,1.0/mfpt)
					mfpt_file.write(line)
			cutoff+=0.5
		label='{0:d} clusters, MSM lag time {1:.0f} ns'.format(ncluster,msm_lag_time)
		mfpt_ax.plot(cutoff_list,mfpt_list,label=label)
		del cluster
		del tram_model

mfpt_file.close()
mfpt_ax.set_yscale('log')
mfpt_ax.set_xlabel('distance cutoff (A)')
mfpt_ax.set_ylabel('MFPT (s)')
mfpt_ax.legend(loc='lower right',fontsize=8)

#print(cutoff_list)
#print(mfpt_list)
outfname='{0}/{1}-mfpt-all.png'.format(plotdir,name)
#-{2:d}-clusters-trial-{3:d}-msm-lag-{4:.3f}-cf-{5:.1f}.png'.format(
#        plotdir,name,ncluster,itrial,msm_lag_time,connectivity_factor)
mfpt_fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
plt.close(mfpt_fig)


sys.exit(-1)


#try to color the 
cm = matplotlib.cm.get_cmap('nipy_spectral') #'gist_rainbow')
state_colors=cm(distfrac) #will this work?
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
_,pos=pyemma.plots.plot_markov_model(unbiased_model.P,state_colors=state_colors,ax=ax)
if (tica_lag_time>=1.0):
        plotdir='plots'
else:
        plotdir='plots-short-tica-lag'
outfname='{0}/{1}-msm-flux-{2:d}-clusters-trial-{3:d}-msm-lag-{4:.3f}-cf-{5:.1f}.png'.format(
	plotdir,name,ncluster,itrial,msm_lag_time,connectivity_factor)
fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
plt.close(fig)





sys.exit(-1)

