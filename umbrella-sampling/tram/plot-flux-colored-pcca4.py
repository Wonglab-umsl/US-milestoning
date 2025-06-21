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
import random
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
nwindow=int(sys.argv[2])
tica_lag_time=float(sys.argv[3])
ncluster=int(sys.argv[4])
#itrial=int(sys.argv[4])
itrial=1
msm_lag_time=float(sys.argv[5])
#connectivity_factor=float(sys.argv[6])
connectivity_factor=1.0
npcca=int(sys.argv[6])
#nframes_pick=int(sys.argv[6])
seed=None
if (len(sys.argv)>7):
        seed=int(sys.argv[7])

#neigvec=int(sys.argv[7]) #number of eigenvectors
#n_plot_dim=int(sys.argv[8])
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
psf='../anal3/{0}.prmtop'.format(name)
#nwindow=421
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



print('setting up featurizer',flush=True)
featurizer = pyemma.coordinates.featurizer(psf)
#the below adds too many coordinates (a total of 9485)
#drug_indices=featurizer.select("resname CPD and (mass>=2.0)")
#print(drug_indices)
#featurizer.add_distances_ca()
#featurizer.add_distances(indices=featurizer.select_Ca(),indices2=drug_indices)
startres=1
endres=270
drugres=271
pairlist=[[ires,drugres] for ires in range(startres,endres+1)]
pairlist=np.array(pairlist)
print(pairlist)
print(pairlist.shape)
#sys.exit(-1)
featurizer.add_residue_mindist(residue_pairs=pairlist,ignore_nonprotein=False,periodic=False)
#featurizer.add_selection(featurizer.select_Ca())#featurizer.add_selection(drug_indices)
#featurizer.add_backbone_torsions(cossin=True, periodic=False)
#featurizer.add_sidechain_torsions(cossin=True,periodic=False)
#featurizer.add_custom_func(get_com,3,"resname CPD")
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
features_data = pyemma.coordinates.source(dcdlist, features=featurizer, dtype=np.float64)
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



#ncluster_list=[10,20,50,100,200,500,1000]
print('loading msm model {0:d} clusters trial {1:d}'.format(ncluster,itrial),flush=True)

if (tica_lag_time>=1.0):
	objfname='{0}/umbrella-amber3/tram-anal2b/{1}-tram-objects-{2:d}-clusters.dat'.format(
		scratch,name,ncluster)
else:
	objfname='{0}/umbrella-amber3/tram-anal2b/{1}-tram-objects-short-tica-lag-{2:d}-clusters.dat'.format(
		scratch,name,ncluster,msm_lag_time)


objname='clusters-trial-{0:d}'.format(itrial)
cluster=pyemma.load(objfname,objname)

objname='tram-{0:d}-clusters-trial-{1:d}-msm-lag-{2:.3f}-cf-{3:.1f}'.format(
	ncluster,itrial,msm_lag_time,connectivity_factor)
tram_model=pyemma.load(objfname,objname)




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

print(unbiased_model.eigenvalues())
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


#sys.exit(-1)

outfname='{0}-lcom-{1:d}-clusters-trial-{2:d}-msm-lag-{3:.3f}-cf-{4:.1f}'.format(name,ncluster,itrial,msm_lag_time,connectivity_factor)

output=open(outfname,'w')
for i in range(0,len(mapped_dtrajs)):
	line='{0:d} {1:f} {2:f} {3:f} {4:d}\n'.format(i,ligand_com_masked[i,0],ligand_com_masked[i,1],ligand_com_masked[i,2],
		mapped_dtrajs[i])
	output.write(line)
output.close()



#fig=plt.figure()
#ax=fig.add_subplot(projection='3d')
#cm = matplotlib.cm.get_cmap('nipy_spectral') #'gist_rainbow')
#for iclus in range(0,nactive):
#        indexlist=np.where(mapped_dtrajs==iclus)
#        ax.scatter(ligand_com_masked[indexlist,0],ligand_com_masked[indexlist,1],ligand_com_masked[indexlist,2],
#		s=1.0,marker=',',edgecolors='none') #,color=cm(frac))
#if (tica_lag_time>=1.0):
#        plotdir='plots'
#else:
#        plotdir='plots-short-tica-lag'
#outfname='{0}/{1}-lcom-{2:d}-clusters-trial-{3:d}-msm-lag-{4:.3f}-cf-{5:.1f}.png'.format(
#        plotdir,name,ncluster,itrial,msm_lag_time,connectivity_factor)
#fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
#plt.close(fig)




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
	cutoff+=0.5

#print(cutoff_list)
#print(mfpt_list)
#fig=plt.figure()
#ax=fig.add_subplot(1,1,1)
#if (tica_lag_time>=1.0):
#        plotdir='plots'
#else:
#        plotdir='plots-short-tica-lag'
#ax.plot(cutoff_list,mfpt_list)
#ax.set_yscale('log')
#ax.set_xlabel('distance cutoff (A)')
#ax.set_ylabel('MFPT (s)')
#outfname='{0}/{1}-mfpt-{2:d}-clusters-trial-{3:d}-msm-lag-{4:.3f}-cf-{5:.1f}.png'.format(
#        plotdir,name,ncluster,itrial,msm_lag_time,connectivity_factor)
#fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
#plt.close(fig)


pcca_model = unbiased_model.pcca(npcca)
print(pcca_model.metastable_assignment)
print(pcca_model.metastable_assignment.shape)
print(pcca_model.metastable_sets)
print(unbiased_model.transition_matrix)
print(unbiased_model.pi) #stationary distribution
print(np.sum(unbiased_model.pi))
minchi=np.min(pcca_model.memberships)
print("Minchi: ",minchi)



#sys.exit(-1)
#print(pcca_model)
#cg_trans=np.copy(pcca_model.coarse_grained_transition_matrix)
#for some reason it doesn't satisfly requirements to be a transition matrix
#instead of trying to "fix" cg_trans from above, recalculate the CG transition matrix using the "crisp" membership assignments.
#it should satisfy the stochastic property directly 
#the correct method is equation 6 from my PCCA notes
cg_trans=np.zeros((npcca,npcca))
cg_stationary=np.zeros(npcca)
cg_frame_count=np.zeros(npcca,dtype=np.int64)
for i in range(0,nactive):
	icg=pcca_model.metastable_assignment[i]
	cg_stationary[icg]+=unbiased_model.pi[i]
	for j in range(0,nactive):
		jcg=pcca_model.metastable_assignment[j]
		cg_trans[icg,jcg]+=unbiased_model.transition_matrix[i,j]*unbiased_model.pi[i]
for icg in range(0,npcca):
	cg_trans[icg,:]=cg_trans[icg,:]/cg_stationary[icg]
	s=np.sum(cg_trans[icg,:])
	print(icg,s)
	#cg_trans[i,:]=cg_trans[i,:]/s
print(cg_trans)
print(cg_trans.shape)
print(cg_stationary)
#sys.exit(-1)

#output transition matrix in latex format
outfname='{0}-trans-{1:d}-clusters-{2:d}-pcca-{3:.3f}-msm-lag-time'.format(name,ncluster,npcca,msm_lag_time)
output=open(outfname,'w')
line=''
for jcg in range(0,npcca):
	line=line+'& {0:d} '.format(jcg)
line=line+'\\\\\n'
output.write(line)
output.write('\\hline\n')
for icg in range(0,npcca):
	line='{0:d}'.format(icg)
	for jcg in range(0,npcca):
		s='{0:.1E}'.format(cg_trans[icg,jcg])
		ss='$'+s.replace('E',' \\times 10^{')+'}$'
		line=line+' & ' +ss
	line=line+'\\\\\n'
	output.write(line)
output.close()
#print(pcca_model.memberships)
#print(pcca_model.n_metastable)
#sys.exit(-1)

#try to color the 
cm = matplotlib.cm.get_cmap('nipy_spectral') #'gist_rainbow')
state_colors=cm(distfrac) #will this work?
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
#labels=[str(unbiased_model.active_set[i]) for i in range(0,len(unbiased_model.active_set))]
_,pos=pyemma.plots.plot_markov_model(cg_trans,ax=ax,minflux=1.0e-12) #,state_labels=labels)
#flux=unbiased_model.
#_,pos=pyemma.plots.plot_flux
if (tica_lag_time>=1.0):
        plotdir='plots'
else:
        plotdir='plots-short-tica-lag'
outfname='{0}/{1}-msm-cgflux-{2:d}-clusters-{3:d}-pcca-msm-lag-{4:.3f}.png'.format(
	plotdir,name,ncluster,npcca,msm_lag_time,connectivity_factor)
fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
plt.close(fig)

#create a mapping from the original clusters, through the active set (sampled in unbiased ensemble) 
orig_cluster_to_pcca=[-1]*ncluster
for iactive in range(0,nactive):
	iclus=unbiased_model.active_set[iactive]
	ipcca=pcca_model.metastable_assignment[iactive]
	orig_cluster_to_pcca[iclus]=ipcca
orig_cluster_to_pcca=np.array(orig_cluster_to_pcca)
print(orig_cluster_to_pcca)
#sys.exit(-1)
#reform dtrajs into a list of frames by cluster
frames_by_cluster=[[] for iclus in range(0,npcca)]
for itraj,traj in enumerate(cluster.dtrajs):
	mapped_traj=orig_cluster_to_pcca[traj]
	#print(itraj,mapped_traj)
	for iclus in range(0,npcca):
		frames=np.where(mapped_traj==iclus)[0]
		#print(iclus,frames)
		#sys.exit(-1)
		if (len(frames)>0):
			#print(frames.shape)
			frames=np.vstack([np.full(frames.shape,itraj),frames]).transpose()
			frames_by_cluster[iclus].append(frames)
			#print(itraj,iclus,frames.shape,frames)

#write out data re stationary distribution,

outfname='{0}-fe-{1:d}-clusters-{2:d}-pcca-{3:.3f}-msm-lag-time'.format(name,ncluster,npcca,msm_lag_time)
output=open(outfname,'w')
maxprob=np.max(cg_stationary)
for icg in range(0,npcca):
        frames_by_cluster[icg]=np.vstack(frames_by_cluster[icg])
        fe=-1.987e-3*300*np.log(cg_stationary[icg]/maxprob)
        line='{0:d} & {1:d} & {2:.1E} & {3:.1f} \\\\\n'.format(icg,frames_by_cluster[icg].shape[0],
		cg_stationary[icg],fe)
        output.write(line)
output.close()


#print(frames_by_cluster)
#sys.exit(-1)



#write the picked frames from each cluster out to a dcd file
for iclus in range(0,npcca):
	indices=frames_by_cluster[iclus]
	outfname='{0}/umbrella-amber3/tram-anal2b/clusters/{1}-{2:d}-cg-clusters-msm-lag-time-{3:.3f}-cluster-{4:d}-of-{5:d}-all.dcd'.format(
		scratch,name,ncluster,msm_lag_time,iclus,npcca)
	print('writing cluster {0:d} to file {1} {2:d} frames'.format(iclus,outfname,indices.shape[0]))
	pyemma.coordinates.save_traj(traj_inp=features_data,indexes=indices,outfile=outfname)




sys.exit(-1)

