#!/usr/bin/env python3.6

#must use python 2.7.16
import sys
import pyemma
import pyemma.coordinates.transform
from pyemma._base.progress import ProgressReporterMixin, ProgressReporter
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import pickle
import mdtraj as md

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
nunbiased=int(sys.argv[2])
nwindow=int(sys.argv[3])
nsim=nunbiased+nwindow
#tica_lag_time=float(sys.argv[4])
#ndim=int(sys.argv[5])
short=False
scratch='/mnt/stor/ceph/scratch/jmsc87/pyk2/umbrella-amber3/'
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



tica_lag_times=[0.001,0.002,0.005,0.01,0.02,0.05,1.0,2.0,5.0,10.0]
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
for tica_lag_time in tica_lag_times:
	print('loading TICA object for lag time {0:.3f}'.format(tica_lag_time),flush=True)
	infname='{0}/tram-anal2b/{1}-tica-lag-time-{2:.3f}.dat'.format(scratch,name,tica_lag_time)
	tica=pyemma.coordinates.transform.TICA.load(infname,model_name='tica')
	ndim_tot=len(tica.cumvar)+1
	ndim_50=np.argmax(tica.cumvar > 0.5)+1
	ndim_75=np.argmax(tica.cumvar > 0.75) + 1
	ndim_90=np.argmax(tica.cumvar > 0.9)+1
	ndim_95=np.argmax(tica.cumvar > 0.95)+1
	print("Dimensions for 50/75/90/95% of variance:",ndim_50,ndim_75,ndim_90,ndim_95)
	x=range(1,len(tica.cumvar)+1)
	ax.plot(x,tica.cumvar,label='lag time {0:.3f} ns'.format(tica_lag_time))
	del tica

ax.set_xlabel('TICA dimension')
ax.set_ylabel('cumulative variance')
ax.set_ylim(bottom=0,top=1)
ax.set_xticks(range(0,ndim_tot,50))
ax.set_yticks(np.linspace(0,1,num=11))
ax.grid()
ax.legend(loc='lower right',fontsize=11)
outfname='tica-plots/{0}-tica-cumvar.png'.format(name)
fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
plt.close(fig)



sys.exit(-1)






#print(tica_output)
cm = matplotlib.cm.get_cmap('gist_rainbow')
do_density_plots=True
#do_time_plots=True
#do_autocorr_plots=False
if (do_density_plots):
	for idim in range(0,ndim-1):
		#print(x)
		#print(x.shape)
		#t=np.linspace((start+1)*dt,nframes*dt,num=(nframes-start))
		#print(t)
		#if (do_time_plots):
		#	print('time plot for dimension {0:d}'.format(idim+1),flush=True)
		#	fig=plt.figure()
		#	ax=fig.add_subplot(1,1,1)
		#	for iwindow in range(0,nwindow+1):
		#		if (iwindow==0):
		#		frac=float(iwindow)/float(nwindow)
		#		ax.plot(t,tica_output[iwindow][:,idim],c=cm(frac))
		#	ax.set_xlabel('time (ns)')
		#	ax.set_ylabel('dimension {0:d}'.format(idim+1))
		#	outfname='{0}-tica-{1:d}-time.png'.format(name,idim+1)
		#	fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
		#	plt.close(fig)
		#sys.exit(-1)
		for jdim in range(idim+1,idim+2):
			print('scatterplot for dimensions {0:d} and {1:d}'.format(idim+1,jdim+1),flush=True)
			fig=plt.figure()
			ax=fig.add_subplot(1,1,1)
			for iwindow in (list(range(1,nwindow+1))+[0]):
				if (iwindow==0):
					color='black'
				else:
					frac=float(iwindow)/float(nwindow)
					color=cm(frac)
				x=tica_output[iwindow][:,idim]
				#print(np.min(x),np.max(x))
				y=tica_output[iwindow][:,jdim]
				#print(np.min(y),np.max(y))
				#print(x.shape)
				#print(y.shape)
				ax.scatter(x,y,s=0.1,marker=',',edgecolors='none',color=color)
				#pyemma.plots.plot_density(x,y,ax,logscale=True)
			ax.set_xlabel('dimension {0:d}'.format(idim+1))
			ax.set_ylabel('dimension {0:d}'.format(jdim+1))
			outfname='tica-plots/{0}-tica-lag-time-{1:.3f}-{2:d}-{3:d}.png'.format(name,tica_lag_time,idim+1,jdim+1)
			fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
			plt.close(fig)


sys.exit(-1)
#fig, axes = plt.subplots(4, 1, figsize=(12, 5), sharex=True)
#x = 0.1 * np.arange(tica_output[0].shape[0])
#for i, (ax, tic) in enumerate(zip(axes.flat, tica_output[0].T)):
#    ax.plot(x, tic)
#    ax.set_ylabel('IC {}'.format(i + 1))
#axes[-1].set_xlabel('time / ns')
#fig.tight_layout()

print(np.array(tica_output).shape)
#n_clustercenters = [5, 10, 30, 75, 200]
n_clustercenters=[5,10,20,30,40,50,75,100,200,300]
scores = np.zeros((len(n_clustercenters), 5))
for n, k in enumerate(n_clustercenters):
    for m in range(5):
        #with pyemma.util.contexts.settings(show_progress_bars=False):
            print('estimating markov model with {0:d} clusters iteration {1:d}'.format(k,m))
            _cl = pyemma.coordinates.cluster_kmeans(
                tica_output, k=k, max_iter=50, stride=50)
            _msm = pyemma.msm.estimate_markov_model(_cl.dtrajs, 5)
            scores[n, m] = _msm.score_cv(
                _cl.dtrajs, n=1, score_method='VAMP2', score_k=min(10, k))

fig, ax = plt.subplots()
lower, upper = pyemma.util.statistics.confidence_interval(scores.T.tolist(), conf=0.9)
ax.fill_between(n_clustercenters, lower, upper, alpha=0.3)
ax.plot(n_clustercenters, np.mean(scores, axis=1), '-o')
ax.semilogx()
ax.set_xlabel('number of cluster centers')
ax.set_ylabel('VAMP-2 score')
outfname='tica-vamp-scores.png'
fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
plt.close(fig)
#fig.tight_layout()


