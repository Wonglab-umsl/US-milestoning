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
from pyemma._base.progress import ProgressReporterMixin, ProgressReporter
from datetime import datetime

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

#The class tram_its contains a single field called dict, which is a dictionary
#associating timescales to MEMM objects.
#class tram_its:
#	def __init__(self):
#		self.dict={}
#
#	def add_model(self,timescale,model):
#		self.dict[timescale]=model

	

name=sys.argv[1]
#nwindow=int(sys.argv[2])
tica_lag_time=float(sys.argv[2])
ncluster=int(sys.argv[3])
itrial=1
nits=int(sys.argv[4])
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
psf='../anal3/{0}.prmtop'.format(name)







#ncluster_list=[50,100]
msm_lag_time_list=[1.0,2.0,5.0]
#cf_list=[0.1,1.0,10.0]
itrial=1
connectivity_factor=1.0
if (tica_lag_time>=1.0):
        objfname='{0}/umbrella-amber3/tram-anal2b/{1}-tram-objects-{2:d}-clusters.dat'.format(
                scratch,name,ncluster)
else:
        objfname='{0}/umbrella-amber3/tram-anal2b/{1}-tram-objects-short-tica-lag-{2:d}-clusters.dat'.format(
                scratch,name,ncluster,msm_lag_time)
tram_models=[]
for msm_lag_time in msm_lag_time_list:
	print('loading msm model {0:d} clusters msm lag time {1:.3f}'.format(ncluster,msm_lag_time),flush=True)
	objname='tram-{0:d}-clusters-trial-{1:d}-msm-lag-{2:.3f}-cf-{3:.1f}'.format(
		ncluster,itrial,msm_lag_time,connectivity_factor)
	tram_model=pyemma.load(objfname,objname)
	tram_models.append(tram_model)
	print('timescales: ',tram_model.models[tram_model.nthermo-1].timescales(k=nits),flush=True)

print('plotting implied timescales -- {0:d} timescales'.format(nits),flush=True)
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
#unbiased model is the last one in the list
pyemma.plots.plot_memm_implied_timescales(tram_models,ax=ax,units='s',dt=dt*1e-9,therm_state=tram_models[0].nthermo-1,
	nits=nits,annotate=False,marker='o')
for fmt in ['png','pdf']:
	outfname='plots/{0}-msm-implied-timescales-{1:d}-clusters-tica-lag-time-{2:.3f}.{3}'.format(name,ncluster,tica_lag_time,fmt)
	fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)

plt.close(fig)



sys.exit(0)
