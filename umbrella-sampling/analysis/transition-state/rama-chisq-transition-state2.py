#!/usr/bin/env python

import sys
import math
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import chisquare,chi2
import pdb
scratch='/lustre/scratch/jmsc87/pyk2/umbrella-amber/'
dt=0.001 #in ns
binsize=float(sys.argv[1])
nbins=int(360/binsize)
binsize=360.0/nbins
print("binsize, nbins = ",binsize,nbins,flush=True)
#stride=int(sys.argv[2])
res_start=int(sys.argv[2])
res_end=int(sys.argv[3])
nres=res_end-res_start+1
#autocorr=int(sys.argv[4])

def block_average(x,s):
	work=np.full(s,np.nan,dtype=np.float64)
	n=len(x)//s
	if (len(x)%s>0): n+=1
	avg=np.full(n,np.nan,dtype=np.float64)
	for i in range(0,n):
		work[:]=np.nan
		start=s*i
		end=s*(i+1)
		if (end>len(x)): end=len(x)
		work[0:end-start]=x[start:end]
		idx=np.where(work-work[0]>180)
		work[idx]-=360
		idx=np.where(work-work[0]<-180)
		work[idx]+=360
		avg[i]=np.mean(work[0:end-start])
	idx=np.where(avg>180)
	avg[idx]-=360
	idx=np.where(avg<-180)
	avg[idx]+=360
	return avg


names=['pyk2-cpd1-amber','pyk2-cpd2-amber','pyk2-cpd8-amber']
labels={'pyk2-cpd1-amber':'ligand 1','pyk2-cpd2-amber':'ligand 2','pyk2-cpd8-amber':'ligand 8'}
strides=[1,2,5,10]
#fig=plt.figure(constrained_layout=False,figsize=(nres,len(names)))
#spec=gridspec.GridSpec(nrows=len(names),ncols=nres,figure=fig)
ticks=np.linspace(start=-180,stop=180,num=7)

do_rama_plots=True
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
xres=range(res_start,res_end+1)
for iname,name in enumerate(names):

	pvalues=np.full((nres,len(strides)),np.nan,dtype=np.float64)
	#load everything
	psffname='prmtop/{0}.prmtop'.format(name)
	top=md.load_prmtop(psffname)

	dcdfname='{0}/transition-state/{1}-bound-state.dcd'.format(scratch,name)
	print('reading trajectory ',dcdfname,flush=True)
	bound_traj=md.load(dcdfname,top=top,stride=1)

	dcdfname='{0}/transition-state/{1}-transition-state.dcd'.format(scratch,name)
	print('reading trajectory ',dcdfname,flush=True)
	ts_traj=md.load(dcdfname,top=top,stride=1)

	#use mdtraj to calculate the phi and psi angles, convert from radians to degrees
	phi=md.compute_phi(bound_traj)
	phi_atoms=phi[0]
	#print(iatom_by_res)
	bound_phi=phi[1]*(180/np.pi)
	psi=md.compute_psi(bound_traj)
	psi_atoms=psi[0]
	bound_psi=psi[1]*(180/np.pi)
	print(bound_phi.shape,bound_psi.shape)

	phi=md.compute_phi(ts_traj)
	ts_phi=phi[1]*(180/np.pi)
	psi=md.compute_psi(ts_traj)
	ts_psi=psi[1]*(180/np.pi)
	print(ts_phi.shape,ts_psi.shape)

	#construct a dictionary relating residue numbers to indices
	#pdb.set_trace()
	phi_resnum_to_index={}
	psi_resnum_to_index={}
	resnames={}
	for ires,iatom in enumerate(phi_atoms[:,2]):
		residue=top.atom(iatom).residue.resSeq+420
		phi_resnum_to_index[residue]=ires
		resnames[residue]=top.atom(iatom).residue.name
	for ires,iatom in enumerate(psi_atoms[:,1]):
		residue=top.atom(iatom).residue.resSeq+420
		psi_resnum_to_index[residue]=ires
	#print(resnum_to_index)	
	#print(phi)
	#print(psi)
	#sys.exit(0)
	fname='{0}-rama-chisq'.format(name)
	output=open(fname,'w')
	for res in xres:
		output.write('{0:d}'.format(res))
		if ((res not in phi_resnum_to_index) or (res not in psi_resnum_to_index)):
			continue
		phi_ires=phi_resnum_to_index[res]
		psi_ires=psi_resnum_to_index[res]
		
		for istride,stride in enumerate(strides):
			print('residue ',res,' stride ',stride,flush=True)
			if (stride>1):
				bound_phi_avg=block_average(bound_phi[:,phi_ires],stride)
				bound_psi_avg=block_average(bound_psi[:,psi_ires],stride)
				ts_phi_avg=block_average(ts_phi[:,phi_ires],stride)
				ts_psi_avg=block_average(ts_psi[:,psi_ires],stride)
			else:
				#optimization -- don't bother block averaging if stride=1
				bound_phi_avg=bound_phi[:,phi_ires]
				bound_psi_avg=bound_psi[:,psi_ires]
				ts_phi_avg=ts_phi[:,phi_ires]
				ts_psi_avg=ts_psi[:,psi_ires]

			bound_hist=np.histogram2d(bound_phi_avg,bound_psi_avg,bins=nbins,range=[[-180.0,180.0],[-180.0,180.0]])[0]
			ts_hist=np.histogram2d(ts_phi_avg,ts_psi_avg,bins=nbins,range=[[-180.0,180.0],[-180.0,180.0]])[0]
			bound_nframes=np.sum(bound_hist)
			ts_nframes=np.sum(ts_hist)
			#formula 14.3.3 in Numerical Recipes -- chi^2 = sum_i [sqrt(S/R)*R_i - sqrt(R/S)*S_i]^2/(R_i+S_i)
			#sum does not include bins where R_i+S_i=0, R=sum_i R_i, S=sum_i S_i (total number of data points for each)
			#R = transition state, S = bound_state
			bound_hist=np.ndarray.flatten(bound_hist)
			ts_hist=np.ndarray.flatten(ts_hist)
			sum_hist=bound_hist[:]+ts_hist[:]
			good_bins=np.where(sum_hist>0)[0]
			ndf=len(good_bins)-1 #I think, need to check this
			sum_hist=sum_hist[good_bins] #eliminate empty bins
			sqrtrs=np.sqrt(ts_nframes/bound_nframes)
			sqrtsr=1/sqrtrs
			diff_hist=(sqrtsr*ts_hist[good_bins]-sqrtrs*bound_hist[good_bins])**2
			chisq=np.sum(diff_hist/sum_hist)
			#Lucy correction for small counts, eqs. 14.3.14 and eq. 14.3.10 Numerical Recipes
			invsum=((ts_nframes+bound_nframes)**2/(ts_nframes*bound_nframes)-6)*np.sum(1.0/sum_hist)
			var_corr_factor=np.sqrt(2*ndf/(2*ndf+invsum))
			chisq=ndf+var_corr_factor*(chisq-ndf)
			#the p-value is 1 - cdf for chi square with ndf degrees of freedom
			#pvalue=chi2.sf(chisq,ndf)
			logpvalue=chi2.logsf(chisq,ndf)/np.log(10) #to convert to common logarthim
			#try using one term of continued fraction expansion for incomplete gamma function
			#eq. 6.2.7 Numerical Recipes, with log taken
			if (logpvalue<=-np.inf): 
				logpvalue=-0.5*chisq+0.5*ndf*np.log(0.5*chisq)-np.log(0.5*chisq+1-0.5*ndf)-math.lgamma(0.5*ndf)
				logpvalue/=np.log(10)
			pvalues[res-res_start,istride]=logpvalue
			output.write(' {0:.2f} {1:d} {2:.2f}'.format(chisq,ndf,logpvalue))
		#if (pvalue<1e-50):
			#pdb.set_trace()
		#	pass
		#line='{0:d} {1:.2f} {2:d} {3:.2f}\n'.format(res,chisq,ndf,logpvalue)
		#output.write(line)
		output.write('\n')
		if (do_rama_plots and (np.max(pvalues[res-res_start,:])<np.log10(0.05/nres))):
			#statistically-significant
			rama_fig=plt.figure()
			rama_ax=rama_fig.add_subplot(1,1,1)
			rama_ax.scatter(bound_phi[:,phi_ires],bound_psi[:,psi_ires],marker='.',s=0.01,c='C0',label='bound state')
			rama_ax.scatter(ts_phi[:,phi_ires],ts_psi[:,psi_ires],marker='.',s=0.01,c='C1',label='transition state')
			rama_ax.set_xlim(left=-180,right=180)
			rama_ax.set_ylim(bottom=-180,top=180)
			rama_ax.set_xticks(ticks)
			rama_ax.set_yticks(ticks)
			rama_ax.legend(loc='lower right',fontsize=8)
			rama_ax.tick_params(labelsize=8)
			rama_ax.grid()
			outfname='rama-plots/{0}-rama-{1:d}.png'.format(name,res)
			rama_fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
			plt.close(rama_fig)

	#pdb.set_trace() 
	#choose from among the various strides the one that gives the maximum p-value
	#corresponding to the largest standard error
	pvalues=np.max(pvalues,axis=1)
	ax.plot(xres,-pvalues,label=labels[name])
	output.close()
#finish and save the plot
#ax.set_yscale('log')
ax.set_xlabel('residue number')
ax.set_ylabel('-$\\log_{10}(p)$')
ax.legend(loc='upper right',fontsize=11)
for fmt in ['png','pdf']:
	outfname='pyk2-allcpds-rama-chisq-{0:d}-{1:d}.{2}'.format(res_start,res_end,fmt)
	fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)
plt.close(fig)

