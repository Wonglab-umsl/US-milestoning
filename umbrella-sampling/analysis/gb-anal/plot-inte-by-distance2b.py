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

windowcounts={'pyk2-cpd1-amber':421,'pyk2-cpd2-amber':353,'pyk2-cpd8-amber':187}
ligandlabels={'pyk2-cpd1-amber':'ligand 1','pyk2-cpd2-amber':'ligand 2','pyk2-cpd8-amber':'ligand 8'}
#windowcounts={'pyk2-cpd1-amber':300,'pyk2-cpd2-amber':50}
colors={'pyk2-cpd1-amber':'C0','pyk2-cpd2-amber':'C1','pyk2-cpd8-amber':'C4'}
linestyles={'pyk2-cpd1-amber':'-','pyk2-cpd2-amber':'--','pyk2-cpd8-amber':':'}
#these field numbers are in the an2/*-inte-anal file, under "
terms={'vdw':1,'elec':2,'gb':3,'surf':4,'total':7}
terms2=list(terms.keys())+['elec+gb','vdw+surf']
titles={'vdw':'van der Waals','elec':'electrostatic','gb':'polar solvation (GB)','surf':'hydrophobic solvation (SA)',
	'elec+gb':'electrostatic + GB','vdw+surf':'van der Waals + hydrophobic','total':'total interaction'}
interval=int(sys.argv[1])
binsize=float(sys.argv[2])
do_bw=((len(sys.argv)>3) and (sys.argv[3]=='bw'))
#storage='/mnt/stor/ceph/scratch/jmsc87/pyk2'
scratch='/lustre/scratch/jmsc87/pyk2/umbrella-amber/'
#scratch='/storage/hpc/group/wong/scratch/Justin/pyk2/umbrella-amber/'
igroup=1 #for now

figs={}
axes={}
for term in terms2:
	figs[term]=plt.figure()
	axes[term]=figs[term].add_subplot(1,1,1)
for name in windowcounts.keys():
	nwindow=windowcounts[name]
	print('processing data for ',name,flush=True)
	data={term:[] for term in terms2}
	good_windows=[] #ones where data is available
	for iwindow in range(1,nwindow+1):
		fname='{0}/mmpbsa-temp/an/{1}-window-{2:d}-inte-anal'.format(scratch,name,iwindow)
		print('reading data from ',fname,flush=True)
		try:
			with open(fname,'r') as input:
				flag=False #have we reached the "Delta" section yet?
				for line in input:
					#if (line[0:1]=='#'):
					#	continue
					#look for a line that says "DELTA Energy Terms", then skip over it and the next line
					if (line.find('DELTA')>=0):
						flag=True
						input.readline()
						continue
					#if (line.find('Frame')):
					#	continue
					if (not flag):
						continue
					words=line.split(',')
					if (len(words)<2):
						#skip blank and newlines
						continue
					try:
						for term in terms.keys():
							field=terms[term]
							data[term].append(float(words[field]))
					except:
						pdb.set_trace()
			good_windows.append(iwindow)
		except IOError as err:
			#data not available
			print('interaction energy data not available for window ',iwindow,flush=True)
			print('error message ',err,flush=True)
			continue


	ligand_com_data=[]
	#need to read the COM data, and slice by interval
	for iwindow in good_windows:
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
		ligand_com_data.append(this_lig_com_data)
	ligand_com_data=np.stack(ligand_com_data)[:,::interval,:]
	print(ligand_com_data.shape) #(window, frame number (decimated), dimension)
	#need to calculate distances 
	ref_com=ligand_com_data[0,0,:]
	ligand_com_dist=np.sqrt(np.sum((ligand_com_data[:,:,:]-ref_com[:])**2,axis=2))
	#concatenate by window
	ligand_com_dist=np.concatenate(ligand_com_dist,axis=0)


	#this creates one array for each term containing all hte windows -- a one-dimensional array
	for term in terms.keys():
		data[term]=np.array(data[term])
		if (len(data[term])!=len(ligand_com_dist)):
			print('error: inconsistent data lengths for ',term,len(data[term]),len(ligand_com_dist),flush=True)
		print(term,data[term].shape)

	data['elec+gb']=data['elec']+data['gb']
	data['vdw+surf']=data['vdw']+data['surf']

	#sys.exit(-1)
	#try using np.digitize to divide the distances into bins
	maxdist=np.max(ligand_com_dist)
	nbins=int(np.ceil(maxdist/binsize))
	binlimits=np.linspace(start=0,stop=nbins*binsize,num=nbins+1)
	print(binlimits)
	bincenters=(binlimits[0:-1]+binlimits[1:])/2
	print(bincenters)
	bins=np.digitize(ligand_com_dist,binlimits)
	print(bins)

	actual_bincenters=[]
	#average and std by bin 
	mean_by_bin={term:[] for term in terms2}
	std_by_bin={term:[] for term in terms2}
	for ibin in range(0,nbins):
		selected=np.where(bins==ibin)
		if (len(selected)>0):
			actual_bincenters.append(bincenters[ibin])
			for term in terms2:
				mean_by_bin[term].append(np.mean(data[term][selected]))
				std_by_bin[term].append(np.std(data[term][selected])) #/np.sqrt(len(selected[0])))


	#make plots
	if (do_bw):
		print("doing black & white plots",flush=True)
		color='black'
	else:
		color=colors[name]
	ls=linestyles[name]
	for term in terms2:
		mean_by_bin[term]=np.array(mean_by_bin[term])
		std_by_bin[term]=np.array(std_by_bin[term])
		axes[term].plot(actual_bincenters,mean_by_bin[term],label=ligandlabels[name],color=color,ls=ls)
		axes[term].fill_between(actual_bincenters,mean_by_bin[term]-std_by_bin[term],mean_by_bin[term]+std_by_bin[term],
			color=color,alpha=0.2)



#output plots
for term in terms2:
	axes[term].legend(loc='lower right',fontsize=17)
	axes[term].set_xlim(left=0,right=None)
	#axes[term].set_xlabel('ligand center-of-mass distance (A)')
	#if (term=='total'):
	#	axes[term].set_ylabel('total interaction energy (kcal/mol)')
	#else:
	#	axes[term].set_ylabel('interaction energy component (kcal/mol)')
	axes[term].text(0.05,0.9,titles[term],transform=axes[term].transAxes,fontsize=22)
	for fmt in ['png','pdf']:
		outfname='plots/pyk2-allcpds-{0}-by-distance.{1}'.format(term,fmt)
		if (fmt=='pdf'):
			axes[term].set_xlabel('ligand center-of-mass distance (A)')
			if (term=='total'):
				axes[term].set_ylabel('total interaction energy (kcal/mol)')
			else:
				axes[term].set_ylabel('interaction energy component (kcal/mol)')
		figs[term].savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)

