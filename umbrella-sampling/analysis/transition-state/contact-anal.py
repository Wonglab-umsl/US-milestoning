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
import itertools as it
import mdtraj as md

#aadict={'ALA':'A','CYS':'C','ASP':'D','GLU':'E','GLH':'E','PHE':'F','GLY':'G','HID':'H','HIE':'H','HIP':'H','ILE':'I','LYS':'K',
#        'LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y'}


#windowcounts={'pyk2-cpd1-amber':421,'pyk2-cpd2-amber':353,'pyk2-cpd8-amber':187}
ligandlabels={'pyk2-cpd1-amber':'ligand 1','pyk2-cpd2-amber':'ligand 2','pyk2-cpd8-amber':'ligand 8'}
names=ligandlabels.keys()
#windowcounts={'pyk2-cpd1-amber':300,'pyk2-cpd2-amber':50}
colors={'pyk2-cpd1-amber':'C0','pyk2-cpd2-amber':'C1','pyk2-cpd8-amber':'C2'}
dir='/mnt/stor/ceph/scratch/jmsc87/pyk2/umbrella-amber3/dcd-joined/'
scratch='/lustre/scratch/jmsc87/pyk2/umbrella-amber/'
#name=sys.argv[1]
cutoff=float(sys.argv[2])

nres=690-420+1
frac=np.zeros((nres,2*len(names)),dtype=np.float64)

for iname,name in enumerate(names):
	psffname='prmtop/{0}.prmtop'.format(name)
	top=md.load_prmtop(psffname)
	ca_atoms=top.select('name CA')
	ligand=top.select('resname CPD')
	#construct a list of the residue numbers
	resnums=[]
	res_labels=[]
	for ires,idx in enumerate(ca_atoms):
		resnums.append(top.atom(idx).residue.resSeq)
		label=top.atom(idx).residue.code+str(ires+420)
		res_labels.append(label)
	nres=len(resnums)
	bound_mindist_fname='{0}/transition-state/{1}-bound-state-mindist.npy'.format(scratch,name)
	print('loading bound state distance data from ',bound_mindist_fname,flush=True)
	bound_mindist=np.load(bound_mindist_fname)
	ts_mindist_fname='{0}/transition-state/{1}-transition-state-mindist.npy'.format(scratch,name)
	print('loading transition state distance data from ',ts_mindist_fname,flush=True)
	ts_mindist=np.load(ts_mindist_fname) #,mmap_mode='r',allow_pickle=True)
	#pdb.set_trace()
	bound_counts=np.count_nonzero(bound_mindist[:,:]<=cutoff,axis=0)
	ts_counts=np.count_nonzero(ts_mindist[:,:]<=cutoff,axis=0)
	frac[:,2*iname]=bound_counts/float(bound_mindist.shape[0])
	frac[:,2*iname+1]=ts_counts/float(ts_mindist.shape[0])

#figure out which residues make at least one contact, restrict the plot to those
#pdb.set_trace()
maxfrac=np.max(frac,axis=1)
contact_res=np.where(maxfrac>0.01)[0]
limited_frac=frac[contact_res,:]
ncontact=len(contact_res)
print("number of residues with at least 1% contact",ncontact)

limited_res_labels=[res_labels[r] for r in contact_res]

#write out a table in latex format
outfname='pyk2-allcpds-ts-contacts-{0:.1f}'.format(cutoff)
output=open(outfname,'w')
output.write('Residue & ')
for iname,name in enumerate(names):
	output.write('\\multicolumn{2}{c}{'+name+'}')
	output.write('\\\\\n' if (iname==len(names)-1) else ' & ')
output.write('& & ')
for iname in range(0,len(names)):
	output.write(' bound & transition')
	output.write('\\\\\n' if (iname==len(names)-1) else ' & ')
output.write('\\hline\n')
for ires,res_label in enumerate(limited_res_labels):
	output.write(res_label+' & ')
	for iname,name in enumerate(names):
		f_bound=limited_frac[ires,2*iname]*100.0
		f_bound='' if (f_bound==0) else '{0:.1f}'.format(f_bound)
		f_ts=limited_frac[ires,2*iname+1]*100.0
		f_ts='' if (f_ts==0) else '{0:.1f}'.format(f_ts)
		output.write(f_bound+' & '+f_ts)
		output.write('\\\\\n' if (iname==len(names)-1) else ' & ')

#try and make a plot
cm = matplotlib.cm.get_cmap('binary') #,vmin=-m,vmax=m)
heatmap_fig=plt.figure()
heatmap_ax=heatmap_fig.add_subplot(1,1,1)
_ = heatmap_ax.imshow(limited_frac*100.0,cmap=cm,aspect='equal')
#heatmap_ax.set_xticks(range(1,len(bincenters)+1,4))
#heatmap_ax.set_xticklabels(bincenters[::4],fontsize=8)
#heatmap_ax.set_xticklabels(bincenters)

xticklabels=[]
for name in names:
	lig_label=ligandlabels[name]
	xticklabels.append('{0}, bound state'.format(lig_label))
	xticklabels.append('{0}, transition state'.format(lig_label))
heatmap_ax.set_xticks(ticks=range(0,2*len(names)))
heatmap_ax.set_xticklabels(labels=xticklabels,fontsize=8,rotation=90)
heatmap_ax.set_yticks(range(0,ncontact))
heatmap_ax.set_yticklabels(limited_res_labels,fontsize=8)
#heatmap_ax.set_xlabel('ligand COM distance (A)')
heatmap_ax.set_ylabel('residue')
heatmap_ax.set_aspect(1/5) #'equal')

#heatmap_ax.text(0.05,0.9,ligandlabels[name],transform=axes[term].transAxes,fontsize=22)
cb=heatmap_fig.colorbar(_)
cb.set_label(label='% contact') #,fontsize=11) #,fontsize=6) #,pad=0.2,shrink=0.5,aspect=10)
for fmt in ['png','pdf']:
	outfname='plots/pyk2-allcpds-ts-contacts-{0:.1f}'.format(cutoff)
	outfname=outfname.replace('.','_')+'.'+fmt
	heatmap_fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)
plt.close(heatmap_fig)
