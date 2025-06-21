#!/usr/bin/env python
import sys
import os
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import pdb
import warnings
#this program sometimes gives warnings due to taking logarithms of zero probability
#warnings.filterwarnings("error")
kt=1.987e-3*300
#this is the area under a standard normal distribution from -1 to 1 SD
confidence=0.6826894809

name=sys.argv[1]
milestonefile=sys.argv[2]
product_milestone_file=sys.argv[3]
nbootstrap=int(sys.argv[4])

#milestonefile='{0}-milestones-3'.format(name)
initial_milestones=[]
input=open(milestonefile,'r')
for line in input:
	words=line.split()
	ianchor=int(words[1])
	janchor=int(words[2])
	initial_milestones.append((ianchor,janchor))
input.close()
print(initial_milestones)
nmilestone=len(initial_milestones)

#read list of product milestones
print("reading product milestones from file ",product_milestone_file,flush=True)
input=open(product_milestone_file,'r')
product_milestones=[]
for line in input:
	words=line.split()
	iprod_anchor=int(words[0])
	jprod_anchor=int(words[1])
	try:
		product_milestone=initial_milestones.index((iprod_anchor,jprod_anchor))
	except ValueError:
		print("product anchors {0:d} {1:d} are not one of the milestones".format(iprod_anchor,jprod_anchor))
		sys.exit(-1)
	product_milestones.append(product_milestone)

nprod_milestone=len(product_milestones)
print("{0:d} product milestones detected".format(nprod_milestone),flush=True)




#sys.exit(-1)




data=[]

logfile='{0}-log'.format(name)
input=open(logfile,'r')
for line in input:
	#format of log file: name ianchor janchor id inewanchor jnewanchor newid stoptime
	#field number        0    1       2       3  4          5          6     7
	words=line.split()
	ianchor=int(words[1])
	janchor=int(words[2])
	t1=(ianchor,janchor)
	if (t1 in initial_milestones):
		imilestone=initial_milestones.index(t1)
	else:
		print("log file has initial milestone not in list ",line)
		continue
		
	inewanchor=int(words[4])
	jnewanchor=int(words[5])
	t2=(inewanchor,jnewanchor)
	stoptime=float(words[7])
	t=(t1,t2,stoptime)
	data.append(t)
input.close()

ndata=len(data)
print(ndata)
print(data[0:5])

#instead of calculating the free energy profile (which gives errors due to zero p_norm)
#we will stack p_norms
free_energy_samples=[]
#p_norm_samples=[]
mfpt_samples=[] #np.full(nbootstrap,np.nan,dtype=np.float64)

for ibootstrap in range(0,nbootstrap+1):
	if (ibootstrap==0):
		print("using all data")
		bootstrap_sample=np.array(range(0,ndata))
	else:
		print("beginning bootstrap sample ",ibootstrap)
		bootstrap_sample=np.sort(np.random.choice(ndata,ndata,replace=True))
	print(bootstrap_sample)
	#accumulate n_ij and n_i
	n_ij=np.zeros((nmilestone,nmilestone),dtype=np.int64)
	n_i_incl_new=np.zeros(nmilestone,dtype=np.int64)
	milestone_lifetime=np.zeros(nmilestone,dtype=np.float64)
	for i in range(0,ndata):
		t=data[bootstrap_sample[i]]
		#format of t: ((ianchor, janchor), (inewanchor, jnewanchor), stoptime)
		imilestone=initial_milestones.index(t[0])
		n_i_incl_new[imilestone]+=1
		milestone_lifetime[imilestone]+=t[2]
		if (t[1] in initial_milestones):
			jmilestone=initial_milestones.index(t[1])
			#n_ij only counts trajectories that land on original milestone
			n_ij[imilestone,jmilestone]+=1
	milestone_lifetime[:]=milestone_lifetime[:]/n_i_incl_new[:]
	#print(milestone_lifetime)
	#print(n_ij)
	#k_ij = n_ij / n_i -- we do not use n_i_incl_new here because it includes the trajectories
	#that reach new milestones, we don't want that
	n_i_old_only=np.sum(n_ij,axis=1)
	print(n_i_old_only)
	kmatrix = n_ij[:,:]/np.expand_dims(n_i_old_only,axis=1)
	kmatrix_test=np.sum(kmatrix,axis=1)
	if (np.any(np.isnan(kmatrix)) or not (np.allclose(kmatrix_test,1))):
		print("error in kmatrix")
		continue


	#to calculate free energies, from Brajesh's script -- calculate eigenvalues and eigenvectors
	#find the eigenvalue closest to 1, and ensure it is positive. then multiply elementwise by the milestone lifetime
	k_trans = np.transpose(kmatrix)
	e_v, e_f = np.linalg.eig(k_trans)
	#print(e_v)
	idx = np.abs(e_v - 1).argmin()
	#print('eigenvalue: ',e_v[idx])
	if np.all(e_f[:, idx] > 0):
		q = np.abs(e_f[:, idx])
	else:
		q = np.abs(-1 * e_f[:, idx])
	p = np.transpose(q) * milestone_lifetime
	p_norm = p / np.sum(p)
	#print(p)
	#free_energy_samples[ibootstrap,:] = -1.0 * kt * np.log(p_norm)  
	#p_samples.append(p_norm)

	kmatrix_copy=np.copy(kmatrix)
	kmatrix_copy[product_milestones,:]=0.0
	milestone_lifetime_copy=np.copy(milestone_lifetime)
	milestone_lifetime_copy[product_milestones]=0.0
	mat=np.eye(nmilestone)-kmatrix_copy
	condnum=np.linalg.cond(mat)
	#make sure the initial estimate is always included
	if ((ibootstrap==0) or (condnum<1e15)):
		#the linear system is well conditioned
		aux=np.linalg.solve(mat,milestone_lifetime_copy)
		if (aux[0]<0):
			print("negative mfpt")
		#	pdb.set_trace()
		mfpt_samples.append(aux[0])
		#try:
		free_energy_samples.append(-1.0*kt*np.log(p_norm))
		#except:
		#	pdb.set_trace()
		#	pass
		#p_norm_samples.append(p_norm)
		if (ibootstrap==0):
			#calculate the committors
			kmatrix_copy2=np.copy(kmatrix)
			##### Condition on Reactant milestone  :: assuming reactant milestone is at 0 th position in kmatrix
			kmatrix_copy2[0,:]=0
			####    condition on product/LAST MILESTONE 
			kmatrix_copy2[product_milestones,:]=0
			kmatrix_copy2[product_milestones,product_milestones]=1
			##### Committor Function
			c=np.linalg.matrix_power(kmatrix_copy2,4294967296)
			c_val=c[:, product_milestones]
			#This matrix now gives, for each milestone, the probability of reaching each product milestone 
			#before reaching the reactant milestone.  We want the probability of reaching any product milestone
			#before reaching the reactant milestone.  Since once we reach a product milestone we can't reach any other,
			#they are mutually exclusive.  We can add probabilities.
			c_val=np.sum(c_val,axis=1)
			print(c_val)

			#write out committors to a file
			suffix=product_milestone_file.replace(name+'-','')
			outfname='{0}-committors-{1}'.format(name,suffix)
			output=open(outfname,'w')
			for i,m in enumerate(initial_milestones):
				line='{0:d} {1:d} {2:d} {3:f}\n'.format(i+1,m[0],m[1],c_val[i])
				output.write(line)
			output.close()

			#try to plot it
			np_milestones=np.array(initial_milestones)
			cm = matplotlib.cm.get_cmap('viridis')
			fig=plt.figure()
			ax=fig.add_subplot(1,1,1)
			_ = ax.scatter(np_milestones[:,0],np_milestones[:,1],c=c_val,cmap=cm,s=20.0)
			_ = ax.scatter(np_milestones[:,1],np_milestones[:,0],c=c_val,cmap=cm,s=20.0)
			ax.set_aspect('equal')
			ax.set_xlabel('anchor number')
			ax.set_ylabel('anchor number')
			fig.colorbar(_)
			for fmt in ['png','pdf']:
				outfname='{0}-committors-{1}.{2}'.format(name,suffix,fmt)
				fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)
			plt.close(fig)	

			committors_by_anchor=np.full(np.max(np_milestones),np.nan,dtype=np.float64)
			for i,m in enumerate(initial_milestones):
				if ((m[1]-m[0])==1):
					committors_by_anchor[m[0]]=c_val[i]
			fig=plt.figure()
			ax=fig.add_subplot(1,1,1)
			ax.plot(range(1,np.max(np_milestones)),committors_by_anchor[1:])
			ax.set_xlabel('anchor number')
			ax.set_ylabel('committor value')
			for fmt in ['png','pdf']:
				outfname='{0}-committors1d-{1}.{2}'.format(name,suffix,fmt)
				fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)
			plt.close(fig)
			#sys.exit(-1)
	else:
		d=np.linalg.det(mat)
		print("ill conditioned: condnum = ",condnum," det = ",d)
		#print(mat)
		#u, s, vh = np.linalg.svd(mat)
		#u2, s2, vh2 = np.linalg.svd(kmatrix)
		#print(s)
		#print(s2)
		#pdb.set_trace()
		#pass
	del n_i_incl_new
	del n_ij
	del milestone_lifetime

#The distribution of the free energies and especially the MFPT may not be Gaussian.
#therefore, use 68% confidence intervals instead of the mean and SD to characterize it.
n_successful_bootstrap=len(free_energy_samples)
idx_ci_lo=int(np.floor(n_successful_bootstrap*(0.5-confidence/2)))
idx_ci_hi=int(np.ceil(n_successful_bootstrap*(0.5+confidence/2)))
#working directly with the probability distributions avoids errors due to zero probability
#the free energy is a monotonically decreasing functino of the probability
free_energy_samples=np.vstack(free_energy_samples)
#p_norm_samples=np.vstack(p_norm_samples)
#print(free_energy_samples.shape)
sorted_fe_samples=np.sort(free_energy_samples,axis=0)
#print(sorted_fe_samples)
#hopefully the confidence interval will not include zeros
fe_ci_lo=sorted_fe_samples[idx_ci_lo,:]
fe_ci_hi=sorted_fe_samples[idx_ci_hi,:]

#put this out to a file, for further plotting
outfname='{0}-fe-from-milestoning'.format(name)
output=open(outfname,'w')
for i,m in enumerate(initial_milestones):
        line='{0:d} {1:d} {2:d} {3:.6f} {4:.6f} {5:.6f}\n'.format(i,m[0],m[1],free_energy_samples[0,i],fe_ci_lo[i],fe_ci_hi[i])
        output.write(line)   
output.close()


#convert MFPT from ns to s
mfpt_samples=np.array(mfpt_samples)*1e-9
mfpt_all_data=mfpt_samples[0]


sorted_mfpt_samples=np.sort(mfpt_samples)
if (n_successful_bootstrap%2==0):
	#python uses zero based indices, so if (say) n==6, then n//2=3 but we want indices 2 and 3
	mfpt_median=(sorted_mfpt_samples[n_successful_bootstrap//2-1]+sorted_mfpt_samples[n_successful_bootstrap//2])/2
else:
	mfpt_median=sorted_mfpt_samples[n_successful_bootstrap//2]
mfpt_ci_lo=sorted_mfpt_samples[idx_ci_lo]
mfpt_ci_hi=sorted_mfpt_samples[idx_ci_hi]

#plot the sampling distribution of the mfpt estimates
nbins=int(np.ceil(len(mfpt_samples)/50))
pmin=np.floor(np.log10(np.min(mfpt_samples)))
pmax=np.ceil(np.log10(np.max(mfpt_samples)))
bins=np.logspace(start=pmin,stop=pmax,num=nbins)
#print(bins)
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.hist(mfpt_samples,bins=bins)
ax.set_xscale('log')
ax.set_xlabel('MFPT (s)')
ax.set_ylabel('count of MFPT values')
ax.axvline(x=mfpt_all_data,color='k',linestyle='-')
ax.axvline(x=mfpt_median,color='k',linestyle='--')
ax.axvspan(xmin=mfpt_ci_lo,xmax=mfpt_ci_hi,color='k',alpha=0.2)
major_ticks=[10**i for i in range(int(pmin),int(pmax)+1)]
ax.set_xticks(major_ticks)
outfname='{0}-mfpt-dist-{1:d}-{2:d}.png'.format(name,product_ianchor,product_janchor)
fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
plt.close()



#there is actually one less "successful" bootstrap than n_successful_bootstrap
#because the "all data" estimate is automatically included
print("total number of bootstraps: ",nbootstrap)
print("fraction successful bootstraps: {0:.3f}".format((n_successful_bootstrap-1)/nbootstrap))
print("total number of data points: ",ndata)
print("MFPT estimate from all data (s): {0:.3g}".format(mfpt_all_data))
print("rate estimate from all data (s^-1): {0:.3g}".format(1.0/mfpt_all_data))
print("median MFPT, 68% CI (s): {0:.3g} {1:.3g} {2:.3g}".format(mfpt_median,mfpt_ci_lo,mfpt_ci_hi))
print("median rate, 68% CI (s^-1): {0:.3g} {1:.3g} {2:.3g}".format(1.0/mfpt_median,1.0/mfpt_ci_hi,1.0/mfpt_ci_lo))
halfwidth=np.log10(mfpt_ci_hi/mfpt_ci_lo)/2
print("half 68% CI width in orders of magnitude: {0:.3f}".format(halfwidth))
sys.exit(-1)
