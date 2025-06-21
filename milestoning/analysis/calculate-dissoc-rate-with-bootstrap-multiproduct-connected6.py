#!/usr/bin/env python
import sys
import os
import numpy as np
import numpy.ma as ma
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import networkx as nx
import pdb
import warnings
import pickle
#this program sometimes gives warnings due to taking logarithms of zero probability
#warnings.filterwarnings("error")
kt=1.987e-3*300
#this is the area under a standard normal distribution from -1 to 1 SD
confidence=0.6826894809

def portray(center,low,high):
        expo=int(np.floor(np.log10(center)))
        mcenter=center/10**expo
        #dhi=(high-center)/10**expo
        #dlo=(center-low)/10**expo
        mhi=high/10**expo
        mlo=low/10**expo
        s='{0:.3f} {1:.3f} {2:.3f} x 10^{3:d}'.format(mcenter,mlo,mhi,expo)
        return s


name=sys.argv[1]
logfile=sys.argv[2]
milestonefile=sys.argv[3]
start_ianchor=int(sys.argv[4])
start_janchor=int(sys.argv[5])
product_milestone_file=sys.argv[6]
#product_ianchor=int(sys.argv[3])
#product_janchor=int(sys.argv[4])
nbootstrap=int(sys.argv[7])

#if (name=='pyk2-cpd1-amber'):
#	milestonefile='{0}-milestones-3'.format(name)
#elif (name=='pyk2-cpd2-amber'):
#	milestonefile='{0}-milestones-2'.format(name)
#else:
#	print('unrecognized name')
#	sys.exit(-1)

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
#pdb.set_trace()

#sys.exit(-1)




data=[]

#logfile='{0}-log'.format(name)
input=open(logfile,'r')
for line in input:
	#format of log file: name ianchor janchor id inewanchor jnewanchor newid stoptime
	#field number        0    1       2       3  4          5          6     7
	words=line.split()
	ianchor=int(words[1])
	janchor=int(words[2])
	id=int(words[3])
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
	t=(t1,t2,stoptime,id)
	data.append(t)
input.close()

#try reading the incomplete file for censored data
nincomplete=0
incompletelogfile='{0}-incomplete'.format(name)
try:
	input=open(incompletelogfile,'r')
except OSError:
	#not found, skip it
	pass
else:
	for line in input:
		#format: name ianchor janchor id stoptime
		#fields  0    1       2       3  4
		words=line.split()
		ianchor=int(words[1])
		janchor=int(words[2])
		id=int(words[3])
		t1=(ianchor,janchor)
		if (t1 in initial_milestones):
			imilestone=initial_milestones.index(t1)
		else:
			#print("log file has initial milestone not in list ",line)
			continue
		stoptime=float(words[4])
		t=(t1,None,stoptime,id)
		#check to make sure it is not already included 
		l=[tt for tt in data if ((tt[0]==t1) and (tt[3]==id))]
		if (len(l)<=0):
			data.append(t)
			nincomplete+=1

ndata=len(data)
print(ndata)
#print(data[0:5])
#sys.exit(-1)

#generate the milestone graph, test it for connectivity, and if necessary 
starting_milestone=0
milestone_graph=nx.DiGraph()
milestone_graph.add_nodes_from(range(0,nmilestone))
orig_n_ij=np.zeros((nmilestone,nmilestone))
for i in range(0,ndata):
	t=data[i]
	#format of t: ((ianchor, janchor), (inewanchor, jnewanchor), stoptime)
	imilestone=initial_milestones.index(t[0])
	if ((t[1] is not None) and (t[1] in initial_milestones)):
		jmilestone=initial_milestones.index(t[1])
		#milestone_graph.add_edge(imilestone,jmilestone)
		orig_n_ij[imilestone,jmilestone]+=1


#breakpoint()
n_ij_mask=(orig_n_ij<=0)
masked_n_ij=ma.masked_array(orig_n_ij,n_ij_mask)
min_n_ij=int(np.min(masked_n_ij))
print(f'minimum transition count: {min_n_ij:d}',flush=True)
count_one_transition=np.count_nonzero(orig_n_ij==1)
count_nonzero_transition=np.count_nonzero(orig_n_ij>0)
print(f'exactly one/at least one transition milestone pairs: {count_one_transition:d}/{count_nonzero_transition:d}',flush=True)
 

#require at least two simulations to connect milestones
for i in range(0,nmilestone):
	for j in range(0,nmilestone):
		if (orig_n_ij[i,j]>0):
			milestone_graph.add_edge(i,j)

good_milestones=set(product_milestones)
for prod_milestone in product_milestones:
	good_milestones=good_milestones.union(nx.ancestors(milestone_graph,prod_milestone))


#test it for connectivity.  If necessary, identify the good milestones 
#pdb.set_trace()
if (not nx.is_strongly_connected(milestone_graph)):
	print("warning: milestone network is not strongly connected",flush=True)
	#each component is a set of milestone indices
	comps=list(nx.strongly_connected_components(milestone_graph))
	#get information on each of the strongly connected components
	output=open('{0}-connected-comps'.format(name),'w')
	for i,c in enumerate(comps):
		l=[initial_milestones[m] for m in c]
		output.write('{0:d}: {1}\n'.format(i,str(l)))
	output.close()
	cond_graph=nx.condensation(milestone_graph)
	fig=plt.figure()
	ax=fig.add_subplot(1,1,1)
	nx.draw_networkx(cond_graph,ax=ax,font_size=4,node_size=150)
	ax.set_aspect('equal')
	outfname='{0}-milestones-condensed.png'.format(name)
	fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
	plt.close(fig)

#sys.exit(-1)

#we will not choose a starting milestone now. rather calculate the mfpt for all the good milestones
#the starting milestone is the one 

#if (starting_milestone not in good_milestones):
#	print("warning: starting milestone not connected to product milestones",flush=True)
#	starting_milestone=min(good_milestones)
#	print("new starting milestone: ",initial_milestones[starting_milestone],flush=True)
ngoodmilestone=len(good_milestones)
if (ngoodmilestone<nmilestone):
	print("reduced matrix to ",ngoodmilestone," good milestones",flush=True)
good_milestones=list(good_milestones)
#starting_good_milestone=good_milestones.index(starting_milestone)
good_product_milestones=[]
for m in product_milestones:
	if (m in good_milestones):
		good_product_milestones.append(good_milestones.index(m))

#instead of calculating the free energy profile (which gives errors due to zero p_norm)
#we will stack p_norms
free_energy_samples=[]
#p_norm_samples=[]
mfpt_samples=[] #np.full(nbootstrap,np.nan,dtype=np.float64)
starting_milestone=0
#min_trans_counts=[]
count_one_transition=[]
for ibootstrap in range(0,nbootstrap+1):
	if (ibootstrap==0):
		print("using all data",flush=True)
		bootstrap_sample=np.array(range(0,ndata))
	else:
		print("beginning bootstrap sample ",ibootstrap,flush=True)
		bootstrap_sample=np.sort(np.random.choice(ndata,ndata,replace=True))
	#print(bootstrap_sample)
	#accumulate n_ij and n_i -- we could rewrite these loops to use only the good milestones. but it is easier just to reduce them afterwards
	n_ij=np.zeros((nmilestone,nmilestone),dtype=np.int64)
	n_i_incl_new=np.zeros(nmilestone,dtype=np.int64)
	milestone_lifetime=np.zeros(nmilestone,dtype=np.float64)
	for i in range(0,ndata):
		t=data[bootstrap_sample[i]]
		#format of t: ((ianchor, janchor), (inewanchor, jnewanchor), stoptime)
		imilestone=initial_milestones.index(t[0])
		#always include in numerator of lifetime average, even if incomplete
		milestone_lifetime[imilestone]+=t[2]
		#but do not count in denominator if not complete (censored sample MLE estimate)
		if (t[1] is not None):
			n_i_incl_new[imilestone]+=1
			if (t[1] in initial_milestones):
				jmilestone=initial_milestones.index(t[1])
				#n_ij only counts trajectories that land on original milestone
				n_ij[imilestone,jmilestone]+=1
	#try to reduce the set of milestones so that there are no all-zero columns in n_ij 
	#(that is, only those milestones that have been reached at least once
	orig_n_ij=n_ij
	#statistics on the minimum number of 
	#min_trans_counts.append(np.min(ma.masked_array(orig_n_ij,n_ij_mask)))
	count_one_transition.append(np.count_nonzero(orig_n_ij==1))
	#use only the rows and columns of n_ij that correspond to the largest strongly connected component
	n_ij=n_ij[good_milestones,:][:,good_milestones]
	#good_milestones indexing should be in terms of the original milestones
	#need to reduce all the other arrays
	#pdb.set_trace()
	milestone_lifetime=milestone_lifetime[good_milestones]
	n_i_incl_new=n_i_incl_new[good_milestones]
	#convert the set of product milestones to be in the reduced list of milestones
	#can keep going as before
	milestone_lifetime[:]=milestone_lifetime[:]/n_i_incl_new[:]
	#print(milestone_lifetime)
	if (ibootstrap==0):
		#print(n_ij)
		n_in_k_matrix=np.sum(orig_n_ij)
	#k_ij = n_ij / n_i -- we do not use n_i_incl_new here because it includes the trajectories
	#that reach new milestones, we don't want that
	n_i_old_only=np.sum(n_ij,axis=1)
	#print(n_i_old_only)
	kmatrix = n_ij[:,:]/np.expand_dims(n_i_old_only,axis=1)
	kmatrix_test=np.sum(kmatrix,axis=1)
	if (np.any(np.isnan(kmatrix)) or not (np.allclose(kmatrix_test,1))):
		print("error in kmatrix")
		#pdb.set_trace()
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
	#this needs to be converted back to original milestones!
	p_norm = p / np.sum(p)
	#print(p)
	#free_energy_samples[ibootstrap,:] = -1.0 * kt * np.log(p_norm)  
	#p_samples.append(p_norm)
	#pdb.set_trace()
	kmatrix_copy=np.copy(kmatrix)
	kmatrix_copy[good_product_milestones,:]=0.0
	milestone_lifetime_copy=np.copy(milestone_lifetime)*1e-9
	milestone_lifetime_copy[good_product_milestones]=0.0
	mat=np.eye(ngoodmilestone)-kmatrix_copy
	condnum=np.linalg.cond(mat)
	#make sure the initial estimate is always included
	if ((ibootstrap==0) or (condnum<1e15)):
		if (condnum>1e15):
			d=np.linalg.det(mat)
			print("ill conditioned: condnum = ",condnum," det = ",d)
		#the linear system is well conditioned, factor of 1e-9 converts to s
		#pdb.set_trace()
		if (ibootstrap==0):
			#plot contributions to mfpt
			#mfpt = e_0 . (I-K')^-1 . <t> = sum_j [(I-K')^-1]_{0,j}  * <t>_j
			#pdb.set_trace()
			fig=plt.figure()
			ax=fig.add_subplot(1,1,1)
			x=list(range(ngoodmilestone))
			ax.step(x,milestone_lifetime_copy,where='mid',label='$\\langle t \\rangle_j$ (s)')
			matinv=np.linalg.inv(mat)
			ax.step(x,matinv[starting_milestone,:],where='mid',label='$(I-K\')^{-1}_{0j}$')
			prod=matinv[starting_milestone,:]*milestone_lifetime_copy[:]
			ax.step(x,prod,where='mid',label='$(I-K\')^{-1}_{0j} \\langle t \\rangle_j$ (s)')
			milestone_labels=[initial_milestones[m] for m in good_milestones]
			ax.set_yscale('log')
			ax.set_xticks(x,milestone_labels)
			ax.tick_params(axis='x',labelsize=8,labelrotation=90)
			ax.set_xlabel('milestone $j$')
			ax.set_ylabel('MFPT contribution (s)')
			ax.legend(loc='best',fontsize=11)
			outfname='{0}-mfpt-components.png'.format(name)
			fig.savefig(outfname,format='png',bbox_inches='tight',dpi=600)
			plt.close(fig)
			outfname='{0}-mfpt-components'.format(name)
			output=open(outfname,'w')
			for im,m in enumerate(milestone_labels):
				line='{0:d} {1:d} {2} {3} {4}\n'.format(m[0],m[1],
					matinv[starting_milestone,im],milestone_lifetime_copy[im],prod[im])
				output.write(line)
			output.close()
			#matinv_sorted=np.sort(np.ravel(np.linalg.inv(mat)))
		mfpt_all=np.linalg.solve(mat,milestone_lifetime_copy)
		if (np.any(mfpt_all<0)):
			print("warning: negative mfpt")
		#	pdb.set_trace()
		mfpt_samples.append(mfpt_all)
		free_energy_samples.append(-1.0*kt*np.log(p_norm))
		#except:
		#	pdb.set_trace()
		#	pass
		#p_norm_samples.append(p_norm)
	else:
		d=np.linalg.det(mat)
		print("ill conditioned: condnum = ",condnum," det = ",d)
		u, s, vt=np.linalg.svd(mat)
		#pdb.set_trace()
		print("highest and lowest sing. values: ",s[0],s[-1])
	del n_i_incl_new
	del n_ij
	del milestone_lifetime

pdb.set_trace()
fname='{0}-bootstrap-distributions.p'.format(name)
output=open(fname,'wb')
actual_good_milestones=[initial_milestones[i] for i in good_milestones]
pickle.dump(actual_good_milestones,output)
pickle.dump(free_energy_samples,output)
rate_samples=1.0/np.vstack(mfpt_samples)
pickle.dump(rate_samples,output)
output.close()


n_successful_bootstrap=len(free_energy_samples)
if (n_successful_bootstrap==0):
	print('no successful bootstraps')
	sys.exit(-1)


idx_ci_lo=int(np.ceil(n_successful_bootstrap*(0.5-confidence/2)))-1
idx_ci_hi=int(np.ceil(n_successful_bootstrap*(0.5+confidence/2)))-1
#working directly with the probability distributions avoids errors due to zero probability
#the free energy is a monotonically decreasing functino of the probability
free_energy_samples=np.vstack(free_energy_samples)
#choose the starting milestone to be the good milestone with the lowest free energy

#p_norm_samples=np.vstack(p_norm_samples)
#print(free_energy_samples.shape)
sorted_fe_samples=np.sort(free_energy_samples,axis=0)
#print(sorted_fe_samples)
#hopefully the confidence interval will not include zeros
fe_ci_lo=sorted_fe_samples[idx_ci_lo,:]
fe_ci_hi=sorted_fe_samples[idx_ci_hi,:]
#consider calculating the median fe estimate
#put this out to a file, for further plotting
outfname='{0}-fe-from-milestoning-bootstrap'.format(name)
output=open(outfname,'w')
#the indexing is in terms of the good milesotnes
for igood,iactual in enumerate(good_milestones):
	m=initial_milestones[iactual]
	line='{0:d} {1:d} {2:d} {3:.6f} {4:.6f} {5:.6f}\n'.format(iactual,m[0],m[1],
		free_energy_samples[0,igood],fe_ci_lo[igood],fe_ci_hi[igood])
	output.write(line)
output.close()



try:
	starting_actual_milestone=initial_milestones.index((start_ianchor,start_janchor))
	starting_good_milestone=good_milestones.index(starting_actual_milestone)
except ValueError:
	print('{0} is not sufficiently connected to the product milestone, substituting'.format((start_ianchor,start_janchor)),flush=True)
	#choose the milestone with lowest free energy that is not a product milestone)
	good_non_product_milestones=[m for m in range(0,ngoodmilestone) if good_milestones[m] not in product_milestones]
	starting_good_milestone=good_non_product_milestones[np.argmin(free_energy_samples[0,good_non_product_milestones])]
	#starting_good_milestone=np.argmin(free_energy_samples[0,:])
	starting_actual_milestone=good_milestones[starting_good_milestone]
	print("milestone with lowest free energy: {0:d} {1} {2:.3f}".format(
		starting_actual_milestone,str(initial_milestones[starting_actual_milestone]),
		free_energy_samples[0,starting_good_milestone]),flush=True)





#mfpt_all_data=mfpt_samples[0]
#The distribution of the free energies and especially the MFPT may not be Gaussian.
#therefore, use 68% confidence intervals instead of the mean and SD to characterize it.
if (n_successful_bootstrap==1):
	#we cannot calculate confidence intervals
	print("total number of bootstraps: ",nbootstrap)
	#print("fraction successful bootstraps: {0:.3f}".format((n_successful_bootstrap-1)/nbootstrap))
	#print("total number of data points: ",ndata)
	print("total number of data points, complete, incomplete sims: ",ndata,ndata-nincomplete,nincomplete)
	print("total sims in k matrix (not hitting new milestone): ",n_in_k_matrix)
	#convert MFPT from ns to s
	mfpt_samples=np.array(mfpt_samples)
	mfpt_all_data=mfpt_samples[0]
	print("MFPT estimate from all data (s): {0:.3g}".format(mfpt_all_data[starting_good_milestone]))
	print("rate estimate from all data (s^-1): {0:.3g}".format(1.0/mfpt_all_data[starting_good_milestone]))
	sys.exit(0)

#breakpoint()
#extract confidence intervals for 
rate_samples=1.0/np.vstack(mfpt_samples)
sorted_rate_samples=np.sort(rate_samples,axis=0)
#convert MFPT from ns to s, convert to rates
rate_all_data=rate_samples[0,:]
if (n_successful_bootstrap%2==0):
	#python uses zero based indices, so if (say) n==6, then n//2=3 but we want indices 2 and 3
	rate_median=(sorted_rate_samples[n_successful_bootstrap//2-1,:]+sorted_rate_samples[n_successful_bootstrap//2,:])/2
else:
	rate_median=sorted_rate_samples[n_successful_bootstrap//2,:]
rate_ci_lo=sorted_rate_samples[idx_ci_lo,:]
rate_ci_hi=sorted_rate_samples[idx_ci_hi,:]


#put the rate data out to a file
outfname='{0}-rate-data-bootstrap'.format(name)
output=open(outfname,'w')
#the indexing is in terms of the good milesotnes
for igood,iactual in enumerate(good_milestones):
	m=initial_milestones[iactual]
	halfwidth=np.log10(rate_ci_hi[igood]/rate_ci_lo[igood])/2
	line='{0:d} {1:d} {2:d} {3:.3g} {4:.3g} {5:.3g} {6:.3g} {7:.3f}\n'.format(iactual,m[0],m[1],
		rate_all_data[igood],rate_median[igood],rate_ci_lo[igood],rate_ci_hi[igood],halfwidth)
	output.write(line)   
output.close()

rate_samples_from_starting=rate_samples[:,starting_good_milestone]

nbins=int(np.ceil(2.0*np.sqrt(float(n_successful_bootstrap))))
#there is actually one less "successful" bootstrap than n_successful_bootstrap
#because the "all data" estimate is automatically included
print('total number of bootstraps: ',nbootstrap)
print('fraction successful bootstraps: {0:.3f}'.format((n_successful_bootstrap-1)/nbootstrap))
print('total number of data points, complete, incomplete sims: ',ndata,ndata-nincomplete,nincomplete)
print('total sims in k matrix (not hitting new milestone): ',n_in_k_matrix)
print('number of histogram bins: ',nbins)
print('MFPT estimate from all data (s): {0:.3g}'.format(1.0/rate_all_data[starting_good_milestone]))
print('rate estimate from all data (s^-1): {0:.3g}'.format(rate_all_data[starting_good_milestone]))
print('median MFPT, 68% CI (s): {0:.3g} {1:.3g} {2:.3g}'.format(
	1.0/rate_median[starting_good_milestone],1.0/rate_ci_hi[starting_good_milestone],1.0/rate_ci_lo[starting_good_milestone]))
print('median rate, 68% CI (s^-1): {0:.3g} {1:.3g} {2:.3g}'.format(
	rate_median[starting_good_milestone],rate_ci_lo[starting_good_milestone],rate_ci_hi[starting_good_milestone]))
halfwidth=np.log10(rate_ci_hi[starting_good_milestone]/rate_ci_lo[starting_good_milestone])/2
print('half 68% CI width in orders of magnitude: {0:.3f}'.format(halfwidth))



pmin=np.floor(np.log10(np.nanmin(rate_samples_from_starting)))
pmax=np.ceil(np.log10(np.nanmax(rate_samples_from_starting)))
bins=np.logspace(start=pmin,stop=pmax,num=nbins)
#print(bins)
#try:
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
ax.hist(rate_samples_from_starting,bins=bins)
ax.set_xscale('log')
ax.set_xlabel('rate (s)')
ax.set_ylabel('count of rate values')
ax.axvline(x=rate_all_data[starting_good_milestone],color='k',linestyle='-')
ax.axvline(x=rate_median[starting_good_milestone],color='k',linestyle='--')
ax.axvspan(xmin=rate_ci_lo[starting_good_milestone],xmax=rate_ci_hi[starting_good_milestone],color='k',alpha=0.2)
major_ticks=[10**i for i in range(int(pmin),int(pmax)+1)]
ax.set_xticks(major_ticks)
for fmt in ['pdf','png']:
        outfname='{0}-rate-dist.{1}'.format(name,fmt)
        fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)
plt.close()

fig=plt.figure()
ax=fig.add_subplot(1,1,1)
y=np.linspace(start=0,stop=1,num=n_successful_bootstrap)
ax.plot(np.sort(rate_samples_from_starting),y)
ax.set_xscale('log')
ax.set_xlabel('rate (s)')
ax.set_ylabel('CDF')
ax.axhspan(ymin=0.025,ymax=0.975,color='k',alpha=0.1)
ax.axhspan(ymin=0.5-confidence/2,ymax=0.5+confidence/2,color='k',alpha=0.1)
major_ticks=[10**i for i in range(int(pmin),int(pmax)+1)]
ax.set_xticks(major_ticks)
for fmt in ['pdf','png']:
        outfname='{0}-rate-dist-cdf.{1}'.format(name,fmt)
        fig.savefig(outfname,format=fmt,bbox_inches='tight',dpi=600)
plt.close()

