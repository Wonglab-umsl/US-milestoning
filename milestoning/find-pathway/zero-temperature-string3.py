#!/usr/bin/env python
import sys
import numpy as np
import math
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 22})
import matplotlib.pyplot as plt
import pdb
#import fpectl

#from scipy import interpolate
#Note: this requires the 1.9.0 version of scipy.  Earlier versions do not have the cubic interpolator,
#and version 1.9.1 gives an error on trying to import scipy.special.
#To install, use pip install scipy==1.9.0 --user
from scipy.interpolate import RegularGridInterpolator as rgi
from scipy.optimize import minimize_scalar #we will need this later
from scipy.interpolate import CubicSpline

#read the WHAM info file and construct the grid
def read_and_construct_interpolators(infname,res,maxfe):
	grid_points=[]
	grid_values=[]
	input=open(infname,'r')
	for line in input:
		#coordinates of the bin center
		words=line.split()
		bincenter=(float(words[0]),float(words[1]),float(words[2]))
		grid_points.append(bincenter)
		fe=float(words[-1])
		grid_values.append(fe)
	input.close()
	grid_points=np.array(grid_points)
	grid_values=np.array(grid_values)
	imin=np.argmin(grid_values)
	min_coords=grid_points[imin]
	#better to fill in points with the highest free energy -- this produces fewer numerical errors on interpolation
	if (maxfe is None):
		maxfe=np.max(grid_values)
		print("maxfe: ",maxfe)
	npoints=grid_points.shape[0]
	print("total number of grid values: ",npoints)
	print("global minimum: ",min_coords,grid_values[imin])
	bbox_min=np.min(grid_points,axis=0)
	bbox_max=np.max(grid_points,axis=0)
	print("bounding box minimum: ",bbox_min)
	print("bounding box maximum: ",bbox_max)
	bounds=[(bbox_min[k],bbox_max[k]) for k in range(0,3)]
	dims=(bbox_max[:]-bbox_min[:])/res
	dims=np.ndarray.astype(dims,np.int32)+1
	dims=tuple(dims)
	print("number of points: ",dims)
	grid=np.full(dims,maxfe,dtype=np.float64)
	for i in range(0,npoints):
		idx=(grid_points[i,:]-bbox_min[:])/res
		idx=tuple(np.ndarray.astype(idx,np.int32))
		#print(idx,grid_values[i])
		grid[idx]=grid_values[i]
	pointlist=[]
	grad_pointlist=[]
	for k in range(0,3):
		#points are bbox_min[k], bbox_min[k]+res, bbox_min[k]+2*res, etc.
		l=np.linspace(start=bbox_min[k],stop=bbox_max[k],num=dims[k])
		#points for gradient grid, halfway in between
		l2=np.linspace(start=bbox_min[k]+0.5*res,stop=bbox_max[k]-0.5*res,num=dims[k]-1)
		#print(l)
		pointlist.append(l)
		grad_pointlist.append(l2)
		#print(l2)
	pointlist=tuple(pointlist)
	interp=rgi(points=pointlist,values=grid,method='cubic',bounds_error=False,fill_value=maxfe)
	#construct separate interpolators for the gradient
	#they are on points halfway in between those for the main grid, 
	#see "l2" above
	grad_grid_x=np.full((dims[0]-1,dims[1],dims[2]),np.nan,dtype=np.float64)
	grad_grid_x[:,:,:]=(grid[1:,:,:]-grid[:-1,:,:])/res
	grad_interp_x=rgi(points=[grad_pointlist[0],pointlist[1],pointlist[2]],values=grad_grid_x,
		method='linear',bounds_error=False,fill_value=0)
	grad_grid_y=np.full((dims[0],dims[1]-1,dims[2]),np.nan,dtype=np.float64)
	grad_grid_y[:,:,:]=(grid[:,1:,:]-grid[:,:-1,:])/res
	grad_interp_y=rgi(points=[pointlist[0],grad_pointlist[1],pointlist[2]],values=grad_grid_y,
		method='linear',bounds_error=False,fill_value=0)
	grad_grid_z=np.full((dims[0],dims[1],dims[2]-1),np.nan,dtype=np.float64)
	grad_grid_z[:,:,:]=(grid[:,:,1:]-grid[:,:,:-1])/res
	grad_interp_z=rgi(points=[pointlist[0],pointlist[1],grad_pointlist[2]],values=grad_grid_z,
		method='linear',bounds_error=False,fill_value=0)
	grad_interp=[grad_interp_x,grad_interp_y,grad_interp_z]
	return min_coords, interp, grad_interp,maxfe



def total_energy(interp,maxfe,krepel,rrepel,n,x):
	fe=0.0
	for i in range(0,n):
		en=interp(x[i,:])[0]
		fe+=en
	diff=np.zeros(3,dtype=np.float64)
	aux=np.zeros(3,dtype=np.float64)
	d2=rrepel*rrepel
	erepel=0.0
	#pdb.set_trace()
	for i in range(0,npoints): 
		for j in range(i+krepel,npoints):
			diff[:]=points[j,:]-points[i,:]
			dist2=np.dot(diff,diff)
			if (dist2<d2):
				#pdb.set_trace()
				erepel+=maxfe*(1-dist2/d2)*(1-dist2/d2)
	etot=fe+erepel
	return etot, fe, erepel


def print_points(interp,n,x,origin,spline):
	#print(x)
	etot=0.0
	en=interp(origin)[0]
	line='{0:4d} {1:6.3f} {2:6.3f} {3:6.3f} {4:6.3f}'.format(1,origin[0],origin[1],origin[2],en)
	print(line,flush=True)
	for i in range(0,n):
		en=interp(x[i,:])[0]
		etot+=en
		if (i==0):
			diff=x[i,:]-origin[:]
		else:
			diff=x[i,:]-x[i-1,:]
		l=np.sqrt(np.dot(diff,diff))
		#if (i<(n-1)):
			#diff=x[i+1,:]-x[i,:]
			#l=np.sqrt(np.dot(diff,diff))
		grad=np.zeros(3,dtype=np.float64)
		for k in range(0,3):
			grad[k]=grad_interp[k](x[i-1,:])[0]
		mgrad=np.sqrt(np.dot(grad,grad))
		if (spline is not None):
			t=spline(i+1,1)
			mt=np.sqrt(np.dot(t,t))
			#calculate the cosine of the angle between tangent and gradient
			if (mgrad>0.0):
				ctheta=np.dot(t,grad)/(mt*mgrad)
			else:
				ctheta=0.0
			if (ctheta<1.0):
				theta=math.acos(ctheta)*(180.0/math.pi)
			else:
				theta=0.0
			#print(i+1,x[3*i],x[3*i+1],x[3*i+2],en)
			#the origin is anchor #1; the first movable point is anchor 2
			line='{0:4d} {1:12.6f} {2:12.6f} {3:12.6f} {4:6.3f} {5:6.3f} {6:6.3f} {7:6.3f}'.format(
				i+2,x[i,0],x[i,1],x[i,2],en,l,mgrad,theta)
		else:
			line='{0:4d} {1:12.6f} {2:12.6f} {3:12.6f} {4:6.3f} {5:6.3f}'.format(
				i+2,x[i,0],x[i,1],x[i,2],en,l)
		print(line,flush=True)


def calc_grad(interp,grad_interp,maxfe,krepel,rrepel,npoints,points,grad):
	#pdb.set_trace()
	for i in range(0,npoints):
		for k in range(0,3):
			grad[i,k]=grad_interp[k](points[i,:])[0]
   	#self repulsion term: sum_{ij} maxfe*(1-(d_ij/d)^2)
   	#gradient: -maxfe*(2/d^2)*(x_j-x_i)
	diff=np.zeros(3,dtype=np.float64)
	aux=np.zeros(3,dtype=np.float64)
	d2=rrepel*rrepel
	for i in range(0,npoints): 
		for j in range(i+krepel,npoints):
			diff[:]=points[j,:]-points[i,:]
			dist2=np.dot(diff,diff)
			if (dist2<d2):
				#pdb.set_trace()
				grad[j,:]-=(4.0*maxfe/d2)*(1-dist2/d2)*diff[:]
				grad[i,:]+=(4.0*maxfe/d2)*(1-dist2/d2)*diff[:]
	#return grad

def calc_dx(interp,grad_interp,maxfe,krepel,rrepel,npoints,points,dt,dx):
	#fourth-order r-k method
	#pdb.set_trace()
	grad=np.zeros((4,npoints,3),dtype=np.float64)
	p=np.zeros((3,npoints,3),dtype=np.float64)
	#dx=np.zeros((n,3),dtype=np.float64)
	#grad[0,:,:]=
	calc_grad(interp,grad_interp,maxfe,krepel,rrepel,npoints,points,grad[0,:,:])
	p[0,:,:]=points[:,:]-0.5*dt*grad[0,:,:]  #for k(1), minus sign is because we have x'=-del(U)
	#grad[1,:,:]=
	calc_grad(interp,grad_interp,maxfe,krepel,rrepel,npoints,p[0,:,:],grad[1,:,:])
	p[1,:,:]=points[:,:]-0.5*dt*grad[1,:,:]
	#grad[3,:,:]=
	calc_grad(interp,grad_interp,maxfe,krepel,rrepel,npoints,p[1,:,:],grad[2,:,:])
	p[2,:,:]=points[:,:]-dt*grad[2,:,:]
	#grad[4,:,:]=
	calc_grad(interp,grad_interp,maxfe,krepel,rrepel,npoints,p[2,:,:],grad[3,:,:])
	dx[:,:]=-dt*(grad[0,:,:]+2*grad[1,:,:]+2*grad[2,:,:]+grad[3,:,:])/6.0

#respace the points so that they are spaced equal to d
#the spline should be set up so that 0 = origin, etc.
#the positions will be offset by 1 from indices into the points array
def respace_points(npoints,origin,d,points,newpoints):
	#positions along spline where the new points will be constrcuted
	s=np.full(npoints+1,np.nan,dtype=np.float64)
	s[0]=0
	#new_points=np.full((npoints,3),np.nan,dtype=np.float64)
	prev=np.full(3,np.nan,dtype=np.float64)
	diff=np.full(3,np.nan,dtype=np.float64)
	for i in range(0,npoints):
		if (i==0):
			prev[:]=origin[:]
		else:
			prev[:]=points[i-1,:]
		diff[:]=points[i,:]-prev[:]
		s[i+1]=s[i]+np.sqrt(np.dot(diff,diff))
	#this is a slight variation from the procedure in the zts paper
	s[:]=s[:]/d
	spoints=np.vstack([origin,new_points])
	spline=CubicSpline(x=s,y=spoints,axis=0,
		bc_type='natural',extrapolate=True)
	for i in range(1,npoints+1):
		newpoints[i-1,:]=spline(i)
	return spline #print_points needs it for tangent vectors


infname=sys.argv[1]
res=float(sys.argv[2])
print('reading free energy surface from ',infname,flush=True)
print('table resolution ',res,' A',flush=True)
origin, interp, grad_interp, maxfe=read_and_construct_interpolators(infname,res,maxfe=None)
#pdb.set_trace()
#print(grad_interp)
#sys.exit(-1)
npoints=int(sys.argv[3])
length=float(sys.argv[4])
guess_direction=[]
for iarg in range(5,8):
	guess_direction.append(float(sys.argv[iarg]))
guess_direction=np.array(guess_direction)
guess_direction=guess_direction/np.sqrt(np.dot(guess_direction,guess_direction))
print("guess direction: ",guess_direction)
#parameters for repulsive potential: only apply to positions that are at least k apart, and use r as the repulsion radius
krepel=int(sys.argv[8])
rrepel=float(sys.argv[9])
print("parameters for repulsive potential: k, r = ",krepel,rrepel)
orig_stepsize=float(sys.argv[10])
stepsize=orig_stepsize
maxiter=int(sys.argv[11])
#construct the initial pathway
init_coords=np.full((npoints,3),np.nan,dtype=np.float64)
#tangent=np.full((npoints,3),np.nan,dtype=np.float64)

d=length/npoints
print("number of points, overall length, distance between points: ",npoints,length,d)
for i in range(0,npoints):
	init_coords[i,:]=origin[:]+float(i+1)*d*guess_direction[:]
	#tangent[i,:]=guess_direction[:]
#print(init_coords)
print("step size, max iterations: ",stepsize,maxiter)
print("initial coordinates: ")
print_points(interp,npoints,init_coords,origin,None)

iter=1
points=np.copy(init_coords)
new_points=np.full((npoints,3),np.nan,dtype=np.float64)
new_points2=np.full((npoints,3),np.nan,dtype=np.float64)
#print_points(interp,npoints,x)
#pdb.set_trace()
etot=total_energy(interp,maxfe,krepel,rrepel,npoints,points)
print("total energy: ",etot)
dx=np.full((npoints,3),np.nan,dtype=np.float64)
diff=np.full(3,np.nan,dtype=np.float64)
dev=np.full(npoints,np.nan,dtype=np.float64)
#new_tangent=np.full((npoints,3),np.nan,dtype=np.float64)
tol=0.001
spline=None
while (iter<=maxiter):
	#calculate optimizing adjustments and move the points
	calc_dx(interp,grad_interp,maxfe,krepel,rrepel,npoints,points,stepsize,dx)
	new_points[:,:]=points[:,:]+dx[:,:]
	#respace the points
	del spline
	spline=respace_points(npoints,origin,d,new_points,new_points2)
	#print(spline)
	#newe=total_energy(interp,npoints,new_points)
	#print("newe: ",newe)
	#check deviation
	for i in range(0,npoints):
		diff[:]=new_points2[i,:]-points[i,:]
		dev[i]=np.sqrt(np.dot(diff,diff))
	maxdev=np.max(dev)
	points[:,:]=new_points2[:,:]

	if ((iter%10==0) or (maxdev<tol)):
		etot, fe, erepel=total_energy(interp,maxfe,krepel,rrepel,npoints,points)
		print("iteration {0:d} pathway energy {1:.3f} repulsion energy {2:.3f} total energy {3:3f} deviation {4:.6f}".format(
			iter,fe,erepel,etot,maxdev),flush=True)
		if (maxdev<tol):
			break
		#print_points(interp,npoints,x)
	iter=iter+1

final_coords=np.copy(points)
print("final coordinates: ")
print_points(interp,npoints,final_coords,origin,spline)
sys.exit(0)

