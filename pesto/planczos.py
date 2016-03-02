#!/usr/bin/python
import numpy
from pio import *
from putil import *
import copy
import random
from dstev.dstev import *
import math
from scipy.optimize import *
from scipy.linalg import *
import pglobals
import pminimise
import pdefect
lanczos_accuracy = 0
counter = 0

			
def optimise_saddle_direct(lattice, min_dir = [], min_tol=1E-2, bfgs = True):
	# uses a the lanczos within the lbfgs to convert saddle into positive definite lattice
	global q, min_mode, eigen_old, counter
	counter = 0
	eigen_old = 1000
	print len(min_dir), lattice.NAtoms*3
	if(len(min_dir) != lattice.NAtoms*3):
		print "Randomly guessing initial dir"
		min_mode = numpy.random.rand(lattice.NAtoms*3)*2 - 1
	else:
		print "Using passed min_mode"
		min_mode = min_dir
		
	print len(min_mode)
	if(bfgs):
		lattice = lanczos_wrapped_bfgs(lattice, min_tol)
	else:
		lattice = lanczos_wrapped_sd(lattice,min_tol)
	
	print "Total force calls", pio.forcecalls
	
	if(eigen_old>0):
		lattice.Failed = True
	
	return lattice, min_mode
	
def reconverge_saddle(lattice):
	global q, min_mode, eigen_old
	eigen_old = 1000
	if(len(min_dir) != lattice.NAtoms*3):
		min_mode = numpy.random.rand(lattice.NAtoms*3)*2 - 1
	else:
		print "Using passed min_mode"
		min_mode = min_dir
		
	print len(min_mode)
	
	bounds = []
	for i in xrange(0, 3*lattice.NAtoms):
		bounds.append([lattice.Pos[i]-0.2, lattice.Pos[i]+0.2])
		
	lattice = lanczos_wrapped_bfgs(lattice, bounds=bounds)
	
	print "Total force calls", pio.forcecalls
	return lattice, min_mode
	
	
def optimise_saddle_decomposed(lattice, min_dir = [], min_tol=1E-4):
	# estimates distance along softest mode, uses constrained lbfgs to minimise perpendicular to mode
	global q, min_mode
	if(len(min_dir) != lattice.NAtoms*3):
		min_mode = numpy.random.rand(lattice.NAtoms*3)*2 - 1
	else:
		min_mode = min_dir
	old_energy = -1E6
	while(math.fabs(lattice.TPE - old_energy)> min_tol):
		old_energy = lattice.TPE
		lattice = eval_lattice(lattice)
		q, min_mode, eigenvals = lanczos_method(lattice, min_mode)
		FdN = numpy.dot(min_mode, lattice.Force)
		if(eigenvals[0] > 0.1):
			lattice.Failed = True
			print "Stumbled into +ve curvature"
			return lattice
		deltax = -FdN/eigenvals[0]
		
		print "DeltaX:", deltax
		if(math.fabs(deltax) > 0.1):
			deltax = 0.1 * numpy.sign(deltax)
		
		print "Displacing:", deltax, "along negative mode", magnitude(min_mode)
		lattice.Pos -= deltax * min_mode
		lattice = constrained_lbfgs(lattice, min_mode, min_tol)
		print "Force is:", lattice.MaxForce
		print "Lanczos delta:", math.fabs(lattice.TPE - old_energy)

	print "Total force calls", pio.forcecalls	
	lattice.min_mode = min_mode
	lattice.Failed = False
	return lattice
	
		

def lanczos_method(lattice, r, min_tol=1E-2,max_steps=20):
	Orig_Force = copy.copy(lattice.Force)	
	Orig_TPE = copy.copy(lattice.TPE)
	Orig_PE = copy.copy(lattice.PE)	
	# generate initial test vector
	
	sigma = magnitude(r)
	q = 0
	#start lanczos iteration
	alpha_vec = []
	sigma_vec = []
	q_vec = []
	steps = 0
	info = 0
	eigen_old = 1000.0
	eigen_diff = 0
	eigen_diff_old = 0
	while(steps < max_steps):
		
		if(steps == 1):
			q_orig = q
		q_old = copy.copy(q)
		q = r/sigma
		if(len(q_vec)>0):
			pass
			nq = numpy.array(q_vec,dtype=float).T
			q = q - numpy.dot(nq, lstsq(nq,q)[0])

		u = hessian_vec_prod(lattice, q, Orig_Force)
		r = u - sigma*q_old
		alpha = numpy.dot(q,r)
		r = r - alpha*q
		alpha_vec.append(alpha)
		q_vec.append(q)
		sigma = magnitude(r)
		sigma_vec.append(sigma)
		
		
		if(steps>1):
			# lets try and get the eigenvalues and vectors using lapacks dstev
			alpha_size = len(alpha_vec)
			workspace = numpy.zeros(alpha_size*2, dtype=float)
			z = numpy.zeros((alpha_size, alpha_size),dtype=float)
			temp_alpha=numpy.array(alpha_vec, dtype=float)		
			temp_sigma =numpy.array(sigma_vec[:-1], dtype=float)		
			pyz = dstev(jobz="V", n=alpha_size, d=temp_alpha, e=temp_sigma, ldz=alpha_size, z=z, work=workspace, info=info)
			eigenvals = temp_alpha
			eigen_diff_old = eigen_diff
			eigen_diff = math.fabs((eigen_old - eigenvals[0])/eigen_old)
				
			print "Eigen_diff:", eigen_diff
			eigen_old = eigenvals[0]
			if(eigen_diff < min_tol):
				break
			elif(eigen_diff>1): steps = 0
		

		steps += 1
	print steps	
	deco = numpy.dot(q,q_orig)	
	print "Decoherence Angle is:", math.fabs(90 - math.degrees(math.acos(deco)))

	nq = numpy.array(q_vec,dtype=float).T

	lattice.TPE = Orig_TPE
	lattice.Force = Orig_Force
	lattice.PE = Orig_PE
	
	min_eig = 0
	for i in xrange(0, 3):
		if(eigenvals[i] < 0 and eigenvals[i] > -1.0): 
			min_eig = i
			break
	
	

	#else:
	print eigenvals
	#print nq.shape, pyz.shape, len(alpha_vec)
	#lattice.Pos = Orig_Pos
	return q, numpy.dot(nq, pyz.T[min_eig]), eigenvals

	

def lanczos_wrapped_sd(lattice, mintol=1E-8):
	Null, lattice.Force = lanczos_inverter(lattice.Pos, lattice, None)
	while(lattice.MagForce>mintol):
		lattice.Pos += normalise(lattice.Force)*0.01
		Null, lattice.Force = lanczos_inverter(lattice.Pos, lattice, None)
	
def lanczos_wrapped_bfgs(lattice, mintol=1E-8, bounds=None):
	global lanczos_accuracy,  min_mode
	print "mintol is:", mintol
	try:
		(lattice.Pos, energy, data) = fmin_l_bfgs_b(lanczos_inverter, lattice.Pos, args=(lattice,None), m=5, maxfun=50, iprint=0, pgtol=mintol, factr=0)# factr=mintol)
		#lattice = sheppards_lbfgs(lattice, tolerance=mintol)
		#energy = lattice.TPE
		print energy
		lattice.eval()
		print "Lanczos Accuracy:", lanczos_accuracy
		if(math.fabs(lanczos_accuracy) > 1E-1):
			lattice.Failed = True
		else:
			lattice.Failed = False
		return lattice	
	except:
		"Saddle convergence failed"
		lattice.Failed = True
		return lattice
	
def lanczos_inverter(Pos, lattice, null):
	global q, min_mode, eigen_old, lanczos_accuracy, counter

	# this is the lancoz approx wrapped in the lbfgs method
	lattice.Pos = Pos
	lattice = eval_lattice(lattice)
	Pos = copy.copy(lattice.Pos)
	Force = copy.copy(lattice.Force)
	PE = copy.copy(lattice.PE)
	#q = numpy.random.rand(lattice.NAtoms*3)*2 - 1
	q, min_mode, eigenvals = lanczos_method(lattice, min_mode)
	lattice.Pos = Pos
	lattice.Force = Force
	lattice.PE = PE
	
	if(eigenvals[0] > 0 and eigen_old > 0):
		print "Positive eigenvalue detected so terminating convergence"
		return 0, 0
	
	eigen_old = eigenvals[0]
	FdI = numpy.dot(lattice.Force, min_mode)
	lattice.Force -= 2*FdI*min_mode
	lattice.min_mode = min_mode
	try:
		sep = separation_mod(lattice, lattice.Start_Pos)
		print "Displacement from activation stage:", sep
		barrier = lattice.TPE - lattice.Start_TPE
		print "Barrier:", barrier
		if(barrier < 1E-2 and sep[0] < 0.5):
			return False
	except: pass
	
	# Assume harmonic, use force and curvature to determine position along 1D function
	# this is the newton method which provides an exact deltax for quadratic functions
	#deltax = FdI*eigenvals[0]
	#deltax = - FdI/eigenvals[0]
	eigenvalue = eigenvals[0]

	deltax = -FdI/eigenvalue
	lanczos_accuracy = FdI

	# then use the approximate quadratic expression, using first and 2nd derivatives
	xenergy = 0.5*eigenvalue*(deltax**2)
	
	#modify the total energy to subtract this modes contribution
	
	print "Magnitude of Force along -ve Eigenvector is:", FdI, xenergy	
	print "Real Energy", lattice.TPE
	# subtract 2* energy to mirror the energy surface
	lattice.TPE -= 2*xenergy
	
	pglobals.status = "Converging with LBFGS wrapped Lanczos. MaxForce is:" + str(lattice.MaxForce())

	return lattice.TPE, -lattice.Force
	

		
def sheppards_lbfgs(lattice, tolerance):
	alpha = 0.05
	tmemory = 25
	offset = 10
	ialpha = 0.1/(alpha+tmemory-offset)
	max_move = 0.2
	# initialize
	s=[]
	y=[]
	rho=[]
	lattice.TPE, lattice.Force = lanczos_inverter(lattice.Pos, lattice, None)
	lattice.Force *= -1
	a = numpy.zeros(tmemory)
	while(lattice.MaxForce()>tolerance):
		Ho = alpha+max(0,(len(s)-offset-1)*ialpha)
		print "Alpha:", Ho
		q = -lattice.Force
		for i in xrange(len(s)-1,-1,-1):
			a[i] = rho[i]*numpy.dot(s[i],q)
			q -= a[i]*y[i]
		
		d = Ho*q
		
		for i in xrange(len(s)):
			b = rho[i] * numpy.dot(y[i],d)
			d += s[i] * (a[i] - b)
			
		d *= -1.0
		
		Pos_Old = lattice.Pos.copy()
		Force_Old = lattice.Force.copy()
		magd = putil.magnitude(d)
		if(magd>max_move):
			lattice.Pos += max_move*d/magd
		else:
			lattice.Pos += d
		lattice.TPE, lattice.Force = lanczos_inverter(lattice.Pos, lattice, None)
		lattice.Force *= -1
		
		print "Energy:", lattice.TPE, "Force:", lattice.MagForce()
		

		s.append(lattice.Pos - Pos_Old)
		# wrap boundaries
		y.append(Force_Old - lattice.Force)
		rho.append(1.0/numpy.vdot(y[-1],s[-1]))
		
		if(len(s)+1>tmemory):
			s.pop(0)
			y.pop(0)
			rho.pop(0)
		
		a2 = putil.magnitude(Force_Old)
		a1 = math.fabs(numpy.vdot(lattice.Force, Force_Old))
		if(a1<0 and a2!=0):
			print "Resetting memory.", a1,a2
			s=[]
			y=[]
			rho=[]
		
			
	return lattice


def hessian_vec_prod(lattice, q, Orig_Force):
	step = 0.001
	#new_lattice = lattice
	lattice.Pos += step*q
	lattice = eval_lattice(lattice)
	lattice.Pos -= step*q
	return -((lattice.Force - Orig_Force)/step)
	print "Total force calls:", io_module.forcecalls
	
	
	
if __name__ == '__main__':
	main()	
			
	
