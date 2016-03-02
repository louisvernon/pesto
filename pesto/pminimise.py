#!/usr/bin/python
import os,sys,math,random,copy,time
import numpy
from numpy import dot
from math import sqrt
import pio,putil
import pglobals
from scipy.optimize import *

# test
exclusion=0
dv =[]
dv_total=0
ref=0
variable_tol = True

def minimise_lattice(lattice,mintype,mintol):
	global exclusion,dv,dv_total,ref
	exitstatus =0
	totalsteps =0
	funcevalcount = 0
	step = 0.05


	ref = lattice
	lattice = lattice.eval()
	old_lattice = putil.copylattice(lattice)
	orig_energy = old_lattice.TPE
	
	putil.printl("Energy:", orig_energy, "Max Force",lattice.MaxForce(), "FuncEvalCount", pio.forcecalls)

	if(lattice.MaxForce() < 0.002 ):
		putil.printl("We've started too close to the minimum (or saddle point) quitting")
		exitstatus =1;
		return lattice
	
	if(mintype == "SD"):
		echange_old = 1	
		putil.printl("MINIMISING USING SD")	
		echange = 2
		while((mintol<1 and lattice.MaxForce() > mintol) or (mintol>=1 and totalsteps < mintol)):

			search_direction = putil.normalise(lattice.Force)
			prevV = lattice.TotalPE()
			lattice.Pos += step*search_direction
			lattice = lattice.eval()					
			echange =  prevV - lattice.TPE

			if(numpy.dot(lattice.Force,old_lattice.Force)<0):
				step /= 2
			else:
				step *= 1.5
				
			
			
			# new bit of code, so long as the stepsize is increasing we are far away from the minimum				
			if(echange > echange_old and echange_old > 1e-5):
				totalsteps-=1
			echange_old = echange
			old_lattice.Force = copy.copy(lattice.Force)
			totalsteps+=1


	
	elif(mintype == "CG"):
		putil.printl("MINIMISING USING CG")
		echange = 10
		search_direction = lattice.Force
		method = 1
		fstarttime = time.time()
		cg_count = 0
	

		while(lattice.MaxForce() > mintol or cg_count < 5): 
				
			cg_count += 1
			starttime = time.time()
			prevV = lattice.TotalPE()		
			search_direction = putil.normalise(search_direction)		
			fstart = copy.copy(lattice.Force)
			#calculate first step length in line minimiser
			#calculate_first_step(lattice,search_direction,exclusion,dv,dv_total,ref)
			#or fixed step length
			step = 0.01
			
			#use step length in line minimiser
			
			
			lattice = secant_line_min(lattice, search_direction, step, mintol*10)
			#lattice = quadratic_line_min(lattice, search_direction, step, 0.005)
			
			if(method==1):	
				#weird version in tsse
				
				yk = lattice.Force - fstart #determine the new search direction based 
				beta = max(0, dot(yk, lattice.Force)/dot(fstart, fstart)) #on scipy's method
				search_direction = lattice.Force + beta * search_direction  
			else:	
				#polak-rib...
				
				gamma = dot((lattice.Force-fstart),lattice.Force)
				fstartmag  = putil.magnitude(fstart)
				gamma = gamma/(fstartmag**2)
				#print "GAMMA",gamma
				search_direction = lattice.Force + gamma*search_direction
				
			echange =  prevV -lattice.TotalPE()


			totalsteps+=1
	
	elif(mintype=="SLBFGS"):
		print "Attempting to use SCIPYs bundled BFGS mimiser"
		(lattice.Pos, energy, data) = fmin_l_bfgs_b(putil.scipy_eval_lattice, lattice.Pos, args=(lattice,1), m=10, factr=0, pgtol=mintol, iprint=0)
		print energy
	
	elif(mintype=="LBFGS"):
		print "Using Local LBFGS implementation"
		lattice = lbfgs(lattice, mintol)
	
	
	elif(mintype=="Quickmin"):
		dt = 0.1
		dR_Max = 1E8
		V = numpy.zeros(3*lattice.NAtoms)
		while(lattice.MaxForce() > mintol):
			VdotF = numpy.vdot(V,putil.normalise(lattice.Force))
			if(VdotF > 0):
				V = putil.normalise(lattice.Force) * VdotF
			else:
				V = numpy.zeros(3*lattice.NAtoms)
			dR = V * dt
			V += lattice.Force * dt
			dR += V * dt
			dR /= 2
			if(putil.magnitude(dR)>dR_Max):
				dR = putil.normalise(dR) * dR_Max
			
			lattice.Pos += dR
			lattice.eval()
			totalsteps += 1
			
				
	else:
		print "You haven't selected a valid minimiser"
		sys.exit()
			
	
	
        print "Energy change", lattice.TPE-orig_energy
	return lattice	



def calculate_first_step(lattice,step,search_direction,exclusion,dv,dv_total,ref):	
	step = 0.001
	F = dot(lattice.Force,search_direction)
	lattice.Pos += step * search_direction #take the step
					
	lattice = lattice.eval()

	Vold = lattice.TotalPE()

	F2 = dot(lattice.Force,search_direction)
	
	printd("INITIAL STEP", step)
	#Calculate step length using :  search_dir.force(r)/( search_dir.force(r+step) - search_dir.force(r) )
	step = -step*(F/(F2-F))
	printd("CALCULATED STEP", step)


def quadratic_line_min(lattice, direction, step, tol):
	#if(pio.verbose):putil.printl("USING REGULA FALSI LINEMIN")
	a_lattice = copy.copy(lattice)
	b_lattice = copy.copy(lattice)
	c_lattice = copy.copy(lattice)	
	a = 0
	b = step
	dfa = numpy.dot(a_lattice.Force, direction)
	b_lattice.Pos = lattice.Pos + (direction * b)
	b_lattice = b_lattice.eval()
	
	dfb = numpy.dot(b_lattice.Force, direction)

	loopcount=0
	dfc = tol+1
	while(abs(dfc) > tol and loopcount < 50):
		c = (dfa*b - dfb*a) / (dfa - dfb)
		c_lattice.Pos = lattice.Pos + (direction * c)
		c_lattice = c_lattice.eval()
		dfc = numpy.dot(c_lattice.Force, direction)
		#if(pio.verbose): print a, b, c, dfa, dfb, dfc
		if(dfa*dfc>0):
			a = c
			dfa = dfc
		else:
			b = c
			dfb = dfc	
	return c_lattice
	
def secant_line_min(lattice, direction, step, tol):
	#if(pio.verbose):putil.printl("USING SECANT LINEMIN")
	a_lattice = copy.copy(lattice)
	b_lattice = copy.copy(lattice)
	c_lattice = copy.copy(lattice)	
	a = 0
	b = step
	dfa = numpy.dot(a_lattice.Force, direction)
	b_lattice.Pos = lattice.Pos + (direction * b)
	b_lattice = b_lattice.eval()
	dfb = numpy.dot(b_lattice.Force, direction)
	#while(dfa*dfb >0):
	#	#if(pio.verbose):print "NOT BRACKETED",step,dfb,dfa
	#	bold = b
	#	b = ((dfa*b - dfb*a) / (dfa - dfb))*1.5
	#	a = bold
	#	dfa=dfb
	#	b_lattice.Pos = lattice.Pos + (direction * b)
	#	b_lattice = eval_lattice(b_lattice)
	#	dfb = numpy.dot(b_lattice.Force, direction)
		
	#if(pio.verbose):print "BRACKETED",a,b,dfb,dfa	
	loopcount=0
	dfc = tol+1
	while(abs(dfc) > tol and loopcount < 50):
		c = b - ((b-a)/(dfb - dfa))*dfb
		c_lattice.Pos = lattice.Pos + (direction * c)
		c_lattice = c_lattice.eval()
		dfc = numpy.dot(c_lattice.Force, direction)
		#if(pio.verbose): print a, b, c, dfa, dfb, dfc
		a = b
		dfa = dfb
		b = c
		dfb = dfc
		
	return c_lattice

def invquad_line_min(lattice, direction, step, tol):
	#if(pio.verbose):putil.printl("USING INVERSE QUADRATIC LINEMIN")
	a_lattice = copy.copy(lattice)
	b_lattice = copy.copy(lattice)
	c_lattice = copy.copy(lattice)
	d_lattice = copy.copy(lattice)	
	a = 0
	b = a+step/2.0
	c = step
	dfa = numpy.dot(a_lattice.Force, direction)
	b_lattice.Pos = lattice.Pos + (direction * b)
	b_lattice = b_lattice.eval()
	
	dfb = numpy.dot(b_lattice.Force, direction)
	c_lattice.Pos = lattice.Pos + (direction * c)
	c_lattice = c_lattice.eval()
	dfc = numpy.dot(c_lattice.Force, direction)
	
	
	while(dfa*dfc >0):
		#if(pio.verbose):print "NOT BRACKETED",step,dfa,dfb
		cold=c
		c = ((dfa*c - dfc*a) / (dfa - dfc))*1.5
		a = cold
		dfa=dfc
		c_lattice.Pos = lattice.Pos + (direction * c)
		c_lattice = c_lattice.eval()
		dfc = numpy.dot(c_lattice.Force, direction)
		cold=c
	#if(pio.verbose):print "BRACKETED",a,b,c,dfa,dfb,dfc
	
	
	loopcount=0
	dfd = tol+1
	
	while(abs(dfd) > tol and loopcount < 50):
		R = dfb/dfc
		S = dfb/dfa
		T = dfa/dfc
		P = S*(T*(R-T)*(c-b)-(1-R)*(b-a))
		Q = (T-1)*(R-1)*(S-1)
		d = b + P/Q
		d_lattice.Pos = lattice.Pos + (direction * d)
		d_lattice = d_lattice.eval()
		dfd = numpy.dot(d_lattice.Force, direction)
		#if(pio.verbose): print a, b, c, dfa, dfb, dfc
		a = b
		dfa = dfb
		b = c
		dfb = dfc
		c = d
		dfc = dfd		
	return d_lattice		
		
		
		
def lbfgs(lattice, tolerance=1e-2, max_move = 0.2, memory = 50, alpha = 0.05, scale=3):
	Ho = alpha
	ialpha = scale*alpha/(memory)
	s=[]
	y=[]
	rho=[]
	lattice.eval()
	a = numpy.zeros(memory)

	while(lattice.MaxForce()>tolerance):
		Ho = alpha+max(0,(len(s)-1)*ialpha)
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
		lattice.eval()
		print "Energy:", lattice.TPE, "Force:", lattice.MaxForce()
		

		s.append(lattice.Pos - Pos_Old)
		# wrap boundaries
		y.append(Force_Old - lattice.Force)
		rho.append(1.0/numpy.vdot(y[-1],s[-1]))
		
		if(len(s)+1>memory):
			s.pop(0)
			y.pop(0)
			rho.pop(0)
		

		a1 = numpy.vdot(lattice.Force, Force_Old)
		if(a1<-0.5*putil.magnitude(Force_Old)):
			print "Resetting memory.", a1
			s=[]
			y=[]
			rho=[]
			
	return lattice
		

					
