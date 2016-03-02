"""
pneb is a parallelised nudged elastic band module.
"""

__author__ =  'Louis Vernon <louis.vernon@gmail.com>'


import os,sys,math,copy,random,time, multiprocessing
import numpy, scipy
from pesto import pio, putil
from multiprocessing import Queue

	


class Band():
	pass



def worker(lattice, inputq, outputq):
	while(True):
		work = inputq.get()
		if(not work): return
		lattice.Pos = work[1]
		lattice.eval()
		#print work[0], lattice.TPE
		outputq.put((work[0], lattice.TPE, lattice.Force))
	
		



def run_neb(initial, final, images=11, springk=1.0, minimiser="BFGS", tolerance=1e-3, saddle=False, cpus=0, trajectory=False):
	""" this takes an initial and final lattice configuration spanning a saddle (or series of saddles) and returns an array of lattices converged along the intermediate minimum energy pathway. """
	
	# some initialisation sanity checking
	# first we verify that the initial and final images are relaxed
	
	initial.eval()
	final.eval()
	
	if(initial.MagForce() > 0.1 or final.MagForce() > 0.1):
		print "The initial and final images don't seem to be relaxed."
		return False
	# verify the images span an actual transition
	
	if(putil.separation(initial, final)<0.5):
		print "The initial and final images seem to lie within the same minimum"
		return False
	

	# ok we are good to go, let's spawn the parallel workers
	
	if(cpus==0): cpus = multiprocessing.cpu_count()
	inputq = Queue()
	outputq = Queue()
	
	workers = []
	for i in xrange(cpus):
		workers.append(multiprocessing.Process(target=worker, args=(putil.copylattice(initial),inputq, outputq)))
		workers[i].start()
		
	# now we generate the interpolated band
	band = interpolate_band(initial, final, images, saddle, trajectory)
	if(not band): return False
	band.count = 0
	band.Dim = initial.Dim
	band.PBC = initial.PBC
	band.nimages = images
	band.inputq = inputq
	band.outputq = outputq
	band.springk = springk
	
	band = converge_band(band,minimiser, tolerance)

	for i in xrange(cpus):
		band.inputq.put((False))
	for i in xrange(cpus):
		workers[i].join()
	
	return band
		
def interpolate_band(initial, final, images, saddle, trajectory):
	if(saddle):
		mid_point = int((images-1)/2)
		out_vector = putil.separation_vector(saddle, initial)
		out_vector /= mid_point
		
		in_vector = putil.separation_vector(final, saddle)
		in_vector /= (images-1)-mid_point
		
		
		print "Interpolating through provided midpoint:", mid_point
	elif(trajectory):
		print "Initisialising using predefined trajectory"
		if(len(trajectory) != images):
			print "The trajectory list does not map to the number of images requested."
			return False
	else:
		print "Linearly interpolating between states."
		sep_vector = putil.separation_vector(final, initial)
		sep_vector /= (images-1)
	
	band = Band()
	band.NAtoms = initial.NAtoms
	band.Pos = numpy.zeros(initial.NAtoms*3*images)
	band.Force = numpy.zeros(initial.NAtoms*3*images)
	band.Energy = numpy.zeros(images)
	band.Energy[0] = initial.TPE
	band.Energy[-1] = final.TPE
	band.TPE = 0
	band.Force[0:3*band.NAtoms] = initial.Force
	band.Force[(images-1)*3*band.NAtoms:images*3*band.NAtoms] = final.Force
	if(saddle):
		for i in xrange(mid_point):
			band.Pos[i*3*band.NAtoms:(i+1)*3*band.NAtoms] = initial.Pos + out_vector*i
		for i in xrange(mid_point, images):
			band.Pos[i*3*band.NAtoms:(i+1)*3*band.NAtoms] = saddle.Pos + in_vector*(i - mid_point)
	elif(trajectory):
		for i in xrange(images):
			band.Pos[i*3*band.NAtoms:(i+1)*3*band.NAtoms] = trajectory[i].Pos
	else:	
		for i in xrange(images):
			band.Pos[i*3*band.NAtoms:(i+1)*3*band.NAtoms] = initial.Pos + sep_vector*i
	return band



	
def eval_band(Pos, band, null):
	band.Pos = Pos
	band.count += 1
	for i in xrange(1, band.nimages-1):
		#print "eval image", i
		band.inputq.put((i,band.Pos[i*3*band.NAtoms:(i+1)*3*band.NAtoms]))
		
	for i in xrange(1, band.nimages-1):
		work = band.outputq.get()
		band.Energy[work[0]] = work[1]
		band.Force[work[0]*3*band.NAtoms:(work[0]+1)*3*band.NAtoms] = work[2]
		
	band.Barrier =  numpy.max(band.Energy) - band.Energy[0]
	band.TPE = numpy.sum(band.Energy)
	band = neb_modify_force(band)
	
	print "Barrier:", band.Barrier, "Force:", band.mForce
	return band.mForce, - band.Force
	

def neb_modify_force(band):
	peak_image = numpy.argmax(band.Energy)
	for i in xrange(1, band.nimages-1):
		# "improved" tangent calculation
		Forwards = putil.separation_vector_raw(band.Pos[(i+1)*3*band.NAtoms:(i+2)*3*band.NAtoms], band.Pos[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms], band.PBC, band.Dim)
		Backwards = putil.separation_vector_raw(band.Pos[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms], band.Pos[(i-1)*3*band.NAtoms:(i)*3*band.NAtoms], band.PBC, band.Dim)
		mForwards = putil.magnitude(Forwards)
		mBackwards = putil.magnitude(Backwards)
	
		if(band.Energy[i]>band.Energy[i-1] and band.Energy[i+1] > band.Energy[i]):
			norm = putil.normalise(Forwards)
		elif(band.Energy[i] < band.Energy[i-1] and band.Energy[i+1] < band.Energy[i]):
			norm = putil.normalise(Backwards)
		else:
			norm = putil.normalise(Backwards * math.fabs(band.Energy[i] - band.Energy[i-1]) + Forwards * math.fabs(band.Energy[i+1] - band.Energy[i]))

			
		# if climbing then move uphill independent of the rest of the band
		if(i==peak_image):
			band.Force[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms] -= 2 * norm * numpy.vdot(band.Force[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms], norm)
			
		else:
			# remove parallel component of real force
			band.Force[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms] -= norm * numpy.vdot(band.Force[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms], norm)
			# add in force due to the springs
			band.Force[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms] += norm*(mForwards - mBackwards)*band.springk
			
		#band.mForce = numpy.absolute(band.Force).max()
		band.mForce = putil.magnitude(band.Force)
		band.peak_image = peak_image
		

		
	return band	
		
		
	

def converge_band(band,minimiser, tolerance):
	
	if(minimiser=="SLBFGS"):
		(band.Pos, energy, data) = scipy.optimize.fmin_l_bfgs_b(eval_band, band.Pos, args=(band,1), m=1, factr=0, pgtol=tolerance, iprint=1)
	elif(minimiser=="LBFGS"):
		band = lbfgs(band, tolerance)
	else:
		print "Unknown NEB optimizer", minimiser
	
	return band

def lbfgs(band, tolerance=1e-2, max_move = 0.2, memory = 5, alpha = 0.01, scale=0):
	#Ho = alpha
	ialpha = scale*alpha/(memory)
	s=[]
	y=[]
	rho=[]

	a = numpy.zeros(memory)
	max_move *= math.sqrt(band.nimages)

	work = 0
	# initialize
	s=[]
	y=[]
	rho=[]
	band.mForce, band.Force = eval_band(band.Pos, band, None)
	band.Force *= -1
	a = numpy.zeros(memory)
	while(band.mForce>tolerance):
		Ho = alpha+max(0,(len(s)-1)*ialpha)
		q = -band.Force
		for i in xrange(len(s)-1,-1,-1):
			a[i] = rho[i]*numpy.dot(s[i],q)
			q -= a[i]*y[i]
		
		d = Ho*q
		
		for i in xrange(len(s)):
			b = rho[i] * numpy.dot(y[i],d)
			d += s[i] * (a[i] - b)
			
		d *= -1.0
		
		Pos_Old = band.Pos.copy()
		Force_Old = band.Force.copy()
		magd = putil.magnitude(d)
		if(magd>max_move):
			band.Pos += max_move*d/magd
		else:
			band.Pos += d
		band.mForce, band.Force = eval_band(band.Pos, band, None)
		band.Force = -band.Force
		#print "Energy:", band.TPE, "Force:", band.MagForce()
		

		s.append(band.Pos - Pos_Old)
		# wrap boundaries
		y.append(Force_Old - band.Force)
		rho.append(1.0/numpy.vdot(y[-1],s[-1]))
		
		if(len(s)+1>memory):
			s.pop(0)
			y.pop(0)
			rho.pop(0)
			
		a1 = numpy.vdot(band.Force, Force_Old)
		if(a1<-0.5*putil.magnitude(Force_Old)):
			print "Resetting memory.", a1
			s=[]
			y=[]
			rho=[]	
		work+=1
	print "Evaluations:", work
	return band
		
	


		
		
				
	
	
	
	
	
	
	
	
	
