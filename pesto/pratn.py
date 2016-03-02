#!/usr/bin/python

"""


"""


import numpy
from pio import *
from putil import *
import copy
import random
from pesto.dstev.dstev import *
import math
from scipy.optimize import *
from scipy.linalg import *
import planczos
import pdefect
import pminimise
import pneb


def failed_return(lattice):
	lattice.Failed = True
	if(initial_lattice.local):
		write_lattice(lattice, "failed.dat")
	return lattice


def ratn_method(initial_lattice, translation_step=0.1, min_eig = -0.1, activation_length = 0.5, mode=1, min_tol=1E-3, persistant = False, find_intermediate=True):
	"""	# RAT Modes:
	# 0 - Activate until ridge crossed using RAT, converge using the Lanczos
	# 1 - Activate to negative curvature using RAT, converge using the Lanczos
	# 2 - Activate using Lanczos to climb, converge once in negative curvature using inverted Lanczos 
	"""
	try:
		if(initial_lattice.local):
			pass
	except Exception:
		initial_lattice.local=False
		
		
	# permute step size and activation length a little to avoid convergence bias
	
	translation_step += numpy.random.random() * (translation_step/2) - translation_step/4
	activation_length += numpy.random.random() * (activation_length/2) - activation_length/4
	
	print "Translation_Step:", translation_step, "Activation_Length:", activation_length
	
	lattice = copylattice(initial_lattice)

	print "RATN mobile atoms:", numpy.sum(lattice.DefectMask)/3
	lattice.ForceCalls = pio.forcecalls

	# I've got probably unfounded concerns about the random seed during MPI runs. I'm going to get numpy to seed from /dev/urandom just to be sure
	numpy.random.seed()
	# get initial energy
	lattice = eval_lattice(lattice)
	orig_lattice = copylattice(lattice)
	orig_energy = lattice.TotalPE()
	old_lattice = copylattice(lattice)
	
	# here I am returning a copy of the original lattice position in order to let the server verify the saddle is synchronised with the kmc
	lattice.Start_Pos = copy.copy(orig_lattice.Pos)
	lattice.Start_TPE = lattice.TPE
	
	
	#generate random initial displacement vector
	try: n = lattice.n * lattice.DefectMask
	except: 
		n = numpy.random.normal(0, 1, lattice.NAtoms*3)
		n *= lattice.DefectMask
#		print n
	
	#alpha controls the mixing between the relaxed displacement vector and the previous translation vector
	alpha = 0.1
	
	
	# generate random initial lanczos vector
	min_mode = normalise(numpy.random.rand(lattice.NAtoms*3)*2 - 1)

	# normalise
	n = normalise(n)

	DefectMask = copy.copy(lattice.DefectMask)
	
	print "Translation step for:", numpy.sum(DefectMask)/3, "atoms is:", translation_step
	print "Activation length is:", activation_length
	
	print separation_mod(lattice, lattice.Start_Pos)
	eigen_old = 1000
	

	count = 0
	while(1):

		count += 1
		if(count > 500):
			print "Exiting search early due to excessive number of internal steps"
			failed_return(lattice)

		old_lattice.Pos = copy.copy(lattice.Pos)
		old_lattice.TPE = lattice.TPE		
		disp_old = separation_mod(lattice, lattice.Start_Pos)
		#print magnitude(n)

		lattice.Pos += (n*translation_step)
		lattice.eval()

		disp_temp = separation_mod(lattice, lattice.Start_Pos)
		
		forward_energy_diff = lattice.TPE - old_lattice.TPE
		barrier = lattice.TPE - orig_lattice.TPE

		print "Energy change", barrier

		print "Displacement from origin:", count, disp_old

		disp = separation_mod(lattice, lattice.Start_Pos)
		
		curvature = numpy.dot(-n, lattice.Force)
		print "Curvature:", curvature
		if(disp[2] > activation_length):	
			dot_direction = numpy.dot((lattice.Pos-lattice.Start_Pos), n)
			print "Dot direction, forward energy diff", dot_direction, forward_energy_diff
			if(mode > 0):
				lattice.DefectMask = numpy.ones(3*lattice.NAtoms,dtype=bool)
				lattice.eval()

				q, min_mode, eigenvals = planczos.lanczos_method(lattice, n, 1E-3, 8)


				if(mode==2):		
					n = normalise(min_mode) 
					if(numpy.dot(min_mode, lattice.Force) > 0):
						n *= -1.0

				lattice.DefectMask = copy.copy(DefectMask)
				
				#if(eigenvals[0] < eigen_old and eigen_old < min_eig):
				if(eigenvals[0] < min_eig):				 
					break
				eigen_old = eigenvals[0]
			else:
				print dot_direction, forward_energy_diff
				if(dot_direction > 0 and (forward_energy_diff < 0 or curvature < 0)):

					lattice.DefectMask = numpy.ones(3*lattice.NAtoms,dtype=bool)
					q, min_mode, eigenvals = planczos.lanczos_method(lattice, n, 1E-3, 8)
				
					lattice.DefectMask = copy.copy(DefectMask)
					
					print "Eigenvalue is:", eigenvals[0]
					if(eigenvals[0] < eigen_old and eigen_old < min_eig):				 
						break
					eigen_old = eigenvals[0]
						
						
		
		# the RAT part, we construct a new displacement vector bases on n_old and 
		if(mode<2):
			lattice = constrained_sd(lattice, n, step_orig=0.1)
			
			n_prime = normalise(separation_vector(lattice,old_lattice.Pos))
			n = normalise((n+n_prime*alpha)*lattice.DefectMask)
			
			
		#sign = numpy.sign(numpy.dot(n, normalise(lattice.Pos-orig_lattice.Pos)))

		
		
	print "Displacement from activation stage:", separation_mod(lattice, orig_lattice)
	print "Energy barrier:", lattice.TPE - orig_lattice.TPE
	lattice.Barrier = lattice.TPE - orig_lattice.TPE

		
	Activated = copylattice(lattice)
			
	lattice.DefectMask[:] = True
	#lattice.eval()
	# clear the defectmask so we can converge on the correct saddle, will the negative eigenvalue remain?
	lattice, min_mode = planczos.optimise_saddle_direct(lattice, min_mode, min_tol)

	
	lattice.Saddle = copy.copy(lattice.Pos)


	print lattice.ForceCalls, pio.forcecalls
	
	lattice.ForceCalls = pio.forcecalls - lattice.ForceCalls
	print lattice.ForceCalls
	lattice.Barrier = lattice.TPE - orig_lattice.TPE
		
	print "Converged? Saddle barrier is:", lattice.Barrier

	sep_before = separation_mod(lattice, old_lattice)
	sep_after = separation_mod(lattice, Activated)
	sep_orig = separation_mod(lattice, orig_lattice)
	print "Distance between activated and converged states is:", sep_after, sep_before
	print "Distance from the origin is:", sep_orig
	

	# check we've found something that could be a real saddle
	if(sep_orig[0] < 0.5):
		print "Converged too close to origin"
		failed_return(lattice)




	# saddle check verifies we are bracketing a saddle
	lattice = saddle_check(lattice)
	if(find_intermediate):
		end_lattice = copylattice(lattice)
		end_lattice.Pos = lattice.Approx_End_Pos
		end_lattice = pminimise.minimise_lattice(lattice, "LBFGS", 1E-3)
		displacement = separation_mod(end_lattice, old_lattice)
		if(displacement[0]>1):
			band = pneb.run_neb(old_lattice, end_lattice, images=9)
	
	
	lattice = local_saddle_check(lattice)
	
		
	if(lattice.Failed == True):
		failed_return(lattice)
	else:
		print "Search completed succesfully."
	

	return lattice

