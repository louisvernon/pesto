import os,sys,math
import numpy
import random
import time
import pglobals
import pio
import copy
import inspect
import timeit
try: from collections import defaultdict
except: pass

from scipy.optimize import *
from scipy.linalg import *
import pminimise


def printl(*args):
    if(pio.verbose):
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        string= mod.__name__+" : "
    
        for arg in args:
            string+=str(arg)+" "    
     
        print string
 
def printd(*args):
    if(pio.debug):
        frm = inspect.stack()[1]
        mod = inspect.getmodule(frm[0])
        string= mod.__name__+" : "
    
        for arg in args:
            string+=str(arg)+" "    
     
        print string        


def check_numerical_forces(lattice,h):
    HVec = h*normalise(randvec(lattice.Pos))
    
    forwardlattice = copylattice(lattice)
    backlattice = copylattice(lattice)
    
    forwardlattice.Pos +=HVec
    backlattice.Pos -=HVec
    
    forwardlattice = pio.eval_lattice(forwardlattice)
    backlattice = pio.eval_lattice(backlattice)
    lattice = pio.eval_lattice(lattice)
    
    numerical = -1*(forwardlattice.TotalPE()-backlattice.TotalPE())/(2.0*h)
    analytical = numpy.dot(lattice.Force, normalise(HVec))
    error = numerical - analytical
    ratio = analytical/numerical
    print "NUMERICAL FORCE ", numerical," ANALYTICAL FORCE ",analytical, " ERROR ",error," RATIO ",ratio
    
    

def magnitude(vector):
	mag = math.sqrt(numpy.dot(vector,vector))	
	return mag

	
	
def normalise(vector):
	
	mag = magnitude(vector)
	if mag == 0:
		print "Can't normalise a zero vector" 
		return vector
	vector = vector/mag
	return vector
	
def separation(lattice1, lattice2):
	vector = separation_vector(lattice1, lattice2)
	return math.sqrt(numpy.dot(vector, vector))

def separation_vector(lattice1, lattice2):
	try:
		sep_vector = lattice1.Pos - lattice2.Pos
	except:
		sep_vector = lattice1.Pos - lattice2
		
	if(lattice1.PBC[0]==1):
		sep_vector[0::3] -= numpy.rint((sep_vector[0::3])/lattice1.Dim[0])*lattice1.Dim[0]
	if(lattice1.PBC[1]==1):
		sep_vector[1::3] -= numpy.rint((sep_vector[1::3])/lattice1.Dim[1])*lattice1.Dim[1]
	if(lattice1.PBC[2]==1):
		sep_vector[2::3] -= numpy.rint((sep_vector[2::3])/lattice1.Dim[2])*lattice1.Dim[2]
	return sep_vector



def separation_vector_raw(lattice1, lattice2, PBC, Dim):
	sep_vector = lattice1 - lattice2
		
	if(PBC[0]==1):
		sep_vector[0::3] -= numpy.rint((sep_vector[0::3])/Dim[0])*Dim[0]
	if(PBC[1]==1):
		sep_vector[1::3] -= numpy.rint((sep_vector[1::3])/Dim[1])*Dim[1]
	if(PBC[2]==1):
		sep_vector[2::3] -= numpy.rint((sep_vector[2::3])/Dim[2])*Dim[2]
	return sep_vector

def separation_vector_point_point(lattice1, point1, point2):
	sep_vector = point1-point2
	if(lattice1.PBC[0]==1):
		sep_vector[0::3] -= numpy.rint((sep_vector[0::3])/lattice1.Dim[0])*lattice1.Dim[0]
	if(lattice1.PBC[1]==1):
		sep_vector[1::3] -= numpy.rint((sep_vector[1::3])/lattice1.Dim[1])*lattice1.Dim[1]
	if(lattice1.PBC[2]==1):
		sep_vector[2::3] -= numpy.rint((sep_vector[2::3])/lattice1.Dim[2])*lattice1.Dim[2]
	return sep_vector	


def get_max_move_atom(lattice1,lattice2):
    #returns the PTA 
    dimx = lattice1.Dim[0]
    dimy = lattice1.Dim[1]
    dimz = lattice1.Dim[2]
    vector3 = numpy.zeros(len(lattice1.Pos),float)
    rmax = 0
    for i in numpy.arange(len(lattice1.Pos)/3):
        x_sep = lattice1.Pos[i*3] - lattice2.Pos[i*3]
        if(lattice1.PBC[0]==1):
            x_sep -= int(x_sep/dimx)*dimx
            x_sep -= int(2*x_sep/dimx)*dimx
        y_sep = lattice1.Pos[i*3+1] - lattice2.Pos[i*3+1]
        if(lattice1.PBC[1]==1):
            y_sep -= int(y_sep/dimy)*dimy
            y_sep -= int(2*y_sep/dimy)*dimy
        z_sep = lattice1.Pos[i*3+2] - lattice2.Pos[i*3+2]
        if(lattice1.PBC[2]==1):
            z_sep -= int(z_sep/dimz)*dimz
            z_sep -= int(2*z_sep/dimz)*dimz    
        
        r = math.sqrt(x_sep**2 + y_sep**2 + z_sep**2)
        if(r>rmax):
            maxmove = i
            rmax = r
        vector3[i*3]=x_sep
        vector3[i*3+1]=y_sep
        vector3[i*3+2]=z_sep
    return maxmove    


def atom_separation(atom,lattice1,lattice2):
    #returns atom separation between images
    #vector3 = []
    dimx = lattice1.Dim[0]
    dimy = lattice1.Dim[1]
    dimz = lattice1.Dim[2]
    x_sep = lattice1.Pos[atom*3] - lattice2.Pos[atom*3]
    if(lattice1.PBC[0]==1):
        x_sep -= numpy.rint(x_sep/dimx)*dimx
        
    y_sep = lattice1.Pos[atom*3+1] - lattice2.Pos[atom*3+1]
    if(lattice1.PBC[1]==1):
        y_sep -= numpy.rint(x_sep/dimx)*dimx
        
    z_sep = lattice1.Pos[atom*3+2] - lattice2.Pos[atom*3+2]
    if(lattice1.PBC[2]==1):
        z_sep -= numpy.rint(x_sep/dimx)*dimx  
        
    r = math.sqrt(x_sep**2 + y_sep**2 + z_sep**2)    
    return r
    
def atom_separation_pos(atom,lattice1,Pos):
    dimx = lattice1.Dim[0]
    dimy = lattice1.Dim[1]
    dimz = lattice1.Dim[2]
    x_sep = lattice1.Pos[atom*3] - Pos[atom*3]
    if(lattice1.PBC[0]==1):
        x_sep -= numpy.rint(x_sep/dimx)*dimx
        
    y_sep = lattice1.Pos[atom*3+1] - Pos[atom*3+1]
    if(lattice1.PBC[1]==1):
        y_sep -= numpy.rint(y_sep/dimy)*dimy
        
    z_sep = lattice1.Pos[atom*3+2] - Pos[atom*3+2]
    if(lattice1.PBC[2]==1):
        z_sep -= numpy.rint(z_sep/dimz)*dimz
        
    r = math.sqrt(x_sep**2 + y_sep**2 + z_sep**2)    
    return r
    
def atom_pair_separation_local(atom1, atom2, lattice):
    dimx = lattice.Dim[0]
    dimy = lattice.Dim[1]
    dimz = lattice.Dim[2]
    
    x_sep = lattice.Pos[atom1*3] - lattice.Pos[atom2*3]
    if(lattice.PBC[0]==1):
        x_sep -= numpy.rint(x_sep/dimx)*dimx
        
    y_sep = lattice.Pos[atom1*3+1] - lattice.Pos[atom2*3+1]
    if(lattice.PBC[1]==1):
        y_sep -= numpy.rint(y_sep/dimy)*dimy
        
    z_sep = lattice.Pos[atom1*3+2] - lattice.Pos[atom2*3+2]
    if(lattice.PBC[2]==1):
        z_sep -= numpy.rint(z_sep/dimz)*dimz 
        
    r = math.sqrt(x_sep**2 + y_sep**2 + z_sep**2)    
    return r	

def atom_pair_separation(atom,lattice1,atom2,lattice2):
    #returns r, xsep,ysep,zsep between 2 atoms in 2 vectors
    dimx = lattice1.Dim[0]
    dimy = lattice1.Dim[1]
    dimz = lattice1.Dim[2]
    x_sep = lattice1.Pos[atom*3] - lattice2.Pos[atom2*3]
    if(lattice1.PBC[0]==1):
        x_sep -= numpy.rint(x_sep/dimx)*dimx
    y_sep = lattice1.Pos[atom*3+1] - lattice2.Pos[atom2*3+1]
    if(lattice1.PBC[1]==1):
        y_sep -= numpy.rint(y_sep/dimy)*dimy
    z_sep = lattice1.Pos[atom*3+2] - lattice2.Pos[atom2*3+2]
    if(lattice1.PBC[2]==1):
        z_sep -= numpy.rint(z_sep/dimz)*dimz 
          
    r = math.sqrt(x_sep**2 + y_sep**2 + z_sep**2)    
    return r,x_sep,y_sep,z_sep    
 
 
def get_displaced_atom_array(lattice1,lattice2,threshold):
    #returns the PTA 
    vector3 = numpy.zeros(len(lattice1.Pos),float)
    atoms = []
    sep = []
    rmax = 0
    for i in numpy.arange(len(lattice1.Pos)/3):
        x_sep = lattice1.Pos[i*3] - lattice2.Pos[i*3]
        if(lattice1.PBC[0]):
            x_sep -= int(x_sep/dimx)*dimx
            x_sep -= int(2*x_sep/dimx)*dimx
        y_sep = lattice1.Pos[i*3+1] - lattice2.Pos[i*3+1]
        if(lattice1.PBC[1]):
            y_sep -= int(y_sep/dimy)*dimy
            y_sep -= int(2*y_sep/dimy)*dimy
        z_sep = lattice1.Pos[i*3+2] - lattice2.Pos[i*3+2]
        if(lattice1.PBC[2]):
            z_sep -= int(z_sep/dimz)*dimz
            z_sep -= int(2*z_sep/dimz)*dimz    
        
        r = math.sqrt(x_sep**2 + y_sep**2 + z_sep**2)
        if(r>threshold):
            atoms.append(i)
            sep.append(r)

    return atoms,sep 
    
def separation_mod(lattice1, lattice2):
    
    vector = separation_vector(lattice1, lattice2)
    sep = math.sqrt(numpy.dot(vector, vector))
    sep_vector = vector[0::3]**2 + vector[1::3]**2 + vector[2::3]**2
    primary_atom = sep_vector.argmax()
    max_atom_sep = sep_vector[primary_atom]
    
    return sep, primary_atom, math.sqrt(max_atom_sep)
    	
    
def saddle_check(lattice):
  """ Displaces along lattices min_mode vector in both directions and relaxes using SD to ensure the points role away from each other. """
  displacement_step = 0.1
  try: print "Rank", lattice.rank, "Saddle check"
  except: pass
  orig_pos = copy.copy(lattice.Pos)
  displaced_lattice = copylattice(lattice)	
  displaced_lattice.Pos = orig_pos - lattice.min_mode * displacement_step
  displaced_lattice = pminimise.minimise_lattice(displaced_lattice, "SD", 10)
  end1 = copy.copy(displaced_lattice.Pos)
  end1F = copy.copy(displaced_lattice.Force)
  end1disp = separation(displaced_lattice, lattice.Start_Pos)
  displaced_lattice.Pos = orig_pos + lattice.min_mode * displacement_step
  displaced_lattice = pminimise.minimise_lattice(displaced_lattice, "SD", 10)
  end2 = copy.copy(displaced_lattice.Pos)
  end2disp = separation(displaced_lattice, lattice.Start_Pos)
  end2F = copy.copy(displaced_lattice.Force)
  
  if(end1disp < end2disp):
    lattice.Approx_Start_Pos = end1
    lattice.Approx_End_Pos = end2
  else:
    lattice.Approx_Start_Pos = end2
    lattice.Approx_End_Pos = end1

  saddle_product = numpy.vdot(end1F, lattice.min_mode) * numpy.vdot(end2F, lattice.min_mode)
  displacement = magnitude(end2-end1)
  print "Relative Displacement:", displacement
  print "Saddle Product:", saddle_product

  if(saddle_product > 0 or displacement < displacement_step):
    lattice.Failed = True
  else:
    lattice.Failed = False
  return lattice
    
def local_saddle_check(lattice):
  try: print "Rank", lattice.rank, "LBFGS local minima check"
  except: pass
  min_lattice = copylattice(lattice)
  min_lattice.Pos = copy.copy(lattice.Start_Pos)
  print "Approx Start", separation_mod(min_lattice, lattice.Approx_Start_Pos)
  print "Approx End", separation_mod(min_lattice, lattice.Approx_End_Pos)
  (min_lattice.Pos, energy, data) = fmin_l_bfgs_b(scipy_eval_lattice, lattice.Approx_Start_Pos, fprime=None, args=(min_lattice,1), m=10, factr=1E-8, maxfun=100, iprint=0)
  sep_mod_towards = separation_mod(min_lattice, lattice.Start_Pos)
  print "Displacement from Origin:", sep_mod_towards
  lattice.Approx_Start_Pos = min_lattice.Pos
  if sep_mod_towards[2] > 1:
    print "Converged on non-local saddle"
    lattice.Failed = True
  else:
    print "Converged on local saddle"
    lattice.Failed = False
    
  try: print "Rank", lattice.rank, "LBFGS local minima check complete"
  except: pass
  return lattice
    
	
def roll_check(lattice):
	
	sep_mod = separation_mod(lattice, lattice.Start_Pos)
	orig_pos = copy.copy(lattice.Pos)
	
	dispdotmin = numpy.dot(separation_vector(lattice, lattice.Start_Pos), lattice.min_mode)
	disp_vector = normalise(lattice.min_mode) * 0.1 * numpy.sign(dispdotmin)
		
	displaced_lattice = copylattice(lattice)	
	displaced_lattice.Pos -= disp_vector
	
	
	displaced_lattice = minimise.minimise_lattice(displaced_lattice, "SD", 0.01)
	
	sep_mod_towards = separation_mod(displaced_lattice, lattice.Start_Pos)
	print "Roll Check Trial Towards:", sep_mod, sep_mod_towards
	if(sep_mod_towards[0] > sep_mod[0]):
		print "Roll Check Trial Towards Failed"
		lattice.Failed = True
		lattice.Approx_Start_Pos = copy.copy(displaced_lattice.Pos)
		return lattice

	
	(displaced_lattice.Pos, energy, data) = fmin_l_bfgs_b(scipy_eval_lattice, displaced_lattice.Pos, fprime=None, args=(displaced_lattice,1), m=10, factr=1E-8, maxfun=1000, iprint=0)
	
	sep_mod_towards = separation_mod(displaced_lattice, lattice.Start_Pos)
	print "Roll Check Complete Towards:", sep_mod, sep_mod_towards
	if(sep_mod_towards[2] > 0.5):
		print "Roll Check Complete Towards Failed"
		lattice.Failed = True
		lattice.Approx_Start_Pos = copy.copy(displaced_lattice.Pos)
		return lattice

	lattice.Approx_Start_Pos = copy.copy(displaced_lattice.Pos)
	
	displaced_lattice.Pos = copy.copy(lattice.Pos) + disp_vector
	

	displaced_lattice = minimise.minimise_lattice(displaced_lattice, "SD", 0.01)
	
	sep_mod_away = separation_mod(displaced_lattice, lattice.Start_Pos)
	print "Roll Check Trial Away:", sep_mod, sep_mod_away
	if(sep_mod_away[0] < sep_mod[0]):
		print "Roll Check Trial Away Failed"
		lattice.Failed = True
		lattice.Approx_End_Pos = copy.copy(displaced_lattice.Pos)
		return lattice
	
	(displaced_lattice.Pos, energy, data) = fmin_l_bfgs_b(scipy_eval_lattice, displaced_lattice.Pos, fprime=None, args=(displaced_lattice,1), m=10, factr=1E-12, maxfun=1000, iprint=0)
	
	sep_mod_away = separation_mod(displaced_lattice, lattice.Start_Pos)
	print "Roll Check Complete away:", sep_mod, sep_mod_away
	if(sep_mod_away[2] < 0.5):
		print "Roll Check Complete Away Failed"
		lattice.Failed = True
		lattice.Approx_End_Pos = copy.copy(displaced_lattice.Pos)		
		return lattice	
	
	lattice.End_Pos = displaced_lattice.Pos
	lattice.Approx_End_Pos = lattice.End_Pos

	print "We converged on a local saddle"
	return lattice	

    
def correct_drift_non_corresponding(input_lattice, ref_lattice):
	
	# as a result of recieving 2 lattices for the input and ref, with the input having undregone uncorrelated drift (2 atom sets dont match) i need a correction mechanism
	# this constructs x,y,z rdfs for the input and ref and then seeks to minimise the error by lining the plots up		

	num_bins = 40
	x_width = input_lattice.Dim[0]/num_bins
	y_width = input_lattice.Dim[1]/num_bins
	z_width = input_lattice.Dim[2]/num_bins
	input_rdf = numpy.zeros((3, num_bins), dtype=float)
	ref_rdf = numpy.zeros((3, num_bins), dtype=float)
	
	#loop over all the input atoms
	for i in range (0, input_lattice.NAtoms):
		x = input_lattice.Pos[3*i]
		y = input_lattice.Pos[3*i+1]
		z = input_lattice.Pos[3*i+2]
		x_bin = int(x / x_width)
		y_bin = int(x / x_width)
		z_bin = int(x / x_width)
		input_rdf[0][x_bin] += 1.0 / input_lattice.NAtoms
		input_rdf[1][y_bin] += 1.0 / input_lattice.NAtoms
		input_rdf[2][z_bin] += 1.0 / input_lattice.NAtoms
		
	for i in range (0, ref_lattice.NAtoms):
		x = ref_lattice.Pos[3*i]
		y = ref_lattice.Pos[3*i+1]
		z = ref_lattice.Pos[3*i+2]
		x_bin = int(x / x_width)
		y_bin = int(x / x_width)
		z_bin = int(x / x_width)
		ref_rdf[0][x_bin] += 1.0 / ref_lattice.NAtoms
		ref_rdf[1][y_bin] += 1.0 / ref_lattice.NAtoms
		ref_rdf[2][z_bin] += 1.0 / ref_lattice.NAtoms
		
	
	f = open("input_rdf.csv", "w")
	writer = csv.writer(f)
	writer.writerows(input_rdf)
	f.close()
	
	f = open("ref_rdf.csv", "w")
	writer = csv.writer(f)
	writer.writerows(ref_rdf)
	f.close()



def image_separation(lattice1,lattice2):
    dimx = lattice1.Dim[0]
    dimy = lattice1.Dim[1]
    dimz = lattice1.Dim[2]
    if(import_c==1):
        #print "RUNNING CLIB"
        vector3 = numpy.zeros(len(lattice1.Pos),float)
        clibs.c_util.image_separation(lattice1.Pos,lattice2.Pos,dimx,dimy,dimz,vector3,lattice1.PBC[0],lattice1.PBC[1],lattice1.PBC[2])    
    else:    
        #returns vector separation between images
	vector3 = separation_vector(lattice1, lattice2)
    return vector3    


def wrap_periodic(lattice):

	if(lattice.PBC[0]==1):
		lattice.Pos[0::3] -= numpy.floor(lattice.Pos[0::3] / lattice.Dim[0]) * lattice.Dim[0]
	if(lattice.PBC[1]==1):
		lattice.Pos[1::3] -= numpy.floor(lattice.Pos[1::3] / lattice.Dim[1]) * lattice.Dim[1]
	if(lattice.PBC[2]==1):
		lattice.Pos[2::3] -= numpy.floor(lattice.Pos[2::3] / lattice.Dim[2]) * lattice.Dim[2]
		

	return lattice
		
	
	
def proj(vector1, vector2):
	mag = magnitude(vector2)
	if mag == 0:
		print "Can't project onto a zero vector" 
		return vector2
	vector2 = vector2/mag
	return vector2 * (numpy.dot(vector1, vector2))
	#return vector2 * (numpy.dot(vector1, vector2)/mag)

		
def nearest_neighbours_sorted(lattice, index):
	"""Calculates the seperation of atoms in lattice from atom index, or position tuple sorted by proximity"""
#	print index
	try:
		len(index)
		ref_pos = index
	except:
			ref_pos = (lattice.Pos[index*3],  lattice.Pos[index*3+1], lattice.Pos[index*3+2])


	sep_list = []
	sep_vector = numpy.zeros(lattice.NAtoms*3)
	sep_vector[0::3] = lattice.Pos[0::3] - ref_pos[0]
	sep_vector[1::3] = lattice.Pos[1::3] - ref_pos[1]
	sep_vector[2::3] = lattice.Pos[2::3] - ref_pos[2]
	
	if(lattice.PBC[0]==1):
		sep_vector[0::3] -= numpy.rint((sep_vector[0::3])/lattice.Dim[0])*lattice.Dim[0]	
	if(lattice.PBC[1]==1):
		sep_vector[1::3] -= numpy.rint((sep_vector[1::3])/lattice.Dim[1])*lattice.Dim[1]	
	if(lattice.PBC[2]==1):
		sep_vector[2::3] -= numpy.rint((sep_vector[2::3])/lattice.Dim[2])*lattice.Dim[2]			

	for i in xrange(0,lattice.NAtoms):
		r = math.sqrt(numpy.vdot(sep_vector[3*i:3*i+3],sep_vector[3*i:3*i+3]))
		sep_list.append([i, r])
		
	return sorted(sep_list, key=lambda sep: sep[1])
	
	
def neighbours_displacement_vectors(lattice, index):
	"""Calculates the seperation of atoms in lattice from point"""
#	print index
	try:
		len(index)
		ref_pos = index
	except:
			ref_pos = (lattice.Pos[index*3],  lattice.Pos[index*3+1], lattice.Pos[index*3+2])


	sep_list = []
	sep_vector = numpy.zeros(lattice.NAtoms*3)
	sep_vector[0::3] = lattice.Pos[0::3] - ref_pos[0]
	sep_vector[1::3] = lattice.Pos[1::3] - ref_pos[1]
	sep_vector[2::3] = lattice.Pos[2::3] - ref_pos[2]
	
	if(lattice.PBC[0]==1):
		sep_vector[0::3] -= numpy.rint((sep_vector[0::3])/lattice.Dim[0])*lattice.Dim[0]	
	if(lattice.PBC[1]==1):
		sep_vector[1::3] -= numpy.rint((sep_vector[1::3])/lattice.Dim[1])*lattice.Dim[1]	
	if(lattice.PBC[2]==1):
		sep_vector[2::3] -= numpy.rint((sep_vector[2::3])/lattice.Dim[2])*lattice.Dim[2]			

	for i in xrange(0,lattice.NAtoms):
		r = math.sqrt(numpy.vdot(sep_vector[3*i:3*i+3],sep_vector[3*i:3*i+3]))
		sep_list.append([i, r, normalise(sep_vector[3*i:3*i+3])])
		
	return sorted(sep_list, key=lambda sep: sep[1])
	
	
   
def purge_lattice(lattice):
	"""Removes some unnessential information from lattice for faster transmission"""
	try: del(lattice.Force)
	except: pass
	try: del(lattice.PE)
	except: pass
	try: del(lattice.lmp)
	except: pass
	return lattice
	

def re_init_lattice(lattice):
	lattice.Force = numpy.zeros(3*lattice.NAtoms, dtype=float)
	lattice.PE = numpy.zeros(lattice.NAtoms, dtype=float)
	return lattice
	
    	
def atomic_displacement_sorted(input_lattice, end_lattice):
	sep_list = []
	for i in xrange(0, input_lattice.NAtoms):
		
		x_sep = input_lattice.Pos[i*3] - end_lattice.Pos[i*3]
		if(input_lattice.PBC[0]):
			x_sep -= int(x_sep/input_lattice.Dim[0])*input_lattice.Dim[0]
			x_sep -= int(2*x_sep/input_lattice.Dim[0])*input_lattice.Dim[0]

		y_sep = input_lattice.Pos[i*3+1] - end_lattice.Pos[i*3+1]
		if(input_lattice.PBC[1]):
			y_sep -= int(y_sep/input_lattice.Dim[1])*input_lattice.Dim[1]
			y_sep -= int(2*y_sep/input_lattice.Dim[1])*input_lattice.Dim[1]

		z_sep = input_lattice.Pos[i*3+2] - end_lattice.Pos[i*3+2]
		if(input_lattice.PBC[2]):
			z_sep -= int(z_sep/input_lattice.Dim[2])*input_lattice.Dim[2]
			z_sep -= int(2*z_sep/input_lattice.Dim[2])*input_lattice.Dim[2]
		r = math.sqrt(x_sep**2 + y_sep**2 + z_sep**2)
		sep_list.append((i, r))
		
	return sorted(sep_list, key=lambda sep: -sep[1])

def scipy_eval_lattice(Pos, lattice, null):
	''' This is an energy calculator designed to work with the scipy modules, the null was required due to a series versus instance error '''
	lattice.Pos = copy.copy(Pos)
	lattice = pio.eval_lattice(lattice)
	return lattice.TPE, -lattice.Force
	
      
	
	
def scipy_eval_lattice_constrained(Pos, lattice, constraint):
	lattice.Pos = copy.copy(Pos)
	lattice = pio.eval_lattice(lattice)
	dp = numpy.dot(lattice.Force, constraint)
	lattice.Force -= dp*constraint
	return lattice.TotalPE(), -lattice.Force

def scipy_calc_energy(Pos, lattice, null):
	lattice.Pos = copy.copy(Pos)
	lattice = pio.eval_lattice(lattice)
	return lattice.TotalPE()
	
def scipy_calc_forces(Pos, lattice, null):
	lattice.Pos = copy.copy(Pos)
	lattice = pio.eval_lattice(lattice)
	return -lattice.Force
	
def check_potential(lattice):
	check_grad(scipy_calc_energy, scipy_calc_forces, lattice.Pos,lattice, 1)	
	
def constrained_lbfgs(lattice, n, min_tol=1E-12):
	orig_energy = lattice.TPE
	mintol = min_tol / (numpy.finfo(float).eps)
	(lattice.Pos, lattice.TPE, data) = fmin_l_bfgs_b(scipy_eval_lattice_constrained, lattice.Pos, fprime=None, args=(lattice,n), m=10, factr=mintol, maxfun=20, iprint=-1)
#	lattice = pio.eval_lattice(lattice)
	print "Delta Energy:", orig_energy-lattice.TPE
	return lattice	
	
def black_hole_lbfgs(lattice, point, alpha, min_tol=1E-12, cutoff=100):
	orig_energy = lattice.TPE
	mintol = min_tol / (numpy.finfo(float).eps)
	(lattice.Pos, lattice.TPE, data) = fmin_l_bfgs_b(scipy_eval_lattice_black_hole, lattice.Pos, fprime=None, args=(lattice,point,alpha, cutoff), m=10, factr=mintol, maxfun=1000, iprint=0)
	lattice = pio.eval_lattice(lattice)
	return lattice	
	
	
	
def constrained_sd(lattice, constraint, step_orig=0.02):

	step = step_orig
	#lattice.eval()
	dp = numpy.dot(lattice.Force, constraint)
	print "Perp Force", dp
	if(math.fabs(dp)<1E-4):
		return lattice
	lattice.Force -= dp*constraint
	lattice.Pos += normalise(lattice.Force) * step
	Force_Old = copy.copy(lattice.Force)
	fdp_old = 0
	while(1):
		lattice.eval()
		dp = numpy.dot(lattice.Force, constraint)
		lattice.Force -= dp*constraint
		lattice.Pos += normalise(lattice.Force) * step
		fdp = numpy.dot(lattice.Force, Force_Old)
		if(fdp > fdp_old):
			step*=1.5
		else:
			step/=2
		fdp_old = fdp
			
		if(step < step_orig/4):
			step = step_orig/4
		elif(step > step_orig*8):
			step = step_orig*8
			
		print fdp, step
		if(fdp<0):
			return lattice
		else:
			Force_Old = copy.copy(lattice.Force)

	
def copylattice(lattice):
	# been using deepcopy because copy.copy uses references for internal structures (like lattice.Pos)
	# well this doesn't work with the lmp interface object, will create a soft copy and then deep copy the important stuff
	try:
	#	del(lattice.lammps_types)
		obj = copy.copy(lattice)
		#obj = pio.Lattice()
		for key in lattice.__dict__:
			#print key
			if key == "lmp":
				#obj.__dict__[key] = copy.copy(lattice.__dict__[key])
				#print "Deleting lmp"
				del(obj.lmp)
				#pass
				
			else:
				try:
					# try a deep copy on attribute
					#print "Copying key:", key
					obj.__dict__[key] = copy.deepcopy(lattice.__dict__[key])
				except:
					# fail and make a soft copy
					print "Making a soft copy of:", key
					
					try:
						obj.__dict__[key] = copy.copy(lattice.__dict__[key])
					except:
						print "Not copying:", key
		
		try:del(obj.lmp)
		except: pass
		return obj
	except:
		#del obj.key
		try:
			print "Failed to copy lattice key:", key
		except: 
			print "Something disconcerting is going on with this object. Who knows what is being copied."
		return obj
		

def copyobject(object):
	try:
		class Cobject:
			pass
		obj = Cobject()
		for key in object.__dict__:
			try:
				# try a deep copy on attribute
				print "Copying key:", key
				obj.__dict__[key] = copy.deepcopy(object.__dict__[key])
			except:
				# failed
				print "Skipping copy of:", key
		return obj
	except:
		print "Failed to copy object."
		
		
		
def benchforcecall(lattice):
	t = time.time() 
	for i in xrange(0,1000):
		lattice = lattice.eval()
	elapsed = (time.time() - t)/1000.0
	print elapsed, "seconds per force call"
	return		
	

def vineyard_rate(barrier, temp=300, prefactor=1e13):
	prob = prefactor * math.exp(-barrier / ( 8.617342295e-5 * temp))
	return prob



