"""
pio handles all of the input and output performed by pesto. This includes interfacing with MD and Visualisation packages.
"""

__author__ =  'Louis Vernon <louis.vernon@gmail.com>'




import os,sys
import numpy
import random
import commands
import math
import time
import putil
import pglobals
import copy
import gzip

verbose = True

try:
	import ctypes
except:	pass

class Lattice:
	""" The base lattice object contains standardised lattice attributes in addition to custom attributes required for specific MD implementations """
	NDefects = False
	workingdir = False
	def __init__(self):
		pass
	def MagForce(self):
		""" Returns the magnitude of the lattice force vector. Calculates when called. """
		return putil.magnitude(self.Force)
	def MaxForce(self):
		""" Return the maximum force on any single atom. """
		return math.sqrt((self.Force[0::3]**2 + self.Force[1::3]**2 + self.Force[2::3]**2).max())
	def TotalPE(self):
		""" Returns the total potential energy of the lattice. Lookup, no cost. """
		return self.TPE
	def TotalKE(self):
		""" Returns the kinetic energy of the system. Not really used. Calculates when called. """
		return numpy.sum(self.KE)
	def EnableVerbose(self, verbosity=True):
		""" Enable verbose output when using LAMMPS"""
		try:
			self.lmp.free()
			self.lmp.close()
			del(self.lmp)
		except: pass
		if(verbosity):
			self.lmp_args = []
		else:
			self.lmp_args = ["-echo","none", "-screen", "none", "-log", "none"]
			
	
	
	def RemoveAtom(self,index):  
		""" Allows deletion of atoms by index. You can provide a single index or an array or list of indices. """
		
		index = numpy.array(index)
		self.deletion_key = numpy.delete(numpy.arange(self.NAtoms), index)
		self.NAtoms-=index.size
		Nindex = numpy.concatenate((3*index, 3*index+1, 3*index+2), None)
		Nindex = numpy.sort(Nindex)
		self.Pos = numpy.delete(self.Pos, Nindex, None)
		self.Force = numpy.delete(self.Force, Nindex, None)
		self.PE = numpy.delete(self.PE, index, None)
		
		
		try:
			self.KE =numpy.delete(self.KE,index)
		except:
			pass
		self.Charge =numpy.delete(self.Charge,index,None)
		
		# treats species as numpy array by default, assumes list otherwise
		

		try:
			self.Specie = numpy.delete(self.Specie, index, None)
		except:
			self.Specie = numpy.delete(numpy.array(self.Specie), index, None).tolist()
	
		try:
			self.DefectMask = numpy.delete(self.DefectMask, Nindex, None)
		except: pass
		
	
		
		try:
			self.lmp.close()
			del(self.lmp)
			
		except:
			pass
			
	def AddAtom(self, Pos, Specie, Charge=0.0):
		""" Allows insertion of a single atom. Provide a 3 coordinate position, the specie and the charge (optional). """
		self.NAtoms+=1
		self.Pos = numpy.append(self.Pos, Pos)
		self.Specie = numpy.append(self.Specie, Specie)
		self.Force = numpy.zeros(3*self.NAtoms)
		self.DefectMask = numpy.append(self.DefectMask, [True,True,True])
		self.PE = numpy.zeros(self.NAtoms)
		try: self.Charge = numpy.append(self.Charge, Charge)
		except: pass
		try:
			self.lmp.free()
			self.lmp.close()
			del(self.lmp)
		except:
			pass		
				
	def eval(self):	
		""" Performs in place calculation of the forces and energies on the current lattice. Returns reference to lattice. """	
		return eval_lattice(self)





class Species_Info:
	""" Stores custom species information beyond label """
	Charge=0.0
	Mass = 0.0
	AtomsPerSpec = 0
	Label = ""


def _create_ram_disk(sizeMB):
	""" Creates a RAM disc on OSX. Should not be called by end user. """
	blocks = (sizeMB*2048)
	print "hdiutil eject /Volumes/ramdisk"
	print "diskutil erasevolume HFS+ \"ramdisk\" `hdiutil attach -nomount ram://"+str(blocks)+"` "
	os.system("hdiutil eject /Volumes/ramdisk")
	os.system("diskutil erasevolume HFS+ \"ramdisk\" `hdiutil attach -nomount ram://"+str(blocks)+"` ")


mdtype = os.environ.get("MDTYPE")
vistype = os.environ.get("VISTYPE")

forcecalls=0

print "mdtype",mdtype
print "vistype",vistype

valid_mdtypes=["LAMMPY","LAMMPS","LBOMD","DLPOLY","ARTMD"]
if(mdtype not in valid_mdtypes):
	print "No valid MDTYPE defined."
	print "Valid MTYPE's:", valid_mdtypes

lmp_args = ["-echo","none", "-screen", "none", "-log", "none"]
#lmp_args = ["-log", "none"]

# autodetect os
def _get_temp_path():
	if(os.uname()[0] == "Linux"):
		#rootfolder = "/dev/shm/"
		rootfolder = "/tmp/"
	elif(os.uname()[0] == "Darwin"):
		if(not os.path.exists("/Volumes/ramdisk/")):
			_create_ram_disk(256)
		rootfolder = "/Volumes/ramdisk/"
	else:
		rootfolder = "/tmp/"
	
	rootfolder += "PyKMC/"
		
	
			
	rootfolder += str(os.getpid()) + "/"
	
		
	if(not os.path.exists(rootfolder)):
		os.makedirs(rootfolder)
	return rootfolder	


pid = str(os.getpid())
	


def prepare_lattice_for_MD(lattice):
	""" Allocate arrays ready for a local MD simulation. """
	lattice.V = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.KE = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.Time = 0.0
	lattice.Temp = 0.0
	return lattice

	
def read_lattice(filename):
	""" Read in a lattice file depending on specified mdtype. """
	if(mdtype == "LBOMD" or mdtype == "LBOMDPIPING" ):
		lattice = read_lattice_lbomd(filename)
	elif(mdtype == "LAMMPS"):
		lattice = read_lattice_lammps(filename)
	elif(mdtype == "ARTMD"):
		lattice = read_lattice_artmd_buck(filename)
	elif(mdtype == "DLPOLY"):
		lattice = read_lattice_dlpoly(filename)
	elif(mdtype == "LAMMPY"):
		lattice = read_lattice_lammpy(filename)
	else:
		print "READ_LATTICE: Environmental variable MDTYPE not recognised or set."
		sys.exit()
		
	lattice.workingdir = os.getcwd()
	lattice = putil.wrap_periodic(lattice)
	# shift the lattice so the cell is centered at 0, 0, 0
	try:
		lattice.Dim[0] -= lattice.Dim[3]
		lattice.Dim[1] -= lattice.Dim[4]
		lattice.Dim[2] -= lattice.Dim[5]
		lattice.Pos[0::3] -= lattice.Dim[3]
		lattice.Pos[1::3] -= lattice.Dim[4]
		lattice.Pos[2::3] -= lattice.Dim[5]
	except: pass
	return lattice



	
def read_lattice_lbomd(filename):
	lattice = Lattice()
	lattice.Dim = [0]*6

	f = open(filename, "r")
	lattice.NAtoms = int(f.readline())

	lattice.Specie = []*lattice.NAtoms
	lattice.Pos = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.Charge = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.Force = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.PE = numpy.zeros((lattice.NAtoms), dtype=float)
#	lattice.KE = numpy.zeros((lattice.NAtoms), dtype=float)
	
	(lattice.Dim[0], lattice.Dim[1], lattice.Dim[2]) = f.readline().split()
	lattice.Dim[0] = float(lattice.Dim[0])
	lattice.Dim[1] = float(lattice.Dim[1])
	lattice.Dim[2] = float(lattice.Dim[2])
	
	for i in range(0,lattice.NAtoms):
		temparray = f.readline().split()
		lattice.Specie.append(temparray[0])
		lattice.Pos[3*i] = float(temparray[1])
		lattice.Pos[3*i+1] = float(temparray[2])
		lattice.Pos[3*i+2] = float(temparray[3])	
		lattice.Charge[i] = float(temparray[4])
		
	f.close()
	
	try:
		f = open("lbomd.IN", "r")
		f.readline()
		f.readline()
		f.readline()
		lattice.JobName = f.readline().split()
		temparray = f.readline().split()
		lattice.PBC = [0]*3
		lattice.PBC[0] = int(temparray[0])
		lattice.PBC[1] = int(temparray[1])
		lattice.PBC[2] = int(temparray[2])
		f.close()
		lattice = putil.wrap_periodic(lattice)
	except:
		print "lbomd.IN not found"
	
	
	return lattice





def read_lattice_dlpoly(filename):
	#ok dlpoly input file is the CONFIG file for pos
	#needs to read the flag to determine if reading pos/vel/forces
	lattice = Lattice()
	lattice.Dim = [0]*6
	lattice.PBC = [0]*3

	

	  
	try:
	  f = open("FIELD", "r")
	except:
	  print "Error: read_lattice_dlpoly: FIELD not found, quitting"
	  sys.exit() 
	
	lattice.NAtoms=0
	lattice.NSpecies = 0
	
#	lattice.NSpecies = int(f.readline())
#	lattice.species_info = {}
#	for i in range(0, lattice.NSpecies):
#		temparray = f.readline().split()
#		lattice.species_info[temparray[0]] = Species_Info()
#		lattice.species_info[temparray[0]].Charge=float(temparray[1])
#		lattice.species_info[temparray[0]].Mass = float(temparray[2])
#		lattice.species_info[temparray[0]].Label = temparray[-1].rstrip()
#	
	
	lattice.species_info = {}
	
	for line in f:
		array = line.split()
		#print line, array
		if(len(array)>0):
			if(array[0] == "NUMMOLS"):
			   lattice.NAtoms+=int(array[1])
			   lastlabel=str(prevarray[0])
			   nmols = int(array[1])
			if(len(str(array[0]))==2 or len(str(array[0]))==1 ):
				try:
					if(lattice.species_info[str(array[0])].Charge):
						pass
				except:
					lattice.species_info[str(array[0])] = Species_Info()
					lattice.species_info[str(array[0])].Mass = float(array[1])
					lattice.species_info[str(array[0])].Charge = float(array[2])
					lattice.species_info[str(array[0])].Label=lastlabel
					lattice.species_info[str(array[0])].AtomsPerSpec=nmols
					lattice.NSpecies+=1
					#lattice.MassPerSpec[str(array[0])] = float(array[1])
					#lattice.ChargePerSpec[str(array[0])] = float(array[2])
					#print "SPEC",  str(array[0]),lattice.MassPerSpec[str(array[0])],lattice.ChargePerSpec[str(array[0])]
		prevarray = array		
	  
	f.close()
	
	#print "IN DLPOLY READ"
	for key in lattice.species_info:
			print key, lattice.species_info[key].Mass,lattice.species_info[key].Charge,lattice_info[key].Label
	
	try:
		f = open(filename, "r")
	except:
		print filename + " not found, quitting"
		sys.exit()  
	  
	  
	lattice.JobName = f.readline().split("\n")[0] 
	pbc = f.readline().split()
	lattice.LevCfg = int(pbc[0])
	lattice.PBCtype = int(pbc[1])
	if(int(pbc[1])>0):
		lattice.PBC[0] = 1
		lattice.PBC[1] = 1
		lattice.PBC[2] = 1
	else:
		print "Error: read_lattice_dlpoly: pbc flag unrecognised "
		sys.exit() 
		
	lattice.avec = f.readline().split()
	lattice.bvec = f.readline().split()
	lattice.cvec = f.readline().split()
	xdims = [0]*3
	ydims = [0]*3
	zdims = [0]*3
	xdims[0] = float(lattice.avec[0])
	xdims[1] = float(lattice.bvec[0])
	xdims[2] = float(lattice.cvec[0])
	ydims[0] = float(lattice.avec[1])
	ydims[1] = float(lattice.bvec[1])
	ydims[2] = float(lattice.cvec[1])
	zdims[0] = float(lattice.avec[2])
	zdims[1] = float(lattice.bvec[2])
	zdims[2] = float(lattice.cvec[2])
	lattice.Dim[0] = numpy.max(xdims)
	lattice.Dim[1] = numpy.max(ydims)
	lattice.Dim[2] = numpy.max(zdims)
	lattice.Dim[3] = numpy.min(xdims)
	lattice.Dim[4] = numpy.min(ydims)
	lattice.Dim[5] = numpy.min(zdims)
	lattice.X0 = 0.5*(lattice.Dim[0]-lattice.Dim[3])
	lattice.Y0 = 0.5*(lattice.Dim[1]-lattice.Dim[4])
	lattice.Z0 = 0.5*(lattice.Dim[2]-lattice.Dim[5])
	
	lattice.Specie = []*lattice.NAtoms
	lattice.Pos = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.Charge = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.Force = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.PE = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.KE = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.V = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.DefectMask = numpy.ones((3*lattice.NAtoms), dtype=bool)	
	
	
	for i in range(0,lattice.NAtoms):
		symb = f.readline().split()
		lattice.Specie.append(symb[0])
		atomnum = int(symb[1])
		index = atomnum-1	
		if(lattice.LevCfg ==0):
			temparray = f.readline().split()
			lattice.Pos[3*index] = float(temparray[0])+lattice.X0
			lattice.Pos[3*index+1] = float(temparray[1])+lattice.Y0
			lattice.Pos[3*index+2] = float(temparray[2])+lattice.Z0	
			lattice.Charge[index] = float(lattice.species_info[lattice.Specie[index]].Charge)
			#print lattice.Specie[i],lattice.Pos[3*i],lattice.Pos[3*i+1],lattice.Pos[3*i+2],lattice.Charge[i]
		if(lattice.LevCfg ==1):
			temparray = f.readline().split()
			lattice.Pos[3*index] = float(temparray[0])+lattice.X0
			lattice.Pos[3*index+1] = float(temparray[1])+lattice.Y0
			lattice.Pos[3*index+2] = float(temparray[2])+lattice.Z0	
			lattice.Charge[index] = float(lattice.species_info[lattice.Specie[index]].Charge)
			temparray = f.readline().split()
			lattice.V[3*index] = float(temparray[0])
			lattice.V[3*index+1] = float(temparray[1])
			lattice.V[3*index+2] = float(temparray[2])
		
		if(lattice.LevCfg ==2):
			temparray = f.readline().split()
			lattice.Pos[3*index] = float(temparray[0])+lattice.X0
			lattice.Pos[3*index+1] = float(temparray[1])+lattice.Y0
			lattice.Pos[3*index+2] = float(temparray[2])+lattice.Z0	
			lattice.Charge[index] = float(lattice.species_info[lattice.Specie[index]].Charge)
			temparray = f.readline().split()
			lattice.V[3*index] = float(temparray[0])
			lattice.V[3*index+1] = float(temparray[1])
			lattice.V[3*index+2] = float(temparray[2])
			temparray = f.readline().split()
			lattice.Force[3*index] = float(temparray[0])
			lattice.Force[3*index+1] = float(temparray[1])
			lattice.Force[3*index+2] = float(temparray[2])
		
	lattice = putil.wrap_periodic(lattice)	
	return lattice




	
	
def read_lattice_artmd_buck(filename):	
	lattice = Lattice()
	lattice.Dim = [0]*6
	lattice.PBC = [0]*3
	try:
	  f = open(filename, "r")
	except:
	  print filename + " not found, quitting"
	  sys.exit()


	
	lattice.NAtoms = int(f.readline().split()[0])
#	print "test"
	f.readline();	

	lattice.Specie = []*lattice.NAtoms
	lattice.Pos = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.Charge = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.Force = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.PE = numpy.zeros((lattice.NAtoms), dtype=float)
#	lattice.KE = numpy.zeros((lattice.NAtoms), dtype=float)
#	lattice.V = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.DefectMask = numpy.ones((3*lattice.NAtoms), dtype=bool)
	  
	(lattice.Dim[0], lattice.Dim[1], lattice.Dim[2]) = f.readline().split()
	lattice.Dim[0] = float(lattice.Dim[0])
	lattice.Dim[1] = float(lattice.Dim[1])
	lattice.Dim[2] = float(lattice.Dim[2])
	
	
	for i in range(0, lattice.NAtoms):
		
		temparray = f.readline().split()
		
		lattice.Specie.append(temparray[3])
		lattice.Pos[3*i] = float(temparray[0])
		lattice.Pos[3*i+1] = float(temparray[1])
		lattice.Pos[3*i+2] = float(temparray[2])
	
	f.close()
	
	# artmd doesnt support per atom charges so we have to apply them via the buckingham parameters file (will this get used?)
	
#	filename = "buck.parms"
#	try:
#		f = open(filename, "r")
#		skip = float(f.readline().rstrip())
#		skip = int(skip)
#		for i in range(0,skip):
#			f.readline()
#		lattice.NSpecies = int(f.readline())
#		lattice.species_info = {}
#		for i in range(0, lattice.NSpecies):
#			temparray = f.readline().split()
#			lattice.species_info[temparray[0]] = Species_Info()
#			lattice.species_info[temparray[0]].Charge=float(temparray[1])
#			lattice.species_info[temparray[0]].Mass = float(temparray[2])
#			lattice.species_info[temparray[0]].Label = temparray[-1].rstrip()

	lattice.species_info = {}
	filename = "atomic.attribs"
	try:
		f=open(filename,"r")
		attribs = f.readlines()
#		print attribs
		# dangerous peice of code, an attacker could construct an atomic.attribs file that could damage your system
		# but if they are able to do that you are in trouble already!
		try:
			for line in attribs:
				exec(line)
		except:
			print "When attempting to apply the atomic attributes python had a read error, please check atomic.attribs"
			sys.exit()

		for i in range(0, lattice.NAtoms):

			lattice.Charge[i] = lattice.species_info[lattice.Specie[i]].Charge
	except:
	#  print filename + " not found, hope you don't want to perform any MD"
		pass
	  

		
	lattice.PBC[0] = 1
	lattice.PBC[1] = 1
	lattice.PBC[2] = 1
	lattice = putil.wrap_periodic(lattice)
	return lattice
	
	
def read_lattice_artmd(filename):	
	lattice = Lattice()
	lattice.Dim = [0]*6
	lattice.PBC = [0]*3
	try:
	  f = open(filename, "r")
	except:
	  print filename + " not found, quitting"
	  sys.exit()


	
	lattice.NAtoms = int(f.readline().split()[0])
#	print "test"
	f.readline();	

	lattice.Specie = []*lattice.NAtoms
	lattice.Pos = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.Charge = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.Force = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.PE = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.DefectMask = numpy.ones((3*lattice.NAtoms), dtype=bool)
#	lattice.KE = numpy.zeros((lattice.NAtoms), dtype=float)
#	lattice.V = numpy.zeros((3*lattice.NAtoms), dtype=float)


	  
	(lattice.Dim[0], lattice.Dim[1], lattice.Dim[2]) = f.readline().split()
	lattice.Dim[0] = float(lattice.Dim[0])
	lattice.Dim[1] = float(lattice.Dim[1])
	lattice.Dim[2] = float(lattice.Dim[2])
	
	
	for i in range(0, lattice.NAtoms):
		temparray = f.readline().split()
		lattice.Specie.append(temparray[3])
		lattice.Pos[3*i] = float(temparray[0])
		lattice.Pos[3*i+1] = float(temparray[1])
		lattice.Pos[3*i+2] = float(temparray[2])

	
	f.close()
	print "wag"
	
	# artmd doesnt support per atom charges so we have to apply them via the buckingham parameters file (will this get used?)
	
	
	try:
		try:
			f = open("buck.parms", "r")
		except:
			f = open("buckmorse.parms", "r")
	except:
	  print filename + " not found, quitting"
	  sys.exit()
	  
	skip = float(f.readline().rstrip())
	skip = int(skip)
	for i in range(0,skip):
		f.readline()
	lattice.NSpecies = int(f.readline())
	lattice.species_info = []
	for i in range(0, lattice.NSpecies):
		lattice.species_info.append(Species_Info())
		temparray = f.readline().split()
		lattice.species_info[i].Charge=float(temparray[1])
		lattice.species_info[i].Mass = float(temparray[2])
		lattice.species_info[i].Label = temparray[-1].rstrip()
		
	lattice.PBC[0] = 1
	lattice.PBC[1] = 1
	lattice.PBC[2] = 1
	lattice = putil.wrap_periodic(lattice)
	return lattice	
		
	#
	
def read_lattice_lammps(filename):
	lattice = Lattice()
	lattice.Dim = [0]*6
	lattice.PBC = [0]*3
	try:
	  f = open(filename, "r")
	except:
	  print filename + " not found, quitting"
	  sys.exit()
	line = f.readline()
	line_array = line.rstrip().split()
	lattice.lammps_jobfile = line_array[-1]
	try:
		g = open(line_array[-1])
	
		line = g.readline().rstrip()
		lattice.JobName = line
		lattice.lammps_input_files = line_array
		
		while True:
		  line = g.readline()
		  if not line:
		    break
		  line_array = line.rstrip().split()
		  if(line_array):
		    if(line_array[0] == "atom_style"):
		      lattice.lammps_atom_style = line_array[-1]
		    elif(line_array[0] == "boundary"):
		      if(line_array[-1] == "p"):
			lattice.PBC[2] = 1
		      if(line_array[-2] == "p"):
			lattice.PBC[1] = 1	    
		      if(line_array[-3] == "p"):
			lattice.PBC[0] = 1	    
		    elif(line_array[0] == "read_data"):
		      lattice.lammps_lattice_name = line_array[-1]
		g.close()
	except:
	  print "The first line of the lattice needs to end with the path to the lammps input file"
	  
	
	while True:
	  line = f.readline()
	  if not line:
	    break
	  line_array = line.rstrip().split()
	  if(line_array):
	    if(line_array[-1] == "atoms"):
	      lattice.NAtoms = int(line_array[0])
	      lattice.Specie = []
	      lattice.Pos = numpy.zeros((3*lattice.NAtoms), dtype=float)
	      lattice.Charge = numpy.zeros((lattice.NAtoms), dtype=float)
	      lattice.Force = numpy.zeros((3*lattice.NAtoms), dtype=float)
	      lattice.PE = numpy.zeros((lattice.NAtoms), dtype=float)
	      lattice.DefectMask = numpy.ones((3*lattice.NAtoms), dtype=bool)
#	      lattice.KE = numpy.zeros((lattice.NAtoms), dtype=float)
#	      lattice.V = numpy.zeros((3*lattice.NAtoms), dtype=float)
	    elif(line_array[-1] == "types"):
	      lattice.lammps_types = int(line_array[0])
	    elif(line_array[-1] == "xhi"):
	      lattice.Dim[3] = float(line_array[0])
	      lattice.Dim[0] = float(line_array[1])
	    elif(line_array[-1] == "yhi"):
	      lattice.Dim[4] = float(line_array[0])
	      lattice.Dim[1] = float(line_array[1])
	    elif(line_array[-1] == "zhi"):
	      lattice.Dim[5] = float(line_array[0])
	      lattice.Dim[2] = float(line_array[1])
	    elif(line_array[-1] == "Atoms"):
	      f.readline()
	      for i in range(0,lattice.NAtoms):
		      temparray = f.readline().split()
		      lattice.Specie.append(temparray[1])
		      lattice.Pos[3*i] = float(temparray[2])
		      lattice.Pos[3*i+1] = float(temparray[3])
		      lattice.Pos[3*i+2] = float(temparray[4])
		      try:
			lattice.Charge[i] = float(temparray[5])
		      except:
			lattice.Charge[i] = 0
			    
	   
	f.close() 
	lattice.Specie = numpy.array(lattice.Specie)
	
	return lattice    
	
	
def read_lattice_lammpy(filename):
#	global lmp
	# This is a little complicated. We read in the lattice file like usual, mainly to initialise the variables
	# an lmp object will be created that links to lammps using cython.
	
	

	lattice = Lattice()
	lattice.lmp_args = lmp_args
	lattice.Dim = [0]*6
	lattice.PBC = [0]*3
	try:
		try:
			f = gzip.open(filename, "r")
			line = f.readline()
		except:
			f = open(filename, "r")
			line = f.readline()
	except:
	  print filename + " not found, quitting"
	  sys.exit()
	#line = f.readline()
	#print line
	line_array = line.rstrip().split()
	lattice.lammps_jobfile = line_array[-1]
	try:
		g = open(line_array[-1])	  
		line = g.readline().rstrip()
		lattice.lammps_jobfile_array = []
		lattice.JobName = line
		lattice.lammps_jobfile_array.append(line)
		lattice.lammps_input_files = line_array
	except:
	  print "The first line of the lattice needs to end with the path to the lammps input file"
	  sys.exit()		
	#print lattice.JobName
	# At this point we will parse the job file, we can now actually execute it using the lammps object
#	lattice.lmp = lammps(copy.copy(lmp_args))

	
	#lattice.lammps_jobfile_array = []
	
	while True:
	  line = g.readline()
	  #print line
	  lattice.lammps_jobfile_array.append(line)
#	  lattice.lmp.command(line)
	  if not line:
	    break
	  line_array = line.rstrip().split()
	  if(line_array):
	    if(line_array[0] == "atom_style"):
	      lattice.lammps_atom_style = line_array[-1]
	    elif(line_array[0] == "boundary"):
	      if(line_array[-1] == "p"):
		lattice.PBC[2] = 1
	      if(line_array[-2] == "p"):
		lattice.PBC[1] = 1	    
	      if(line_array[-3] == "p"):
		lattice.PBC[0] = 1
	    # After reading in the periodic boundary conditions it is time to get lammps to process the lattice file:
#	      lattice.lmp.command("read_data " + filename)
#	      lattice.lammps_jobfile_array.append("read_data " + filename)
	      lattice.lammps_jobfile_array.append("read_data")
		
	
	lattice.lammps_lattice_name = filename
	
	
#	lattice.lammps_jobfile_array.append("compute total_pe all pe")
	lattice.lammps_jobfile_array.append("compute atom_pe all pe/atom")
	#lattice.lammps_jobfile_array.append("variable fx atom fx")
	#lattice.lammps_jobfile_array.append("variable fy atom fy")
	#lattice.lammps_jobfile_array.append("variable fz atom fz")
#	lattice.lammps_jobfile_array.append("variable e equal pe")


	g.close()
	
	
	
	while True:
	  line = f.readline() 
	  if not line:
	    break
	  line_array = line.rstrip().split()
	  if(line_array):
	    if(line_array[-1] == "atoms"):
	      lattice.NAtoms = int(line_array[0])
	      lattice.Specie = []
	      
	
	      for i in xrange(0, lattice.NAtoms):
	      	lattice.Specie.append("0")
	      lattice.Pos = numpy.zeros((3*lattice.NAtoms), dtype=float)
	      lattice.Force = numpy.zeros((3*lattice.NAtoms), dtype=float)
	      lattice.PE = numpy.zeros((lattice.NAtoms), dtype=float)
	      lattice.DefectMask = numpy.ones((3*lattice.NAtoms), dtype=bool)
		     
	      lattice.Charge = numpy.zeros((lattice.NAtoms), dtype=float)
	      lattice.index = numpy.zeros((lattice.NAtoms), dtype=int)
#	      lattice.KE = numpy.zeros((lattice.NAtoms), dtype=float)
#	      lattice.V = numpy.zeros((3*lattice.NAtoms), dtype=float)
	    elif(line_array[-1] == "types"):
	      lattice.lammps_types = int(line_array[0])
	    elif(line_array[-1] == "xhi"):
	      lattice.Dim[3] = float(line_array[0])
	      lattice.Dim[0] = float(line_array[1])
	    elif(line_array[-1] == "yhi"):
	      lattice.Dim[4] = float(line_array[0])
	      lattice.Dim[1] = float(line_array[1]) 
	    elif(line_array[-1] == "zhi"):
	      lattice.Dim[5] = float(line_array[0])
	      lattice.Dim[2] = float(line_array[1])
	    elif(line_array[-1] == "Atoms"):
	      f.readline()
	      for i in range(0,lattice.NAtoms):
			temparray = f.readline().split()
			id = int(temparray[0])-1
			#lattice.index[i] = id
			
			if(lattice.lammps_atom_style=="atomic" or "gc" in lattice.lammps_atom_style):
			      lattice.Specie[id] = temparray[1]
			      lattice.Pos[3*id] = float(temparray[2])
			      lattice.Pos[3*id+1] = float(temparray[3])
			      lattice.Pos[3*id+2] = float(temparray[4])
			      try:
					lattice.Charge[id] = float(temparray[5])
			      except:
					lattice.Charge[id] = 0
			elif(lattice.lammps_atom_style=="full"):
			      lattice.Specie[id] = temparray[1]
			      lattice.Pos[3*id] = float(temparray[4])
			      lattice.Pos[3*id+1] = float(temparray[5])
			      lattice.Pos[3*id+2] = float(temparray[6])
			      try:
					lattice.Charge[id] = float(temparray[3])
			      except:
				lattice.Charge[id] = 0	
                        elif("charge" in lattice.lammps_atom_style):
			  try:
			      lattice.Specie[id] = temparray[1]
                              lattice.Pos[3*id] = float(temparray[3])
                              lattice.Pos[3*id+1] = float(temparray[4])
                              lattice.Pos[3*id+2] = float(temparray[5])
			      lattice.Charge[id] = float(temparray[2])
			  except:
			      lattice.Specie[id] = temparray[1]
                              lattice.Pos[3*id] = float(temparray[2])
                              lattice.Pos[3*id+1] = float(temparray[3])
                              lattice.Pos[3*id+2] = float(temparray[4])
			      lattice.Charge[id] = 0				
				
                              

		
			else:
				print "Lammps input atom style not supported:", lattice.lammps_atom_style
				sys.exit()
				
	
		
	f.close()
	
	
	#lattice = wrap_periodic(lattice)
	
	lattice.species_info = {}
	filename = "atomic.attribs"
	try:
		f=open(filename,"r")
		attribs = f.readlines()
		#print attribs
		# dangerous peice of code, an attacker could construct an atomic.attribs file that could damage your system
		# but if they are able to do that you are in trouble already!
		try:
			for line in attribs:
				exec(line)
			#	print line
		except:
			print "When attempting to apply the atomic attributes python had a read error, please check atomic.attribs"
			sys.exit()

		try:
			for i in range(0, lattice.NAtoms):

				lattice.Charge[i] = lattice.species_info[lattice.Specie[i]].Charge
		except: pass
	except:
	#  print filename + " not found, hope you don't want to perform any MD" 
	  pass
	  
#	lattice.index = numpy.arange(lattice.NAtoms)
	lattice.Specie = numpy.array(lattice.Specie)
	
	return lattice   
	
def read_lattice_lammps_dump(filename):
	lattice = Lattice()
	lattice.Dim = [0]*6
	lattice.PBC = [0]*3
	lattice.Specie = []
	try:
	  f = open(filename, "r")
	except:
	  print filename + " not found, quitting"
	  sys.exit()

	while(True):
		line = f.readline()
		if(len(line) == 0): break
		if("NUMBER OF ATOMS" in line):
			lattice.NAtoms = int(f.readline())
		elif("BOX BOUNDS" in line):			
			line_array = line.rstrip().split()
			if(line_array[-1] == "pp"):
				lattice.PBC[-1] = 1
			else:
				lattice.PBC[-1] = 0
				
			if(line_array[-2] == "pp"):
				lattice.PBC[-2] = 1
			else:
				lattice.PBC[-2] = 0
			if(line_array[-3] == "pp"):
				lattice.PBC[-3] = 1
			else:
				lattice.PBC[-3] = 0	
			
			for i in xrange(0,3):
				line_array = f.readline().rstrip().split()
				lattice.Dim[i] = float(line_array[1])
				lattice.Dim[i+3] = float(line_array[0])
		elif("ITEM: ATOMS" in line):	      	
			lattice.Pos = numpy.zeros((3*lattice.NAtoms), dtype=float)
			lattice.Charge = numpy.zeros((lattice.NAtoms), dtype=float)
			lattice.Force = numpy.zeros((3*lattice.NAtoms), dtype=float)
			lattice.PE = numpy.zeros((lattice.NAtoms), dtype=float)
			lattice.DefectMask = numpy.ones((3*lattice.NAtoms), dtype=bool)
			#lattice.index = numpy.arange((lattice.NAtoms), dtype=int)
			
			for i in xrange(0, lattice.NAtoms):
				line_array = f.readline().rstrip().split()
				lattice.Specie.append(line_array[1])
				lattice.Pos[3*i:3*i+3] = line_array[2:5]
			
		
	f.close()
	return lattice
			
	
			
def read_lattice_eon(filename):	
	lattice = Lattice()
	lattice.Dim = [0]*6
	lattice.PBC = [0]*3
	try:
	  f = open(filename, "r")
	except:
	  print filename + " not found, quitting"
	  sys.exit()


	f.readline()
	f.readline()
	(lattice.Dim[0], lattice.Dim[1], lattice.Dim[2]) = f.readline().split()
	lattice.Dim[0] = float(lattice.Dim[0])
	lattice.Dim[1] = float(lattice.Dim[1])
	lattice.Dim[2] = float(lattice.Dim[2])
	
	f.readline()
	f.readline()
	f.readline()
	f.readline()
	
	lattice.NAtoms = int(f.readline())
	f.readline()
	f.readline()

	f.readline()

	lattice.Specie = []*lattice.NAtoms
	lattice.Pos = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.Charge = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.Force = numpy.zeros((3*lattice.NAtoms), dtype=float)
	lattice.PE = numpy.zeros((lattice.NAtoms), dtype=float)
	lattice.DefectMask = numpy.ones((3*lattice.NAtoms), dtype=bool)
	  

	for i in range(0, lattice.NAtoms):
		
		temparray = f.readline().split()
		
		lattice.Specie.append("1")
		lattice.Pos[3*i] = float(temparray[0])
		lattice.Pos[3*i+1] = float(temparray[1])
		lattice.Pos[3*i+2] = float(temparray[2])
	
	f.close()
	

		
	lattice.PBC[0] = 1
	lattice.PBC[1] = 1
	lattice.PBC[2] = 1

	return lattice		 				
		 								
	
	
def write_lattice(orig_lattice, filename):
#	print filename
#	print os.getcwd()
	lattice = putil.copylattice(orig_lattice)
	try:
		lattice.Dim[0] += lattice.Dim[3]
		lattice.Dim[1] += lattice.Dim[4]
		lattice.Dim[2] += lattice.Dim[5]
		lattice.Pos[0::3] += lattice.Dim[3]
		lattice.Pos[1::3] += lattice.Dim[4]
		lattice.Pos[2::3] += lattice.Dim[5]
	except: pass
	
	if(mdtype == "LBOMD" or mdtype == "LBOMDPIPING" ):
		write_lattice_lbomd(lattice, filename)
	elif(mdtype == "LAMMPS" or mdtype == "LAMMPY"):
		write_lattice_lammps(lattice, filename)
	elif(mdtype == "ARTMD_BUCK"):
		write_lattice_artmd_buck(lattice, filename)	
	elif(mdtype == "DLPOLY"):
		write_lattice_dlpoly(lattice, filename)		
	else:
		print "write_lattice: Environmental variable MDTYPE not recognised or set."
		sys.exit()
	return


def write_lattice_dlpoly(lattice, filename):
	f = open(filename, 'w')
	f.write(str(lattice.JobName) + "\n")
	#print "PBC",lattice.PBCtype
	f.write(str(lattice.LevCfg) +" "+str(lattice.PBCtype)+"\n")
	f.write(str(lattice.avec[0])+" "+str(lattice.avec[1])+" "+str(lattice.avec[2])+" "+"\n")
	f.write(str(lattice.bvec[0])+" "+str(lattice.bvec[1])+" "+str(lattice.bvec[2])+" "+"\n")
	f.write(str(lattice.cvec[0])+" "+str(lattice.cvec[1])+" "+str(lattice.cvec[2])+" "+"\n")
	for i in range(0, lattice.NAtoms):
		if(lattice.LevCfg ==0):
			index = i+1	
			f.write(str(lattice.Specie[i])+" "+str(index)+"\n")
			string = str(lattice.Pos[3*i]-lattice.X0) + " " + str(lattice.Pos[3*i+1]-lattice.Y0) + " "  + str(lattice.Pos[3*i+2]-lattice.Z0) +"\n"
			#print i,lattice.NAtoms,string
			f.write(string)
		if(lattice.LevCfg ==1):
			index = i+1	
			f.write(str(lattice.Specie[i])+" "+str(index)+"\n")
			string = str(lattice.Pos[3*i]-lattice.X0) + " " + str(lattice.Pos[3*i+1]-lattice.Y0) + " "  + str(lattice.Pos[3*i+2]-lattice.Z0) +"\n"
			f.write(string)	
			string = str(lattice.V[3*i]) + " " + str(lattice.V[3*i+1]) + " "  + str(lattice.V[3*i+2]) +"\n"
			f.write(string)
		if(lattice.LevCfg ==2):
			index = i+1	
			f.write(str(lattice.Specie[i])+" "+str(index)+"\n")
			string = str(lattice.Pos[3*i]-lattice.X0) + " " + str(lattice.Pos[3*i+1]-lattice.Y0) + " "  + str(lattice.Pos[3*i+2]-lattice.Z0) +"\n"
			f.write(string)	
			string = str(lattice.V[3*i]) + " " + str(lattice.V[3*i+1]) + " "  + str(lattice.V[3*i+2]) +"\n"
			f.write(string)
			string = str(lattice.Force[3*i]) + " " + str(lattice.Force[3*i+1]) + " "  + str(lattice.Force[3*i+2]) +"\n"
			f.write(string)		
	f.close()
	return

	f = open(str(dir)+"/animation-reference.xyz", 'w')
	f.write(str(lattice.NAtoms)+'\n')
	string = str(lattice.Dim[0]) + ' ' + str(lattice.Dim[1])  + ' ' + str(lattice.Dim[2]) + '\n'
	f.write(string)
	for i in range(0,lattice.NAtoms):
		num=i+1
		#print i,lattice.Specie[i]
		string = str(i)+" "+lattice.Specie[i] + " %14.7g" % lattice.Pos[3*i] + " %14.7g" % lattice.Pos[3*i+1] + " %14.7g" % lattice.Pos[3*i+2] + " %14.7g" % lattice.Charge[i] + " %14.7g" % lattice.KE[i] + " %14.7g" % lattice.PE[i] + " %14.7g" % lattice.Force[3*i] + " %14.7g" % lattice.Force[3*i+1] + " %14.7g" % lattice.Force[3*i+2] + "\n"
		f.write(string)
	f.close()	
	return	


def write_lattice_lbomd(lattice, filename):
	f = open(filename, 'w')
	f.write(str(lattice.NAtoms)+'\n')
	string = str(lattice.Dim[0]) + ' ' + str(lattice.Dim[1])  + ' ' + str(lattice.Dim[2]) + '\n'
	f.write(string)
	for i in range(0,lattice.NAtoms):
		string = lattice.Specie[i] + " %14.10g" % lattice.Pos[3*i] + " %14.10g" % lattice.Pos[3*i+1] + " %14.10g" % lattice.Pos[3*i+2] + " %14.10g" % lattice.Charge[i] + "\n"
		f.write(string)
	f.close()	
	return
	
def write_lattice_artmd_buck(lattice, filename):
	f = open(filename, 'w')
	string = "%5i" % lattice.NAtoms + "%3i" % 1 + "\nAutogenerated lattice using PyKMC\n"
	f.write(string)
	string = "%20.10f" %  lattice.Dim[0] + "%20.10f" % lattice.Dim[1] + "%20.10f" % lattice.Dim[2] + '\n'
	f.write(string)
#	f.write("(3f20.10,i5)\n")
	for i in range(0,lattice.NAtoms):
		string = "%20.10f" % lattice.Pos[3*i] + "%20.10f" % lattice.Pos[3*i+1] + "%20.10f" % lattice.Pos[3*i+2] + "%5i" % int(lattice.Specie[i]) + "\n"
		f.write(string)
	f.close()	
	return


def write_lattice_lammps(lattice, filename):
	if(lattice.NAtoms > 1000000):
		f = gzip.open(filename, 'w')
	else:
		f = open(filename, 'w')	
	
	
	
	
		NAtoms = lattice.NAtoms
#	f.write(qlattice.lammps_jobfile + "\n")
	try:
		f.write(" ".join(lattice.lammps_input_files) + "\n")
	except:
		pass
	#print " ".join(lattice.lammps_input_files)
	f.write(str(NAtoms) + " atoms\n")
	try:	f.write(str(lattice.lammps_types) + " atom types\n")
	except: f.write(str(len(numpy.unique(lattice.Specie))) + " atom types\n")
	f.write(str(lattice.Dim[3]) + " " + str(lattice.Dim[0]) + " xlo xhi\n")
	f.write(str(lattice.Dim[4]) + " " + str(lattice.Dim[1]) + " ylo yhi\n")
	f.write(str(lattice.Dim[5]) + " " + str(lattice.Dim[2]) + " zlo zhi\n")	
	f.write("Atoms\n\n")
	for i in range(0, NAtoms):
	  string = "%5i" % (i+1) + " " + str(lattice.Specie[i]) + " "
	  try:
		  if(lattice.Charge[0] != 0 or "charge" in lattice.lammps_atom_style ):
		    string +=  str(lattice.Charge[i]) + " "
	  except:
	  	pass
	  string +=  "%20.15f" % lattice.Pos[3*i] + " " + "%20.15f" % lattice.Pos[3*i+1] + " "  + "%20.15f" %  lattice.Pos[3*i+2]
	  string += "\n"
	  
	  f.write(string)
	f.flush()
	f.close()
	return
	  
def write_anim_ref_lbomd(lattice, dir):
	f = open(str(dir)+"/animation-reference.xyz", 'w')
	f.write(str(lattice.NAtoms)+'\n')
	string = str(lattice.Dim[0]) + ' ' + str(lattice.Dim[1])  + ' ' + str(lattice.Dim[2]) + '\n'
	f.write(string)
	for i in range(0,lattice.NAtoms):
		num=i+1
		#print i,lattice.Specie[i],lattice.Specie[i],lattice.Pos[3*i],lattice.Pos[3*i+1],lattice.Pos[3*i+2] ,lattice.Charge[i] ,lattice.KE[i] ,lattice.PE[i] ,lattice.Force[3*i] , lattice.Force[3*i+1] , lattice.Force[3*i+2] 
		string = str(i)+" "+lattice.Specie[i] + " %14.7g" % lattice.Pos[3*i] + " %14.7g" % lattice.Pos[3*i+1] + " %14.7g" % lattice.Pos[3*i+2] + " %14.7g" % lattice.Charge[i] + " %14.7g" % lattice.KE[i] + " %14.7g" % lattice.PE[i] + " %14.7g" % lattice.Force[3*i] + " %14.7g" % lattice.Force[3*i+1] + " %14.7g" % lattice.Force[3*i+2] + "\n"
		f.write(string)
	f.close()	
	return	

def write_xyz_lbomd(lattice, dir, prefix, number,time):
	number = "%04d" % (number)
	filename = str(dir)+"/"+str(prefix)+str(number)+".xyz"
	f = open(filename, 'w')
	f.write(str(lattice.NAtoms)+'\n')
	string = str(time) + '\n'
	f.write(string)
	for i in range(0,lattice.NAtoms):
		num=i+1	
		string = str(i)+" %14.7g" % lattice.Pos[3*i] + " %14.7g" % lattice.Pos[3*i+1] + " %14.7g" % lattice.Pos[3*i+2] + " %14.7g" % lattice.KE[i] +  "%14.7g" % lattice.PE[i] + "\n"
		f.write(string)
	f.close()	
	return
	


def eval_lattice(lattice):
#	global lmp

	
	
	if(mdtype == "LBOMD"):
		lattice = eval_lattice_lbomd(lattice)
	elif(mdtype == "LBOMDPIPING"):
		lattice = eval_lattice_lbomd_piping(lattice)
	elif(mdtype == "LAMMPS"):
		lattice = eval_lattice_lammps(lattice)
	elif(mdtype == "LAMMPY"):
		lattice = eval_lattice_lammpy(lattice)		
	elif(mdtype == "ARTMD_BUCK"):
		lattice = eval_lattice_artmd_buck(lattice)	
	elif(mdtype == "DLPOLY"):
		lattice = eval_lattice_dlpoly(lattice)			
		
	else:
		print "eval_lattice: Environmental variable MDTYPE not recognised or set."
		sys.exit()
	
	try:
		lattice.global_DefectMask
	except:
		lattice.global_DefectMask = numpy.ones(lattice.NAtoms*3, dtype=bool)
		lattice.global_DefectMask[numpy.where(lattice.Force==0)[0]]=False

	if(pglobals.status_queue.qsize()):
		print "Raising exception due to message in status queue."
		raise Exception("non-empty status queue")
		return 0
	

	try:
		lattice.Force *= lattice.DefectMask
		
	except:
		pass	
	return lattice
	
	
def eval_lattice_lbomd(lattice):
	global forcecalls
	#not sure of the best way to proceed here, but for now will work in a random directory :S bit unhappy because of lots of io, but it makes the code cleaner
	random.seed(time.time())
	tempfolder = "Temp" + str(random.randrange(0,100000))
	
	while(os.path.isdir(tempfolder)):
		tempfolder = "Temp" + str(random.randrange(0,100000))		
	os.makedirs(tempfolder)
	write_lattice_lbomd(lattice, tempfolder + "/lattice.dat")
	child = os.fork()
	if child:
		os.waitpid(child,0)
	else:	
		os.chdir(tempfolder)
		#make symlinks to all the input files
	#	os.symlink("LBOMD.exe",  tempfolder + "/LBOMD.exe")
		os.symlink("../Puphi.txt",  "Puphi.txt")
		#os.symlink("../08_08.spl",  "08_08.spl")	
		#os.symlink("../22_08.spl",  "22_08.spl")	
		#os.symlink("../22_22.spl",  "22_22.spl")
		os.symlink("../collisions.IN",  "collisions.IN")
		os.symlink("../equilibration.IN",  "equilibration.IN")
		os.symlink("../meam_parameters.IN",  "meam_parameters.IN")
		os.symlink("../lbomd.IN", "lbomd.IN")	
		os.symlink("../potfor.IN", "potfor.IN")
		os.symlink("../lattice.IN", "lattice.IN")
		os.symlink("../dpmta.IN", "dpmta.IN")
		os.system(".././LBOMD.exe")
		os.chdir("../")
		os._exit(0)	
		
	lattice = read_lattice_output(lattice,tempfolder)
	#lattice = read_lattice_output(lattice)
	#os.system("mv *.xyz ../")
	#os.chdir("../")
	os.system("rm -rf " + tempfolder)
	forcecalls+=1
	return lattice


def get_unique_folder_name():
	random.seed(time.time()*float(os.getpid()))
	#random.seed(time.time())
	tempfolder = "Temp" + str(random.randrange(0,100000))
	tempfolder = rootfolder + tempfolder
	
	while(os.path.isdir(tempfolder)):
		tempfolder = rootfolder+"Temp" + str(random.randrange(0,100000))
				
	os.makedirs(tempfolder)	
	
	return tempfolder

def run_dlpoly(folder,lattice):

	
	cwd = os.getcwd()
#	child = os.fork()
#	if child:
#		os.waitpid(child,0)
#	else:	
#		os.chdir(tempfolder)
#		write_dlpoly_FIELD(lattice)
#		printl("MAKING SYMLINKS", os.getcwd())
#		#make symlinks to all the input files
#		os.symlink(cwd+"/CONTROL",  "CONTROL")
#		#os.symlink(cwd+"/FIELD",  "FIELD")
#		os.symlink(cwd+"/TABLE",  "TABLE")
#		os.symlink(cwd+"/DLPOLY.Y",  "DLPOLY.Y")
#		os.system("./DLPOLY.Y")
#		#os.system(cwd+"/./DLPOLY.Y")
#		#os.chdir(cwd+"/../")
#		os._exit(0)	
		
	os.chdir(folder)
	try:os.symlink(cwd+"/CONTROL",  "CONTROL")
	except:pass#if CONTROL exists must have mean manually created, do not symlink
	try:os.symlink(cwd+"/TABLE",  "TABLE")
	except:pass
	try:os.symlink(cwd+"/DLPOLY.Y",  "DLPOLY.Y")
	except:pass
	
	printd("RUNNING DLPOLY IN CWD",os.getcwd())
	write_dlpoly_FIELD(lattice)
	write_lattice(lattice, "CONFIG")	
	os.system("./DLPOLY.Y")
	
	os.chdir(cwd)
	return lattice
	
def eval_lattice_dlpoly(lattice):
	global forcecalls

	tempfolder = get_unique_folder_name()
	lattice = run_dlpoly(tempfolder,lattice)
	lattice = read_lattice_output(lattice,tempfolder)	
	os.system("rm -rf " + tempfolder)
	
	
	
	forcecalls+=1
	return lattice

def eval_lattice_artmd_buck(lattice):
	global forcecalls, pid
	#not sure of the best way to proceed here, but for now will work in a random directory :S bit unhappy because of lots of io, but it makes the code cleaner
#	random.seed(time.time())
#	tempfolder = "Temp" + str(random.randrange(0,100000))
#	tempfolder = rootfolder + tempfolder + pid

	
#	while(os.path.isdir(tempfolder)):
#		tempfolder = rootfolder + "Temp" + str(random.randrange(0,100000) + pid)		

	tempfolder = _get_temp_path()

#	os.makedirs(tempfolder)
	write_lattice_artmd_buck(lattice, tempfolder + "/tadmd1.before.dat")
	cwd = os.getcwd()
	child = os.fork()
	if child:
		os.waitpid(child,0)
	else:	
		os.chdir(tempfolder)
#		make symlinks to all the input files
		os.symlink(cwd+"/buck.parms",  "buck.parms")
		os.symlink(cwd+"/buckmorse.parms",  "buckmorse.parms")
		os.system(cwd + "/./tadmd1.exe >> /dev/null")
		os._exit(0)	


		
	lattice = read_lattice_output(lattice,tempfolder)
	#lattice = read_lattice_output(lattice)
	#os.system("mv *.xyz ../")
	#os.chdir("../")
	os.system("rm -rf " + tempfolder)
	forcecalls+=1
	return lattice


def eval_lattice_lammps(lattice):
	global forcecalls
	random.seed(time.time())
	tempfolder = "Temp" + str(random.randrange(0,100000))
	
	while(os.path.isdir(tempfolder)):
		tempfolder = "Temp" + str(random.randrange(0,100000))		
	os.makedirs(tempfolder)
	child = os.fork()
	if child:
		os.waitpid(child,0)
	else:	
		os.chdir(tempfolder)
		#os.symlink("../Ag_u3.eam", "Ag_u3.eam")
		for filename in lattice.lammps_input_files:
		  os.symlink("../"+filename, filename)
		write_lattice_lammps(lattice,lattice.lammps_lattice_name)
		os.system(".././lmp < " + lattice.lammps_jobfile + " >/dev/null")
		os._exit(0)		
	
	lattice = read_lattice_output(lattice,tempfolder)
	os.system("rm -rf " + tempfolder)
	forcecalls+=1
	return lattice
	
def eval_lattice_lammpy(lattice):
	global forcecalls
	
	#lattice = putil.wrap_periodic(lattice)
	try: lattice.lmp
	except:
		import lammps
		lattice.lmp = lammps.lammps(cmdargs=copy.copy(lattice.lmp_args))
		for command in lattice.lammps_jobfile_array:
			if(command == "read_data"):
				path = _get_temp_path()
				write_lattice_lammps(lattice, path+"lattice.lammps")
				lattice.lmp.command("read_data " + path+"lattice.lammps")
			else:
				lattice.lmp.command(command)
	


	if(lattice.lmp.get_natoms() != lattice.NAtoms):
		print "Raising exception due to lammps losing atoms."
		raise NameError('LAMMPS has lost atoms')		

	x = lattice.lmp.extract_atom("x", 3)	
	lattice.lmp.scatter_atoms("x",1,3,lattice.Pos.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))
	
	#print "Performing Force Call"
	lattice.lmp.command("run 0")
	f = lattice.lmp.extract_atom("f", 3)
	lattice.Force[0:lattice.NAtoms*3] = f[0].__getslice__(0, lattice.NAtoms*3)		
	lattice.TPE = lattice.lmp.extract_compute("thermo_pe",0,0)
	lattice.PE[0:lattice.NAtoms] = numpy.array(lattice.lmp.extract_compute("atom_pe",1,1).__getslice__(0, lattice.NAtoms))
	#x = lattice.lmp.extract_atom("x", 3)
	#lattice.Pos[0:lattice.NAtoms*3] = x[0].__getslice__(0, lattice.NAtoms*3)
	
	forcecalls+=1
	return lattice	
	
	
	
	
	
	
	
	

def eval_lattice_lbomd_piping(lattice):
	global forcecalls
	DimString = str(lattice.Dim[0]) + ' ' + str(lattice.Dim[1])  + ' ' + str(lattice.Dim[2]) + '\n'
		
	dir=os.getcwd()
	random.seed(time.time())
	id = random.randint(0,99999)
	fullid = "%04d" % (id)
	fifoname = "latticepipe"+str(fullid)
	#print "FIFO",dir+str(fifoname)
	while (os.path.exists(dir+"/"+str(fifoname))):
		print "FIFO",dir+"/"+str(fifoname ), "EXISTS"
		id = random.randint(0,99999)
		fullid = "%04d" % (id)
		fifoname = "latticepipe"+str(fullid)       
		
	os.mkfifo(fifoname)

	if(verbose):print "FORKING TO PIPES",fifoname
	child = os.fork()
	if child:
		#Father
		command  = "echo "+str(fifoname)+" | ./LBOMD.exe"
		output = commands.getoutput(command)
	
		array = output.split('ENERGYXX')
		emergyarray= array[1].split('\n')
		TOTALENERGY = emergyarray[0];
		arrayb = array[0].split('PIPEXX')
		#print "ARRAY 0 @arrayb\n";

		count=0
		Atomcount = 0
		magF = 0
		MaxForce = -10000
	
		for line in arrayb:
			if(count>0):
				arrayc = line.split()				
		#		lattice.Pos[3*Atomcount] = float(arrayc[2])
		#		lattice.Pos[3*Atomcount+1] = float(arrayc[3])
		#		lattice.Pos[3*Atomcount+2] = float(arrayc[4])	
				lattice.KE[Atomcount] = float(arrayc[5])
				try:	
					lattice.PE[Atomcount] = float(arrayc[6])
				except:
					lattice.PE[Atomcount] = 99999
				lattice.Force[3 * Atomcount + 0] = float(arrayc[7])
				lattice.Force[3 * Atomcount + 1] = float(arrayc[8])
				lattice.Force[3 * Atomcount + 2] = float(arrayc[9])
				lattice.Charge[Atomcount] = float(arrayc[10])
				
#				magF = putil.magnitude((lattice.Force[3 * Atomcount + 0],lattice.Force[3 * Atomcount + 1],lattice.Force[3 * Atomcount + 2]))
#				if(magF > MaxForce):
#					MaxForce = magF
				Atomcount+=1
			count+=1
#		lattice.MaxForce=MaxForce	
		os.unlink(str(fifoname))
		#os.system("rm -rf "+str(fifoname))
		os.waitpid(child,0)
			
	else:
		#Child
		#fifo = open(fifoname, 'w')
		#NEED os.open and os.O_WRONLY for pipe write only when there's a reader.
		fifo=os.open(fifoname, os.O_WRONLY)
		os.write(fifo,str(lattice.NAtoms)+"\n")
		os.write(fifo,str(DimString))
		for i in numpy.arange(0,lattice.NAtoms):
			string = str(lattice.Specie[i])+" "+str(lattice.Pos[i*3])+" "+str(lattice.Pos[i*3+1])+" "+str(lattice.Pos[i*3+2])+" "+str(lattice.Charge[i])+" \n"
			os.write(fifo,string)
		os.close(fifo)
		os._exit(0)
		
	
	forcecalls+=1			
	return lattice


def read_lattice_output(lattice,dir='./'):
	if(mdtype == "LBOMD"):
		lattice = read_lattice_output_lbomd(lattice, dir)
	elif(mdtype == "LAMMPS"):
		lattice = read_lattice_output_lammps(lattice, dir)
	elif(mdtype == "ARTMD_BUCK"):
		lattice = read_lattice_output_artmd_buck(lattice, dir)
	elif(mdtype == "DLPOLY"):
		lattice = read_lattice_output_dlpoly(lattice, dir)	
	else:
		print "read_lattice_output: Environmental variable MDTYPE not recognised or set."
		sys.exit()
		
	return lattice  

def read_lattice_output_dlpoly(lattice,dir='./'):
	try:
		f = open(str(dir)+"/STATIS", "r")
	except:
		printl("Error: read_lattice_output_dlpoly: must output energy to STATIS",str(dir)+"/STATIS")
	f.readline()
	f.readline()
	f.readline()
	temparray = f.readline().split()
	totalpe = float(temparray[0])
	lattice.TPE = totalpe
	#print temparray
	#print "ENERGY",totalpe
	f.close()
	

	f = open(str(dir)+"/HISTORY", "r")
	f.readline()
	temparray = f.readline().split()
	forcecheck = int(temparray[0])
	if(forcecheck!=2):
		print "Error: read_lattice_output_dlpoly: must output forces to HISTORY"
		sys.exit()
		
	f.readline()
	f.readline()
	f.readline()
	f.readline()
#	MaxForce = -1000
	for i in range(0,lattice.NAtoms):
		temparray = f.readline().split()
		lattice.Charge[i] = float(temparray[3])
		f.readline()
		f.readline()
		temparray = f.readline().split()
		lattice.PE[i] = (totalpe)/float(lattice.NAtoms)
		lattice.Force[3*i] = float(temparray[0])/9648.530821
		lattice.Force[3*i+1] = float(temparray[1])/9648.530821
		lattice.Force[3*i+2] = float(temparray[2])/9648.530821
		#print i,lattice.PE[i],lattice.Force[3*i],lattice.Force[3*i+1],lattice.Force[3*i+2]
#		magF = putil.magnitude((lattice.Force[3 * i + 0],lattice.Force[3 * i + 1],lattice.Force[3 * i + 2]))
#		if(magF > MaxForce):
#			MaxForce = magF
#	lattice.MaxForce=MaxForce
	

	return lattice	
	



	
def read_lattice_output_lbomd(lattice,dir='./'):
	f = open(str(dir)+"/animation-reference.xyz", "r")
	f.readline()
	f.readline()
#	MaxForce = -10000
	for i in range(0,lattice.NAtoms):
		temparray = f.readline().split()
	#	lattice.Pos[3*i] = float(temparray[2])
	#	lattice.Pos[3*i+1] = float(temparray[3])
	#	lattice.Pos[3*i+2] = float(temparray[4])
	#	lattice.KE[i] = float(temparray[5])
		try:	
			lattice.PE[i] = float(temparray[6])
		except:
			lattice.PE[i] = 99999	
		lattice.Force[3*i] = float(temparray[7])
		lattice.Force[3*i+1] = float(temparray[8])
		lattice.Force[3*i+2] = float(temparray[9])	
		lattice.Charge[i] = float(temparray[10])
#		magF = magnitude((lattice.Force[3 * i + 0],lattice.Force[3 * i + 1],lattice.Force[3 * i + 2]))
#		if(magF > MaxForce):
#			MaxForce = magF
#	lattice.MaxForce=MaxForce
	

	return lattice
	
	
def read_lattice_output_artmd_buck(lattice,dir='./'):
	f = open(str(dir)+"/dump.forces", "r")
	f.readline()
	f.readline()
	MaxForce = -10000
	lattice.TPE = 0
	for i in range(0,lattice.NAtoms):
		temparray = f.readline().split()
	#	lattice.Pos[3*i] = float(temparray[2])
	#	lattice.Pos[3*i+1] = float(temparray[3])
	#	lattice.Pos[3*i+2] = float(temparray[4])
	#	lattice.KE[i] = 0
		lattice.PE[i] = float(temparray[8])
		lattice.TPE += lattice.PE[i]

		lattice.Force[3*i] = -float(temparray[5])
		lattice.Force[3*i+1] = -float(temparray[6])
		lattice.Force[3*i+2] = -float(temparray[7])	

#		magF = magnitude((lattice.Force[3 * i + 0],lattice.Force[3 * i + 1],lattice.Force[3 * i + 2]))
#		if(magF > MaxForce):
#			MaxForce = magF
#	lattice.MaxForce=MaxForce


	return lattice	
	
	
def read_lattice_output_lammps(lattice,dir='./'):
	f = open(str(dir)+"/dump.forces", "r")
#	MaxForce = -10000
	while(True):
	  line = f.readline()
	
	  if not line:
	    break
	  line_array = line.rstrip().split()
	  if len(line_array) > 1:
	    if line_array[1] == "ATOMS":
	      for j in range(0,lattice.NAtoms):
		temparray = f.readline().rstrip().split()
		i = int(temparray[0]) - 1
	#	lattice.Pos[3*i] = float(temparray[2])
	#	lattice.Pos[3*i+1] = float(temparray[3])
	#	lattice.Pos[3*i+2] = float(temparray[4])
	#	lattice.KE[i] = 0
		lattice.PE[i] = float(temparray[5])
		
		lattice.Force[3*i] = float(temparray[6])
		lattice.Force[3*i+1] = float(temparray[7])
		lattice.Force[3*i+2] = float(temparray[8])	
		if("charge" in lattice.lammps_atom_style):
		  lattice.Charge[i] = float(temparray[9])
		
#		magF = magnitude((lattice.Force[3 * i + 0],lattice.Force[3 * i + 1],lattice.Force[3 * i + 2]))
#		if(magF > MaxForce):
#			MaxForce = magF
	
			
		  
#	lattice.MaxForce=MaxForce

			
	return lattice
	    
	
	
def read_lattice_output_old(lattice):
	f = open("animation-reference.xyz", "r")
	f.readline()
	f.readline()
	MaxForce = -10000
	for i in range(0,lattice.NAtoms):
		temparray = f.readline().split()
		lattice.Pos[3*i] = float(temparray[2])
		lattice.Pos[3*i+1] = float(temparray[3])
		lattice.Pos[3*i+2] = float(temparray[4])
		lattice.KE[i] = float(temparray[5])
		try:	
			lattice.PE[i] = float(temparray[6])
		except:
			lattice.PE[i] = 99999	
		lattice.Force[3*i] = float(temparray[7])
		lattice.Force[3*i+1] = float(temparray[8])
		lattice.Force[3*i+2] = float(temparray[9])	
		lattice.Charge[i] = float(temparray[10])
		magF = putil.magnitude((lattice.Force[3 * i + 0],lattice.Force[3 * i + 1],lattice.Force[3 * i + 2]))
		if(magF > MaxForce):
			MaxForce = magF
	lattice.MaxForce=MaxForce		
	return lattice
		
def vis_lattice2(lattice, prefix,folder=None,suffix=None,num=None):
	if(folder):folder = folder+"/"
	
	if(vistype == "ATOMEYE"):
		vis_lattice_atomeye(lattice, folder+prefix+str(num)+suffix)
	elif(vistype == "XYZ"):
		vis_lattice_xyz(lattice, folder+prefix+str(num)+suffix)		
	elif(vistype == "LBOMD"):
		vis_lattice_lbomd_xyz(lattice,folder,prefix,suffix,num)	
	else:
		print "vis_lattice_output: Environmental variable VIS not recognised or set."
		return

def vis_lattice(orig_lattice, filename, group=0):
	
	lattice = putil.copylattice(orig_lattice)
	#try:
	#	lattice.Dim[0] += lattice.Dim[3]
	#	lattice.Dim[1] += lattice.Dim[4]
	#	lattice.Dim[2] += lattice.Dim[5]
	#	lattice.Pos[0::3] += lattice.Dim[3]
	#	lattice.Pos[1::3] += lattice.Dim[4]
	#	lattice.Pos[2::3] += lattice.Dim[5]
	#except: pass
	
	if(vistype == "ATOMEYE"):
		vis_lattice_atomeye(lattice, filename)
	elif(vistype == "XYZ"):
		vis_lattice_xyz(lattice, filename)		
	elif(vistype == "LBOMD"):
		#can only pass as lbomd needs prefix/number otherwise it gets messy (see vis_lattice2)
		pass
	elif(vistype == "LAMMPS"):
		vis_lattice_lammps(lattice, filename, group)
	else:
		print "vis_lattice_output: Environmental variable VIS not recognised or set."
		return
	 
	
def vis_lattice_atomeye(lattice, filename):
	f=open(filename, "w")
	f.write("Number of particles = " + str(lattice.NAtoms) + "\n")
	f.write("A = 5 Angstrom\n")
	f.write("H0(1,1) = " + str(lattice.Dim[0]) + " A\n")
	f.write("H0(1,2) = 0 A\nH0(1,3) = 0 A\nH0(2,1) = 0 A\n")
	f.write("H0(2,2) = " + str(lattice.Dim[1]) + " A\n")	
	f.write("H0(2,3) = 0 A\nH0(3,1) = 0 A\nH0(3,2) = 0 A\n")
	f.write("H0(3,3) = " + str(lattice.Dim[2]) + " A\n")	
	
	f.write(".NO_VELOCITY.\n")
	
	f.write("entry_count = 7\n")
	
	f.write("auxiliary[0] = fx [eV/A]\n")
	f.write("auxiliary[1] = fy [eV/A]\n")
	f.write("auxiliary[2] = fz [eV/A]\n")		
	f.write("auxiliary[0] = pe [eV]\n")	
	
	for i in range(0, lattice.NAtoms):
		f.write("1.00\n")
		f.write(str(lattice.species_info[int(lattice.Specie[i])-1].Label) + "\n")
		f.write(str(lattice.Pos[3*i]/lattice.Dim[0]) + " " + str(lattice.Pos[3*i+1]/lattice.Dim[1]) + " " + str(lattice.Pos[3*i+2]/lattice.Dim[2]) + " " + str(lattice.PE[i]) + " " + str(lattice.Force[3*i+0]) + " " + str(lattice.Force[3*i+1]) + " " + str(lattice.Force[3*i+2]) + "\n")
	
	f.close()
	
def vis_lattice_xyz(lattice, filename):
	f=open(filename, "w")
	f.write(str(lattice.NAtoms)+"\n")
	f.write("PyKMC Automatically Generated XYZ File\n")
	for i in range(0, lattice.NAtoms):
		f.write(str(lattice.Specie[i]) + " " + str(lattice.Pos[3*i])+ " " + str(lattice.Pos[3*i+1])+ " " + str(lattice.Pos[3*i+2])+"\n")
	f.close()

	
def vis_lattice_lbomd_xyz(lattice,folder,prefix,suffix,num,time=0.0):	
	#print "IN VIS LBOMD"
	number = "%04d" % (num)
	filename = folder+str(prefix[0:5])+str(number)+suffix
	f = open(filename, 'w')
	f.write(str(lattice.NAtoms)+'\n')
	string = str(time) + '\n'
	f.write(string)
	for i in range(0,lattice.NAtoms):
		string = str(i)+" %14.7g" % lattice.Pos[3*i] + " %14.7g" % lattice.Pos[3*i+1] + " %14.7g" % lattice.Pos[3*i+2] + " %14.7g" % lattice.KE[i] +  "%14.7g" % lattice.PE[i] + "\n"
		f.write(string)
	f.close()	
	return
	
def vis_lattice_lammps(lattice, filename, group=0):
	f = open(filename, 'w')
	lattice = putil.copylattice(lattice)
	f.write("ITEM: NUMBER OF ATOMS\n")
	if(group!=0):
		f.write(str(len(lattice.defect_list_group[group]))+"\n")
	else:
		f.write(str(lattice.NAtoms)+"\n")
	f.write("ITEM: BOX BOUNDS ")
	string = ""
	if(lattice.PBC[0] == 1):
		string += "pp "
	else:
		string += "ff "
	if(lattice.PBC[1] == 1):
		string += "pp "
	else:
		string += "ff "
	if(lattice.PBC[2] == 1):
		string += "pp "
	else:
		string += "ff "
	f.write(string + "\n")
	f.write(str(lattice.Dim[3]) + " " + str(lattice.Dim[0]) + "\n")
	f.write(str(lattice.Dim[4]) + " " + str(lattice.Dim[1]) + "\n")
	f.write(str(lattice.Dim[5]) + " " + str(lattice.Dim[2]) + "\n")
	
	
	
	f.write("ITEM: ATOMS id type x y z pe\n")
	if(group!=0):
		count = 0
		print lattice.defect_list_group[group]
		for i in lattice.defect_list_group[group]:
			count += 1
			string = str(count)+ " " + lattice.Specie[i] + " %14.7g" % lattice.Pos[3*i] + " %14.7g" % lattice.Pos[3*i+1] + " %14.7g" % lattice.Pos[3*i+2] + "%14.7g" % lattice.PE[i] + "\n"
			f.write(string)
	else:
		for i in range(0,lattice.NAtoms):
			string = str(i+1)+ " " + lattice.Specie[i] + " %14.7g" % lattice.Pos[3*i] + " %14.7g" % lattice.Pos[3*i+1] + " %14.7g" % lattice.Pos[3*i+2] + "%14.7g" % lattice.PE[i] + "\n"
			#print string
			f.write(string)
	f.close()
		


def write_dlpoly_FIELD(lattice):
	#print "DLPOLY FIELD",lattice.species_info
	
	f = open("FIELD","w")
	
	f.write(lattice.JobName+"\n")
	f.write("UNITS eV \n")
	f.write("MOLECULES "+str(lattice.NSpecies)+"\n")
	#for m in range(0,lattice.NSpecies):
	for key in lattice.species_info:
		#print "FIELD",key
		if(lattice.species_info[key].AtomsPerSpec >0):
			f.write(lattice.species_info[key].Label+" \n")	
			f.write("NUMMOLS "+str(lattice.species_info[key].AtomsPerSpec)+" \n")	
			f.write("ATOMS 1 \n")	
			f.write(str(key)+" "+str(lattice.species_info[key].Mass)+" "+str(lattice.species_info[key].Charge)+" 1 0 \n")	
			f.write("FINISH \n")	
			f.write("\n")
	
	for line in lattice.InteractionLine:
		f.write(line)		

	f.close()	
	
def write_dlpoly_FIELD(lattice):
	#print "DLPOLY FIELD",lattice.species_info
	
	f = open("FIELD","w")
	
	f.write(lattice.JobName+"\n")
	f.write("UNITS eV \n")
	f.write("MOLECULES "+str(lattice.NSpecies)+"\n")
	#for m in range(0,lattice.NSpecies):
	for key in lattice.species_info:
		#print "FIELD",key
		if(lattice.species_info[key].AtomsPerSpec >0):
			f.write(lattice.species_info[key].Label+" \n")	
			f.write("NUMMOLS "+str(lattice.species_info[key].AtomsPerSpec)+" \n")	
			f.write("ATOMS 1 \n")	
			f.write(str(key)+" "+str(lattice.species_info[key].Mass)+" "+str(lattice.species_info[key].Charge)+" 1 0 \n")	
			f.write("FINISH \n")	
			f.write("\n")
	
	for line in lattice.InteractionLine:
		f.write(line)		

	f.close()	
	
def write_dlpoly_CONTROL(jobname,temp,press,ensemble,steps,equilsteps,traj,timestep,cutoff,
						pka_specie=None,pka_loc=None,pka_id=None,lattice=None,pka_energy=None,
						pka_direction=None,thermallayer=None,folder=None,print_stats=50,
						restart=False):
	collision = False
	if(pka_specie or pka_loc or pka_id):
		collision = True
		if(pka_id):
			pass
		else:
			find_pka(lattice,pka_specie,pka_loc)
	if(folder):	f = open(str(folder)+"/CONTROL","w")
	else:f = open("CONTROL","w")
	
	f.write(jobname+"\n")
	if(restart):f.write("restart \n")
	f.write("integrator velocity \n")
	f.write("temperature      "+str(temp)+"\n")
	f.write("pressure         "+str(press)+"\n")
	f.write("ensemble "+str(ensemble)+"   \n")
	if(thermallayer):f.write(str(thermallayer)+"\n")
	if(collision):
		f.write("impact "+str(pka_id+1)+" 0 "+str(float(pka_energy)/1000.0)+" "+str(pka_direction[0])+" "+str(pka_direction[1])+" "+str(pka_direction[2])+"  \n")
	
	f.write("steps               "+str(steps)+"\n")
	f.write("equilibration       "+str(equilsteps)+"\n")
	
	f.write("traj "+str(traj)+"\n")
	f.write("print "+str(print_stats)+"\n")
	f.write("stack "+str(print_stats)+"\n")
	f.write("stats "+str(print_stats)+"\n")
	
	f.write("variable timestep "+str(timestep)+"\n")
	f.write("cutoff "+str(cutoff)+"\n")
	f.write("ewald precision  1.0E-6\n")
	
	f.write("job time   62000.0\n")
	f.write("close time    2000.00\n")
	
	f.write("finish \n")	
	
def get_info_from_STATIS(step,rundir=None):
	printl("GETTING INFO FROM STATIS",step,rundir)
	if(rundir):f = open(str(rundir)+"/STATIS","r")
	else:f = open("STATIS","r")
	f.readline()
	f.readline()
    
#    for s in range(0,step-1):
#    		print
#            f.readline()
	line = "blah"
	while(1 and len(line)>0):        
		line = f.readline().split()
		#printl("STATIS",line)
		try:
			stepin = int(line[0])
			time = float(line[1])
			if(stepin == step):
				line = f.readline().split()
				energy = float(line[0])
				break
		except:
			pass	
	#        for i in range(0,8):
	#            f.readline()
	
	printl("STEP:",step,"TIME:",time, "DIR:",rundir)
	    
	return time,energy	





	
	
