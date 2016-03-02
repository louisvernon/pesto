#!/usr/bin/python


import os,sys,math, copy
import subprocess
import numpy
import time
from pesto.pio import *	
from pesto.putil import *
import k3match
#import pynauty
try:
	import networkx
	import pygraphviz
	import matplotlib.pyplot as plt
except:	pass
import scipy

class Defect:
	pass
		

def identify_defects_classical(input_lattice, reference_lattice, vacancy_radius = 1, defect_radius=5):
	""" a faster version of the original defect module, without spatial decom in the traditional sense
	 uses a tree algorithm to rapidly find atomic neighbours (k3match) with some  skinning in order to handle pbcs
	 unfortunately this algorithm doesn't guarentee 1-1 mapping between the input and reference lattice, strictly using a cutoff to determine defect status
	 as such there likely won't be a correct total number of vacancies / interstitials
	 however it will produce reasonable points for the defect cluster """
	
	t = time.time()

	
	ix = input_lattice.Pos[0::3]
	iy = input_lattice.Pos[1::3]
	iz = input_lattice.Pos[2::3]
	
	rx = reference_lattice.Pos[0::3]
	ry = reference_lattice.Pos[1::3]
	rz = reference_lattice.Pos[2::3]
	
	if(input_lattice.PBC[0]):
		
		key_index = numpy.arange(0, input_lattice.NAtoms)
		
		x_skin_l = numpy.where(ix < defect_radius)[0]
		x_skin_u = numpy.where(ix > input_lattice.Dim[0] - defect_radius)[0]
		ix = numpy.append(ix, ix[x_skin_l]+input_lattice.Dim[0])
		iy = numpy.append(iy, iy[x_skin_l])
		iz = numpy.append(iz, iz[x_skin_l])
		key_index = numpy.append(key_index, key_index[x_skin_l])
		
		#x_skin = numpy.where(ix > input_lattice.Dim[0] - defect_radius)[0]
		ix = numpy.append(ix, ix[x_skin_u]-input_lattice.Dim[0])
		iy = numpy.append(iy, iy[x_skin_u])
		iz = numpy.append(iz, iz[x_skin_u])
		key_index = numpy.append(key_index, key_index[x_skin_u])
		
		
		
		y_skin_l = numpy.where(iy < defect_radius)[0]
		y_skin_u = numpy.where(iy > input_lattice.Dim[1] - defect_radius)[0]
		iy = numpy.append(iy, iy[y_skin_l]+input_lattice.Dim[1])
		ix = numpy.append(ix, ix[y_skin_l])
		iz = numpy.append(iz, iz[y_skin_l])
		key_index = numpy.append(key_index, key_index[y_skin_l])
		
		#y_skin = numpy.where(iy > input_lattice.Dim[1] - defect_radius)[0]
		iy = numpy.append(iy, iy[y_skin_u]-input_lattice.Dim[1])
		ix = numpy.append(ix, ix[y_skin_u])
		iz = numpy.append(iz, iz[y_skin_u])
		key_index = numpy.append(key_index, key_index[y_skin_u])

					
		z_skin_l = numpy.where(iz < defect_radius)[0]
		z_skin_u = numpy.where(iz > input_lattice.Dim[2] - defect_radius)[0]
		iz = numpy.append(iz, iz[z_skin_l]+input_lattice.Dim[2])
		ix = numpy.append(ix, ix[z_skin_l])
		iy = numpy.append(iy, iy[z_skin_l])
		key_index = numpy.append(key_index, key_index[z_skin_l])
		
		#z_skin = numpy.where(iz > input_lattice.Dim[2] - defect_radius)[0]
		iz = numpy.append(iz, iz[z_skin_u]-input_lattice.Dim[2])
		ix = numpy.append(ix, ix[z_skin_u])
		iy = numpy.append(iy, iy[z_skin_u])
		key_index = numpy.append(key_index, key_index[z_skin_u])	
			
					
					
		#print "PBC Skin. Total added atoms:", len(key_index) - input_lattice.NAtoms
						
		
	result = k3match.cartesian(ix,iy,iz, rx, ry,rz, vacancy_radius)		
	
	input_map = result[0]
	ref_map = result[1]
	sep = result[2]
	
	
	# map the results back to their original ids
	input_map = key_index[input_map]
	
	#print "K3match time was:", time.time()-t
	t = time.time()
	
	interstitials = numpy.setdiff1d( numpy.arange(input_lattice.NAtoms),input_map)
	vacancies = numpy.setdiff1d(numpy.arange(reference_lattice.NAtoms), ref_map)
	print "Interstitials:", len(interstitials), interstitials
	print "Vacancies:", len(vacancies), vacancies
	
	input_lattice.Interstitials = len(interstitials)
	input_lattice.Vacancies = len(vacancies)
	
	input_lattice.defect_points = []
	rx = []
	ry = []
	rz = []
	for i in interstitials:
		#if(mask_filter and input_lattice.global_DefectMask[3*i]==0):
		#	continue
		input_lattice.defect_points.append(input_lattice.Pos[3*i:3*i+3])
		rx.append(input_lattice.defect_points[-1][0])
		ry.append(input_lattice.defect_points[-1][1])
		rz.append(input_lattice.defect_points[-1][2])
	for i in vacancies:
		input_lattice.defect_points.append(reference_lattice.Pos[3*i:3*i+3])
		rx.append(input_lattice.defect_points[-1][0])
		ry.append(input_lattice.defect_points[-1][1])
		rz.append(input_lattice.defect_points[-1][2])
	
	
	# build the list of point list neighbours using k3match again.
	result = k3match.cartesian(ix,iy,iz, rx, ry,rz, defect_radius)		
	
	input_map = result[0]
	ref_map = result[1]
	sep = result[2]
	
	
	# map the results back to their original ids
	input_map = key_index[input_map]
	
	sort_key = numpy.argsort(input_map)
	input_map = input_map[sort_key]
	ref_map = ref_map[sort_key]
	
	
	# find common atoms:
	# create initial group assignments (each source defect atom is its own group)
	groups = numpy.arange(len(input_lattice.defect_points))
	
	for i in xrange(0, len(input_map)-1):
		if input_map[i] == input_map[i+1]:
			# groups intersect
			a = ref_map[i]
			b = ref_map[i+1]
			
			if(groups[a]<groups[b]):
				#print a, b, numpy.where(groups==groups[a]), numpy.where(groups==groups[b])
				old_group = groups[b]
				groups[numpy.where(groups==groups[b])[0]] = groups[a]
				groups[numpy.where(groups>old_group)[0]] -= 1
				#groups[b] = groups[a]
			elif(groups[b]<groups[a]):
				#print a, b, numpy.where(groups==groups[a]), numpy.where(groups==groups[b])					
				old_group = groups[a]
				groups[numpy.where(groups==groups[a])[0]] = groups[b]
				groups[numpy.where(groups>old_group)[0]] -= 1
				#groups[a] = groups[b]
				
					
		#print groups
		
	
	groups += 1

	
	input_lattice.defect_list, index = numpy.unique(input_map, return_index=True)
	input_map = input_map[index]
	ref_map = ref_map[index]
	input_lattice.NDefects = len(input_lattice.defect_list)
	
	chull = numpy.zeros((input_lattice.NDefects, 3))
	for i in xrange(0, len(input_lattice.defect_list)):
		chull[i] = input_lattice.Pos[3*i:3*(i+1)]
	#delaunay = scipy.spatial.Delaunay(chull)
	
	
	# set group_id for each defect atom
	input_lattice.group_id = groups[ref_map]
	input_lattice.NGroups = max(groups)
	print "Number of defect groups:", input_lattice.NGroups
	
	
	
	input_lattice.defect_list_group = []
	input_lattice.defect_list_group.append(input_lattice.defect_list)
	for i in xrange(1, input_lattice.NGroups+1):
		input_lattice.defect_list_group.append(input_map[numpy.where(groups[ref_map] == i)[0]])

	
	#print input_map[numpy.where(input_lattice.group_id == 4)[0]]
	defect_lattice = get_defect_lattice(input_lattice,0)
	#vis_lattice(defect_lattice, "Defect0.xyz")
	



	input_lattice.DefectMask = numpy.zeros((input_lattice.NAtoms*3), dtype=float)
	
	for i in xrange(0, input_lattice.NDefects):
		input_lattice.DefectMask[3*input_lattice.defect_list[i]] = 1
		input_lattice.DefectMask[3*input_lattice.defect_list[i]+1] = 1
		input_lattice.DefectMask[3*input_lattice.defect_list[i]+2] = 1			
	
	print "Number of atoms interacting with defects:", len(input_lattice.defect_list)
	
	
	#print "Python time was: ", time.time()-t
	return input_lattice				
		
		
		
		
		
def get_defect_lattice(lattice, group=0):
  defect_lattice = copylattice(lattice)
  
  if(group == 0):
	defect_lattice.Pos.resize(lattice.NDefects*3)
	try:
		defect_lattice.Charge.resize(lattice.NDefects)
	except: pass
	defect_lattice.Specie.resize(lattice.NDefects)
	for j in xrange(0, lattice.NDefects):
		i = lattice.defect_list[j]
		#print j, i
		defect_lattice.Pos[3*j] = lattice.Pos[3*i]
		defect_lattice.Pos[3*j+1	] = lattice.Pos[3*i+1]
		defect_lattice.Pos[3*j+2] = lattice.Pos[3*i+2]
		defect_lattice.Charge[j] = lattice.Charge[i]
		defect_lattice.Specie[j] = lattice.Specie[i]
		defect_lattice.NAtoms = lattice.NDefects
		try: defect_lattice.PE = lattice.group_id
		except: pass
  else:
  	defect_lattice.NAtoms = 0
	defects = numpy.sum(lattice.group_id==group)
	defect_lattice.Pos.resize(defects*3)
	try:
		defect_lattice.Charge.resize(defects)
	except: pass
	defect_lattice.Specie.resize(defects)	
	for j in xrange(0, lattice.NDefects):
		#print lattice.group_id[j]
		i = lattice.defect_list[j]
		if(lattice.group_id[j] == group):
			defect_lattice.Pos[3*defect_lattice.NAtoms] = lattice.Pos[3*i]
			defect_lattice.Pos[3*defect_lattice.NAtoms+1] = lattice.Pos[3*i+1]
			defect_lattice.Pos[3*defect_lattice.NAtoms+2] = lattice.Pos[3*i+2]
			defect_lattice.Charge[defect_lattice.NAtoms] = lattice.Charge[i]
			defect_lattice.Specie[defect_lattice.NAtoms] = lattice.Specie[i]
		
			defect_lattice.NAtoms+=1
	#print "Returning Group:", group, defect_lattice.NAtoms
	#print lattice.group_id
  return defect_lattice		
		
		
def get_defect_lattice_from_energy(input_lattice, defect_variation = 20, defect_radius = 0, sub_lattice=None):

	# try and establish bulk energies from a given lattice
	# this is complex because the lattice could have multiple species
	# bin the energies of each atom
	# I guess bins that have % or more atoms will be marked as bulk
	# deviations of energy by a certain percentage of those bins will be 
	# marked as defective
	
	print "Charge before:", numpy.sum(input_lattice.Charge)

	input_lattice.eval()
	print input_lattice.TPE
	threshold = 0.001  # threshold of atom numbers for marking bin as bulk
	target_bins = 100
	threshold = input_lattice.NAtoms*threshold
	
	max_energy = input_lattice.PE.max()
	min_energy = input_lattice.PE.min()
	bin_range = max_energy - min_energy
	bin_width = bin_range/target_bins # width in ev for the binning of atoms
 	print bin_width, max_energy, min_energy, threshold
	num_bins = int(bin_range/bin_width)+1
	
 	if (num_bins < 10):
 		print "Energy distribution not well defined, try a different bin width"
 		return input_lattice
	else: print num_bins, "bins"
 		
 	
 	bin_ids = numpy.zeros(input_lattice.NAtoms) - 1
 	
 	bins = numpy.zeros(num_bins, dtype=int)
 	for i in xrange(0, input_lattice.NAtoms):
 		if(sub_lattice):
 		#	print check_pos_in_sub_lattice(input_lattice.Pos[3*i:3*i+3], sub_lattice)
 			if(not check_pos_in_sub_lattice(input_lattice.Pos[3*i:3*i+3], sub_lattice)):
 				continue
 		index = int((input_lattice.PE[i] - min_energy) / bin_width)
 		bins[index] += 1
 		bin_ids[i] = index
 		
	#print threshold, bins 

	bad_bins = numpy.where(bins<threshold)[0]
	#print bad_bins, numpy.sum(bins[bad_bins])
 	defect_points = []
 	for bin in bad_bins:
 		#print bin, numpy.where(bin_ids==bin)[0]
 		defect_points.extend(numpy.where(bin_ids==bin)[0])

 	temp_lattice = copylattice(input_lattice)
 	temp_lattice.RemoveAtom(defect_points)
 	temp_lattice.eval()

 	defect_list = []
 	defect_list.extend(defect_points)
 	
 	print "Original Defects:", len(defect_list)
 	for i in xrange(0, temp_lattice.NAtoms):
 		delta_energy = math.fabs(temp_lattice.PE[i] - input_lattice.PE[temp_lattice.deletion_key[i]])
		#print delta_energy
 		if (delta_energy > defect_variation):
 			defect_list.append(temp_lattice.deletion_key[i])
			#input_lattice.Specie[temp_lattice.deletion_key[i]]="3"
			
 		
 	print "Total Defects:", len(defect_list)	
 	input_lattice.defect_list = defect_list
 	input_lattice.NDefects = len(input_lattice.defect_list)
 	
 	#input_lattice = get_defect_groups(input_lattice, defect_radius)
 	input_lattice.DefectMask = numpy.zeros((input_lattice.NAtoms*3), dtype=bool)

	for i in xrange(0, input_lattice.NDefects):
		input_lattice.DefectMask[3*input_lattice.defect_list[i]] = True
		input_lattice.DefectMask[3*input_lattice.defect_list[i]+1] = True
		input_lattice.DefectMask[3*input_lattice.defect_list[i]+2] = True	
	
	#input_lattice = group_defects(input_lattice, defect_radius)

	defect_lattice = copylattice(input_lattice)
	defect_lattice.Specie[input_lattice.defect_list]="d"
	vis_lattice(defect_lattice, "defectlattice.xyz")
	print "Charge after:", numpy.sum(input_lattice.Charge)
	return input_lattice		

def get_defect_lattice_from_specie(input_lattice, specie, defect_radius=3):
	input_lattice.defect_points = []
	specie_list = numpy.where(input_lattice.Specie==specie)[0]
	point_list = []
	for atom in specie_list:
		point_list.append([input_lattice.Pos[3*atom], input_lattice.Pos[3*atom+1], input_lattice.Pos[3*atom+2]])
	
	defect_atoms = get_list_nearest_neighbours(point_list, input_lattice, defect_radius)
	print defect_atoms
	input_lattice.defect_points = []
	for atom in defect_atoms:
		input_lattice.defect_points.append([input_lattice.Pos[3*atom], input_lattice.Pos[3*atom+1], input_lattice.Pos[3*atom+2]])
	
	input_lattice.DefectMask[:] = False
		
	input_lattice.DefectMask[3*defect_atoms] = True
	input_lattice.DefectMask[3*defect_atoms+1] = True
	input_lattice.DefectMask[3*defect_atoms+2] = True
	
	print numpy.sum(input_lattice.DefectMask)/3, "Defects."
	
	input_lattice.NGroups = 1
	
	input_lattice.group_id = numpy.zeros(input_lattice.NAtoms)
	input_lattice.group_id[defect_atoms] = 1
	input_lattice.defect_list_group = []
	input_lattice.defect_list_group.append(defect_atoms)
	input_lattice.defect_list_group.append(defect_atoms)
	
	return input_lattice

	
	

		

def get_list_nearest_neighbours(point_list, input_lattice, defect_radius=3):
	rx = input_lattice.Pos[0::3]
	ry = input_lattice.Pos[1::3]
	rz = input_lattice.Pos[2::3]
	ix = []
	iy = []
	iz = []
	
	for point in point_list:
		ix.append(point[0])
		iy.append(point[1])
		iz.append(point[2])
		
	ix = numpy.array(ix)
	iy = numpy.array(iy)
	iz = numpy.array(iz)
		

	if(input_lattice.PBC[0]):
		
		key_index = numpy.arange(0, len(ix))
		
		x_skin_l = numpy.where(ix < defect_radius)[0]
		x_skin_u = numpy.where(ix > input_lattice.Dim[0] - defect_radius)[0]
		ix = numpy.append(ix, ix[x_skin_l]+input_lattice.Dim[0])
		iy = numpy.append(iy, iy[x_skin_l])
		iz = numpy.append(iz, iz[x_skin_l])
		key_index = numpy.append(key_index, key_index[x_skin_l])
		
		#x_skin = numpy.where(ix > input_lattice.Dim[0] - defect_radius)[0]
		ix = numpy.append(ix, ix[x_skin_u]-input_lattice.Dim[0])
		iy = numpy.append(iy, iy[x_skin_u])
		iz = numpy.append(iz, iz[x_skin_u])
		key_index = numpy.append(key_index, key_index[x_skin_u])
		
		
		
		y_skin_l = numpy.where(iy < defect_radius)[0]
		y_skin_u = numpy.where(iy > input_lattice.Dim[1] - defect_radius)[0]
		iy = numpy.append(iy, iy[y_skin_l]+input_lattice.Dim[1])
		ix = numpy.append(ix, ix[y_skin_l])
		iz = numpy.append(iz, iz[y_skin_l])
		key_index = numpy.append(key_index, key_index[y_skin_l])
		
		#y_skin = numpy.where(iy > input_lattice.Dim[1] - defect_radius)[0]
		iy = numpy.append(iy, iy[y_skin_u]-input_lattice.Dim[1])
		ix = numpy.append(ix, ix[y_skin_u])
		iz = numpy.append(iz, iz[y_skin_u])
		key_index = numpy.append(key_index, key_index[y_skin_u])

					
		z_skin_l = numpy.where(iz < defect_radius)[0]
		z_skin_u = numpy.where(iz > input_lattice.Dim[2] - defect_radius)[0]
		iz = numpy.append(iz, iz[z_skin_l]+input_lattice.Dim[2])
		ix = numpy.append(ix, ix[z_skin_l])
		iy = numpy.append(iy, iy[z_skin_l])
		key_index = numpy.append(key_index, key_index[z_skin_l])
		
		#z_skin = numpy.where(iz > input_lattice.Dim[2] - defect_radius)[0]
		iz = numpy.append(iz, iz[z_skin_u]-input_lattice.Dim[2])
		ix = numpy.append(ix, ix[z_skin_u])
		iy = numpy.append(iy, iy[z_skin_u])
		key_index = numpy.append(key_index, key_index[z_skin_u])
		
	result = k3match.cartesian(ix,iy,iz, rx, ry,rz, defect_radius)		
	
	input_map = result[0]
	ref_map = result[1]
	sep = result[2]
	
	return ref_map
	

def get_defect_lattice_from_mask(defect_lattice):
	lattice = copylattice(defect_lattice)
	remove = []
	for i in xrange(0, lattice.NAtoms):
		if (not lattice.DefectMask[3*i]): remove.append(i)
	remove.reverse()
	#print "Deleting:", remove
	lattice.RemoveAtom(remove)
	lattice.Dim[0] = lattice.Pos[0::3].max()
	lattice.Dim[1] = lattice.Pos[1::3].max()
	lattice.Dim[2] = lattice.Pos[2::3].max()
	lattice.Dim[3] = lattice.Pos[0::3].min()
	lattice.Dim[4] = lattice.Pos[1::3].min()
	lattice.Dim[5] = lattice.Pos[2::3].min()
	return lattice

def classify_defects(lattice, group):
	# this needs to identify whether the passed defect is new or previously recognised
	# first we need to box our atoms (min, max, x, y, z)
	# we use the first atom in the defect to construct a relative set of displacement vectors
	# then we translate the positions so that the center is 0,0
	tc = time.time()
	

	try:
		graph_cutoff  = lattice.graph_cutoff
		partition_min = lattice.partition_cutoff
	except:
		graph_cutoff  = 5.415
		partition_cutoff = 3.61
		
	print "Graph Cutoff:", graph_cutoff
	
	defect = Defect()	
	ref_atom = 0
	defect.S = []
	defect.NDefectAtoms = len(lattice.defect_list_group[group])
	print "Group:", group, "NDefectAtoms:", defect.NDefectAtoms
	defect.center = numpy.array([0.,0.,0.])
	t1 = time.time()
	#initialise defect arrays
	defect.x = numpy.empty(defect.NDefectAtoms, dtype=float)
	defect.y = numpy.empty(defect.NDefectAtoms, dtype=float)
	defect.z = numpy.empty(defect.NDefectAtoms, dtype=float)
	defect.atom = numpy.empty(defect.NDefectAtoms, dtype=int)
	
	for i in xrange(0, defect.NDefectAtoms):
		if(i == 0):
			defect.x[0] = 0;
			defect.y[0] = 0;
			defect.z[0] = 0;
			atom_mapped = lattice.defect_list_group[group][i]
			defect.S.append(lattice.Specie[atom_mapped])
			defect.atom[0] = atom_mapped
		else:
			# get displacement relative to reference atom
			i_mapped = lattice.defect_list_group[group][i]
			#print i, i_mapped
		
			x_sep = lattice.Pos[3*i_mapped] - lattice.Pos[3*atom_mapped]
			if(lattice.PBC[0]==1):
		 	  x_sep -= numpy.rint(x_sep/lattice.Dim[0])*lattice.Dim[0]
				  
			y_sep = lattice.Pos[3*i_mapped+1] - lattice.Pos[3*atom_mapped+1]
			if(lattice.PBC[1]==1):
		 	  y_sep -= numpy.rint(y_sep/lattice.Dim[1])*lattice.Dim[1]			  
			  
			  
			z_sep = lattice.Pos[3*i_mapped+2] - lattice.Pos[3*atom_mapped+2]
			if(lattice.PBC[2]==1):
		 	  z_sep -= numpy.rint(z_sep/lattice.Dim[2])*lattice.Dim[2]
		 	  		  				  
			defect.x[i] = x_sep;
			defect.y[i] = y_sep;
			defect.z[i] = z_sep;
			defect.S.append(lattice.Specie[i_mapped])
		
			defect.center[0]+=x_sep
			defect.center[1]+=y_sep
			defect.center[2]+=z_sep								
		
			defect.atom[i] = i_mapped
			
	
	

	# ok so we've built a defect object with positions relative to the reference atom, now we shift so that it starts from 0,0,0
#	correction, we will shift it so that its COM sits at 0,0,0	

	defect.center[0] /= defect.NDefectAtoms
	defect.center[1] /= defect.NDefectAtoms
	defect.center[2] /= defect.NDefectAtoms	

	#print defect.center
	defect.real_center = defect.center

	
	for i in xrange(0, defect.NDefectAtoms):
		#print i

		defect.x[i] -= defect.center[0]
		defect.y[i] -= defect.center[1]
		defect.z[i] -= defect.center[2]

	defect.center = numpy.array([0.,0.,0.])
	

	# get sample defect center translated back into real space 
	
	defect.real_center = [lattice.Pos[3*defect.atom[0]] - defect.x[0], lattice.Pos[3*defect.atom[0]+1] - defect.y[0], lattice.Pos[3*defect.atom[0]+2] - defect.z[0]];

	# now we construct the nauty graph
	graph = "l=10000 \n c\n n=" +str(defect.NDefectAtoms)+ " g\n";

	
	size = []
	color = []
	label = {}
	
	# new idea to try and avoid an issue art described with equivalent defects that are not symettrically equivelant!
	# we now partition by summation of integer distances
	partition = {}
	
	t2 = time.time()		
	for i in xrange(0, defect.NDefectAtoms):
		#G.add_nodes_from(numpy.arange(defect.NDefects))
		#P.add_node(str(i),)
		neb_dist_count = 0
		for j in xrange(0, defect.NDefectAtoms):
			if(i==j): continue
			# get displacement relative to reference atom, no need for pbc correction as i stored the relative displacement, right?
		
			x_sep = defect.x[i] - defect.x[j]
			y_sep = defect.y[i] - defect.y[j]
			z_sep = defect.z[i] - defect.z[j]
			  
			sep = math.sqrt(x_sep*x_sep + y_sep*y_sep + z_sep*z_sep)
			#neb_dist_count += numpy.rint(sep/partition_cutoff)
			if(sep < graph_cutoff):
				graph += str(j) + " "
				#neb_dist_count += numpy.rint(sep)
				if(sep < partition_cutoff):
					neb_dist_count += 1
			
		if neb_dist_count in partition:
			partition[neb_dist_count].append(i)
		else:
			partition[neb_dist_count] = [i]
				
		graph += ";\n"
	# new technique to partition by total displacement

	print partition
#	print "Second loop time", time.time() - t2
	keys = sorted(partition.keys())
	colors = ["red","green","blue"]
	print "Number of partitions:", len(keys)
	print keys
	graph += "f=["
	for j in xrange(0, len(keys)):
		for k in xrange(0, len(partition[keys[j]])):
			graph += str(partition[keys[j]][k]) + ","
		graph += "|"
	graph += "]"
	graph += "\nx";
	
	defect.graph = graph
	graph += "\nz\n"
	
	
	fpopen = subprocess.Popen(["./dreadnaut"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	output = fpopen.communicate(graph)	
	result = output[0]
	
	result_lines = result.rstrip().split("\n")
	if(len(result_lines) <= 1):
		print "Can't run without dreadnaut!"
		sys.exit()
	defect.hash = result_lines[-1]
	print defect.hash
	print "Time to classify defect:", time.time() - tc
	defect.group_id = group

	return defect
	
	

def sub_lattice(lattice):
	pass
	
	
def create_sub_lattice(lattice, group = 0, force_thresh = 1E-2):
	# this function has been added to deal with larger lattices
	
	# First we delete the atoms that belong to the defect cluster
	print "Number of atoms before extraction:", lattice.NAtoms
	sub_lattice = copylattice(lattice)
	sub_lattice.DefectMask[:] = True
	sub_lattice.eval()
	
	
	# we make a reference of the original lattice forces
	Reference_Force = copy.copy(sub_lattice.Force)

	
	remove_id = []
	for i in xrange(0, len(lattice.defect_list_group[group])):
		index = lattice.defect_list_group[group][-1-i]
		remove_id.append(index)
	

	print remove_id
	local_index = numpy.array(remove_id)
	index = numpy.arange(lattice.NAtoms)
	index = numpy.delete(index, local_index)
	vis_sub_lattice = copylattice(lattice)
	vis_sub_lattice.RemoveAtom(index)
	vis_lattice(vis_sub_lattice, "sub_lattice0.xyz")


	# delete the atoms from the defect cluster
	sub_lattice.index = numpy.arange(lattice.NAtoms)
	sub_lattice.RemoveAtom(remove_id)
	print "Number of defect atoms extracted:", len(remove_id)
	vis_lattice(sub_lattice, "sub_lattice0.xyz")
	sys.exit()
	Sub_Reference_Force = copy.copy(sub_lattice.Force)
	
	# calculate the forces on the system sans the defect
	sub_lattice.eval()
	#sub_lattice.Force -= Sub_Reference_Force
	Force_Vec = (sub_lattice.Force[0::3])**2 + (sub_lattice.Force[1::3])**2 + (sub_lattice.Force[2::3])**2
	
	# now we subtract the shell that interacted with the defect atoms
	interaction_count = 0
	for i in xrange(0, sub_lattice.NAtoms):
		if(Force_Vec[i] > force_thresh):
			index = sub_lattice.deletion_key[i]
			remove_id.append(index)
			interaction_count += 1

	print remove_id
	local_index = numpy.array(remove_id)
	index = numpy.arange(lattice.NAtoms)
	index = numpy.delete(index, local_index)
	vis_sub_lattice = copylattice(lattice)
	vis_sub_lattice.RemoveAtom(index)
	vis_lattice(vis_sub_lattice, "sub_lattice1.xyz")

	
	print "Number of atoms directly interacting with defect:", interaction_count
	
	sys.exit()
	
	sub_lattice = copylattice(lattice)
	sub_lattice.DefectMask[:] = True

	id_mapping = numpy.arange(lattice.NAtoms)
	
	# delete the defect and its interaction shell
	
	sub_lattice.id_mapping = numpy.delete(numpy.arange(sub_lattice.NAtoms),remove_id,None)
	sub_lattice.RemoveAtom(remove_id)
	
	remove_id_n = numpy.array(remove_id)
	Nindex = numpy.concatenate((3*remove_id_n, 3*remove_id_n+1, 3*remove_id_n+2), None)
	Sub_Reference_Force = numpy.delete(Reference_Force, Nindex)
	
	#print "Interaction extraction:", sub_lattice.NAtoms
	
	# calculate the interacting radius on the original defect and its surrounding shell (this is used for the fixed outer region)
	
	sub_lattice.eval()
	sub_lattice.Force -= Sub_Reference_Force
	Force_Vec = (sub_lattice.Force[0::3])**2 + (sub_lattice.Force[1::3])**2 + (sub_lattice.Force[2::3])**2
	
	shell_count = 0
	for i in xrange(0, sub_lattice.NAtoms):
		if(Force_Vec[i] > force_thresh):
			remove_id.append(id_mapping[i])
			shell_count += 1
			
	print "Number of atoms indirectly interacting with defect:", shell_count
	# so now we have a list "remove id" that corresponds to the atoms that we need in our simulation.
	

	print "Total number of retained atoms:", len(remove_id)
	# Now we use what we know to translate
	local_index = numpy.array(remove_id)
	index = numpy.arange(lattice.NAtoms)
	index = numpy.delete(index, local_index)
	vis_sub_lattice = copylattice(lattice)
	vis_sub_lattice.RemoveAtom(index)
	vis_lattice(vis_sub_lattice, "sub_lattice3.xyz")
	
	print "Number of atoms after extraction:", sub_lattice.NAtoms

	
	sys.exit()
	return sub_lattice

def main():
	if(len(sys.argv)>3 or len(sys.argv)<2):
		print "defect_module.py input_lattice reference_latice - Identify defects through XOR" 
		print "defect_module.py input_lattice - Identify defects through energy binning"
		
	if(len(sys.argv)==3):
		reference_lattice = read_lattice(sys.argv[2])
		input_lattice = read_lattice(sys.argv[1])
		lattice = identify_defects_classical(copylattice(input_lattice), reference_lattice, 1, 3)
		for i in xrange(0, lattice.NDefects):
			lattice.Specie[lattice.defect_list[i]] = lattice.group_id[i] + 1
		vis_lattice(lattice, "defect_lattice.xyz")
	else:
		input_lattice = read_lattice(sys.argv[1])
		lattice = get_defect_lattice_from_energy(input_lattice)
		
	defect_lattice = copylattice(input_lattice)
	for i in xrange(0, lattice.Interstitials):
		defect_lattice.AddAtom(Pos=lattice.defect_points[i], Specie="i")
	for i in xrange(lattice.Interstitials, lattice.Interstitials+lattice.Vacancies):
		defect_lattice.AddAtom(Pos=lattice.defect_points[i], Specie="v")
		
	vis_lattice(defect_lattice, "literal_defect_lattice.xyz")
	
	
	
	
	#for i in xrange(0, len(lattice.chull_list)):
	#	lattice.Specie[lattice.chull_list[i]] = "2"

	
	
	

if __name__ == '__main__':
	lattice = main()
