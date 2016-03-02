import os
MDTYPE = os.environ.get("MDTYPE")
if(MDTYPE != "LAMMPY"):
    print "These examples require that the environmental variable MDTYPE is set to LAMMPY - to use LAMMPS via the python interface."
    print "For the examples the variable will be set for you, but this will not happen under normal execution."
    os.environ["MDTYPE"] = "LAMMPY"
    
    
from pesto.pio import *
from pesto.putil import *
from pesto import pminimise 

lattice = read_lattice("../perfect.lammps")
lattice.eval()

# also eval_lattice(lattice)


# save the energy of the perfect crystal
reference_energy = lattice.TPE

# identify center of lattice
print "Lattice Dimensions (max, min):", lattice.Dim
midpoint = [lattice.Dim[0]/2, lattice.Dim[1]/2, lattice.Dim[2]/2]
print "Lattice MidPoint", midpoint

# get the sorted list of neighbours to the midpoint
neighbors = nearest_neighbours_sorted(lattice, midpoint)

print "The first 5 nearest atoms to the midpoint are [id, distance]:\n", neighbors[0:5]

# create monovacancy
lattice.RemoveAtom(neighbors[0][0])

# evaluate the lattice
lattice.eval()

print "Unrelaxed monovacancy energy:", lattice.TPE
lattice = pminimise.minimise_lattice(lattice, "LBFGS", 1E-4)
work = pio.forcecalls

# insert some noise and try all optimisers

print "Relax using Scipy LBFGS"
lattice.Pos += numpy.random.random(lattice.NAtoms*3) / 100
test_lattice = pminimise.minimise_lattice(copylattice(lattice), "SLBFGS", 1E-5)
print "Force calls:", pio.forcecalls - work
print "Vacancy formation energy (after LBFGS):", lattice.TPE - reference_energy, lattice.MagForce()
work = pio.forcecalls

print "Relax using Locally implemented LBFGS"
test_lattice = pminimise.minimise_lattice(copylattice(lattice), "LBFGS", 1E-8)
print "Force calls:", pio.forcecalls - work
work = pio.forcecalls

print "Vacancy formation energy (after SBFGS):", lattice.TPE - reference_energy, lattice.MagForce()

print "Relax using special SD"
test_lattice = pminimise.minimise_lattice(copylattice(lattice), "SD", 1E-5)
print "Force calls:", pio.forcecalls - work
work = pio.forcecalls

print "Vacancy formation energy (after SD):", lattice.TPE - reference_energy, lattice.MagForce()

print "Relax using Quickmin"
test_lattice = pminimise.minimise_lattice(copylattice(lattice), "Quickmin", 1E-5)
print "Force calls:", pio.forcecalls - work
work = pio.forcecalls

print "Vacancy formation energy (after Quickmin):", lattice.TPE - reference_energy, lattice.MagForce()

print "Relax using Conjugate Gradient"
test_lattice = pminimise.minimise_lattice(copylattice(lattice), "CG", 1E-5)
print "Force calls:", pio.forcecalls - work
work = pio.forcecalls

print "Vacancy formation energy (after Conjugate Gradient):", lattice.TPE - reference_energy, lattice.MagForce()


print "Saving the monovacancy lattice"

write_lattice(test_lattice, "../mono_vacancy.lammps")
