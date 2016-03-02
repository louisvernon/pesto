import os
MDTYPE = os.environ.get("MDTYPE")
if(MDTYPE != "LAMMPY"):
    print "These examples require that the environmental variable MDTYPE is set to LAMMPY - to use LAMMPS via the python interface."
    print "For the examples the variable will be set for you, but this will not happen under normal execution."
    os.environ["MDTYPE"] = "LAMMPY"

from pesto.pio import *
from pesto.putil import *
from pesto import pminimise, pdefect, pratn

# read in mono_vacancy generated in 2-example
defect_lattice = read_lattice("../mono_vacancy.lammps")

lattice = read_lattice("../perfect.lammps")

# perform defect analysis against reference lattice

defect_lattice = pdefect.identify_defects_classical(defect_lattice, lattice)

# do a simple transition search, trying a variety of accuracies

calls = pio.forcecalls

hop_lattice = pratn.ratn_method(copylattice(defect_lattice), min_tol=0.1)

pe_low_tol = hop_lattice.Barrier
calls_low_tol = pio.forcecalls - calls

calls = pio.forcecalls

hop_lattice = pratn.ratn_method(copylattice(defect_lattice), min_tol=0.005)

pe_medium_tol = hop_lattice.Barrier
calls_medium_tol = pio.forcecalls - calls

calls = pio.forcecalls

hop_lattice = pratn.ratn_method(copylattice(defect_lattice), min_tol=0.0001)

pe_high_tol = hop_lattice.Barrier
calls_high_tol = pio.forcecalls - calls

print "Energy Barriers, force calls (low, medium, high accuracies):", pe_low_tol, calls_low_tol, pe_medium_tol, calls_medium_tol, pe_high_tol, calls_high_tol

print "Perhaps contrary to intuition the RAT method converges down the PES towards the saddle. Lower accuracy saddles tend to overestimate the transition barrier."
print "The calculated barriers and the required forcecalls will always change due the random initial displacement of the search."

# relax the saddle into the new minimum by projecting along the minimum mode

hop_lattice.Pos += hop_lattice.min_mode*0.2

end_lattice = pminimise.minimise_lattice(hop_lattice, "LBFGS", 1E-4)

print "Separation of relaxed end point and origin is:", separation(defect_lattice, end_lattice)

write_lattice(end_lattice, "../mono_vacancy_hop.lammps")

saddle_lattice = copylattice(hop_lattice)
saddle_lattice.Pos = saddle_lattice.Saddle
write_lattice(saddle_lattice, "../mono_vacancy_saddle.lammps")

