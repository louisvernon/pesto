import sys, math
import os

# test whether lammps is correctly built and installed
import lammps
lmp = lammps.lammps()

# check whether they have set the environmental variable.

MDTYPE = os.environ.get("MDTYPE")
if(MDTYPE != "LAMMPY"):
    print "These examples require that the environmental variable MDTYPE is set to LAMMPY - to use LAMMPS via the python interface."
    print "For the examples the variable will be set for you, but this will not happen under normal execution."
    os.environ["MDTYPE"] = "LAMMPY"
    


from pesto.pio import *

lattice = read_lattice("../perfect.lammps")
lattice.EnableVerbose(True)
lattice.eval()


print "Lattice Number of Atoms:", lattice.NAtoms,  (lattice.NAtoms == 5324)
test = (lattice.NAtoms == 5324)
print "Lattice Energy:", lattice.TPE, math.fabs(-18848.1222851-lattice.TPE) < 0.1
test += math.fabs(-18848.1222851-lattice.TPE) < 0.1
test += lattice.MaxForce() < 0.1
print "Lattice Force Maximum:", lattice.MaxForce(), lattice.MaxForce() < 0.1
test += lattice.MagForce() < 0.1
print "Lattice Force Magnitude:", lattice.MagForce(), lattice.MagForce()< 0.1

if (test == 4):
    print "All tests completed without error, everything seems to be working correctly"
else:
    print "At least one test came back with an unexpected result. This suggests something is wrong with your configuration."


