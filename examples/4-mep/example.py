import os
MDTYPE = os.environ.get("MDTYPE")
if(MDTYPE != "LAMMPY"):
    print "These examples require that the environmental variable MDTYPE is set to LAMMPY - to use LAMMPS via the python interface."
    print "For the examples the variable will be set for you, but this will not happen under normal execution."
    os.environ["MDTYPE"] = "LAMMPY"

from pesto.pio import *
from pesto.putil import *
from pesto import pneb

initial = read_lattice("../mono_vacancy.lammps")
final = read_lattice("../mono_vacancy_hop.lammps")
saddle = read_lattice("../mono_vacancy_saddle.lammps")
print "Converging NEB without interpolated saddle"
pneb.run_neb(initial, final, images=8, springk=1.0, minimiser="LBFGS")

print "Converging NEB via interpolated saddle:"
pneb.run_neb(initial, final, images=8, springk=1.0, minimiser="LBFGS", saddle=saddle)
