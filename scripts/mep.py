
from pesto.pio import *
from pesto.putil import *
from pesto import pneb
import sys
try:
    import matplotlib.pyplot as plt
    
except:
    pass


initial = read_lattice(sys.argv[1])
final = read_lattice(sys.argv[2])
images=int(sys.argv[3])
if(len(sys.argv)==5):    
    saddle = read_lattice(sys.argv[4])
else:
    saddle = False

band = pneb.run_neb(initial, final, images=int(sys.argv[3]), springk=1.0, minimiser="SBFGS", saddle=saddle)


initial.Pos=band.Pos[0:3*band.NAtoms]
vis_lattice(initial, "Band0.xyz")

energy = []
distance = []
energy.append(0)
distance.append(0)



for i in xrange(1,band.nimages):
    initial.Pos=band.Pos[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms]
    energy.append(band.Energy[i] - band.Energy[0])
    distance.append(magnitude(separation_vector_raw(band.Pos[(i)*3*band.NAtoms:(i+1)*3*band.NAtoms], band.Pos[0:3*band.NAtoms], band.PBC, band.Dim)))
    vis_lattice(initial, "Band"+str(i)+".xyz")
    

csv = open("Band.csv","w")
csv.write("Displacement (A), Energy (eV)\n")
for i in xrange(len(energy)):
    csv.write(str(distance[i]) + "," + str(energy[i]) + "\n")
csv.flush()
csv.close()
    
try:
    plt.figure()
    plt.plot(distance, energy)
    plt.xlabel("Displacement Along Band (A)")
    plt.ylabel("Energy (eV)")
    plt.savefig("Band.pdf")
except:
    print "Install matplotlib for figure generation."