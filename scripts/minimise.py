import pesto.pio
import pesto.pminimise

import sys

def main():
#	print os.getcwd()
	if len(sys.argv) != 4:
		print "lattice, minimisation type, tolerance. E.g.\nminimise.py lattice.dat SD 200, where SD is an accelerated steepest descent and 200 is the maximum number of steps"
		print "minimise.py lattice.dat CG 0.0001, where CG is conjugate gradient and 0.0001 is the minimum energy change per line minimisation."
		print "minimise.py latice.dat TA, where TA implies thermal annealing"
		quit()
	
        filename = str(sys.argv[1])
	type=str(sys.argv[2])
	tol=float(sys.argv[3])

	lattice = pesto.pio.read_lattice(filename)
	lattice = lattice.eval()

	lattice = pesto.pminimise.minimise_lattice(lattice,type,tol)
	print "Force Calls", pesto.pio.forcecalls

	pesto.pio.write_lattice(lattice, "relaxedlattice.dat")
					
if __name__ == '__main__':
	main()					
		