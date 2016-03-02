import pesto.pio
import pesto.putil
import pesto.pratn
import pesto.pdefect

import sys



def main():
#	print len(sys.argv)	
	if len(sys.argv) < 2:
		print "Arguments: Input_lattice. \nOpt: Reference_lattice"
		quit()
	initial = pesto.pio.read_lattice(sys.argv[1])

	if(len(sys.argv)==2):
		#initial = defect_module.get_defect_lattice_from_energy(initial)
		initial = pesto.pdefect.get_defect_lattice_from_specie(initial, "3", 3)
	else:
		reference = pesto.pio.read_lattice(sys.argv[2])
		initial = pesto.pdefect.identify_defects_classical(initial, reference,1, 4)
		
		
	#initial.EnableVerbose(True)
	initial.eval()
	
	pesto.pio.vis_lattice(pesto.pdefect.get_defect_lattice_from_mask(initial), "defect.xyz")
	
	pesto.pio.vis_lattice(initial, "RAT0.xyz")
	#lattice = ratn_method(initial)
	totalpe = 0
	succeeded = 0
	failed_old_min = 0
	failed_new_min = 0
	failed_didnt_converge = 0
	loops = 100
	Failed = True
	
	initial.local=True
	saddle_lattice = pesto.pratn.ratn_method(pesto.putil.copylattice(initial))
	if(saddle_lattice.Failed):
		print "Transition search failed."
		return


	pesto.pio.write_lattice(saddle_lattice, "saddle.dat")
	pesto.pio.vis_lattice(saddle_lattice, "RAT1.xyz")


	
if __name__ == '__main__':
	main()	