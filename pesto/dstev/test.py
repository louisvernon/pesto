
import numpy
from dstev import dstev
from scipy import linalg

info = 0
alpha = numpy.array([1,2,3,4,5,6], dtype=float)
sigma = numpy.array([5,4,3,2,1], dtype=float)

workspace = numpy.zeros(len(alpha)*2, dtype=float)
z = numpy.zeros((len(alpha), len(alpha)),dtype=float)
#print workspace
print alpha
print sigma
pyz = dstev(jobz="V", n=len(alpha), d=alpha, e=sigma, ldz=len(alpha), z=z, work=workspace, info=info)
#pyz = pyz.transpose()

print "Completion Status:", info

print "DSTEV Eigenvalues:" 
print alpha
print "DSTEV Eigenvectors:"
print pyz

test_matrix = numpy.array(([1,5,0,0,0,0],[5,2,4,0,0,0],[0,4,3,3,0,0],[0,0,3,4,2,0],[0,0,0,2,5,1],[0,0,0,0,1,6]), dtype=float)
print test_matrix
la, v = linalg.eig(test_matrix)
print "scipy.linalg Eigenvalues"
print la
print "scipy.linalg Eigenvectors"
print v
