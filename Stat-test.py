from re import L
from matplotlib.pyplot import ylabel
import numpy as np
from numpy.linalg import eig , det
import math

covmat=np.matrix([[100,55],[55,36]])
n= 50

# correlation matrix
d = np.sqrt(covmat.diagonal())
corr_mat = ((covmat.T/d).T)/d

print(corr_mat)

L,v=eig(covmat)


print(math.log(8,np.exp))

