import numpy as np
import matplotlib.pyplot as plt

from sympy import Matrix, init_printing
init_printing()

def printM(M):
    return display(Matrix(M))

from IPython.display import Image

########################################

NSITES = 3
N_electrons = 2
t = 1.0
U = 4.0