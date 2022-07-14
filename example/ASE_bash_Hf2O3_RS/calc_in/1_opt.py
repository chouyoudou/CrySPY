from ase import Atoms
from ase.optimize import FIRE
from ase.optimize import BFGS
import numpy as np
from deepmd.calculator import DP
from ase.io import (read, write)
calculator=DP(model="frozen_model.pb")

struc = read('POSCAR',format="vasp")
struc.calc = calculator
dyn = BFGS(struc)
dyn.run(fmax=0.5)
write('CONTCAR', struc, format='vasp')

