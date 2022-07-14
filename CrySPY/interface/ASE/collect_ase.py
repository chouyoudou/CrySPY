'''
Collect results in ase
'''

import numpy as np

from . import structure as ase_structure
from ...IO import read_input as rin
from pymatgen.core import Structure

def collect_ase(current_id, work_path):
    # ---------- check optimization in current stage & obtain energy
    energy = np.nan
    check_opt = ' Done'
    try:
        with open(work_path+rin.ase_outfile, 'r') as fout:
            lines = fout.readlines()
            energy = float(lines[-1].split()[3])  # in eV (units is metal)
            energy = energy/float(rin.natot)    # eV/cell --> eV/atom
            check_opt = 'done'
    except Exception as e:
        print(e)
        print('    Structure ID {0}, could not obtain energy from {1}'.format(
            current_id, rin.ase_outfile))
        energy = np.nan    # error
        check_opt = 'no_file'

    # ---------- obtain magmom
    magmom = np.nan    # magnetic moment is not calculated

    # ---------- collect the last structure
    opt_struc = get_opt_struc_vasp(work_path+'CONTCAR')

    # ---------- return
    return opt_struc, energy, magmom, check_opt

def get_opt_struc_vasp(file_name):
    try:
        opt_struc = Structure.from_file(file_name)
    except Exception:
        opt_struc = None
    return opt_struc