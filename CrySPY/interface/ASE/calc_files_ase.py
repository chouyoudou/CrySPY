'''
Calculation files in ase
'''

import os

from ...IO import read_input as rin


def check_input_ase():
    # ---------- prepare rin.jobfile, rin.ase_potential, rin.rammps_infile
    if rin.ase_potential is None:
        calc_inputs = [rin.jobfile, rin.ase_infile]
    else:
        calc_inputs = [rin.jobfile, rin.ase_infile] + rin.ase_potential

    # ----- check required files
    for f in calc_inputs:
        if f == rin.ase_infile:
            finfiles = ['{}_'.format(i) + rin.ase_infile  for i in range(
                1, rin.nstage+1)]
            for ff in finfiles:
                if not os.path.isfile('./calc_in/'+ff):
                    raise IOError('Could not find ./calc_in/'+ff)
        else:
            if not os.path.isfile('./calc_in/'+f):
                raise IOError('Could not find ./calc_in/'+f)
