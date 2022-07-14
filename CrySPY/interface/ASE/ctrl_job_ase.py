'''
Control jobs in ase
'''

import os
import shutil

from . import structure as ase_structure
from ...IO import read_input as rin


def next_stage_ase(stage, work_path):
    # ---------- skip_flag
    skip_flag = False

    # ---------- rename ase files at the current stage
    ase_files = [rin.ase_infile, rin.ase_outfile,
                 'POSCAR', 'CONTCAR']
    for f in ase_files:
        if not os.path.isfile(work_path+f):
            raise IOError('Not found '+work_path+f)
        os.rename(work_path+f, work_path+'stage{}_'.format(stage)+f)

   # ---------- cp CONTCAR POSCAR
    shutil.copyfile(work_path+'stage{}_CONTCAR'.format(stage),
                    work_path+'POSCAR')

    # ---------- copy the input file from ./calc_in for the next stage
    finfile = './calc_in/'+'{}_'.format(stage + 1)+rin.ase_infile
    shutil.copyfile(finfile, work_path+rin.ase_infile)

    # ---------- return
    return skip_flag


def next_struc_ase(structure, current_id, work_path):
    # ---------- copy files
    if rin.ase_potential is None:
        calc_inputs = [rin.ase_infile]
    else:
        calc_inputs = [rin.ase_infile] + rin.ase_potential
    for f in calc_inputs:
        ff = '1_'+f if f == rin.ase_infile else f
        if not os.path.isfile('./calc_in/'+ff):
            raise IOError('Could not find ./calc_in/'+ff)
        shutil.copyfile('./calc_in/'+ff, work_path+f)

    # ---------- generate POSCAR
    structure.to(fmt='poscar', filename=work_path+'POSCAR')
    if not os.path.isfile(work_path+'POSCAR'):
        raise IOError('Could not find {}POSCAR'.format(work_path))

    # ---------- Change the title of POSCAR
    with open(work_path+'POSCAR', 'r') as f:
        lines = f.readlines()
    lines[0] = 'ID_{}\n'.format(current_id)
    with open(work_path+'POSCAR', 'w') as f:
        for line in lines:
            f.write(line)