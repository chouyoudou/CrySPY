from ase.atoms import Atoms
from ase.io import (read, write)
from pyxtal.lattice import Lattice
import numpy as np
import os
from ase.constraints import FixAtoms
from pyxtal import pyxtal
from random import randint 
from ase.calculators.lj import LennardJones



def update_pos(up_struc, sub_struc, buffer):
    '''Automatically match substrate and superstructure'''
    up_pos = up_struc.positions
    sub_cellpar = sub_struc.cell.cellpar()
    displace = sub_cellpar[2]   ### c
    for i in range(len(up_pos)):
        up_pos[i][2] = up_pos[i][2] + displace + buffer
    return up_pos

def translate_pos(structure, direction, displace):
    '''Automatically match substrate and superstructure
    direction 1, 2 ,3 = x, y, z'''

    new_structure = structure
    #new_structure.positions[:][direction-1] = new_structure.positions[:][direction-1] + displace

    for i in range(len(new_structure.positions)):
        new_structure.positions[i][direction-1] = new_structure.positions[i][direction-1] + displace
    return new_structure

def create_inter_lattice(up_struc, sub_struc, buffer, vacum):
    sub_cellpar = sub_struc.cell.cellpar()
    up_cellpar = up_struc.cell.cellpar()
    
    ### update cell
    inter_cellpar = sub_cellpar
    inter_cellpar[2] = inter_cellpar[2] + up_cellpar[2] + vacum + buffer
    return inter_cellpar

def match_cell(up_struc, sub_struc, buffer, vacum):
    sub_pos = sub_struc.positions
    sub_symbol = sub_struc.get_chemical_formula()
    sub_cellpar = sub_struc.cell.cellpar()
    
    up_symbol = up_struc.get_chemical_formula()
    new_up_pos = update_pos(up_struc, sub_struc, buffer)
    
    inter_pos= np.vstack((new_up_pos , sub_pos))
    inter_cellpar = create_inter_lattice(up_struc, sub_struc, buffer, vacum)
    inter_struc = Atoms(up_symbol+sub_symbol, inter_pos, cell=inter_cellpar)
    #new_up_struc = Atoms(up_symbol, new_up_pos, cell=inter_cellpar)
    fix = FixAtoms( mask = inter_struc.positions[:,2] < sub_cellpar[2] ) ##Fix all atoms whose z-coordinate is less than the lattice constant c of the base, that is, fix all atoms of the base
    inter_struc.set_constraint(fix)
    
    return inter_struc

def separate(interface_struc, thickness):
    sub_thickness = thickness
    interface_parcell = interface_struc.cell.cellpar()
    up_thickness = interface_parcell[2] - thickness


def generate_up_structure(spice, spice_num, cell_par, up_thickness, sub_file='SUB_POSCAR'):
    '''To generate a structure with fixed lattice parameters, this function needs to obtain the element species and lattice parameters. 
    After entering 6 lattice parameters, the Lattice class of pyxtal is generated, and then Lattice is passed to generate_fix_lattice and from_random to generate the structure'''

    lattice = Lattice.from_para(cell_par[0], cell_par[1], up_thickness, cell_par[3], cell_par[4], cell_par[5])
    up_struc = generate_fix_lattice(spice, spice_num, lattice)
    #up_struc.to_file('UP_POSCAR', fmt='poscar')

    return up_struc

def generate_fix_lattice(spice, spice_num, lat_par):
    '''Call pyxtal's functions to generate structures'''
    structure=pyxtal()
    while True:
        sg = randint(2, 230)
        try:
            structure.from_random(3, sg, spice, spice_num, lattice=lat_par)
            structure = structure.to_ase()
            cell_par = structure.cell.cellpar()
            if True not in (structure.positions[:,2] < 0):
                if True not in (structure.positions[:,2] > cell_par[2]):
                    print('finish')
                    break  
                else:
                    print('Some atoms > c, try again')
                    structure=pyxtal()
            else:
                print('Some atoms < 0, try again')
                structure=pyxtal()

        except: 
            print('can not generate up structure, try again')
            pass
    return structure

def generate_interface(sub_struc_path, spice, spice_num, up_thickness, buffer=0, vacum=0, sub_struc_format='vasp'):
    sub_struc = read(sub_struc_path,format=sub_struc_format)
    cell_par = sub_struc.cell.cellpar()
    up_struc = generate_up_structure(spice, spice_num, cell_par, up_thickness)
    inter = match_cell(up_struc, sub_struc, buffer, vacum)
    return inter


if __name__ == '__main__':
    generate_interface(sub_struc_path, spice, spice_num, up_thickness, buffer=0, vacum=0, sub_struc_format='vasp')
