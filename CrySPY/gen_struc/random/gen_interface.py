from ase.atoms import Atoms
from ase.io import (read, write)
from pyxtal.lattice import Lattice
import numpy as np
import os
from ase.constraints import FixAtoms
from pyxtal import pyxtal
from random import randint 
from ase.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes
from ase.optimize import BFGS

def update_pos(up_struc, sub_struc, buffer):
    '''Automatically match substrate and superstructure'''
    up_pos = up_struc.positions
    sub_cellpar = sub_struc.cell.cellpar()
    displace = sub_cellpar[2] + buffer

    # 直接将 z 坐标加上位移量
    up_pos[:, 2] += displace

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
    sub_symbol = sub_struc.get_chemical_symbols()
    sub_cellpar = sub_struc.cell.cellpar()

    up_symbol = up_struc.get_chemical_symbols()
    new_up_pos = update_pos(up_struc, sub_struc, buffer)

    inter_cellpar = create_inter_lattice(up_struc, sub_struc, buffer, vacum)

    # 创建一个新的 Atoms 对象，用于存储界面结构
    inter_struc = Atoms(symbols=sub_symbol, positions=sub_pos, cell=inter_cellpar, pbc=True)

    # 使用 ase.Atoms.extend 方法将更新后的顶层结构添加到基底结构中
    inter_struc.extend(Atoms(symbols=up_symbol, positions=new_up_pos, cell=inter_cellpar, pbc=True))

    # 设置约束条件
    fix = FixAtoms(mask=inter_struc.positions[:, 2] < sub_cellpar[2])
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



class LinearInverseForce(Calculator):

    implemented_properties = ['energy', 'energies', 'forces', 'free_energy']
    implemented_properties += ['stress', 'stresses']  # bulk properties
    default_parameters = {
        'k': 1.0,
        'rc': None,
    }
    nolabel = True

    def __init__(self, **kwargs):

        Calculator.__init__(self, **kwargs)

        if self.parameters.rc is None:
            self.parameters.rc = 2

        self.nl = None

    def calculate(
        self,
        atoms=None,
        properties=None,
        system_changes=all_changes,
    ):
        if properties is None:
            properties = self.implemented_properties

        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)

        k = self.parameters.k
        rc = self.parameters.rc

        if self.nl is None or 'numbers' in system_changes:
            self.nl = NeighborList(
                [rc / 2] * natoms, self_interaction=False, bothways=True
            )

        self.nl.update(self.atoms)

        positions = self.atoms.positions
        cell = self.atoms.cell

        energies = np.zeros(natoms)
        forces = np.zeros((natoms, 3))

        for ii in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(ii)
            cells = np.dot(offsets, cell)

            # pointing *towards* neighbours
            distance_vectors = positions[neighbors] + cells - positions[ii]

            r = np.sqrt((distance_vectors ** 2).sum(1))
            within_cutoff = r < rc

            pairwise_forces = -k * (1.0 / r) * within_cutoff

            pairwise_forces = pairwise_forces[:, np.newaxis] * distance_vectors

            forces[ii] += pairwise_forces.sum(axis=0)



        energy = energies.sum()
        self.results['energy'] = energy
        self.results['energies'] = energies

        self.results['free_energy'] = energy

        self.results['forces'] = forces

def pre_optimize(atoms,rcut=1.0):
    calc = LinearInverseForce(rc=rcut)
    atoms = Atoms(atoms, calculator=calc)
    opt = BFGS(atoms)
    opt.run(fmax=0.05, steps=1000)
    
    return atoms



if __name__ == '__main__':
    generate_interface(sub_struc_path, spice, spice_num, up_thickness, buffer=0, vacum=0, sub_struc_format='vasp')
