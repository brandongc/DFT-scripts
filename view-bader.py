#!/usr/bin/env python

"""This script is for visualizing the output of a bader charge analysis
of a VASP calculation using the ASE package.

The default version of attach_charges provided in ASE associates the charge 
field of each atom as Atomic number - Bader result. This is undesirable in 
many cases.

There is also a units (Angstrom vs Borh) mismatch between the bader code and 
VASP output which I have added an option to correct.

I'm hoping these changes will be accepted into a new version of the ASE
package.
"""
import numpy as np
from ase.io import read, write
from ase.units import Bohr
from ase.visualize import view

def attach_charges(atoms, fileobj='ACF.dat', displacement=1e-4, use_diff=True,
                   use_bohr=True):
    """Attach the charges from the fileobj to the Atoms."""
    if isinstance(fileobj, str):
        fileobj = open(fileobj)

    sep = '---------------'
    i = 0 # Counter for the lines
    k = 0 # Counter of sep
    assume6columns = False
    for line in fileobj:
        if line[0] == '\n': # check if there is an empty line in the 
            i -= 1          # head of ACF.dat file
        if i == 0:
            headings = line
            if 'BADER' in headings.split():
                j = headings.split().index('BADER')
            elif 'CHARGE' in headings.split():
                j = headings.split().index('CHARGE')
            else:
                print 'Can\'t find keyword "BADER" or "CHARGE".' \
                +' Assuming the ACF.dat file has 6 columns.'
                j = 4
                assume6columns = True
        if sep in line: # Stop at last seperator line
            if k == 1:
                break
            k += 1
        if not i > 1:
            pass
        else:
            words = line.split()
            if assume6columns is True:
                if len(words) != 6:
                    raise IOError('Number of columns in ACF file incorrect!\n'
                                  'Check that Bader program version >= 0.25')
                
            atom = atoms[int(words[0]) - 1]
            if use_diff:
                atom.charge = atom.number - float(words[j])
            else:
                atom.charge = float(words[j])
            if displacement is not None: # check if the atom positions match
                if use_bohr:
                    xyz = np.array([float(w) for w in words[1:4]]) * Bohr
                else:
                    xyz = np.array([float(w) for w in words[1:4]])
                assert np.linalg.norm(atom.position - xyz) < displacement
        i += 1

atoms = read('POSCAR')
attach_charges(atoms, 'ACF.dat', use_bohr=False, use_diff=False)
# example of a cheap trick to get spin-charges to stand out
# chgsplit.pl CHGCAR
# chgsum.pl AECCAR0 AECCAR2
# after e.g. bader CHGCAR_mag -ref CHGCAR_sum
#nums = atoms.get_atomic_numbers()
#for a in atoms:
#    if abs(a.charge) > 0.05:
#        print a.symbol, a.charge, a.position
#        nums[a.index] += 1
#atoms.set_atomic_numbers(nums)
view(atoms)
        
