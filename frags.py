# frags.py
from table import *
import numpy as np
from numpy import linalg as LA
from collections import deque
import logging
import itertools
import os


def read_geom_from_file(input):
    f = open(input, "r")
    natom = int(f.readline().strip().split()[0])
    xyz_in = np.zeros(natom * 3)
    xyz_in.shape = (natom, 3)
    atnum = []

    ### Read in coordinates
    for i in range(natom):
        line = f.readline().strip()
        atnum.append(int(line.split()[0]))
        for j in range(3):
            xyz_in[i][j] = line.split()[j + 1]
    ### Add symbol logic
    atsym = []
    for i in atnum:
        atsym.append(get_symbol(i))

    ### Make a tuple out of atsym and xyz
    coords = (atsym, xyz_in)
    return coords


def find_mons(input):
    fraglogger = logging.getLogger("frags")

    ### Slightly modified vdW radii, especially for HF clusters
    ### alkali metal values for cations =
    cov = {
        "H": 1.266,
        "Li": 1.8,
        "C": 1.829,
        "N": 1.757,
        "O": 1.682,
        "F": 1.287,
        "Ne": 1.4,
        "Na": 1.5,
        "P": 1.9,
        "S": 1.8,
        "Cl": 1.7,
        "Ar": 1.8,
        "Br": 1.8,
        "Kr": 1.8,
        "K": 1.8,
        "Rb": 2.0,
        "I": 1.9,
        "Cs": 2.1,
    }

    atsym, xyz_in = read_geom_from_file(input)
    natom = len(atsym)

    ### Make dist. matrix
    dist = np.zeros(natom * natom)
    dist.shape = (natom, natom)
    for i in range(natom):
        for k in range(natom):
            if i > k:
                for j in range(3):
                    dist[i][k] = LA.norm(xyz_in[i] - xyz_in[k])
                    dist[k][i] = dist[i][k]
    ### Set up queue
    queue = deque(list(range(natom)))
    ## Lol (list of lists) holds frags
    lol = []
    #### Put first atom from queue into lol
    while queue:
        lol.append([])
        lol[-1].append(queue.popleft())
        ### Loop over frags
        for i in range(len(lol)):
            ###     Loop over atoms
            j = 0
            while j < len(lol[i]):
                ###     Loop over queue
                n_no = 0
                # Count the number of "no"s
                while n_no < len(queue):
                    k = 0
                    while k < len(queue):
                        an1 = lol[i][j]
                        as1 = atsym[an1]
                        an2 = queue[k]
                        as2 = atsym[an2]
                        thr = cov[as1] + cov[as2]
                        if dist[an1][an2] < thr:
                            ### Add queue[k] to lol[i] (atom to frag)
                            lol[i].append(an2)
                            queue.remove(an2)
                            break
                        else:
                            n_no = n_no + 1
                        k = k + 1
                j = j + 1

    mons = lol

    fraglogger.info("\nIdentified these groups as monomers\n%s\n" % mons)
    coords = (atsym, xyz_in)
    frag_back = (coords, mons)
    return frag_back


def read_mons(input, frag_file):
    ## Manually specified list of fragments if frags.in provided.
    f = open(frag_file, "r")
    lines = f.readlines()
    f.close()
    mons = []
    for line in lines:
        line = line.strip()  # throw whitespace away
        # Ignore comments
        try:
            if line[0] == "#":
                continue
        except IndexError:  # line was empty
            continue

        # Lines look like  0 1 2 for monomer of atoms 0, 1, and 2 in a frag
        atoms_in_fragment = line.split()
        # Make ints not strings
        atoms_in_fragment = [int(a) for a in atoms_in_fragment]
        mons.append(atoms_in_fragment)

    # coords is the tuple of (atomic_symbols, xyz_geometry)
    coords = read_geom_from_file(input)

    # Do a couple of QA checks and raise fatal errors
    flat_atoms_list = list(
        itertools.chain(*mons)
    )  # standard way of flattening a list of lists into 1d array
    natom = len(coords[0])
    n_atoms_mapped = len(flat_atoms_list)

    # Fail if n_atoms_mapped != natom
    assert n_atoms_mapped == natom, (
        "Mapped " + str(n_atoms_mapped) + " atoms vs natom==" + str(natom) + "!"
    )

    # Checks for missing / repeat atomic indices
    missing_atoms = [a for a in range(natom) if a not in flat_atoms_list]
    assert len(missing_atoms) == 0, (
        "Uhoh missing some atoms in the fragment specifications "
        + str(missing_atoms)
        + " . I bet you repeated one of the atom indices"
    )

    return coords, mons
