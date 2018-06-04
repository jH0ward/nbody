# frags.py
from table import *
import numpy as np
from numpy import linalg as LA
from collections import deque
import logging

def find_mons(input):
    fraglogger=logging.getLogger('frags')

### Slightly modified vdW radii, especially for HF clusters
### alkali metal values for cations = 
    cov = {'H' : 1.266, 'Li': 1.8, 'C' : 1.829, 'N' : 1.757, 'O' : 1.682, 
           'F' : 1.287, 'Na': 1.5, 'P' : 1.9,   'S' : 1.8,   'Cl': 1.7,
           'Br': 1.8  , 'K' : 1.8, 'Rb': 2.0,   'I' : 1.9,   'Cs': 2.1 }

    f = open(input, "r")
    natom = int(f.readline().strip().split()[0])
    xyz_in = np.zeros(natom*3)
    xyz_in.shape=(natom,3)
    atnum=[]


### Read in coordinates
    for i in range (natom):
       line = f.readline().strip()
       atnum.append(int(line.split()[0]))
       for j in range (3):
          xyz_in[i][j]= line.split()[j+1]
### Add symbol logic
    atsym=[]
    for i in atnum:
       atsym.append(get_symbol(i))

### Make a tuple out of atsym and xyz
    coords = (atsym, xyz_in)

### Make dist. matrix
    dist = np.zeros(natom*natom)
    dist.shape=(natom,natom)
    for i in range (natom):
        for k in range(natom):
            if i > k:
                for j in range (3):
                    dist[i][k] = LA.norm(xyz_in[i]-xyz_in[k])
                    dist[k][i] = dist[i][k]
### Set up queue
    queue = deque(range(natom))
## Lol (list of lists) holds frags
    lol=[]
#### Put first atom from queue into lol
    while queue:
        lol.append([])
        lol[-1].append(queue.popleft())
### Loop over frags
        for i in range (len(lol)):
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
                            n_no=n_no+1
                        k=k+1 
                j=j+1
    
    mons=lol
    fraglogger.info('\nIdentified these groups as monomers\n%s\n'%mons)
    frag_back = (coords, mons)
    return frag_back
