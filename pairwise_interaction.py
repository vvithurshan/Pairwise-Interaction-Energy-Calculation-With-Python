import MDAnalysis as mda
import itertools
from numba import jit
import warnings
import time
import numpy as np
import pickle

warnings.filterwarnings("ignore")
warnings.filterwarnings("default")

psf = "ionized_dry.psf"
# dcd = "only_10_frames.dcd"
dcd = "protein.pdb"
# psf = "PSF_12A21_231965_C1.psf"
# dcd = "PDB_12A21_231965_C1.pdb"
# Load the trajectory and topology files
u = mda.Universe(psf, dcd)

sel1 = u.select_atoms('protein')
sel2 = u.select_atoms('protein')

start = time.time()

def distance(atm_1, atm_2):
    dist = np.linalg.norm(np.array(atm_1) - np.array(atm_2))
    return dist
    
def generate_pairs(sel1, sel2):
    target = [f"{res.segid}_{res.resname}_{res.resid}" for res in sel1.residues]
    source = [f"{res.segid}_{res.resname}_{res.resid}" for res in sel2.residues]
    pairProduct = [tuple(sorted(pairs)) for pairs in itertools.product(source, target) if pairs[0] != pairs[1]]
    pairProduct = list(set(pairProduct))
    return pairProduct

pairProduct = generate_pairs(sel1, sel2)
# pairProduct = [list(x) for x in set(tuple(x) for x in pairProduct)]
end = time.time()
print(end - start)
print(len(pairProduct))
print(pairProduct)



from parmed.charmm import CharmmPsfFile, CharmmCrdFile, CharmmParameterSet

psf = "ionized_dry.psf"
dcd = "only_10_frames.dcd"
# dcd = "protein.pdb"
u = mda.Universe(psf, dcd)
params = CharmmParameterSet('par_all36_prot.prm')
switchdist = 10  # switching distance
cutoff = 12  # cutoff distance

## Dict for Sigma
Dict_Sigma = {}

## Energy functions
def LJ_Energy(epsilon, sigma, r):
    r_6 = r**6
    r_12 = r_6**2
    # Calculate switching function
    if r > switchdist:
        switching = 1 - np.exp(-(r-switchdist)**2/(2*(cutoff-switchdist)**2))
    else:
        switching = 1
        
    if Dict_Sigma.get(sigma):
        sigma_6, sigma_12 = Dict_Sigma[sigma]
    else:
        sigma_6 = sigma**6
        sigma_12 = sigma_6**2
        Dict_Sigma[sigma] = (sigma_6, sigma_12)

    LJ_energy = 4*epsilon*((1/r*(sigma_12/r_12)) - 2.0/r*(sigma_6/r_6)) ## corr 0.99
    LJ_potentials = LJ_energy * switching 
    return LJ_potentials

def Ele_Energy(charge_a1, charge_a2, r):
    return charge_a1 * charge_a2/r * 1390 * 0.239



## Calling pairs
def Energy_Calculation(pairProduct):
    Dict_Energy = {}
    for I, pairs in enumerate(pairProduct):
        print(f'pair {I} ',end = '')
        ## continue if pair 1 == pair 2
        if pairs[0] == pairs[1]:
            continue

        key1 = f"{pairs[0].split('_')[0]}{pairs[0].split('_')[1]}{pairs[0].split('_')[2]}"
        key2 = f"{pairs[1].split('_')[0]}{pairs[1].split('_')[1]}{pairs[1].split('_')[2]}"

        if Dict_Energy.get(f"{key2}-{key1}"):
            continue

        ### extracting pair 1 and pair 2 
        pair1 = f"segid {pairs[0].split('_')[0]} and resid {pairs[0].split('_')[2]}"
        pair2 = f"segid {pairs[1].split('_')[0]} and resid {pairs[1].split('_')[2]}"

        ### CA atoms
        pair1_CA = u.select_atoms(f"segid {pairs[0].split('_')[0]} and resid {pairs[0].split('_')[2]} and name CA")
        pair2_CA = u.select_atoms(f"segid {pairs[1].split('_')[0]} and resid {pairs[1].split('_')[2]} and name CA")

        if distance(pair1_CA.positions, pair2_CA.positions) > 12:
            continue
            
        ###
        atom_1 = u.select_atoms(pair1)
        atom_2 = u.select_atoms(pair2)

        ###
        n_frames = u.trajectory.n_frames

        lst_LJ = np.zeros(n_frames)
        lst_Ele = np.zeros(n_frames)
        lst_Total = np.zeros(n_frames)

        for frame, ts in enumerate(u.trajectory):
            # print('tarj')
            LJ_potentials = 0
            Electrostatics = 0

            for a1 in atom_1:
                name_a1 = a1.type
                charge_a1 = a1.charge
                epsilon_a1 = params.atom_types[name_a1].epsilon
                sigma_a1 = params.atom_types[name_a1].rmin

                for a2 in atom_2:
                    name_a2 = a2.type
                    charge_a2 = a2.charge
                    epsilon_a2 = params.atom_types[name_a2].epsilon
                    sigma_a2 = params.atom_types[name_a2].rmin 
                    r = distance(tuple(a1.position), tuple(a2.position))
                
                    # Calculate Lennard-Jones energy with switching function
                    epsilon = np.sqrt(epsilon_a1 * epsilon_a2)
                    sigma = (sigma_a1 + sigma_a2)

                    LJ_potentials += LJ_Energy(epsilon, sigma, r)

                    Electrostatics += Ele_Energy(charge_a1, charge_a2, r)
            
            ### loading energies into array
            lst_LJ[frame] = LJ_potentials
            lst_Ele[frame] = Electrostatics
            lst_Total[frame] = LJ_potentials + Electrostatics


        Dict_Energy[f"{key1}-{key2}"] = {'LJ':lst_LJ, "Ele": lst_Ele, "Total": lst_Total}

    return Dict_Energy


import multiprocessing
import numpy as np

# Create a Pool object.
pool = multiprocessing.Pool()

# Split the data into chunks.
pairChunks = np.array_split(pairProduct, 2)

# Calculate the energy for each chunk in parallel.
results = pool.map(Energy_Calculation, pairChunks)

Dict_Energy = {}
# Print the results.
for result in results:
    Dict_Energy.update(result)
    
with open('Dict_Energy.pickle', 'wb') as file:
    pickle.dump(Dict_Energy, file) 
   
   

# Close the pool.
pool.close()
pool.join()

