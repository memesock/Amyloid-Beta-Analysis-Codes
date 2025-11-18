#FIRST IM INCLUDING CONTACT BASED CLUSTERING.
#THE WAY IT WORKS IS FOLLOWING:

# WE FIRST NOTE DOWN SIDE CHAIN ATOMS FROM ALL RESIDUES. IT MUST BE NOTED SOME OF THE RESIDUES DONE HAVE ANY SIDE CHAINS AT ALL
# FOR THE RESIDUES HAVING NO SIDE CHAINS, I BELEIVE A GOOD WAY OF DEPICTION IS JUST ASSIGNING IT INFINITE DISTANCE.
# FOR ANY OTHER RESIDUES, MINIMUM OF ANY OF THE SIDE CHAIN IS TAKEN AS THE DISTANCE BETWEEN TWO RESIDUES.



import os
import math
import numpy as np
import mdtraj as md
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
print(f"Current working directory: {os.getcwd()}")
plt.rcParams["font.family"] = "DejaVu Serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"

# ------------------------ User parameters ------------------------
sc_names = {'SC', 'SC1', 'SC2', 'SC3', 'SC4'}

# ONE CAN CHAGE THIS TO THEIR OWN SYSTEM AND DIMER SYSTEM USING (0,40) AND (40,80) FOR DIFFERENT CHAINS

# MAKE SURE TO CHANGE +4 TO +1 AT LINE 83
chainA_range = range(0, 40)   
chainB_range = range(0, 40)



#######################################################################################################
########################################## PARAMETERS #################################################
# THIS IS THE INFINITE DISTANCE IM TAKING FOR ANY RESIDUES NOT HAVING THE SIDE CHAINS
large_r = 100.0
# BETA FOR ASSIGNING THE PROBABILITES
beta = 500.0   # nm^-1
# r0 USED FOR DISTANCE BASED PROBABILITY CALCULATION
r0 = 0.8       # nm
# NUMBER OF CLUSTERS I WANTED FOR MY SYSTEM
n_clusters = 10
# THIS FOLDER WILL SAVE THE CONTACT MAP OF THE CLUSTERED AVERAGE, AND THE VALUES USER CAN USE AT LATER STAGE FOR CHANGING FORMATING AND COLORS.
# IT ALSO SAVES A cluster_representatives.txt FILE WHICH HAS NUMBER OF FRAMES INSIDE EACH CLUSTER AND ALSO THE REPRESENTATIVE FRAME OF THE CLUSTER FOR VISUALIZATION PRUPOSES
output_dir = 'cluster_output'
xtc_files = [
   '/scratch/ykv210/SPRINGS/BBBB/75/40/monomer/md/noPBC.xtc',
   ]

topology = '../md/eq.pdb'

#######################################################################################################

#-LETS LOAD TRAJECTOIRY USING MDTRAJ
print('Loading trajectory...')
traj = md.load(xtc_files, top=topology)
print(f'Total frames loaded: {len(traj)}')

###########################################################################
###########################################################################
#  Filter frames based on COM FOR DIMER TRAJECTORY - UNCOMMENT THIS BLOCK IF DIFFERENT CHAINS, AND HAVE ALREADY SAVED COM DISTANCES LIKE HERE
data = []
filez = ['../md/real_com.xvg', '../md1/real_com.xvg', '../md2/real_com.xvg']

for filename in filez:
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('@', '#')):
                parts = line.strip().split()
                if len(parts) >= 2:
                    data.append(float(parts[1]))  # second column

com_array = np.array(data)
frame_mask = com_array < 3.7   # only frames where COM < 3.7
print(f"Frames satisfying COM<3.7: {np.sum(frame_mask)} / {len(com_array)}")

# Apply mask to trajectory
traj = traj[frame_mask]
print(f"Trajectory after COM filtering: {traj.n_frames} frames")

###########################################################################
###########################################################################










residues = list(traj.topology.residues)
n_res = len(residues)
print(f'Topology residues found: {n_res if (n_res := n_res if False else n_res) else len(residues)}')

# Ensure chain ranges are inside available residues
max_res_index = n_res - 1
A_start, A_stop = chainA_range.start, chainA_range.stop
B_start, B_stop = chainB_range.start, chainB_range.stop
if A_stop - 1 > max_res_index or B_stop - 1 > max_res_index:
    raise ValueError('Specified chain residue ranges exceed the number of residues in topology.')

# collect sidechain atom indices per residue
sc_atoms_per_res = []
for res in residues:
    sc_atoms = [atom.index for atom in res.atoms if atom.name in sc_names]
    sc_atoms_per_res.append(sc_atoms)



# ------------------------ Build pair list (within-chain) ------------------------
pairs = []
pair_res = []

# IVE CONSIDERED RESIDUE PAIRS WHICH ARE MORE THAN 3 RESIDUES AWAY THAN THE INDEX i RESIDUE, SO AS TO REMOVE ANY CLOSE RESIDUE CONTATCS, LIMITING THE NOISE.
# ONE SHOULD CHANGE THIS +4 TO +1 IF CONSIDERING DIFFERENT CHAINS


for i in chainA_range:
    for j in range(i + 4, chainA_range.stop):  # FOR SINGLE CHAIN
    #for j in range(i + 1, chainA_range.stop):  #FOR DIFFERENT CHAINS
        atoms_i = sc_atoms_per_res[i]
        atoms_j = sc_atoms_per_res[j]
        if atoms_i and atoms_j:
            for a in atoms_i:
                for b in atoms_j:
                    pairs.append([a, b])
            pair_res.append((i, j))

pairs = np.array(pairs, dtype=int)
print(f'Total atom pairs (within-chain): {len(pairs)}')

#  Compute distances - THIS IS THE MOST RESOURCE HEAVY LINES OF THE CODE
print('Computing distances for all atom pairs...')
distances = md.compute_distances(traj, pairs)  # shape (n_frames, n_pairs)

# ATOM PAIR DISTANCES TO RESIDUE RESIDUE DISTANCES, BY TAKING THE MINIMUM DISTANCES OUT OF ALL THE SIDE CHAIN BEADS 

n_A = len(chainA_range)
n_B = len(chainB_range)
n_frames = traj.n_frames
min_dist_matrix = np.full((n_frames, n_A, n_B), np.inf)

start = 0
for (i_res, j_res) in pair_res:
    ni = len(sc_atoms_per_res[i_res])
    nj = len(sc_atoms_per_res[j_res])
    block = distances[:, start:start + ni * nj]
    min_block = np.min(block, axis=1)
    iA = i_res - chainA_range.start
    jB = j_res - chainB_range.start
    min_dist_matrix[:, iA, jB] = min_block
    start += ni * nj

min_dist_matrix[~np.isfinite(min_dist_matrix)] = large_r

#  Convert distances -> probabilities using sigmoidal function 
print('Converting distances to probabilities using sigmoidal function...')
# p = 1 / (1 + exp(beta * (r - r0)))
prob_matrix = 1.0 / (1.0 + np.exp(beta * (min_dist_matrix - r0)))


# THIS FLATTENING IS BEING DONE SO AS TO DO KMEANS CLUSTERING EASILY.

Y = prob_matrix.reshape(n_frames, n_A * n_B)
print('Feature matrix shape for clustering:', Y.shape)
n_features = Y.shape[1]

Y = np.array(Y, dtype=float)
Y[~np.isfinite(Y)] = 0.0

## K MEANS CLUSTERING DOING HERE
print(f'Running KMeans with n_clusters={n_clusters} ...')
kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init='auto')
labels = kmeans.fit_predict(Y)

#  Cluster averages (in probability space)
cluster_avgs = np.zeros((n_clusters, n_features))
for k in range(n_clusters):
    idx = np.where(labels == k)[0]
    if len(idx) == 0:
        print(f'Warning: cluster {k} is empty. Filling with zeros.')
        cluster_avgs[k] = np.zeros(n_features)
    else:
        cluster_avgs[k] = Y[idx].mean(axis=0)

# Save cluster average maps 

for k in range(n_clusters):
    avg_vec = cluster_avgs[k]
    avg_mat = avg_vec.reshape(n_A, n_B)
 
    outpath = os.path.join(output_dir, f'cluster_avg_{k}.txt')
    np.savetxt(outpath, avg_mat, fmt='%.6e')
    print(f'Saved cluster average {k} → {outpath}')


# Representative frames
centers = kmeans.cluster_centers_  
rep_frames = []
for k in range(n_clusters):
    center = centers[k]
    dists = np.linalg.norm(Y - center, axis=1)
    rep_idx = int(np.argmin(dists))
    rep_frames.append(rep_idx)

# Save representative frames
repfile = os.path.join(output_dir, 'cluster_representatives.txt')
with open(repfile, 'w') as fh:
    fh.write('# cluster_id  representative_frame_index  n_members\n')
    for k in range(n_clusters):
        members = int(np.sum(labels == k))
        fh.write(f'{k:2d}  {rep_frames[k]:8d}  {members:6d}\n')
print(f'Saved representative frames → {repfile}')

#  Save frame numbers for each cluster 
print("Saving frame numbers for each cluster...")

for k in range(n_clusters):
    frames = np.where(labels == k)[0]  # all frame indices for this cluster
    outpath = os.path.join(output_dir, f"cluster_{k}_frames.txt")
    np.savetxt(outpath, frames, fmt='%d')
    print(f"  Saved {len(frames)} frames → {outpath}")



#  Plotting: separate heatmaps for each cluster 
print('Plotting and saving individual cluster-average heatmaps...')

# Use perceptually uniform light→dark colormaps
cmaps = ['Greens', 'Blues', 'Oranges', 'Purples', 'Reds',
         'YlGn', 'YlOrBr', 'PuBuGn', 'BuPu', 'GnBu']


for k in range(n_clusters):
    avg_vec = cluster_avgs[k]
    avg_mat = avg_vec.reshape(n_A, n_B)
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        avg_mat, vmin=0.0, vmax=1.0,
        origin='lower', interpolation='nearest', cmap=cmaps[k % len(cmaps)]
    )
    
    ax.set_title(f'Cluster {k} (n={int(np.sum(labels == k))})', fontsize=14)
    ax.set_xlabel('Residue index', fontsize=12)
    ax.set_ylabel('Residue index', fontsize=12)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Contact Probability')
    
    fig.tight_layout()
    outpath = os.path.join(output_dir, f'cluster_avg_heatmap_{k}.png')
    plt.savefig(outpath, dpi=300)
    plt.close(fig)
    print(f'Saved heatmap for cluster {k} → {outpath}')


print('All individual cluster heatmaps saved successfully.')



#  Print summary 
print('\nSummary:')
for k in range(n_clusters):
    members = int(np.sum(labels == k))
    print(f' - Cluster {k:2d}: members={members:4d}, representative_frame={rep_frames[k]}')

print('\nDone. All outputs are in the folder:', output_dir)



#FEEL FREE TO CONTACT ME AT Y.K.Verma@sms.ed.ac.uk for any query




