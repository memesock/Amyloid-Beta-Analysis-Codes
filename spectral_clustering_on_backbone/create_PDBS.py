#!/usr/bin/env python3

import os
import numpy as np
import MDAnalysis as mda

# ======================
# Input files
# ======================

TRAJ_FILES = [
    '../md/pbc_emma.xtc',
    '../md1/pbc_emma1.xtc',
    '../md2/pbc_emma2.xtc'
]

TOPO_PDB = "../md/eq.pdb"

COM_FILES = [
    '../md/real_com.xvg',
    '../md1/real_com.xvg',
    '../md2/real_com.xvg'
]

OUT_DIR = "pdb_10_stride_com_2.5"
STRIDE = 10
COM_CUTOFF = 2.5

# ======================
# Load trajectory
# ======================

print("Loading trajectory...")
u = mda.Universe(TOPO_PDB, TRAJ_FILES)
traj = u.trajectory
n_frames = traj.n_frames
print(f"Total frames loaded: {n_frames}")

# ======================
# Load COM data
# ======================

print("Loading COM data...")
data = []

for filename in COM_FILES:
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(('@', '#')):
                parts = line.strip().split()
                if len(parts) >= 2:
                    data.append(float(parts[1]))

com_array = np.array(data)

if len(com_array) != n_frames:
    raise RuntimeError(
        f"ERROR: COM distances ({len(com_array)}) do not match loaded frames ({n_frames})"
    )

# ======================
# Frame filtering
# ======================

frame_mask = com_array < COM_CUTOFF
valid_frames = np.where(frame_mask)[0]

print(f"Frames satisfying COM < {COM_CUTOFF}: {len(valid_frames)} / {n_frames}")

# ======================
# Apply stride (on filtered frames)
# ======================

strided_frames = valid_frames[::STRIDE]

print(f"Frames after stride {STRIDE}: {len(strided_frames)}")

# ======================
# Backbone atom selection
# ======================

# Standard protein backbone atoms
bb = u.select_atoms("name BB")

if len(bb) == 0:
    raise RuntimeError("ERROR: No backbone atoms selected. Check atom names.")

# ======================
# Output directory
# ======================

os.makedirs(OUT_DIR, exist_ok=True)

# ======================
# Write PDBs
# ======================

print("Saving PDB files...")

for frame_idx in strided_frames:
    traj[frame_idx]  # move to frame

    pdb_name = f"frame_{frame_idx}.pdb"
    pdb_path = os.path.join(OUT_DIR, pdb_name)

    bb.write(pdb_path)

print("Done.")
print(f"PDBs saved in: {OUT_DIR}")

