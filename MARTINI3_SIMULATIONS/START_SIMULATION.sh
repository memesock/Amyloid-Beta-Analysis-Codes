
# ONE MUST HAVE AN ALL ATOMISTIC PDB STRUCTURE TO START THE SIMULATION.IN THIS SCRIPT, ITS DENOTED BY dimer.pdb

#######################################################################
##############################################################
######################################################
# FOR MARTINI 3 SIMULATION
martinize2 -f dimer.pdb -o topol.top -x protein_cg.pdb -p backbone -ff martini3001 -cys auto -merge all

#FOR MARTINI 3 PROTEIN-WATER SCALED SIMULATION, ONE MUST CHANGE THE EPSILON VALUE FOR LJ POTENTIAL BETWEEN PROTEIN AND WATER ATOMS. THE martini_v3.0.0/martini_v3.0.0.itp FILE CONTAINS THE PROTEIN-WATER VALUES, WHICH CAN BE MODIFIED USING PYTHON FILE INSIDE THE FOLDER "LAMBDA".
martinize2 -f dimer.pdb -o topol.top -x protein_cg.pdb -p backbone -ff martini3001 -cys auto -merge all

#FOR MARTINI 3 IDP SIMULATION, RUN FOLLOWWING COMMAND
martinize2 -f dimer.pdb -o topol.top -x protein_cg.pdb -ff martini3IDP -cys auto -idr-tune -id-regions 1:42 -merge all

######################################################
##############################################################
#######################################################################

# UPTO NOW ONE MUST HAVE A TOPOL.TOP FILE AND A COARSE GRAINED PDB FILE, AND AN .ITP FILE READY ACCORDING TO FORCE FIELD SELECTED READY WITHIN THE WORKING DIRECTORY.
# ONE MUST HAVE A MDP FOLDER AND martini_v3 FOLDER READY TOO. I HAVE INCLUDED THE DEFUALT LAMBDA OF 1 FOR PROTEIN-WATER SCALING IN martini_v3.0.0/martini_v3.0.0.itp, ONE CAN CHANGE IT USING THE LAMBDA DIRECTORY PYTHON FILE.

cat > topol.top <<EOF
#include "martini_v3.0.0/martini_v3.0.0.itp"
#include "martini_v3.0.0/martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0/martini_v3.0.0_solvents_v1.itp"
#include "protein0.itp"
EOF

#CHANGE THE BOX SIZE ACCORDINGLY
gmx editconf -f protein_cg.pdb -box 13.6 13.6 13.6 -bt cubic -o protein_cg.gro

#ENERGY MINIMIZATION IN VACCUM
mkdir -p em_vac
gmx grompp -f mdp_files/em.mdp -p topol.top -c protein_cg.gro -r protein_cg.gro -o em_vac/em_vac.tpr
gmx mdrun -v -deffnm em_vac/em_vac

#SOLVATION IN THE BOX 
mkdir -p sol
gmx solvate -cp em_vac/em_vac.gro -cs mdp_files/water.gro -radius 0.21 -o sol/sol.gro -p topol.top -maxsol 21000 

#ADDYING IONS IN THE SOLVATED BOX
gmx grompp -f mdp_files/ions.mdp -c sol/sol.gro -p topol.top -o sol/ions.tpr -maxwarn 1
echo W | gmx genion -s sol/ions.tpr -p topol.top -neutral -conc 0.15 -o sol/sol_neutral.gro

#ENERGY MINIMIZATION
mkdir -p em
gmx grompp -f mdp_files/em.mdp -p topol.top -c sol/sol_neutral.gro -r sol/sol.gro -o em/em.tpr
gmx mdrun -v -deffnm em/em

#EQUILLIBRIATION
mkdir -p eq
gmx grompp -f mdp_files/eq.mdp -c em/em.gro -p topol.top -o eq/eq.tpr -maxwarn 2
srun gmx mdrun -v -deffnm eq/eq -nb gpu -update gpu

#MD RUN
mkdir -p md
gmx grompp -f mdp_files/md.mdp -c eq/eq.gro -p topol.top -o md/md.tpr
srun gmx mdrun -v -deffnm md/md -nb gpu -update gpu
















