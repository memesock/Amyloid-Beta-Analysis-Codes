
# ONE MUST HAVE AN ALL ATOMISTIC PDB STRUCTURE TO START THE SIMULATION.IN THIS SCRIPT, ITS DENOTED BY dimer.pdb

# FOR MARTINI 3 SIMULATION
martinize2 -f dimer.pdb -o topol.top -x protein_cg.pdb -p backbone -ff martini3001 -cys auto 

#FOR MARTINI 3 PROTEIN-WATER SCALED SIMULATION, ONE MUST CHANGE THE EPSILON VALUE FOR LJ POTENTIAL BETWEEN PROTEIN AND WATER ATOMS.
# THE martini_v3.0.0/martini_v3.0.0.itp FILE CONTAINS THE 




cat > topol.top <<EOF
#include "martini_v3.0.0/martini_v3.0.0.itp"
#include "martini_v3.0.0/martini_v3.0.0_ions_v1.itp"
#include "martini_v3.0.0/martini_v3.0.0_solvents_v1.itp"
#include "protein0.itp"
#include "protein1.itp"
