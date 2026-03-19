#!/bin/bash
# edt5min.sh: Minimization of edt.* on Rockfish
#SBATCH -J min107      # job name
#SBATCH -o %j.out      # output and error file name (%j expands to jobID)
#SBATCH -p defq        # Queue name: defq or a100
#SBATCH -N 1           # number of nodes
#SBATCH --ntasks=8     # CPU cores
#SBATCH -A mcb140208_gpu
#SBATCH -t 24:00:00    # run time (hh:mm:ss) - 48 hour (max)

# Load environment variables
# This is for sander.MPI on defq nodes
source /home/cannonj/amber18/amber.sh 

# Range of models to process
#for ((i=186; i<=200;i++)); do
for i in 107 109 110 111 112 113 114 115 116 117 118 119 120; do
name=edt.$i

# Make *.top and *.crd files
cat<<eof>parmin
chamber -top $HOME/charm/toppar/top_all36_prot.rtf \
-param $HOME/charm/toppar/par_all36m_prot.prm \
-top $HOME/charm/toppar/top_all36_lipid.rtf \
-param $HOME/charm/toppar/par_all36_cgenff.prm \
-param $HOME/charm/toppar/par_all36_carb.prm \
-param $HOME/charm/toppar/par_all36_lipid.prm \
-str $HOME/charm/toppar/toppar_all36_lipid_cholesterol.str \
-str $HOME/charm/toppar/toppar_all36_carb_glycolipid.str \
-str $HOME/charm/toppar/toppar_all36_lipid_inositol.str \
-str $HOME/charm/toppar/toppar_water_ions.str \
-psf  $name.psf \
-box bounding \
-crd $name.pdb
parmout $name.top $name.crd
# The -box bounding ensures *.crd has box information like it would from tleap.
#outPDB edt4.$mod.pdb
# Note the TER cards increment ATOM number in this context!
go
eof
parmed -n -i parmin

# Minimize solvent first by restraining protein and membrane
cat<<eof>mdin1
Minimization with protein and lipid restraints
 &cntrl
    ! Minimization options
    imin=1,       ! Turn on minimization
    maxcyc=1000,  ! Maximum number of minimization cycles
    ncyc=500,     ! shift to conjugant gradient halfway
    
    ! Potential energy function options
    cut=12.0,      ! nonbonded cutoff, in angstroms
    fswitch=10.0,  ! Force-based switching

    ! Control how often information is printed to the output file
    ntpr=100,      ! Print energies every 100 steps
    ntxo=2,        ! Write NetCDF format

    ! Restraint options
    ntr=1,         ! Positional restraints for proteins and lipid head groups

    ! Set water atom/residue names for SETTLE recognition
    watnam='WAT',  ! Water residues are named WAT
    owtnm='O',     ! Water oxygens are named O
 /
 &ewald
    vdwmeth = 0,
 /
Protein position restraint
10.0
RES 1 74 75 157
END
Membrane position restraint
2.5
FIND
O3 * * CHL1
P * * POPS
P * * DPPC
P * * POPC
P * * DOPE
P * * POPI
P * * DOPC
SEARCH
RES 1 88273
END
END
eof
# First minimization
mpirun -n 8 sander.MPI -O -i mdin1 -p $name.top -c $name.crd \
-ref $name.crd \
-r $name.min.rst -o $name.min.out -inf $name.inf 
#
done

