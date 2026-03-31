#!/bin/bash
# tcr39NPT3a.sh: Run 40 ns NPT of tcr39.*.** models on Anvil
#SBATCH -J 39.5.02       # job name
#SBATCH -o %j.out        # output and error file name (%j expands to jobID)
#SBATCH -p gpu           # Queue (partition) name: gpu
#SBATCH --nodes=1             # Total # of nodes 
#SBATCH --ntasks-per-node=1   # Number of MPI ranks per node (one rank per GPU)
#SBATCH --gpus-per-node=1     # Number of GPUs per node
#SBATCH -A mcb140208-gpu  # allocation name
#SBATCH -t 48:00:00        # run time (hh:mm:ss) - 48 hour (max)
# Load environment variables
module --force purge
module load modtree/gpu
module purge
module load gcc/8.4.1 openmpi/4.0.6-cu11.0.3
module load amber

#name=tcr39.0.00
name=tcr$SLURM_JOB_NAME
topology=$name.top
# Starting eq*, increase by 7 from last one
i=45
((j=i+1))

# 5 ns production NPT 
cat<<eof>mdin5
5 ns NPT 310 K without restraints, 10 ps frames
 &cntrl
  imin=0, irest=1, ntx=5, 
  ntt=3, gamma_ln=1.0, temp0=310.0, ig=-1,
  cut=12.0, fswitch=10.0,
  nstlim=2500000, dt=0.002,
  ntc=2, ntf=2,     ! SHAKE
  ntpr=5000, ntwr=5000, ntwx=5000,
  iwrap=1,          ! Keep atoms in periodic box
  ntxo=2, ioutfm=1, ! Save in NetCDF format 
  ntwprt=11037      ! Save only 11037 protein, CHL1 519, CHL1 201 atoms
  barostat=2,       ! MC barostat
  ntp=3,            ! semi-isotropic scaling
  pres0=1.0,        ! 1 bar pressure
  csurften=3,       ! Interfaces in xy plane
  gamma_ten=0.0,    ! No surface tension for pure semi-iso scaling
  ninterface=2,     ! Number of interfaces (2 for bilayer)
  ! Set water atom/residue names for SETTLE recognition
  watnam='WAT',     ! Water residues are named WAT
  owtnm='O',        ! Water oxygens are named O
 /
 &ewald
    vdwmeth = 0,
 /
  &wt
    type='END'
 /
END
eof
########################################################### 
# Eight 5 ns blocks
# Removing the ">& /dev/null" to catch potential errors.
date
echo -n "Running on "; hostname
echo "Output $name.eq$j.cdf"
pmemd.cuda -O -i mdin5 -p $topology -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf 
((i=i+1))
((j=i+1))
echo "Output $name.eq$j.cdf"
pmemd.cuda -O -i mdin5 -p $topology -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf 
((i=i+1))
((j=i+1))
echo "Output $name.eq$j.cdf"
pmemd.cuda -O -i mdin5 -p $topology -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf 
((i=i+1))
((j=i+1))
echo "Output $name.eq$j.cdf"
pmemd.cuda -O -i mdin5 -p $topology -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf 
((i=i+1))
((j=i+1))
echo "Output $name.eq$j.cdf"
pmemd.cuda -O -i mdin5 -p $topology -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf 
((i=i+1))
((j=i+1))
echo "Output $name.eq$j.cdf"
pmemd.cuda -O -i mdin5 -p $topology -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf 
((i=i+1))
((j=i+1))
echo "Output $name.eq$j.cdf"
pmemd.cuda -O -i mdin5 -p $topology -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf
((i=i+1))
((j=i+1))
echo "Output $name.eq$j.cdf"
pmemd.cuda -O -i mdin5 -p $topology -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf -inf $name.inf 
date

 
