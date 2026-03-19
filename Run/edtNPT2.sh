#!/bin/bash
# edtNPT2.sh: continue 5 ns NPT of edt models on Rockfish
#SBATCH -J edt91-100    # job name
#SBATCH -o %j.out         # output and error file name (%j expands to jobID)
#SBATCH -p a100           # Queue name: defq or a100
#SBATCH -N 1
#SBATCH --cpus-per-gpu=1
#SBATCH --gres=gpu:1
#SBATCH -x gpu01          # Exclude the V100 node, gpu01
#SBATCH -A mcb140208_gpu
#SBATCH -t 48:00:00       # run time (hh:mm:ss) - 48 hour (max)
# Load environment variables
module load cuda/11.1.0
source /home/cannonj/amber20_src/amber.sh

# Range of models to process
for ((i=91;i<=100;i++)); do
name=edt.$i
echo "Running continuation of 5 ns NPT MD for $name"
# 5 ns production NPT continuation
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
  ntwprt=2448       ! Save only 2448 protein atoms
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
END
eof
#
pmemd.cuda -O -i mdin5 -p $name.top -c $name.eq2.rst \
-r $name.eq3.rst -o $name.eq3.out -x $name.eq3.cdf -inf $name.inf 

done

