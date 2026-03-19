#!/bin/bash
# MakeEpsilon.sh: Create random CD3epsilon cytoplamsic domain models.
# This is a self-contained bash job script.
# These have the humanized 57-residue CD3e CT sequence from PDB 2K4F.
# These  models will use the ff14IDPSFF force field and 
# the PRS domain will be restrained to the Nck-bound conformation.
# 250 models will be equilibrated for 5 ns with chiral, trans omega, and PRS
# restraints. 
# Save only final frame after 5 ns MD in binary restart format.      
# Before this job script was run, getPRSrestraints.sh 
# made eps3.RST with chiral, omega, and PRS restraints and
# eps4.top topology using the ff14IDPSFF force field.
# It also made the PRS impose commands.
# temp.top=ff14SB force field, eps4.top=ff14IDPSFF force field.
topology=eps4.top
# eps3.RST=chiral, omega, and PRS restraints; eps.RST=chiral and omega
restraints=eps3.RST
# Output table
table=epsTable
rm -f $table

# CD3epsilon CT 57-residues in 2k4f, converted mouse to human
# Start eps at CD3e Trp151, end at C-terminus Ile207.      
cat<<eof>eps.seq
hI2seg = sequence { NTRP SER LYS ASN ARG LYS ALA LYS 
ALA LYS PRO VAL THR ARG GLY ALA GLY ALA GLY GLY ARG 
GLN ARG GLY GLN ASN LYS GLU ARG PRO PRO PRO VAL PRO 
ASN PRO ASP TYR GLU PRO ILE ARG LYS GLY GLN ARG ASP
LEU TYR SER GLY LEU ASN GLN ARG ARG CILE }
eof
# tleap made eps.pdb from above and then makeCHIR_RST eps.pdb 
# made eps.RST, which has chiral and omega restraints.
# Need to add PRS dihedral restraints.
# Need to make tleap impose commands to set the PRS dihedrals.

## There are three pmemd input files necessary: one for minimization, one for SA,
# one for equilibration.
cat <<eof >mdmin
Initial 1000 steepest descent step minimization
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  cut=999., rgbmax=999.,igb=1, ntb=0,
  ntpr=100
 /
eof
####################################################################
cat <<eof >SAin
150 ps simulated annealing protocol
 &cntrl
  nstlim=150000, ntt=1,
  ntb=0, igb=1,
  cut=999.,rgbmax=999.
  ntpr=500,
  nmropt=1,
 &end
#
# Simulated annealing algorithm:
#
#   from steps       0 to  1000: raise target temperature 10->1500K
#   from steps    1000 to 50000: leave at 1500K
#   from steps  50000 to 150000: cool to 300K
#
 &wt type='TEMP0',istep1=0, istep2=1000,value1=10.,value2=1500.,&end
 &wt type='TEMP0',istep1=1001, istep2=50000, value1=1500.,value2=1500.0,&end
 &wt type='TEMP0',istep1=50001, istep2=150000, value1=300,value2=300,&end
#
# Strength of temperature coupling:
#
#    steps      0 to  50000: tight coupling for heating and equilibration
#    steps  50000 to 130000: slow cooling phase
#    steps 130000 to 150000: somewhat faster cooling
#
 &wt type='TAUTP',istep1=0,istep2=3000,value1=0.2,value2=0.2,&end
 &wt type='TAUTP',istep1=3001,istep2=50000,value1=4.0,value2=2.0,&end
 &wt type='TAUTP',istep1=50001,istep2=130000,value1=1.0,value2=1.0,&end
 &wt type='TAUTP',istep1=130001,istep2=150000,value1=0.5,value2=0.5,&end

 &wt type='END',  &end
LISTOUT=POUT
DISANG=$restraints 
 /
eof
####################################################################
cat <<eof >Eqin
 5 ns NPT at 300K, restraints, only last restart saved
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=2500000, dt=0.002,
  ntc=2, ntf=2,
  ntt=3, gamma_ln=1, ig=-1,
  ntpr=500, ntwx=0,
  cut=999.,rgbmax=999. 
  igb=5, saltcon=0.2, ntb=0,
  ntxo=2,ioutfm=1,
  nmropt=1,pencut=0.1,
 &end
 &wt type='END',  &end
LISTOUT=POUT
DISANG=$restraints
 /
eof
####################################################################
for ((i=1; i<=250;i++)); do
# Make tleap script, the HGBLE string will be in a comment
echo "source leaprc.protein.ff14SB" > leapin
echo "source eps.seq" >> leapin
echo "set default PBradii mbondi2" >> leapin
# Note backticks!
HGBLEstring=`makeHGBLE eps.seq HGBLEfreq2 `
makeTorsions $HGBLEstring >> leapin
# Add PRS impose commands
cat epsStart.txt>>leapin
echo "saveAmberParm hI2seg temp.top temp.crd" >> leapin
echo "quit" >> leapin
# Run tleap
tleap -f leapin >leapout
# Run sander to minimize
mpirun -np 8 sander.MPI -O -i mdmin -p $topology -c temp.crd \
-r min.rst -o min.out
# Check if NaN energy resulting from steric clashes in mdinfo.
if grep -q  NaN mdinfo
then continue
fi
# Run simulated annealing
pmemd.cuda -O -i SAin -p $topology -c min.rst \
-r SA.rst -o SA.out
# Get potential energy from mdinfo. Note backticks!
high=`grep -c "*" mdinfo `
if [ $high -gt 0 ]
then continue
fi
energy=`grep Etot mdinfo | awk '{print $3}' `
# Do not consider models with high energy
# Pipe statement to floating point calculator, bc, to compare. Special chars need quotes.
highE=`echo $energy '>' -500 | bc -l `
# The spaces are necessary!
if [ $highE -eq 1 ] 
then continue
fi
# Write to table using enabled blackslash escapes
# Model# penalty energy HGBLEstring
echo -e $i "\t" $energy "\t" $HGBLEstring >> $table
# Now start 5 ns equilibration.
pmemd.cuda -O -i Eqin -p $topology -c SA.rst \
-r traj/eps.$i.rst -o eq.out 
# Only save the final restart.
rm leap.log    # Otherwise it gets big!
done

