#!/bin/bash
# runLck.sh: Run 100 ns of explicit solvent Lck-peptide model at 310 K.
# This has the Lck bound CD3ζ Tyr111 15-residue ITAM from Neel Shah.
cd lck
name=lckPep
# Trim away offending hydrogens
sed '/HH33/d;/HH32/d;/HH31/d;' lck.pdb >lck2.pdb
sed '/10  H   LYS/d' peptide.pdb>peptide2.pdb
# Step 1: build system
cat<<eof >leapin
# Use recommended force fields
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.DNA.OL15
loadamberparams atp.frcmod
loadamberprep atp2.prep
lck=loadpdb lck2.pdb
pep=loadpdb peptide2.pdb
atp=loadpdb ATP.pdb
lckPep=combine {lck pep atp}
check lckPep
# ^That gave charge of -7.
solvateOct lckPep TIP3PBOX 12.0 iso
# ^That added 13,041 waters
# Add counterions with solvent replacement. 
# Number of ions = 0.01867 * (Desired Concentration) * (Number of waters)
# =37
# 7 Na+ will neutralize charge, 44 extra bring NaCl to 150 mM.
addIons lckPep Na+ 44 Cl- 37
saveamberparm lckPep $name.top $name.crd
# Just to check
savepdb lckPep $name.out.pdb
quit
eof
tleap -f leapin
# Last protein atom 4654.
# Step 2: Minimize solvent
cat << eof >mdin1
Initial 1000 step solvent minimization
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  cut=9, ntb=1,
  ntr=1,
  ntpr=100,
 /
Hold the protein fixed
500.0
RES 1 4654
END
END
eof
pmemd.cuda -O -i mdin1 -p $name.top -c $name.crd \
-ref $name.crd \
-r $name.min.rst -o $name.min.out -inf $name.inf 
# Step 3: Minimize everything
cat << eof >mdin2
Total 1000 step  minimization
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  cut=9, ntb=1,
  ntpr=100,
 /
END
eof
pmemd.cuda -i mdin2 -p $name.top -c $name.min.rst \
-r $name.min2.rst -o $name.min2.out -inf $name.inf 
# Step 4: Position-restrained 20 ps NVT dynamic heating
cat << eof >mdin3
NVT position-restrained 20 ps MD heating
 &cntrl
  imin=0, irest=0, ntx=1,
  cut=9, ntb=1,
  ntr=1,
  ntc=2, ntf=2,
  ntt=1, tempi=0.0,
  nstlim = 10000, dt = 0.002,
  ntpr=100,
  nmropt=1, restraint_wt=10.0, restraintmask=':1-4654',
 &end
 &wt type='TEMP0', istep1=0, istep2=10000, value1=0.0, value2=310.0, /
 &wt type = 'END', /
 LISTOUT=POUT
&wt type='END', &end
 /
eof
pmemd.cuda -O -i mdin3 -p $name.top -c $name.min2.rst \
-ref $name.min2.rst \
-r $name.eq1.rst -o $name.eq1.out -inf $name.inf >& /dev/null 
# Step 5: 10 ns NPT unrestrained
cat << eof >mdin4
Equilibration for 10 ns at 310K, save every 10 ps.
 &cntrl
  imin=0, irest=1, ntx=7,
  ntb=2, pres0=1.0, ntp=1, taup=2.0,
  cut=9, ntr=0,
  ntc=2, ntf=2,
  tempi=310.0, temp0=310.0,
  ntt=1,
  nstlim=5000000, dt=0.002, ioutfm=1, ntwprt=4654
  ntpr=5000, ntwx=5000,
 /
END
eof
# Run 100 ns NPT MD in 10 ns intervals
for ((i=1;i<=10;i++)); do
((j=i+1))
pmemd.cuda -O -i mdin4 -p $name.top -c $name.eq$i.rst \
-r $name.eq$j.rst -o $name.eq$j.out -x $name.eq$j.cdf \
-inf $name.inf </dev/null
done


