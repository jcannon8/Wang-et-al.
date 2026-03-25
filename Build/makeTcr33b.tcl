# makeTcr33b.tcl: Make variation of tcr33 with S1 CHL from fgt.185.
# vmd -dispdev none on Eire
cd /t/tcr33
set outMod tcr33b
# tcr30.psf is the same as /t/tcr30/charmm-gui-7001280423/amber/step5_input.psf
# made by CHARMM-GUI from 7FJD.
set tcr [mol new /t/tcr30/tcr30.psf]
mol addfile /t/tcr30/tcr30.pdb waitfor all
# Assign chain names, 8 protein subunits, 2 structural CHL, 1600 membrane lipids, 
# CLA and SOD ions, and water
[atomselect $tcr "segname PROA"] set chain A
[atomselect $tcr "segname PROB"] set chain B
[atomselect $tcr "segname PROC"] set chain D
[atomselect $tcr "segname PROD"] set chain E
[atomselect $tcr "segname PROE"] set chain F
[atomselect $tcr "segname PROF"] set chain G
[atomselect $tcr "segname PROG"] set chain M
[atomselect $tcr "segname PROH"] set chain N
[atomselect $tcr "segname MEMB"] set chain J
# Now chain M is just TCRalpha protein atoms
set tcr12 [mol new /home/jcannon/charm/charmm-gui-2876761245/amber/tcr12.psf]
mol addfile /home/jcannon/charm/charmm-gui-2876761245/amber/tcr12.pdb
[atomselect $tcr12 "segname PROD"] set chain F
# Confirm missing four chain F residues
#[atomselect $tcr "name CA and chain F"] get {resname resid}
#[atomselect $tcr12 "name CA and chain F"] get {resname resid}
set F30 [atomselect $tcr "name CA and chain F and (resid 33 to 69 or resid 74 to 156)"] 
set F12 [atomselect $tcr12 "name CA and chain F and (resid 33 to 69 or resid 74 to 156)"] 
# Move tcr12 chain F to fit tcr30 
set M [measure fit $F12 $F30]
set F12b [atomselect $tcr12 "chain F"]
$F12b move $M
measure rmsd $F12 $F30; # =0.9731
# Save parts of chain F
[atomselect $tcr "chain F and resid 22 to 69"] writepdb fchainN.pdb
[atomselect $tcr12 "chain F and resid 70 to 74"] writepdb fchainX.pdb
[atomselect $tcr "chain F and resid 75 to 156"] writepdb fchainC.pdb
# Fuse chain F together using bash commands.
set out newFchain.pdb
exec cat fchainN.pdb > $out
exec cat fchainX.pdb >> $out
exec cat fchainC.pdb >> $out
# Remove END and CRYST records.
exec sed -i /END/d $out
exec sed -i /CRYST/d $out
#############################
# Extend chain A by three residues from 22-54 to 22-57
[atomselect $tcr12 "segname PROA"] set chain A
# Use five CA to position
set A30 [atomselect $tcr "name CA and chain A and resid 50 to 54"]
set A12 [atomselect $tcr12 "name CA and chain A and resid 50 to 54"] 
# Move tcr12 chain A to fit tcr30 
set M [measure fit $A12 $A30]
set A12b [atomselect $tcr12 "chain A"]
$A12b move $M
measure rmsd $A12 $A30; # =0.3844 
# Save parts of chain A
[atomselect $tcr "chain A and resid 22 to 54"] writepdb AchainN.pdb
[atomselect $tcr12 "chain A and resid 55 to 57"] writepdb AchainC.pdb
# Fuse chain A together using bash commands.
set out newAchain.pdb
exec cat AchainN.pdb > $out
exec cat AchainC.pdb >> $out
# Remove END and CRYST records.
exec sed -i /END/d $out
exec sed -i /CRYST/d $out
#############################
# Extend chain B by N- and C-terminal residues from 26-54 to 24-55
[atomselect $tcr12 "segname PROB"] set chain B
# Use all of 7FJD CA to position
set B30 [atomselect $tcr "name CA and chain B and resid 26 to 54"]
set B12 [atomselect $tcr12 "name CA and chain B and resid 26 to 54"] 
# Move tcr12 chain B to fit tcr30 
set M [measure fit $B12 $B30]
set B12b [atomselect $tcr12 "chain B"]
$B12b move $M
measure rmsd $B12 $B30; # =1.488
# Save parts of chain B
[atomselect $tcr12 "chain B and resid 24 to 26"] writepdb BchainN.pdb
[atomselect $tcr "chain B and resid 27 to 54"] writepdb BchainX.pdb
[atomselect $tcr12 "chain B and resid 55"] writepdb BchainC.pdb
# Fuse chain B together using bash commands.
set out newBchain.pdb
exec cat BchainN.pdb > $out
exec cat BchainX.pdb >> $out
exec cat BchainC.pdb >> $out
# Remove END and CRYST records.
exec sed -i /END/d $out
exec sed -i /CRYST/d $out
#############################
# Extend chain D by one residue
[atomselect $tcr12 "segname PROC"] set chain D
# Use five CA to position
set D30 [atomselect $tcr "name CA and chain D and resid 124 to 128"]
set D12 [atomselect $tcr12 "name CA and chain D and resid 124 to 128"] 
# Move tcr12 chain D to fit tcr30 
set M [measure fit $D12 $D30]
set D12b [atomselect $tcr12 "chain D"]
$D12b move $M
measure rmsd $D12 $D30; # =0.365
# Save parts of chain D
[atomselect $tcr "chain D and resid 22 to 128"] writepdb DchainN.pdb
[atomselect $tcr12 "chain D and resid 129"] writepdb DchainC.pdb
# Fuse chain D together using bash commands.
set out newDchain.pdb
exec cat DchainN.pdb > $out
exec cat DchainC.pdb >> $out
# Remove END and CRYST records.
exec sed -i /END/d $out
exec sed -i /CRYST/d $out
#############################
# Save all parts not changed
[atomselect $tcr "chain E"] writepdb eChain.pdb
[atomselect $tcr "chain G"] writepdb gChain.pdb
[atomselect $tcr "chain M"] writepdb mChain.pdb
[atomselect $tcr "chain N"] writepdb nChain.pdb
# Save two structural cholesterols separately
#[atomselect $tcr "segid HETA"] writepdb s1Chain.pdb
# Use the S1 CHL519.pdb from anTcr33.sh, which came from fgt.185.
set S2 [atomselect $tcr "segid HETB"]
$S2 set chain S
$S2 writepdb s2Chain.pdb
[atomselect $tcr "chain J"] writepdb mem.pdb
# Ion chain names already assigned
#[atomselect $tcr "chain I"] writepdb ions.pdb
# Remove 2 CLA to neutralize
set resCLA [[atomselect $tcr "name CLA"] get resid]
set saveCLA [lrange $resCLA 0 end-2]
# Now save fewer CLA
[atomselect $tcr "name SOD or (name CLA and resid $saveCLA)"] writepdb ions.pdb
# The waters need to be saved in smaller chunks for psfgen.
set saveWat [atomselect $tcr water]
$saveWat set chain z
set Res [lsort -unique -integer [$saveWat get residue]]
set last [lrange $Res end end]
# A PDB residue can only have 4-characters
# Save water in 9999 residue chunks
set chunk 9999
for {set n 1; set k 0} {$k<=$last} {incr n} {
  set j [expr ($n-1) * $chunk +1]
  set k [expr $j + $chunk -1]
  puts "$n $j $k"
  set toSave [atomselect $tcr "water and residue $j to $k"]
  # Need to adjust residue# so that they go from 1 to $chunk in each file.
  # That ensures they will be unique within each psfgen segment.
  unset -nocomplain z
  set nRes [expr [$toSave num] / 3]
  for {set i 1} {$i<=$nRes} {incr i} {lappend z $i $i $i}
  $toSave set resid $z
  $toSave writepdb water$n.pdb
}
# water1.pdb to water[$n-1].pdb saved
###########################################################
# Use psfgen to build psf and pdb outputs.
package require psfgen
# Load the Charmm32m topology and parameters necessary for this system.
topology /home/jcannon/charm/toppar/top_all36_prot.rtf
topology /home/jcannon/charm/toppar/par_all36m_prot.prm
topology /home/jcannon/charm/toppar/top_all36_lipid.rtf
topology /home/jcannon/charm/toppar/par_all36_cgenff.prm
topology /home/jcannon/charm/toppar/par_all36_carb.prm
topology /home/jcannon/charm/toppar/par_all36_lipid.prm
topology /home/jcannon/charm/toppar/toppar_all36_lipid_cholesterol.str
topology /home/jcannon/charm/toppar/toppar_all36_carb_glycolipid.str
topology /home/jcannon/charm/toppar/toppar_all36_lipid_inositol.str
topology /home/jcannon/charm/toppar/toppar_water_ions.str
# The above reports many errors but they are OK (I hope so).
# Load segments with coordinates. 
# This defines the order in this model (same as 7FJD and tcr30 order)
# Even though these chain names are lower case, they are upper case in the output.
segment a {pdb newAchain.pdb}; coordpdb newAchain.pdb a
segment b {pdb newBchain.pdb}; coordpdb newBchain.pdb b
segment d {pdb newDchain.pdb}; coordpdb newDchain.pdb d
segment e {pdb eChain.pdb}; coordpdb eChain.pdb e
segment f {pdb newFchain.pdb}; coordpdb newFchain.pdb f
segment g {pdb gChain.pdb}; coordpdb gChain.pdb g
segment m {pdb mChain.pdb}; coordpdb mChain.pdb m
segment n {pdb nChain.pdb}; coordpdb nChain.pdb n
# Two structural cholesterols, membrane, ions
#segment s1 {pdb s1Chain.pdb}; coordpdb s1Chain.pdb s1
# Use the S1 CHL519.pdb from anTcr33.tcl, which came from fgt.185.
segment s1 {pdb CHL519.pdb}; coordpdb CHL519.pdb s1
segment s2 {pdb s2Chain.pdb}; coordpdb s2Chain.pdb s2
segment j {pdb mem.pdb}; coordpdb mem.pdb j
segment k {pdb ions.pdb}; coordpdb ions.pdb k
# Load in water
pdbalias residue WAT TIP3
pdbalias atom WAT O OH2
# Water is in multiple pdb files
for {set i 1} {$i<$n} {incr i} {
  segment w$i {pdb water$i.pdb}; coordpdb water$i.pdb w$i
  puts $i
}
# 14 disulfides in TCR counting the missing E:119 E:122 from 7fjd.pdb
patch DISU A:32 B:32
patch DISU D:37 D:73
patch DISU D:93 D:96
patch DISU E:49 E:98
# Add disulfide with missing CONECT record for E:119 E:122
patch DISU E:119 E:122
patch DISU F:49 F:98
patch DISU F:119 F:122
patch DISU G:46 G:87
patch DISU G:104 G:107
patch DISU M:45 M:111
patch DISU M:155 M:205
patch DISU N:42 N:110
patch DISU N:164 N:229
patch DISU N:264 M:227
# Very important to add missing H atom coordinates for mutant residues.
guesscoord
writepsf $outMod.psf
writepdb $outMod.pdb
# Check charge to see if ions need alteration.
set tcr3 [mol new $outMod.psf]
mol addfile $outMod.pdb waitfor all
set sel [atomselect $tcr3 all]
measure sumweights $sel weight charge
# ^The charge was zero.
###########################################################
# Clean up
file delete water.pdb; file delete ions.pdb; file delete mem.pdb;
file delete newAchain.pdb; file delete newBchain.pdb; file delete newDchain.pdb;
file delete AchainN.pdb; file delete AchainC.pdb;
file delete BchainN.pdb; file delete BchainX.pdb; file delete BchainC.pdb;
file delete DchainN.pdb; file delete DchainC.pdb;
file delete eChain.pdb; file delete fChain.pdb; file delete gChain.pdb;
file delete mChain.pdb; file delete nChain.pdb; 
file delete s1Chain.pdb; file delete s2Chain.pdb;
file delete fchainN.pdb; file delete fchainX.pdb; file delete fchainC.pdb;
file delete newFchain.pdb;
foreach f [glob water*.pdb] {file delete $f}
# Get zero-based protein residue and index numbers. 
# Change to one-based, Amber numbers for restraining during initial MD.
lindex [[atomselect $tcr3 "name CA and protein"] get residue] end; # max zero-based residue=1068 
lindex [[atomselect $tcr3 protein] get index] end; # max zero-based index=16728
lindex [[atomselect $tcr3 lipid] get residue] end; # max zero-based lidid residue=2670
