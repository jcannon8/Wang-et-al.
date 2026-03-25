# makeTcr35.tcl: Make TCRβ L281A/L285A mutant from PDB 7FJE
# with added DOPC 190 from fgr185.  
# It has one extracellular side structural cholesterol and DOPC in S2.
# Derived from makeTcr34.tcl:
# chain A add three C-terminal residues
# chain B add residues to N- and C-termini
# chain D add one C-termial residue
# chain E is missing E:119 E:122 disulfide.
# chain F is missing Glu70 Asp71 Asp72 Lys73 from chain F.
# ^No, Those missing residues added by CHARRM_GUI.
# vmd -dispdev none on Eire
cd /t/tcr35
# There are two potential DPPC240 orientations for tcr35a and tcr35b.
# There is one DOPC190 for tcr35c
set outMod tcr35c
# Start with CHARMM-GUI output from 7FJE.
set tcr [mol new /t/tcr34/charmm-gui-7295924549/amber/step5_input.psf]
mol addfile /t/tcr34/charmm-gui-7295924549/amber/step5_input.pdb waitfor all
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
#############################
# Use tcr12 from 6JXR to extend some chains.
set tcr12 [mol new /home/jcannon/charm/charmm-gui-2876761245/amber/tcr12.psf]
mol addfile /home/jcannon/charm/charmm-gui-2876761245/amber/tcr12.pdb
# Extend chain A by three residues from 22-54 to 22-57
[atomselect $tcr12 "segname PROA"] set chain A
# Use five CA to position
set A30 [atomselect $tcr "name CA and chain A and resid 50 to 54"]
set A12 [atomselect $tcr12 "name CA and chain A and resid 50 to 54"] 
# Move tcr12 chain A to fit tcr30 
set M [measure fit $A12 $A30]
set A12b [atomselect $tcr12 "chain A"]
$A12b move $M
measure rmsd $A12 $A30; # =0.2912 
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
# Use all of 7FJE CA to position
set B30 [atomselect $tcr "name CA and chain B and resid 26 to 54"]
set B12 [atomselect $tcr12 "name CA and chain B and resid 26 to 54"] 
# Move tcr12 chain B to fit tcr33 
set M [measure fit $B12 $B30]
set B12b [atomselect $tcr12 "chain B"]
$B12b move $M
measure rmsd $B12 $B30; # =2.1315
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
measure rmsd $D12 $D30; # =0.701
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
[atomselect $tcr "chain F"] writepdb fChain.pdb
[atomselect $tcr "chain G"] writepdb gChain.pdb
[atomselect $tcr "chain M"] writepdb mChain.pdb
[atomselect $tcr "chain N"] writepdb nChain.pdb
# Save one structural cholesterol separately
set chl [atomselect $tcr "segid HETA"]
$chl set chain S
$chl writepdb s1Chain.pdb
#set unk [mol new /t/mn34/DPPC240b.pdb]
set unk [mol new /t/mn34/DOPC190a.pdb]
set unkSel [atomselect $unk all]
$unkSel set chain S
$unkSel writepdb DOPC190a.pdb
[atomselect $tcr "chain J"] writepdb mem.pdb
# Ion chain names already assigned
[atomselect $tcr "chain I"] writepdb ions.pdb
# No need to alter ions because CHARMM-GUI got it right.
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
segment f {pdb fChain.pdb}; coordpdb fChain.pdb f
segment g {pdb gChain.pdb}; coordpdb gChain.pdb g
segment m {pdb mChain.pdb}; coordpdb mChain.pdb m
segment n {pdb nChain.pdb}; coordpdb nChain.pdb n
# One structural cholesterols, DPPC240, membrane, ions
segment s1 {pdb s1Chain.pdb}; coordpdb s1Chain.pdb s1
segment s2 {pdb DOPC190a.pdb}; coordpdb DOPC190a.pdb s2
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
# 14 disulfides in TCR counting the missing E:119 E:122 from 7fje.pdb
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
# ^The charge was zero, the mutations did not change charge.
###########################################################
# Clean up
file delete water.pdb; file delete ions.pdb; file delete mem.pdb;
file delete newAchain.pdb; file delete newBchain.pdb; file delete newDchain.pdb;
file delete AchainN.pdb; file delete AchainC.pdb;
file delete BchainN.pdb; file delete BchainX.pdb; file delete BchainC.pdb;
file delete DchainN.pdb; file delete DchainC.pdb;
file delete eChain.pdb; file delete fChain.pdb; file delete gChain.pdb;
file delete mChain.pdb; file delete nChain.pdb; 
file delete s1Chain.pdb;
foreach f [glob water*.pdb] {file delete $f}
# Get zero-based protein residue and index numbers. 
# Change to one-based, Amber numbers for restraining during initial MD.
lindex [[atomselect $tcr3 "name CA and protein"] get residue] end; # max zero-based residue=1068 
lindex [[atomselect $tcr3 protein] get index] end; # max zero-based index=16710
lindex [[atomselect $tcr3 lipid] get residue] end; # max zero-based lidid residue=2670

