# removeExt4.tcl: Remove most extracellular atoms from tcr35c to make tcr38.
# tcr35c.eq70a.rst (322 ns) makes tcr38.[psf, pdb].
# Derived from removeExt3.tcl
# vmd -dispdev none on Eire.
cd /t/tcr38
set outMod tcr38
set tcr [mol new /t/tcr35/tcr35c.psf]
# Use restart file autoimaged on molecule 1
mol addfile /t/tcr35/tcr35c.eq70a.rst type netcdf; # 322 ns
display projection Orthographic
color Display Background white
rotate x by 90
#
set lipids  [atomselect $tcr "chain J"]
llength [lsort -unique [$lipids get resid]]; # 1600 residues in membrane
measure minmax $lipids
# ^ This shows lipids are from z=65.24 to 135.68
# tcr33 has the extracellular domain at z=0-65. Remove down to z=zCap.
# If I remove water with z<53, that would leave at least 12 angstroms above membrane.
set zCap 53; 
###########################################################
# Remove extracellular protein.
# Mark all atoms to save with beta=1.
[atomselect $tcr all] set beta 1
# Remove protein residues up to first transmembrane domain residue.
# Keep all of CD3 zeta residues, chains a and b. They are so small anyway.
[atomselect $tcr "chain A"] writepdb aChain.pdb
[atomselect $tcr "chain B"] writepdb bChain.pdb
### delta chain d TM helix is Asp100 to His128
# Remove Phe22 to Val97, save 2 extracellular residues, Glu98 & Leu99. 
[atomselect $tcr "chain D and resid 22 to 97"] set beta 0
[atomselect $tcr "chain D and beta 1"] writepdb dChain.pdb
### epsilon chain f TM helix Met128 to Arg155
# Remove Gln33 to Cys122
[atomselect $tcr "chain F and resid 33 to 122"] set beta 0
[atomselect $tcr "chain F and beta 1"] writepdb fChain.pdb
### gamma chain g TM helix is Asn111 to Ala137
# Remove Ser24 to Cys107
[atomselect $tcr "chain G and resid 24 to 107"] set beta 0
[atomselect $tcr "chain G and beta 1"] writepdb gChain.pdb
### alpha chain m TM helix is Asp239 to Ser273
# Remove Val26 to Phe236
[atomselect $tcr "chain M and resid 26 to 236"] set beta 0
[atomselect $tcr "chain M and beta 1"] writepdb mChain.pdb
### beta chain n TM helix is Glu269 to Lys308
# Remove Gly22 to Ser268
[atomselect $tcr "chain N and resid 22 to 268"] set beta 0
[atomselect $tcr "chain N and beta 1"] writepdb nChain.pdb
### epsilon chain e TM helix is Asp126 to Asn154
# Remove Gln33 to Glu124
[atomselect $tcr "chain E and resid 33 to 124"] set beta 0
[atomselect $tcr "chain E and beta 1"] writepdb eChain.pdb
# ^They all checked out via inspection.
# Save S1 structural cholesterol and S2 DOPC. 
set S1 [atomselect $tcr "chain S and resname CHL1"] 
$S1 set segid S1; $S1 writepdb S1.pdb
set S2 [atomselect $tcr "chain S and resname DOPC"]
$S2 set segid S2; $S2 writepdb S2.pdb
# Save the whole membrane, 1600 residues.
$lipids writepdb lipids.pdb
###########################################################
# Remove extra water by setting beta to 2 or 3
set rmWater [atomselect $tcr "water and z<$zCap"]
$rmWater set beta 2
[atomselect $tcr "beta 2"] num
# ^ Identifies 173,208 water atoms. 191,633 in tcr15 construction
[atomselect $tcr "same residue as beta 2"] set beta 3
[atomselect $tcr "same residue as beta 3"] num
# ^ Identifies 174,630 water atoms. 193,158  tcr15 construction
# Now save whole waters with beta 1.
# The waters need to be saved in smaller chunks for psfgen.
set saveWat [atomselect $tcr "water and beta 1"]
$saveWat set chain z
# Use "residue" here instead of "resid" because tcr33 has duplicate resid for each segment.
set Res [lsort -unique -integer [$saveWat get residue]]
set last [lrange $Res end end]
# A PDB residue can only have 4-characters
# Save water in 9999 residue chunks
set chunk 9999
for {set n 1; set k 0} {$k<=$last} {incr n} {
  set j [expr ($n-1) * $chunk +1]
  set k [expr $j + $chunk -1]
  puts "$n $j $k"
  set toSave [atomselect $tcr "water and beta 1 and residue $j to $k"]
  # Need to adjust residue# so that they go from 1 to $chunk in each file.
  # That ensures they will be unique within each psfgen segment.
  unset -nocomplain z
  set nRes [expr [$toSave num] / 3]
  for {set i 1} {$i<=$nRes} {incr i} {lappend z $i $i $i} 
  $toSave set resid $z
  $toSave writepdb water$n.pdb
}
# water1.pdb to water[$n-1].pdb saved
# Remove extra extracelllar ions 
set rmIons [atomselect $tcr "ions and z<$zCap"]
$rmIons set beta 4
[atomselect $tcr "beta 4"] num
# ^ removing 364 ions.
[atomselect $tcr "ions and beta 1"] num
# ^ Saving 878 ions
# Adjust charge, initial charge is potentially not zero
[atomselect $tcr "name SOD and beta 1"] num; # 558 SOD
[atomselect $tcr "name CLA and beta 1"] num; # 320 CLA
# Rather than adjust charge of tcr38 made here, worry about
# charge neutrality in models after CD3 CTs added.
[atomselect $tcr "ions and beta 1"] writepdb ions.pdb
###########
# The system above has a void where protein was removed and there is no water.
# Add water into space vacated by protein removal.
# removeExt.tcl built a big water box to fill in removed protein.
set wat [mol new /home/jcannon/tcr2/tcr10/bigWater3.pdb]
[atomselect $wat all] set beta 0
set watSel [atomselect $wat "z>$zCap and z<70"]
# ^ 131,920 water atoms in range. 98320 for tcr15 construction
set rmAtoms [atomselect $tcr "beta 0"]
# ^ 12322 removed atoms, 12383 for tcr15 construction
# Find waters near the removed protein atoms. 
set con [measure contacts 1.0 $rmAtoms $watSel]
llength [lindex $con 1]
# ^That gives 2121 contacts. Set beta 1 for all near water atoms.
set i [lindex $con 1]
[atomselect $wat "index $i"] set beta 1
set addWater [atomselect $wat "beta 1"]
$addWater num
# Now this gives 1708 unique atoms in contact. 1387 for tcr15 construction 
# Water resid is not unique, but residue is.
# Set beta 1 for all near water residues.
set i [$addWater get residue]
[atomselect $wat "residue $i"] set beta 1
set addWater [atomselect $wat "beta 1"]
$addWater num
# Now 3156 atoms, to add. 2586 for tcr15 construction
# Do not add waters that collide with existing waters.
set opWater [atomselect $tcr "water and z<80 and beta 1"]
# Find waters near current waters. 
set con [measure contacts 1.0 $opWater $addWater]
# Contacts with existing waters.
set i [lindex $con 1]
[atomselect $wat "index $i"] set beta 3
set rmWater [atomselect $wat "beta 3"]
set i [$rmWater get residue]
[atomselect $wat "residue $i"] set beta 3
# Potential added waters that contacted exisiting water.
[atomselect $wat "beta 1"] num
# Now only 3051 water atoms will be added. 2499 for tcr15 constructio
[atomselect $wat "beta 1"] writepdb addWater.pdb
###########################################################
# Use psfgen to build from parts
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
segment a {pdb aChain.pdb}; coordpdb aChain.pdb a
segment b {pdb bChain.pdb}; coordpdb bChain.pdb b
segment d {pdb dChain.pdb}; coordpdb dChain.pdb d
segment e {pdb eChain.pdb}; coordpdb eChain.pdb e
segment f {pdb fChain.pdb}; coordpdb fChain.pdb f
segment g {pdb gChain.pdb}; coordpdb gChain.pdb g
segment m {pdb mChain.pdb}; coordpdb mChain.pdb m
segment n {pdb nChain.pdb}; coordpdb nChain.pdb n
# Structural cholesterol CHL1 201, DOPC 190, membrane, ions
segment s1 {pdb S1.pdb}; coordpdb S1.pdb s1
segment s2 {pdb S2.pdb}; coordpdb S2.pdb s2
guesscoord; # Important because of protein truncations.
segment j {pdb lipids.pdb}; coordpdb lipids.pdb j
segment k {pdb ions.pdb}; coordpdb ions.pdb k
# Water is in multiple pdb files
pdbalias residue WAT TIP3
pdbalias atom WAT O OH2
for {set i 1} {$i<$n} {incr i} {
  segment w$i {pdb water$i.pdb}; coordpdb water$i.pdb w$i
  puts $i
}
segment z {pdb addWater.pdb}; coordpdb addWater.pdb z
# Only one disulfide remains from 6jxr.pdb
patch DISU A:32 B:32
# TODO Check on charge and ion balance.
writepsf $outMod.psf
writepdb $outMod.pdb
###########################################################
# clean up
file delete aChain.pdb; file delete bChain.pdb; file delete dChain.pdb;
file delete eChain.pdb; file delete fChain.pdb; file delete gChain.pdb;
file delete mChain.pdb; file delete nChain.pdb; 
file delete S1.pdb; file delete S2.pdb;
file delete addWater.pdb; file delete ions.pdb; file delete lipids.pdb;
foreach f [glob water*.pdb] {file delete $f}
###########################################################
# Check charge to see if ions need alteration.
set tcr3 [mol new $outMod.psf]
mol addfile $outMod.pdb waitfor all
set sel [atomselect $tcr3 all]
measure sumweights $sel weight charge
# ^ shows initial charge is +16
