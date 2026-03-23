# makeNckref.tcl: Make a trajectory of Nck-PRS references
# from 11 edr.*.eq54.rst, 14 fgr.*.eq53.rst frames in NckPRS.dcd, 
# and 2jxb3 frame 9803 in NckPRS2.dcd.
# See notes at end.
# vmd -dispdev none 
cd /t/tcr39an
# Make NckPRS.pdb: 7 PRS CA atoms, chain E and Nck heavy atoms resid 31-86, chain K.
if {[file exists NckPRS.pdb]==0} {
  # Use edm.147, edr.147 as an example
  set psfEx "/t/tcr16/edm/edm.147.psf"
  set pdbEx "/t/tcr16/edr/edr.147.eq54.rst"
  mol load psf $psfEx netcdf $pdbEx
  # Selection of atoms to delete, save only chain K and PRS CA atoms.
  [atomselect top all] set beta 0
  [atomselect top "name CA and chain E and resid 180 to 186"] set beta 1
  [atomselect top "chain K and not hydrogen"] set beta 1
  set del [atomselect top "beta 0"]
  [atomselect top "beta 1"] writepdb NckPRS.pdb
}
###########################################################
# Sorted collection of 11 Nck references from edr.*.eq54.rst frames
# on the basis of frequency of Nck binding in edtPRS3.dat and fgtPRS3.dat.
# In anTcr39d.tcl, all atoms of these models were loaded. 
set refTraj [mol new NckPRS.pdb]
set trajXYZ [atomselect $refTraj all]
set edrModels "147 16 176 154 107 201 20 62 174 75 197"
foreach ref $edrModels {
  set nckRef [mol new /t/tcr16/edm/edm.$ref.psf]
  mol addfile /t/tcr16/edr/edr.$ref.eq54.rst type netcdf waitfor all
  # Nck chain K, resid 31-86. PRS chain E, resid 180-188.
  set PRSref [atomselect $nckRef "name CA and chain E and resid 180 to 186"]
  set nckSH [atomselect $nckRef "chain K and not hydrogen"]
  set allSel [atomselect $nckRef all]
  # Both Nck in chain K and the PRS CA atoms in chain E need to move.
  $allSel set beta 0
  $PRSref set beta 1
  $nckSH set beta 1
  set nckPRSref [atomselect $nckRef "beta 1"]
  # Save nckPRSref as one frame of refTraj
  # Get nckPRSref coordinates and save into trajectory frame.
  set xyz [$nckPRSref get {x y z}]
  $trajXYZ set {x y z} $xyz
  # Add new frame to trajectory (initially duplicate of last).
  animate dup $refTraj
  # Delete atomselections and mol.
  $PRSref delete 
  $nckSH delete
  $nckPRSref delete
  $allSel delete
  mol delete $nckRef
}
# Load collection of 14 Nck references from fgr.*.eq53.rst frames.
set fgrModels "105 177 199 93 10 83 168 142 30 55 154 117 86 201"
foreach ref $fgrModels {
  set nckRef [mol new /t/tcr16/fgm/fgm.$ref.psf]
  mol addfile /t/tcr16/fgr/fgr.$ref.eq53.rst type netcdf waitfor all
  # Nck chain K, resid 31-86. PRS chain F, resid 180-188.
  set PRSref [atomselect $nckRef "name CA and chain F and resid 180 to 186"]
  set nckSH [atomselect $nckRef "chain K and not hydrogen"]
  set allSel [atomselect $nckRef all]
  # Both Nck in chain K and the PRS CA atoms in chain E need to move.
  $allSel set beta 0
  $PRSref set beta 1
  $nckSH set beta 1
  set nckPRSref [atomselect $nckRef "beta 1"]
  # Save nckPRSref as one frame of refTraj
  # Get nckPRSref coordinates and save into trajectory frame.
  set xyz [$nckPRSref get {x y z}]
  $trajXYZ set {x y z} $xyz
  # Add new frame to trajectory (initially duplicate of last).
  animate dup $refTraj
  # Delete atomselections and mol.
  $PRSref delete 
  $nckSH delete
  $nckPRSref delete
  $allSel delete
  mol delete $nckRef
}
molinfo $refTraj get {numframes frame}; # =26, 25
set nFrames [molinfo $refTraj get frame]; # =25
# Save trajectory (frames zero-based), don't repeat the first frame
animate write dcd NckPRS.dcd beg 1 end $nFrames waitfor all $refTraj
# NckPRS.dcd has 7 PRS CA and all heavy Nck atoms, in that order.
###########################################################
# Cannot add first frame from 2JXB because it has a different atom order.
# Add Nck-PRS reference 2jxb3 frame 9803
# Protein-only topology of 2jxb3.top made by clusterNck2.sh 
set nckRef [mol new /home/jcannon/tcr2/cd3e/nowat.top type parm7]
for {set i 2} {$i<=11} {incr i} {
  mol addfile /home/jcannon/tcr2/cd3e/2jxb3.eq$i.cdf type netcdf waitfor all
}
# The 2jxb3 trajectory has 10,000 frames, use frame 9803 for Nck and PRS.
# PRS here is: 6-PRO PRO PRO VAL PRO ASN PRO-12
set PRSref [atomselect $nckRef "name CA and resid 6 to 12" frame 9803]
set nckSH [atomselect $nckRef "resid 31 to 86 and mass>2" frame 9803]
# ^Only use nckSH heavy atoms
# Get nckPRSref coordinates and save into trajectory frame.
set allSel [atomselect $nckRef all]
# Both Nck in chain K and the PRS CA atoms in chain E need to move.
$allSel set beta 0
$PRSref set beta 1
$nckSH set beta 1
set nckPRSref [atomselect $nckRef "beta 1" frame 9803]
set xyz [$nckPRSref get {x y z}]
$trajXYZ set {x y z} $xyz
# Add new frame to trajectory (initially duplicate of last).
animate dup $refTraj
# Delete atomselections and mol.
$PRSref delete 
$nckSH delete
mol delete $nckRef
molinfo $refTraj get {numframes frame}; # =27, 26
set nFrames [molinfo $refTraj get frame]; # =25
# Save trajectory (frames zero-based), don't repeat the first frame
animate write dcd NckPRS2.dcd beg 1 end $nFrames waitfor all $refTraj
###########################################################
# The edm and edr models have chains D,E,K,J,... in that order.
# The fgm and fgr models have chains F,G,K,J,... in that order.
# So even though the DE and FG orders are different, NckPRS.dcd
# has PRS and Nck atom coordinates in the right order because
# chains D and G were not used.
# The output NckPRS.dcd has 7 PRS CA atoms, chain E and Nck heavy atoms resid 31-86, chain K.
# NckPRS.dcd has 26 frames, use frames 0-24 (frames 24 and 25 are the same).

