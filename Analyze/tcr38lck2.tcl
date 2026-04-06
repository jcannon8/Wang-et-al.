# tcr38lck2.tcl: Report Lck binding data for tcr38.*.** full trajectories.
# Derived from tcr38lck1.tcl
# Collisions are heavy atoms less than or equal to 1.0 angtroms apart.
# Contacts are any atoms less than or equal to 2.0 angtroms apart.
# Call from an external script RunanTcr39r.sh using
# export mod=0.00
# vmd -dispdev none -e tcr38lck2.tcl
set mod $env(mod);  # tcr38.[0-11].[00-22] model number
set name tcr38.$mod  
cd /t/tcr38an
set rmdsMax 3.69;  # RMSD limit for Lck substrate fit.
set out [open temp w]; # Save Lck bindings that pass all criteria
# Concatenated into tcr38lck2.dat by Runtcr38lck2.sh
#
proc memZ {tcr} {
# Average cytoplasmic side P Z-dimension.
set P [atomselect $tcr "name P"]
set span [measure minmax $P]
# Cytoplasm has greater Z dimension
set MemTop [lindex $span 0 2]
set MemBottom [lindex $span 1 2]
set ave [expr ($MemBottom + $MemTop)/2]
# The cytoplasmic side P atoms have z>ave.
set cP [atomselect $tcr "name P and z>$ave"]
# Now get the average of cytoplasmic P
set spanP [measure minmax $cP]
set minP [lindex $spanP 0 2]
set maxP [lindex $spanP 1 2]
set aveP [expr ($minP + $maxP)/2]
$P delete; $cP delete
return $aveP
}
###########################################################
# ITAM list with {chain resid}
set itams {{A 72} {A 83} {A 111} {A 123} {A 142} {A 153} 
{B 72} {B 83} {B 111} {B 123} {B 142} {B 153}  
{D 149} {D 160} {E 188} {E 199} {F 188} {F 199} {G 160} {G 171} }
# Setup Lck-ITAM atomselections.
set lckpep [mol new /home/jcannon/lck/Lck-Peptide.pdb]
set lck [atomselect $lckpep "chain A"]
set lcksub [atomselect $lckpep all]
set lckH [atomselect $lckpep "chain A and mass>2"]; # Lck heavy atoms
# Choose residues for fitting
set t 4; # The 9 residues around pTyr for fitting.
set c [expr 285- $t]; # first Lck substrate residue to fit
set d [expr 285+ $t]; # last Lck substrate residue to fit
set pSel [atomselect $lckpep "name CA and chain B and resid $c to $d"]
# Start testing model
set tcr [mol new /t/tcr38b/$name.psf]
# Report for eq4-eq85
for {set i 4} {$i<=85} {incr i} {
mol addfile /t/tcr38b/$name.eq$i.rst type netcdf waitfor all
#
set mem [atomselect $tcr "chain J and mass>2"]
set memC [memZ $tcr]
# Start testing each ITAM
foreach s $itams {
set chn [lindex $s 0]
set n [lindex $s 1]
set a [expr $n- $t]; # first ITAM residue to fit
set b [expr $n+ $t]; # last ITAM residue to fit
set iSel [atomselect $tcr "chain $chn and name CA and resid $a to $b"]
# Fitting: fit substrate to ITAM and move lcksub
set M [measure fit $pSel $iSel]
$lcksub move $M
set rms [measure rmsd $pSel $iSel]
$iSel delete
# Use third quartile of RMSD to disqualify Lck binding.
if {$rms > $rmdsMax} {continue}
# Now check if ITAM containing subunit collides with Lck
# beyond the ITAM residues in the active site.
set w 7; # How far away from pTyr are Lck collisions allowed? 
# With w=7, 15-residues of substrate in active site.
set e [expr $n- $w];set f [expr $n+ $w]
set iNot [atomselect $tcr "chain $chn and mass>2 and not (resid $e to $f)"]
set jNot [atomselect $tcr "protein and mass>2 and not chain $chn"]
set conSub [measure contacts 1.0 $lckH $iNot]
# Check Lck collision with ITAM containing subunit (proximal collision).
if {[lsearch $conSub {}] ==-1} {$iNot delete;$jNot delete;continue}  
# Check Lck collision with membrane.
set conMem [measure contacts 1.0 $lckH $mem]
if {[lsearch $conMem {}] ==-1} {$iNot delete;$jNot delete;continue}
# Check Lck collision with CD3 CTs
set conCT [measure contacts 1.0 $lckH $jNot]
if {[lsearch $conCT {}] ==-1} {$iNot delete;$jNot delete;continue}
# No collisions
# Lck to mean membrane distance
set minLck [lindex [measure minmax $lck] 0 2]
set meanMem [expr $minLck - $memC]
#
# Get shortest Lck distance to membrane.
# Increment contacts cutoff, tt, until atom pairs are found.
set tt 10;set LckNum 0 
while {$LckNum==0} {
  set LckMem [measure contacts $tt $lckH $mem]
  set LckNum [llen [lindex $LckMem 1]]
  incr tt 2
}        
unset -nocomplain dl     
#
foreach a2 [lindex $LckMem 0] b2 [lindex $LckMem 1] {
  set ax [atomselect $lckpep "index $a2"]
  set bx [atomselect $tcr "index $b2" frame $i]
  set aa [$ax get {x y z}]
  set bb [$bx get {x y z}]
  set d [vecdist [lindex $aa 0] [lindex $bb 0]]; # distance between atom pair
  lappend dl $d
  set minMem [lindex [lsort $dl] 0]; # Shortest Lck distance to membrane           
  $ax delete; $bx delete
}
# Output for ITAMs that can bind Lck with RMSD<3.67 and 
# no collisions with CT protein or membrane.
# Example of tcr38lck2.tcl output:
#   tcr38.0.00   8 A  153  3.12    9.57    9.45
# Output: <model><MD frame eq><chain><resid><RMSD><to mem mean><to mem>
# Post-processing can count Lck binding in each MD frame by screening column 3.
# The mean membrane distance can be collected from columns 8 and 9.
puts $out [format "%12s %3d %s %4d %5.2f %7.2f %7.2f" \
$name $i $chn $n $rms $meanMem $minMem] 
# Done with this ITAM, delete atomelect's
$iNot delete; $jNot delete
}; # next ITAM
$mem delete
}; # next MD frame
quit
