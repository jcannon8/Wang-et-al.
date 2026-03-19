# makeEdn3.tcl: Make CD3 epsilon-delta dimers with random CTs fused to truncated subunits
# with Nck bound to CD3 epsilon PRS.
# This uses Nck and PRS coordinates in the input cd3edm model.
# Output models, edm.*.[pdb, psf], have protein, membrane, and water.
# This is derived from makeEdn2.tcl.
# vmd -dispdev none on Eire
# Environment variables for model range desired
set first $env(first)
set last $env(last)
cd /t/tcr16/edm
set out [open log.dat w]
# Use psfgen to build edm.$mod.[psf,pdb]
package require psfgen
# Use 6.875 ns equilibrated cd3edm.eq7 autoimaged with chain e for reference. 
# Autoimage using chain E anchor to bring the most water into cytoplasm.
set tcr2 [mol new /t/tcr16/charm-gui8/amber/step5_input.psf ]
mol addfile /t/tcr16/charm-gui8/amber/cd3edm.ref1.pdb waitfor all
# That psf marks chains by segid. E=CD3 epsilon, D=CD3 delta, K=Nck
# chain E: segid PROA, Met600-Ile682; chain D: segid PROB, Glu285-Lys358
# chain K: segid PROC, Glu31-Lys86
set nckAll [atomselect $tcr2 "segid PROC"]
# Save reference Nck coordinates to use before repositioning.
set refKxyz [$nckAll get {x y z}]
# chain E Met125-Ile171 and Glu98-Lys171 in 6JXR resid
# Save TM coordinates of reference, tcr2, so that chains can be returned to original position.
set tmE0 [atomselect $tcr2 "name CA and segid PROA and resid 600 to 629"]
set tmD0 [atomselect $tcr2 "name CA and segid PROB and resid 285 to 314"]
# CT residues: chain e Lys156-Ile207; chain d Thr130-Lys171 (in 6JXR)
set eCT [atomselect $tcr2 "segid PROA and resid 631 to 682"]
set dCT [atomselect $tcr2 "segid PROB and resid 317 to 358"]
# Reference PRS residues
set rPRS [atomselect $tcr2 "name CA and segid PROA and resid 655 to 663"]
set rPRSall [atomselect $tcr2 "segid PROA and resid 655 to 663"]
set all2 [atomselect $tcr2 all]
set limitAll [measure minmax $all2]
# Add back waters that current chain E and D CTs contact beause those chains are removed.
# removeExt.tcl built a big water box to fill in removed protein.
# Use method in removeExt2.tcl
set wat [mol new /home/jcannon/tcr2/tcr10/bigWater3.pdb]
set water0 [atomselect $wat all]
$water0 set beta 0;  # Mark all waters to not add, beta=0
# Add cytoplasmic water that collides with old chainD and E CT and chain K. 
set conD [measure contacts 1.0 $dCT $water0]
set wats [lindex $conD 1]
set selD [atomselect $wat "same residue as index $wats"]
$selD set beta 2
set conE [measure contacts 1.0 $eCT $water0]
set wats [lindex $conE 1]
set selE [atomselect $wat "same residue as index $wats"]
$selE set beta 3
set chainK [atomselect $tcr2 "segname PROC"]
set conK [measure contacts 1.0 $chainK $water0]
set wats [lindex $conK 1]
set selK [atomselect $wat "same residue as index $wats"]
$selK set beta 4
set wat1 [atomselect $wat "beta 2 3 4"]
# Water diagnostics
set initWat [[atomselect $tcr2 "name O and water"] num]
set wat1Num [[atomselect $wat "beta 2 3 4 and name OH2"] num]
puts $out "Initial waters=$initWat, wat1 contact waters=$wat1Num"
# Later, add waters with beta 2 3 4 add to final output if they do not colide with new CTs.
# Save lipids and ions
set lipids [atomselect $tcr2 "segname MEMB"]
$lipids set chain J
$lipids writepdb lipids.pdb
set ions [atomselect $tcr2 "name SOD CLA"]
$ions writepdb ions.pdb
##########
# A second cd3edm copy to change and save.
set tcr [mol new /t/tcr16/charm-gui8/amber/step5_input.psf ]
mol addfile /t/tcr16/charm-gui8/amber/cd3edm.ref1.pdb waitfor all
# Transmembrane residues
set tmE1 [atomselect $tcr "name CA and segid PROA and resid 600 to 629"]
set tmD1 [atomselect $tcr "name CA and segid PROB and resid 285 to 314"]
# CT residues: chain e Lys156-Ile207; chain d Thr130-Lys171
set eCT [atomselect $tcr "segid PROA and resid 631 to 682"]
set dCT [atomselect $tcr "segid PROB and resid 317 to 358"]
set eCT2 [atomselect $tcr "segid PROA and resid 631 to 682 and mass>2"]
set dCT2 [atomselect $tcr "segid PROB and resid 317 to 358 and mass>2"]
# Rename chains
set chainE [atomselect $tcr "segname PROA"]; $chainE set chain E
set chainD [atomselect $tcr "segname PROB"]; $chainD set chain D
set chainK [atomselect $tcr "segname PROC"]; $chainK set chain K
# Renumber resid's for chainE to those in 6JXR.
unset -nocomplain newR
set oldR [$chainE get resid]
foreach r $oldR {lappend newR [expr $r - 475]}
$chainE set resid $newR
# Renumber resid's for chainD to those in 6JXR.
unset -nocomplain newR
set oldR [$chainD get resid]
foreach r $oldR {lappend newR [expr $r - 187]}
$chainD set resid $newR
# chainK is already numbered according to 6JXR.
set ePRS [atomselect $tcr "name CA and chain E and resid 180 to 188"]
set ePRSall [atomselect $tcr "chain E and resid 180 to 188"]
set notPRS [atomselect $tcr "protein and not (chain E and resid 178 to 188)"]
set water [atomselect $tcr water]
set mem [atomselect $tcr "segname MEMB"] 
# atomselections above moved out of loop so they are not created more than once.
#
proc setAng {molA chn ang res} {
# Set backbone angles: molA=molid; chn=segid; ang=HBGLE; res=resid.
# Use HGBLE angles defined in http://eire.umh.edu/mediawiki/index.php/MakeTorsions
# phi, psi H: -57, -47; G: -49, -26; B: -119, 113; L 57, 47; E: -150, 155
  if {$ang == "H"} {
    [atomselect $molA "name CA and segid $chn and resid $res"] set phi -57
    [atomselect $molA "name CA and segid $chn and resid $res"] set psi -47
  }
  if {$ang == "G"} {
    [atomselect $molA "name CA and segid $chn and resid $res"] set phi -49
    [atomselect $molA "name CA and segid $chn and resid $res"] set psi -26
  }
  if {$ang == "B"} {
    [atomselect $molA "name CA and segid $chn and resid $res"] set phi -119
    [atomselect $molA "name CA and segid $chn and resid $res"] set psi 133
  }
  if {$ang == "L"} {
    [atomselect $molA "name CA and segid $chn and resid $res"] set phi 57
    [atomselect $molA "name CA and segid $chn and resid $res"] set psi 47
  }
  if {$ang == "E"} {
    [atomselect $molA "name CA and segid $chn and resid $res"] set phi -150
    [atomselect $molA "name CA and segid $chn and resid $res"] set psi 155
  }
}
proc checkThread {molA resN loopName} {
  # Check loop penetration: 
  # molA=molid; resN=residue numbers; loopName=loop atom names
  set OK 0
  foreach p $resN {
    set others [atomselect $molA "beta 1 and not residue $p"]
    set this [atomselect $molA "residue $p and sidechain"]
    set loop [atomselect $molA "residue $p and name $loopName"]
    set loopCnt [measure center $loop]
    set oops [measure contacts 2.1 $others $this]
    set nCon [llength [lindex $oops 0]]; # nCon=number of contact atom pairs
    if {$nCon>2} {
      foreach i [lsort -unique [lindex $oops 0]] {
        set left [atomselect $molA "index $i"]
        # Get distance of left atom to loopCnt
        set leftPt [$left get {x y z}]
        set d [veclength [vecsub $loopCnt [lindex $leftPt 0]]]               
        if {$d<1.0} {set OK 1}
        $left delete
      }
    }
    $others delete; $this delete; $loop delete
  }  
  return $OK 
}
# Start of loop to make models
for {set mod $first; set t 1} {$mod<=$last} {incr t} {
if {[expr $t % 10]==0} {puts "$t $mod" }
# The makeHGBLE program makes a random HGBLE string using the CT sequence. 
# and the table of frequencies I have used before.
# HGBLEstring=`makeHGBLE eps.seq HGBLEfreq2 ` was used in MakeEpsilon.sh
set strE [exec makeHGBLE /home/jcannon/tcr2/cd3f/eps.seq /home/jcannon/tcr2/cd3f/HGBLEfreq2]
# That ^ CD3 epsilon started at Trp 151. I only need Lys156 to Ile207
# Changing phi and psi in the CT will move the whole subunit.
# Use strE to set new phi, psi in CD3 chain E.
for {set r 156} {$r<=207} {incr r} {
  set j [expr $r-151]
  set ang [string range $strE $j $j]
  if {$r<180 || $r>188} {setAng $tcr "PROA" $ang $r}
}
# Check for chain E self-contacts
set conE [measure contacts 1.0 $eCT2 $eCT2]
set nConE [llength [lindex $conE 0]]
if {$nConE>0} {continue}; # No heavy atom contacts allowed!
# Move chain E back to original location using reference, tcr2, TM CA coordinates.
set M [measure fit $tmE1 $tmE0]
$chainE move $M 
set RMSD [measure rmsd $tmE1 $tmE0]
# ^ The RMSD should be very low
# chainE cannot extend into cytoplasm more than 110.
set limit [measure minmax $chainE]
set deltaZ [expr abs([lindex $limit 0 2] - [lindex $limit 1 2])]
if {$deltaZ>110} {continue}
# Check chain E is within periodic box.
if {[lindex $limit 0 0] < [lindex $limitAll 0 0]} {continue}
if {[lindex $limit 0 1] < [lindex $limitAll 0 1]} {continue}
if {[lindex $limit 1 0] > [lindex $limitAll 1 0]} {continue}
if {[lindex $limit 1 1] > [lindex $limitAll 1 1]} {continue}
# Check chain E collision with membrane
set con [measure contacts 1.0 $eCT $mem]
set flCon [llength [lindex $con 0]]
if {$flCon >0} {continue}
puts "found a chain E"
##########
# Set coordinates of chainK to that of reference
$chainK set {x y z} $refKxyz
# Move Nck to bind to chain E PRS
# measure fit <reference><moveable>
set M [measure fit $rPRS $ePRS]
$chainK move $M
# Compare the Nck to PRS contacts, output below.
set conStart [measure contacts 2.3 $nckAll $rPRSall]
set conFinal [measure contacts 2.3 $chainK $ePRSall]
# Check Nck is within the periodic box 
set limit [measure minmax $chainK]
if {[lindex $limit 0 0] < [lindex $limitAll 0 0]} {continue}
if {[lindex $limit 0 1] < [lindex $limitAll 0 1]} {continue}
if {[lindex $limit 1 0] > [lindex $limitAll 1 0]} {continue}
if {[lindex $limit 1 1] > [lindex $limitAll 1 1]} {continue}
# Check Nck contacts with membrane
set con [measure contacts 1.0 $nckAll $mem]
set flCon [llength [lindex $con 0]]
if {$flCon >0} {continue}
puts "Nck bound"
#########
# Make CD3 delta CT similarly. Try 100 times with a discovered chain E
set notfound 1
for {set k 0} {$k<100 && $notfound} {incr k} {
set strD [exec makeHGBLE /home/jcannon/tcr2/cd3d/delta.seq /home/jcannon/tcr2/cd3f/HGBLEfreq2]
# That ^ CT delta starts at Ala126. Last 6JXR delta residue Glu129.
# I need Thr130 to Lys171 
for {set r 130} {$r<=171} {incr r} {
  set j [expr $r-126]
  set ang [string range $strD $j $j]
  setAng $tcr "PROB" $ang $r
}
# Check for chain D self-contacts
set conD [measure contacts 1.0 $dCT2 $dCT2]
set nConD [llength [lindex $conD 0]]
if {$nConD>0} {incr t; continue}; # No heavy atom contacts allowed!
# Move chain D back to original location using reference, tcr2, TM CA coordinates.
set M [measure fit $tmD1 $tmD0]
$chainD move $M 
#set RMSD [measure rmsd $tmD1 $tmD0]
# ^ The RMSD should be very low
# chainD cannot extend into cytoplasm more than 110.
set limit [measure minmax $chainD]
set deltaZ [expr abs([lindex $limit 0 2] - [lindex $limit 1 2])]
if {$deltaZ>110} {incr t; continue}
# Check chain D is within periodic box.
if {[lindex $limit 0 0] < [lindex $limitAll 0 0]} {continue}
if {[lindex $limit 0 1] < [lindex $limitAll 0 1]} {continue}
if {[lindex $limit 1 0] > [lindex $limitAll 1 0]} {continue}
if {[lindex $limit 1 1] > [lindex $limitAll 1 1]} {continue}
# Check chain D collision with membrane
set con [measure contacts 1.0 $dCT $mem]
set dlCon [llength [lindex $con 0]]
if {$dlCon >0} {incr t; continue}
puts "found a chain D"
set notfound 0
}
# End of k loop for chain D
if {$notfound == 1} {puts "try another chain E"; continue}
#
# Check collison of eCT and dCT
set con [measure contacts 1.0 $eCT $dCT]
set edCon [llength [lindex $con 0]]
if {$edCon >0} {incr t; continue}
# Check if Nck binds to anything besides the PRS in chain E
set con [measure contacts 1.0 $chainK $notPRS]
set flCon [llength [lindex $con 0]]
if {$flCon >0} {continue}
#######################################
# Ready to create and save the model.
# Save protein chains separately for psfgen, remove water new CTs and Nck contact,
# and add back water old CTs contacted.
# The resid's were altered above to match 6JXR and 2JXB before saving.
$chainD writepdb chainD.pdb
$chainE writepdb chainE.pdb
$chainK writepdb chainK.pdb
# Remove cytoplasmic water that collides with new chain D E CT, and chain K.
# Set beta=1 for water to keep.
$water set beta 1
set conD [measure contacts 1.0 $dCT $water]
set wats [lindex $conD 1]
if {[llen $wats] >0} {
  set selD [atomselect $tcr "same residue as index $wats"]
  $selD set beta 8
  $selD delete
}
set conE [measure contacts 1.0 $eCT $water]
set wats [lindex $conE 1]
if {[llen $wats] >0} {
  set selE [atomselect $tcr "same residue as index $wats"]
  $selE set beta 7
  $selE delete
  }  
set conK [measure contacts 1.0 $chainK $water]
set wats [lindex $conK 1]
if {[llen $wats] >0} {
  set selK [atomselect $tcr "same residue as index $wats"]
  $selK set beta 9
  $selK delete
}
# Water diagnostics
set subWat [atomselect $tcr "beta 7 8 9 and name O OH2"]
set conNew [$subWat num]
puts $out "For model edn.$mod"
puts $out "Subtract waters contacting new chain ED and K = $conNew"
$subWat delete
# Now save all waters in mol tcr with beta=1
# The waters need to be saved in smaller chunks for psfgen.
set saveWat [atomselect $tcr "water and beta 1"]
$saveWat set chain z
set Res [lsort -unique -integer [$saveWat get residue]]
set lastWat [lrange $Res end end]
# A PDB residue can only have 4-characters
# Save water in 9999 residue chunks
set chunk 9999
for {set n 1; set k 0} {$k<=$lastWat} {incr n} {
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
  $toSave delete
}
$saveWat delete
# water1.pdb to water[$n-1].pdb saved above
# Save waters added from removed old CTs, wat1, which are not in contact with new CTs.
set conD2 [measure contacts 1.0 $dCT $wat1]
set conE2 [measure contacts 1.0 $eCT $wat1]
set conK2 [measure contacts 1.0 $chainK $wat1]
# Indices of wat2 waters to not add
unset -nocomplain notWat2
append notWat2 [lindex $conE2 1] [lindex $conD2 1] [lindex $conK2 1]
# This addWat is different from the above used in makeEdn.tcl
set addWat [atomselect $wat "same residue as beta 2 3 4 and not index $notWat2"]
# Water diagnostics
set addWatNum [expr [$addWat num] / 3]
puts $out "Add waters from initial CTs = $addWatNum"
# Compare the Nck to PRS contacts
set nStart [llength [lindex $conStart 1]]
set nFin [llength [lindex $conFinal 1]]
puts $out "Nck-PRS: start=$nStart, final=$nFin"
flush $out
# Now save addWat
$addWat writepdb waterExt.pdb 
$addWat delete
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
segment d {pdb chainD.pdb}; coordpdb chainD.pdb d
segment e {pdb chainE.pdb}; coordpdb chainE.pdb e
segment k {pdb chainK.pdb}; coordpdb chainK.pdb k
segment j {pdb lipids.pdb}; coordpdb lipids.pdb j
segment q {pdb ions.pdb}; coordpdb ions.pdb q
# Load in water
pdbalias residue WAT TIP3
pdbalias atom WAT O OH2
# Water is in multiple pdb files
for {set i 1} {$i<$n} {incr i} {
  segment w$i {pdb water$i.pdb}; coordpdb water$i.pdb w$i
  puts $i
}
# Add extra water
incr i
segment w$i {pdb waterExt.pdb}; coordpdb waterExt.pdb w$i
# Only one disulfide remains from 6jxr.pdb, between zeta chains, not here.
# Hydrogens are in the cd3edm input for all protein, so guesscoord is unnecessary.
guesscoord
writepsf edm.$mod.psf
writepdb edm.$mod.pdb
# Important to reset the psf for next model.
resetpsf
#######################################
# Check for loop penetration
set tcr3 [mol new edm.$mod.psf]
mol addfile edm.$mod.pdb waitfor all
# Set the CD3 CTs beta 1 to screen for threading
set all3 [atomselect $tcr3 all]; $all3 set beta 0
set eCT3 [atomselect $tcr3 "chain E and resid 156 to 207"]
set dCT3 [atomselect $tcr3 "chain D and resid 130 to 171"]
set nck [atomselect $tcr3 "chain K"]
$eCT3 set beta 1; $dCT3 set beta 1; $nck set beta 1
$all3 delete; $eCT3 delete; $dCT3 delete; $nck delete
# Get residue numbers. These are global residue numbers (Amber numbering), 
# not like resid, which are chain-specific.
set ProRes [atomselect $tcr3 "name CA and resname PRO and beta 1"]
set pro [$ProRes get residue]; $ProRes delete
# No Phe in CD3 CTs or in Nck SH3.1
set TyrRes [atomselect $tcr3 "name CA and resname TYR and beta 1"] 
set tyr [$TyrRes get residue]; $TyrRes delete
set HisRes [atomselect $tcr3 "name CA and resname HSD and beta 1"]
set his [$HisRes get residue]; $HisRes delete
set TrpRes [atomselect $tcr3 "name CA and resname TRP and beta 1"]
set trp [$TrpRes get residue]; $TrpRes delete
# Start checking
if {[checkThread $tcr3 $pro "N CA CB CD CG"]} {mol delete $tcr3; continue}
if {[checkThread $tcr3 $tyr "CZ CE1 CD1 CG CD2 CE2"]} {mol delete $tcr3; continue}
if {[checkThread $tcr3 $his "CG ND1 CE1 NE2 CD2"]} {mol delete $tcr3; continue}
if {[checkThread $tcr3 $trp "CG CD1 NE1 CE2 CD2"]} {mol delete $tcr3; continue}
if {[checkThread $tcr3 $trp "CE2 CD2 CE3 CZ3 CH2 CZ2"]} {mol delete $tcr3; continue}
puts "edn.$mod.pdb passed the loop penetration test"
mol delete $tcr3
incr mod
}
# End of model construction loop
###########################################################
# Clean up
file delete waterExt.pdb; file delete ions.pdb; file delete lipids.pdb;
file delete chainD.pdb; file delete chainE.pdb; file delete chainK2.pdb
foreach f [glob water*.pdb] {file delete $f}
# I used 1.0 angstrom as the distance to consider atoms colliding. This might
# be too short because Van der Waals radii of H 1.0-1.4, C 1.9, O 1.66, N 1.8
# are much larger and contacts would occur at the sum of VdW radii. Nevertheless,
# contacts are discovered in equilibrated systems with the 1.0 angstrom distance.
close $out
