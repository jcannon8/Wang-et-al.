# thread2.tcl: Detect threading through Pro, Tyr, His, and Trp sidechains in edt models.
# Use vmd -dispdev none
# Screening models made by makeCDed3.tcl, which are not minimized will
# have many false positives because atoms from neighboring residues are near.
# That is if the thread1.tcl scheme is used that merely tests distances to 
# side chains.
set out [open thread2.dat w]
for {set mod 100} {$mod<=200} {incr mod} {
#set mod 29
set tcr [mol new /home/jcannon/tcr16/edt/edt.$mod.psf]
mol addfile /home/jcannon/tcr16/edt/edt.$mod.pdb waitfor all
puts $out "testing edt.$mod.pdb"
# Set the CD3 CTs beta 1 to screen for threading
set all [atomselect $tcr all]; $all set beta 0
set dCT [atomselect $tcr "chain D and resid 130 to 171"]; $dCT set beta 1
set eCT [atomselect $tcr "chain E and resid 156 to 207"]; $eCT set beta 1
$all delete; $dCT delete; $eCT delete
# Get residue numbers. These are global residue numbers (Amber numbering), 
# not like resid, which are chain-specific.
set ProRes [atomselect $tcr "name CA and resname PRO and beta 1"]
set pro [$ProRes get residue]; $ProRes delete
# No Phe!
set TyrRes [atomselect $tcr "name CA and resname TYR and beta 1"] 
set tyr [$TyrRes get residue]; $TyrRes delete
set HisRes [atomselect $tcr "name CA and resname HSD and beta 1"]
set his [$HisRes get residue]; $HisRes delete
set TrpRes [atomselect $tcr "name CA and resname TRP and beta 1"]
set trp [$TrpRes get residue]; $TrpRes delete
#
foreach p $pro {
  set others [atomselect $tcr "beta 1 and not residue $p"]
  set this [atomselect $tcr "residue $p and sidechain"]
  # Loop heavy atoms for Pro: name N CA CB CD CG
  set loop [atomselect $tcr "residue $p and name N CA CB CD CG"]
  set loopCnt [measure center $loop]
  set oops [measure contacts 2.1 $others $this]
  set nCon [llength [lindex $oops 0]]; # nCon=number of contact atom pairs
  if {$nCon>2} {
    set right [atomselect $tcr "name CA and residue $p"] 
    set rightName [$right get {resname resid chain}]
    foreach i [lsort -unique [lindex $oops 0]] {
      set left [atomselect $tcr "index $i"]
      set leftName [$left get {resname resid chain name}]
      # Get distance of left atom to loopCnt
      set leftPt [$left get {x y z}]
      set d [veclength [vecsub $loopCnt [lindex $leftPt 0]]]               
      if {$d<1.0} {puts $out "Warning: threaded Pro $rightName-$leftName $d"}
      $left delete
    }
    $right delete
  }
  $others delete; $this delete; $loop delete
}
foreach p $tyr {
  set others [atomselect $tcr "beta 1 and not residue $p"]
  set this [atomselect $tcr "residue $p and sidechain"]
  # Loop heavy atoms for Tyr: CZ CE1 CD1 CG CD2 CE2
  set loop [atomselect $tcr "residue $p and name CZ CE1 CD1 CG CD2 CE2"]
  set loopCnt [measure center $loop]
  set oops [measure contacts 2.1 $others $this]
  set nCon [llength [lindex $oops 0]]; # nCon=number of contact atom pairs
  if {$nCon>2} {
    set right [atomselect $tcr "name CA and residue $p"] 
    set rightName [$right get {resname resid chain}]
    foreach i [lsort -unique [lindex $oops 0]] {
      set left [atomselect $tcr "index $i"]
      set leftName [$left get {resname resid chain name}]
      # Get distance of left atom to loopCnt
      set leftPt [$left get {x y z}]
      set d [veclength [vecsub $loopCnt [lindex $leftPt 0]]]               
      if {$d<1.0} {puts $out "Warning: threaded Tyr $rightName-$leftName $d"}
      $left delete
    }
    $right delete
  }
  $others delete; $this delete; $loop delete
}
foreach p $his {
  set others [atomselect $tcr "beta 1 and not residue $p"]
  set this [atomselect $tcr "residue $p and sidechain"]
  # Loop heavy atoms for HIS: CG ND1 CE1 NE2 CD2
  set loop [atomselect $tcr "residue $p and name CG ND1 CE1 NE2 CD2"]
  set loopCnt [measure center $loop]
  set oops [measure contacts 2.1 $others $this]
  set nCon [llength [lindex $oops 0]]; # nCon=number of contact atom pairs
  if {$nCon>2} {
    set right [atomselect $tcr "name CA and residue $p"] 
    set rightName [$right get {resname resid chain}]
    foreach i [lsort -unique [lindex $oops 0]] {
      set left [atomselect $tcr "index $i"]
      set leftName [$left get {resname resid chain name}]
      # Get distance of left atom to loopCnt
      set leftPt [$left get {x y z}]
      set d [veclength [vecsub $loopCnt [lindex $leftPt 0]]]               
      if {$d<1.0} {puts $out "Warning: threaded His $rightName-$leftName $d"}
      $left delete
    }
    $right delete
  }
  $others delete; $this delete; $loop delete
}
foreach p $trp {
  set others [atomselect $tcr "beta 1 and not residue $p"]
  set this [atomselect $tcr "residue $p and sidechain"]
  # Loop heavy atoms for Trp: CG CD1 NE1 CE2 CD2 and  CE2 CD2 CE3 CZ3 CH2 CZ2  
  set loop [atomselect $tcr "residue $p and name CG CD1 NE1 CE2 CD2"]
  set loopCnt [measure center $loop]
  set oops [measure contacts 2.1 $others $this]
  set nCon [llength [lindex $oops 0]]; # nCon=number of contact atom pairs
  if {$nCon>2} {
    set right [atomselect $tcr "name CA and residue $p"] 
    set rightName [$right get {resname resid chain}]
    foreach i [lsort -unique [lindex $oops 0]] {
      set left [atomselect $tcr "index $i"]
      set leftName [$left get {resname resid chain name}]
      # Get distance of left atom to loopCnt
      set leftPt [$left get {x y z}]
      set d [veclength [vecsub $loopCnt [lindex $leftPt 0]]]               
      if {$d<1.0} {puts $out "Warning: threaded Trp $rightName-$leftName $d"}
      $left delete
    }
    $right delete
  }
  $others delete; $this delete; $loop delete
}
foreach p $trp {
  set others [atomselect $tcr "beta 1 and not residue $p"]
  set this [atomselect $tcr "residue $p and sidechain"]
  # Loop heavy atoms for Trp: CG CD1 NE1 CE2 CD2 and  CE2 CD2 CE3 CZ3 CH2 CZ2  
  set loop [atomselect $tcr "residue $p and name CE2 CD2 CE3 CZ3 CH2 CZ2"]
  set loopCnt [measure center $loop]
  set oops [measure contacts 2.1 $others $this]
  set nCon [llength [lindex $oops 0]]; # nCon=number of contact atom pairs
  if {$nCon>2} {
    set right [atomselect $tcr "name CA and residue $p"] 
    set rightName [$right get {resname resid chain}]
    foreach i [lsort -unique [lindex $oops 0]] {
      set left [atomselect $tcr "index $i"]
      set leftName [$left get {resname resid chain name}]
      # Get distance of left atom to loopCnt
      set leftPt [$left get {x y z}]
      set d [veclength [vecsub $loopCnt [lindex $leftPt 0]]]               
      if {$d<1.0} {puts $out "Warning: threaded Trp $rightName-$leftName $d"}
      $left delete
    }
    $right delete
  }
  $others delete; $this delete; $loop delete
}
# Now test heavy atom contacts within a chain
set eCT2 [atomselect $tcr "chain E and resid 156 to 207 and mass>2"]
set dCT2 [atomselect $tcr "chain D and resid 130 to 171 and mass>2"]
set conE [measure contacts 1.0 $eCT2 $eCT2]
set conD [measure contacts 1.0 $dCT2 $dCT2]
set nConE [llength [lindex $conE 0]]
set nConD [llength [lindex $conD 0]]
$eCT2 delete; $dCT2 delete
puts $out "$nConD chain D contacts, $nConE chain E contacts"
mol delete $tcr
}
close $out
###########################################################
