# anTcr38r.tcl: Nck analysis of tcr38b models every 5 ns from *.rst files. 
# This uses multiple references for PRS-bound Nck from [edr,fgr].*.eq54.rst frames.
# Those are in NckPRS.dcd made by makeNckref.tcl
# Collisons are heavy atoms less than or equal to 1.0 angstroms.
# Output: $inMod.Nckt3.dat
# vmd -dispdev none
# Call from an external script RunanTcr38r.sh using
# export mod=0.00
# vmd -dispdev none -e anTcr38r.tcl
set mod $env(mod);  # tcr38.[0-11].[00-22] model number
set name tcr38.$mod 
cd /t/tcr38an
# Using wrapped angles (Hollingsworth and Karplus 2010).
proc wrapPhi {a} {
  return [expr $a<0?$a+360:$a]
}
proc wrapPsi {a} {
  return [expr $a<-120?$a+360:$a]
}
# Average cytoplasmic side P Z-dimension.
proc memZ2 {tcr} {
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
$P delete
return $aveP
}
#
set dist 1.0;  # Nck collision distance for heavy atoms
# Load the Nck-PRS references made by makeNckref.tcl
set nckRef [mol new /t/tcr39an/NckPRS.pdb]
mol addfile /t/tcr39an/NckPRS.dcd waitfor all
# Define reference atomselections (use only heavy atoms)
set PRSref [atomselect $nckRef "name CA and chain E and resid 180 to 186"]
set nckSH [atomselect $nckRef "chain K and mass>2"]
set refAll [atomselect $nckRef "mass>2"]
set numRef [expr [molinfo $nckRef get numframes]-1]
###########################################################
# For all tcr38b models, load *.eq[4-85].rst
set out [open $name.Nckt3.dat w]
set tcr [mol new /t/tcr38b/$name.psf]
# Report for eq4-eq85
for {set i 4} {$i<=85} {incr i} {
  mol addfile /t/tcr38b/$name.eq$i.rst type netcdf waitfor all
}
set lastframe [molinfo $tcr get numframes]
# After loading all 82 rst files, only 5% Eire memory is used.
###########################################################
# PRS is Pro180-Pro186 in 6JXR chains E and F (see anEdt2.sh).
set ePRS [atomselect $tcr "name CA and chain E and resid 180 to 186"]
set fPRS [atomselect $tcr "name CA and chain F and resid 180 to 186"]
# Set beta=1 for heavy cytoplasmic atoms
[atomselect $tcr all] set beta 0
[atomselect $tcr "chain A and mass>2 and resid>57"] set beta 1
[atomselect $tcr "chain B and mass>2 and resid>57"] set beta 1
[atomselect $tcr "chain D and mass>2 and resid>128"] set beta 1
[atomselect $tcr "chain F and mass>2 and resid>155"] set beta 1
[atomselect $tcr "chain G and mass>2 and resid>137"] set beta 1
[atomselect $tcr "chain E and mass>2 and resid>154"] set beta 1
# notPRS atoms are more than two residues away from PRS
set eNotPRS [atomselect $tcr "beta 1 and not (chain E and resid 178 to 188)"]
set fNotPRS [atomselect $tcr "beta 1 and not (chain F and resid 178 to 188)"] 
set mem [atomselect $tcr "chain J and mass>2"]
# Look at critical PRS angles P182psi and P184psi.
set P182e [atomselect $tcr "chain E and name CA and resid 182"]
set P184e [atomselect $tcr "chain E and name CA and resid 184"]
set P182f [atomselect $tcr "chain F and name CA and resid 182"]
set P184f [atomselect $tcr "chain F and name CA and resid 184"]
# ^ The 16 atomselections above should be deleted after each tcr model.
# Analyze the rst MD frames.S
for {set f 0} {$f<$lastframe} {incr f} {
  $ePRS frame $f; $fPRS frame $f;
  $eNotPRS frame $f; $fNotPRS frame $f
  $mem frame $f
  $P182e frame $f;$P184e frame $f
  $P182f frame $f;$P184f frame $f
  # Z-axis points towards cytoplasm.
  set memZ [memZ2 $tcr]; # Average of cytoplasmic leaflet phosphates.
  # Check contacts for Nck-PRS references
  set nRefE -1; set nRefF -1;  # Output first Nck-PRS reference with no contacts. 
  set eNckMem 99; set fNckMem 99;  # Default Nck to membrane distance
  set eNckMem2 99; set fNckMem2 99;  # Default Nck to membrane distance
  set t 10;set NckNum 0
  # For each Nck reference, check ePRS binding
  for {set ref 0} {$ref<$numRef} {incr ref} {
    $PRSref frame $ref; $nckSH frame $ref; $refAll frame $ref;
    # Move refAll to bind to ePRS
    set M [measure fit $PRSref $ePRS]; # Fit the 7 PRS CA atoms
    $refAll move $M
    # Check protein contacts with Nck SH3.1
    set nckCon [measure contacts $dist $nckSH $eNotPRS] 
    set conNum [llen [lindex $nckCon 1]]
    if {$conNum==0} {
      set nRefE $ref; # Output Nck-PRS reference# with no protein contacts
      # Get Nck to mean membrane distance
      set loc [measure minmax $nckSH]
      set minZ [lindex $loc 0 2]
      set eNckMem [expr $minZ- $memZ]; # minimum distance to mean membrane.
      if {$eNckMem<-6} {break}; # Previous criterion of Nck membrane collision.
      # Get shortest Nck distance to membrane.
      # Increment contacts cutoff, t, until atom pairs are found. 
      while {$NckNum==0} {
        set NckMem [measure contacts $t $nckSH $mem]
        set NckNum [llen [lindex $NckMem 1]]
        incr t 2
      }        
      unset -nocomplain dl     
      foreach a [lindex $NckMem 0] b [lindex $NckMem 1] {
        set ax [atomselect $nckRef "index $a" frame $ref]
        set bx [atomselect $tcr "index $b" frame $f]
        set aa [$ax get {x y z}]
        set bb [$bx get {x y z}]
        set d [vecdist [lindex $aa 0] [lindex $bb 0]]; # distance between atom pair
        lappend dl $d
        set eNckMem2 [lindex [lsort $dl] 0]; # Shortest Nck distance to membrane           
        $ax delete; $bx delete
      }
      break; # Stop searching once the first Nck binding without contact found.
    }
  }
  puts "Done with eNck, t=$t"
  # If no Nck-PRS binds without protein collision, membrane distance or collision
  # was not tested.
  set t 10;set NckNum 0
  for {set ref 0} {$ref<$numRef} {incr ref} {
    $PRSref frame $ref; $nckSH frame $ref; $refAll frame $ref;
    # Move refNck to bind to fPRS
    set M [measure fit $PRSref $fPRS]; # Fit the 7 PRS CA atoms
    $refAll move $M  
    # Check protein contacts with Nck SH3.1
    set nckCon [measure contacts $dist $nckSH $fNotPRS]
    set conNum [llen [lindex $nckCon 1]]
    if {$conNum==0} {
      set nRefF $ref; # Output Nck-PRS reference# with no protein contacts
      # Get Nck to mean membrane distance
      set loc [measure minmax $nckSH]
      set minZ [lindex $loc 0 2]
      set fNckMem [expr $minZ- $memZ]; # minimum distance to mean membrane.
      if {$eNckMem<-6} {break}; # Previous criterion of Nck membrane collision.
      # Get shortest Nck distance to membrane.
      # Increment contacts cutoff, t, until atom pairs are found. 
      while {$NckNum==0} {
        set NckMem [measure contacts $t $nckSH $mem]
        set NckNum [llen [lindex $NckMem 1]]
        incr t 2
      }        
      unset -nocomplain dl   
      foreach a [lindex $NckMem 0] b [lindex $NckMem 1] {
        set ax [atomselect $nckRef "index $a" frame $ref]
        set bx [atomselect $tcr "index $b" frame $f]
        set aa [$ax get {x y z}]
        set bb [$bx get {x y z}]
        set d [vecdist [lindex $aa 0] [lindex $bb 0]]; # distance between atom pair
        lappend dl $d             
        set fNckMem2 [lindex [lsort $dl] 0]; # Shortest Nck distance to membrane           
        $ax delete; $bx delete
      }
      break; # Stop searching once the first Nck binding without contact found.
    }
  }
  puts "Done with fNck, t=$t"
  puts [format "%s.eq%d atomselect#=%d" $name [expr $f+4] [llen [atomselect list]]]
  # Output PRS dihedrals
  set eP182psi [wrapPsi [$P182e get psi]]  
  set eP184psi [wrapPsi [$P184e get psi]]
  set fP182psi [wrapPsi [$P182f get psi]] 
  set fP184psi [wrapPsi [$P184f get psi]]
  # Output: <eq frame><eNck binding><fNck binding><eMem1><eMem2><fMem1><fMem2> |
  # <eP182psi><eP184psi> | <fP182psi><fP184psi>
  puts $out [format "%6d%4d%4d %7.2f%7.2f%7.2f%7.2f | %7.2f %7.2f | %7.2f %7.2f" \
    [expr $f+4] $nRefE $nRefF $eNckMem $eNckMem2 $fNckMem $fNckMem2 \
    $eP182psi $eP184psi $fP182psi $fP184psi]

}
close $out
quit
# below not used if called with environment variable.
mol delete $tcr
}
}

