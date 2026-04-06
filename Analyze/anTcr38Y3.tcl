# anTcr38Y3.tcl: Investigate ITAM Tyr membrane contacts for tcr38.*.** models.
# vmd -dispdev none
# Single 202 ns (eq45) or 405 ns (eq85) frames analyzed. 
# Consider contacts of heavy atoms less than or equal to 3.0 angstroms apart.
cd /t/tcr38an
# This outputs a comprehensive list, *.memy85a.dat
set m tcr39; # =tcr38 for tcr38b; =tcr39 for tcr39b
set out [open $m.memy85a.dat w]; # For 425 ns (eq85) frame analysis
###########################################################
set itams {{A 72} {A 83} {A 111} {A 123} {A 142} {A 153} 
{B 72} {B 83} {B 111} {B 123} {B 142} {B 153}  
{D 149} {D 160} {E 188} {E 199} {F 188} {F 199} {G 160} {G 171} }
#
for {set j 0} {$j<=11} {incr j} {
foreach k {00 01 02 10 11 12 20 21 22} {
set fName $m.$j.$k  
# Skip missing models. These are 425 ns frames.
if {[file exist /t/${m}b/$fName.psf]==0} continue
if {[file exist /t/${m}b/$fName.eq85.rst ]==0} continue
set tcr [mol new /t/${m}b/$fName.psf]
mol addfile /t/${m}b/$fName.eq85.rst type netcdf waitfor all
set mem [atomselect $tcr "chain J and mass>2"];  # Only membrane heavy atoms
# Find shortest Tyr-Mem distances in this model
foreach s $itams {
set chn [lindex $s 0]
set n [lindex $s 1]
set ySel [atomselect $tcr "chain $chn and resid $n and mass>2"]; # Only Tyr heavy atoms
set cMem [measure contacts 3.0 $mem $ySel]
set conNum [llength [lindex $cMem 1]]; # Number of membrane contacts
if {$conNum>0} {
  # If mem contacts found, report number
  puts $out "$fName $chn $n $conNum"
  # Report each atom name pair
  foreach a2 [lindex $cMem 0] b2 [lindex $cMem 1] {
    set ax [atomselect $tcr "index $a2"]
    set bx [atomselect $tcr "index $b2"]
    set aa [$ax get {x y z}]
    set bb [$bx get {x y z}]
    set d [vecdist [lindex $aa 0] [lindex $bb 0]]; # distance between atom pair
    puts $out [format "  %s %s - %s %s %4.2f"\
      [$ax get resname] [$ax get name] [$bx get resname] [$bx get name] $d]
    # These all have a TYR residue as the second in these pairs.     
    $ax delete
    $bx delete
  }  
}
$ySel delete
} 
}
}
$mem delete
close $out
###########################################################
# This outputs a list focusing ITAM substrate membane atom, *.memy85b.dat
# vmd -dispdev none
set m tcr39; # =tcr38 for tcr38b; =tcr39 for tcr39b
set out [open $m.memy85b.dat w]; # For 425 ns (eq85) frame analysis
set itams {{A 72} {A 83} {A 111} {A 123} {A 142} {A 153} 
{B 72} {B 83} {B 111} {B 123} {B 142} {B 153}  
{D 149} {D 160} {E 188} {E 199} {F 188} {F 199} {G 160} {G 171} }
#
for {set j 0} {$j<=11} {incr j} {
foreach k {00 01 02 10 11 12 20 21 22} {
set fName $m.$j.$k  
# Skip missing models. These are 425 ns frames.
if {[file exist /t/${m}b/$fName.psf]==0} continue
if {[file exist /t/${m}b/$fName.eq85.rst ]==0} continue
set tcr [mol new /t/${m}b/$fName.psf]
mol addfile /t/${m}b/$fName.eq85.rst type netcdf waitfor all
set mem [atomselect $tcr "chain J and mass>2"];  # Only membrane heavy atoms
# Find shortest Tyr-Mem distances in this model
foreach s $itams {
set chn [lindex $s 0]; # Subunit chain
set n [lindex $s 1];   # Tyr residue
set ySel [atomselect $tcr "chain $chn and resid $n and mass>2"]; # Only Tyr heavy atoms
set cMem [measure contacts 3.0 $mem $ySel]
set conNum [llength [lindex $cMem 1]]; # Number of membrane contacts
if {$conNum>0} {
  # If mem contacts found, report ITAM TYR, lipid resname, and atom name
  # Report each atom name pair
  foreach a2 [lindex $cMem 0] b2 [lindex $cMem 1] {
    set ax [atomselect $tcr "index $a2"];  # membrane residue index
    set memRes [$ax get resname];  # membrane residue name 
    set memAtom [$ax get name];     # membrane residue atom
    puts $out [format "%s %s%s-%s %s" $fName $chn $n $memRes $memAtom]     
    $ax delete
  }  
}
$ySel delete
} 
}
}
$mem delete
close $out
###########################################################
# Example *.memy85b.dat output
tcr38.11.22 B83-DOPC O13
tcr38.11.22 E199-POPC O14
tcr38.11.22 D149-POPC O14
tcr38.11.22 F188-CHL1 O3 
