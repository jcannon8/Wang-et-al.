# tcr33restraint.tcl: Change the tcr30 lipid restraints for tcr33.
# After makeTcr33.tcl fixed the errors in the CHARM_GUI output for tcr30,
# the atom count increased and the CHARM_GUI restraints for lipids need adjustment.
# Increase the atom indices by 176 for all restraints.
cd /t/tcr33
set delta 176
set in [open /t/tcr30/charmm-gui-7001280423/amber/dihe.restraint r]
set out [open dihe.restraint w]
set inData [read $in]; close $in
foreach line [split $inData "\n"] {
  if {[lsearch $line "iat*"]==0} {
    set s [split $line "=, "]
    puts $out [format "    iat=%d, %d, %d, %d," \
      [expr [lindex $s 5]+$delta] \
      [expr [lindex $s 7]+$delta] \
      [expr [lindex $s 9]+$delta] \
      [expr [lindex $s 11]+$delta] ]
    
  } else {puts $out $line}
} 
close $out
 
