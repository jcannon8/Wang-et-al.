# anPRS4.tcl: Analyze PRS dihedral angles in edt and fgt.
cd /t/tcr16/eds
# Get backbone dihedral angles with VMD.
# vmd -dispdev none
# Wrap angles as decribed by Hollingsworth and Karplus 2010
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3061398/
# Wrapped Ramachandran plots reset phi<0 to phi+360 and psi< -120 to psi+360.
# Then plot -120 to 240 psi on y-axis, 0 to 360 phi on x-axis.
proc phiOut {phi outFile last} {
# Output one line of phi angles with wrapping.
  for {set j 0} {$j<=8} {incr j} {
    set a [lindex $phi $j]
    set b [expr $a<0?$a+360:$a]
    puts -nonewline $outFile [format "%6.3f " $b]
  }
  puts $outFile $last
}
proc psiOut {psi outFile last} {
# Output one line of psi angles with wrapping.
  for {set j 0} {$j<=8} { incr j} {
    set a [lindex $psi $j]
    set b [expr $a<-120?$a+360:$a]
    puts -nonewline $outFile [format "%6.3f " $b]
  }
  puts $outFile $last
}
# anEdt11.tcl used cpptraj to get the dihedral angles but did not wrap them.
set outPhi [open "edtfgt.phi.dat" w]
set outPsi [open "edtfgt.psi.dat" w]
# 11 edt references. 
set Edtmod {84 171 30 115 35 90 13 41 32 153 198}
foreach mod $Edtmod {
  set edt [mol new /t/tcr16/edt/edt.$mod.top type parm7 waitfor all]
  # Use the 205 ns frame
  mol addfile /t/tcr16/edt/edt.$mod.eq52.rst type netcdf waitfor all
  set phi [[atomselect $edt "name CA and resid 130 to 138"] get phi]
  phiOut $phi $outPhi 1
  set psi [[atomselect $edt "name CA and resid 130 to 138"] get psi]
  psiOut $psi $outPsi 1
  mol delete $edt
}
# 14 fgt references. 
set Fgtmod {185 121 191 180 137 75 14 153 93 28 109 122 42 45}
foreach mod $Fgtmod {
  set fgt [mol new /t/tcr16/fgt/fgt.$mod.top type parm7 waitfor all]
  # Use the 205 ns frame
  mol addfile /t/tcr16/fgt/fgt.$mod.eq52.rst type netcdf waitfor all
  set phi [[atomselect $fgt "name CA and resid 58 to 66"] get phi]
  phiOut $phi $outPhi 2
  set psi [[atomselect $fgt "name CA and resid 58 to 66"] get psi]
  psiOut $psi $outPsi 2
  mol delete $fgt
}
# Original NckPRS9803 reference
set 2jxb3 [mol new /home/jcannon/tcr2/cd3e/2jxb3.c0.pdb waitfor all]
set phi [[atomselect $2jxb3 "name CA and resid 6 to 14"] get phi]
phiOut $phi $outPhi 3
set psi [[atomselect $2jxb3 "name CA and resid 6 to 14"] get psi]
psiOut $psi $outPsi 3
#
close $outPhi
close $outPsi
###########################################################
# Plot together with edr/fgr in anPRS3.tcl

