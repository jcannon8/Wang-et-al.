# viewTcrEx1.tcl: Render example tcr38b and tcr39b models (Fig 14).
# use vmd on Eire
cd /t/tcr39an
set m tcr39; # ensemble name, =tcr38 for tcr38b; =tcr39 for tcr39b
set j 0;  # model name for EFDG CD3 CTs
set k 21; # model name for AB CD3 CTs
# Scripts anTcr38rep1.tcl and anTcr38=9rep1.tcl make tcr38b85.dat and tcr39b85.dat.
# Use data from tcr38b85.dat or tcr39b85.dat to choose interesting models.
# Example line from tcr39b85.dat:
#  0.10 |  0  8.06 |  0  7.09 | A142 10.34 | A83-DPPC O14 | A123-DOPE N | 
# A123-DOPE O14 | A123-DPPC O13 | B153-POPC O11 | F199-POPS C13 | F199-POPS O13B
#
# Get the line with data for this model.
set t "$j.$k"
set sumFile ${m}b85.dat
set a "awk {{if (\$1==$t) {print \$0}}} $sumFile"; # Build awk command to get model line.
# Notice the necessary use of "eval" below
# https://stackoverflow.com/questions/28095820/tcl-expect-exec-how-to-execute-program-with-parameters
set line [eval exec $a]
# Extract features for model.
set eNck [lindex $line 2]
if {[string compare $eNck "no"]==0} {set eNck -1}
set fNck [lindex $line 5]
if {[string compare $fNck "no"]==0} {set fNck -1}
set Lck [lindex $line 8]
if {[string compare $Lck "no"]==0} {set Lck -1}
# There could be multiple deep ITAM Tyr.
set lip {}; set lip2 {}; set lip3 {} 
set c 11; # column for next deep lipid
set n 1
while {$n} {
  set lip0 [lindex $line $c]
  set n [string compare $lip0 ""]
  if {$n} {lappend lip2 $lip0}
  incr c 3
}
foreach v [split $lip2 " "] {
  lappend lip3 [lindex [split $v "-"] 0]
}
set lip [lsort -unique $lip3]
puts "$eNck $fNck $Lck $lip"
# Load TCR model.
set fName $m.$j.$k 
set tcr [mol new /t/${m}b/$fName.psf]
mol addfile /t/${m}b/$fName.eq85.rst type netcdf waitfor all
# Bind Nck to ePRS.
if {$eNck>-1} {
# Load the Nck-PRS references made by makeNckref.tcl
set eRef [mol new /t/tcr39an/NckPRS.pdb]
mol addfile /t/tcr39an/NckPRS.dcd waitfor all
# Define Nck-PRS reference atomselections. Chain E is PRS, chain K is Nck.
set eNckPRS [atomselect $eRef "name CA and chain E and resid 180 to 186" frame eNck]
set eNckSH [atomselect $eRef "chain K" frame eNck]
# PRS is Pro180-Pro186 in 6JXR chains E and F (see anEdt2.sh).
set ePRS [atomselect $tcr "name CA and chain E and resid 180 to 186"]
# Move eNckSH to bind to ePRS
set M [measure fit $eNckPRS $ePRS]; # Fit the 7 PRS CA atoms
$eNckSH move $M
# Save eNckSH.pdb
set Nckname $fName.eNck.pdb
$eNckSH writepdb $Nckname
}
# fPRS Nck binding
if {$fNck>-1} {
# Load the Nck-PRS references made by makeNckref.tcl
set fRef [mol new /t/tcr39an/NckPRS.pdb]
mol addfile /t/tcr39an/NckPRS.dcd waitfor all
# Define Nck-PRS reference atomselections. Chain E is PRS, chain K is Nck.
set fNckPRS [atomselect $fRef "name CA and chain E and resid 180 to 186" frame fNck]
set fNckSH [atomselect $fRef "chain K" frame fNck]
# PRS is Pro180-Pro186 in 6JXR chains E and F (see anEdt2.sh).
set fPRS [atomselect $tcr "name CA and chain F and resid 180 to 186"]
# Move fNckSH to bind to fPRS
set M [measure fit $fNckPRS $fPRS]; # Fit the 7 PRS CA atoms
$fNckSH move $M
# Save fNckSH.pdb
set Nckname $fName.fNck.pdb
$fNckSH writepdb $Nckname
}
animate goto 0; # Otherwise the rendering will have the wrong Nck references
# Bind Lck to ONE ITAM Tyr substrate
if {[string compare $Lck -1]} { 
# Define Lck-ITAM reference atomselections. Chain A is Lck, chain B is ITAM substrate
set lckpep [mol new /home/jcannon/lck/Lck-Peptide.pdb]
set lck [atomselect $lckpep "chain A"]
set lcksub [atomselect $lckpep all]
set lckH [atomselect $lckpep "chain A and mass>2"]; # Lck heavy atoms
# Choose residues for fitting
set t 4; # The 9 residues around pTyr for fitting.
set c [expr 285- $t]; # first Lck substrate residue to fit
set d [expr 285+ $t]; # last Lck substrate residue to fit
set pSel [atomselect $lckpep "name CA and chain B and resid $c to $d"]
# Define ITAM atomselections based on Lck above.
set chn [lindex [split $Lck ""] 0]
set n [lindex [split $Lck $chn] 1]
set a [expr $n- $t]; # first ITAM residue to fit
set b [expr $n+ $t]; # last ITAM residue to fit
set iSel [atomselect $tcr "chain $chn and name CA and resid $a to $b"]
# Fitting: fit substrate to ITAM and move lcksub
set M [measure fit $pSel $iSel]
$lcksub move $M
# Save lck.pdb
set Lckname $fName.lck.pdb
$lcksub writepdb $Lckname
}
##################
# Render the TCR complex with bound Nck and Lck, and deep Tyr.
display resize 560 300; # This size (width, height) allows four models across full page width.
display projection Orthographic
display depthcue off
color Display Background white
axes location Off
rotate x by 90
# Adjust the representations
mol delrep 0 $tcr
# Change default representation to tube otherwise they are lines
mol representation Tube 0.50 12.0; # thicker than default.
# CD3zeta: chain a green
mol color ColorID 7
mol selection chain A
mol addrep $tcr
# CD3zeta: chain b purple
mol color ColorID 11
mol selection chain B
mol addrep $tcr
# CD3delta: chain d red
mol color ColorID 1
mol selection chain D
mol addrep $tcr
# CD3epsilon: chain f grey
mol color ColorID 2
mol selection chain F
mol addrep $tcr
# CD3gamma: chain g orange
mol color ColorID 3
mol selection chain G
mol addrep $tcr
# TCRalpha: chain m pink
mol color ColorID 9
mol selection chain M
mol addrep $tcr
# TCRbeta: chain n yellow3
mol color ColorID 18
mol selection chain N
mol addrep $tcr
# CD3epsilon: chain e blue
mol color ColorID 0
mol selection chain E
mol addrep $tcr
# Membrane P atoms
mol color Name
mol representation VDW 1.0 12.0
mol selection name P
mol material Transparent
mol addrep $tcr
# Structural lipids black Licorice
mol color ColorID 16
mol material Opaque
mol selection chain S
mol representation Licorice
mol addrep $tcr
##############
# Render eNck magenta2 NewCartoon
if {$eNck>-1} {
mol delrep 0 $eRef
mol selection all
mol material Opaque
mol color ColorID 28
mol representation NewCartoon
mol addrep $eRef
}
# Render fNck magenta2 NewCartoon
if {$fNck>-1} {
mol delrep 0 $fRef
mol selection all
mol material Opaque
mol color ColorID 28
mol representation NewCartoon
mol addrep $fRef
}
##############
# Render Lck ochre NewCartoon
if {[string compare $Lck -1]} {
mol delrep 0 $lckpep
mol selection all
mol material Opaque
mol color ColorID 14
mol representation NewCartoon
mol addrep $lckpep
}
##############
# If deep ITAM Tyr black VDW
if {[string compare $lip -1]} {
# Parse ITAM Tyr from lip list above.
foreach d $lip {
  set chn2 [lindex [split $d ""] 0]
  set res [lindex [split $d $chn2] 1]
  mol color ColorID 16
  mol material Opaque
  mol selection chain $chn2 and resid $res
  mol representation VDW
  mol addrep $tcr
}
}
###########################################################
# Position system to fit into display.
# The rotation makes the display have X-axis, -Z-axis in the XY display plane.
set S [molinfo top get scale_matrix]
# For tcr39.5.20 with Nck and Lck, scale factor is 0.022 and nothing is visible.
#
# Center the system by changing coordinates.
set sel [atomselect $tcr all]
set com [measure center $sel]
set negcom [vecscale -1.0 $com]
# Change coordinates for all molecules.
foreach mm [molinfo list] { 
set selM [atomselect $mm all]
$selM moveby $negcom
}
# ^ That leaves the extracellular membrane leaflet at Z=0.
# That leaflet needs to move up to the top of the display expose everything.
set ext [measure minmax [atomselect $tcr "protein or lipid"]]
set dZ [expr -1 * ([lindex $ext 1 2] + [lindex $ext 0 2] -20)]
foreach mm [molinfo list] {
set selM [atomselect $mm all]
$selM moveby [list 0 0 $dZ]; # perfect for tcr39.5.20
}
# Above works for tcr39.5.20, tcr39.5.11, tcr39.0.21, tcr38.3.21
