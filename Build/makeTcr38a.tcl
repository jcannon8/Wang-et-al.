# makeTcr38a.tcl: Test adding four CD3 cytoplasmic tails from CT ensembles to tcr38
# made by removeExt4.tcl. 
# This script outputs potential framesets to make models with reported collisions.
# This is derived from makeTcr36a.tcl.
# Instead of using only 2jxb3 frame 9803 to bind Nck, this uses many references 
# Contacts are heavy atoms separated less than or equal to 2.0 angtroms. 
# This outputs a table of fusion data in tcr38.$model.dat.
# Use vmd -dispdev none with exported nLast and model variables. 
cd /t/tcr38
# Number of requested framesets from nlast environment variable.
set nlast $env(nlast)
set model $env(model)
proc base3 {m} {
  # Convert decimal m into a base 3 number
  set r27 [expr int($m/27)]
  set n [expr $m - $r27 * 27]
  set r9 [expr int($n/9)] 
  set n [expr $n - $r9 *9]
  set r3 [expr int($n/3)]
  set r1 [expr $n - $r3 *3]
  return "$r27 $r9 $r3 $r1"
}
set mb3 [base3 $model]
# The four-digit mb3 is (a)(b)(fe)(dg)
set fRef [lindex $mb3 2]; set eRef $fRef
set dRef [lindex $mb3 3]; set gRef $dRef
# tcr38 model number comes from treating refs in base3. As in (a)(b)(fe)(dg)
set model [expr 3*$fRef + $dRef]
puts "model 38.$model framesets, refs=$fRef.$dRef.$gRef.$eRef"
set out [open tcr38.$model.dat w]
# tcr38: removeExt4.tcl removed extracellular protein and water from tcr35c.eq70a.rst (322 ns)
set tcr [mol new tcr38.psf]
mol addfile tcr38.pdb
# ^ That gets the whole system with 6JXR chain names and residue numbering.
#
# This script adds tails from trajectories, adds Nck to chains EF,
# and tests whether they are feasible based on contacts with membrane and other CTs.
# This samples from four ff14IDPSFF CD3 CT reference trajectories 
# in the /home/jcannon/tcr15/ directory.
# CT trajectories: eps.c[0-11] (12), delta3.c[0-8] (9), gamma3.c[0-8] (9), zeta3.c[0-8] (9)
######################################################
# This is where multiple Nck-PRS references are loaded.
# Load collection of 11 Nck references from edr.*.eq54.rst frames.
set edrModels "154 16 75 174 147 20 176 197 62 107 201"
set j 0; # j=index for Nck references
foreach ref $edrModels {
  set nckRef($j) [mol new /t/tcr16/edm/edm.$ref.psf]
  mol addfile /t/tcr16/edr/edr.$ref.eq54.rst type netcdf waitfor all
  # Nck chain K, resid 31-86. PRS chain E, resid 180-188.
  set PRSref($j) [atomselect $nckRef($j) "name CA and chain E and resid 180 to 186"]
  set nckSH($j) [atomselect $nckRef($j) "chain K and not hydrogen"]
  # Both Nck in chain K and the PRS CA atoms in chain E need to move.
  [atomselect $nckRef($j) all] set beta 0
  $PRSref($j) set beta 1
  [atomselect $nckRef($j) "chain K"] set beta 1
  set nckPRSref($j) [atomselect $nckRef($j) "beta 1"]
  incr j
}
set lastEdrJ $j
# Load collection of 14 Nck references from fgr.*.eq53.rst frames.
set fgrModels "105 177 199 93 10 83 168 142 30 55 154 117 86 201"
foreach ref $fgrModels {
  set nckRef($j) [mol new /t/tcr16/fgm/fgm.$ref.psf]
  mol addfile /t/tcr16/fgr/fgr.$ref.eq53.rst type netcdf waitfor all
  # Nck chain K, resid 31-86. PRS chain F, resid 180-186.
  set PRSref($j) [atomselect $nckRef($j) "name CA and chain F and resid 180 to 186"]
  set nckSH($j) [atomselect $nckRef($j) "chain K and not hydrogen"]
  # Both Nck in chain K and the PRS CA atoms in chain F need to move.
  [atomselect $nckRef($j) all] set beta 0
  $PRSref($j) set beta 1
  [atomselect $nckRef($j) "chain K"] set beta 1
  set nckPRSref($j) [atomselect $nckRef($j) "beta 1"]
  incr j
}
set lastFgrJ $j
####################
# Testing membrane collision of heavy atoms.
set lipids [atomselect $tcr "chain J"]
llength [lsort -unique [$lipids get resid]]
# ^ There are 1600 lipid residues
# Get average membrane z-dimension for inner leaflet.
set memDim [measure minmax $lipids]
set min [lindex $memDim 0 2]
set max [lindex $memDim 1 2]
set aveMem [expr ($min+$max)/2]
# Get inner phosphate z-dimension
set phos [atomselect $tcr "name P"]
set memDim2 [measure minmax $phos]
set maxP31 [lindex $memDim2 1 2]
# cytoplasmic leaflet heavy atoms
set mem [atomselect $tcr "chain J and z>$aveMem and mass>2"]
set conX 1.0;  # Collisions with membrane heavy atoms are <= 1.0 angstrom
set maxCon 0;  # maximum number of contacts acceptable
set maxCon2 0;  # maximum number of contacts acceptable with Nck
# Finding membrane collisions is tedious. Save non-collision frames after first pass
# in files called chain[A,B,D,E,F,G].$ref.
# Also check that CTs do not stick out of box, particularly chain a and b in z-dimension.
set all [atomselect $tcr all]
set box [measure minmax $all]
set maxZ [lindex $box 1 2]
###########################################################
# Add eps CT to chain f.
set m [file exists chainF.$fRef]
if {$m==0} {set outF [open chainF.$fRef w]}
set cd3f [mol new /home/jcannon/tcr15/eps.top type parm7]
mol addfile /home/jcannon/tcr15/eps.c$fRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/eps.c$fRef.eq3.cdf type netcdf waitfor all
set lastframe [molinfo $cd3f get numframes]
# Cytoplasmic tail is all of eps. New addition to tcr36 are residues 2-57.
set CTf [atomselect $cd3f "resid 2 to 57 and mass>2"]
set CTfa [atomselect $cd3f "resid 4 to 57"]
[atomselect $tcr "chain F and name N"] get {resname resid}
# CD3 epsilon: 6JXR chain f resid 33-156,  tcr36 resid 123-156.
# Change fLim from 153 to 152 to consider more chain f fusions.
set fLim [[atomselect $tcr "chain F and name N and resid 152"] get z]
# Fit between tcr and eps uses SER LYS ASN ARG residues in each.
# Align resid 2-5 of cd3f and resid 152 to 155 of tcr36.
set End [atomselect $tcr "chain F and name CA and resid 152 to 155"]
set Start [atomselect $cd3f "resid 2 to 5 and name CA"]
set cd3CD [atomselect $cd3f all]
# For testing Nck binding 
set fPRS [atomselect $cd3f "name CA and resid 30 to 36"]
# notPRS atoms are more than two residues away from PRS
set fNotPRS [atomselect $cd3f "not hydrogen and not resid 28 to 38"] 
for {set f 0 } {$f < $lastframe} {incr f} {
	$Start frame $f; $cd3CD frame $f
  $CTf frame $f; $CTfa frame $f 
  $fPRS frame $f; $fNotPRS frame $f
	# RMS fit Start to End
	set M [measure fit $Start $End]
	# Move CD3 CT to fit tcr 
	$cd3CD move $M
 	if {$m==0} {
    # Determine membrane collision
    $CTf frame $f; $CTfa frame $f
    set limit [measure minmax $CTfa]
    set min [lindex $limit 0 2]
    if {$min <$fLim} {continue}
    set con [measure contacts $conX $CTf $mem]
  	set flCon [llength [lindex $con 0]]
    if {$flCon<=$maxCon} {
      puts [format "%d Fframe, min=%5.2f, fLim=%5.2f" $f $min $fLim]
      lappend fFL1 $f
    } 
  }
}
if {$m==0} {puts $outF $fFL1; close $outF}
######################################################
# File chainF.0 already has frames from eps that clear the membrane.
# Make nckF.$fRef if necessary with frames that bind Nck without collisions.
set m [file exists nckF.$fRef]
if {$m==0} {
  set outF [open nckF.$fRef w]
  # Get the frames that cleared the membrane
  set inF [open chainF.$fRef r]
  set fFL1 [read $inF]; close $inF
  # Check only frames that clear the membrane for collisions with bound Nck.
  foreach f $fFL1 {
  	$Start frame $f; $cd3CD frame $f
    $fPRS frame $f; $fNotPRS frame $f
  	# RMS fit eps Start to tcr36 End
  	set M [measure fit $Start $End]
  	# Move eps CD3 CT to fit tcr 
  	$cd3CD move $M
    # Bind Nck references and check collisions with chain F.
    # For each Nck reference, collect fNckRMSD and conNumF.
    unset -nocomplain nRMSD nFcon
    for {set ref 0} {$ref<$lastFgrJ} {incr ref} {
      # Move refNck to bind to fPRS
      set M [measure fit $PRSref($ref) $fPRS]
      # This moves chain F PRS CA atoms and all of Nck so that the RMSD of 
      # the moved atoms can calculate the PRS RMSD fit.
      $nckPRSref($ref) move $M
      set RMSD [measure rmsd $PRSref($ref) $fPRS]  
      # Check contacts with Nck SH3.1 
      set con [measure contacts 2.0 $nckSH($ref) $fNotPRS]
      set conNumF [llen [lindex $con 1]]
      # Append list of Nck binding parameters 
      lappend nRMSD $RMSD; lappend nFcon $conNumF
    }
    # Sort Nck binding parameters, get the one with lowest chain F contacts
    set bestFcon [lindex [lsort $nFcon] 0]
    set bestFconi [lindex [lsort -indices $nFcon] 0]
    set bRMSD [lindex $nRMSD $bestFconi]
    # Get the membrane contacts to best
    set con [measure contacts 1.0 $nckSH($bestFconi) $mem]
    set bMem [llength [lindex $con 0]]
    # Only save frames that have zero chain F or membrane contacts.
    # Low bRMSD could also be used as a criteria because some are above 3.0. 
    if {[expr $bestFcon + $bMem] ==0 && $bRMSD < 3.0} {
      puts $outF [format "%4d %3d %3d %3d %5.3f" $f $bestFcon $bMem $bestFconi $bRMSD]
    }
  }
  close $outF
}
# Fill fFL with viable chain F frames and fNck with corresponding Nck refs.
unset -nocomplain fFL fNck
set inF [open nckF.$fRef r]
set inData [read $inF]
close $inF
foreach line [split $inData "\n"] {
  if {$line == ""} break
  lappend fFL [lindex $line 0] 
  lappend fNck [lindex $line 3]
}  
set fLen [llen $fFL]
###########################################################
# New Add eps CT to chain e.
set m [file exists chainE.$eRef]
if {$m==0} {set outE [open chainE.$eRef w]}
set cd3e [mol new /home/jcannon/tcr15/eps.top type parm7]
mol addfile /home/jcannon/tcr15/eps.c$eRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/eps.c$eRef.eq3.cdf type netcdf waitfor all
set lastframe [molinfo $cd3e get numframes]
# Cytoplasmic tail is all of eps. New addition to tcr36 are residues 2-57.
set CTe [atomselect $cd3e "resid 2 to 57 and mass>2"]
set CTea [atomselect $cd3e "resid 4 to 57"]
[atomselect $tcr "chain E and name N"] get {resname resid}
# CD3 epsilon: 6JXR chain e resid 125-155 in tcr36
# No CT atoms should go deeper into the membrane than the last-3 N atom in 6JXR.
set eLim [[atomselect $tcr "chain E and name N and resid 152"] get z]
# Fit between tcr and eps uses SER LYS ASN ARG residues in each.
# Align resid 2-5 of cd3e and resid 152 to 155 of tcr36
set End [atomselect $tcr "chain E and name CA and resid 152 to 155"]
set Start [atomselect $cd3e "resid 2 to 5 and name CA"]
set cd3CD [atomselect $cd3e all]
# For testing Nck binding
set ePRS [atomselect $cd3e "name CA and resid 30 to 36"]
# notPRS atoms are more than two residues away from PRS
set eNotPRS [atomselect $cd3e "not hydrogen and not resid 28 to 38"] 
for {set f 0 } {$f < $lastframe} {incr f} {
  $Start frame $f; $cd3CD frame $f
  $CTe frame $f; $CTea frame $f
  $ePRS frame $f; $eNotPRS frame $f
	# RMS fit Start to End
	set M [measure fit $Start $End]
	# Move CD3 CT to fit tcr 
	$cd3CD move $M
	if {$m==0} {
    # Determine membrane collision	
    set limit [measure minmax $CTea]
    set min [lindex $limit 0 2]
    if {$min <$eLim} {continue}
    set con [measure contacts $conX $CTe $mem]
  	set flCon [llength [lindex $con 0]]
    if {$flCon<=$maxCon} { 
      puts [format "%d Eframe, min=%5.2f, eLim=%5.2f" $f $min $eLim]
      lappend eFL1 $f
    }
  } 			
}
if {$m==0} {puts $outE $eFL1; close $outE}
###########################################################
# File chainE.0 already has frames from eps that clear the membrane.
# Make nckE.$eRef if necessary with frames that bind Nck without collisions.
set m [file exists nckE.$eRef]
if {$m==0} {
  set outE [open nckE.$eRef w]
  # Get the frames that cleared the membrane
  set inE [open chainE.$eRef r]
  set eFL1 [read $inE]; close $inE
  # Check only frames that clear the membrane for collisions with bound Nck.
  foreach f $eFL1 {
    $Start frame $f; $cd3CD frame $f
    $ePRS frame $f; $eNotPRS frame $f
  	# RMS fit eps Start to tcr36 End
  	set M [measure fit $Start $End]
  	# Move eps CD3 CT to fit tcr 
  	$cd3CD move $M
    # Bind Nck references and check collisions with chain E.
    # For each Nck reference, collect eNckRMSD and conNumE.
    unset -nocomplain nRMSD nEcon
    for {set ref 0} {$ref<$lastFgrJ} {incr ref} {
      # Move refNck to bind to ePRS
      set M [measure fit $PRSref($ref) $ePRS]
      # This moves chain E PRS CA atoms and all of Nck so that the RMSD of 
      # the moved atoms can calculate the PRS RMSD fit.
      $nckPRSref($ref) move $M
      set RMSD [measure rmsd $PRSref($ref) $ePRS]  
      # Check contacts with Nck SH3.1 
      set con [measure contacts 2.0 $nckSH($ref) $eNotPRS]
      set conNumE [llen [lindex $con 1]]
      # Append list of Nck binding parameters 
      lappend nRMSD $RMSD; lappend nEcon $conNumE
    }
    # Sort Nck binding parameters, get the one with lowest chain E contacts
    set bestEcon [lindex [lsort $nEcon] 0]
    set bestEconi [lindex [lsort -indices $nEcon] 0]
    set bRMSD [lindex $nRMSD $bestEconi]
    # Get the membrane contacts to best
    set con [measure contacts 1.0 $nckSH($bestEconi) $mem]
    set bMem [llength [lindex $con 0]]
    # Only save frames that have zero chain E or membrane contacts.
    # Low bRMSD could also be used as a criteria 
    if {[expr $bestEcon + $bMem] ==0 && $bRMSD < 3.0} {
      # Output: <eps frame><E contacts><membrane contacts><Nck reference #><PRS RMSD>
      puts $outE [format "%4d %3d %3d %3d %5.3f" $f $bestEcon $bMem $bestEconi $bRMSD]
    }
  }
  close $outE
}
# Fill eFL with viable chain E frames and eNck with corresponding Nck refs.
unset -nocomplain eFL eNck
set inF [open nckE.$eRef r]
set inData [read $inE]
close $inE
foreach line [split $inData "\n"] {
  if {$line == ""} break
  lappend eFL [lindex $line 0] 
  lappend eNck [lindex $line 3]
}  
set eLen [llen $eFL]
###########################################################
# Add delta3 CT to chain d in tcr model.
set m [file exists chainD.$dRef]
if {$m==0} {set outD [open chainD.$dRef w]}
set delta [mol new /home/jcannon/tcr15/delta.top type parm7]
mol addfile /home/jcannon/tcr15/delta3.c$dRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/delta3.c$dRef.eq3.cdf type netcdf waitfor all
set lastframe [molinfo $delta get numframes]
# Cytoplasmic tail is residues 5-46 in delta.
set CTd [atomselect $delta "resid 3 to 46 and mass>2"]
set CTda [atomselect $delta "resid 2 to 46"]
[atomselect $tcr "chain D and name N"] get {resname resid}
# CD3 delta: 6JXR chain d resid 98 to 129 tcr36 
# No CT atoms should go deeper into the membrane than the last-3 N atom in 6JXR.
set dLim [[atomselect $tcr "chain D and name N and resid 126"] get z]
# delta3 has delta residues 100-171
# Align resid 1-4 of delta3 and resid 126 to 129 of tcr36
# Fit between tcr and delta uses ALA GLY HIS GLU residues in each.
set End [atomselect $tcr "chain D and name CA and resid 126 to 129"]
set Start [atomselect $delta "resid 1 to 4 and name CA"]
set cd3CD [atomselect $delta all]
for {set f 0 } {$f < $lastframe} {incr f} {
	$Start  frame $f; $cd3CD  frame $f
	# RMS fit Start to End
	set M [measure fit $Start $End]
	# Move CD3 CT to fit tcr 
	$cd3CD move $M
	if {$m==0} {
    # Determine membrane collision	
    $CTd frame $f; $CTda frame $f
    set limit [measure minmax $CTda]
    set min [lindex $limit 0 2]
    if {$min <$dLim} {continue}
    set con [measure contacts $conX $CTd $mem]
  	set flCon [llength [lindex $con 0]]
    if {$flCon<=$maxCon} {
      puts [format "%d Dframe, min=%5.2f, dLim=%5.2f" $f $min $dLim]
      lappend dFL $f
    }
  } 
}
if {$m==0} {puts $outD $dFL; close $outD;}
set inD [open chainD.$dRef r]
set dFL [read $inD]; close $inD
set dLen [llen $dFL]
###########################################################
# Add gamma3 CT to chain g in tcr model.
set m [file exists chainG.$gRef]
if {$m==0} {set outG [open chainG.$gRef w]}
set gamma [mol new /home/jcannon/tcr15/gamma.top type parm7]
mol addfile /home/jcannon/tcr15/gamma3.c$gRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/gamma3.c$gRef.eq3.cdf type netcdf waitfor all
set lastframe [molinfo $gamma get numframes]
# Cytoplasmic tail is residues 5-48 in gamma.
set CTg [atomselect $gamma "resid 3 to 48 and mass>2"]
set CTga [atomselect $gamma "resid 3 to 48"]
# CD3 gamma: 6JXR chain g resid 108-138
[atomselect $tcr "chain G and name N"] get {resname resid}
# No CT atoms should go deeper into the membrane than the last-3 N atom in 6JXR.
set gLim [[atomselect $tcr "chain G and name N and resid 135"] get z]
# gamma3 has gamma residues 111-182
# Align resid 1-4 of gamma and resid 135 to 138 of tcr36
# Fit between tcr and gamma uses PHE ILE ALA GLY residues in each.
set End [atomselect $tcr "chain G and name CA and resid 135 to 138"]
set Start [atomselect $gamma "resid 1 to 4 and name CA"]
set cd3CD [atomselect $gamma all]
for {set f 0 } {$f < $lastframe} {incr f} {
	$Start frame $f; $cd3CD frame $f
	# RMS fit Start to End
	set M [measure fit $Start $End]
	# Move CD3 CT to fit tcr 
	$cd3CD move $M
	if {$m==0} {
    # Determine membrane collision	
    $CTg frame $f; $CTga frame $f
    set limit [measure minmax $CTga]
    set min [lindex $limit 0 2]
    if {$min<$gLim} {continue}
    set con [measure contacts $conX $CTg $mem]
  	set flCon [llength [lindex $con 0]]
    if {$flCon<=$maxCon} {
      puts [format "%d Gframe, min=%5.2f, gLim=%5.2f" $f $min $gLim]
      lappend gFL $f
    }
  } 	
}
if {$m==0} {puts $outG $gFL; close $outG}
set inG [open chainG.$gRef r]
set gFL [read $inG]; close $inG
set gLen [llen $gFL]
###########################################################
# Now test collisions CTs that do not collide with membrane.
# Load lists of frames that did not collide with membrane.
set n 0;		  # number of successful frameset, found
# nlast=requested number of framesets from environment variable.
set i 0;		  # number of frame sets tested, tries
set conX 1.0;  # Collisions with all CT atoms are <= 1.0 angstrom
set maxCon 5;  # Maximum tolerable contacts between CTs
set maxCon2 0;  # Maximum tolerable contacts to Nck
set memAll [atomselect $tcr "chain J and z>$aveMem"]
while {$n<$nlast } {
	incr i
  # Output tries and found framesets every 100 tries.
  if {[expr $i % 100]==0} {puts "$i $n"}
  # set f d g e a b to random frame that does not collide with membrane.
	set if [expr round(rand()*($fLen-1))]; set f [lindex $fFL $if]
	set id [expr round(rand()*($dLen-1))]; set d [lindex $dFL $id]
	set ig [expr round(rand()*($gLen-1))]; set g [lindex $gFL $ig]
	set ie [expr round(rand()*($eLen-1))]; set e [lindex $eFL $ie] 
  $CTfa frame $f; $CTda frame $d; $CTga frame $g; $CTea frame $e;
  $fPRS frame $f; $ePRS frame $e
	# Test f and d collision 
	set con [measure contacts $conX $CTfa $CTda]
  # con now has { {list of CTf atoms} {list of CTd atoms} }
  set fd [llength [lindex $con 0]]
  if {$fd>$maxCon} {continue}
	# Test f and g collision 
	set con [measure contacts $conX $CTfa $CTga]
	set fg [llength [lindex $con 0]]
  if {$fg>$maxCon} {continue}
	# Test f and e collision 
	set con [measure contacts $conX $CTfa $CTea]
	set fe [llength [lindex $con 0]]
  if {$fe>$maxCon} {continue}
 	# Test d and g collision 
	set con [measure contacts $conX $CTda $CTga] 
	set dg [llength [lindex $con 0]]
  if {$dg>$maxCon} {continue}
	# Test d and e collision 
	set con [measure contacts $conX $CTda $CTea]
	set de [llength [lindex $con 0]]
  if {$de>$maxCon} {continue}
	# Test g and e collision 
	set con [measure contacts $conX $CTga $CTea]
	set ge [llength [lindex $con 0]]
  if {$ge>$maxCon} {continue}
  puts "no CT to CT collision, test Nck collision"
  ###############################################
  # A frameset was found that met the CT to CT collision criteria.
  # Now check bound Nck collision with CTs
  # Bind Nck to chain F eps fPRS
  # Get Nck reference from fNck using index if.
  set refF [lindex $fNck $if]
  # Bind Nck to chain F PRS using fPRS set above.
  set M [measure fit $PRSref($refF) $fPRS]
  # This moves chain F PRS CA atoms and all of Nck so that the RMSD of 
  # the moved atoms can calculate the PRS RMSD fit.
  $nckPRSref($refF) move $M
  set RMSD [measure rmsd $PRSref($refF) $fPRS]   
  # Check if chain F-bound Nck contacts other CTs (d,g,e)
  set con [measure contacts $conX $nckSH($refF) $CTda]
  set nFD [llength [lindex $con 0]]
  if {$nFD>$maxCon2} {continue}
  set con [measure contacts $conX $nckSH($refF) $CTga]
  set nFG [llength [lindex $con 0]]
  if {$nFG>$maxCon2} {continue}
  set con [measure contacts $conX $nckSH($refF) $CTea]
  set nFE [llength [lindex $con 0]]
  if {$nFE>$maxCon2} {continue}
  puts "Chain F-bound Nck OK, index=$refF, RMSD=$RMSD"
  # Bind Nck to chain E eps ePRS 
  # Get Nck reference from eNck using index ie.
  set refE [lindex $eNck $ie]
  # Bind Nck to chain E PRS using ePRS set above.
  set M [measure fit $PRSref($refE) $ePRS]
  # This moves chain E PRS CA atoms and all of Nck so that the RMSD of 
  # the moved atoms can calculate the PRS RMSD fit.
  $nckPRSref($refE) move $M
  set RMSD [measure rmsd $PRSref($refE) $ePRS]   
  # Check if chain E-bound Nck contacts other CTs (f,d,g)
  set con [measure contacts $conX $PRSref($refE) $CTfa]
  set nEF [llength [lindex $con 0]]
  if {$nEF>$maxCon2} {continue}
  set con [measure contacts $conX $PRSref($refE) $CTda]
  set nED [llength [lindex $con 0]]
  if {$nED>$maxCon2} {continue}
  set con [measure contacts $conX $PRSref($refE) $CTga]
  set nEG [llength [lindex $con 0]]
  if {$nEG>$maxCon2} {continue}
  puts "Chain E-bound Nck OK, index=$refE, RMSD=$RMSD"
  # All collisions tested except for Nck-Nck collision.
  # Now check all atom contacts with membrane
  # Above and saved in the files are only heavy atom membrane contacts.
  # Additional max all-atom contacts.
  set strict 1
  set con [measure contacts $conX $CTfa $memAll]
 	set fmem [llength [lindex $con 0]]
  if {$strict==1 && $fmem >$maxCon} {continue}  
  set con [measure contacts $conX $CTda $memAll]
 	set dmem [llength [lindex $con 0]] 
   if {$strict==1 && $dmem >$maxCon} {continue}   
  set con [measure contacts $conX $CTga $memAll]
 	set gmem [llength [lindex $con 0]] 
  if {$strict==1 && $gmem >$maxCon} {continue} 
  set con [measure contacts $conX $CTea $memAll]
 	set emem [llength [lindex $con 0]]
  if {$strict==1 && $emem >$maxCon} {continue}       
  incr n
  set memSum [expr $fmem+$dmem+$gmem+$emem]
  set icSum [expr $fd+$fg+$fe+$dg+$de+$ge]
  set sum [expr $memSum+$icSum]
  puts -nonewline $out [format "%3d%6d%6d%6d%6d |%3d%3d%3d%3d |" \
    $sum $f $d $g $e $fmem $dmem $gmem $emem]  
  puts $out " $fd $fg $fe $dg $de $ge" 
#  puts $out " | $nFD $nFG $nFE $nEF $nED $nEG"
  # Output: <Total contacts><f frame><d frame><g frame><e frame> \
  # <fmem><dmem><gmem><emem>
# Con   f     d     g     e  mem f  d  g  e | fdfgfedgdege 
  flush $out		
}
puts $out "Total tries: $i"
close $out


