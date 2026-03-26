# makeTcr38b.tcl: Add four CD3 cytoplasmic tails to tcr38 to make tcr38.[0-8].
# and extend chain n to end by four residues.
# The resulting models are called tcr38.[0-8].
# Quads of reference numbers and frame numbers come from environment variables
# model and frames. Use best framesets from makeTcr38a.tcl.
#
# This version is used as a VMD script as in 
# export model=0
# export frames="8778   413  2820  1310"
# vmd -e makeTcr38b.tcl or vmd -dispdev none 
set model $env(model)
set frames $env(frames)
proc base3 {m} {
  # Convert m into a base 3 number
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
# =27*B3+9*C3+3*D3+E3 in tcr15s.xlsx
set outFile tcr38.$model
# A sextet of reference numbers {f, d, g, e, a, b} for addTail14b.tcl input.
puts "$fRef $dRef $gRef $eRef"
# Parse frames from frameset
set fFrame [lindex $frames 0]; set dFrame [lindex $frames 1]
set gFrame [lindex $frames 2]; set eFrame [lindex $frames 3]
# tcr38: removeExt4.tcl removed extracellular protein and water from tcr35c.eq70a.rst (322 ns)
set tcr [mol new tcr38.psf]
mol addfile tcr38.pdb
# No Nck addition for tcr38 models
###########################################################
# Add eps CT to chain f.
set cd3f [mol new /home/jcannon/tcr15/eps.top type parm7]
mol addfile /home/jcannon/tcr15/eps.c$fRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/eps.c$fRef.eq3.cdf type netcdf waitfor all
# Fit between tcr and eps uses SER LYS ASN ARG residues in each.
# Align resid 2-5 of cd3f and resid 152 to 155 of tcr15.
set End [atomselect $tcr "chain F and name CA and resid 152 to 155"]
set Start [atomselect $cd3f "resid 2 to 5 and name CA" frame $fFrame]
set cd3CD [atomselect $cd3f all frame $fFrame]
# RMS fit Start to End
set M [measure fit $Start $End]
# Move cd3ex to fit tcr 
$cd3CD move $M
# Check alignment
set fit [measure rmsd $Start $End]
# Cross over in the center of the shared residues instead of the ends.
set cd3fN [atomselect $tcr "chain F and resid 123 to 153"]
set cd3fC [atomselect $cd3f "resid 4 to 57" frame $fFrame]
$cd3fN set chain f
$cd3fC set chain f
# Residues in CT need to continue the numbering from the N-terminus
foreach i [[atomselect $cd3f "resid 4 to 57"] get index] {
  set r [[atomselect $cd3f "index $i"] get resid]
  [atomselect $cd3f "index $i"] set resid [expr $r + 150]
}
# Use bash commands to fuse N- and C-termini
$cd3fN writepdb cd3fN.pdb
$cd3fC writepdb cd3fC.pdb
exec cat cd3fN.pdb > fChain.pdb
exec cat cd3fC.pdb >> fChain.pdb
# Remove END and CRYST records.
exec sed -i /END/d fChain.pdb
exec sed -i /CRYST/d fChain.pdb
file delete cd3fN.pdb; file delete cd3fC.pdb;
###########################################################
# Remove cytoplasmic water that collides with new chain f CT.
set water [atomselect $tcr water]
# Set beta=1 for water to keep.
$water set beta 1
# Test fChain collision with water
set fChain [mol new fChain.pdb]
set fProt [atomselect $fChain all]
# measure contacts works for atoms in different molecules
set con [measure contacts 1.0 $fProt $water]
set wats [lindex $con 1]
# This is more elegant than the proc tagWater used in addTail7.tcl et al.
# Use atomselection "same <keyword> as <selection>"
[atomselect $tcr "same residue as index $wats"] set beta 8
[atomselect $tcr "beta 8"] num
# 528 water atoms are tagged for removal by chain f.
###########################################################
# Add eps CT to chain e.
set cd3e [mol new /home/jcannon/tcr15/eps.top type parm7]
mol addfile /home/jcannon/tcr15/eps.c$fRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/eps.c$fRef.eq3.cdf type netcdf waitfor all
# Align resid 2-5 of cd3ex and resid 152 to 155 of tcr15
# Fit between tcr and eps uses SER LYS ASN ARG residues in each.
set End [atomselect $tcr "chain E and name CA and resid 152 to 155"]
set Start [atomselect $cd3e "resid 2 to 5 and name CA" frame $eFrame]
set cd3CD [atomselect $cd3e all frame $eFrame]
# RMS fit Start to End
set M [measure fit $Start $End]
# Move cd3ex to fit tcr 
$cd3CD move $M
# Check alignment
set fit [measure rmsd $Start $End]
# Cross over in the center of the shared residues instead of the ends.
set cd3fN [atomselect $tcr "chain E and resid 125 to 153"]
set cd3fC [atomselect $cd3e "resid 4 to 57" frame $eFrame]
$cd3fN set chain e
$cd3fC set chain e
# Residues in CT need to continue the numbering from the N-terminus
foreach i [[atomselect $cd3e "resid 4 to 57"] get index] {
  set r [[atomselect $cd3e "index $i"] get resid]
  [atomselect $cd3e "index $i"] set resid [expr $r + 150]
}
$cd3fN writepdb cd3fN.pdb
$cd3fC writepdb cd3fC.pdb
exec cat cd3fN.pdb > eChain.pdb
exec cat cd3fC.pdb >> eChain.pdb
# Remove END and CRYST records.
exec sed -i /END/d eChain.pdb
exec sed -i /CRYST/d eChain.pdb
file delete cd3fN.pdb; file delete cd3fC.pdb;
###########################################################
# Test eChain collision with water
set eChain [mol new eChain.pdb]
set eProt [atomselect $eChain all]
# measure contacts works for atoms in different molecules
set con [measure contacts 1.0 $eProt $water]
set wats [lindex $con 1]
if {$wats>0} {[atomselect $tcr "same residue as index $wats"] set beta 5}
[atomselect $tcr "beta 5"] num
# 510 water atoms are tagged for removal by chain e.
###########################################################
# Add delta3 CT to chain d in tcr model.
set delta [mol new /home/jcannon/tcr15/delta.top type parm7]
mol addfile /home/jcannon/tcr15/delta3.c$dRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/delta3.c$dRef.eq3.cdf type netcdf waitfor all
# delta3 has delta residues 100-171
# Align resid 1-4 of delta and resid 126 to 129 of tcr15
# Fit between tcr and delta uses ALA GLY HIS GLU residues in each.
set End [atomselect $tcr "chain D and name CA and resid 126 to 129"]
set Start [atomselect $delta "resid 1 to 4 and name CA" frame $dFrame]
set cd3CD [atomselect $delta all frame $dFrame]
# RMS fit Start to End
set M [measure fit $Start $End]
# Move cd3ex to fit tcr 
$cd3CD move $M
# Check alignment
set fit [measure rmsd $Start $End]
# Cross over in the center of the shared residues instead of the ends.
set cd3fN [atomselect $tcr "chain D and resid 98 to 126"]
set cd3fC [atomselect $delta "resid 2 to 46" frame $dFrame]
$cd3fN set chain d
$cd3fC set chain d
# Residues in CT need to continue the numbering from the N-terminus
foreach i [[atomselect $delta "resid 2 to 46"] get index] {
  set r [[atomselect $delta "index $i"] get resid]
  [atomselect $delta "index $i"] set resid [expr $r + 125]
}
# Use bash commands to fuse N- and C-termini
$cd3fN writepdb cd3fN.pdb
$cd3fC writepdb cd3fC.pdb
exec cat cd3fN.pdb > dChain.pdb
exec cat cd3fC.pdb >> dChain.pdb
# Remove END and CRYST records.
exec sed -i /END/d dChain.pdb
exec sed -i /CRYST/d dChain.pdb
file delete cd3fN.pdb; file delete cd3fC.pdb;
###########################################################
# Test dChain collision with water
set dChain [mol new dChain.pdb]
set dProt [atomselect $dChain all]
# measure contacts works for atoms in different molecules
set con [measure contacts 1.0 $dProt $water]
set wats [lindex $con 1]
[atomselect $tcr "same residue as index $wats"] set beta 7
[atomselect $tcr "beta 7"] num
# 381 water atoms are tagged for removal by chain d.
###########################################################
# Add gamma3 CT to chain g in tcr model.
set gamma [mol new /home/jcannon/tcr15/gamma.top type parm7]
mol addfile /home/jcannon/tcr15/gamma3.c$gRef.eq2.mdcrd type crd waitfor all
mol addfile /home/jcannon/tcr15/gamma3.c$gRef.eq3.cdf type netcdf waitfor all
# gamma3 has gamma residues 111-182
# Align resid 1-4 of gamma and resid 135 to 138 of tcr15
# Fit between tcr and gamma uses PHE ILE ALA GLY residues in each.
set End [atomselect $tcr "chain G and name CA and resid 135 to 138"]
set Start [atomselect $gamma "resid 1 to 4 and name CA" frame $gFrame]
set cd3CD [atomselect $gamma all frame $gFrame]
# RMS fit Start to End
set M [measure fit $Start $End]
# Move cd3ex to fit tcr 
$cd3CD move $M
# Check alignment
set fit [measure rmsd $Start $End]
# Cross over in the center of the shared residues instead of the ends.
set cd3fN [atomselect $tcr "chain G and resid 108 to 135"]
set cd3fC [atomselect $gamma "resid 2 to 48" frame $gFrame]
$cd3fN set chain g
$cd3fC set chain g
# Residues in CT need to continue the numbering from the N-terminus
foreach i [[atomselect $gamma "resid 2 to 48"] get index] {
  set r [[atomselect $gamma "index $i"] get resid]
  [atomselect $gamma "index $i"] set resid [expr $r + 134]
}
# Use bash commands to fuse N- and C-termini
$cd3fN writepdb cd3fN.pdb
$cd3fC writepdb cd3fC.pdb
exec cat cd3fN.pdb > gChain.pdb
exec cat cd3fC.pdb >> gChain.pdb
# Remove END and CRYST records.
exec sed -i /END/d gChain.pdb
exec sed -i /CRYST/d gChain.pdb
file delete cd3fN.pdb; file delete cd3fC.pdb;
###########################################################
# Test gChain collision with water
set gChain [mol new gChain.pdb]
set gProt [atomselect $gChain all]
# measure contacts works for atoms in different molecules
set con [measure contacts 1.0 $gProt $water]
set wats [lindex $con 1]
[atomselect $tcr "same residue as index $wats"] set beta 6
[atomselect $tcr "beta 6"] num
# 414 water atoms are tagged for removal by chain g.
###########################################################
# Save chains ABM and structural cholesterol and DOPC
# Here just save the chain A, B TMs.
[atomselect $tcr "chain A"] writepdb aChain.pdb
[atomselect $tcr "chain B"] writepdb bChain.pdb
[atomselect $tcr "chain M"] writepdb mChain.pdb
# Save S1 structural cholesterol and S2 DOPC. 
set S1 [atomselect $tcr "chain S and resname CHL1"] 
$S1 set segid S1; $S1 writepdb S1.pdb
set S2 [atomselect $tcr "chain S and resname DOPC"]
$S2 set segid S2; $S2 writepdb S2.pdb
###########################################################
# Extend chain N by 4 residues
if {[file exists nChain3.pdb]==0} {
# Use chain n extended by four C-terminal residues by extChainN.sh, nChain2.pdb
set nExt [mol new /home/jcannon/tcr15/nChain2.pdb]
set nChain [atomselect $nExt all]
# Align resid 37 to 40 of nExt with chain N resid 305 to 308 of tcr38
set End [atomselect $tcr "chain N and name CA and resid 305 to 308"]
set Start [atomselect $nExt "name CA and resid 37 to 40"]
# RMS fit Start to End
set M [measure fit $Start $End]
# Move cd3ex to fit tcr 
$nChain move $M	
# Check alignment
set fit [measure rmsd $Start $End]
# Cross over in the center of the shared residues instead of the ends.
set cd3fN [atomselect $tcr "chain N and resid 269 to 306"]
set cd3fC [atomselect $nExt "resid 39 to 44"]
#$cd3fN set chain n
#$cd3fC set chain n
# Residues in nChain extension need renumbering 
foreach i [[atomselect $nExt "resid 39 to 44"] get index] {
  set r [[atomselect $nExt "index $i"] get resid]
  [atomselect $nExt "index $i"] set resid [expr $r + 268]
}
# Use bash commands to fuse N- and C-termini
$cd3fN writepdb cd3fN.pdb
$cd3fC writepdb cd3fC.pdb
exec cat cd3fN.pdb > nChain3.pdb
exec cat cd3fC.pdb >> nChain3.pdb
# Remove END and CRYST records.
exec sed -i /END/d nChain3.pdb
exec sed -i /CRYST/d nChain3.pdb
file delete cd3fN.pdb; file delete cd3fC.pdb;
}
##############################
# Save lipids and ions
set lipids [atomselect $tcr "chain J"]
$lipids set chain j
$lipids writepdb lipids.pdb
set ions [atomselect $tcr "name SOD CLA"]
$ions set chain i
#$ions writepdb ions.pdb
# Initially save all ions. That gave netCharge=+34.
[atomselect $tcr "name SOD"] num; # 556
[atomselect $tcr "name CLA"] num; # 296
# Remove 34 SOD to neutralize
set resSOD [[atomselect $tcr "name SOD"] get resid]
set saveSOD [lrange $resSOD 0 end-34]
# Now save fewer SOD
[atomselect $tcr "name CLA or (name SOD and resid $saveSOD)"] writepdb ions.pdb
#### Water
# Saved waters already marked above.
# The waters need to be saved in smaller chunks for psfgen.
set saveWat [atomselect $tcr "water and beta 1"]
$saveWat set chain z
set Res [lsort -unique -integer [$saveWat get residue]]
set last [lrange $Res end end]
# A PDB residue can only have 4-characters
# Save water in 9999 residue chunks
set chunk 9999
for {set n 1; set k 0} {$k<=$last} {incr n} {
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
}
# water1.pdb to water[$n-1].pdb saved
###########################################################
# Use psfgen to build tcr38.$model.[psf,pdb]
# https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf
# http://www.ks.uiuc.edu/Training/Tutorials/namd/namd-tutorial-html/node6.html
package require psfgen
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
# The CD3 CTs were made with Amber HIE -> Charmm HSE
# However, change to delta protonation so all protein protonated the same.
pdbalias residue HIE HSD
# Load segments with coordinates.
# This defines the order in this model (same as 7FJD and tcr30 order)
segment a {pdb aChain.pdb}; coordpdb aChain.pdb a
segment b {pdb bChain.pdb}; coordpdb bChain.pdb b
segment d {pdb dChain.pdb}; coordpdb dChain.pdb d
segment e {pdb eChain.pdb}; coordpdb eChain.pdb e
segment f {pdb fChain.pdb}; coordpdb fChain.pdb f
segment g {pdb gChain.pdb}; coordpdb gChain.pdb g
segment m {pdb mChain.pdb}; coordpdb mChain.pdb m
segment n {pdb nChain3.pdb}; coordpdb nChain3.pdb n
# Structural cholesterol CHL1 201, DOPC 190, membrane, ions
segment s1 {pdb S1.pdb}; coordpdb S1.pdb s1
segment s2 {pdb S2.pdb}; coordpdb S2.pdb s2
guesscoord
# The hydrogen coordinates were poorly guessed above
# These hydrogen warnings could be resolved by eliminating them when pdb chains are saved.
segment j {pdb lipids.pdb}; coordpdb lipids.pdb j
segment k {pdb ions.pdb}; coordpdb ions.pdb k
# Load in water
pdbalias residue WAT TIP3
pdbalias atom WAT O OH2
# Water is in multiple pdb files
for {set i 1} {$i<$n} {incr i} {
  segment w$i {pdb water$i.pdb}; coordpdb water$i.pdb w$i
  puts $i
}
# Only one disulfide remains from 6jxr.pdb, between zeta chains.
#patch DISU A:32 B:32
# TODO Check on charge and ion balance.
writepsf $outFile.psf
writepdb $outFile.pdb
###########################################################
# Clean up
file delete water.pdb; file delete ions.pdb; file delete lipids.pdb;
file delete aChain.pdb; file delete bChain.pdb; 
file delete dChain.pdb;
file delete eChain.pdb; file delete fChain.pdb; file delete gChain.pdb;
file delete mChain.pdb; file delete S1.pdb; file delete S2.pdb;
foreach f [glob water*.pdb] {file delete $f}
###########################################################
# Testing
# Check charge to see if ions need alteration.
set tcr3 [mol new $outFile.psf]
mol addfile $outFile.pdb waitfor all
set sel [atomselect $tcr3 all]
measure sumweights $sel weight charge
# ^The charge was originally +34, now zero.
################################ 
# Check for loop penetration
proc checkThread {molA resN loopName} {
  # Check loop penetration: 
  # molA=molid; resN=residue numbers; loopName=loop atom names
  set OK 0; # Start with OK=0, no penetration.
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
        if {$d<1.0} {
          set OK 1
          # Output penetration information
          puts -nonewline [lindex [$loop get {chain resname resid}] 0]
          puts -nonewline " penetrated by "
          puts [lindex [$left get {chain resname resid}] 0]
        }
        $left delete
      }
    }
    $others delete; $this delete; $loop delete
  }  
#  return $OK 
}
# Set the CD3 CTs beta 1 to screen for loop penetration
set all3 [atomselect $tcr3 all]; $all3 set beta 0
set eCT3 [atomselect $tcr3 "chain E and resid 156 to 207"]
set dCT3 [atomselect $tcr3 "chain D and resid 130 to 171"]
set gCT3 [atomselect $tcr3 "chain G and resid 139 to 182"]
set fCT3 [atomselect $tcr3 "chain F and resid 156 to 207"]
set nCT3 [atomselect $tcr3 "chain N and resid 309 to 312"]
$eCT3 set beta 1; $dCT3 set beta 1;$gCT3 set beta 1; $fCT3 set beta 1
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
# Check loop penetrations
puts "Passed the loop penetration test?"
checkThread $tcr3 $pro "N CA CB CD CG"
checkThread $tcr3 $tyr "CZ CE1 CD1 CG CD2 CE2"
checkThread $tcr3 $his "CG ND1 CE1 NE2 CD2"
checkThread $tcr3 $trp "CG CD1 NE1 CE2 CD2"
checkThread $tcr3 $trp "CE2 CD2 CE3 CZ3 CH2 CZ2"
##########################
# Get zero-based protein residue and index numbers. 
# Change to one-based, Amber numbers for restraining during initial MD.
lindex [[atomselect $tcr3 "name CA and protein"] get residue] end; # max zero-based residue=465 
lindex [[atomselect $tcr3 protein] get index] end; # max zero-based index=7455
lindex [[atomselect $tcr3 "not water and not ion"] get residue] end; # residues=2067
lindex [[atomselect $tcr3 "protein or lipid"] get residue] end; # residues=2067
# VMD viewing
display projection Orthographic
display depthcue off
color Display Background white
rotate x by 90
# CD3delta: chain d red
mol color ColorID 1
mol selection chain D
mol addrep $tcr3
# CD3epsilon: chain f grey
mol color ColorID 2
mol selection chain F
mol addrep $tcr3
# CD3gamma: chain g orange
mol color ColorID 3
mol selection chain G
mol addrep $tcr3
# TCRalpha: chain m pink
mol color ColorID 9
mol selection chain M
mol addrep $tcr3
# TCRbeta: chain n yellow3
mol color ColorID 18
mol selection chain N
mol addrep $tcr3
# CD3epsilon: chain e blue
mol color ColorID 0
mol selection chain E
mol addrep $tcr3

