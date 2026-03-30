# lck4.tcl: Compare A111, E188, and E199 ITAM substrates binding to Lck (S8 fig).
# vmd -dispdev none on Eire.
cd /home/jcannon/lck
# Change display
display projection Orthographic
display depthcue off
color Display Background white
axes location Off
########################################
# Load Amber format files for lckPep, starting and 100 ns frames
set a111 [mol new lckPep.top type parm7]
mol addfile lckPep.crd type rst7
mol addfile lckPep.eq11.rst type netcdf $a111
# Adjust representations
# Lck NewCartoon
mol delrep 0 $a111
mol color ColorID 5
mol material Transparent
mol selection "resid 1 to 273"
mol representation NewCartoon
mol addrep $a111
# Substrate Licorice
mol color Name
mol material Opaque
mol selection "resid 274 to 288 and mass>2"
mol representation Licorice
mol addrep $a111
# ATP Licorice
mol selection "resid 289"
mol representation Licorice
mol addrep $a111
# MG VDW
mol selection "name MG"
mol representation VDW 0.5 12
mol addrep $a111
# Reposition models by fitting Lck CA
# ^ Lck resid 1-273; A111 substrate 274-288; ATP 289; two MG
set lckA111 [atomselect $a111 "name CA and resid 1 to 273" frame 0]
set n 281; # This is the pTyr of A111.
set t 4; # The 9 residues around pTyr.
set a [expr $n- $t]; # first ITAM residue to fit
set b [expr $n+ $t]; # last ITAM residue to fit
set A111sub [atomselect $a111 "name CA and resid $a to $b" frame 0]
# Second frame
set lckA111b [atomselect $a111 "name CA and resid 1 to 273" frame 1]
set A111 [atomselect $a111 all frame 1]
set A111subb [atomselect $a111 "name CA and resid $a to $b" frame 1]
# ^ Other models will be fit to Lck frame 0 of A111
set M [measure fit $lckA111b $lckA111]
$A111 move $M
## Take renderings here for frames 0,1 comparison
# Rotate around to get good view of active site, get rotate_matrix to use for other models.
set RM [molinfo top get rotate_matrix]
molinfo top set rotate_matrix $RM
set rmsLck [measure rmsd $lckA111b $lckA111]; # Lck CA RMSD after fitting Lck =1.93
# Fit 9-residue substrate
set M [measure fit $A111subb $A111sub]
$A111 move $M
set rmsSub [measure rmsd $A111subb $A111sub]; # 9-residue substrate RMSD =1.37
########################################
# Load Amber format files for E188
set e188 [mol new E188.top type parm7]
mol addfile E188.crd type rst7 $e188
mol addfile E188.eq11.rst type netcdf $e188
# Adjust representations
# Lck NewCartoon
mol delrep 0 $e188
mol color ColorID 5
mol material Transparent
mol selection "resid 1 to 273"
mol representation NewCartoon
mol addrep $e188
# Substrate Licorice
mol color Name
mol material Opaque
mol selection "resid 274 to 288 and mass>2"
mol representation Licorice
mol addrep $e188
# ATP Licorice
mol selection "resid 289"
mol representation Licorice
mol addrep $e188
# MG VDW
mol selection "name MG"
mol representation VDW 0.5 12
mol addrep $e188
# Adjust E188 orientation by fitting to A111 Lck position.
set lckE188 [atomselect $e188 "name CA and resid 1 to 273" frame 0]
set E188 [atomselect $e188 all frame 0]
# Fitting: fit to lckA111 and move E188prot
set M [measure fit $lckE188 $lckA111]
$E188 move $M
set rms [measure rmsd $lckE188 $lckA111]; # Lck CA RMSD after fitting Lck =0
# Rotate to match A111
molinfo top set rotate_matrix $RM
# Render frame 0
# Fit 9-residue substrate
set E188sub [atomselect $e188 "name CA and resid $a to $b" frame 0]
set M [measure fit $A111sub $E188sub]
$E188 move $M
set rmsSub [measure rmsd $E188sub $A111sub]; # 9-residue substrate RMSD =1.61
# frame 1
# Adjust E188 orientation by fitting to A111 Lck position.
set lckE188 [atomselect $e188 "name CA and resid 1 to 273" frame 1]
set E188 [atomselect $e188 all frame 1]
# Fitting: fit to lckA111 and move E188prot
set M [measure fit $lckE188 $lckA111]
$E188 move $M
set rms [measure rmsd $lckE188 $lckA111]; # Lck CA RMSD 0 to 100 ns comparison =1.75
# Fit 9-residue substrate
set E188sub1 [atomselect $e188 "name CA and resid $a to $b" frame 1]
set M [measure fit $E188sub1 $E188sub]
$E188 move $M
set rmsSub [measure rmsd $E188sub1 $E188sub]; # substrate RMSD 0 to 100 ns comparison = 1.69
# Compare E188 substrate to reference A111 substrate
set M [measure fit $E188sub1 $A111sub]
$E188 move $M
set rmsSub [measure rmsd $E188sub1 $A111sub]; # = 1.71
########################################
# Load Amber format files for E199
set e199 [mol new E199.top type parm7]
mol addfile E199.crd type rst7
mol addfile E199.eq11.rst type netcdf $e199
# Adjust representations
# Lck NewCartoon
mol delrep 0 $e199
mol color ColorID 5
mol material Transparent
mol selection "resid 1 to 273"
mol representation NewCartoon
mol addrep $e199
# Substrate Licorice
mol color Name
mol material Opaque
mol selection "resid 274 to 288 and mass>2"
mol representation Licorice
mol addrep $e199
# ATP Licorice
mol selection "resid 289"
mol representation Licorice
mol addrep $e199
# MG VDW
mol selection "name MG"
mol representation VDW 0.5 12
mol addrep $e199
# Adjust E199 orientation by fitting to A111 Lck position.
set lckE199 [atomselect $e199 "name CA and resid 1 to 273" frame 0]
set E199 [atomselect $e199 all frame 0]
# Fitting: fit to lckA111 and move E188prot
set M [measure fit $lckE199 $lckA111]
$E199 move $M
set rms [measure rmsd $lckE199 $lckA111]; # Lck CA RMSD after fitting Lck =0
# Rotate to match A111
molinfo top set rotate_matrix $RM
# Render frame 0
# Fit 9-residue substrate
set E199sub [atomselect $e199 "name CA and resid $a to $b" frame 0]
set M [measure fit $A111sub $E199sub]
$E199 move $M
set rmsSub [measure rmsd $E199sub $A111sub]; # 9-residue substrate RMSD =1.47
# frame 1
# Adjust E199 orientation by fitting to A111 Lck position.
set lckE199 [atomselect $e199 "name CA and resid 1 to 273" frame 1]
set E199 [atomselect $e199 all frame 1]
# Fitting: fit to lckA111 and move E188prot
set M [measure fit $lckE199 $lckA111]
$E199 move $M
set rms [measure rmsd $lckE199 $lckA111]; # Lck CA RMSD after fitting Lck =1.98
# Fit 9-residue substrate
set E199sub1 [atomselect $e199 "name CA and resid $a to $b" frame 1]
set M [measure fit $E199sub1 $E199sub]
$E199 move $M
set rmsSub [measure rmsd $E199sub1 $E199sub]; # 9-residue substrate RMSD =2.45
# Compare E188 substrate to reference A111 substrate
set M [measure fit $E199sub1 $A111sub]
$E199 move $M
set rmsSub [measure rmsd $E199sub1 $A111sub]; # = 2.02
