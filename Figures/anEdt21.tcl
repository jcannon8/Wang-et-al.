# anEdt21.tcl: Lck binding and proximal collisions to edt 200 ns MD (S9 fig).
# Output edtLck2.$mod.dat
# Derived from anEdt20.tcl 
# vmd -dispdev none on Eire.
# Contacts are heavy atoms less than or equal to 1.0 angstroms apart. 
set rmdsMax 3.69;  # Limit for Lck substrate fit.
# ^ Not used here
cd /t/tcr16/edt
###########################################################
# Lck substrate atomselection
set lckpep [mol new /home/jcannon/lck/Lck-Peptide.pdb]
set lck [atomselect $lckpep "chain A"]
set lcksub [atomselect $lckpep all]
set lckH [atomselect $lckpep "chain A and mass>2"]; # Lck heavy atoms
# Choose residues for fitting
set t 4; # The 9 residues around pTyr for fitting.
set w 7; # How far away from pTyr are Lck collisions allowed? 
set c [expr 285- $t]; # first Lck substrate residue to fit
set d [expr 285+ $t]; # last Lck substrate residue to fit
set pSel [atomselect $lckpep "name CA and chain B and resid $c to $d"]
###########################################################
# ITAM list with {chain resid}
set itamsAll {{A 72} {A 83} {A 111} {A 123} {A 142} {A 153} 
{B 72} {B 83} {B 111} {B 123} {B 142} {B 153}  
{D 149} {D 160} {E 188} {E 199} {F 188} {F 199} {G 160} {G 171} }
set itams { {D 149} {D 160} {E 188} {E 199} }
###########################################################
# Load edt 5.4-205 ns MD trajectories for 11 edt models.
set edtModels "84 171 30 115 35 90 13 41 32 153 198"
foreach mod $edtModels {
# Use stripped topology with just 2448 protein atoms that were saved in trajectory.
# made by anEdt1.tcl
set tcr [mol new /t/tcr16/edt.prot.psf]
# The NPT MD started with edt.*.eq2.
for {set i 2} {$i<=53} {incr i} {
  set mdFile /t/tcr16/edt/edt.$mod.eq$i.cdf
  if { [file exists $mdFile] } {
    mol addfile $mdFile type netcdf waitfor all
  }
}
set out [open edtLck2.$mod.dat w]
set lastframe [molinfo $tcr get numframes]
###########################################################
# ITAM atomselection's
set i 0
foreach s $itams {
  set chn [lindex $s 0]; # ITAM chain
  set n [lindex $s 1];   # ITAM resid
  set a [expr $n- $t]; # first ITAM residue to fit
  set b [expr $n+ $t]; # last ITAM residue to fit
  set iSel($i) [atomselect $tcr "chain $chn and name CA and resid $a to $b"]
  # Check proximal iNot (same as ITAM chain) or other protein colissions jNot.
  set e [expr $n- $w];set g [expr $n+ $w]
  set iNot($i) [atomselect $tcr "chain $chn and mass>2 and not (resid $e to $g)"]
  set jNot($i) [atomselect $tcr "protein and mass>2 and not chain $chn"]
  # Each subunit will collect a list of contacts per residue
  set chnSel [atomselect $tcr "name CA and chain $chn"]
  set resN($i) [$chnSel get resid]
  $chnSel delete
  set lenRes [llength $resN($i)]; # Zero fill list 
  for {set x 0} {$x<$lenRes} {incr x} {lappend conP($i) 0}
  incr i
}
# Go through every 10 ps frame.
for {set f 0} {$f<$lastframe} {incr f} {
  for {set k 0} {$k<$i} {incr k} {
    $iSel($k) frame $f; $iNot($k) frame $f; $jNot($k) frame $f
    # Fitting: fit substrate to ITAM and move lcksub
    set M [measure fit $pSel $iSel($k)]
    $lcksub move $M
    set rms [measure rmsd $pSel $iSel($k)]
    # Check proximal Lck collisions
    set conSub1 [measure contacts 1.0 $lckH $iNot($k)]
    set prox [llen [lindex $conSub1 1]]
    if {$prox==0} continue 
    # Get unique resid of proximal contacts   
    set conInd [lindex $conSub1 1]
    set cons [atomselect $tcr "index $conInd"]
    set conRes [lsort -unique [$cons get resid]]
    foreach r $conRes {
      set iRes [lsearch $resN($k) $r]; # index for resid in resN
      set oCon [lindex $conP($k) $iRes]; # original contact number
      incr oCon; # increment and save new contact number
      lset conP($k) $iRes $oCon   
    }
    $cons delete
  }
  if {[expr $f % 10]==0} {puts "$f edt.$mod"}; # progress report
}
# Output <ITAM chn resid><resid's><collisions> for each ITAM Tyr, followed by frame numbers
for {set k 0} {$k<$i} {incr k} {
  puts $out [lindex $itams $k]
  puts $out $resN($k)
  puts $out $conP($k)
}
puts $out $lastframe;  # output frame numbers for model
close $out
mol delete $tcr
for {set k 0} {$k<$i} {incr k} {
  # Delete atomselections
  $iSel($k) delete; $iNot($k) delete; $jNot($k) delete; 
  # Empty collision vectors, conP($i)
  unset conP($k)
}
}
###########################################################
# Analysis for frequency of Lck proximal collisions per residue.
# Get edtLck2.$mod.dat data for all models and
# output collision data for each model in columns (sum in last).
tclsh<<"eof"
# Load data for all models into line($nMod).
set edtModels "84 171 30 115 35 90 13 41 32 153 198"
set nMod 0
foreach mod $edtModels {
  if {[file exist edtLck2.$mod.dat]==1} {
    set in [open edtLck2.$mod.dat r]
    set inData [read $in]
    close $in
    set line($nMod) [split $inData "\n"]
    incr nMod
  } 
}
# Use first edtLck2.$mod.dat to get resid's, ITAM names, and open output files.
for {set i 0; set x 0} {$i <4} {incr i; set x [expr $x+3]} {
  # i=index to ITAM output file; x=index to data in input file.
  set itam($i) [lindex $line(0) $x]
  # Remove space from ITAM name for output file names.
  set iName [string map {" " ""} $itam($i)]
  set out($i) [open edtRes.$iName.dat w]
  # Get resid's
  set resN($i) [lindex $line(0) [expr $x+1]]     
}
# Save resN(), collision data, and close output files
for {set i 0} {$i <4} {incr i} {
  # i=index to ITAM output file
  set len [llength $resN($i)]; # Number of output lines
  for {set j 0} {$j<$len} {incr j} {
    set sum 0
    # j=output line number, per line; first output resid
    puts -nonewline $out($i) [format "%4d" [lindex $resN($i) $j]]
    # Then append the line with collision data from each edtLck2.$mod.dat file 
    set x [expr $i*3 +2]; # x=index to data in input file. 
    for {set k 0} {$k < $nMod} {incr k} {
      if {[llength $line($k)]==0} break; # End line if file not parsed.
      puts -nonewline $out($i) [format "%7d" [lindex $line($k) $x $j]]
      set sum [expr $sum + [lindex $line($k) $x $j]]
    }
    # Terminate line with line sum
    puts $out($i) [format "%8d" $sum]  
  } 
  close $out($i)
} 
eof
# This places data columns together like bash paste by using tcl puts.
###########################################################
cd /t/tcr16/edtAn
# Plot edt Lck proximal collision averages for ensemble
# S9 fig panel A
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 1000,750
set out 'edt.lck5.jpg'
# layout: rows, columns
set multiplot layout 2,2 title "{/:Bold Proximal Lck collisions with CD3εδ}"
set border 3
set lmargin 4
set rmargin 2
set xtics out nomirror 
set ytics out nomirror offset 1,0
set grid ytics
set xlabel "{/:Bold Residue}"  
#set ylabel "{/:Bold Collision frequency (%)}" offset 2,0
set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
# These use ps 7 circle, red for ITAM, blue for other, magenta for TM.
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "magenta" lt 1 lw 2 pt 7 pi -1 ps 1.0
# Use ITAM specific label and file name. Subunit specific xrange
# These trajectories have 20,500 frames 
set label "{/:Bold CD3δ Tyr149}" at graph 0.05,1.0
# CD3δ: 171 residues, TM 100-129
set xrange [*:171]
set xtics 90,20,170 
set xtics add ("" 100, "" 120, "" 140, "" 160)
plot "edtRes.D149.dat" using ($1<143||$1>155?$1:1/0):($13/(11*205)) with linespoints ls 2 t "",\
'' using ($1>141&&$1<157?$1:1/0):($13/(11*205)) with linespoints ls 1 t "",\
'' using ($1<129?$1:1/0):($13/(11*205)) with linespoints ls 3 t ""
unset label
set label "{/:Bold CD3δ Tyr160}" at graph 0.05,1.0
plot "edtRes.D160.dat" using ($1<154||$1>166?$1:1/0):($13/(11*205)) with linespoints ls 2 t "",\
'' using ($1>152&&$1<168?$1:1/0):($13/(11*205)) with linespoints ls 1 t "",\
'' using ($1<129?$1:1/0):($13/(11*205)) with linespoints ls 3 t ""
unset label
set xrange [*:207]
set xtics 120,20,200 
set xtics add ("" 130, "" 150, "" 170, "" 190)
set label "{/:Bold CD3ε Tyr188}" at graph 0.05,1.0
# CD3ε: 207 residues, TM 127-155 
plot "edtRes.E188.dat" using ($1<182||$1>194?$1:1/0):($13/(11*205)) with linespoints ls 2 t "",\
'' using ($1>180&&$1<196?$1:1/0):($13/(11*205)) with linespoints ls 1 t "",\
'' using ($1<155?$1:1/0):($13/(11*205)) with linespoints ls 3 t ""
unset label
set label "{/:Bold CD3ε Tyr199}" at graph 0.05,1.0
plot "edtRes.E199.dat" using ($1<193||$1>205?$1:1/0):($13/(11*205)) with linespoints ls 2 t "",\
'' using ($1>191&&$1<207?$1:1/0):($13/(11*205)) with linespoints ls 1 t "",\
'' using ($1<155?$1:1/0):($13/(11*205)) with linespoints ls 3 t ""  
quit
eof
