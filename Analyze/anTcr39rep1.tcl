# anTcr39rep1.tcl: Report Nck, Lck binding, and ITAM Tyr penetration in tcr39b models.
# Use data from 
# tcr39.*.**.Nckt3.dat, tcr39lck2.dat, and tcr39.memy85b.dat. 
# tcr39.*.**.Nckt3.dat files from anTcr39r.tcl, one line per model.
# tcr39lck2.dat from Runtcr39lck2.sh and tcr39lck2.tcl, one line per Lck binding.
# ITAM TYR penetration tcr39.memy85b.dat from anTcr38Y3.tcl
# Use tclsh because no VMD extensions  necessary.
cd /t/tcr39an
set out [open tcr39b85.dat w]
# anTcr39r.tcl output are 62 tcr39.*.**.Nckt3.dat files.
# Example *.Nckt3.dat output:
#     76   0  24    3.10   1.08  -8.48   0.57 |  125.67  129.85 |  131.20  132.49
# <eq frame><eNck binding><fNck binding><eMem1><eMem2><fMem1><fMem2> |
# Columns 1 & 2 have Nck-PRS reference # or -1 if no Nck binding without protein collision.
# Columns 3 & 5 have Nck distance to average membrane.
# Columns 3 & 6 have shortest Nck distance to membrane.
for {set j 0} {$j<=11} {incr j} {
foreach k {00 01 02 10 11 12 20 21 22} {
set fname tcr39.$j.$k 
# Skip missing  models.
if {[file exist /t/tcr39b/$fname.psf]==0} continue
# Get Nck binding from *.Nckt3.dat files.
set in [open $fname.Nckt3.dat r]
set inData1 [read $in]
close $in
foreach line [split $inData1 "\n"] {
  if {$line == ""} break
  set tm [lindex $line 0];    # eq frame (5 ns).
  if {$tm!=85} continue;  # Only consider last, 425 ns eq85 frame
  set eProt [lindex $line 1]; # eNck protein collision 
  set fProt [lindex $line 2]; # fNck protein collision
  set eMem2 [lindex $line 4]; # eNck to membrane distance
  set fMem2 [lindex $line 6]; # fNck to membrane distance
  # CD3εδ PRS angles
  set eP182 [lindex $line 8] 
  set eP184 [lindex $line 9] 
  # CD3εγ PRS angles
  set fP182 [lindex $line 11]   
  set fP184 [lindex $line 12] 
  # Nck can bind if no collision with other proteins, 
  # no Nck collision with membrane, and no bad angles.
  if {$eProt>-1 && $eProt<24 && $eMem2>1.0 && $eMem2 !=99 && $eP182>100 && $eP184>100} {
    puts -nonewline $out [format "%2d.%s | %2d  %4.2f |" $j $k $eProt $eMem2]
  } else {
    puts -nonewline $out [format "%2d.%s |  no bind |" $j $k ]
  }
  if {$fProt>-1 && $fProt<24 && $fMem2>1.0 && $fMem2 !=99 && $fP182>100 && $fP184>100} {
    puts -nonewline $out [format " %2d %5.2f |" $fProt $fMem2]
  } else {
    puts -nonewline $out [format "  no bind |"]
  }
  # Get potential Lck binding from tcr39lck2.dat.
  # tcr39lck2.dat has Lck binding and no model binds more than one,
  # at least in the 425 ns frame, eq85.
  # See if this model has any Lck binding.
  set lck1 [exec awk {{if ($2==85) {print $0}}} /t/tcr39an/tcr39lck2.dat]
  set bind [lsearch $lck1 $fname]; # returns -1 if no binding
  if {$bind==-1} {puts -nonewline $out "    no bind "} else {
    # Model binds Lck, get substrate and Lck distance to membrane.
    set lck0 [exec awk {{if ($2==85) {print $0}}} /t/tcr39an/tcr39lck2.dat | grep $fname]
    # Example lck0:
    #   tcr39.0.21  85 B  142  3.24   15.01   11.42
    # Column 2 has chain, column 3 has residue, column 6 has Lck to membrane.
    set chn [lindex $lck0 2]
    set res [lindex $lck0 3]
    set lMem [lindex $lck0 6]
    puts -nonewline $out [format " %s%s %5.2f " $chn $res $lMem]
  }
}    
  # Next ITAM Tyr membrane penetration from [tcr38,tcr39].memy85b.dat
  set in [open /t/tcr38an/tcr39.memy85b.dat r]
  set inData2 [read $in]
  close $in
  foreach line [split $inData2 "\n"] {
    if {$line == ""} break
    set mod [lindex $line 0];   # TCR model
    set pair [lindex $line 1]; # ITAM substrate-lipid resname
    set atom [lindex $line 2]; # lipid atom name
    if {[string compare $mod $fname] ==0} {puts -nonewline $out "| $pair $atom "}   
  }
  puts $out " "
}
}
close $out  
###########################################################
# Example output:
#<model>|<eRef><eNck mem>|<fRef><fNck mem>|<Lck bind><Lck mem>|
# <ITAM TYR chain+residue-lipid atom>
# 0.00 |  no bind |  no bind | B153 15.39 | B142-DOPC O13 | E188-POPC O31
# 0.01 |  0  1.10 |  no bind |    no bind | F199-POPS O13A | G160-DOPC O14
