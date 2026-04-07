# anTcr38k.tcl: PLot tcr38.*.** & tcr39.*.** lipid distribution (Fig 18).
# Single 405 ns (eq85) frame analyzed.
# vmd -dispdev none
cd /t/tcr38an
###########################################################
set dist 20; # Distance to count lipids
set m tcr39; # =tcr38 for tcr38b; =tcr39 for tcr39b
set out1 [open ${m}k85a.dat w]; # lipid distribution near protein
set out2 [open ${m}k85b.dat w]; # lipid distribution near CTs
#set out3 [open ${m}k85c.dat w]; # lipid distribution 20 near anchors
set out3 [open ${m}k85d.dat w]; # lipid distribution 40 near anchors
for {set j 0} {$j<=11} {incr j} {
foreach k {00 01 02 10 11 12 20 21 22} {
set fName $m.$j.$k 
# Skip missing  models. These are 405 ns frames.
if {[file exist /t/${m}b/$fName.psf]==0} continue
if {[file exist /t/${m}b/$fName.eq85.rst ]==0} continue
set tcr [mol new /t/${m}b/$fName.psf]
mol addfile /t/${m}b/$fName.eq85.rst type netcdf waitfor all
# Study atomselection of cytoplasmic side lipids
set mem [atomselect $tcr "mass>2 and chain J"]; 
set memMin [lindex [measure minmax $mem] 0 2]; # minimum Z-dimension of membrane.
set memMax [lindex [measure minmax $mem] 1 2]; # maximum Z-dimension of membrane.
set aveMem [expr ($memMin+$memMax)/2]
if {0} {
# Below studies how to specifically select cytoplasmic lipids. 
set cytMem [atomselect $tcr "chain J and mass>2 and z>$aveMem"]
set resCyt [lsort -unique [$cytMem get resid]]; # resid of cytoplasmic lipids
set nameCyt [lsort -unique [$cytMem get resname]]; # Seven lipids in membrane
llength $resCyt; # 937 cytoplasmic lipids, one-half the lipids would be 800
[atomselect $tcr "chain J and name P"] num; # 880 total phospholipids
[atomselect $tcr "chain J and name P and z>$aveMem"] num; # 440 cytoplasmic phospholipids
set cytChl [atomselect $tcr "chain J and resname CHL1 and z>$aveMem"]
set resCytChl [lsort -unique [$cytChl get resid]]
llength $resCytChl; # 432 cytoplasmic cholesterols
set allChl [atomselect $tcr "chain J and resname CHL1"]
set resAllChl [lsort -unique [$allChl get resid]]
llength $resAllChl; # 720 total cholesterols
set extChl [atomselect $tcr "chain J and resname CHL1 and z<$aveMem"]
set resExtChl [lsort -unique [$extChl get resid]]
llength $resExtChl; # 550 extracellular cholesterols
# ^ Those data show that some CHL1 heavy atoms cross the average membrane Z-dimension.
set nameAllChl [lsort -unique [$allChl get name]]; # CHL1 atom names, O3 is shallowest
[atomselect $tcr "chain J and z>$aveMem and name P"] num; # 440 phospholipids
[atomselect $tcr "chain J and resname CHL1 and z>$aveMem and name O3"] num; # 360 CHL
}
# Make atomselection of cytoplasmic side lipids
set allMem [atomselect $tcr "chain J"]
$allMem set beta 0
# https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node86.html#SECTION00838200000000000000
# Set heavy cytoplasmic atoms to beta=1. Find phospholipids via P, CHL1 via O3 atoms.
set cytPL [atomselect $tcr "mass>2 and same residue as (chain J and z>$aveMem and name P)"]
$cytPL set beta 1
set CH "mass>2 and same residue as (chain J and resname CHL1 and z>$aveMem and name O3)"
set cytCH [atomselect $tcr $CH]
$cytCH set beta 1
#############
# Get lipid atoms within dist of any TCR protein
set nearL [atomselect $tcr "beta 1 and within $dist of protein"]
set resNear [lsort -unique [$nearL get residue]]
llength $resNear; # 210 residues within 20 of any TCR protein
set tot 0
foreach r {CHL1 POPC DPPC DOPC POPS DOPE POPI} {
  set res [atomselect $tcr "residue $resNear and resname $r"] 
  set n($r) [llength [lsort -unique [$res get residue]]]
  $res delete
  set tot [expr $tot + $n($r) ]
}
puts -nonewline $out1 [format "%6s" $fName]
foreach r {CHL1 POPC DPPC DOPC POPS DOPE POPI} {
puts -nonewline $out1 [format "%4d" $n($r)]
}
puts $out1 [format "%4d" $tot]
$allMem delete; $cytPL delete; $cytCH delete; $nearL delete
#############
# Get lipid atoms within dist of CD3 CTs
set eCT [atomselect $tcr "chain E and mass>2 and resid>155"]
set fCT [atomselect $tcr "chain F and mass>2 and resid>155"]
set dCT [atomselect $tcr "chain D and mass>2 and resid>129"]
set gCT [atomselect $tcr "chain G and mass>2 and resid>138"]
set aCT [atomselect $tcr "chain A and mass>2 and resid>57"]
set bCT [atomselect $tcr "chain B and mass>2 and resid>57"]
$eCT set beta 2
$fCT set beta 2
$dCT set beta 2
$gCT set beta 2
$aCT set beta 2
$bCT set beta 2
set nearL2 [atomselect $tcr "beta 1 and within $dist of beta 2"]
set resNear2 [lsort -unique [$nearL2 get residue]]
llength $resNear2; # 206 residues within 20 of CD3 CTs
set tot 0
foreach r {CHL1 POPC DPPC DOPC POPS DOPE POPI} {
  set res [atomselect $tcr "residue $resNear2 and resname $r"] 
  set n($r) [llength [lsort -unique [$res get residue]]]
  $res delete
  set tot [expr $tot + $n($r) ]
}
puts -nonewline $out2 [format "%6s" $fName]
foreach r {CHL1 POPC DPPC DOPC POPS DOPE POPI} {
puts -nonewline $out2 [format "%4d" $n($r)]
}
puts $out2 [format "%4d" $tot]
$eCT delete; $fCT delete; $dCT delete; $gCT delete; $aCT delete; $bCT delete; 
$nearL2 delete 
#############
# Get lipids near TCR anchors (last 6JXR residue)
set eAn [atomselect $tcr "chain E and mass>2 and resid 155"]
set fAn [atomselect $tcr "chain F and mass>2 and resid 155"]
set dAn [atomselect $tcr "chain D and mass>2 and resid 129"]
set gAn [atomselect $tcr "chain G and mass>2 and resid 138"]
set aAn [atomselect $tcr "chain A and mass>2 and resid 57"]
set bAn [atomselect $tcr "chain B and mass>2 and resid 57"]
$eAn set beta 3
$fAn set beta 3
$dAn set beta 3
$gAn set beta 3
$aAn set beta 3
$bAn set beta 3
set nearL3 [atomselect $tcr "beta 1 and within 40 of beta 3"]
set resNear3 [lsort -unique [$nearL3 get residue]]
llength $resNear3; # 77 residues within 20 of anchors
set tot 0
foreach r {CHL1 POPC DPPC DOPC POPS DOPE POPI} {
  set res [atomselect $tcr "residue $resNear3 and resname $r"] 
  set n($r) [llength [lsort -unique [$res get residue]]]
  $res delete
  set tot [expr $tot + $n($r) ]
}
puts -nonewline $out3 [format "%6s" $fName]
foreach r {CHL1 POPC DPPC DOPC POPS DOPE POPI} {
puts -nonewline $out3 [format "%4d" $n($r)]
}
puts $out3 [format "%4d" $tot]
$eAn delete; $fAn delete; $dAn delete; $gAn delete; $aAn delete;$bAn delete
# Done with this model
mol delete $tcr
}
}
close $out1; close $out2; close $out3
#
###########################################################
# Four boxplot graphs of lipid distribution 
# Figure 18A-D
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr38lipid.jpg'
# layout: rows, columns
set multiplot layout 2,2
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 8
set xtic out nomirror rotate by -45 offset -1 
set ytic out nomirror offset 1,0
set xrange [0:21]
set grid ytics
set ylabel "{/:Bold Lipid number}" offset 2,0
set style boxplot nooutliers pointtype 7
set style data boxplot
set boxwidth 0.5
set key horizontal textcolor variable samplen -1
set xtics ("CHL" 1.5, "POPC" 4.5, "DPPC" 7.5, "DOPC" 10.5, "POPS" 13.5, "DOPE" 16.5, "POPI" 19.5)
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 7 pi -1 ps 1.0
#
set title "{/:Bold Within 20 Å of all protein}"
f38='tcr38k85a.dat'; f39='tcr39k85a.dat'
#
plot f39 u (1):2 with boxplot ls 1 t "",\
f38 u (2):2 with boxplot ls 2 t "",\
f39 u (4):3 with boxplot ls 1 t "",\
f38 u (5):3 with boxplot ls 2 t "",\
f39 u (7):4 with boxplot ls 1 t "",\
f38 u (8):4 with boxplot ls 2 t "",\
f39 u (10):5 with boxplot ls 1 t "",\
f38 u (11):5 with boxplot ls 2 t "",\
f39 u (13):6 with boxplot ls 1 t "",\
f38 u (14):6 with boxplot ls 2 t "",\
f39 u (16):7 with boxplot ls 1 t "",\
f38 u (17):7 with boxplot ls 2 t "",\
f39 u (19):8 with boxplot ls 1 t "TCR",\
f38 u (20):8 with boxplot ls 2 t "TCR-GOF"
#
set title "{/:Bold Within 20 Å of CD3 CTs}"
f38='tcr38k85b.dat'; f39='tcr39k85b.dat'
replot
#
set title "{/:Bold Within 20 Å of CD3 anchors}"
f38='tcr38k85c.dat'; f39='tcr39k85c.dat'
replot
#
set title "{/:Bold Within 40 Å of CD3 anchors}"
f38='tcr38k85d.dat'; f39='tcr39k85d.dat'
replot
eof
###########################################################
# Use R to get t-tests for lipids
read.table("tcr38k85d.dat")->t38
read.table("tcr39k85d.dat")->t39
p <- numeric()
m38 <- numeric()
m39 <- numeric()
dis <- numeric()
for(i in seq(2,8,1)) {
  # two-sided t-test of columns 2-8
  t.test(t39[,i],t38[,i])->q  
  p[i-1] <- format(q$p.value, digits =3)
  m39[i-1] <- round(mean(t39[, i]),digits=3) 
  m38[i-1] <- round(mean(t38[, i]),digits=3)
  dis[i-1] <- round((mean(t38[, i]) - mean(t39[, i])),digits=3) 
}
n<-c("CHL","POPC","DPPC","DOPC","POPS","DOPE","POPI")
d<-data.frame(n,m39,m38,p,dis)
write.table(d,"tcr38kd.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
###########################################################
# All protein
cat tcr38ka.dat
CHL 120.146 114.5 0.233 -5.646
POPC 52.049 45.387 0.000835 -6.662
DPPC 30.902 32.306 0.277 1.404
DOPC 33.024 29.71 0.00773 -3.315
POPS 36.098 26.984 1.13e-11 -9.114
DOPE 17.488 16.677 0.209 -0.81
POPI 16.415 13.032 1.22e-05 -3.382
# CD3 CTs
cat tcr38kb.dat
CHL 112.073 107.903 0.381 -4.17
POPC 49.78 43.5 0.00204 -6.28
DPPC 29.829 31.29 0.273 1.461
DOPC 31.561 28.371 0.0142 -3.19
POPS 35.098 26.129 1.9e-11 -8.969
DOPE 16.902 16.226 0.299 -0.677
POPI 15.976 12.742 2.55e-05 -3.234
# 20 from anchors
cat tcr38kc.dat
CHL 26.537 33.952 2.41e-16 7.415
POPC 14.659 11.919 4.25e-09 -2.739
DPPC 6.78 8.839 1.26e-07 2.058
DOPC 9.439 9.968 0.169 0.529
POPS 13.22 9.935 1.97e-12 -3.284
DOPE 4.268 6.323 4.42e-12 2.054
POPI 6.195 4.048 1.5e-13 -2.147
# 40 from anchors
cat tcr38kd.dat
CHL 93.415 100.887 8.83e-11 7.472
POPC 40.024 38.79 0.0903 -1.234
DPPC 21.195 27.177 3.62e-17 5.982
DOPC 24.659 25.306 0.203 0.648
POPS 27.415 21.774 9.41e-21 -5.64
DOPE 13.488 14.387 0.0476 0.899
POPI 12.341 10.048 6.34e-09 -2.293
###########################################################
#
