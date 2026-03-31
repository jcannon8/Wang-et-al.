# anTcr38v4.tcl: Report anchor radii for several tcr38b & tcr39b CT residues (Fig 19).
# Single 425 ns (eq85) frame analyzed.
# The goal is to get max radius and mean Z dimension for each anchor.
# Derived from anTcr38v3.tcl
# vmd -dispdev none
cd /t/tcr39an
###########################################################
set m tcr39; # =tcr38 for tcr38b; =tcr39 for tcr39b
set out1 [open ${m}v85h.dat w]
for {set j 0} {$j<=11} {incr j} {
foreach k {00 01 02 10 11 12 20 21 22} {
set fName $m.$j.$k 
# Skip missing  models. These are 425 ns frames.
if {[file exist /t/${m}b/$fName.psf]==0} continue
if {[file exist /t/${m}b/$fName.eq85.rst ]==0} continue
set tcr [mol new /t/${m}b/$fName.psf]
mol addfile /t/${m}b/$fName.eq85.rst type netcdf waitfor all
puts -nonewline $out1 [format "%6s" $fName]
foreach v {0 2 4 8 12 20} {
# Anchor points for six CD3 subunits, v'th CT residue CA.
set w [expr 155+$v]
set eV [atomselect $tcr "name CA and chain E and resid $w"]
set fV [atomselect $tcr "name CA and chain F and resid $w"]
set w [expr 129+$v]
set dV [atomselect $tcr "name CA and chain D and resid $w"]
set w [expr 138+$v]
set gV [atomselect $tcr "name CA and chain G and resid $w"]
set w [expr 57+$v]
set aV [atomselect $tcr "name CA and chain A and resid $w"]
set bV [atomselect $tcr "name CA and chain B and resid $w"]
# https://www.ks.uiuc.edu/Research/vmd/vmd-1.7.1/ug/node163.html
# Get anchor 2D coordinates for max XY radius.
set eVc [list [$eV get x] [$eV get y]]
set fVc [list [$fV get x] [$fV get y]]
set dVc [list [$dV get x] [$dV get y]]
set gVc [list [$gV get x] [$gV get y]]
set aVc [list [$aV get x] [$aV get y]] 
set bVc [list [$bV get x] [$bV get y]]
# Get anchor Z coordinates for average Z dimension.
set eVz [$eV get z]
set fVz [$fV get z]
set dVz [$dV get z]
set gVz [$gV get z]
set aVz [$aV get z] 
set bVz [$bV get z]
# Average Z dimension
set aveZ [expr (1.0/6) * ($eVz+$fVz+$dVz+$gVz+$aVz+$bVz)]
# Calculate anchor XY center
set cen [vecscale [expr 1.0/6] [vecadd $eVc $fVc $dVc $gVc $aVc $bVc]]
# anchor distances from center
set eDis [vecdist $eVc $cen]
set fDis [vecdist $fVc $cen]
set dDis [vecdist $dVc $cen]
set gDis [vecdist $gVc $cen]
set aDis [vecdist $aVc $cen]
set bDis [vecdist $bVc $cen]
# Get maximum distance for anchor circle radius, max.
set max 0; set min 100
# Get anchor nearest to center, near; farthest from center, far.
foreach d "$eDis $fDis $dDis $gDis $aDis $bDis" chn {E F D G A B} {
if {$d > $max} {set max $d; set far $chn}
if {$d < $min} {set min $d; set near $chn}
}
# ^ max is now the radius of the circle of v'th CT anchor
puts -nonewline $out1 [format "%7.2f%7.2f" $max $aveZ]
$eV delete; $fV delete; $dV delete; $gV delete; $aV delete; $bV delete
}
puts $out1 " "; # Terminate output line
mol delete $tcr
}
}
# Output: <model><CT2 max radius><CT2 aveZ><CT4 max radius><CT4 aveZ>...
close $out1
#
###########################################################
# Use R to get mean and SD for anchor circle radii and Z dimension
read.table("tcr38v85h.dat")->t38
read.table("tcr39v85h.dat")->t39
v<- c("0","2","4","8","12","20") 
r38 <- numeric()
r39 <- numeric()
z38 <- numeric()
z39 <- numeric()
rs38 <- numeric()
rs39 <- numeric()
zs38 <- numeric()
zs39 <- numeric()
# Number of anchors
numCT = (length(t38[1,])-1)/2
# mean Z dimension of TCR CT2
mZ38=mean(t38[,3])
mZ39=mean(t39[,3])
for(i in seq(1,numCT,1)) {
  r39[i] <- round(mean(t39[,i*2]),digits=2) 
  r38[i] <- round(mean(t38[,i*2]),digits=2)
  rs39[i] <- round(sd(t39[,i*2]*0.5),digits=2)
  rs38[i] <- round(sd(t38[,i*2]*0.5),digits=2)
# Z dimension is relative to mean CT0 of the ensemble
  z39[i] <- round(mean(t39[,i*2+1])- mZ39,digits=2) 
  z38[i] <- round(mean(t38[,i*2+1])- mZ38,digits=2)
  zs39[i] <- round(sd(t39[,i*2+1]*0.5),digits=2) 
  zs38[i] <- round(sd(t38[,i*2+1]*0.5),digits=2)
}
# gnuplot Xyerrorbars needs
# 4 columns: x y xdelta ydelta
d38<-data.frame(v,z38,r38,zs38,rs38)
write.table(d38,"tcr38xy.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
d39<-data.frame(v,z39,r39,zs39,rs39)
write.table(d39,"tcr39xy.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
###########################################################
# Plot anchor point positions with linespoints (best) or xyerrorbars
# Z dimension versus Radii
# Fig 19A
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 500, 375
set out 'tcr38vxy2.jpg'
set title "{/:Bold CD3 CT Cα positions}"
set border 3
set xtic nomirror 
set ytic out nomirror
set grid y
set ylabel "{/:Bold Z dimension (Å)}" offset 1,0
set xlabel "{/:Bold Radii (Å)}" offset 1,0
set key bottom right Left reverse  textcolor variable samplen -1
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 7 pi -1 ps 1.0
# Z dimension versus radii
plot 'tcr39xy.dat' u 3:2 with linespoints ls 1 t "TCR",\
'tcr38xy.dat' u 3:2  with linespoints ls 2 t "TCR-GOF",\
'tcr39xy.dat' u 3:2:1 with labels point tc ls 1 offset char 0.5,-0.7 t ""
eof
#
