# anTcr39nck2.tcl: Plot Nck binding to tcr39*.** (Fig 15).
# Similar to anTcr39e.tcl
# vmd -dispdev none
# Collisons are heavy atoms less than or equal to 1.0 angstroms.
cd /t/tcr39an
# Using data from anTcr39r.tcl
# anTcr39r.tcl output are 62 tcr39.*.**.Nckt3.dat files.
# Example *.Nckt3.dat output:
#     76   0  24    3.10   1.08  -8.48   0.57 |  125.67  129.85 |  131.20  132.49
# <eq frame><eNck binding><fNck binding><eMem1><eMem2><fMem1><fMem2> |
# Columns 1 & 2 have Nck-PRS reference # or -1 if no Nck binding without protein collision.
# Columns 3 & 5 have Nck distance to average membrane.
# Columns 3 & 6 have shortest Nck distance to membrane.
set out [open tcr39nck2.dat w]; # Output for Nck binding criteria
set outE [open tcr39nck3e.dat w]; # Output for eNck to membrane distances
set outF [open tcr39nck3f.dat w]; # Output for fNck to membrane distances
set files [glob *.Nckt3.dat]; # Output files from anTcr39r.tcl for each tcr39b model.
set nMod 0
foreach f $files {
  puts $f
  incr nMod
  set in [open $f r]
  set inData [read $in]
  close $in
  foreach line [split $inData "\n"] {
    if {$line == ""} break
    set tm [lindex $line 0];    # eq frame (5 ns).
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
    # Testing binding without collisions
    # No protein collision exlcuding fgr.201 if >-1
    if {$eProt>-1 && $eProt<24} {incr pConE($tm)}
    if {$fProt>-1 && $fProt<24} {incr pConF($tm)}
    # no membrane collision
    if {$eMem2 >1.0} {incr mConE($tm)}
    if {$fMem2 >1.0} {incr mConF($tm)}
    # 7-27-24 criteria: P182 ψ >100°, P184 ψ >100° permissive for binding.
    if {$eP182>100 && $eP184>100} {incr aConE($tm)}
    if {$fP182>100 && $fP184>100} {incr aConF($tm)} 
    # Nck can bind if no collision with other proteins, 
    # no Nck collision with membrane, and no bad angles.
    if {$eProt>-1 && $eProt<24 && $eMem2>1.0 && $eMem2 !=99 && $eP182>100 && $eP184>100} {
      incr nBindE($tm); # Count how many Nck bindings at time, tm
      # Save eMem2 for average and SD.
      lappend eMemSum($tm) $eMem2
    }
    if {$fProt>-1 && $fProt<24 && $fMem2>1.0 && $fMem2 !=99 && $fP182>100 && $fP184>100} {
      incr nBindF($tm)
      # Save fMem2 for average and SD.
      lappend fMemSum($tm) $fMem2
    }
}
}
for {set tm 4} {$tm<=85} {incr tm} {
# Output, tcr38nck2.dat, from this script in percents :
# <eP182psi><eP184psi> | <fP182psi><fP184psi>
# <time><Pe><Pf><Me><Mf><Ae><Af><Te><Tf>
    puts $out [format "%4d %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f" \
        [expr $tm*5] \
        [expr $pConE($tm)*100./$nMod] [expr $pConF($tm)*100./$nMod] \
        [expr $mConE($tm)*100./$nMod] [expr $mConF($tm)*100./$nMod] \
        [expr $aConE($tm)*100./$nMod] [expr $aConF($tm)*100./$nMod] \
        [expr $nBindE($tm)*100./$nMod] [expr $nBindF($tm)*100./$nMod]]
# Output: list of Nck to membrane distances for each time.
  puts $outE $eMemSum($tm); puts $outF $fMemSum($tm)
}
close $out; close $outE; close $outF
############################################################
# To check one MD frame in all models
# awk '{if ($1==45) {print $0}}' *.Nckt3.dat
###########################################################
# Plot tcr39b and tcr38b Nck binding averages 
# With all criteria for each dimer, four panels.
# Figure 15A,B
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr3938b.sum9.jpg'
# layout: rows, columns
set multiplot layout 2,2 
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 3.5
set xtics out nomirror 0,100,430 add ("" 50,"" 150,"" 250,"" 350) offset 0,0.2
set ytics out nomirror offset 1,0
set xrange [0:430]
set yrange [0:100]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
#set ylabel "{/:Bold Potential Nck binding (%)}" offset 3,0
#set key Left reverse textcolor variable samplen -1 font "arial.ttf,12"
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
#
set title "{/:Bold TCR CD3εδ}" offset 0,-1
plot '/t/tcr39an/tcr39nck2.dat' u 1:2 with linespoints ls 1 t "",\
'' using 1:4 with linespoints ls 2 t "" ,\
'' using 1:6 with linespoints ls 3 t "" ,\
'' using 1:8 with linespoints ls 4 t ""  
#
set title "{/:Bold TCR CD3εγ}" offset 0,-1
plot '/t/tcr39an/tcr39nck2.dat' u 1:3 with linespoints ls 1 t "" ,\
'' using 1:5 with linespoints ls 2 t "" ,\
'' using 1:7 with linespoints ls 3 t "" ,\
'' using 1:9 with linespoints ls 4 t ""  
#
set title "{/:Bold TCR-GOF CD3εδ}" offset 0,-1
set key bottom left Left reverse textcolor variable samplen -1 font "arial.ttf,12"
# Legend split between two plots to accomodate 12 pt font height.
plot '/t/tcr38an/tcr38nck2.dat' u 1:2 with linespoints ls 1 t "No protein collision" ,\
'' using 1:4 with linespoints ls 2 t "No membrane collision" ,\
'' using 1:6 with linespoints ls 3 t "" ,\
'' using 1:8 with linespoints ls 4 t ""
#
set title "{/:Bold TCR-GOF CD3εγ}" offset 0,-1
plot '/t/tcr38an/tcr38nck2.dat' u 1:3 with linespoints ls 1 t "",\
'' using 1:5 with linespoints ls 2 t "",\
'' using 1:7 with linespoints ls 3 t "Favorable PRS angles",\
'' using 1:9 with linespoints ls 4 t "Potential Nck binding"  
quit
eof
###########################################################
# Plot tcr39b and tcr38b Nck binding averages
# With all criteria for each dimer, four panels.
# Figure 15C
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,375
set out 'tcr3938c.sum9.jpg'
# layout: rows, columns
set multiplot layout 1,2 
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 3.5
set xtics out nomirror 0,100,430 add ("" 50,"" 150,"" 250,"" 350) offset 0,0.2
set ytics out nomirror offset 1,0
set xrange [0:430]
set yrange [0:60]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
#set ylabel "{/:Bold Potential Nck binding (%)}" offset 2,0
set key Left top reverse textcolor variable samplen -1 
set style line 1 lc rgb "red" lt 1 lw 2 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "blue" lt 1 lw 2 pt 5 pi -1 ps 1.0
set style line 3 lc rgb "dark-green" lt 1 lw 2 pt 9 pi -1 ps 1.0
set style line 4 lc rgb "magenta" lt 1 lw 2 pt 11 pi -1 ps 1.0
#
set title "{/:Bold CD3εδ}" offset 0,-1
plot '/t/tcr39an/tcr39nck2.dat' u 1:8 with linespoints ls 1 t "TCR",\
'/t/tcr38an/tcr38nck2.dat' u 1:8 with linespoints ls 2 t "TCR-GOF"
#
set title "{/:Bold CD3εγ}" offset 0,-1
plot '/t/tcr39an/tcr39nck2.dat' u 1:9 with linespoints ls 1 t "TCR",\
'/t/tcr38an/tcr38nck2.dat' u 1:9 with linespoints ls 2 t "TCR-GOF"
quit
eof
###########################################################
# Use R to get Nck to membrane mean and SD.
# tcr38nck3e
# Need to accomdate possible 41 columns and fill blanks with NA.
read.table("/t/tcr38an/tcr38nck3e.dat",header=FALSE,col.names=seq(1,82),fill=TRUE)->e38
# Treat data as numeric and strip NA before computation.
nRows=length(e38[,1])
out<-matrix(,nrow=nRows,ncol=3) 
for(i in seq(1,nRows,1)) {
  (i+3)*5-> out[i,1]
  round(mean(as.numeric(e38[i,]),na.rm=TRUE),3)-> out[i,2]
  round((sd(as.numeric(e38[i,]),na.rm=TRUE)/2),3)-> out[i,3]
}
write.table(out,file="/t/tcr38an/tcr38.nck4e.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
#
# tcr38nck3f
# Need to accomdate possible 41 columns and fill blanks with NA.
read.table("/t/tcr38an/tcr38nck3f.dat",header=FALSE,col.names=seq(1,82),fill=TRUE)->f38
# Treat data as numeric and strip NA before computation.
nRows=length(f38[,1])
out<-matrix(,nrow=nRows,ncol=3) 
for(i in seq(1,nRows,1)) {
  (i+3)*5-> out[i,1]
  round(mean(as.numeric(f38[i,]),na.rm=TRUE),3)-> out[i,2]
  round((sd(as.numeric(f38[i,]),na.rm=TRUE)/2),3)-> out[i,3]
}
write.table(out,file="/t/tcr38an/tcr38.nck4f.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
####
# tcr39nck3e
# Need to accomdate possible 41 columns and fill blanks with NA.
read.table("/t/tcr39an/tcr39nck3e.dat",header=FALSE,col.names=seq(1,41),fill=TRUE)->e39
# Treat data as numeric and strip NA before computation.
nRows=length(e39[,1])
out<-matrix(,nrow=nRows,ncol=3) 
for(i in seq(1,nRows,1)) {
  (i+3)*5-> out[i,1]
  round(mean(as.numeric(e39[i,]),na.rm=TRUE),3)-> out[i,2]
  round((sd(as.numeric(e39[i,]),na.rm=TRUE)/2),3)-> out[i,3]
}
write.table(out,file="/t/tcr39an/tcr39.nck4e.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
# tcr39nck3f
# Need to accomdate possible 41 columns and fill blanks with NA.
read.table("/t/tcr39an/tcr39nck3f.dat",header=FALSE,col.names=seq(1,41),fill=TRUE)->f39
# Treat data as numeric and strip NA before computation.
nRows=length(f39[,1])
out2<-matrix(,nrow=nRows,ncol=3) 
for(i in seq(1,nRows,1)) {
  (i+3)*5-> out2[i,1]
  round(mean(as.numeric(f39[i,]),na.rm=TRUE),3)-> out2[i,2]
  round((sd(as.numeric(f39[i,]),na.rm=TRUE)/2),3)-> out2[i,3]
}
write.table(out2,file="/t/tcr39an/tcr39.nck4f.dat",quote=FALSE,row.names=FALSE,col.names=FALSE)
# Comparison of Nck to membrane differences
e<- numeric()
d<- numeric()
for(i in seq(1,82,1)) {
t.test(e38[i,],e39[i,])$p.value-> e[i]
m<- mean(as.numeric(e38[i,]),na.rm=TRUE ) - mean(as.numeric(e39[i,]),na.rm=TRUE)
round(m,3)-> d[i]
}
# ^ 40 frames (40/82 = 49%) have a p<0.05
summary(d)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
 -1.857   2.901   5.158   4.666   6.277  10.650
# ^ Median e38-e39 difference = 5.158
f<- numeric()
for(i in seq(1,82,1)) {
t.test(f38[i,],f39[i,])$p.value-> f[i]
}
# ^ One time (1%) has a p<0.05 
# Get means at one time for E chain
t.test(e38[i,],e39[i,])
###########################################################
# Plot tcr39b and tcr38b Nck to membrane averages
# Figure 15D,E
gnuplot<<eof
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr.nck4.sum9.jpg'
# layout: rows, columns
set multiplot layout 2,2 
set border 3
# margin units are character heights or widths
set rmargin 1
set lmargin 3.5
set xtics out nomirror 0,100,430 add ("" 50,"" 150,"" 250,"" 350) offset 0,0.2
set ytics out nomirror 0,5,25 offset 1,0
set xrange [0:430]
set yrange [0:25]
set grid y
set xlabel "{/:Bold Time (ns)}" offset 0,0.8
#set ylabel "{/:Bold Nck to membrane (\305)}" offset 2,0
set style line 1 lc rgb "blue" lt 1 lw 3 pt 7 pi -1 ps 1.0
set style line 2 lc rgb "red" lt 1 lw 3 pt 7 pi -1 ps 1.0
# tcr39.nck4e
set title "{/:Bold TCR CD3εδ}" offset 0,-1
plot '/t/tcr39an/tcr39.nck4e.dat' u 1:2:3 with yerrorbars ls 1 t "",\
'' u 1:2 with lines ls 2 t ""
# tcr39.nck4f
set title "{/:Bold TCR CD3εγ}" offset 0,-1
plot '/t/tcr39an/tcr39.nck4f.dat' u 1:2:3 with yerrorbars ls 1 t "",\
'' u 1:2 with lines ls 2 t ""
# tcr38.nck4e
set title "{/:Bold TCR-GOF CD3εδ}" offset 0,-1
plot '/t/tcr38an/tcr38.nck4e.dat' u 1:2:3 with yerrorbars ls 1 t "",\
'' u 1:2 with lines ls 2 t ""
# tcr38.nck4f
set title "{/:Bold TCR-GOF CD3εγ}" offset 0,-1
plot '/t/tcr38an/tcr38.nck4f.dat' u 1:2:3 with yerrorbars ls 1 t "",\
'' u 1:2 with lines ls 2 t ""
quit
eof


