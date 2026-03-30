# anEdt12.sh: PRS dihedral cluster of 200-member edt and fgt ensembles (S5 Fig).
cd /t/tcr16/eds
# Use protein-only edt.prot.pdb and edtEns.dcd
# Use protein-only fgt.prot.pdb and fgtEns.dcd
###########################################################
# Edt ensemble clustering without intermediate dihedral files.
cpptraj<<eof
parm /t/tcr16/edt.prot.pdb
trajin /t/tcr16/edtEns.dcd
# Pro182ψ and Pro184ψ show angles that disfavor Nck binding
dihedral P182psi :132@N :132@CA :132@C :133@N
dihedral P184psi :134@N :134@CA :134@C :135@N
# This "go" is necessary before the write command.
go
write edtP182psi.dat P182psi
write edtP184psi.dat P184psi
# Make COORDS set crd1
createcrd crd1
# It is unclear why a COORDS set is needed for the clustering because it is the dihedral data
# that is clustered, which is irrelevant to coordinates at this point.
cluster crdset crd1 c1 kmeans data P182psi,P184psi clusters 5 out edt5.out \
summary edt5.summary.dat info edt5.info.dat 
eof
###########################################################
# Fgt ensemble clustering without intermediate dihedral files.
cpptraj<<eof
parm /t/tcr16/fgt.prot.pdb
trajin /t/tcr16/fgtEns.dcd
# Pro182ψ and Pro184ψ show angles that disfavor Nck binding
dihedral P182psi :60@N :60@CA :60@C :61@N
dihedral P184psi :62@N :62@CA :62@C :63@N
# This "go" is necessary before the write command.
go
write fgtP182psi.dat P182psi
write fgtP184psi.dat P184psi
# Make COORDS set crd1
createcrd crd1
cluster crdset crd1 c1 kmeans data P182psi,P184psi clusters 5 out fgt5.out \
summary fgt5.summary.dat info fgt5.info.dat 
eof
###########################################################
# Combine angle and cluster assignments
paste edtP182psi.dat edtP184psi.dat edt5.out > edtAng2.dat
paste fgtP182psi.dat fgtP184psi.dat fgn5.out > fgtAng2.dat
###########################################################
# Plot angles color-coded by cluster
# S5 Fig panel A 2x3 panels
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 2250, 1120
# ^ size is width, height
set out 'dimer.ang2.jpg'
# layout: rows, columns
set multiplot layout 2,3
set border 3
set xtic in nomirror 
set xrange [-100:200]
set yrange [-100:200]
set ytic in nomirror
set xzeroaxis
set yzeroaxis
set xlabel "{/:Bold P182ψ}"
set ylabel "{/:Bold P184ψ}"
# For the "lc var" to work in the plot command, use "linetypes", not "style line".
# Colors: https://i.stack.imgur.com/x6yLm.png
set linetype 1 lc rgb "blue" lw 3 pt 7 pi -1 ps 1.0
set linetype 2 lc rgb "red" lw 3 pt 7 pi -1 ps 1.0
set linetype 3 lc rgb "forest-green" lw 3 pt 7 pi -1 ps 1.0
set linetype 4 lc rgb "purple" lw 3 pt 7 pi -1 ps 1.0
set linetype 5 lc rgb "dark-cyan" lw 3 pt 7 pi -1 ps 1.0
fname="edtAng2.dat"
set title "{/:Bold CD3εδ angles}"
# P182psi in column 2, P184psi in column 4, cluster# (0-4) in column 6
plot fname using 2:4:($6+1) with points ls 1 lc var t ""
#
set title "{/:Bold CD3εδ^N angles}"
fname="eAng2.dat"
replot
#
set title "{/:Bold CD3εδ^R angles}"
fname="edrAng2.dat"
replot
#
set title "{/:Bold CD3εγ angles}"
fname="fgtAng2.dat"
replot
#
set title "{/:Bold CD3εγ^N angles}"
fname="fAng2.dat"
replot
#
set title "{/:Bold CD3εγ^R angles}"
fname="fgrAng2.dat"
replot
quit
eof

###########################################################
gnuplot<<"eof"
set term jpeg font "arial.ttf,18"
set out 'edn.ang2.jpg'
#set out 'fgn.ang2.jpg'
set title "{/:Bold CD3εδ^N angles}"
#set title "{/:Bold CD3εγ^N angles}"
set border 3
set xtic in nomirror 
set xrange [-100:200]
set yrange [-100:200]
set ytic in nomirror
set xzeroaxis
set yzeroaxis
set xlabel "{/:Bold P182ψ}"
set ylabel "{/:Bold P184ψ}"
# For the "lc var" to work in the plot command, use "linetypes", not "style line".
# Colors: https://i.stack.imgur.com/x6yLm.png
set linetype 1 lc rgb "blue" lw 3 pt 7 pi -1 ps 1.0
set linetype 2 lc rgb "red" lw 3 pt 7 pi -1 ps 1.0
set linetype 3 lc rgb "forest-green" lw 3 pt 7 pi -1 ps 1.0
set linetype 4 lc rgb "purple" lw 3 pt 7 pi -1 ps 1.0
set linetype 5 lc rgb "dark-cyan" lw 3 pt 7 pi -1 ps 1.0
# P182psi in column 2, P184psi in column 4, cluster# (0-4) in column 6
plot 'eAng2.dat' using 2:4:($6+1) with points ls 1 lc var t ""
#plot 'fAng2.dat' using 2:4:($6+1) with points ls 1 lc var t ""
quit
eof
