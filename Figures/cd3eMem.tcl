# cd3eMem.tcl: Plot membrane contacts in CD3e subunit residues (S13 fig).
# Use data made by /home/wangly/TCR_full_analysis/Scripts
# anTcr38t_LW.tcl, MakeupZero.sh,TCRtoMEMplot.tcl 
cd /t/tcr38an 
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 1600,600
set out 'cd3Mem.jpg'
# layout: rows, columns
set multiplot layout 1,2
set xlabel "{/:Bold Residue}"
set ylabel "{/:Bold Percent membrane contact}" 
set xrange [157:207]
set yrange [0:85]
set border 3
set xtic nomirror 
set ytic out nomirror
set grid y
set style line 1 lt 1 lw 2 pt 7 ps 1.5 lc rgb "red"
set style line 2 lt 1 lw 2 pt 7 ps 1.5 lc rgb "blue"
set key top right textcolor variable samplen -1
# Two plots
plot '/home/wangly/TCR_full_analysis/eq85/fixed_tcr39e85.dat' \
using 1:($2*100/42) with linespoints ls 1 t "TCR CD3εδ",\
'/home/wangly/TCR_full_analysis/eq85/fixed_tcr38e85.dat' \
using 1:($2*100/61) with linespoints ls 2 t "TCR-GOF CD3εδ"
#
plot '/home/wangly/TCR_full_analysis/eq85/fixed_tcr39f85.dat' \
using 1:($2*100/42) with linespoints ls 1 t "TCR CD3εγ",\
'/home/wangly/TCR_full_analysis/eq85/fixed_tcr38f85.dat' \
using 1:($2*100/61) with linespoints ls 2 t "TCR-GOF CD3εγ"
eof
###########################################################
# Four plots
set title "{/:Bold TCR CD3εδ}" offset 0,-3
plot '/home/wangly/TCR_full_analysis/eq85/fixed_tcr39e85.dat' \
using 1:($2*100/42) with linespoints ls 1 t ""
set title "{/:Bold TCR CD3εγ}"
plot '/home/wangly/TCR_full_analysis/eq85/fixed_tcr39f85.dat' \
using 1:($2*100/42) with linespoints ls 1 t ""
set title "{/:Bold TCR-GOF CD3εδ}"
plot '/home/wangly/TCR_full_analysis/eq85/fixed_tcr38e85.dat' \
using 1:($2*100/61) with linespoints ls 2 t ""
set title "{/:Bold TCR-GOF CD3εγ}"
plot '/home/wangly/TCR_full_analysis/eq85/fixed_tcr38f85.dat' \
using 1:($2*100/61) with linespoints ls 2 t ""
eof
