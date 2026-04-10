# tcr39.tyr4.sh: Plot tcr39.*.** Tyr to membrane distances (Fig 12 A-D).
# Derived from Gnuplot_LWv2_tcr38Box.sh Written by Liangyu Wang
cd /home/wangly/tyr_analysis/NoImage_Test
gnuplot<<eof
#set terminal pngcairo size 1280,960 enhanced font 'arial.ttf,18'
set term jpeg font "arial.ttf,18" size 2000,1500
set output 'tcr39.tyr85.jpg'
# layout: rows, columns
set multiplot layout 2,2
set border 3
set rmargin 1
set lmargin 3
set xtic nomirror rotate by -45 offset -1.2,0
set ytic out nomirror offset 1,0
set grid ytics
#set ylabel "{/:Bold Distance to membrane (Å)}"
set xrange [0.5:20.5]
set yrange [0:45]
set grid y
set style boxplot nooutliers pointtype 7
set style data boxplot
set boxwidth 0.5
set key off
set xtics ( \
    "AY72" 1, "AY83" 2, "AY111" 3, "AY123" 4, "AY142" 5, "AY153" 6, "BY72" 7, \
    "BY83" 8, "BY111" 9, "BY123" 10, "BY142" 11, "BY153" 12, "DY149" 13, "DY160" 14, \
    "EY188" 15, "EY199" 16, "FY188" 17, "FY199" 18, "GY160" 19, "GY171" 20)

set title "{/:Bold TCR at 225 ns}" offset 10,-2
file="/home/wangly/tyr_analysis/NoImage_Test/ITAMTYRtoMEM_NoImage45.dat"
# Color coded by chain
plot for [i=1:6] file u (i):(column(i)) with boxplot lc rgb "purple" lw 3 t "", \
     for [i=7:12] file u (i):(column(i)) with boxplot lc rgb "blue" lw 3 t "", \
     for [i=13:14] file u (i):(column(i)) with boxplot lc rgb "red" lw 3 t "", \
     for [i=15:16] file u (i):(column(i)) with boxplot lc rgb "orange" lw 3 t "", \
     for [i=17:18] file u (i):(column(i)) with boxplot lc rgb "dark-green" lw 3 t "", \
     for [i=19:20] file u (i):(column(i)) with boxplot lc rgb "light-green" lw 3 t ""
#     
set title "{/:Bold TCR at 425 ns}" offset 10,-2
file="/home/wangly/tyr_analysis/NoImage_Test/ITAMTYRtoMEM_NoImage85.dat"
replot
#
set title "{/:Bold TCR-GOF at 225 ns}" offset 10,-2
file="/home/wangly/tyr_analysis/TCR38_TYRAnlyz/tcr38_eq45/ITAMTYRtoMEM_TCR38eq45.dat"
replot
#
set title "{/:Bold TCR-GOF at 425 ns}" offset 10,-2
file="/home/wangly/tyr_analysis/TCR38_TYRAnlyz/tcr38_eq45/ITAMTYRtoMEM_TCR38eq85.dat"
replot
quit
eof
