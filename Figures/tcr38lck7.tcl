# tcr38lck7.tcl: Plot Lck binding criteria for each ITAM (Fig 15C-G).
# Uses data from tcr38lck2.tcl
# derived from tcr38lck6.tcl
# This uses different fill patterns and consistent TCR red, TCR-GOF blue
# Figure 15C: Total Lck binding potential
cd /t/tcr38an
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr.lck7c.jpg'
# layout: rows, columns
set multiplot layout 2,1 
set border 3
set xtic nomirror rotate by -45 offset -1
set ytic out nomirror
set xtics ( \
    "AY72" 0, "AY83" 1, "AY111" 2, "AY123" 3, "AY142" 4, "AY153" 5, "BY72" 6, \
    "BY83" 7, "BY111" 8, "BY123" 9, "BY142" 10, "BY153" 11, "DY149" 12, "DY160" 13, \
    "EY188" 14, "EY199" 15, "FY188" 16, "FY199" 17, "GY160" 18, "GY171" 19)
set grid ytics
set ylabel "{/:Bold Frequency (%)}"
set xrange [-0.5:19.5]
set yrange [0:7]
set style line 1 lc rgb "blue" 
set style line 2 lc rgb "red" 
set boxwidth 0.9 relative
set style data histograms
set style histogram cluster
#set style fill solid 1.0 border lt -1
set title "{/:Bold TCR Lck ITAM binding}" tc rgb "red"
# tcr39lck4.dat data:
# <ITAM chain><ITAM res><RMSDa><RMSDb><Proxa><Proxb><CTa><CTb><Mema><Memb><Binda><Bindb>
# Divide tally by total frames = frames * models
set key left Left reverse tc rgb "red"
plot '/t/tcr39an/tcr39lck4.dat' using (100*$11/1722) ls 2 fs pattern 4 t "0-225 ns",\
'' u (100*$12/1640) ls 2 fs pattern 3 t "225-425 ns"
#
set key left Left tc rgb "blue"
set title "{/:Bold TCR-GOF Lck ITAM binding}" tc rgb "blue"
plot '/t/tcr38an/tcr38lck4.dat' using (100*$11/2604) ls 1 fs pattern 4 t "0-225 ns",\
'' u (100*$12/2480) ls 1 fs pattern 3 t "225-425 ns"
eof
###########################################################
# Lck binding criteria for each ITAM for
# Fig 16D: Favorable RMSD 
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr.lck7d.jpg'
# layout: rows, columns
set multiplot layout 2,1 
set border 3
set xtic nomirror rotate by -45 offset -1
set ytic out nomirror
set xtics ( \
    "AY72" 0, "AY83" 1, "AY111" 2, "AY123" 3, "AY142" 4, "AY153" 5, "BY72" 6, \
    "BY83" 7, "BY111" 8, "BY123" 9, "BY142" 10, "BY153" 11, "DY149" 12, "DY160" 13, \
    "EY188" 14, "EY199" 15, "FY188" 16, "FY199" 17, "GY160" 18, "GY171" 19)
set grid ytics
set ylabel "{/:Bold Frequency (%)}"
set xrange [-0.5:19.5]
set yrange [40:100]
set style line 1 lc rgb "blue" 
set style line 2 lc rgb "red" 
set boxwidth 0.9 relative
set style data histograms
set style histogram cluster
set title "{/:Bold TCR Lck ITAM RMSD}" tc rgb "red" 
# tcr39lck4.dat data:
# <ITAM chain><ITAM res><RMSDa><RMSDb><Proxa><Proxb><CTa><CTb><Mema><Memb><Binda><Bindb>
# Divide tally by total frames = frames * models
set key left Left reverse tc rgb "red" at first 0, first 100 
plot '/t/tcr39an/tcr39lck4.dat' using (100*$3/1722) ls 2 fs pattern 4 t "0-225 ns",\
'' u (100*$4/1640) ls 2 fs pattern 3 t "225-425 ns"
#
set key left Left tc rgb "blue" at first 0, first 100
set title "{/:Bold TCR-GOF Lck ITAM RMSD}" tc rgb "blue"
plot '/t/tcr38an/tcr38lck4.dat' using (100*$3/2604) ls 1 fs pattern 4 t "0-225 ns",\
'' u (100*$4/2480) ls 1 fs pattern 3 t "225-425 ns"
eof
###########################################################
# Fig 16E: No ITAM collision
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr.lck6e.jpg'
# layout: rows, columns
set multiplot layout 2,1 
set border 3
set xtic nomirror rotate by -45 offset -1
set ytic out nomirror
set xtics ( \
    "AY72" 0, "AY83" 1, "AY111" 2, "AY123" 3, "AY142" 4, "AY153" 5, "BY72" 6, \
    "BY83" 7, "BY111" 8, "BY123" 9, "BY142" 10, "BY153" 11, "DY149" 12, "DY160" 13, \
    "EY188" 14, "EY199" 15, "FY188" 16, "FY199" 17, "GY160" 18, "GY171" 19)
set grid ytics
set ylabel "{/:Bold Frequency (%)}"
set xrange [-0.5:19.5]
set yrange [0:40]
set style line 1 lc rgb "blue" 
set style line 2 lc rgb "red" 
set boxwidth 0.9 relative
set style data histograms
set style histogram cluster
set title "{/:Bold TCR Lck ITAM no proximal collision}" tc rgb "red"
# tcr39lck4.dat data:
# <ITAM chain><ITAM res><RMSDa><RMSDb><Proxa><Proxb><CTa><CTb><Mema><Memb><Binda><Bindb>
# Divide tally by total frames = frames * models
set key left Left reverse tc rgb "red" at first 0, first 40 
plot '/t/tcr39an/tcr39lck4.dat' using (100*$5/1722) ls 2 fs pattern 4 t "0-225 ns",\
'' u (100*$6/1640) ls 2 fs pattern 3 t "225-425 ns"
#
set key left Left tc rgb "blue" at first 0, first 40
set title "{/:Bold TCR-GOF Lck ITAM no proximal collision}" tc rgb "blue"
plot '/t/tcr38an/tcr38lck4.dat' using (100*$5/2604) ls 1 fs pattern 4 t "0-225 ns",\
'' u (100*$6/2480) ls 1 fs pattern 3 t "225-425 ns"
eof
###########################################################
# Fig 16F: ITAM no CT collision
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr.lck6f.jpg'
# layout: rows, columns
set multiplot layout 2,1 
set border 3
set xtic nomirror rotate by -45 offset -1
set ytic out nomirror
set xtics ( \
    "AY72" 0, "AY83" 1, "AY111" 2, "AY123" 3, "AY142" 4, "AY153" 5, "BY72" 6, \
    "BY83" 7, "BY111" 8, "BY123" 9, "BY142" 10, "BY153" 11, "DY149" 12, "DY160" 13, \
    "EY188" 14, "EY199" 15, "FY188" 16, "FY199" 17, "GY160" 18, "GY171" 19)
set grid ytics
set ylabel "{/:Bold Frequency (%)}"
set xrange [-0.5:19.5]
#set yrange [0:50]
set style line 1 lc rgb "blue" 
set style line 2 lc rgb "red" 
set boxwidth 0.9 relative
set style data histograms
set style histogram cluster
set title "{/:Bold TCR Lck ITAM no CT collision}" tc rgb "red"
# tcr39lck4.dat data:
# <ITAM chain><ITAM res><RMSDa><RMSDb><Proxa><Proxb><CTa><CTb><Mema><Memb><Binda><Bindb>
# Divide tally by total frames = frames * models
set key left Left reverse tc rgb "red" at first 12, first 70 
plot '/t/tcr39an/tcr39lck4.dat' using (100*$7/1722) ls 2 fs pattern 4 t "0-225 ns",\
'' u (100*$8/1640) ls 2 fs pattern 3 t "225-425 ns"
#
set key left Left tc rgb "blue" at first 12, first 70
set title "{/:Bold TCR-GOF Lck ITAM no CT collision}" tc rgb "blue"
plot '/t/tcr38an/tcr38lck4.dat' using (100*$7/2604) ls 1 fs pattern 4 t "0-225 ns",\
'' u (100*$8/2480) ls 1 fs pattern 3 t "225-425 ns"
eof
###########################################################
# F 16G: No membrane collision
gnuplot<<"eof"
set term jpeg font "arial.ttf,18" size 1000,750
set out 'tcr.lck6g.jpg'
# layout: rows, columns
set multiplot layout 2,1 
set border 3
set xtic nomirror rotate by -45 offset -1
set ytic out nomirror
set xtics ( \
    "AY72" 0, "AY83" 1, "AY111" 2, "AY123" 3, "AY142" 4, "AY153" 5, "BY72" 6, \
    "BY83" 7, "BY111" 8, "BY123" 9, "BY142" 10, "BY153" 11, "DY149" 12, "DY160" 13, \
    "EY188" 14, "EY199" 15, "FY188" 16, "FY199" 17, "GY160" 18, "GY171" 19)
set grid ytics
set ylabel "{/:Bold Frequency (%)}"
set xrange [-0.5:19.5]
set yrange [0:60]
set style line 1 lc rgb "blue" 
set style line 2 lc rgb "red" 
set boxwidth 0.9 relative
set style data histograms
set style histogram cluster
set title "{/:Bold TCR Lck ITAM no membrane collision}" tc rgb "red"
# tcr39lck4.dat data:
# <ITAM chain><ITAM res><RMSDa><RMSDb><Proxa><Proxb><CTa><CTb><Mema><Memb><Binda><Bindb>
# Divide tally by total frames = frames * models
set key left Left reverse tc rgb "red" at first 12, first 60 
plot '/t/tcr39an/tcr39lck4.dat' using (100*$9/1722) ls 2 fs pattern 4 t "0-225 ns",\
'' u (100*$10/1640) ls 2 fs pattern 3 t "225-425 ns"
#
set key left Left tc rgb "blue" at first 12, first 60
set title "{/:Bold TCR-GOF Lck ITAM no membrane collision}" tc rgb "blue"
plot '/t/tcr38an/tcr38lck4.dat' using (100*$9/2604) ls 1 fs pattern 4 t "0-225 ns",\
'' u (100*$10/2480) ls 1 fs pattern 3 t "225-425 ns"
eof
