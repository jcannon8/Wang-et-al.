# viewLegend.tcl: Color legend for TCR model subunits (Fig 14).
# Use graphics and draw VMD commands to make a legend for TCR models.
# This size (width, height) allows four models across full page width.
display resize 560 300; 
color Display Background white
axes location Off
set tcr [mol new]; # plain blank molecule for drawing canvas.
# The size and thickness control font size
# The encoding does not allow Greek letters
# size 3 text is about 24 pt on Adobe Illustrator.
# Three columns at x=-1, 0, and 1.
graphics $tcr delete all
#
set x -1;  # X-dimension, column 
set y 0.3; # Y-dimension, row. Decrease by 0.3
set x2 [expr $x+0.2];set x3 [expr $x+0.4]
puts "x=$x x2=$x2 x3=$x3 y=$y"
# CD3zeta: chain a green
graphics $tcr color 7
draw text [list $x $y 0] "A" size 3 thickness 3
draw line [list $x2 $y 0] [list $x3 $y 0] width 4
set y 0
# CD3zeta: chain b purple
graphics $tcr color 11
draw text [list $x $y 0] "B" size 3 thickness 3
draw line [list $x2 $y 0] [list $x3 $y 0] width 4
set y -0.3
# CD3delta: chain d red
graphics $tcr color 1
draw text [list $x $y 0] "D" size 3 thickness 3
draw line [list $x2 $y 0] [list $x3 $y 0] width 4
set y -0.6
# CD3gamma: chain g orange
graphics $tcr color 3
draw text [list $x $y 0] "G" size 3 thickness 3
draw line [list $x2 $y 0] [list $x3 $y 0] width 4
################
# Second column
set x 0;  # X-dimension, column 
set y 0.3; # Y-dimension, row. Decrease by 0.3
set x2 [expr $x+0.2];set x3 [expr $x+0.4]
# CD3epsilon: chain e blue
graphics $tcr color 0
draw text [list $x $y 0] "E" size 3 thickness 3
draw line [list $x2 $y 0] [list $x3 $y 0] width 4
set y 0
# CD3epsilon: chain f grey
graphics $tcr color 2
draw text [list $x $y 0] "F" size 3 thickness 3
draw line [list $x2 $y 0] [list $x3 $y 0] width 4
set y -0.3
# TCRalpha: chain m pink
graphics $tcr color 9
draw text [list $x $y 0] "M" size 3 thickness 3
draw line [list $x2 $y 0] [list $x3 $y 0] width 4
set y -0.6
# TCRbeta: chain n yellow3
graphics $tcr color 18
draw text [list $x $y 0] "N" size 3 thickness 3
draw line [list $x2 $y 0] [list $x3 $y 0] width 4
################
# Third column
set x 1;  # X-dimension, column 
set y 0.3; # Y-dimension, row. Decrease by 0.3
set x2 [expr $x+0.2];set x3 [expr $x+0.4]
# Structural lipids black Licorice
graphics $tcr color 16
draw text [list $x $y 0] "Lipids" size 3 thickness 3
set y 0
# Nck magenta2
graphics $tcr color 28
draw text [list $x $y 0] "Nck" size 3 thickness 3
set y -0.3
# Lck ochre NewCartoon
graphics $tcr color 14
draw text [list $x $y 0] "Lck" size 3 thickness 3
