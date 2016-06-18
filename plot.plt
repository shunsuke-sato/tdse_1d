set term postscript eps enhanced color
unset key
set nosurface
set contour
set view 0,0


set cntrparam levels discrete 0.1, 0.2, 0.3, 0.4, 0.5, 0.6

set output "wfn.eps"
#splot [-5:5][-5:5] "GS_wfn.out" matrix w l

set xtics ("-5" 151, "0" 200, "5" 251) 
set ytics ("-5" 151, "0" 200, "5" 251) 
unset ztics
set xrange [151:251]
set yrange [151:251]
set size square

splot "GS_wfn_matrix.out" matrix w l linecolor 7 lt 1
unset output