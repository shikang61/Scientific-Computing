set terminal pngcairo size 1400, 700 enhanced font 'Verdana,12'
set output './euler_final2D/'.ARG1.'.png'
unset warnings

set xrange [*:*]    # Automatically scale to x data range
set yrange [*:*]    # Automatically scale to y data range
set cbrange [*:*]   # Automatically scale color bar to u1 data range

set xlabel "X-axis (x)"
set ylabel "Y-axis (y)"
set palette rgb 33,13,10   # Optional: set color palette (modify to suit preference)

datafile = ARG1.'.dat'

stats datafile using 1 nooutput
n_time_steps = STATS_blocks
i = n_time_steps-2
#i = 1

stats datafile index i using 1 nooutput
actual_time = STATS_min

set multiplot layout 1, 2 title sprintf("Final solution (2D) at t = (Time = %.3f)", actual_time)

# First plot: rho
set zlabel "Density"
set cblabel "rho"
unset key           
set view map        
set style data pm3d
set title "Density"
splot datafile index i using 2:3:4

# First plot: p
set zlabel "Pressure"
set cblabel "p"
unset key
set view map
set style data pm3d
set title "Pressure"
splot datafile index i using 2:3:7 


