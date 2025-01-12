set terminal pngcairo size 1400,1000 enhanced font 'Verdana,12'
set output './euler_final3D/'.ARG1.'.png'
unset warnings

set xrange [*:*]    # Automatically scale to x data range
set yrange [*:*]    # Automatically scale to y data range

set xlabel "X-axis (x)"
set ylabel "Y-axis (y)"

datafile = ARG1.'.dat'

stats datafile using 1 nooutput
n_time_steps = STATS_blocks
i = n_time_steps-2

stats datafile index i using 1 nooutput
actual_time = STATS_min

set multiplot layout 1, 2 title sprintf("Final solution (3D Mesh) at t = (Time = %.3f)", actual_time)

# General settings for mesh
unset pm3d               # Turn off pm3d to use lines for a mesh plot
set style data lines     # Use lines for a mesh appearance
set view 70, 40, 1.2, 1  # Set 3D view angles and scale for zoom (closer view)

# First plot: rho
set zlabel "Density"
unset key           
set title "Density"
splot datafile index i using 2:3:4 with lines

# Fourth plot: p
set zlabel "Pressure"
unset key
set title "Pressure"
splot datafile index i using 2:3:7 with lines

unset multiplot
