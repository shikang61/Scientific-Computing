set terminal pngcairo size 1200, 1000 enhanced font 'Verdana,12'
set output './euler_final/'.ARG1.'.png'
unset warnings

RHO = 3
VX = 4
VY = 5
VZ = 6
P = 7
BX = 8
BY = 9
BZ = 10
BMAG = 11

set xlabel "Position (x)"
set grid
set pointsize 0.5

datafile = ARG1.'.dat'  # Replace with your actual data file

stats datafile using 2 nooutput
x_min = STATS_min
x_max = STATS_max
set xrange [x_min:x_max]

stats datafile using 1 nooutput
n_time_steps = STATS_blocks
i = n_time_steps-2

stats datafile index i using 1 nooutput
actual_time = STATS_min

set multiplot layout 4, 2 title "Final MHD solution" title sprintf("Final time (Time = %.6f)", actual_time)

# First plot: rho
stats datafile index i using RHO nooutput
expand_range = 0.05*(STATS_max - STATS_min) + 0.01
rho_min = STATS_min - expand_range
rho_max = STATS_max + expand_range
set yrange [rho_min:rho_max]
set ylabel "rho"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:RHO with points title "rho"
set yrange[*:*]

# Second plot: vx
stats datafile index i using VX nooutput
expand_range = 0.05*(STATS_max - STATS_min) + 0.01
vx_min = STATS_min - expand_range
vx_max = STATS_max + expand_range
set yrange [vx_min:vx_max]
set ylabel "vx"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:VX with points title "vx"
set yrange[*:*]

# Second plot: vy
stats datafile index i using VY nooutput
expand_range = 0.05*(STATS_max - STATS_min) + 0.01
vy_min = STATS_min - expand_range
vy_max = STATS_max + expand_range
set yrange [vy_min:vy_max]
set ylabel "vy"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:VY with points title "vy"
set yrange[*:*]

# Second plot: vz
stats datafile index i using VZ nooutput
expand_range = 0.05*(STATS_max - STATS_min) + 0.01
vz_min = STATS_min - expand_range
vz_max = STATS_max + expand_range
set yrange [vz_min:vz_max]
set ylabel "vz"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:VZ with points title "vz"
set yrange[*:*]

# Second plot: p
stats datafile index i using P nooutput
expand_range = 0.05*(STATS_max - STATS_min) + 0.01
p_min = STATS_min - expand_range
p_max = STATS_max + expand_range
set yrange [p_min:p_max]
set ylabel "P"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:P with points title "P"
set yrange[*:*]

# # Third plot: Bx
# stats datafile index i using BX nooutput 
# expand_range = 0.05*(STATS_max - STATS_min) + 0.01
# Bx_min = STATS_min - expand_range
# Bx_max = STATS_max + expand_range
# set yrange [Bx_min:Bx_max]
# set ylabel "B_x"
# set xzeroaxis linewidth 2.5
# plot datafile index i using 2:BX with points title "Bx"
# set yrange[*:*]

# Third plot: By
stats datafile index i using BY nooutput 
expand_range = 0.05*(STATS_max - STATS_min) + 0.01
By_min = STATS_min - expand_range
By_max = STATS_max + expand_range
set yrange [By_min:By_max]
set ylabel "By"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:BY with points title "By"
set yrange[*:*]

# Third plot: Bz
stats datafile index i using BZ nooutput 
expand_range = 0.05*(STATS_max - STATS_min) + 0.01
Bz_min = STATS_min - expand_range
Bz_max = STATS_max + expand_range
set yrange [Bz_min:Bz_max]
set ylabel "Bz"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:BZ with points title "Bz"
set yrange[*:*]

# Third plot: B_mag
stats datafile index i using BMAG nooutput 
expand_range = 0.05*(STATS_max - STATS_min) + 0.01
B_mag_min = STATS_min - expand_range
B_mag_max = STATS_max + expand_range
set yrange [B_mag_min:B_mag_max]
set ylabel "|B|"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:BMAG with points title "|B|"

unset multiplot
set output
