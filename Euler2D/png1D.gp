set terminal pngcairo size 1200,1000 enhanced font 'Verdana,12'
set output './euler_final1D/'.ARG1.'.png'
unset warnings

set xlabel "Position"
set grid

set xrange [*:*]  # Replace xmin and xmax with appropriate values for your x

datafile = ARG1.'.dat'

stats datafile using 1 nooutput
n_time_steps = STATS_blocks
i = n_time_steps-2

stats datafile index i using 1 nooutput
actual_time = STATS_min

buffer = 0.05

set multiplot layout 2, 2 title sprintf("Final solution at t = Time = %.3f", actual_time)

# First plot: rho
stats datafile index i using 3 nooutput 
expand_range = 0.15*(STATS_max - STATS_min)
rho_min = STATS_min - expand_range - buffer
rho_max = STATS_max + expand_range + buffer
set yrange [rho_min:rho_max]
set ylabel "rho"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:3 with points title "rho"
set yrange[*:*]

# Second plot: vx
stats datafile index i using 4 nooutput 
expand_range = 0.15*(STATS_max - STATS_min)
vx_min = STATS_min - expand_range - buffer
vx_max = STATS_max + expand_range + buffer
set yrange [vx_min:vx_max]
set ylabel "vx"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:4 with points title "vx"
set yrange[*:*]


# Second plot: vy
stats datafile index i using 5 nooutput 
expand_range = 0.15*(STATS_max - STATS_min)
vy_min = STATS_min - expand_range - buffer
vy_max = STATS_max + expand_range + buffer
set yrange [vy_min:vy_max]
set ylabel "vy"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:5 with points title "vy"
set yrange[*:*]


# Third plot: p
stats datafile index i using 6 nooutput  # Retrieve min/max for "v"
expand_range = 0.15*(STATS_max - STATS_min)
p_min = STATS_min - expand_range - buffer
p_max = STATS_max + expand_range + buffer
set yrange [p_min:p_max]
set ylabel "p"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:6 with points title "p"

unset multiplot
set output
