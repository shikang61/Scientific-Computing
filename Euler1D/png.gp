set terminal pngcairo size 1500,600 enhanced font 'Verdana,12'
set output './euler_final/'.ARG1.'.png'
unset warnings

set xlabel "Position (x)"
set grid

set xrange [0:1]  # Replace xmin and xmax with appropriate values for your x

datafile = ARG1.'.dat'  # Replace with your actual data file

stats datafile using 1 nooutput
n_time_steps = STATS_blocks
i = n_time_steps-2

stats datafile index i using 1 nooutput
actual_time = STATS_min

set multiplot layout 1, 3 title "Final solution" title sprintf("Final time (Time = %.6f)", actual_time)

# First plot: rho
stats datafile index i using 3 nooutput  # Retrieve min/max for "rho"
expand_range = 0.15*(STATS_max - STATS_min) + 0.1
# rho_min = STATS_min - expand_range
# rho_max = STATS_max + expand_range
#set yrange [rho_min:rho_max]
set ylabel "rho"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:3 with points title "rho", datafile index i using 2:6 with lines title "rho exact" linecolor "red"

set yrange[*:*]

# Second plot: v
stats datafile index i using 4 nooutput  # Retrieve min/max for "v"
expand_range = 0.15*(STATS_max - STATS_min) + 0.1
# v_min = STATS_min - expand_range
# v_max = STATS_max + expand_range
# set yrange [v_min:v_max]
set ylabel "v"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:4 with points title "v", datafile index i using 2:7 with lines title "v exact" linecolor "red"

set yrange[*:*]

# Third plot: p
stats datafile index i using 5 nooutput  # Retrieve min/max for "v"
expand_range = 0.15*(STATS_max - STATS_min) + 0.1
# p_min = STATS_min - expand_range
# p_max = STATS_max + expand_range
# set yrange [p_min:p_max]
set ylabel "p"
set xzeroaxis linewidth 2.5
plot datafile index i using 2:5 with points title "p", datafile index i using 2:8 with lines title "p exact" linecolor "red"


unset multiplot
set output
