set terminal gif animate delay 8 size 1500,600 enhanced font 'Verdana,12'
set output './euler_animation/'.ARG1.'.gif'
unset warnings

set xlabel "Position (x)"
set grid

set xrange [*:*]  # Replace xmin and xmax with appropriate values for your x

datafile = ARG1.'.dat'  # Replace with your actual data file

stats datafile using 1 nooutput
n_time_steps = STATS_blocks
 
do for [i=0:n_time_steps-2] {

    stats datafile index i using 1 nooutput
    actual_time = STATS_min

    set multiplot layout 1, 3 title sprintf("Time Evolution of rho, v, and p (Time = %.6f)", actual_time)
    
    # First plot: rho
    set ylabel "rho"
    plot datafile index i using 2:3 with points title "rho", datafile index i using 2:6 with lines title "rho exact" linecolor "red"

    # Second plot: vx
    set ylabel "vx"
    plot datafile index i using 2:4 with points title "v", datafile index i using 2:7 with lines title "v exact" linecolor "red"

    # Third plot: p
    set ylabel "p"
    plot datafile index i using 2:5 with points title "p", datafile index i using 2:8 with lines title "p exact" linecolor "red"
    
    unset multiplot
}

