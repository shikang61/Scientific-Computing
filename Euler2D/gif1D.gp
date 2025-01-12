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

    set multiplot layout 2, 2 title sprintf("Time Evolution of rho, v, and p (Time = %.3f)", actual_time)
    
    # First plot: rho
    set ylabel "rho"
    plot datafile index i using 2:3 with points title "rho"

    # Second plot: vx
    set ylabel "vx"
    plot datafile index i using 2:4 with points title "v"

    # Second plot: vy
    set ylabel "vy"
    plot datafile index i using 2:5 with points title "vy"

    # Third plot: p
    set ylabel "p"
    plot datafile index i using 2:6 with points title "p"
    
    unset multiplot
}

