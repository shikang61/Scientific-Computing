set terminal gif animate delay 8 size 1200, 900 enhanced font 'Verdana,12'
set output './euler_animation/'.ARG1.'.gif'
unset warnings

set xlabel "Position (x)"
set grid

RHO = 3
VX = 4
VY = 5
VZ = 6
P = 7
BX = 8
BY = 9
BZ = 10
BMAG = 11

set xrange [*:*]  # Replace xmin and xmax with appropriate values for your x

datafile = ARG1.'.dat'  # Replace with your actual data file

stats datafile using 1 nooutput
n_time_steps = STATS_blocks
 
do for [i=0:n_time_steps-2] {

    stats datafile index i using 1 nooutput
    actual_time = STATS_min

    set multiplot layout 3, 2 title sprintf("Time Evolution of rho, vx, Bz and |B| (Time = %.6f)", actual_time)
    
    # First plot: rho
    set ylabel "rho"
    plot datafile index i using 2:RHO with points title "rho"

    # Second plot: vx
    set ylabel "vx"
    plot datafile index i using 2:VX with points title "vx"

    # Third plot: P
    set ylabel "P"    
    plot datafile index i using 2:P with points title "P"

    # Third plot: Bx
    # set ylabel "Bx"    
    # plot datafile index i using 2:BX with points title "Bx"

    # Third plot: By
    set ylabel "By"
    plot datafile index i using 2:BY with points title "By"

    # Third plot: Bz
    set ylabel "Bz"
    plot datafile index i using 2:BZ with points title "Bz"

    # Third plot: B_mag
    set ylabel "|B|"
    plot datafile index i using 2:BMAG with points title "|B|"
    
    unset multiplot
}

