set terminal gif animate delay 12 size 1400, 1000 enhanced font 'Verdana,12'
set output './euler_animation/'.ARG1.'.gif'

set xrange [*:*]    # Automatically scale to x data range
set yrange [*:*]    # Automatically scale to y data range

set xlabel "X-axis (x)"
set ylabel "Y-axis (y)"
set palette rgb 33,13,10   # Optional: set color palette (modify to suit preference)

datafile = ARG1.'.dat'  # Replace with your actual data file

buffer = 0.05
stats datafile using 4 nooutput 
expand_range = 0.15*(STATS_max - STATS_min)
rho_min = STATS_min - expand_range - buffer
rho_max = STATS_max + expand_range + buffer

stats datafile using 5 nooutput 
expand_range = 0.15*(STATS_max - STATS_min)
vx_min = STATS_min - expand_range - buffer
vx_max = STATS_max + expand_range + buffer


stats datafile using 6 nooutput 
expand_range = 0.15*(STATS_max - STATS_min)
vy_min = STATS_min - expand_range - buffer
vy_max = STATS_max + expand_range + buffer


stats datafile using 7 nooutput 
expand_range = 0.15*(STATS_max - STATS_min)
p_min = STATS_min - expand_range - buffer
p_max = STATS_max + expand_range + buffer


stats datafile using 1 nooutput
n_time_steps = STATS_blocks
 
do for [i=0:n_time_steps-2] {

    stats datafile index i using 1 nooutput
    actual_time = STATS_min

    set multiplot layout 2, 2 title sprintf("Time Evolution of rho, v, and p (Time = %.3f)", actual_time)
    
    # First plot: rho
    set zlabel "Density"
    set cblabel "rho"
    unset key           
    set view map        
    set style data pm3d
    set title "Density"
    set cbrange [rho_min:rho_max]
    splot datafile index i using 2:3:4 

    # First plot: vx
    set zlabel "vx"
    set cblabel "vx"
    unset key
    set view map        
    set style data pm3d
    # set hidden3d
    set title "vx"
    set cbrange [vx_min:vx_max]
    splot datafile index i using 2:3:5 

    # First plot: vy
    set zlabel "vy"
    set cblabel "vy"
    unset key
    set view map        
    set style data pm3d
    # set hidden3d
    set title "vy"
    set cbrange [vy_min:vy_max]
    splot datafile index i using 2:3:6

    # First plot: p
    set zlabel "Pressure"
    set cblabel "p"
    unset key
    set view map
    set style data pm3d
    set cbrange [p_min:p_max]
    set title "Pressure"
    splot datafile index i using 2:3:7 

    unset multiplot
}

