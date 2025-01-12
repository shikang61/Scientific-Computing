set terminal gif animate delay 8 size 1600,600 enhanced font 'Verdana,12'
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
    set arrow from graph 0, first 0 to graph 1, first 0 nohead lc "black" lw 1
    set arrow from graph 0, second 0 to graph 1, second 0 nohead lc "black" lw 1
    set arrow from 0.5, graph 0 to 0.5, graph 1 nohead lc "black" lw 1 dt 2
    set yrange[*:*]
    # First plot: rho
    set ylabel "rho"
    set y2label "phi"
    set y2tics
    set ytics nomirror
    set y2range[-0.5:0.5]
    # stats datafile using 3 nooutput  # Retrieve min/max for "rho"
    # expand_range = 0.15*(STATS_max - STATS_min) + 0.1
    # rho_min = STATS_min - expand_range
    # rho_max = STATS_max + expand_range
    # set yrange [rho_min:rho_max]
    # plot datafile index i using 2:3 with points title "rho", datafile index i using 2:6 with lines title "rho exact" linecolor "red"
    plot datafile index i using 2:($9 < 0 ? $3 : 1/0) with points title "rho", \
         datafile index i using 2:($9 > 0 ? $6 : 1/0) with points title "rho2", \
         datafile index i using 2:10 with lines title "exact" linecolor "red",\
         datafile index i using 2:9 axes x1y2 with lines title "phi" linecolor "green"
    set yrange[*:*]

    # Second plot: vx
    set ylabel "vx"
    set y2label "phi"
    set y2tics
    set ytics nomirror
    set y2range[-0.5:0.5]
    # stats datafile using 4 nooutput  # Retrieve min/max for "v"
    # expand_range = 0.15*(STATS_max - STATS_min) + 0.1 
    # v_min = STATS_min - expand_range
    # v_max = STATS_max + expand_range
    # set yrange [v_min:v_max]
    # plot datafile index i using 2:4 with points title "v", datafile index i using 2:7 with lines title "v exact" linecolor "red"
    plot datafile index i using 2:($9 < 0 ? $4 : 1/0) with points title "v", \
         datafile index i using 2:($9 > 0 ? $7 : 1/0) with points title "v2", \
         datafile index i using 2:11 with lines title "exact" linecolor "red",\
         datafile index i using 2:9 axes x1y2 with lines title "phi" linecolor "green" 
     #     datafile index i using 2:11 with points title "v2", \
    set yrange[*:*]

    # Third plot: p
    set ylabel "p"
    set y2label "phi"
    set y2tics
    set ytics nomirror
    set y2range[-0.5:0.5]
    # stats datafile using 5 nooutput  # Retrieve min/max for "v"
    # expand_range = 0.15*(STATS_max - STATS_min) + 0.1
    # p_min = STATS_min - expand_range
    # p_max = STATS_max + expand_range
    # set yrange [p_min:p_max]
    # plot datafile index i using 2:5 with points title "p", datafile index i using 2:8 with lines title "p exact" linecolor "red"
    plot datafile index i using 2:($9 < 0 ? $5 : 1/0) with points title "p", \
         datafile index i using 2:($9 > 0 ? $8 : 1/0) with points title "p2",\
         datafile index i using 2:12 with lines title "exact" linecolor "red",\
         datafile index i using 2:9 axes x1y2 with lines title "phi" linecolor "green" 
     #     datafile index i using 2:12 with points title "p2"

    
    unset multiplot
}

