base_name = "slice2D100_"
file_type = ".dat"
set xlabel "Normalized Position (x/L)"
set ylabel "Ux_{LBM}"
set yrange [-0.00075:0.00066]
set xrange [0:1]
set grid
do for [i=0:50] {
    set terminal png size 800,600
    set output sprintf("frame_%03d.png", i)
    filename = base_name . i . file_type
    set title sprintf("Ux_{LBM} - Timestep %d", i)
    plot filename using 1:3 with lines linewidth 2 title sprintf("t = %d", i)
}
set output
