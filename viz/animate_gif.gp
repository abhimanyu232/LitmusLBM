set terminal gif animate optimize delay 10 size 800,600
in_basename = "slice2D100_"
in_file_type = ".dat"
set output "ux_lbm_animation.gif"
set xlabel "Normalized Position (x/L)"
set ylabel "Ux_{LBM}"
set yrange [-0.00075:0.00066]
set xrange [0:1]
set grid
do for [i=0:50] {
    filename = sprintf("%s%d%s", in_basename, i, in_file_type)
    set title sprintf("Ux_{LBM} — Timestep %d", i)
    plot filename using 1:3 with lines linewidth 2 title sprintf("t = %d", i)
}
set output
