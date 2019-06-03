set xlabel 't'
set ylabel 'r,G/N'
set title 'EPES'

p 'FG_N=500_m=1_a=-19.9_sig=0.05.txt' u 1:2
rep 'FG_N=500_m=1_a=-19.9_sig=0.05.txt' u 1:3:0 lt palette

