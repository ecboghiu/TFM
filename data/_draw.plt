set xlabel 't'
set ylabel 'r,G/N'
set title 'EPES'

p 'FG_N=1001_m=1_a=5.1_sig=0.08.txt' u 1:2
rep 'FG_N=1001_m=1_a=5.1_sig=0.08.txt' u 1:3:0 lt palette

