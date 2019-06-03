set xlabel 't'
set ylabel 'r,G/N'
set title 'EPES'

p 'FG_N=2000_m=1_2_a=20_sig=0.05.txt' u 1:2:0
rep 'FG_N=2000_m=1_2_a=20_sig=0.05.txt' u 1:3:0

