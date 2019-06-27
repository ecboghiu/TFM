set nokey
set terminal png size 500,400
set output '_1_sss.png'
p 'FG_N=2000_m=1_a=5_sig=0.01_EDGELIST.txt' u 4:5:($0/2000.0) lt palette
#rep 'FG_N=2000_m=1_a=5_sig=0.08_EDGELIST.txt' u 1:3:0 lt palette

