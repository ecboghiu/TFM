set xlabel "t"
set ylabel "log(Degree)"
set zlabel "log(Count)"

#set hidden3d
#set xyplane 0
#set dgrid3d splines


#set palette rgb 21,22,23

splot "hist3D.txt" u ($1):($2):($3):($1) palette notitle
#splot "hist3D_2000_20_0.05.txt"