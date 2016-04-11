unset xrange
unset yrange
set xrange[-3.95:3.95]
set yrange[-3.95:3.95]
set nokey
unset xzeroaxis
unset yzeroaxis
set xzeroaxis
set yzeroaxis
unset xtics
unset ytics
set xtics axis 2
set ytics axis 2
unset parametric
set parametric
set term postscript eps enhanced color
set output 'DKA_method.eps'
unset size
set size 0.9,0.9
set size square
unset label
unset xlabel
unset ylabel
set label 'Re z'     at first  3.8 , -0.3 
set label 'Im z'     at first  0.1 , 3.7
set label '( 1.0 , 0.0 )'       at first  1.0 , 0.2
set label '( 1.2 , 2.0 )' at first  0.9 , 2.3
set label '( -1.0 , 1.0 )'   at first -1.9 , 1.2
set label '( 1.5 , - 2.0 )'   at first  1.0 , -1.8
set label '( -1.5 , - 1.5 )'at first -2.0 , -1.2
set label '+' at first 0.17,-0.07
set label 'R = 3.1864'at first -3.5,3.5
unset arrow
R = 3.18640337246646
set arrow from 0.24 , -0.1 to R*cos(3.2*pi/2)+0.24 , R*sin(3.2*pi/2)-0.1  

set style line 8 lt 8
plot [0:2*pi] R*cos(t)+0.24,R*sin(t)-0.1,\
     'z_1.dat'with line,'z_2.dat'with line,'z_3.dat'with line, \
     'z_4.dat'with line,'z_5.dat'with line lt 8
set term x11
