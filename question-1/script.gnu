set key left top
set xlabel 'Temperature (K)'
set ylabel 'C_V (Kg.J/K)'

plot 'temp.dat'  w l lt 7 title "Debye's Heat Capacity"

set term eps
set output 'DebyeHC.eps'
replot

unset term 
