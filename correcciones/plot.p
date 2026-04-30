set terminal pngcairo size 1200,750 enhanced font 'Verdana,23'

set output 'lyapunov-r2.png'
set title "Exponentes Característicos de Lyapunov"
set xlabel "Parámetro r_2"
set ylabel "Exponentes"
set grid
set yrange[-0.8:0.1]
set key left bottom
plot "Espectro1.dat" u 1:2 w l lt 1 lw 2 title "{/Symbol l}_1", \
     "Espectro1.dat" u 1:3 w l lt 6 lw 2 title "{/Symbol l}_2", \
     "Espectro1.dat" u 1:4 w l lt 2 lw 2 title "{/Symbol l}_3"

set output 'lyapunov-r2-positivo.png'
set title "Exponente Positivo de Lyapunov"
set xlabel "Parámetro r_2"
set yrange[-0.05:0.05]
set xrange[0.4:0.8]
set key right top
plot "Espectro1.dat" u 1:($2>0?$2:NaN) \
     w l lt 7 lw 2 title "{/Symbol l}_1 > 0"