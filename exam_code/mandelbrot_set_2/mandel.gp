# png
set title "Mandelbrot set";
set terminal pngcairo size 1366,768 enhanced font 'Verdana,10'
set output 'Mandelbrot.png'

set xrange [-2:2]
set yrange [-2:2]
set size ratio 1

set border linewidth 1.5
# Set line styles to blue (#0060ad) and red (#dd181f) and green (#7cfc00)
set style line 1 \
    linecolor rgb '#0060ad' \
    linetype -1 linewidth 1 \
    #pointtype 7 pointsize 1.5
#notice the following has the same result of above (multi-line vs same-line)
set style line 2 \
    linecolor rgb '#dd181f' \
    linetype -1 linewidth 1 \
    #pointtype 5 pointsize 1.5


#num = system("wc -l points.txt | awk '{ print $1 }'")

set pointsize 0.1

plot "mandelbrot.pgm" using 2:($1==0?$3:NaN) with points pointtype 7 lc rgb "yellow" notitle, \
    "" using 2:($1==1?$3:NaN) with points pointtype 7 lc rgb "black" notitle
