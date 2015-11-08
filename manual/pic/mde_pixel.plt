set terminal pdf color dashed transparent linewidth 2 fontscale 0.5 size 14cm,5cm 
set output 'MDE_pixel.pdf'
 set linetype 10 lt 2 pt 2 ps 0.8 lc rgbcolor "black" lw 1.5; set linetype 11 lt 2 pt 3 ps 0.8 lc rgbcolor "red" lw 1.5; set linetype 12 lt 2 pt 4 ps 0.8 lc rgbcolor "blue" lw 1.5; set linetype 13 lt 2 pt 5 ps 0.8 lc rgbcolor "web-green" lw 1.5; set linetype 14 lt 2 pt 6 ps 0.8 lc rgbcolor "dark-magenta" lw 1.5; set linetype 15 lt 2 pt 7 ps 0.8 lc rgbcolor "orange" lw 1.5; set linetype 16 lt 2 pt 8 ps 0.8 lc rgbcolor "dark-blue" lw 1.5; set linetype 17 lt 2 pt 9 ps 0.8 lc rgbcolor "seagreen" lw 1.5; set linetype 18 lt 2 pt 10 ps 0.8 lc rgbcolor "steelblue" lw 1.5; set linetype 19 lt 2 pt 11 ps 0.8 lc rgbcolor "purple" lw 1; set linetype 1 lt 1 pt 2 ps 0.8 lc rgbcolor "red" lw 1.5; set linetype 2 lt 1 pt 3 ps 0.8 lc rgbcolor "blue" lw 1.5; set linetype 3 lt 1 pt 4 ps 0.8 lc rgbcolor "web-green" lw 1.5; set linetype 4 lt 1 pt 5 ps 0.8 lc rgbcolor "dark-magenta" lw 1.5; set linetype 5 lt 1 pt 6 ps 0.8 lc rgbcolor "orange" lw 1.5; set linetype 6 lt 1 pt 7 ps 0.8 lc rgbcolor "dark-blue" lw 1.5; set linetype 7 lt 1 pt 8 ps 0.8 lc rgbcolor "seagreen" lw 1.5; set linetype 8 lt 1 pt 9 ps 0.8 lc rgbcolor "purple" lw 1.5; set style line 11 lt 1 lc rgb "gray80" lw 0.7; set style line 12 lt 1 lc rgb "gray80" lw 0.4; set border lw 1.2; set grid xtics ytics ls 11, ls 12; set grid mxtics mytics ls 11, ls 12; set key inside left top Left reverse samplen 2 spacing 1.35 font ',7'; set xtics border nomirror center out scale 0.5,0.2 offset 0,0.2; set ytics border nomirror right out scale 0.5,0.2 offset 0.2,0; set palette file './data/default.pal' using 1:2:3:4; unset colorbox; set xzeroaxis lt -1 lw 1.2; set yzeroaxis lt -1 lw 1.2; histogramBin(x, width, Min) = width*(floor((x-Min)/width)+0.5) + Min; set format x '%g'; set format y '%g'; set ylabel offset 1.5,0; histogramBinMin0(x,width)=histogramBin(x,width,0); 
      set xlabel 'positions x [micrometers]'
      set ylabel 'MDE(x) [A.U.]'
      set samples 500
      #set palette rgb 23,28,3;
      unset colorbox
      set key outside 
      MDE(x,w,a)=(erf((a-2.0*x)/(sqrt(2.0)*w))+erf((a+2.0*x)/(sqrt(2.0)*w)))/(2.0*erf(a/(sqrt(2)*w)))
      G(x,w)=exp(-2.0*x*x/w/w)
      plot [-3.5:3.5][0:1.1] MDE(x,0.5,0.5) title 'SPIM MDE: a=500nm' with lines lc palette frac 0.0, \
                  MDE(x,0.5,1) title 'SPIM MDE: a=1000nm' with lines lc palette frac 0.3, \
                  MDE(x,0.5,2) title 'SPIM MDE: a=2000nm' with lines lc palette frac 0.6, \
                  MDE(x,0.5,4) title 'SPIM MDE: a=4000nm' with lines lc palette frac 1, \
                  G(x,0.5) title 'Gaussian' with lines lc rgbcolor "red"
    
