set terminal cairolatex pdf color dashed transparent linewidth 2 fontscale 0.5 size 8cm,6cm 
set output 'manual-gnuplottex-fig7.tex'
 set linetype 10 lt 2 pt 2 ps 0.8 lc rgbcolor "black" lw 1.5; set linetype 11 lt 2 pt 3 ps 0.8 lc rgbcolor "red" lw 1.5; set linetype 12 lt 2 pt 4 ps 0.8 lc rgbcolor "blue" lw 1.5; set linetype 13 lt 2 pt 5 ps 0.8 lc rgbcolor "web-green" lw 1.5; set linetype 14 lt 2 pt 6 ps 0.8 lc rgbcolor "dark-magenta" lw 1.5; set linetype 15 lt 2 pt 7 ps 0.8 lc rgbcolor "orange" lw 1.5; set linetype 16 lt 2 pt 8 ps 0.8 lc rgbcolor "dark-blue" lw 1.5; set linetype 17 lt 2 pt 9 ps 0.8 lc rgbcolor "seagreen" lw 1.5; set linetype 18 lt 2 pt 10 ps 0.8 lc rgbcolor "steelblue" lw 1.5; set linetype 19 lt 2 pt 11 ps 0.8 lc rgbcolor "purple" lw 1; set linetype 1 lt 1 pt 2 ps 0.8 lc rgbcolor "red" lw 1.5; set linetype 2 lt 1 pt 3 ps 0.8 lc rgbcolor "blue" lw 1.5; set linetype 3 lt 1 pt 4 ps 0.8 lc rgbcolor "web-green" lw 1.5; set linetype 4 lt 1 pt 5 ps 0.8 lc rgbcolor "dark-magenta" lw 1.5; set linetype 5 lt 1 pt 6 ps 0.8 lc rgbcolor "orange" lw 1.5; set linetype 6 lt 1 pt 7 ps 0.8 lc rgbcolor "dark-blue" lw 1.5; set linetype 7 lt 1 pt 8 ps 0.8 lc rgbcolor "seagreen" lw 1.5; set linetype 8 lt 1 pt 9 ps 0.8 lc rgbcolor "purple" lw 1.5; set style line 11 lt 1 lc rgb "gray80" lw 0.7; set style line 12 lt 1 lc rgb "gray80" lw 0.4; set border lw 1.2; set palette file 'default.pal' using 1:2:3:4; set grid xtics ytics ls 11, ls 12; set grid mxtics mytics ls 11, ls 12; set key inside left top Left reverse samplen 2 spacing 1.35 font ',7'; set xtics border nomirror center out scale 0.5,0.2 offset 0,0.2; set ytics border nomirror right out scale 0.5,0.2 offset 0.2,0; unset colorbox; set xzeroaxis lt -1 lw 1.2; set yzeroaxis lt -1 lw 1.2; histogramBin(x, width, Min) = width*(floor((x-Min)/width)+0.5) + Min; set format x '\fsfn $%g$'; set format y '\fsfn $%g$'; set ylabel offset 1.5,0; histogramBinMin0(x,width)=histogramBin(x,width,0); 
  set format x '\footnotesize $10^{%T}$';
        N=10
        N0=0
        N1=N*0.2
        N2=N*0.4
        N3=N*0.6
        N4=N*0.8
        N5=N
        w=0.5
        z=1.2
        a=0.4
        c=10
        D2=1
        D1=100
        Veff=sqrt(pi)*a*a*z/(erf(a/w)+w/sqrt(pi)/a*(exp(-a*a/w/w)-1.0))**2
        Aeff=a*a/(erf(a/w)+w/sqrt(pi)/a*(exp(-a*a/w/w)-1.0))**2
        tD(D)=Aeff/4.0/D
        sqr(x)=x*x
        gg(tau,D,N)=N/Veff/(sqrt(pi)*a*a*z)*sqr(erf(a/sqrt(4*D*tau+w*w))+sqrt(4*D*tau+w*w)/a/sqrt(pi)*(exp(-a*a/(4*D*tau+w*w))-1)) /sqrt(1+4.0*D*tau/z/z)
        g(tau, D1, D2, N1, N2)=(gg(tau,D1,N1)+gg(tau,D2,N2))/sqr(N1/Veff+N2/Veff)

        set logscale x
        set xlabel 'lag time $\tau\unitb{s}$'
        set ylabel 'correlation amplitude $g_\gamma(\tau)$'

        set key inside right top vertical  opaque  maxrows 6
        set label 1 '$D_1=100\unit{\mu m^2/s}$' at 3e-5,0.02 left front
        set label 2 '$D_2=1\unit{\mu m^2/s}$' at 3e-5,0.01 left front

        plot [1e-5:10][0:0.105] g(x, D1,   D2, N0, N-N0) title '\footnotesize $\rho_1=0\%$' with lines lc palette frac 0,  \
                                g(x, D1,   D2, N1, N-N1) title '\footnotesize $\rho_1=20\%$' with lines lc palette frac 0.2,  \
                                g(x, D1,   D2, N2, N-N2) title '\footnotesize $\rho_1=40\%$' with lines lc palette frac 0.4,  \
                                g(x, D1,   D2, N3, N-N3) title '\footnotesize $\rho_1=60\%$' with lines lc palette frac 0.6,  \
                                g(x, D1,   D2, N4, N-N4) title '\footnotesize $\rho_1=80\%$' with lines lc palette frac 0.8,  \
                                g(x, D1,   D2, N5, N-N5) title '\footnotesize $\rho_1=100\%$' with lines lc palette frac 0.9
