set ytics  
set y2tics 
set autoscale  y
set autoscale y2
plotname= '.\pippo.txt'
set multiplot layout 3,1
#plot plotname using 1:5 with lines axes x1y1, plotname using 1:8 with lines axes x1y2
plot plotname using 1:6 with lines axes x1y1, plotname using 1:9 with lines axes x1y1, plotname using 1:11 with lines axes x1y2
plot plotname using 1:15 with lines axes x1y1, plotname using 1:16 with lines axes x1y1, plotname using 1:10 with lines axes x1y1
plot plotname using 1:17 with lines axes x1y1
unset multiplot
#set term png             
#set output "chart.png" 
#replot
#set term x11

 
 
 
 