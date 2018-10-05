set ytics  
set y2tics 
set autoscale  y
set autoscale y2
plotname= '.\pippo.txt'
#plot plotname using 1:5 with lines axes x1y1, plotname using 1:8 with lines axes x1y2
plot plotname using 1:6 with lines axes x1y1, plotname using 1:9 with lines axes x1y1, plotname using 1:15 with lines axes x1y2, plotname using 1:11 with lines axes x1y2

#set term png             
#set output "chart.png" 
#replot
#set term x11

 
 
 
 