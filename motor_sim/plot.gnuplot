unset multiplot
set terminal svg
set output "svg.html"
set ytics  
set y2tics 
set autoscale  y
set autoscale y2
set key off
plotname= 'C:\Users\MUNARID\git\motor_sim\motor_sim\pippo.txt'

#set multiplot layout 3,3 spacing -0.2,-0.2

set title 'Romegae, i' offset 0,-1
plot plotname using 't':'Romegae' 			with lines axes x1y1,\
           '' using 't':'Somegae'	 		with lines axes x1y1,\
           '' using 't':'i' 			with lines axes x1y2

set title 'Deltamax, lim, real, phi-\pi/2' offset 0,-1
plot plotname  using 't':'delta_lim.delta_max' 		with lines axes x1y1,\
	'' using 't':'delta_lim.delta_lim' 		with lines axes x1y1,\
	'' using 't':'(Sthetae-Rthetae+2*PI*gve-PI/2)' 	with lines axes x1y1,\
	'' using 't':'gamma.phi'                        with lines axes x1y1

set title 'pow_windmill, windmill.out' offset 0,-1
plot plotname  using 't':'pow_windmill' 		with lines axes x1y1,\
	'' using 't':'windmill.out' 		with lines axes x1y1

set title 'Samplie, Romegae*FLUX' offset 0,-1
plot plotname  using 't':'Samplie' 			with lines axes x1y1,\
	'' using 't':'Romegae*FLUX' 		with lines axes x1y1


#unset multiplot
unset output
#set term png             
#set output "chart.png" 
#replot
#set term x11

 
 
 
 