function [plotgap,nplots,tplot,time]=plot_intervals(dt,tmax)

tplot = tmax/1000;
%tplot = tmax/dt;

plotgap = round(tplot/dt);

nplots = round(tmax/tplot);

time = [0:tplot:(tmax-tplot)]';

