clear all
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set parameters

[beta,gamma,N]=Parameters_SIR();
dt = 10^(-2);   %time-step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial condition

[x,y,z]=IC_SIR();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Integration in time

t = 0;
tmax = 200;
[plotgap,nplots,tplot,time]=plot_intervals(dt,tmax);

%Pre-allocation arrays

XX = zeros(nplots,1);
YY = zeros(nplots,1);
ZZ = zeros(nplots,1);

for i=1:nplots
    
    XX(i,1) = x;
    YY(i,1) = y;
    ZZ(i,1) = z;
    
    X = x;
    Y = y;
    Z = z;
    
    %[M,PD]= evol_matrix_3rd(A0,A1,B1,L,c,epsil,mu0,mu1,x);
    %g = a0 + a1*cos(x) + b1*sin(x);
    %Phi(i,:) = g;

    for n=1:plotgap
        [x,y,z] = rk4_SIR(dt,x,y,z,beta,gamma,N);
        t = t+dt;
    end
end

clear x y z X Y Z i dt n t tmax nplots plotgap tplot

figure
plot(time,XX,'b')
hold on
plot(time,YY,'r')
plot(time,ZZ,'g')
hold off

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Loop on \mu_{1}

Dmu1 = 0.005;

A_zero = zeros(nplots,1);
A_uno = zeros(nplots,1);
B_uno = zeros(nplots,1);
Phi = zeros(nplots,N);
Phi_dot = zeros(nplots,N);
g = zeros(1,N);

for j=1:jmax
    
    j
    
    if (j<=1000)
        tmax = 500;
        [plotgap,nplots,tplot,time]=plot_intervals(dt,tmax);
    elseif ((j>1000)&(j<=1700))
        tmax = 800;
        [plotgap,nplots,tplot,time]=plot_intervals(dt,tmax);
    else 
        tmax = 2000;
        [plotgap,nplots,tplot,time]=plot_intervals(dt,tmax);
    end
    
    t = 0;
    time = [0:tplot:(tmax-tplot)]';
    mu1 = mu1+Dmu1;
    
    for i=1:nplots
        
        A_zero(i,1) = a0;
        A_uno(i,1) = a1;
        B_uno(i,1) = b1;
        
        A0 = a0;
        A1 = a1;
        B1 = b1;
        
        [M,PD]= evol_matrix_SIR(A0,A1,B1,L,c,epsil,mu0,mu1,x);
        g = a0 + a1*cos(x) + b1*sin(x);
        Phi(i,:) = g;
        for n=1:plotgap
            [a0,a1,b1] = rk4_SIR(dt,a0,a1,b1,L,c,epsil,mu0,mu1,x);
            t = t+dt;
        end
    end
    
    %filename = num2str(j,'Data_TW1_%d.mat');
    %if (j>200)
        %save(['/Users/francescosarnari/Documents/MATLAB/ODE_TW1_mu0_10/' filename],'x','A0','A1','B1','A_zero','A_uno','B_uno','time','Phi','L','c','epsil','mu0','mu1');
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return

%to indicize output files
%for K = 1:12
%  filename = num2str(K,'Data%d.txt');
%  save(filename, 'x','y');
%end

%to save output files in a directory different from the default one

%filename = num2str(l,'data_TS_TW1_%d.mat');
%save(['E:\TS_PDE_per_23_04_2008\mu0_0p5\TW1\' filename]);

