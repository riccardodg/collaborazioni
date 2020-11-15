% rk4st.m - function m-file for Fourth order Runge Kutta evaluation of
% steady state integration of the boundary value problem:

function [x,y,z] = rk4_SIR(dt,x,y,z,beta,gamma,N);

% RK coefficient 0:

X = x;
Y = y;
Z = z;

[M]= evol_matrix_SIR(X,Y,Z,beta,gamma,N);

% RK coefficient 1:

k1 = M(1,1);
l1 = M(2,1);
m1 = M(3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = x + (dt/2)*k1;
Y = y + (dt/2)*l1;
Z = z + (dt/2)*m1;

[M]= evol_matrix_SIR(X,Y,Z,beta,gamma,N);

% RK coefficient 2:

k2 = M(1,1);
l2 = M(2,1);
m2 = M(3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = x + (dt/2)*k2;
Y = y + (dt/2)*l2;
Z = z + (dt/2)*m2;

[M]= evol_matrix_SIR(X,Y,Z,beta,gamma,N);

% RK coefficient 3:

k3 = M(1,1);
l3 = M(2,1);
m3 = M(3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = x + (dt/2)*k3;
Y = y + (dt/2)*l3;
Z = z + (dt/2)*m3;

[M]= evol_matrix_SIR(X,Y,Z,beta,gamma,N);

% RK coefficient 4:

k4 = dt*M(1,1);
l4 = dt*M(2,1);
m4 = dt*M(3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

kappa = k1 + 2*k2 + 2*k3 + k4;
elle = l1 + 2*l2 + 2*l3 + l4;
emme = m1 + 2*m2 + 2*m3 + m4;

X = x + (dt/6)*kappa;
Y = y + (dt/6)*elle;
Z = z + (dt/6)*emme;

x = X;
y = Y;
z = Z;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


