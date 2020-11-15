function [M]= evol_matrix_SIR(X,Y,Z,beta,gamma,N)

Linear = zeros(3,1);
Forcing = zeros(3,1);
Nonlinear = zeros(3,1);

Linear(1,1) = 0;
Linear(2,1) = -gamma*Y;
Linear(3,1) = gamma*Y;

Nonlinear(1,1) = -beta*X*Y/N;
Nonlinear(2,1) = beta*X*Y/N;
Nonlinear(3,1) = 0;

M(1,1) = Linear(1,1) + Nonlinear(1,1);
M(2,1) = Linear(2,1) + Nonlinear(2,1);
M(3,1) = Linear(3,1) + Nonlinear(3,1);
