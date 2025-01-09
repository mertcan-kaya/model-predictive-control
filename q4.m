clc, clear all, close all

% State-pace parameters
M = [0.8 0.1;0.1 0.9];
N = [0.2;0];
Q = [0.1 0.9];
D = 0;

N = 30;

yss = 1;

QM30 = Q*M^N;

x = zeros(2,1);

cf = 100;
for t = 1:cf
    
    for k = 1:N
        yf(k) = Q*M^k*x+yB(k)*mu;
    end
    
end
