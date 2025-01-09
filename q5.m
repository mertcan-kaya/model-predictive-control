clc, clear all, close all

% Simulation parameters
ts = 0.1;

% State-pace parameters
A = [-1.5, 0.56];
B = [0.9, -0.6];
T = 1;

N = 10;
Np = 10;
Nu = 2;

sigma = 1;
lambda = 1;
alpha = 0.5;

A_hat = [A(1)-1, A(2)-A(1)];

A_f = [1 0 0; A_hat(1) 1 0; A_hat(2) A_hat(1) 1];
A_b = [-A_hat(1) -A_hat(2) 0; A_hat(2) 0 0; 0 0 0];

B_b = [B(2) 0; 0 0];
B_f = [B(1) 0; B(2) B(1)];

H = -A_f\A_b;
P = A_f(1:Nu,1:Nu)\B_b;
Q = A_f(1:Nu,1:Nu)\B_f;

G = Q;
w = zeros(2,1);
sigma_m = eye(2)*sigma;
lambda_m = eye(2)*lambda;

cf = 10;
r = ones(100,1);

y = 0;
y_1 = 0;
y_2 = 0;
du_1 = 0;
du_2 = 0;
du = [0;0];
u_1 = 0;
u_2 = 0;
u = 0;
for t = 1:cf
    e = 0.001*randn();
    
    f = H(1:2,1:2)*[y;y_1]+P*[du_1;du_2];

    w(1) = y;
    for i = 1:1
        w(i+1) = alpha * w(i) + (1 - alpha) * r(i+1);
    end

    du_2 = du_1;
    du_1 = du(1);
    du = (G'*sigma_m*G+lambda_m)\(G'*sigma_m*(w-f));
    
    u_2 = u_1;
    u_1 = u;
    u = u + du(1);

    y_1 = y;
    y = -A(1)*y_1 - A(2)*y_2 + B(1)*u_1 + B(2)*u_2 + e;

    DU(t) = du(1);
    U(t) = u;
    Y(t) = y;
end

figure
subplot(3,1,1)
plot(Y)
subplot(3,1,2)
plot(U)
subplot(3,1,3)
plot(DU)
