clc, clear all, close all

% Simulation parameters
ts = 0.1;

% State-pace parameters
A = [0.2 1;0 0.7];
B = [1;0];
C = [1 0];
D = 0;
Bd = [0.8;1.5];

r = 1;
d = 0;

Np = 20;
Nu = 3;

z = [C 0;A-eye(length(B)) B]\[r-d;0;0];
xss = z(1:2);
uss = z(3);

Px = zeros(2*Np,2);
P = zeros(Np,2);
Hx = zeros(2*Np,Np);
H = zeros(Np,Np);
for i = 1:Np
    Px((2*i-1):(2*i),:) = A^i;
    P(i,:) = C*A^i;
    Hx((2*i-1):(2*i),1) = A^(i-1)*B;
    H(i,1) = C*A^(i-1)*B;
    for j = 2:Np
        if (i+j) < (Np+2)
            H(i+j-1,j) = H(i,1);
            Hx(2*(i+j-1)-1:2*(i+j-1),j) = Hx((2*i-1):(2*i),1);
        end
    end
end

Q = eye(Np);
R = eye(Nu);
K = (H(:,1:Nu)'*Q*H(:,1:Nu)+R)\(H(:,1:Nu)'*Q*P);

Q2 = eye(2*Np);
R2 = eye(2*Nu);
Kx = (Hx(:,1:2*Nu)'*Q2*Hx(:,1:2*Nu)+R2)\(Hx(:,1:2*Nu)'*Q2*Px);

A_t = [A [0;0];C*A 1];
B_t = [B;C*B];
C_t = [[0 0] 1];

z_t = [C_t 0;A_t-eye(length(B_t)) B_t]\[r-d;0;0;0];
xss_t = z_t(1:2);
uss_t = z_t(3);

l = length(B_t);

Px_t = zeros(l*Np,l);
P_t = zeros(Np,l);
Hx_t = zeros(l*Np,Np);
H_t = zeros(Np,Np);
for i = 1:Np
    Px_t((l*i-2):(l*i),:) = A_t^i;
    P_t(i,:) = C_t*A_t^i;
    Hx_t((l*i-2):(l*i),1) = A_t^(i-1)*B_t;
    H_t(i,1) = C_t*A_t^(i-1)*B_t;
    for j = 2:Np
        if (i+j) < (Np+2)
            H_t(i+j-1,j) = H_t(i,1);
            Hx_t(l*(i+j-1)-2:l*(i+j-1),j) = Hx_t((l*i-2):(l*i),1);
        end
    end
end

Qa = eye(l*Np);
Ra = eye(l*Nu);
Ka = (Hx_t(:,1:l*Nu)'*Qa*Hx_t(:,1:l*Nu)+Ra)\(Hx_t(:,1:l*Nu)'*Qa*Px_t);

%% a)

% Horizon parameters
Np = 20;
Nu = 3;

% Disturbance parameters
d_m = 0;
v = 0;
d = 0.5;

r = 1;
yss = r;

z = [C 0;A-eye(length(B)) B]\[yss-d;0;0];
xss = z(1:2);
uss = z(3);

x = zeros(2,1);

cf = 11;

U = zeros(cf,1);
Y = zeros(cf,1);
for i = 1:cf
%     u = -Ka*x_t;

    u = -K*(x-xss)+uss;

    y = C*x+D*u(1)+d+v;
    x = A*x+B*u(1)+Bd*d_m;

    U(i) = u(1);
    Y(i) = y;
end

ref = ones(cf,1);
time = (0:cf)*ts;

figure
subplot(2,1,1)
stairs(time(1:end-1),ref)
hold on
stairs(time(1:end-1),Y)
axis([0 1 0 1.25])
title('System Output')
xlabel('time (s)')
ylabel('amplitude')
legend('reference','output','Location','southeast')
subplot(2,1,2)
stairs(time(1:end-1),U)
axis([0 1 0 0.5])
title('System Input')
xlabel('time (s)')
ylabel('amplitude')

%% b)

% Horizon parameters
Np = 20;
Nu = 3;

% Disturbance parameters
d_m = 0.8;
v = 0;
d = 0;

r = 1;

Q = eye(Np);
R = eye(Nu);

M = H_t(:,1:Nu)'*Q*H_t(:,1:Nu)+R;

x_t = zeros(3,1);
u = 0;
z = 0;

cf = 11;

U = zeros(cf,1);
Y = zeros(cf,1);
for i = 1:cf    
    f = H_t(:,1:Nu)'*Q*(P_t*x_t-r);
    
    du = -M\f;
    
    z = C_t*x_t;
    x_t = A_t*x_t+B_t*du(1);
    y = z+d+v;

    u = u + du(1);
    
    U(i) = u;
    Y(i) = y;
end

figure
subplot(2,1,1)
stairs(time(1:end-1),ref)
hold on
stairs(time(1:end-1),Y)
axis([0 1 0 1.25])
title('System Output')
xlabel('time (s)')
ylabel('amplitude')
legend('reference','output','Location','southeast')
subplot(2,1,2)
stairs(time(1:end-1),U)
axis([0 1 0 1])
title('System Input')
xlabel('time (s)')
ylabel('amplitude')

%% c)

% Horizon parameters
Np = 20;
Nu = 3;

% Disturbance parameters
d_m = 0;
% v = 0;
% d = 0.5;

a = 0.6;

mu=0;
sigma=0.1;
% dw=sigma*randn(1)+mu;

r = 1;
yss = r;

x = zeros(2,1);

cf = 101;

d_2i = (1+a)*sigma*randn(1)+mu - a*sigma*randn(1)+mu + sigma*randn(1)+mu;
d_i = (1+a)*d(1) - a*sigma*randn(1)+mu + sigma*randn(1)+mu;
U = zeros(cf,1);
Y = zeros(cf,1);
for i = 1:cf
    d = (1+a)*d_i - a*d_2i + sigma*randn(1)+mu;
    d_2i = d_i;
    d_i = d;
    
    z = [C 0;A-eye(length(B)) B]\[yss-d;0;0];
    xss = z(1:2);
    uss = z(3);

    u = -K*(x-xss)+uss;

    v = sigma*randn(1)+mu;

    y = C*x+D*u(1)+d+v;
    x = A*x+B*u(1)+Bd*d_m;

    U(i) = u(1);
    Y(i) = y;
end

ref = ones(cf,1);
time = (0:cf)*ts;

figure
subplot(2,1,1)
stairs(time(1:end-1),ref)
hold on
stairs(time(1:end-1),Y)
title('System Output')
xlabel('time (s)')
ylabel('amplitude')
legend('reference','output','Location','southeast')
subplot(2,1,2)
stairs(time(1:end-1),U)
title('System Input')
xlabel('time (s)')
ylabel('amplitude')