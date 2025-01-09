clc, clear all, close all

% Simulation parameters
ts = 0.1;
cf = 40;

% Horizon parameters
N = 30;
Np = 10;
Nu = 5;

% Control parameters
alpha = [0.7;0];
lambda = [1;0.1];
sigma = 1;

% System parameters
a = 0.8351;
b = 0.2713;

% impulse response model
h_i = zeros(N,1);
H = zeros(N,N);
for ir = 1:N
    if ir < 3
        h_i(ir) = 0;
    else
        h_i(ir) = b*a^((ir-1)-2);
    end
    for ic = 1:N
        if ir >= ic
            H(ir,ic) = h_i(ir - ic + 1);
        end
    end
end

% step response model
g_i = zeros(N,1);
G = zeros(N,N);
for ir = 1:N
    for ic = 1:ir
        g_i(ir) = g_i(ir) + h_i(ic);
    end
    for ic = 1:N
        if ir >= ic
            G(ir,ic) = g_i(ir - ic + 1);
        end
    end
end

dG = zeros(Np,N-Np);
for ir = 1:Np
    for ic = 1:N-Np
        if ic <= N-Np
            dG(ir,ic) = g_i(ir + ic)-g_i(ic);
        else
            dG(ir,ic) = 0;
        end
    end
end

U = zeros(N,length(alpha));
Y = zeros(N,length(alpha));
for j = 1:length(alpha)
    r = ones(100,1);

%     u_fwd = zeros(N,N);
%     u_bck = zeros(N,N);
%     du_fwd = zeros(N,N);
    du_bck = zeros(N,N);

    for t = 1:N
        for k = 1:N
%             u_fwd(t,k) = r(t+k-1);

%             if t-k > 0
%                 u_bck(t,k) = r(t-k);
%             else
%                 u_bck(t,k) = 0;
%             end

%             if t+k-2 > 0            
%                 du_fwd(t,k) = r(t+k-1)-r(t+k-2);
%             elseif t+k-2 == 0 
%                 du_fwd(t,k) = r(t+k)-0;
%             else
%                 du_fwd(t,k) = 0;
%             end

            if t-k-1 > 0
                du_bck(t,k) = r(t-k)-r(t-k-1);
            elseif t-k-1 == 0
                du_bck(t,k) = r(t-k)-0;
            else
                du_bck(t,k) = 0;
            end
        end

    end

    lambda_m = eye(Nu)*lambda(j);
    sigma_m = eye(Np)*sigma;

    L = ones(Np,1);
    G_m = G(1:Np,1:Nu);
    w = zeros(Np,1);

    u_v = zeros(3,1);
    y_c = 0;
    u_c = 0;

    du_v = du_bck(1,1:N-Np)';

    for t = 1:cf
        f = dG*du_v+L*y_c;

        w(1) = y_c;
        for i = 1:Np-1
            w(i+1) = alpha(j) * w(i) + (1 - alpha(j)) * r(i+1);
        end

        du = (G_m'*sigma_m*G_m+lambda_m)\(G_m'*sigma_m*(w-f));

        du_v = [du(1);du_v(1:end-1)];
        u_c = u_c + du(1);    
        u_v = [u_v(2:3);u_c];
        
        % y(t) = 0.8351*y(t-1) + 0.2713*u(t-3)
        y_c = 0.8351*y_c + 0.2713*u_v(1);

        U(t,j) = u_c;
        Y(t,j) = y_c;
    end
end

figure
time = (0:cf)*ts;

subplot(2,1,1)
stairs(time,r(1:cf+1),'--')
% line([0 0], [0 1.25],'Color','black','LineStyle',':')
hold on
stairs(time,[[0 0];Y(1:end,:)])
xline(0);
axis([-0.5 4 0 1.25])
title('System Output')
xlabel('time (s)')
ylabel('amplitude')
legend('unit step input','\alpha = 0.7 & \lambda = 1','\alpha = 0 & \lambda = 0.1','Location','southeast')

subplot(2,1,2)
stairs(time(1:end-1),U(:,1),'Color',[0.8500 0.3250 0.0980])
hold on
stairs(time(1:end-1),U(:,2),'Color',[0.9290 0.6940 0.1250])
xline(0);
axis([-0.5 4 0 2])
title('System Input')
xlabel('time (s)')
ylabel('amplitude')
legend('\alpha = 0.7 & \lambda = 1','\alpha = 0 & \lambda = 0.1','Location','northeast')
