clc
clear
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = -pi/2:0.0001:pi/2;  % AoA/AoD
dist = 1:0.05:80;
M = 32;  % # antennas
c = 3e8; % light speed
fc = 8e9;
B = 1e9;
f_ml = fc-B/2; lambda_l = c/f_ml;
f_mu = fc+B/2; lambda_u = c/f_mu;
lambda_c = c/fc;
lambda_m = [lambda_l,lambda_c,lambda_u];
theta_F = deg2rad(30);
dist_F = 30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
AG = @(M,n,lambda_c,lambda_m,theta,dist,theta_F,dist_F) ...
0.5/M*exp(1j*n*0.5*lambda_c*(2*pi/lambda_m*sin(theta) - 2*pi/lambda_c*sin(theta_F))...
-1j*n^2*0.25*lambda_c^2*(2*pi/lambda_m*((1-sin(theta).^2)./(2*dist)) - ...
2*pi/lambda_c*((1-sin(theta_F)^2)./(2*dist_F))));

AG_theta = zeros(3, length(theta));
AG_dist = zeros(3, length(dist));
%%%%%%

for i=1:3
    for j=-M:M
        AG_theta(i,:)=abs(AG_theta(i,:)+ ...
        AG(M,j,lambda_c,lambda_m(i),theta,dist_F,theta_F,dist_F));
    end
end


for i=1:3
    for j=-M:M
        AG_dist(i,:)=abs(AG_dist(i,:)+ ...
        AG(M,j,lambda_c,lambda_m(i),theta_F,dist,theta_F,dist_F));
    end
end

%%%% plots
figure
hold on
for i=1:3
    plot(theta, AG_theta(i,:))
end
legend('f1', 'fc', 'f2')
grid

figure
hold on
for i=1:3
    plot(dist, AG_dist(i,:))
end
legend('f1', 'fc', 'f2')

grid