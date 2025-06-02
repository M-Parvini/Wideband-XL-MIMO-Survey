clc
clear
%%%%%%%%%%%%%%%%%%%
theta_ULA = -pi/2:0.001:pi/2;  % AoA/AoD
type = 'Tx';
M = 256;  % # Rx antennas

% UPA Tx side horizental and vertical elements
M_h = sqrt(M);
M_v = sqrt(M);

c = 3e8; % light speed
fc = 300e9;
B = 30e9;

lambda = c/fc;
theta_F_Tx = 0.5;
theta_F_Rx = 0.5;

BFc = B/fc;
nN = [-0.5, 0, 0.5];

%%% Function definition
ULA_Response = @(M,n,b,ang_z,ang_zF) ...
    1/M * abs(sin(((b*n+1)*sin(ang_z)-sin(ang_zF))*M*pi/2) ...
    /sin(((b*n+1)*sin(ang_z)-sin(ang_zF))*pi/2));

%% Normalized array gain UPA
theta_UPA = -pi/2:0.005:pi/2;  % AoA/AoD
phi_UPA = -pi/2:0.005:pi/2;  % AoA/AoD
AG_UPA = zeros(length(nN), length(theta_UPA), length(phi_UPA));
for i = 1:length(nN)
    for j=1:length(theta_UPA)
        for k=1:length(phi_UPA)
            AG_UPA(i, j, k) = ULA_Response(M_h, nN(i), BFc, theta_UPA(j), theta_F_Tx)...
            *ULA_Response(M_v, nN(i), BFc, phi_UPA(k), theta_F_Rx);
        end
    end
end


% plots beam squint UPA
figure
for i = 1:length(nN)
    [X,Y] = meshgrid(theta_UPA, phi_UPA);
    p(i) = surf(X,Y, squeeze(AG_UPA(i, :, :)).');
    % p(i) = polarplot(theta, B_g(i, :), 'DisplayName',txt);
    hold on
end
xlabel('$\theta$ [rad]', Interpreter='latex')
ylabel('$\phi$ [rad]', Interpreter='latex')
zlabel('Normalized array gain', Interpreter='latex')
xlim([theta_F_Tx-0.2 theta_F_Tx+0.2])
ylim([theta_F_Rx-0.2 theta_F_Rx+0.2])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i = 1:3
    A = squeeze(AG_UPA(i, :, :));
    [val, ind] = max(A,[], 2);
    plot(phi_UPA, pow2db(A(:, ind(1))));
    hold on
end
xlabel('$\theta$ or $\phi$ [rad]', Interpreter='latex')
ylabel('Normalized array gain', Interpreter='latex')
legend('$f_1$','$f_c$','$f_K$', Interpreter='latex')
xlim([theta_F_Tx-0.1 theta_F_Tx+0.1])
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
yaxisproperties= get(gca, 'ZAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
box on
grid on

