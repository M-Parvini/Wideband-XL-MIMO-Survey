clc
clear
%%%%%%%%%%%%%%%%%%%
Antennas = [16, 64, 256];  % # Antennas
Fc = [3.5e9, 8e9 ,300e9];
BW = [0.1e9, 400e6, 30e9];

theta_F = 0.5;
theta = -pi/2:0.00001:pi/2;  % AoA/AoD

for i=1:3
    M = Antennas(i);
    BFc = BW(i)/Fc(i);
    ArrayGains(i,:,:) = ArrayGainULA(M, BFc, theta_F);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i = 1:3
    plot(theta, pow2db(squeeze(ArrayGains(1, i, :))));
    hold on
end
legend('$f_1$','$f_c$','$f_K$', Interpreter='latex')
xlabel('$\theta$ [rad]', Interpreter='latex')
ylabel('Normalized array gain [dB]', Interpreter='latex')
xlim([theta_F-0.1 theta_F+0.1])
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
yaxisproperties= get(gca, 'ZAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i = 1:3
    plot(theta, pow2db(squeeze(ArrayGains(2, i, :))));
    hold on
end
xlabel('$\theta$ [rad]', Interpreter='latex')
ylabel('Normalized array gain [dB]', Interpreter='latex')
xlim([theta_F-0.1 theta_F+0.1])
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
yaxisproperties= get(gca, 'ZAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
box on
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i = 1:3
    plot(theta, pow2db(squeeze(ArrayGains(3, i, :))));
    hold on
end
xlabel('$\theta$ [rad]', Interpreter='latex')
ylabel('Normalized array gain [dB]', Interpreter='latex')
xlim([theta_F-0.1 theta_F+0.1])
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
yaxisproperties= get(gca, 'ZAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
box on
grid on