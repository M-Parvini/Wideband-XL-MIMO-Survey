clc,
clear,

Nt = 128;
Nr = 1;
K = 512;

fc = 8e9;
B  = 400e6;

freq_range = fc-B*(K-1)/(2*K):B/K:fc+B*(K-1)/(2*K);
f1 = freq_range(1);
fc = freq_range(K/2);
fK = freq_range(K);


theta_range = 0.4:0.001:0.6;
theta_l = 0.5;

%%%%%%% Near-field %%%%%%%%
r0 = 20;
nn = -(Nt-1)/2:1:(Nt-1)/2;
c = 3e8;
d = (c/fc)/2;
r = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sin(theta_l));
nf_arr_res_fc = exp(1j*2*pi*fc*(r - r0)/c)/sqrt(Nt);
%nf_arr_res_fc = 1/sqrt(Nt)*(2*pi*fc*(r - r0)/c);

nf_array_gain_f1 = zeros(length(theta_range),1);
nf_array_gain_fc = zeros(length(theta_range),1);
nf_array_gain_fK = zeros(length(theta_range),1);

%%%% With perfect distance knowledge %%%%
for i=1:length(theta_range)
    theta = theta_range(1,i);
    r_n = sqrt(r0^2 + (nn*d).^2 - 2*r0*nn*d*sin(theta));
    nf_arr_res_fc_H = (exp(1j*2*pi*fc*(r_n - r0)/c)/sqrt(Nt))';
    nf_arr_res_f1_H = (exp(1j*2*pi*f1*(r_n - r0)/c)/sqrt(Nt))';
    nf_arr_res_fK_H = (exp(1j*2*pi*fK*(r_n - r0)/c)/sqrt(Nt))';
    nf_array_gain_f1(i,1) = abs(nf_arr_res_fc*nf_arr_res_f1_H);
    nf_array_gain_fc(i,1) = abs(nf_arr_res_fc*nf_arr_res_fc_H);
    nf_array_gain_fK(i,1) = abs(nf_arr_res_fc*nf_arr_res_fK_H);
end

figure,
plot(theta_range,pow2db(nf_array_gain_f1));
hold on,
plot(theta_range,pow2db(nf_array_gain_fc));
hold on,
plot(theta_range,pow2db(nf_array_gain_fK));
% ylim([0.001 1])
legend('$f_1$','$f_c$','$f_K$', Interpreter="latex");
xlabel("$\theta$ [rad]", Interpreter="latex");
ylabel("Normalized array gain [dB]", Interpreter="latex");
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
yaxisproperties= get(gca, 'ZAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
box on
grid on
%%%% With imperfect distance knowledge %%%%
dist_range = 5:0.1:35;
nf_array_gain_f1_dist = zeros(length(dist_range),1);
nf_array_gain_fc_dist = zeros(length(dist_range),1);
nf_array_gain_fK_dist = zeros(length(dist_range),1);
for i=1:length(dist_range)
    r_prime = dist_range(1,i);
    r_dist = sqrt(r_prime^2 + (nn*d).^2 - 2*r_prime*nn*d*sin(theta_l));
    nf_arr_res_fc_H_dist = (exp(1j*2*pi*fc*(r_dist - r0)/c)/sqrt(Nt))';
    nf_arr_res_f1_H_dist = (exp(1j*2*pi*f1*(r_dist - r0)/c)/sqrt(Nt))';
    nf_arr_res_fK_H_dist = (exp(1j*2*pi*fK*(r_dist - r0)/c)/sqrt(Nt))';
    nf_array_gain_f1_dist(i,1) = abs(nf_arr_res_fc*nf_arr_res_f1_H_dist);
    nf_array_gain_fc_dist(i,1) = abs(nf_arr_res_fc*nf_arr_res_fc_H_dist);
    nf_array_gain_fK_dist(i,1) = abs(nf_arr_res_fc*nf_arr_res_fK_H_dist);
end

figure,

plot(dist_range,pow2db(nf_array_gain_f1_dist));
hold on,
plot(dist_range,pow2db(nf_array_gain_fc_dist));
hold on,
plot(dist_range,pow2db(nf_array_gain_fK_dist));
%ylim([0.1 1])
legend('$f_1$','$f_c$','$f_K$', Interpreter="latex");
xlabel("Distance (m)", Interpreter="latex");
ylabel("Normalized array gain [dB]", Interpreter="latex");

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
yaxisproperties= get(gca, 'ZAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
box on
grid on