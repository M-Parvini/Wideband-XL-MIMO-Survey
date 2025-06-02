clc
clear
% close all
tic
% rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);
%%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
Chan.LSpeed = physconst('LightSpeed');
fc     = 8e9;
Chan.lambda = Chan.LSpeed/fc;
fs     = 2*800e6;
Nt     = 1024;
K      = 1024;
AoD    = deg2rad([-60 -45 -20 20 45 60]);
delays = [2:0.5:4.5]*1e-5;
gains  = ones(1,6);

L      = length(gains);

Atf  = zeros(Nt, L, K);
Atfc = zeros(Nt, L);

H_f=zeros(Nt,K);
H_fc=zeros(Nt,K);
% Wireless Channel
for k=1:K
    fk = fs/(K)*(k-1-(K-1)/2);
    f=fc+fk;
    %%% NLoS paths
    for j=1:L
        % array responses
        Atfc(:, j, k) = ArrayResponse(AoD(j),Nt,Chan,Chan.LSpeed/fc);
        Atf(:, j, k) = ArrayResponse(AoD(j),Nt,Chan,Chan.LSpeed/f);
        % Channel model (non-beams quint channel)
        H_fc(:,k)=H_fc(:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*Atfc(:, j);
        % Channel model (beam squint channel)
        H_f(:,k)=H_f(:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*Atf(:, j, k);
    end        
end

%%%%% plots
AngDelayHfc = dftmtx(Nt)'*H_fc*conj(dftmtx(K));
AngDelayHf = dftmtx(Nt)'*H_f*conj(dftmtx(K));
% normalization
AngDelayHfc = AngDelayHfc./max(AngDelayHfc(:));
AngDelayHfc = fftshift(AngDelayHfc, 1);

AngDelayHf = AngDelayHf./max(AngDelayHf(:));
AngDelayHf = fftshift(AngDelayHf, 1);

figure
grid on

[X,Y] = meshgrid(1:Nt,1:K);
plott = surf(X,Y,abs(AngDelayHfc).');
set(plott,'LineStyle','none')
xlim([1 Nt])
ylim([1 K])
xlabel('Angular-domain index', Interpreter='latex')
ylabel('Delay-domain index', Interpreter='latex')
zlabel('Normalized magnitude', Interpreter='latex')

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
yaxisproperties= get(gca, 'ZAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
box on

figure
grid on

[X,Y] = meshgrid(1:Nt,1:K);
plott = surf(X,Y,abs(AngDelayHf).');
set(plott,'LineStyle','none')
xlim([1 Nt])
ylim([1 K])
xlabel('Angular-domain index', Interpreter='latex')
ylabel('Delay-domain index', Interpreter='latex')
zlabel('Normalized magnitude', Interpreter='latex')

xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
yaxisproperties= get(gca, 'ZAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
box on