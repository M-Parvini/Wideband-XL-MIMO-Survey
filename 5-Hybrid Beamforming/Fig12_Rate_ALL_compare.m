clc
clear
% close all
Addpaths
tic
rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation code for comparing All the proposed algorithms with and
% without beam squint effects
% For Fig.12(a): SNRlist = -20:5:10; percentage = 10/100; RFlist = 4;
% For Fig.12(b): SNRlist = 10; percentage = 10/100,
%                RFlist=[3 5 7 9 11];
% For Fig.12(c): SNRlist = 10; percentage = [0.05 0.1 0.5 1 5 10]./100;    
%                RFlist = 4;
%%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
NumSim = 200;
SNRlist = -20:5:10;                      % SNR values [Fig.12(a)]
% SNRlist = 10;                              % SNR values [Fig.12(b) & (c)]
% percentage = [0.05 0.1 0.5 1 5 10];        % Subcarrier spacing [Fig.12(c)]
percentage = 10;                           % 10 percent BW/fc
SCSlist = (percentage/100)*8e9/64;
% RFlist = [3 5 7 9 11];                     % RF chains values [Fig.12(b)]
RFlist = 4;                                % RF
AntennaConfig = [256 16];                  % [Nt Nr]

%%% Rate variables
% Digital precoding
DgitalRate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% Manifold optimization
MO_Rate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);
BSAMO_Rate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% Orthogonal matching pursuit
OMP_Rate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);
BSAOMP_Rate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% Virtual precoding
VirtualSubarrayRate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% MCCM precoding
MCCMRate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

%%% Main loop
for RF = 1:length(RFlist)
    for SCS = 1:length(SCSlist)
        for SNR = 1:length(SNRlist)
            Nt = AntennaConfig(1, 1);
            Nr = AntennaConfig(1, 2);
            [OFDMParams, ChanParams, BSParams, UEParams] = ...
            InitializeParams_Rate(SNRlist(SNR), SCSlist(SCS), RFlist(RF), Nt, Nr);
            fprintf('::::::::::::::::::\n')
            fprintf(['SNR =', num2str(SNRlist(SNR)), '\n'])
            parfor SimId = 1:NumSim
                %%% Simulation Cycle
                results = ...
                    mMIMO_OFDM_Compare1(OFDMParams, ChanParams, BSParams, UEParams, SNRlist(SNR));

                % Save rate values
                %%% Digital
                DgitalRate(RF, SCS, SNR, SimId) = real(results.DigitalRate);

                %%% Manifold Optimization
                MO_Rate(RF, SCS, SNR, SimId) = real(results.MORate);
                BSAMO_Rate(RF, SCS, SNR, SimId) = real(results.BSAMORate);

                %%% Orthogonal matching pursuit
                OMP_Rate(RF, SCS, SNR, SimId) = real(results.OMPRate);
                BSAOMP_Rate(RF, SCS, SNR, SimId) = real(results.BSAOMPRate);

                %%% Beam broadening
                VirtualSubarrayRate(RF, SCS, SNR, SimId) = real(results.VirtualRate);

                %%% MCCM
                MCCMRate(RF, SCS, SNR, SimId) = real(results.MCCM);

            end

        end
    end
end


toc
%% Rate plot with SNR [Fig. 12(a)]
figure
hold on

for SCS = 1:length(SCSlist)
    %%% Digital
    txt = 'Digital';
    plot(SNRlist, mean(squeeze(DgitalRate(1, SCS, :, :)), 2), 'DisplayName', txt)

    %%% BSA MO and OMP
    txt = 'BSA-MO';
    plot(SNRlist, mean(squeeze(BSAMO_Rate(1, SCS, :, :)), 2), 'DisplayName', txt)
    txt = 'BSA-OMP';
    plot(SNRlist, mean(squeeze(BSAOMP_Rate(1, SCS, :, :)), 2), 'DisplayName', txt)

    %%% MCCM
    txt = 'MCCM';
    plot(SNRlist, mean(squeeze(MCCMRate(1, SCS, :, :)), 2), 'DisplayName', txt)

    %%% MO
    txt = 'MO';
    plot(SNRlist, mean(squeeze(MO_Rate(1, SCS, :, :)), 2), 'DisplayName', txt)

    %%% OMP
    txt = 'OMP';
    plot(SNRlist, mean(squeeze(OMP_Rate(1, SCS, :, :)), 2), 'DisplayName', txt)

    %%% Beam broadening
    txt = 'Beam broadening';
    plot(SNRlist, mean(squeeze(VirtualSubarrayRate(1, SCS, :, :)), 2), 'DisplayName', txt)
end

legend show
grid
box on
xlabel('SNR [dB]', Interpreter='latex')
ylabel('Spectral efficiency [bps/Hz]', Interpreter='latex')
hl = legend;
set(hl, 'Interpreter','latex')
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis


%% Rate plot with SCS [Fig. 12(b)]
% figure
% hold on
% 
% %%%%%% plots
% %%% Digital
% txt = 'Digital';
% plot(1:length(SCSlist), mean(squeeze(DgitalRate(1, :, 1, :)), 2), 'DisplayName', txt)
% 
% %%% BSA-MO
% txt = 'BSA-MO';
% plot(1:length(SCSlist), mean(squeeze(BSAMO_Rate(1, :, 1, :)), 2), 'DisplayName', txt)
% 
% %%% BSA-OMP
% txt = 'BSA-OMP';
% plot(1:length(SCSlist), mean(squeeze(BSAOMP_Rate(1, :, 1, :)), 2), 'DisplayName', txt)
% 
% %%% MCCM
% txt = 'MCCM';
% plot(1:length(percentage), mean(squeeze(MCCMRate(1, :, 1, :)), 2), 'DisplayName', txt)
% 
% %%% MO
% txt = 'MO';
% plot(1:length(SCSlist), mean(squeeze(MO_Rate(1, :, 1, :)), 2), 'DisplayName', txt)
% 
% %%% OMP
% txt = 'OMP';
% plot(1:length(SCSlist), mean(squeeze(OMP_Rate(1, :, 1, :)), 2), 'DisplayName', txt)
% 
% %%% Beam broadening
% txt = 'Beam broadening';
% plot(1:length(percentage), mean(squeeze(VirtualSubarrayRate(1, :, 1, :)), 2), 'DisplayName', txt)
% 
% 
% legend show
% xlabel('$B/f_c [\%]$', 'Interpreter','latex')
% ylabel('Spectral efficiency [bps/Hz]', 'Interpreter','latex')
% grid on
% box on
% hl = legend;
% set(hl, 'Interpreter','latex')
% xaxisproperties= get(gca, 'XAxis');
% xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
% yaxisproperties= get(gca, 'YAxis');
% yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
% xticks([1:6])
% xticklabels(percentage)



% %% Rate plot with RF chain [Fig. 12(c)]
% figure
% hold on
% 
% %%% plots
% %%% Digital
% txt = 'Digital';
% plot(RFlist, mean(squeeze(DgitalRate(:, 1, 1, :)), 2), 'DisplayName', txt)
% 
% %%% BSA-MO
% txt = 'BSA-MO';
% plot(RFlist, mean(squeeze(BSAMO_Rate(:, 1, 1, :)), 2), 'DisplayName', txt)
% 
% %%% BSA-OMP
% txt = 'BSA-OMP';
% plot(RFlist, mean(squeeze(BSAOMP_Rate(:, 1, 1, :)), 2), 'DisplayName', txt)
% 
% %%% MCCM
% txt = 'MCCM';
% plot(RFlist, mean(squeeze(MCCMRate(:, 1, 1, :)), 2), 'DisplayName', txt)
% 
% %%% MO
% txt = 'MO';
% plot(RFlist, mean(squeeze(MO_Rate(:, 1, 1, :)), 2), 'DisplayName', txt)
% 
% %%% OMP
% txt = 'OMP';
% plot(RFlist, mean(squeeze(OMP_Rate(:, 1, 1, :)), 2), 'DisplayName', txt)
% 
% %%% Beam broadening
% txt = 'Beam broadening';
% plot(RFlist, mean(squeeze(VirtualSubarrayRate(:, 1, 1, :)), 2), 'DisplayName', txt)
% 
% legend show
% grid
% box on
% xlabel('RF chain', 'Interpreter','latex')
% ylabel('Spectral efficiency [bps/Hz]', 'Interpreter','latex')
% grid on
% box on
% hl = legend;
% set(hl, 'Interpreter','latex')
% xaxisproperties= get(gca, 'XAxis');
% xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
% yaxisproperties= get(gca, 'YAxis');
% yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis