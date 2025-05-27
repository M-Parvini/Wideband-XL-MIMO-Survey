clc
clear
% close all
Addpaths
tic
rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);

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
BSAMO_Rate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% Orthogonal matching pursuit
BSAOMP_Rate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% Virtual precoding
VirtualSubarrayRate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% MCCM precoding
MCCMRate = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% DPPTDD precoding Fully connected
DPPTTD8_Full = zeros(length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);
DPPTTD4_Full = zeros(length(RFlist), length(SCSlist), ...
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
                    mMIMO_OFDM_Compare2(OFDMParams, ChanParams, BSParams, UEParams, SNRlist(SNR));

                % Save rate values
                %%% Digital
                DgitalRate(RF, SCS, SNR, SimId) = real(results.DigitalRate);

                %%% Manifold Optimization
                BSAMO_Rate(RF, SCS, SNR, SimId) = real(results.BSAMORate);

                %%% Orthogonal matching pursuit
                BSAOMP_Rate(RF, SCS, SNR, SimId) = real(results.BSAOMPRate);

                %%% Beam broadening
                VirtualSubarrayRate(RF, SCS, SNR, SimId) = real(results.VirtualRate);

                %%% MCCM
                MCCMRate(RF, SCS, SNR, SimId) = real(results.MCCM);

                %%% TTD
                DPPTTD8_Full(RF, SCS, SNR, SimId) = real(results.TTD8_FC);
                DPPTTD4_Full(RF, SCS, SNR, SimId) = real(results.TTD4_FC);

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
    
    %%% DPPTDD (FULLY CONNECTED)
    txt = 'TTD-FC ($N_{TTD}$ = 8)';
    plot(SNRlist, mean(squeeze(DPPTTD8_Full(1, SCS, :, :)), 2), 'DisplayName', txt)

    txt = 'TTD-FC ($N_{TTD}$ = 4)';
    plot(SNRlist, mean(squeeze(DPPTTD4_Full(1, SCS, :, :)), 2), 'DisplayName', txt)

    %%% MCCM
    txt = 'MCCM';
    plot(SNRlist, mean(squeeze(MCCMRate(1, SCS, :, :)), 2), 'DisplayName', txt)

    %%% Beam broadening
    txt = 'Beam broadening';
    plot(SNRlist, mean(squeeze(VirtualSubarrayRate(1, SCS, :, :)), 2), 'DisplayName', txt)
end

legend(Interpreter="latex")
grid
box on
xlabel('SNR [dB]', Interpreter='latex')
ylabel('Spectral efficiency [bps/Hz]', Interpreter='latex')
grid on
box on
hl = legend;
set(hl, 'Interpreter','latex')
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
