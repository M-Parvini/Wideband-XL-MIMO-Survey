clc
clear
close all
tic
rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation code for comparing MO and OMP algorithms with and without beam
% squint effect
%%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
NumSim = 200;
SNRlist = -20:5:10;                        % SNR values
percentage = [0.1 10]./100;
% percentage = 10/100;                       % 10 percent BW/fc
SCSlist = percentage*8e9/64;               % 64 -> number of subcarriers
RFlist = 4;                                % RF chains values
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

% Analog precoding
Analog_Rate = zeros(length(RFlist), length(SCSlist), ...
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
                    mMIMO_OFDM_MOvsOMP(OFDMParams, ChanParams, BSParams, UEParams, SNRlist(SNR));
                % Save rate values
                %%% Digital
                DgitalRate(RF, SCS, SNR, SimId) = real(results.DigitalRate);

                %%% Manifold Optimization
                MO_Rate(RF, SCS, SNR, SimId) = real(results.MORate);
                BSAMO_Rate(RF, SCS, SNR, SimId) = real(results.BSAMORate);

                %%% Orthogonal matching pursuit
                OMP_Rate(RF, SCS, SNR, SimId) = real(results.OMPRate);
                BSAOMP_Rate(RF, SCS, SNR, SimId) = real(results.BSAOMPRate);

                %%% Analog
                % AnalogRate(RF, SCS, SNR, SimId) = real(results.AnalogRate);

            end

        end
    end
end

toc


%% Rate plot with SNR
figure
hold on


%%% Digital
txt = 'Digital';
plot(SNRlist, mean(squeeze(DgitalRate(1, 1, :, :)), 2), 'DisplayName', txt)

%%% BSA MO and OMP
% txt = 'BSA-MO';
% plot(SNRlist, mean(squeeze(BSAMO_Rate(1, SCS, :, :)), 2), 'DisplayName', txt)
% 
% txt = 'BSA-OMP';
% plot(SNRlist, mean(squeeze(BSAOMP_Rate(1, SCS, :, :)), 2), 'DisplayName', txt)

%%% Narrowband
%%% MO
txt = 'MO (Narrowband)';
plot(SNRlist, mean(squeeze(MO_Rate(1, 1, :, :)), 2), 'DisplayName', txt)

%%% OMP
txt = 'OMP (Narrowband)';
plot(SNRlist, mean(squeeze(OMP_Rate(1, 1, :, :)), 2), 'DisplayName', txt)

%%% Wideband
%%% MO
txt = 'MO (Wideband)';
plot(SNRlist, mean(squeeze(MO_Rate(1, 2, :, :)), 2), 'DisplayName', txt)

%%% OMP
txt = 'OMP (Wideband)';
plot(SNRlist, mean(squeeze(OMP_Rate(1, 2, :, :)), 2), 'DisplayName', txt)

%%% Analog
% txt = 'Analog';
% plot(SNRlist, mean(squeeze(AnalogRate(1, SCS, :, :)), 2), 'DisplayName', txt)


legend show
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
