clc
clear
% close all
tic
rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Energy Efficiency comparison of TTD networks both for fully and partially
% connected cases.
% Fig 19 (a) -> RFlist = [2 4 8 16]; AntennaConfig = [128, 128];
% Fig 19 (b) -> RFlist = 4, Antennas = [64 128 256 512];
%%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
NumSim = 200;
SNRlist = 10;                         % SNR values
percentage = [1 10];                  % BW [%] [Fig.12(c)]
SCSlist = (percentage/100)*8e9/64;
RFlist = [2 4 8 16];                  % RF chains values Fig 19 (a)
% RFlist = 4;                         % Fig 19 (b)
Antennas = 128;                       % Fig 19 (a)   
% Antennas = [64; 128; 256; 512];     % Fig 19 (b)                  
AntennaConfig = [Antennas Antennas];  % [Nt Nr]

%% SE
% Fully connected TTD
SE_FC = zeros(length(Antennas), length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% Partially connected TTD
SE_PC = zeros(length(Antennas), length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);
%% EE
% Fully connected TTD
EE_FC = zeros(length(Antennas), length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

% Partially connected TTD
EE_PC = zeros(length(Antennas), length(RFlist), length(SCSlist), ...
    length(SNRlist), NumSim);

%%% Main loop
for nt=1:length(Antennas)
    for RF = 1:length(RFlist)
        for SCS = 1:length(SCSlist)
            for SNR = 1:length(SNRlist)
                Nt = AntennaConfig(nt, 1);
                Nr = AntennaConfig(nt, 2);
                [OFDMParams, ChanParams, BSParams, UEParams] = ...
                InitializeParams_TTD(SNRlist(SNR), SCSlist(SCS), RFlist(RF), Nt, Nr);
                fprintf('::::::::::::::::::\n')
                fprintf(['SNR =', num2str(SNRlist(SNR)), '\n'])
                parfor SimId = 1:NumSim
                    %%% Simulation Cycle
                    results = ...
                        mMIMO_OFDM_TDD(OFDMParams, ChanParams, BSParams, UEParams, SNRlist(SNR));

                    % Save rate value                    
                    %%% Fully Connected (SE + EE)
                    SE_FC(nt, RF, SCS, SNR, SimId) = real(results.TDD_FC);
                    EE_FC(nt, RF, SCS, SNR, SimId) = real(results.TDD_FC_EE);

                    %%% Partially Connected (SE + EE)
                    SE_PC(nt, RF, SCS, SNR, SimId) = real(results.TDD_AoSA);
                    EE_PC(nt, RF, SCS, SNR, SimId) = real(results.TDD_AoSA_EE);

                end

            end
        end
    end
end

toc

%% Energy efficiency plot with RF chain
figure
hold on

%%% plots
%%% Partially connected
txt = 'TTD-PC ($B/f_c$ = 1\%)';
plot(1:length(RFlist), mean(squeeze(EE_PC(1, :, 1, 1, :)), 2), 'DisplayName', txt)
txt = 'TTD-PC ($B/f_c$ = 10\%)';
plot(1:length(RFlist), mean(squeeze(EE_PC(1, :, 2, 1, :)), 2), 'DisplayName', txt)
%%% Fully connected
txt = 'TTD-FC ($B/f_c$ = 1\%)';
plot(1:length(RFlist), mean(squeeze(EE_FC(1,:, 1, 1, :)), 2), 'DisplayName', txt)
txt = 'TTD-FC ($B/f_c$ = 10\%)';
plot(1:length(RFlist), mean(squeeze(EE_FC(1,:, 2, 1, :)), 2), 'DisplayName', txt)

legend show
grid
box on

% xticklabels(RFlist)
xlabel('RF chain', Interpreter='latex')
ylabel('Energy efficiency [bps/Hz/W]', Interpreter='latex')
hl = legend;
set(hl, 'Interpreter','latex')
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis


%% Energy efficiency plot with antennas
% figure
% hold on
% 
% %%% plots
% %%% Partially connected
% txt = 'TTD-PC ($B/f_c$ = 1\%)';
% plot(Antennas, mean(squeeze(EE_PC(:,1,1,1,:)), 2), 'DisplayName', txt)
% txt = 'TTD-PC ($B/f_c$ = 10\%)';
% plot(Antennas, mean(squeeze(EE_PC(:,1,2,1,:)), 2), 'DisplayName', txt)
% 
% %%% Fully connected
% txt = 'TTD-FC ($B/f_c$ = 1\%)';
% plot(Antennas, mean(squeeze(EE_FC(:,1,1,1,:)), 2), 'DisplayName', txt)
% txt = 'TTD-FC ($B/f_c$ = 10\%)';
% plot(Antennas, mean(squeeze(EE_FC(:,1,2,1,:)), 2), 'DisplayName', txt)
% 
% 
% legend show
% grid
% box on
% % xticks = 1:length(RFlist);
% % xticklabels(RFlist)
% xlabel('Number of antennas', 'Interpreter','latex')
% ylabel('Energy efficiency [bps/Hz/W]', 'Interpreter','latex')
% hl = legend;
% set(hl, 'Interpreter','latex')
% xaxisproperties= get(gca, 'XAxis');
% xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
% yaxisproperties= get(gca, 'YAxis');
% yaxisproperties.TickLabelInterpreter = 'latex';   % tex for y-axis
