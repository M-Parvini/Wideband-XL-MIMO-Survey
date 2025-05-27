clc
clear
% close all
tic
rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);
%%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
NumSim = 1;
SNR = 30;
RF = 4;
SCS = 7.5/100*8e9/256;
AntennaConfig = [128 128];  % [Nt Nr]
Nt = AntennaConfig(1, 1);
Nr = AntennaConfig(1, 2);
[OFDMParams, ChanParams, BSParams, UEParams] = ...
        InitializeParams_ArrayGain(SNR, SCS, RF, Nt, Nr);
%%% Simulation Cycle
results = ...
    mMIMO_OFDM_ArrayGain(OFDMParams, ChanParams, BSParams, ...
    UEParams, SNR);

toc
%% All the Carriers
figure
grid
title('Array Gain comparison')
%%% Analog
% txt = 'Analog';
% semilogy(1:OFDMParams.nfft, results.AnalogArrayGain, 'DisplayName', txt)
% hold on

%%% Digital
txt = 'Digital';
semilogy(1:OFDMParams.nfft, results.DigitalArrayGain, 'DisplayName', txt)
hold on

%%% DPP-TTD (Fully connected)
% txt = 'full';
% semilogy(1:OFDMParams.nfft, results.DPPTTDArrayGainFull, 'DisplayName', txt)
% hold on

%%% DPP-TTD (Partially connected)
% txt = 'partial';
% semilogy(1:OFDMParams.nfft, results.DPPTTDArrayGainPartial, 'DisplayName', txt)
% hold on

%%% Virtual subarray
txt = 'Beam broadening';
semilogy(1:OFDMParams.nfft, results.VirtualArrayGain, 'DisplayName', txt)
hold on

%%% MCCM
% txt = 'Mean Channel Covariance';
% semilogy(1:OFDMParams.nfft, results.MCCMArrayGain, 'DisplayName', txt)
% hold on

%%% OMP
txt = 'OMP';
semilogy(1:OFDMParams.nfft, results.OMPArrayGain, 'DisplayName', txt)
hold on

% txt = 'BSA-OMP';
% semilogy(1:OFDMParams.nfft, results.BSAOMPArrayGain, 'DisplayName', txt)
% hold on

%%% MO
txt = 'MO';
semilogy(1:OFDMParams.nfft, results.MOArrayGain, 'DisplayName', txt)
hold on

% txt = 'BSA-MO';
% semilogy(1:OFDMParams.nfft, results.BSAMOArrayGain, 'DisplayName', txt)
% hold on


%%% legends
% legend('Analog', 'Digital', 'TTD', 'Virtual subarray', 'OMP', 'BSA-OMP', ...
    % 'MO', 'BSA-MO')
legend show
grid
xlabel('Subcarrier index [n]')
xlim([1 OFDMParams.nfft])
ylabel('Array Gain')

%%%%%%%%%%%%
% rmpath(added_path);
% rmpath(pwd);