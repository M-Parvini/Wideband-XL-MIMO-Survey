function [P_FC, P_AoSA] = PowerConsumption(Chan, OFDM, BS, UE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% Number of RF chains
NRF = OFDM.RFchain;
% Number of OFDM symbols
Ns = OFDM.numOFDMSym;
% central frequency fc
fc = Chan.fc;
% bandwidth
fs = OFDM.BW;
% max time delay
DelaySpread = Chan.delay_spread;
% OFDM subcarrier numbers
N = OFDM.nfft;
% Number of TTD elements
NTTD = Chan.BSTTDsize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% elements power
PPA = 60e-3;        % power amplifier power
PRF = 26e-3;        % RF chain power
PDAC = 110e-3;      % DAC power
PBB = 200e-3;       % Baseband unit power
PPS = 42e-3;        % Phase shifter power
PTTD = 80e-3;       % TTD elements power

%%% calculating the consumed powers (Tx+Rx-->*2)
Pc = PPA*Nt+NRF*(PRF+PDAC)+PBB;                 % Common power
P_FC = (Pc + PTTD*NTTD*NRF+PPS*Nt*NRF)*2;       % fully connected structure power
P_AoSA = (Pc + PTTD*NTTD*NRF+PPS*Nt)*2;         % partially connected power

