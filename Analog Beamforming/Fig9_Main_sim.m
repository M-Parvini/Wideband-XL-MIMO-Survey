clc
clear
% close all
% rng(1997); % For reprodubility
% poolobj=parpool('HPCServerProfile1',200);
tic
%%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
NumSim = 1e3;
SNRlist = 5:10:55;
% SNRlist = 40;
SNRlistBit = SNRlist - 10*log10(4);
AntennaConfig(1,:) = [1 1];  % [Nt Nr]
% AntennaConfig(2,:)  = [2 1];
% AntennaConfig(3,:)  = [1 2];
% AntennaConfig(4,:)  = [2 2];
ResultsBerTot = zeros([size(AntennaConfig, 1), NumSim, length(SNRlist)]);
ResultsBerMid = zeros([size(AntennaConfig, 1), NumSim, length(SNRlist)]);
ResultsBerOut = zeros([size(AntennaConfig, 1), NumSim, length(SNRlist)]);

for antenna = 1:size(AntennaConfig, 1)
    Nt = AntennaConfig(antenna, 1);
    Nr = AntennaConfig(antenna, 2);
    [OFDMParams, ChanParams, BSParams, UEParams] = ...
            InitializeParams(SNRlist, Nt, Nr);
    parfor SimId = 1:NumSim
        %%% Simulation Cycle
        results = ...
            Massive_MIMO_OFDM(OFDMParams, ChanParams, BSParams, UEParams);
        if mod(SimId,100) == 0
           fprintf(':')
        end
        ResultsBerTot(antenna, SimId, :) = results.BerToT;
        ResultsBerMid(antenna, SimId, :) = results.BerMid;
        ResultsBerOut(antenna, SimId, :) = results.BerOut;
    end
end

%% All the Carriers
figure
for antenna = 1:size(AntennaConfig, 1)
    txt = ['[',num2str(AntennaConfig(antenna,:)),']'];
    p(antenna) = semilogy(SNRlistBit, ...
        mean(squeeze(ResultsBerTot(antenna, :,:)),1), 'DisplayName', txt);
    hold on
end
legend(p(1:end))
grid
title('BER of the all the carriers')
%% Middle Carriers
% figure
% for antenna = 1:size(AntennaConfig, 1)
%     txt = ['[',num2str(AntennaConfig(antenna,:)),']'];
%     p(antenna) = semilogy(SNRlistBit, ...
%         mean(squeeze(ResultsBerMid(antenna, :,:)),1), 'DisplayName', txt);
%     hold on
% end
% legend(p(1:end))
% grid
% title('BER of the 50 percent Middle carriers')
%% Outer Carriers
% figure
% for antenna = 1:size(AntennaConfig, 1)
%     txt = ['[',num2str(AntennaConfig(antenna,:)),']'];
%     p(antenna) = semilogy(SNRlistBit, ...
%         mean(squeeze(ResultsBerOut(antenna, :,:)),1), 'DisplayName', txt);
%     hold on
% end
% legend(p(1:end))
% grid
% title('BER of the 50 percent Outer carriers')
toc
% delete(poolobj)