% clc
% clear
% % close all
% % poolobj=parpool('HPCServerProfile1',200);
% rng(1997); % For reprodubility
% %%%%%%%%%%%%%% Parameter Initialization for Simulation %%%%%%%%%%%%%%%%
BW = (0.5:0.5:2)*1e9;
SubList = BW/256;
Ns = [2 4 6];
AntennaType = {'ULA', 'UPA'};
AntennaConfig = [144 16];  % [Nt Nr]
% %%% Rate variables
MCD = zeros(length(AntennaType), length(Ns), length(SubList));

for i = 1:2
    for sub = 1:length(SubList)
        Nt = AntennaConfig(1, 1);
        Nr = AntennaConfig(1, 2);
        SCS = SubList(sub);
        type = AntennaType{i};
        [OFDMParams, ChanParams, BSParams, UEParams] = ...
                InitializeParams(SCS, Nt, Nr, type);

        [Chan, H_f, H_fc] = mmWaveChannelModeling(ChanParams, OFDMParams, ...
            BSParams, UEParams);

        K = OFDMParams.nfft;
        for n=1:length(Ns)
            temp = 0;
            Nsval = Ns(n);
            for k=1:K
                [~,~,V] = svd(H_f(:,:,k));
                [~,~,Vc] = svd(H_fc(:,:,k));
                Vtot = Vc(:,1:Nsval)'*V(:,1:Nsval);
                temp = temp + 1/K*sqrt(Nsval-(norm(Vtot,"fro"))^2);
            end
            MCD(i, n, sub) = temp;
        end
    end
end
%% All the Carriers
figure
hold on
grid
for i=1
    for n=1:length(Ns)
        txt = [AntennaType{i},', $N_s$ = ', num2str(Ns(n))];
        p((i-1)*length(Ns)+n)=plot((0.5:0.5:2),squeeze(MCD(i,n,:)), DisplayName=txt);
    end
end

for i=2
    for n=1:length(Ns)
        txt = [AntennaType{i},', $N_s$ = ', num2str(Ns(n))];
        p((i-1)*length(Ns)+n)=plot((0.5:0.5:2),squeeze(MCD(i,n,:)), '--', DisplayName=txt);
    end
end
legend(p(:), Interpreter="latex")
grid
xlabel('SNR [dB]')
ylabel('Mean chordal distance')
