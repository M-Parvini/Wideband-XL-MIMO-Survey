function [Chan, H_f, H_fc] = ChannelModeling(Nt, Nr, fc, fs, K, L, AoD)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit antenna number N
Nt = BS.nAntenna;
% Received antenna number M
Nr = UE.nAntenna;
% central frequency fc
fc = Chan.fc;
% bandwidth
fs = OFDM.BW;
% OFDM subcarrier numbers
K = OFDM.nfft;

tmax = Chan.delay_spread;
% Chan.numPaths = length(pdpPowAligned);
% L = Chan.numPaths;
L = Chan.numClusters;
%%% Channel Modeling
H_f=zeros(Nr,Nt,K);
H_fc=zeros(Nr,Nt,K);
%%% Random Directions for the channel
Alpha_BS = rand(L,1)*pi-pi/2;
Alpha_UE = rand(L,1)*pi-pi/2;
%%% Predefined directions for checking the correctness of the simulation
Alpha_BS = [deg2rad(-75),deg2rad(-45),deg2rad(0),deg2rad(25),deg2rad(45),deg2rad(75)];
Alpha_BS = [deg2rad(-30)];
% Alpha_UE = ones(L,1)*pi/3+pi;
%%% Random gains (Complex-Guassian)
delays = rand(L,1)*tmax*1000;
delays = (1:L)/L*tmax;
delays = 0;
gains = 1/sqrt(2)*(randn(L,1)+randn(L,1)*1j);
gains = 1*ones(L,1);
% delays = pdpDelayAligned;


%%% Channel Frequency Response creation with subcarrier impact
% frequency dependent array response matrix
Atf  = zeros(Nt, L, K);
Atfc = zeros(Nt, L);
Arf  = zeros(Nr, L, K);
Arfc = zeros(Nr, L);

if Chan.LoS
    for k=1:K
        fk = fs/(K)*(k-1-(K-1)/2);
        f=fc+fk;
        %%% LoS path
        % array responses
        Arf(:, 1, k) = ArrayResponse(UE.angle,UE,Chan,Chan.LSpeed/f);
        Atf(:, 1, k) = ArrayResponse(BS.angle,BS,Chan,Chan.LSpeed/f);
        % channel model
        H_f(:,:,k)=H_f(:,:,k)+sqrt(Chan.LoSKfactor)*gains(1)*...
            exp(-1j*2*pi*delays(1)*f)*Arf(:, 1, k)*Atf(:, 1, k)';
        %%% NLoS paths
        for j=2:L
            % array responses
            Arf(:, j, k) = ArrayResponse(Alpha_UE(j),UE,Chan,Chan.LSpeed/f);
            Atf(:, j, k) = ArrayResponse(Alpha_BS(j),BS,Chan,Chan.LSpeed/f);
            % channel model
            H_f(:,:,k)=H_f(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            Arf(:, j, k)*Atf(:, j, k)';
        end        
    end
else
    for k=1:K
        fk = fs/(K)*(k-1-(K-1)/2);
        f=fc+fk;
        %%% NLoS paths
        for j=1:L
            % array responses
            Arf(:, j, k) = ArrayResponse(Alpha_UE(j),UE,Chan,Chan.LSpeed/f);
            Atf(:, j, k) = ArrayResponse(Alpha_BS(j),BS,Chan,Chan.LSpeed/f);
            % Channel model
            H_f(:,:,k)=H_f(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            Arf(:, j, k)*Atf(:, j, k)';
        end        
    end
end

%%% Channel Frequency Response creation for Center frequency (fc)
if Chan.LoS
    for k=1:K
        fk = fs/(K)*(k-1-(K-1)/2);
        f=fc+fk;
        %%% LoS path
        % array responses
        Arfc(:, 1, k) = ArrayResponse(UE.angle,UE,Chan,Chan.LSpeed/fc);
        Atfc(:, 1, k) = ArrayResponse(BS.angle,BS,Chan,Chan.LSpeed/fc);
        % channel model
        H_fc(:,:,k)=H_fc(:,:,k)+sqrt(Chan.LoSKfactor)*gains(1)*...
            exp(-1j*2*pi*delays(1)*f)*Arfc(:, 1, k)*Atfc(:, 1, k)';
        %%% NLoS paths
        for j=2:L
            % array responses
            Arfc(:, j, k) = ArrayResponse(Alpha_UE(j),UE,Chan,Chan.LSpeed/fc);
            Atfc(:, j, k) = ArrayResponse(Alpha_BS(j),BS,Chan,Chan.LSpeed/fc);
            % Channel model
            H_fc(:,:,k)=H_fc(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            Arfc(:, j, k)*Atfc(:, j, k)';
        end        
    end
else
    for k=1:K
        fk = fs/(K)*(k-1-(K-1)/2);
        f=fc+fk;
        %%% NLoS paths
        for j=1:L
            % array responses
            Arfc(:, j, k) = ArrayResponse(Alpha_UE(j),UE,Chan,Chan.LSpeed/fc);
            Atfc(:, j, k) = ArrayResponse(Alpha_BS(j),BS,Chan,Chan.LSpeed/fc);
            % Channel model
            H_fc(:,:,k)=H_fc(:,:,k)+gains(j)*exp(-1j*2*pi*delays(j)*f)*...
            Arfc(:, j, k)*Atfc(:, j, k)';
        end        
    end
end

H_f = H_f.*sqrt(Nt*Nr/L);
H_fc = H_fc.*sqrt(Nt*Nr/L);

%%% Saving the array response
Chan.Atf  = Atf;
Chan.Atfc = Atfc;
Chan.Arf  = Arf;
Chan.Arfc = Arfc;
Chan.normalizedPathGains = gains;

%%% Saving the BS and UE angles
Chan.BSangle = Alpha_BS;
Chan.UEangle = Alpha_UE;

%%% To observe the angular-delay plots
% H_fc = H_f;
%% Delay-Angular plot
HH = dftmtx(400)'*squeeze(H_f)*conj(dftmtx(1024));
figure
if 1
    [X,Y] = meshgrid(1:400,1:1024);
    plott = surf(X,Y,abs(HH).');
    set(plott,'LineStyle','none')

end

end