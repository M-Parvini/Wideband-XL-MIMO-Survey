function [wT, wR] = Analog_precoder(Chan, BS, UE)

%%%
gains = Chan.normalizedPathGains;
[~, idx] = max(gains);
BSAngle = Chan.BSangle;
UEAngle = Chan.UEangle;

BS.angle = BSAngle(idx);
UE.angle = UEAngle(idx);
%%% Steering vector modeling for BS
wT = ArrayResponse(BS.angle,BS,Chan,Chan.lambda);

%%% Steering vector modeling for UE
wR = ArrayResponse(UE.angle,UE,Chan,Chan.lambda);
